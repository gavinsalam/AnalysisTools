""" module hfile.py

A set of routines to help read histogram files and convert them to
numpy arrays that can then be manipulated before plotting. This module
also includes some helper routines for manipulating the results
(e.g. rebinning) and searching through the file.

To run unit tests, do

  python3 -m hfile

"""
from builtins import range
from builtins import str
from builtins import object
import copy
import os

# This file is based on an earlier version written by Gavin P. Salam
# (2011)

# Comments may appear excessive in places. This a consequence of this
# file being one of the author's early attempts at python programming.



import numpy as np # for arrays
import string    
import re          # for regular expressions
import io   
#import scipy.interpolate
#import sys

# a few variables that users can set(?)
default_lw = 3
default_border_lw = 1

default_encoding='utf-8'

module_path = __file__
module_directory = os.path.dirname(os.path.realpath(module_path))
datafile=module_directory+"/testing/example.dat"

#----------------------------------------------------------------------
class Error(Exception):
    """Base class for exceptions in this module."""
    def __init__(self, msg):
      self.msg = msg

#----------------------------------------------------------------------
def array(file, regexp=None, fortran=False):
  "Same as get_array"
  return get_array(file, regexp, fortran)

#----------------------------------------------------------------------
def open_any(filename, mode=None):
    """Return either a gzipped file or a normal file, based on the extension"""
    if (len(filename)>3 and filename[-3:] == '.gz'):
        import gzip
        return gzip.GzipFile(filename, mode)
    else:
        return open(filename, mode)

def decode(line):
    if (isinstance(line,bytes)): return line.decode(default_encoding)
    else: return line
    
#----------------------------------------------------------------------
def get_2darray(file,regexp=None):
    """
    For data in the form
    
     1  1   1
     1  2   5
    
     2  1   5
     2  2   8
    
    the resulting array will have the following contents
    
         array[0,0,:] = [1,1,1]
         array[0,1,:] = [1,2,5]
         ...
         array[1,1,:] = [2,2,8]
    
    The array is deemed to terminate when there are two blank lines
    """
    # if file is a string, assume it's a filename, otherwise a filehandle
    if (isinstance(file,str)) : file = open_any(file, 'r')
    if (regexp != None)              : search(file,regexp)
    
    # read in the first slice [0th index = 0]
    slices = []
    slices.append(get_array(file))
    # this will leave us in a situation where we've read in the first blank
    # line that follows a block of numbers
    
    # so we try to read in the next line, and if it's the end of the
    # file, or it's a second blank line, then we deem that we have
    # reached the end of the block, and try to go back to the original
    # place in the file.
    while True:
        # record where we are
        pos = file.tell()
        # read in a line
        line = file.readline()
        #print (line)
        if (not line or line.strip() == ""): break
        # go back to the start of the line
        file.seek(pos)
        slices.append(get_array(file))
    
    result = np.ndarray((len(slices), slices[0].shape[0], slices[0].shape[1]))
    for i in range(len(slices)):
        result[i,:,:] = slices[i]
    
    return result

#----------------------------------------------------------------------
class ArrayPlusComments(object):
    '''
    Object for containing an array and associated information, such as a
    header, a footer and names by which to access columns. Main route
    to creation is from get_array_plus_comments.

    By default, the header and footer include their newlines.

    Notes
    -----
    Ideally this would be simply derived from a numpy array, so that any
    numpy operation can be applied straightforwardly. But for now the 
    array is under the "array" member, and the only numpy-like operation
    that works is indexing.

    To work out how to derive cleanly from an np.ndarray, see 
    https://docs.scipy.org/doc/numpy-1.13.0/user/basics.subclassing.html
    
    '''
    def __init__(self, header, array, footer, columns = None, label = None):
        self.header = header
        self.array  = array
        self.footer = footer
        self.columns = columns
        if (type(columns) is list or type(columns) is tuple):
          new_columns = {}
          for i in range(len(columns)):
            new_columns[columns[i]] = i
          columns = new_columns
        if (columns is not None):
            for (key,col) in columns.items():
                if (isinstance(col, list) or isinstance(col,tuple)) and len(col)==2:
                    labelcol = label+str(col) if label else None
                    setattr(self,key,ValueAndError(array[:,col[0]], array[:,col[1]], label = labelcol))
                else:
                    setattr(self,key,array[:,col])

    def __getitem__(self,indices):
        return self.array[indices]

    def __str__(self):
        return self.header + reformat(self.array) + self.footer

    def __repr__(self):
        return "ArrayPlusComments(header={!r},array={!r},footer={!r},columns={!r})".format(
          self.header, self.array, self.footer, self.columns
        )

    def header_quantity(self,label):
      '''Search the header for a string 'label = ...' and return the ...
      - label is a regexp
      - the search is flexible about spaces around the = sign
      - it returns a float if it can successfully be converted, otherwise a string
      - raises a ValueError if the label can't be found
      '''
      m = re.search(fr'{label}\s*=\s*([^\s]+)', self.header)
      if not m: raise ValueError(f"Could not find header quantity '{label}'")
      try:
        return float(m.group(1))
      except ValueError:
        return m.group(1)
        

#----------------------------------------------------------------------    
def get_array_plus_comments(file, regexp=None, fortran=False, columns = None):
  """Returns an ArrayPlusComments object that contains a header, 2d array and footer. 
  The array is in the same format as 2d array.

  Additionally it is possible to specify column names (columns start
  from 0) by passing the columns argument, e.g.

  >>> r = get_array_plus_comments(datafile, "counts-v-time", columns = {'t':0, 'c':1, 'err':2})
  >>> print(r.t)
  [0. 1. 2.]

  If a column is specified as 'y':[1,2] then the result is returned as a
  ValueAndError object
  """
  # NB: various potential improvements may be possible here, e.g.
  # - using # things like "isspace" to detect blank lines,
  # - using np.loadtxt # to do the string to array conversion,
  # - adding option of getting multiple blocks (e.g. as a python array/list)
  # - and adding option of reading gnuplot-style 3d data
  
  # if file is a string, assume it's a filename, otherwise a filehandle
  header = ''
  closeFile = False
  if (isinstance(file,str)) : 
     file = open_any(file, 'r')
     closeFile = True
  if (regexp != None)       : header = search(file,regexp, return_line = True)

  # handle case where we numbers such as 0.4d3 (just replace d -> e)
  fortranRegex = re.compile(r'd', re.IGNORECASE)
  
  lines = []    # temporary store of lines, before conversion to array
  started = False
  footer = ''
  while True:
    line = decode(file.readline())
    if (not line)               : break        # empty line = end-of-file
    line = line.rstrip()                       # strips trailing blanks, \n
    line = line.lstrip()                       # strips leading blanks
    if (not line) :
      if (started) : break                     # empty line = end-of-block
      else         : continue                  
    # lines that are non-numeric are destined for the header / footer
    if (not re.match('[-0-9]', line)):
        if (not started): header += line+"\n"
        else            : footer += line+"\n"
        continue
    if (fortran): line = re.sub(fortranRegex,'e',line) # handle fortran double-prec
    lines.append(line)                         # collect the line
    started = True
  if closeFile: file.close()
  # do some basic error checking
  if (len(lines) < 1):
    raise Error(f"Block in get_array_plus_comments had 0 useful lines (called with file={file}, regexp={regexp})")
  # now we know the size, transfer the information to a numpy ndarray
  ncol = len(lines[0].split())                
  num_array = np.empty( (len(lines), ncol) )
  for i in range(len(lines)):
    num_array[i,:] = lines[i].split()
  return ArrayPlusComments(header, num_array, footer, columns, label = file.name+" "+str(regexp))
    

#----------------------------------------------------------------------
def get_array(file, regexp=None, fortran=False, regexp_transform = None):
  """
  Returns a 2d array that contains the next block of numbers in this
  file (which can be a filehandle or a filename) 
  
  - if a regexp is provided, then first that is searched for
  - Subsequently, any line starting with a non-numeric character is ignored
  - The block ends the first time a blank line is encountered
  
  For data in the form
  
   1  1
   2  4
   3  9
  
  the resulting array will have the following contents
  
       array[:,0] = [1, 2, 3]
       array[:,1] = [1, 4, 9] 
  
  """
  # NB: various potential improvements may be possible here, e.g.
  # - using # things like "isspace" to detect blank lines,
  # - using np.loadtxt # to do the string to array conversion,
  # - adding option of getting multiple blocks (e.g. as a python array/list)
  # - and adding option of reading gnuplot-style 3d data
  
  # if file is a string, assume it's a filename, otherwise a filehandle
  if (isinstance(file,str)) : file = open_any(file, 'r')
  if (regexp != None)              : search(file,regexp)

  # handle case where we numbers such as 0.4d3 (just replace d -> e)
  fortranRegex = re.compile(r'd', re.IGNORECASE)
  
  lines = []    # temporary store of lines, before conversion to array
  started = False
  while True:
    line = decode(file.readline())
    if (not line)               : break        # empty line = end-of-file
    line = line.rstrip()                       # strips trailing blanks, \n
    line = line.lstrip()                       # strips leading blanks
    if (not line) :
      if (started) : break                     # empty line = end-of-block
      else         : continue                  
    if (not re.match('[-0-9]', line)) : continue
    if (fortran): line = re.sub(fortranRegex,'e',line) # handle fortran double-prec
    lines.append(line)                         # collect the line
    started = True
  # do some basic error checking
  if (len(lines) < 1):
    raise Error("Block in get_array had 0 useful lines")

  # carry out transformation if required
  if (regexp_transform != None):
    for i in range(len(lines)):
      lines[i] = re.sub(regexp_transform[0], regexp_transform[1], lines[i])

  # now we know the size, transfer the information to a numpy ndarray
  ncol = len(lines[0].split())                
  num_array = np.empty( (len(lines), ncol) )
  for i in range(len(lines)):
    num_array[i,:] = lines[i].split()
  return num_array

class Histogram(ArrayPlusComments):
    def __init__(self, apc, plot_args):
        super().__init__(apc.header, apc.array, apc.footer)
        self.name = ''
        self.columns = []
        self.plot_args  = plot_args

        header_lines = self.header.split('\n')
        # 
        self.name = header_lines[0].replace("# ","")
        self.name = re.sub(r' \[.*','',self.name)
        #print(header_lines)
        try:
            columns = [line for line in header_lines[1:] if line.startswith("# cols: ")][0].replace("# cols: ","")
        except IndexError:
            raise RuntimeError("No columns found in header:\n" + self.header)
            # raise ("No columns found in header:\n" + self.header)
            # print(self.header)
            # sys.exit(1)
        columns = columns.replace("# cols: ", "")
        columns = re.sub(r'\(.*','',columns)
        self.columns = columns.split()

    def value_array(self):
        return self.array_by_tag('hist', alt_tag = 'avg')

    def value_column_name(self):
        for column in self.columns:
            if 'hist' in column or 'avg' in column: return column
        raise ValueError("No column with tag 'hist' or 'avg'")

    def error_column_name(self):
        for column in self.columns:
            if 'err' in column: return column
        raise ValueError("No column with tag 'err'")
    
    def x_array(self):
        try:
            return self.array_by_tag('vmid')
        except ValueError:
            return self.array_by_tag('v')

    def has_error(self):
        errcols = [column for column in self.columns if 'err' in column]
        return len(errcols) > 0

    def error_array(self):
        return self.array_by_tag('err')

    def array_by_tag(self,tag, alt_tag = None):
        """Return the array for the first column with the given tag"""
        for i,column in enumerate(self.columns):
            if tag in column: return self.array[:,i]
            if alt_tag and alt_tag in column: return self.array[:,i]
        raise ValueError("No column with tag " + tag)

    def value_or_ValueAndError(self):
        """Return the value if there is no error column, otherwise return a ValueAndError object"""
        if self.has_error():
            return ValueAndError(self.value_array(),self.error_array())
        else:
            return self.value_array()

    def __add__(self, other):
        """sum another histogram to this one; NB does not yet handle total_weight, etc."""
        orig = copy.deepcopy(self)

        if orig.name != other.name:
            raise ValueError(f"Can't sum histograms with different names: {orig.name} and {other.name}")
        if orig.columns != other.columns:
            raise ValueError(f"Can't sum histograms with different columns: {orig.columns} and {other.columns}")
        if (orig.x_array() != other.x_array()).any():
            raise ValueError(f"Can't sum histograms with different x arrays for {orig.name}: {orig.x_array()} and {other.x_array()}")

        # now we will sum other into orig
        orig_contents  = 1.0 * self .value_or_ValueAndError() 
        other_contents = 1.0 * other.value_or_ValueAndError() 
        orig_contents  += other_contents
        if orig.has_error():
            orig.array[:,orig.columns.index(orig.value_column_name())] = orig_contents.value
            orig.array[:,orig.columns.index(orig.error_column_name())] = orig_contents.error
        else:
            orig.array[:,orig.columns.index(orig.value_column_name())] = orig_contents
        return orig

    def __mul__(self, scalar):
        """Multiply this histogram by a scalar, returning a new histogram with the result."""
        orig = copy.deepcopy(self)

        # now we multiply scalar into orig
        orig_contents  = scalar * self.value_or_ValueAndError() 
        if orig.has_error():
            orig.array[:,orig.columns.index(orig.value_column_name())] = orig_contents.value
            orig.array[:,orig.columns.index(orig.error_column_name())] = orig_contents.error
        else:
            orig.array[:,orig.columns.index(orig.value_column_name())] = orig_contents
        return orig

    def __rmul__(self, scalar):
        return self.__mul__(scalar)
    
    def __rtruediv__(self, scalar):
        orig = copy.deepcopy(self)

        # now we operate on the value or ValueAndError
        orig_contents  = scalar / self.value_or_ValueAndError() 
        if orig.has_error():
            orig.array[:,orig.columns.index(orig.value_column_name())] = orig_contents.value
            orig.array[:,orig.columns.index(orig.error_column_name())] = orig_contents.error
        else:
            orig.array[:,orig.columns.index(orig.value_column_name())] = orig_contents
        return orig
       

    def plot_to_axes(self, ax, norm = None, **extra):
        """Plot the histogram to the given axes, possibly normalised by 
        the histogram specified by the norm argument (that histogram's error is 
        not included in the error estimate of the ratio)"""
        # multiply by 1.0 to make sure that we take a copy before
        # subsequent normalisation
        contents = 1.0 * self.value_or_ValueAndError() 
        # get arguments, with priority given to extra
        combined_args = self.plot_args.copy()
        combined_args.update(extra)

        # divide by the value, not by the ValueAndError object, to
        # avoid issues when the self and norm are the same thing
        if norm is not None: contents /= norm.value_array()
        if self.has_error():
            line_and_band(ax,self.x_array(), contents, **combined_args)
        else:
            ax.plot(self.x_array(), contents, **combined_args)

    def set_axes_data(self, ax):
        ax.set_title(self.name)
        ax.set_xlabel(re.sub(r'.*:','',self.name))
        ax.set_ylabel(self.value_column_name())

    def plot(self, pdf, others = []):
        import matplotlib.pyplot as plt

        print("Plotting histogram", self.name)
        fig,ax = plt.subplots()

        self.set_axes_data(ax)
        self.plot_to_axes(ax, **styles[0])

        pdf.savefig(fig,bbox_inches='tight')
        plt.close()

        # ax.set_xlim(0,2.5)
        # ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        # ax.yaxis.set_minor_locator(MultipleLocator(0.02))
        # ax.tick_params(top=True,right=True,direction='in',which='both')
        # ax.set_xscale('log')
        # ax.xaxis.set_major_formatter(FuncFormatter(h.log_formatter_fn))
        # ax.grid(True,ls=":")
        # ax.set_title("title")
        # ax.text(x,y,'hello',transform=ax.transAxes)
        # ax.plot(res.x, res.y, label='label', **styles[0])
        #ax.legend(loc='upper left')
        #pdf.savefig(fig,bbox_inches=Bbox.from_extents(0.0,0.0,7.5,4.8))


class HFile(object):
    def __init__(self,filename):
        '''Read in an HFile object from the file with the given name.
        The filename can include labels for the histograms, e.g.
        "file.dat:label=default" will label the histograms
        from this file as "default".
        '''
        has_labels = len(filename.split(':')) > 1
        self.filename = filename.split(':')[0]
        self.plot_args = {}
        if(has_labels):
            for f in filename.split(':')[1].split(" "):
                key, value = f.split("=")
                self.plot_args[key] = value
        else:
            self.plot_args['label'] = self.filename
        self.header = ''
        self.warnings = ''
        self.histograms = []

        f = open_any(self.filename,'r')

        # first get the header
        while True:
            pos = f.tell()
            line = f.readline()
            if line == '': break
            if re.search(r'hist.*:',line): 
                f.seek(pos)
                break
            self.header += line

        # then get the histograms, bailing out if we see a WARNING SUMMARY
        while True:
            pos = f.tell()
            line = f.readline()
            if line == '': break
            #print(line.strip())
            if (line.strip() == ''): continue
            if 'WARNING SUMMARY' in line:
                self.warnings += line
                break
            if 'hist' in line:
                f.seek(pos)
                histogram = Histogram(get_array_plus_comments(f), self.plot_args)

                #print(histogram.name, histogram.columns)
                self.histograms.append(histogram)

        self.map = {}
        for hist in self.histograms:
            self.map[hist.name] = hist

        # finally the warnings
        if self.warnings:
            for line in f: self.warnings += line

    def by_name(self,name):
        """Return the histogram with the exact given name"""
        return self.map[name]
    
    def by_re(self,regexp):
        """Return the histogram with a name that matches the given regexp"""
        hists = []
        for hist in self.histograms:
            if re.search(regexp,hist.name): hists.append(hist)
        if len(hists) == 1: return hists[0]
        if len(hists) == 0: raise ValueError(f"No histogram with name matching '{regexp}'")
        if len(hists) > 1:  
           raise ValueError(f"Multiple histograms with name matching '{regexp}':"
                            +", ".join([hist.name for hist in hists]))

    def __add__(self, other):
        new = copy.deepcopy(self)
        for ih, h in enumerate(new.histograms):
            new.histograms[ih] += other.histograms[ih]
        new.filename += ' + ' + other.filename
        if(other.plot_args['label'] == other.filename):
            new.plot_args['label']   = new.filename
            for h in new.histograms:
                h.plot_args  = new.plot_args
        return new

def line_and_band(ax,x,val_and_err,**extra):
    extra_no_label = copy.copy(extra)
    if ('label' in extra_no_label): del extra_no_label['label']
    ax.fill_between(x,
                    val_and_err.value-val_and_err.error,
                    val_and_err.value+val_and_err.error,
                    alpha=0.2,
                    **extra_no_label
                    )
    ax.plot(x,val_and_err.value, **extra)

#----------------------------------------------------------------------
class XSection(object):
    """\
    Class that contains the result of an "xsc" entry in the file.
    (This is intended to work with files from the CSS framework)
    
    Relevant members are:
    obj.xsc : the cross section
    obj.err : its error
    obj.nentries: the number of entries of which it is built
    obj.name: the name
    obj.line: the full line
    obj.ve(): returns a ValueAndError object
    """
    def __init__(self, line, xsc_label='xsc'):
        self.line = line
        # got this from
        # http://stackoverflow.com/questions/6508043/regular-expression-to-find-any-number-in-a-string
        number_regex = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        # try to get rid of everything to the left of the first equals (with spaces)
        lineRH = re.sub(r"^.*? = ","", line)
        # use the first three numbers of the part beyond the first " = "
        numbers = number_regex.findall(lineRH)
        self.xsc = float(numbers[0])
        self.err = float(numbers[1])
        self.nentries = float(numbers[2])
        # recall that "group" gets the first matching group here...
        if xsc_label != "":
          self.units = re.search(r"([a-z]+) \(n entries",line).group((1))
          self.name = re.search(rf"# (.*): {xsc_label} =",line).group((1))
        else: 
          self.units = None
          self.name = re.search(rf"# (.*) =",line).group((1))

    def ve(self):
       """return the value and error as a tuple as a ValueAndError object"""
       return ValueAndError(self.xsc, self.err)

    def __str__(self):
        return "{0} +- {1} {2} (n entries = {3})".format(self.xsc,self.err,self.units,self.nentries)

#----------------------------------------------------------------------
def get_xsection(file,regexp,xsc_label='xsc'):
    """return a cross section value in the specified file, matching the 
       regexp and ending in xsc_label.
    """
    if (isinstance(file,str)) : file = open_any(file, 'r')
    fullregexp = regexp+rf'.*{xsc_label} ='
    line = search(file, fullregexp, return_line = True)
    #
    result = XSection(line, xsc_label)
    return result

#----------------------------------------------------------------------
def search(filehandle, regexp, return_line=False):
  """ looks through the file described by handle until it finds the regexp
  that's mentioned"""
  regexp_compiled = re.compile(regexp)
  while True:
    line = decode(filehandle.readline())
    if (line == "") : break        # empty line = end-of-file
    #if (re.search(regexp,line)) : 
    if (regexp_compiled.search(line)) : 
        if (return_line): return line
        else            : return filehandle
  return None


#----------------------------------------------------------------------
def reformat(*columns, **keyw):
  """returns a string containing each of the columns placed
      side by side. If instead of columns, 2d numpy arrays
      are supplied, then these are output sensibly too

      For now, it assumes that columns are numpy arrays,
      but this could in principle be relaxed.

      Among the keyword arguments, the only one currently supported is
      "format", which should be a string used to format the output

  """

  ncol = len(columns)
  shapes = []
  ncols  = []
  for i in range(ncol):
    shapes.append(columns[i].shape)
    if    (len(shapes[i]) == 1): ncols.append(0)
    elif  (len(shapes[i]) == 2): ncols.append(shapes[i][1])
    else: raise Error(" a 'column' appears not to be 1 or 2-dimensional" )

  # lazily assume that all lengths are the same
  nlines = shapes[0][0]

  #output = io.BytesIO()
  output = io.StringIO()
  if ("format" in keyw):
    frm=keyw["format"]
    
    for i in range(nlines) :
      for j in range(ncol):
        if (ncols[j] == 0): print(frm.format(columns[j][i]), end=' ', file=output)    # trailing comma kills newline
        else: 
          for k in range (ncols[j]):
            print(frm.format(columns[j][i,k]), end=' ', file=output)
      print(file=output)
  else:
    for i in range(nlines) :
      for j in range(ncol):
        if (ncols[j] == 0): print(columns[j][i], end=' ', file=output)    # trailing comma kills newline
        else: 
          for k in range (ncols[j]):
            print(columns[j][i,k], end=' ', file=output)
      print(file=output)

      
  return output.getvalue()

#----------------------------------------------------------------------
def cumul(array, nx=3, renorm=False, reverse = False):
  """
  Returns the accumulated version of the array. Currently assumes uniform
  spacing.

  @param nx      indicates how many columns of x-axis info there are. If
                 nx = 1, the spacing must be uniform and the x axis is adjusted
                 so as to be shifted by half a bin in the right direction

  @param renorm  if True, it assumes that we are dealing with a 
                 differential distribution (so sum of bins must be divided 
                 by rebin); NB: this is buggy if bin spacings not uniform

  @param reverse accumulates negatively from the high end downwards if
                 true. The high bin will always be at xmax+dx/2,
                 regardless of the value of reverse, which means that
                 if reverse is True, then the high-end y bins will be
                 zero. This facilitates addition of results with and
                 without reverse, since the x bins should always match
  """
  ny = len(array[:,0])
  if (ny < 2) : raise Error("number of bins must be >= 2")

  result = np.empty(array.shape)
  # copy x and y axes
  result[:,:nx] = array[:,:nx]
  if (reverse):
    # some care is needed to make sure that we align the bins properly
    result[ny-1,nx:] = 0
    result[:ny-1,nx:] = -array[ny-1:0:-1,nx:].cumsum(axis=0)[::-1]
  else:
    result[:,nx:] = array[:,nx:].cumsum(axis=0)

  # sort out shifting of bins and renormalisation
  if (nx == 1): 
    dx = array[1,0]-array[0,0]
    result[:,0] += 0.5*dx
    if (renorm): result[:,nx:] *= dx  
  else:
    if (renorm): result[:,nx:] *= (result[:,nx-1]-result[:,0]).reshape(ny,1)
  
  # finally handle the reverse accumulation, keeping identical binning
  # to what we would have had in the forward direction (this
  # effectively loses us one bin of information, but it adds
  # convenience for adding things produced in different directions)
  #if (reverse): result[:,nx:] = -result[-1,nx:] + result[:,nx:]
    
  return result

#----------------------------------------------------------------------
def rebin(array, rebin=2, nx=3, renorm=False, 
                 min=-1e300, max=1e300):
  """ Returns a rebinned version of an array. 
  
  @param rebin   indicates how many bins to combine
  @param nx      indicates how many columns of x-axis info there are
  @param renorm  if True, it assumes that we are dealing with a 
                 differential distribution (so sum of bins must be divided 
                 by rebin)
  @param min     starts only from bins >= min
  @param max     goes only up to bins <= max
  """
  ilo = 0
  ihi = len(array[:,0])-1

  if (ihi < 1) : raise Error("number of bins must be >= 2")
  
  if   (nx == 1) : minslice = array[:,0]; maxslice = array[:,0]
  elif (nx == 2) : minslice = array[:,0]; maxslice = array[:,1]
  elif (nx == 3) : minslice = array[:,0]; maxslice = array[:,2]
  else : raise Error("nx is not 1, 2 or 3")
  
  # get the bounds
  for i in range(len(minslice)):
    if (min <= minslice[i]): ilo = i; break;
  for i in range(len(minslice)-1,-1,-1):
    if (max >= maxslice[i]): ihi = i; break;

  # now figure out how many bins we're going to need
  nbins = int(ihi+1-ilo)//rebin
  #print ilo, ihi, minslice[ilo], maxslice[ihi]
  #print nbins
  rebinned = np.empty( (nbins, len(array[0,:])) )
  
  # and now do the rebinning
  ir = 0
  for i in range(ilo, ilo+nbins*rebin, rebin):
    # first handle the x axes
    if   (nx == 1): rebinned[ir, 0] = sum(array[i:i+rebin,0])/rebin
    elif (nx >= 2): 
      rebinned[ir, 0   ] = array[i        , 0    ]
      rebinned[ir, nx-1] = array[i+rebin-1, nx-1 ]
      if (nx == 3) : rebinned[ir, 1] = 0.5*(rebinned[ir,0]+rebinned[ir,2])
    # then the remaining ones
    rebinned[ir, nx:] = np.sum(array[i:i+rebin, nx:], axis=0)
    ir += 1

  if (renorm) : rebinned[:,nx:] /= rebin
  return rebinned

#----------------------------------------------------------------------
def bins(array):
  """Returns the bin edges of the array, assuming the usual 3-column format of xlo xmid xhi,
  which can be useful with matplotlib's axes.hist. E.g. use as

    #        bin-mid    bin edges      bin contents
    ax.hist(array[:,1], h.bins(array), weights=array[:,3])

  or

    ax.stairs(array[:,3], h.bins(array))

  """
  result = np.empty(len(array[:,0])+1)
  result[0:-1] = array[ :, 0]
  result[-1]   = array[-1, 2]
  return result

#----------------------------------------------------------------------
def extended_binary(dict_or_list_of_arrays, binary_fn, subset = None):
  """
  Return an array each of whose elements is given by the repeated
  pairwise application of binary_fn on the corresponding elements in the
  dict_or_list_of_arrays; if subset is not None, then only the keys in
  the subset are used
  """
  
  if (isinstance(dict_or_list_of_arrays, dict)):
    if (subset is None): subset = list(dict_or_list_of_arrays.keys())
    arrays = [dict_or_list_of_arrays[key] for key in subset]
  else:
    arrays = dict_or_list_of_arrays

  first = True
  for array in arrays:
    if (first): result = array
    else      : result = binary_fn(result,array)
    first = False

  return result

#----------------------------------------------------------------------
def maximum(dict_or_list_of_arrays, subset = None):
  """
  Return an array each of whose elements is given by the maximum of
  the corresponding elements in the dict_or_list_of_arrays; if subset is not
  None, then only the keys in the subset are used
  """
  return extended_binary(dict_or_list_of_arrays, np.maximum, subset = subset)

#----------------------------------------------------------------------
def minimum(dict_or_list_of_arrays, subset = None):
  """
  Return an array each of whose elements is given by the minimum of
  the corresponding elements in the dict_or_list_of_arrays; if subset is not
  None, then only the keys in the subset are used
  """
  return extended_binary(dict_or_list_of_arrays, np.minimum, subset = subset)


#----------------------------------------------------------------------
def minmax(dict_or_list_of_arrays, subset = None):
  """ Return a tuple such that [0][...] is the min
  and [1][...] is the max across the dictionary of arrays
  """
  return (minimum(dict_or_list_of_arrays,subset),
          maximum(dict_or_list_of_arrays,subset))

#------------------------------------------------------------
# should go into a base class?
def opt_or_default(keyw, key, default = None, delete=False):
  """
  If keyw[key] exists, returns it and deletes keyw[key], otherwise
  returns the default value
  """
  if key in keyw:
    result = keyw[key]
    if delete: del keyw[key]
    return result
  else: return default

#-------------------------------------------------------------
def index_of_value(arr, value, tolerance=1e-7):
  '''
  Returns the index at which the array is close to the specified value,
  requiring that the result be within the specified tolerance.
  
  Current implementation takes O(N) time
  '''
  index = np.argmin(np.abs(arr-value))
  if np.abs(arr[index]-value) > tolerance: 
    raise ValueError(f"Could not locate value={value} within array to "
                     f"within tolerance={tolerance}, closest was {arr[index]} at index {index}")
  else:
    return index

#----------------------------------------------------
def log_formatter_fn(value, order):
    """A formatter for matplotlib that uses normal notation for numbers 
    between lower and upper (by default 0.1, 1.0 and 10.0)
    and LaTeX scientific notation otherwise. The "order" is a dummy argument
    required by matplotlib but not used here. 
    
    To use it do

    from matplotlib.ticker import FuncFormatter
    # [...]
    ax.xaxis.set_major_formatter(FuncFormatter(hfile.log_formatter_fn))
    """
    epsilon = 1e-8
    lower=0.09
    upper=11.0

    if value > upper or value < lower:
        exponent = int(np.floor(np.log(value)/np.log(10.0) + epsilon))
        mantissa = value / 10**exponent
        mantissa_is_one = np.abs(mantissa - 1.0) < epsilon
        mantissa_rint = np.rint(mantissa)
        mantissa_is_whole = np.abs(mantissa_rint - mantissa) < epsilon
        if   mantissa_is_one:   return f"$10^{{{exponent}}}$"
        elif mantissa_is_whole: return f"${int(mantissa_rint)} \times 10^{{{exponent}}}$"
        else:                   return f"${mantissa} \times 10^{{{exponent}}}$"
    else: 
        value_rint = np.rint(value)
        value_is_whole = np.abs(value_rint - value) < epsilon
        if value_is_whole: return f"{int(value_rint)}"
        else:              return f"{value}"

#-------------------------------------------------------------
class ValueAndError(object):
    '''
    class that represents a value with its error.

    Arithmetic should be handled smoothly, both with ValueAndError objects
    and with normal arithmetic quantities. 

    The object can also be initialised with numpy arrays and operations 
    then apply element by element.  

    For example
    >>> a = ValueAndError(1.0,0.3)
    >>> b = ValueAndError(2.0,0.4)
    >>> print(a+b)
    3.0 ± 0.5
    '''
    def __init__(self,value,error,label=None):
        '''Create an object from an input value and error, which can be
        plain numbers or numpy arrays; 
        
        Parameters
        ----------

        value : float or numpy array
          the value
        
        error : float or numpy array
          the error on the value

        label        
          if provided this is used as the label (can also be
          a set of labels), if not then the object gets a unique label
          based on id(self); labels are used for tracking correlations
        
        '''      
        self.value = value
        self.error = np.abs(error)
        # if the label is a set, then we assume it's a set of labels
        if isinstance(label,set):
          self.labels = label
        else:
          self.labels = {label} if label else {id(self)}

    def __rmul__(self,fact):
        return self*fact

    def __mul__(self,fact):
        if (isinstance(fact,ValueAndError)):
            self.check_correlation(fact)
            prod = self.value*fact.value
            err  = np.sqrt(self.value**2 * fact.error**2 + fact.value**2 * self.error**2)
            return ValueAndError(prod, err , self.labels | fact.labels)
        else:
            return ValueAndError(self.value*fact, self.error*np.abs(fact), self.labels)
        
    def __truediv__(self,fact):
        if (isinstance(fact,ValueAndError)):
            return self.__mul__(fact.inverse())
        else:
            return ValueAndError(self.value/fact, self.error/np.abs(fact), self.labels)

    def __rtruediv__(self,fact):
        return fact * self.inverse()
        
    def __add__(self,other):
        if (isinstance(other,ValueAndError)):
            self.check_correlation(other)
            return ValueAndError(self.value+other.value, np.sqrt(self.error**2+other.error**2), self.labels | other.labels)
        else:
            return ValueAndError(self.value+other, self.error, self.labels)
        
    def __radd__(self,other):
        if (isinstance(other,ValueAndError)):
            raise TypeError('radd should never have ValueAndError as "other"')
        else:
            return ValueAndError(self.value+other, self.error, self.labels)

    def __sub__(self,other):
        if (isinstance(other,ValueAndError)):
            self.check_correlation(other)
            return ValueAndError(self.value-other.value, np.sqrt(self.error**2+other.error**2), self.labels | other.labels)
        else:
            return ValueAndError(self.value-other, self.error, self.labels)

    def __rsub__(self,other):
        if (isinstance(other,ValueAndError)):
            raise TypeError('rsub should never have ValueAndError as "other"')
        else:
            return ValueAndError(other - self.value, self.error, self.labels)

    def __pow__(self,power):
        return ValueAndError(self.value**power, np.abs(power * self.value**(power-1)) * self.error, self.labels)
        
    def __neg__(self):
        return ValueAndError(-self.value, self.error, self.labels)
    
    def __repr__(self):
        return f'{self.value} ± {self.error}'

    def __format__(self,format_spec):
        return (f"{{:{format_spec}}} ± {{:{format_spec}}}").format(self.value, self.error)

    def __getitem__(self, *args):
        '''Returns a ValueAndError object with the subscripting applied 
        to individual elements of the value and the error.

        For use, e.g., when both are numpy arrays. NB, cannot
        currently handle the case when one is a numpy array while the
        other is a scalar.
        '''
        return ValueAndError(self.value.__getitem__(*args), self.error.__getitem__(*args))
    
    def inverse(self):
        """Return 1.0/self, with the error propagated

        For example
        >>> a = ValueAndError(2.0,0.1)
        >>> print(a.inverse())
        0.5 ± 0.025
        """
        return ValueAndError(1.0/self.value, self.error/self.value**2, self.labels)

    def is_correlated(self, other):
        """
        Return True if the ValueAndError objects have any labels in common

        For example:
        >>> a = ValueAndError(1.0,0.1)
        >>> b = ValueAndError(2.0,0.2)
        >>> a.is_correlated(b)
        False
        >>> c = a*b
        >>> c.is_correlated(a)
        True
        """
        if self.labels.intersection(other.labels): return True
        else: return False
    
    def check_correlation(self, other):
        if self.is_correlated(other):
            import warnings
            warnings.warn(f"Warning, applying binary operator on ValueAndError objects that are correlated, "
                          f"with labels {self.labels} and {other.labels};\n"
                          f"*** Run script as `python3 -W error scriptname.py` to get a traceback ***")
            #import traceback
            #traceback.print_stack()

def unit_tests():
    """Unit tests that cannot easily be included as part of the doctest"""
    import warnings
    # convert warnings to errors, so that we can test for them
    warnings.simplefilter("error")
    # check for warnings on correlated operations
    a = ValueAndError(1.0,0.1,'a')
    b = ValueAndError(2.0,0.2,'b')
    ops = ["*", "/", "+", "-"]
    for op1 in ops:
       for op2 in ops:
           for template in ("c = a {op1} b; d = c {op2} a", #< check c, a correlation
                            "c = a {op1} b; d = c {op2} b", #< check c, b correlation
                            "a {op1}= b; a {op2}= a",       #< check a, a correlation
                            "a {op1}= b; a {op2}= b",       #< check a, b correlation
                            "c = a {op1} b; a {op1}= 0.3; a {op2}= c" #< check correlations aren't lost with scalar operations
                            ):
            try:
                instr = template.format(op1=op1,op2=op2)
                #print(instr, expect_warning)
                exec(instr)
                raise ValueError(f"Operation `{instr}` should have raised a UserWarning")
            except UserWarning:
                pass

    #for operation in ("d = a+c", "d = a/c", "d = a-c", "d = a*c", "a -= c", "a += c", "a *= c", "a /= c"):
    #    try:
    #        a -= c
    #        print(operation)
    #        exec(operation)
    #        raise ValueError(f"Operation {operation} should have raised a UserWarning")
    #    except UserWarning:
    #        pass
        
    #counts_v_time = get_array_plus_comments("testing/example.dat", "histogram:counts-v-time", columns={'time':0, 'count':[1,2]})
    #print(counts_v_time.count, counts_v_time.count.labels)

if __name__ == "__main__": 
  import doctest
  nfailures, ntests = doctest.testmod()
  if nfailures == 0:
      print(f"All {ntests} doctests passed")
  else:
      print(f"{nfailures} doctests failed out of {ntests}")
      exit(1)
  unit_tests()
