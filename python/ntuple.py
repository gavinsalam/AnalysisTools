#!/usr/bin/env python
"""Small library to read the ntuples produced by the generator
framework's SimpleNtuple class

Ntuple is data is in a plain text format (assumed compressed if the filename ends in gz)

  - anything starting with a "#" is a comment
  - anything starting with "!!" indicates a line of ntuple entries (one event)
  - each event has individual entries of the form name#value

Usage:

  - create an NtupleReader

  - do "for event in reader:" to iterator over events

  - then for each event you every ntuple element present as a member
    of event (with some name rewriting, e.g. "." -> "__")

TO DO:
  - get an iterable structure so that one can say
  
    for event in reader:
        ....
  
  - make it command-line executable, e.g. to process ntuple files on the fly,
    extract information from ntuple files (such as list of keys), etc.

  - incorporate other facilities (such as handling of cross sections, etc.)

"""

from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import object
import re
import gzip
import sys

#======================================================================
def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-l","--list-entries", action='store_true', help='list the entries (keys) that are present')
    parser.add_argument("--nev", default=1000, type=float, help="default number of events to process")
    parser.add_argument("--show", metavar="key1[,key2[,...]]", help="show the values of the keys, provided as a comma-separated list")
    parser.add_argument("--select", default=None, help="condition that must be satisfied in order to show/average something; if required elements from the condition are missing, the condition is considered to be false")
    parser.add_argument("--out","-o",default=sys.stdout,type=argparse.FileType('w'),
                        help='output file, default = stdout')
    parser.add_argument("filename")
    args = parser.parse_args()
    
    print("# "+" ".join(sys.argv), file=args.out)
    print("#",args, file=args.out)

    if (args.select is not None):
        select = return_selector(args.select)
    else:
        def select(ev): return True
        
    #----------------------------------------------------------------------
    if (args.list_entries):
        reader = NtupleReader(args.filename, args.nev)
        entries = {}
        for ev in reader:
            for element in ev.__dict__:
                if (element not in entries): entries[element] = 1
                else                       : entries[element] += 1

        sorted_keys = sorted(entries.keys())
        frmt_str="{:30} {:>8}"
        print(frmt_str.format("ntuple key","n occ."), file=args.out)
        print("-----------------------------------------")
        for k in sorted_keys:
            print(frmt_str.format(k, entries[k]), file=args.out)
    elif (args.show):
        reader = NtupleReader(args.filename, args.nev)
        keys = args.show.split(",")
        print("# "+" ".join(keys), file=args.out)
        for ev in reader:
            try:
                if (not select(ev)): continue
            except AttributeError: continue
            out = ""
            for k in keys:
                try:
                    out += str(ev[k]) + " "
                except KeyError:
                    out += "None "
            print(out, file=args.out)
            
#----------------------------------------------------------------------
class NtupleReader(object):
    def __init__(self, filename, nev_max = None):
        if (type(filename) is list or type(filename) is tuple):
            self._n_files = len(filename)
            self._i_file  = 0
            self._filenames = filename
            self.set_handle_to_filename(self._filenames[0])
        else:
            self._filenames = [filename]
            self._n_files = 1
            self._i_file  = 0

            self.set_handle_to_filename(filename)

        self.comments = ""
        self._iev = 0
        self._nev_max = nev_max
        # keep track of the number of generated events through the weight lines
        self._last_iev_with_weight_line = 0
        self._last_iev_with_weight_line_monotonic = True

        # for replica studies it's useful to be able to reset the lumi
        self._last_iev_reset_lumi = 0

    #----------------------------------------------------------------------
    def set_handle_to_filename(self,filename):
        self.filename = filename
        print('Reading ntuples from', self.filename, file=sys.stderr)
        if (re.search(r'.gz$', filename)):
            self.handle = gzip.GzipFile(filename, 'r')
        else:
            self.handle = open(filename, 'r')

    #----------------------------------------------------------------------
    def __iter__(self):
        # needed for iteration to work 
        return self

    #----------------------------------------------------------------------
    # Python 3 compatibility
    def __next__(self):
        return self.next()

    #----------------------------------------------------------------------
    # Python 2 version
    def next(self):
        ev = self.next_event()
        if (ev is None):
            raise StopIteration
        else:
            return ev

            
    #----------------------------------------------------------------------
    def next_event(self):
        '''get the next line with "data"; return None if it's not available; 
        NB, it's probably better to use the iterators'''

        # first check we haven't gone beyond the end of the requested
        # number of events
        if (self._nev_max and self._iev >= self._nev_max): return None
        
        while(True):
            try:
                line = self.handle.readline().decode("utf-8")
            except (IOError, EOFError):
                # if we encounter an error, handle it gracefully and try next file
                print(' --> encountered an IOError or EOFError on file', self.filename, file=sys.stderr)
                if (self._next_file()): continue
                else                  : return None

            # if there is no line (EOF may not always raise an error above), 
            # try the next file
            if (not line):
                if (self._next_file()): continue
                else                  : return None

            if (len(line) > 2):
                if (line[0:2] == '!!'): break
                if (line[0:2] == '%%'):
                    self.register_normalisation(line)
                    continue
            # otherwise assume we have a comment line
            self.comments += "#"+line
            
        self._iev += 1
        self.event = EventNtuple(line[2:].split())
        return self.event

    #----------------------------------------------------------------------
    def _next_file(self): 
        """
        move on to next file; return True if successful, False otherwise
        """
        # when we come to the end of a file, see if there's another file
        self._i_file += 1
        if (self._i_file  >= self._n_files):
            self._i_file -= 1 # correct it
            return False
        else:
            self.set_handle_to_filename(self._filenames[self._i_file])
            return True

    #----------------------------------------------------------------------
    def iev(self): 
        """
        Returns the index of the current event (starts from 1 and equal to
        the total number of events read in)
        """
        return self._iev

    #----------------------------------------------------------------------
    def nev_generated(self):
        return (self._iev * self._iev_ratio_truth_to_recorded)
        
    #----------------------------------------------------------------------
    def register_normalisation(self,line):
        '''take a normalisation line, validate and store its info'''

        # first convert the normalisation line (same format as
        # ntuples, except for leaving %% instead of !!
        norm_info = EventNtuple(line[2:].split())
        
        # the total cross section is
        #
        #  (\sum event_weights) * weight_factor_per_event_nb / iev
        #
        self._weight_factor_per_event_nb = norm_info.weight_factor_per_event_nb
        self._total_cross_section_nb = norm_info.total_xsc_nb
        
        # need to handle cases where only some events are stored in the
        # ntuple; here deduce the ratio of generated events to stored events
        #
        # also aim to be careful in cases where several ntuple files
        # are concatenated; the convention is that a cross section
        # line is always written out for iev=1. So if we detect a
        # cross section line with an iev value that is smaller than
        # the one seen so far, we assume that we we should just bail
        # out of trying to deduce anything more about event ratios
        if (self._last_iev_with_weight_line_monotonic):
            if (norm_info.iev > self._last_iev_with_weight_line):
                self._last_iev_with_weight_line = norm_info.iev
                if (self._iev > 1):
                    self._iev_ratio_truth_to_recorded = (norm_info.iev-1.0)/self._iev
                else:
                    self._iev_ratio_truth_to_recorded = 1.0
            else:
                self._last_iev_with_weight_line_monotonic = False
        #self._iev_ratio_truth_to_recorded = (norm_info.iev-1.0)/self._iev

    #----------------------------------------------------------------------
    def cross_section_nb(self, weight_sum = None):
        """
        Returns the cross section in nb. 

        If weight_sum is not provided, then this is the total cross section.
        Otherwise it's the cross section corresponding to the weight_sum that
        is supplied.
        """
        # use a try/except, just in case the information is missing
        try:
            if (not weight_sum):
                sigma = self._total_cross_section_nb
            else:
                sigma = weight_sum \
                        * self._weight_factor_per_event_nb / \
                        (self._iev * self._iev_ratio_truth_to_recorded)
        except AttributeError:
            sigma = float('nan')
        return sigma

    #----------------------------------------------------------------------
    def lumi_invnb(self):
        """Returns the total lumi used, in invnb"""
        return (self._iev * self._iev_ratio_truth_to_recorded) / self._total_cross_section_nb

    #----------------------------------------------------------------------
    def lumi_invnb_since_reset(self):
        """Returns the total lumi used, in invnb"""
        return ((self._iev - self._last_iev_reset_lumi)
                * self._iev_ratio_truth_to_recorded) / self._total_cross_section_nb
    
    #----------------------------------------------------------------------
    def reset_lumi(self):
        """By calling this, the luminosity counter is reset to zero, which can
        be useful when carrying out Monte Carlo replica studies

        """
        self._last_iev_reset_lumi = self._iev
        
    #----------------------------------------------------------------------
    def __str__(self):
        """Returns a header with information about the ntuple"""
        output = ""
        output += "# NtupleReader from {}\n".format(self._filenames)
        output += "#     processed {} events from {} of {} files\n".format(self._iev, self._i_file+1, self._n_files)
        try:
            output += "#     total_cross_section = {} nb\n".format(self._total_cross_section_nb)
            output += "#     effective_lumi = {} invnb\n".format(self.lumi_invnb())
            output += "#     estimated_ratio_total_ev_to_recorded = {}\n".format(self._iev_ratio_truth_to_recorded)
        except AttributeError:
            pass
        
        return output
    
#----------------------------------------------------------------------        
class EventNtuple(object):
    def __init__(self, items):
        for item in items:
            parts = item.split('#')
            newname = parts[0].replace('.','__')
            #newname = parts[0].replace(':','__')
            self.__setattr__(newname,float(parts[1]))

    def has(self,item):
        return item in self.__dict__

    def __getitem__(self,item):
        return self.__dict__[item]


#----------------------------------------------------------------------
def return_selector(condition):
    """Return a function that selects based on the condition"""
    exec("def _select(ev): return ("+condition+")")
    return _select

    

if __name__ == '__main__': main()
