#!/usr/bin/perl -w
#
# In master version of this code, see test sets in test-sets/combruns-in0x.dat
# 
# NB 2008-12-22: sort out proper weighted recombination in cases where
#                one weight is zero
#
# script to merge files from many runs into a single one, in proportion 
# to the number of events seen in each file
#
# Usage:
#  combine-runs.pl [-err [-noNaN] ] [-ewgt errcol] [-ecol errcol1[,errcol2[...]]\
#                  [-sumfromcol n] [-ignorecomments] [-o outfile] [-nev]\
#                  [-comment char]
#                  file1 file2 [...]
#
# the number of events in each file is determined from the first "nev" tag
# (with the exception of "-nev" which can be the number of events 
# requested rather than the number of events found).
#
# following that, all numbers are added in proportion to the number
# of events in the file (including numbers in comments); note that
# errors are not treated properly...
#
# if a combination includes one or more non-numerical results
# (e.g. NaN, comment), then the result contains the first
# non-numerical string.
#
#  OPTIONS
#  -------
#
#  -err     if present, then each number in the set of input files is
#           replaced with the number and its error as determined from
#           the standard deviation between all files / sqrt(nfiles-1)
#
#  -sumfromcol N 
#           option sums rathers than averages things, starting
#           from column N (files start with col 1). Note: this has not
#           been checked with error calcs (not even clear what should
#           be done...)
#
#  -ewgt errcol
#           column errcol is assumed to contain error and the weight
#           of a given bin of a given file will be taken 
#           proportional to the squared inverse of the error
#
#  -ecol errcol1[,errcol2[,...]]
#           in the combination treat the given columns as if they
#           contain errors; if this is absent then strings in the file
#           "# ecol = i[,j[,...]]" are interpreted to indicate which
#           columns should be combined as errors. 
#           [reset by a blank line -- may be an issue for 2D histograms]
#
#  -comment char
#           set the comment character (default = '#')
#
#  -ignorecomments
#           don't attempt to combine numbers inside comments (lines
#           starting with #; doesn't yet work with non-default comment characters)
#
#  -noNaN
#           This option suppresses writing of NaN as errors for
#           non-numerical words.  warning: this option could mess up
#           the columns if there is a real NaN somewhere instead of a
#           number in the other files
#
#  -nevinfo just return information on the number of events
#
#  -o file  the output is to be sent to this file, rather than standard output
#       

use File::Glob ':glob';
use File::Glob ':globally';
#----------------------------------------------------------------------
# Version history
$version = "2.3.0 (2024-01-30)";
#
# 2.3.0 (2024-01-30)
# - resolved long-standing issue that quantities with +- and ± did not 
#   have any error treatment applied in the combination
#
# 2.2.1 (2020-02-25)
# - fixed issue with sqrt of negative number in error calc
#
# 2.2.0 (2020-02-09)
# - following more discussions with Gregory, reworked the error 
#   calculation to reflect what is used by 
#   https://docs.scipy.org/doc/numpy/reference/generated/numpy.cov.html
#   and gsl_stats_wvariance 
#
# 2.1.1 (2020-01-29)
# - comments out some unused $outerr calculations
#
# 2.0.0 (2020-01-29)
# - fixes issue of wrong error calculation with "-err" option when
#   combining runs with disparate numbers of events
#
#----------------------------------------------------------------------

# for finding numbers
#$numericregex = qw/[+-]?[0-9.]+([EDed][+-]?[0-9]+)?/;
#$numericregex = qw/[+-]?(([0-9]*\.?[0-9]+)|([0-9]+\.?))([EDed][+-]?[0-9]+)?/;
$numericregex = qw/[+-]?((\d*\.?\d+)|(\d+\.?))([EDed][+-]?\d+)?/;

# figure out options
$commandline = "$0 '".join("' '",@ARGV)."'";
$outfile="";
$err="";
$noNaN="";
$wgtcol=-1;
# indicates whether error columns are fixed from the
# command-line (or instead deduced from "ecol=" entries in the file)
$errcolFromCmdLine=0;
@errcol = ();
$sumfromcol=999999999;
$ignorecomments=0;
$nevinfo="";
$commentchar="#";
while ($ARGV[0] =~ /^-/) {
  if ($ARGV[0] eq "-err") {
    $err=1;
    shift @ARGV;
  } elsif ($ARGV[0] eq "-o") {
    shift @ARGV;
    $outfile = shift @ARGV;
  } elsif ($ARGV[0] eq "-nevinfo") {
    shift @ARGV;
    $nevinfo = "yes";
  } elsif ($ARGV[0] eq "-noNaN") {
    # This option suppresses writing of NaN as errors for
    # non-numerical words.
    #
    # warning: this option could mess up the columns if there is a real
    # NaN somewhere instead of a number in the other files
    $noNaN = 1;
    shift @ARGV;
  } elsif ($ARGV[0] eq "-ewgt") {
    # "-ewgt n" indicates that, where available, column n should be taken
    #           as an error and used to determine the weights for the other 
    #           entries in this row; the result of that column will then
    #           by the combined error
    shift @ARGV;
    $wgtcol = shift @ARGV;
    if (! is_numeric($wgtcol)) {
      die "Could not make sense of argument to '-w' option, $wgtcol";}
    # assume user numbers columns from 1; for us it's easier to start from zero
    $wgtcol-- ;
    # register the wgtcol as an error column (i.e. needing special recombination)
    $errcol[$wgtcol] = 1;
  } elsif ($ARGV[0] eq "-ecol") {
    shift @ARGV;
    foreach $ecol (split(",", shift@ARGV)) {
      $errcol[$ecol-1] = 1;
    }
    $errcolFromCmdLine = 1;
  } elsif ($ARGV[0] eq "-sumfromcol") {
    # sumfromcol: 1 is first column
    shift @ARGV;
    $sumfromcol = shift @ARGV;
  } elsif ($ARGV[0] eq "-comment") {
    shift @ARGV;
    $commentchar = shift @ARGV;
  } elsif ($ARGV[0] eq "-ignorecomments") {
    $ignorecomments=1;
    shift @ARGV;
  } else {
    die "Unrecognized option: ".$ARGV[0];
  }
}


@files=();
foreach $argv (@ARGV) {
  push @files, bsd_glob("$argv");
}
$nfile = $#files + 1;


@nev     = (); $nev=0;
@HANDLES = ();


$OUT = STDOUT;
if ($outfile) {open ($OUT, ">$outfile") || die "Could not open $outfile for output"}
if ($nevinfo) {open ($OUT, ">/dev/null") || die "Could not open /dev/null for output"}

print $OUT $commentchar." $commandline\n";
print $OUT $commentchar." until occurrence of 'nev', comments that follow are from first of the ",$nfile," files, $files[0]\n";

for($ifile = 0; $ifile < $nfile; $ifile++) {
  # support for compressed files
  $thisfile = $files[$ifile];
  if    ($thisfile =~ /\.bz2$/) {$open = "bzcat -c $thisfile|";}
  elsif ($thisfile =~ /\.gz$/)  {$open = "gunzip -c $thisfile|";}
  else                          {$open = "< $thisfile";}
  open($HANDLES[$ifile], "$open") || die "Could not read '$open'";
  $HANDLE = $HANDLES[$ifile];
  # figure out the number of events
  while ($line = <$HANDLE>) {
    if ($line =~ /[^a-z_-]nev[^a-z_-]+?($numericregex)/i && $line !~ /\s-nev/i) {
      $nev[$ifile]  = $1; 
      $nev         += $1;
      #print STDERR "file $ifile had $nev[$ifile] events\n";
      #print STDERR $line;
      last;}
    elsif ($ifile == 0) {print $OUT $line;}
  }
}

print $OUT $commentchar." combined nev = ",1.0*$nev,", from ",$#files+1," files (combine-runs.pl version $version)\n";
print STDERR "Total number of events is $nev",sprintf("(%.3e)",$nev*1.0),", from ",$#files+1," files";
if ($outfile) {print STDERR ", being written to $outfile";}
print STDERR "\n";
if ($nevinfo) {exit(0);}

$nLines = 0;
$n_identical_lines = 0;
while () {
  $good = 0;
  @out = ();
  @outerr = ();
  @sum_wgt_yi2     = ();
  @sum_wgt2     = ();
  #@sum_wgt2_yi  = ();
  #@sum_wgt2_yi2 = ();
  @notnumeric=();
  @wgt=();
  $totwgt=0;
  $wascomment = 0;
  $identicalLine = 0;
  for($ifile = 0; $ifile < $nfile; $ifile++) {
    $HANDLE = $HANDLES[$ifile];
    #while ($good = ($line = <$HANDLE>)) {if ($line !~ /^# /) {last;}}
    $good = ($line = <$HANDLE>);
    if (!$good) {
      if ($ifile != 0) {print STDERR "ERROR: file ",$files[$ifile]," was shorter than earlier files\n";}
      last;
    }

    # allow for the file to tell us which columns are error columns from here onwards
    if ($ifile == 0 && $line =~ /^# *ecol *=(.*)/) {
      if ($errcolFromCmdLine) {
        print STDERR "Error: ecol was specified from command line, but is also in file, cf line\n$line\n";
        exit(-1);
      }
      $ecol_list = $1;
      @errcol = ();
      foreach $ecol (split(",", $ecol_list)) {
        $errcol[$ecol-1] = 1;
      }
    }

    $comment_line = ($line =~ /^$commentchar/);
    if ($ignorecomments && $comment_line) {
      # take comment from first instance
      if ($ifile == 0) {
        print $OUT $line;
        $wascomment=1;
      }
      next;
    }
    chomp ($line);

    # a check for identical lines, e.g. in case seeds were accidentally set the same
    if ($ifile == 0) {
      $identicalLine = 0;
      if (!$comment_line) {$nLines++;}
    } else {
      if (!$comment_line && $line eq $lastline) {$identicalLine = 1}
    }
    $lastline = $line;


    # if a line is empty and we are not taking ecol info from the
    # command line, then reset the ecol information, so that it doesn't
    # propagate between histograms that have different structures.
    # (implicit assumption is that blank lines separate histograms)
    if (!$errcolFromCmdLine && $line eq "") {@errcol=();}

    @pieces = split(/[\s,]+/, $line);
    #@pieces = split(/\s+/, $line);
    # lines that start with a space end up having an empty first column.
    # This confuses column numbering for the user; so fix it here.
    if ($#pieces > 0 && $pieces[0] eq "") {shift @pieces;}
    #print STDERR $pieces[0],"\n";

    # support for error-based weighting
    if ($wgtcol >= 0 && $#pieces >= $wgtcol && is_numeric($pieces[$wgtcol])) { 
      $wgt[$ifile] = $pieces[$wgtcol] > 0 ? 1.0/($pieces[$wgtcol]**2) : 0.0;
    } else {
      $wgt[$ifile] = $nev[$ifile];
    }
    $totwgt += $wgt[$ifile];

    # get the information from this row for this file
    # (setting up variables needed to handle temporary error columns)
    $use_tmp_errcol = 0;
    @errcol_save = @errcol;
    for ($icol = 0; $icol <= $#pieces; $icol++) {
      if ($notnumeric[$icol]) {next;}
      $numeric = is_numeric($pieces[$icol]);

      if (!$numeric) {
        $notnumeric[$icol] = 1;
        $out[$icol] = $pieces[$icol];

        # now check for +- or ±, and if so set a tmp error column
        if ($pieces[$icol] eq "+-" || $pieces[$icol] eq "±") {
          if (!use_tmp_errcol) {
            $use_tmp_errcol = 1;
            @errcol = ();
          }
          $errcol[$icol+1] = 1;
        }
      } else {
        $thiswgt = ($icol+1 >= $sumfromcol) ? 1 : $wgt[$ifile];
        $piece = $thiswgt * $pieces[$icol];
        # sum squares if it's an error column
        if ($errcol[$icol]) {$piece = $piece**2;}
        if ($ifile == 0) {
          $out[$icol] = $piece;
          #$outerr[$icol] = $thiswgt * ($pieces[$icol])**2;
          $sum_wgt_yi2[$icol]  = $thiswgt * ($pieces[$icol])**2;
          $sum_wgt2[$icol]     = $thiswgt**2;
          # $sum_wgt2_yi[$icol]  = $thiswgt**2 * $pieces[$icol];
          # $sum_wgt2_yi2[$icol] = $thiswgt**2 * ($pieces[$icol])**2;
        } else {
          $out[$icol] += $piece;
          #$outerr[$icol] += $thiswgt * ($pieces[$icol])**2;
          $sum_wgt_yi2[$icol]  += $thiswgt * ($pieces[$icol])**2;
          $sum_wgt2[$icol]     += $thiswgt**2;
          # $sum_wgt2_yi[$icol]  += $thiswgt**2 * $pieces[$icol];
          # $sum_wgt2_yi2[$icol] += $thiswgt**2 * ($pieces[$icol])**2;
        }
      }
    }
    # reset the error column information if we were using a tmp error column to handle things like +- or ±
    if ($use_tmp_errcol) {
      @errcol = @errcol_save;
    }
  }
  if (!$good) {last;}
  if ($wascomment) {next;}
  # avoid divide by zero
  if ($totwgt == 0) {$totwgt = 1;}
  for ($icol = 0; $icol <= $#out; $icol++) {
    if (!$notnumeric[$icol]) {
      if ($icol+1 < $sumfromcol) {
	      $renorm = $totwgt;
      } else {
      	$renorm = 1;
      }
      $out[$icol] /= $renorm;
      # for an error column, get it to be err = sqrt(sum_i w_i^2 err_i^2)
      if ($errcol[$icol]) {$out[$icol] = sqrt($out[$icol]/$renorm);}
      if ($err) {
        # See also https://en.wikipedia.org/wiki/Sample_mean_and_covariance
        #
        # suppose we have files i=0..N-1, each with nev_i events, summing to nev_tot events;
        # The error on the sum can be worked out as follows
        #
        # - estimate of variance of event set i is var_i = (y_i - y_mean)^2 * nev_i
        # - each event-variance estimate is equally good so we just average them 
        # - estimate of variance as a whole is <var> = sum_i var_i / nfiles
        # - final error is err = sqrt(<var>/nev_tot)
        # - overall that simplifies to err = sqrt[sum_i w_i (y_i - y_mean)^2 / nfiles]
        #
        # Now we need to simplify that:
        # -   sum_i w_i (y_i - y_mean)^2 
        #     = sum_i w_i y_i^2 - 2*sum_i w_i y_i*y_mean + sum_i w_i *y_mean^2
        #     = sum_i w_i y_i^2 - sum_i w_i *y_mean^2
        # 
        # We also include a factor 1/sqrt(1 - sum w_i^2), cf. the
        # documentation of gsl_stats_wvariance (but I haven't derived
        # it, so am not entirely comfortable with it)        
        $mean = $out[$icol];
        $sum_wgt_yi2[$icol] /= $renorm;
        $sum_wgt2[$icol] /= $renorm**2;
        $outerr[$icol] = sqrt(abs(($sum_wgt_yi2[$icol] - $mean**2)/$nfile/(1-$sum_wgt2[$icol])));
        #$outerr[$icol] = sqrt(abs($sum_wgt2_yi2[$icol] -
        #                         2*$sum_wgt2_yi[$icol]*$mean + $sum_wgt2[$icol]*$mean**2)) / $renorm;
        # the following factor adjusts for small numbers of files and assumes file
        # weights are all equal
        #$outerr[$icol] *= sqrt($nfile / ($nfile-1));
        $out[$icol] .= " ".$outerr[$icol];
      }
    } elsif ($err && !$noNaN) {$out[$icol] .= " NaN";}
  };
  if ($identicalLine) {$n_identical_lines++;}
  print $OUT join(' ', @out)."\n";
}
if ($n_identical_lines > 0) {
  print STDERR "WARNING: found $n_identical_lines instances of identical non-comment lines between files, ".
               "out of a total of $nLines non-comment lines. NB some may be legitimate, e.g. histogram bins with zero events.\n";
}


# return true if it's a number
sub is_numeric {
  (my $arg) = @_;
  return ($arg =~ /^$numericregex$/i);
}


# # taken from "perldoc -q determine"
# # but seems marginally slower
# sub getnum {
#   use POSIX qw(strtod);
#   my $str = shift;
#   #$str =~ s/^\s+//;
#   #$str =~ s/\s+$//;
#   $! = 0;
#   my($num, $unparsed) = strtod($str);
#   if (($str eq '') || ($unparsed != 0) || $!) {
#     return undef;
#   }
#   else {
#     return $num;
#   }
# }
# 
# # return true 
# sub is_numeric { defined getnum($_[0]) }


## # open the files, start off with definite number of events
## $nevtot = $nfile;
## for ($ifile=0; $ifile < $nfile; $ifile++) {
##   $nev[$ifile]=1;
##   open ($file[$ifile], "< $ARGV[$ifile]") || die "Could not open $ARGV[$ifile]";
## }
## 
## 
## # run through the files, line-by-line
