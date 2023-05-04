#!/usr/bin/env perl
use warnings;

# usage:
#
#    auto-combine-runs.pl [standard combine-runs.pl options]
#
# it works out which files to combine and creates an output name for
# them, assuming that anything that needs combining has the rseqNNNN
# format.

$do_covariance=0;
$xtraargs = "";
if ($#ARGV >= 0 && $ARGV[0] eq "-cov") {
    # remove the -cov argument
    $xtraargs = shift @ARGV;
    $do_covariance = 1;
}


$args = join(" ",@ARGV);


%fileGroup = ();
$directory=".";
opendir (DIR, $directory) or die $!;
while (my $file = readdir(DIR)) {

  if ($file =~ /(^.*)(rseq[0-9]+)([^~]*)$/) {
      if ($do_covariance) {$outname = "$1allcovrseq$3"; }
      else                {$outname = "$1allrseq$3";}
    #print "$file $1 $2 $3 $outname\n";
    if (!exists $fileGroup{$outname}) {$fileGroup{$outname} = [];}
    push @{$fileGroup{$outname}}, $file;
  }
}

$makefile = "Makefile";
# try to figure out if we have a makefile that
# already exists and with which we may conflict
# If so then we use a different name for our Makefile
if (-e $makefile) {
  open (MAKE, "<", $makefile);
  $firstline = <MAKE>;
  if ($firstline !~ /auto-combine-runs/) {$makefile = "Makefile.autocomb";}
  close MAKE;
}
print "Creating $makefile\n";

$combFiles=join(" ", sort(keys %fileGroup));
# write common structure of Makefile
open (MAKE, ">", $makefile);
print MAKE "# File automatically generated by $0 $args\n";
print MAKE "
all: summary $combFiles

make:
\t$0 $args $xtraargs

summary: $combFiles
\t".'@'."grep combined.nev $combFiles
";

# then write the structure that is specific to each file that
# needs generating.
# (Would it be better to use wildcards here so that if filelist is
# updated we don't need to regenerate the Makefile?)
for $outname (keys %fileGroup) {
    $fileList = join(" ", @{$fileGroup{$outname}});
    if ($do_covariance) {
	print MAKE "
$outname: $fileList
\t../../../../scripts/combine-with-cov.py $args $fileList > $outname 
";
    } else {
	print MAKE "
$outname: $fileList
\tcombine-runs.pl $args -o $outname $fileList

";
    }
}

system("make -f $makefile");
