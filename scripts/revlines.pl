#!/usr/bin/perl -w
# script to reverse the lines coming in on stdin

@lines=();
while ($line=<STDIN>) {
  push @lines , $line;
}


while ($line= pop @lines) {
  print $line;
}


