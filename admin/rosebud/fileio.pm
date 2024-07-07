# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************
#
# Script:     fileio.pm
# Purpose:    Rosebud file operations
# Code owner: UM System Development team

package fileio;

use strict;
use warnings;


sub read_file {

# Read a file

  my $file = shift;  # input filename
  open (my $fh, '<', $file) or die "[FAIL] Cannot read $file: $!\n";
  my @lines = <$fh>;
  close $fh;
  return @lines;
}

sub write_file {

# Write a file

  my $filename = shift; # Output filename
  my @lines = @_;       # Data for writing
  open (my $fh, '>', $filename) or die "[FAIL] Cannot write $filename: $!\n";
  print $fh @lines;
  close $fh;
  return;
}

sub append_file {

# Append to a file

  my $filename = shift; # Output filename
  my @lines = @_;       # Data for writing
  open (my $fh, '>>', $filename) or die "[FAIL] Cannot append to $filename: $!\n";
  print $fh @lines;
  close $fh;
  return;
}

sub create_dirs {

# Create the Rose suite directory structure
  my $hashref = shift;

# Transfer relevant hash values to variables for ease of use
  my $suitedir     = $hashref->{suitedir};     # Suite output directory
  my $makeapp      = $hashref->{makeapp};      # Name of the fcm_make[_um] app
  my $oceanmakeapp = $hashref->{oceanmakeapp}; # Name of the fcm_make_ocean app
  my $runapp       = $hashref->{runapp};       # Name of the atmos/recon app
  my $coupled      = $hashref->{coupled};      # True = Coupled model


  my @suitedirs = ("", "app", "app/$makeapp", "app/$makeapp/file", 
                   "app/$runapp", "app/$runapp/file", "meta");
  push @suitedirs, "app/$oceanmakeapp", "app/$oceanmakeapp/file" if ($coupled);

  foreach my $subdir (@suitedirs) {
    mkdir "$suitedir/$subdir", 0755 unless -z "$suitedir/$subdir";
  }
  return;
}

1;
