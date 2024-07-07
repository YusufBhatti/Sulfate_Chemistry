#!/usr/bin/env perl
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved. 
# For further details please refer to the file COPYRIGHT.txt 
# which you should have received as part of this distribution. 
# *****************************COPYRIGHT******************************* 

# Script to check to see if any new ifdefs have been added in an fcm diff

use strict;
use warnings;
use Data::Dumper;
use 5.010;

# Read branch from command line or assume current directory
my $branch = shift // '.';

# Read text file of retired if-defs from command line or use default file
my $retired_ifdef_file = shift;

unless ($retired_ifdef_file and -f $retired_ifdef_file) {
  die "Retired ifdef filename not provided.\n"
}

# Read in retired if-defs
my @retired = read_file($retired_ifdef_file);

# Exit code - value will be read by the UTF
my $exit = 0;

# File to source to access FCM commands
my $fcm = '/etc/profile';

my $snooze = 120;         # Time to wait before retrying
my $max_snooze = 10;      # Maximum number of retries before aborting
my $suite_mode = 0;

if ($ENV{SOURCE_UM_MIRROR}) {
  print "Detected SOURCE_UM_MIRROR environment variable.\n";
  $branch = $ENV{SOURCE_UM_MIRROR};
  print "Redirecting branch to $branch\n";
  $suite_mode = 1;
}

# Check this is a branch rather than the trunk
my @binfo = `. $fcm; fcm binfo $branch 2>&1`;
my $binfocode = $?;

if (/URL: svn:\/\/fcm\d+\/\w+_svn\/\w+\/trunk/ ~~ @binfo   or
    /URL:.*\/svn\/\w+\/main\/trunk/ ~~ @binfo or
    /URL:..*_svn\/main\/trunk/ ~~ @binfo) {
  print "Detected trunk, nothing to diff against!\n";
  exit 0;
}

unless ($binfocode == 0) {
  if (/svn info --xml/ ~~ @binfo) {
    if ($suite_mode) {
      for (my $i = 1; $i <= $max_snooze; $i++) {
        print "Revision probably doesn't exist yet - waiting $snooze seconds for mirror to update (Snooze $i of $max_snooze).\n";
        sleep $snooze;
        @binfo = `. $fcm; fcm binfo $branch 2>&1`;
        $binfocode = $?;
        last if ($binfocode == 0);
      }
    }
  }
  if ($binfocode != 0) {
    print "Error running fcm binfo:\n";
    print @binfo;
    exit 1;
  }
}

# Results of an 'fcm bdiff' command
my @diff = `. $fcm; fcm bdiff $branch 2>&1`;
my $diffcode = $?;

# Check the bdiff worked correctly
unless ($diffcode == 0 ) {
  die "Error running 'fcm bdiff $branch':\n@diff\n";
}

check_for_wrapping_ifdef(@diff);

# Get lines added
my %added_code = get_added_code(@diff);

# parse the changes for each file separately 
foreach my $filename (keys %added_code) {

  # find the added lines for this file
  my @file_lines = @{$added_code{$filename}};

  # merge cpp continuation lines, removing comments
  my @cpp_lines = get_cpp_lines(remove_comments(merge_cpp_continuation_lines(@file_lines))); 

  # Parse lines for all ifdefs
  my %ifdefs = find_ifdef(@cpp_lines);

  # Compare found ifdefs with those which have been retired
  test_for_retired(\%ifdefs, \@retired, $filename);
}

if ($exit == 0) {
  print "No retired if-defs appear to be used by this branch.\n";
}

exit $exit;

sub test_for_retired {
  my $ifdef_ref = shift;
  my $retired_ref = shift;
  my $filename = shift;
  
  my %ifdefs = %$ifdef_ref;
  my @retired = @$retired_ref;
  
  foreach my $def (keys %ifdefs) {
    if ($def ~~ @retired) {
      print "$def has been retired but is used by $filename !\n";
      $exit = 1;
    }
  }
} 

# Returns added code given an fcm diff, based on the tab checking script
# /home/h01/frum/Unsupported/fortran_tab.pl by Tom Green (frtg)
sub get_added_code {
  my @input = @_;
  my %output = ();

  my $code_file = 0;
  my $filename = "";

  foreach my $line (@input) {
    if ($line =~ /^\+\+\+/) {
      if ($line =~ m/\.[chfF](9[05])?\s/) {
        $code_file = 1;
        my @filenames = split /\s+/, $line;
        $filename = $filenames[1];
      } else {
        $code_file = 0;
      }
    }
    if ($line =~ /^\+[^+]/ and $code_file == 1) {
      my $file_url;
      my $url_revision;
      my $cat_path;

      if ($branch =~/:/) {
        $cat_path = $branch;
        if ($branch =~/@/) {
          $cat_path =~s/(@.*)//;
          $url_revision = $1;
        }
        $cat_path = "$cat_path/";     
      } else {
        $cat_path = '';
      }

      # The $url_revision variable is only present if the URL is a branch
      if ($url_revision) {
        $file_url = "$cat_path$filename$url_revision";
      } else {
        $file_url = "$cat_path$filename";
      }

      $code_file = 0;
      push @{$output{$filename}},cat_file($file_url);
    }
  
  }

  # clean up files with no additions
  foreach my $key (keys %output) {
     delete $output{$key} if !defined($output{$key});
  }

  return %output;
}

sub merge_cpp_continuation_lines {
  my @lines = @_;
  my @newlines;

  # Sort out C continuation lines
  my $entire = join("\n",@lines);
  $entire =~s/\\\s*\n//g;
  @newlines = split/\n+/,$entire;  

  return @newlines;
}


sub get_cpp_lines {
  my @lines = @_;
  my @newlines;
  
  foreach my $line (@lines) {
    if ($line =~/^\s*#/) {
      push @newlines, $line;
    }
  }
  
  return @newlines;
}

sub remove_comments {
  my @lines = @_;
  my @newlines;

  #remove commented sections
  my $entire = join("\n",@lines);
  $entire =~s/\/\*(.|\n)+?(\*\/)//g;
  @newlines = split/\n+/,$entire;  

  return @newlines;
}

sub read_file {
  my $file = shift;
  open (my $fh, '<', $file) or die "Cannot read $file: $!\n";
  chomp (my @lines = <$fh>);
  close $fh;
  return @lines;
}

sub find_ifdef {
  my @lines = @_;
  my @elements;

  my %ifdefs;

  foreach my $line (@lines) {
    # test for #ifdef style tests
    if ($line =~/^\s*#if(n)?def\s+(?<ifdef>\w+)/) {
      my $ifdef = $+{ifdef};
      $ifdefs{$ifdef}++;
    }

    # standardise format for defined(<DEF>) style tests
    $line =~s/defined\s*?\(?\s*?(\w+)[^\S\n]*\)?([|&><*+%^$()\/\-\s])/defined($1) $2/g;
              
    push @elements, split /\s+/, $line;
  }

  # test for defined(<DEF>) style tests
  foreach my $element (@elements) {
    if ($element =~/defined\((?<ifdef>\w+)\)/) {
      my $ifdef = $+{ifdef};
      $ifdefs{$ifdef}++;
    }
  }

  return %ifdefs;
}

sub check_for_wrapping_ifdef {
  my @lines = @_;
  my $fname;
  for (my $i = 0; $i < scalar @lines; $i++) {
    my $line = $lines[$i];
    if ($line =~/^\+\+\+/) {
      $line =~/^\+\+\+\s*(\S+)/;
      $fname = $1;
      next unless ($fname =~/F90$/);
      my $diffline = $lines[$i+1];
      my $firstline = $lines[$i+2];
      if ($diffline =~/\+1,/) {
        if ($firstline =~/^\+#if/) {
          print "File $fname appears to contain a wrapping if-def\n";
          $exit = 1;
        }
      }
    }
  
  }

  return;
}

# Cat a file, either from fcm (if the URL contains a colon) or from disk
sub cat_file {
  my $url = shift;
  my @lines;
  my $error=0;

# If the URL contains a colon treat it as an fcm, else treat as a regular file  
  if ($url =~/:/) {
    @lines = `. $fcm; fcm cat $url 2>&1`;
    $error=$?;
  } else {
    @lines = `cat $url 2>&1`;
    $error=$?;
  }

  if($error) {
    print "Error cating file $url\n";
    exit 1;
  } 

  return @lines;
}