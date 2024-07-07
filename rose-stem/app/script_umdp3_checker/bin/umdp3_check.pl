#!/usr/bin/env perl
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved. 
# For further details please refer to the file COPYRIGHT.txt 
# which you should have received as part of this distribution. 
# *****************************COPYRIGHT******************************* 

# Script to check whether a code change complies with UMDP 003
# Basically the 'guts' of the UMDP3.pm class's perform_task method without 
# the UTF-isms such as status objects and the OO stuff

use strict;
use warnings;
use 5.010;
use Cwd 'abs_path';

# Set the location of the UMDP3 package.
use FindBin;
use lib "$FindBin::Bin";

use UMDP3;

# This is a standalone version of the dispatch tables from UMDP3Job, generated
# by script automatically
use UMDP3DispatchTables;


# Declare variables
my $fcm = '/etc/profile'; # File to source to access 'fcm' commands
my $exit = 0;             # Exit code of this script == number of failing files
my %additions;            # Hash of added code
my %deletions;            # Hash of deleted files
my $filename = '';        # Current filename being tested
my $snooze = 120;         # Time to wait before retrying
my $max_snooze = 10;      # Maximum number of retries before aborting

# Get argument from command line, or assume '.' if not defined
my $branch = shift // '.';

# Cope with UTF-style working copy syntax just in case
$branch =~s/wc://;

# Read text file of whitelisted include files
my $whitelist_includes_file = shift;

unless ($whitelist_includes_file and -f $whitelist_includes_file) {
  die "Whitelist filename not provided.\n"
}

# Read in retired if-defs
my @includes = read_file($whitelist_includes_file);

my $suite_mode = 0;
if ($ENV{SOURCE_UM_MIRROR}) {
  print "Detected SOURCE_UM_MIRROR environment variable.\n";
  $branch = $ENV{SOURCE_UM_MIRROR};
  print "Redirecting branch to $branch\n";
  $suite_mode = 1;
}

my %dispatch_table_diff_fortran = UMDP3DispatchTables::get_diff_dispatch_table_fortran();
my %dispatch_table_file_fortran = UMDP3DispatchTables::get_file_dispatch_table_fortran();
my %dispatch_table_diff_c = UMDP3DispatchTables::get_diff_dispatch_table_c();
my %dispatch_table_file_c = UMDP3DispatchTables::get_file_dispatch_table_c();

my @binfo;
my $binfocode;

start_branch_checking:

# Check this is a branch rather than the trunk
@binfo = `. $fcm; fcm binfo $branch 2>&1`;
$binfocode = $?;

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

foreach my $line (@binfo) {
  if ($line =~/Branch Parent:.*\/trunk@.*/) {
    last;
  } elsif ($line =~/Branch Parent:\s*(.*)/) {
    print "This branch is a branch-of-branch - testing parent ($1)\n";
    $branch = $1;
    goto start_branch_checking;
  }   
}  

my @info;

# Get fcm info for branch
@info = `. $fcm; fcm info $branch 2>&1`;
$binfocode = $?;

if ($binfocode != 0) {
    print "Error running fcm info:\n";
    print @info;
    exit 1;
}

my $repository_branch_path;
my $repository_working_path;
my $repository_relative_path;

foreach my $line (@binfo) {
  if ($line =~/^URL:\s*(.*)/) {
    $repository_branch_path = $1;
    last;
  }
}  

foreach my $line (@info) {
  if ($line =~/^URL:\s*(.*)/) {
    $repository_working_path = $1;
    last;
  }
}  

$repository_relative_path = $repository_working_path;
$repository_relative_path =~ s/$repository_branch_path//;

# replace relative branch paths with absolute paths
if(grep(/Working Copy Root Path:/, @info)) {
  $branch = abs_path($branch);
}

# trim trailing "/"
$branch =~  s/\/$//;

print "Testing branch $branch\n";

if ($repository_relative_path) {
  print "\n[WARN] The relative path between the root of the branch and the script working path ($repository_relative_path) is not empty\n";
  print "       - you are not running from the root of the branch\n\n";
  if ($suite_mode) {
    die "Error - re-run from the root of the branch\n";
  }
}

# Get the diff
my @diff = `. $fcm; fcm bdiff $branch  2>&1`;
my $diffcode = $?;

# Check the bdiff worked correctly
unless ($diffcode == 0 ) {
  die "Error running 'fcm bdiff $branch':\n@diff\n";
}

# We will need to know empty and deleted files - use the bdiff summary to identify these.
my @summary = `. $fcm; fcm bdiff --summarise $branch  2>&1`;
$diffcode = $?;

# Check the second bdiff worked correctly
unless ($diffcode == 0 ) {
  die "Error running 'fcm bdiff --summarise $branch':\n@summary\n";
}

foreach my $line (@summary) {

  # Reset captures to undefined with a trivial successful match.
  "a" =~ /a/;  

  # Add hash entries for added or modified files:
  # These are files which are newly added; or which add or remove lines.
  $line =~ /^(A|M)\s*(?<filename>\S+)$/;
  my $modified_file = $+{filename};
  if($modified_file) {
    #normalise the path
    $modified_file =~ s/$repository_working_path\///;
    $modified_file =~ s/.*trunk$repository_relative_path\///;

    $additions{$modified_file}=[];
  }

  # Reset captures to undefined with a trivial successful match.
  "a" =~ /a/;  
  
  # Add has entries for deleted files
  $line =~ /^D\s*(?<filename>\S+)$/;
  my $deleted_file = $+{filename};
  if($deleted_file) {
    #normalise the path
    $deleted_file =~ s/$repository_working_path\///;
    $deleted_file =~ s/.*trunk$repository_relative_path\///;
    $deletions{$deleted_file}=[];
  }
}

my $store_line = 0;

# Store the lines added in a hash with the filename as the key, i.e.
# %additions =  ( 'filename' => [ 'added line 1', 'added line 2'] )
foreach my $line (@diff) {

  if ($line =~/^\+\+\+/) {
    # Find if the filename is in our additions hash, 
    # and set the subsequent lines to be stored if it is.
    $line =~ /^\+\+\+\s+(?<filename>\S+)/;
    $filename = $+{filename};
    unless (($branch eq ".") || ($filename eq $branch)) {
      $filename =~ s/.*$branch\///;
    }
    $store_line = exists($additions{$filename});

    if($store_line == 0) {
      # if we don't recognise the file as deleted, 
      # or as marking an SVN property change on the branch root, 
      # something has gone wrong.
      if (!exists($deletions{$filename})) {
        my $bdiff_branch_root = $branch;
        if ($filename eq ".") {
          $bdiff_branch_root = ".";
        }
        if ( !(($filename eq $bdiff_branch_root) && grep(/^Property changes on: $bdiff_branch_root$/, @diff)) ) {
          print "Something has failed parsing line '$line'\n";
          print "Filename '$filename' is not contained in the output from fcm bdiff --summarise!\n";
          exit 1;
        }
      }
    }

  } elsif ($line =~/^\+/) {
    if($store_line) {
      # Add the diff to %additions hash
      $line =~ s/^\+//;
      push @{$additions{$filename}}, $line;
    }
  }
}

# This section prints failure messages for each file after each file is tested
print "The following files have failed the UMDP3 compliance tests:\n";

# Set up the error message string to empty
my $message = '';

# set up known includes whitelist
my %includes_hash;
@includes_hash{@includes}=();

# Loop over modified files
foreach my $modified_file (keys %additions) {
  # Initialise variables
  my $failed = 0;
  my @failed_tests;
  my $is_c_file = 0;
  my $is_fortran_include_file = 0;

  # If it's an include file, fail unless it's on the include whitelist 
  # (e.g. its a C header or a Fortran include for reducing code duplication).
  if ($modified_file =~/\.h$/) {

    if (exists($includes_hash{$modified_file})) {
      my @components = split("/", $modified_file);
      if ($components[0] =~ /src/ and $components[-2] =~/include/ and not $components[1] =~ /include/) { 
        $is_fortran_include_file = 1;
      } elsif ($components[0] =~ /src/ and $components[1] =~/include/) {
        $is_c_file = 1;
      } else {
        push @failed_tests, "Added an include file outside of a recognised 'include' directory";
      }
    } else {
      push @failed_tests, "Modified or created non-whitelisted include file rather than using a module";
      $failed++;
    }
  }

  if ($modified_file =~/\.c$/) {
      $is_c_file = 1;
  }

  # if it's Fortran or C apply all the tests
  if ($modified_file =~/\.F90$/ or $modified_file =~/\.f90$/ or $is_c_file or $is_fortran_include_file) {

    my $dispatch_table_diff;
    my $dispatch_table_file;

    if ($is_c_file) {
      $dispatch_table_diff = \%dispatch_table_diff_c;
      $dispatch_table_file = \%dispatch_table_file_c;
    } else {
      $dispatch_table_diff = \%dispatch_table_diff_fortran;
      $dispatch_table_file = \%dispatch_table_file_fortran;
    }

    # Get the diff for this file out of the hash
    my $added_lines_ref = $additions{$modified_file};
    my @added_lines = @$added_lines_ref;

    # Loop over each test which works on a diff
    foreach my $testname (keys %$dispatch_table_diff) {
      UMDP3::reset_extra_error_information();

      # Get the subroutine reference from the tables at the top of this file
      my $subroutine_ref = ${$dispatch_table_diff}{$testname};

      # Run the test
      my $answer = &$subroutine_ref(@added_lines);

      my %extra_error = UMDP3::get_extra_error_information();
      if (scalar keys %extra_error > 0) {
        my @extra_error = keys %extra_error;
        my $extra_text = join(", ",@extra_error);
        $testname .= ": $extra_text";
      }

      # If the test fails, increase the number of failures and add the testname
      # to the array containing the list of problems with this file
      if ($answer) {
        $failed++;
        push @failed_tests, $testname;
      }

    }

    # Analyse the command line argument to work out how to access the whole file
    my $file_url;
    my $url_revision;
    if ($branch =~/@/) {
      $branch =~s/(@.*)//;
      $url_revision = $1;
    }

    # The $url_revision variable is only present if the URL is a branch
    if ($url_revision) {
      $file_url = "$branch/$modified_file$url_revision";
    } else {
      $file_url = "$branch/$modified_file";
    }

    # Get the whole file contents
    my @file_lines = cat_file($file_url);

    # Perform each test which checks the whole file in a similar method to 
    # tests which work on a diff
    foreach my $testname (keys %$dispatch_table_file) {
      UMDP3::reset_extra_error_information();
      my $subroutine_ref = ${$dispatch_table_file}{$testname};
      my $answer = &$subroutine_ref(@file_lines);
      my %extra_error = UMDP3::get_extra_error_information();
      if (scalar keys %extra_error > 0) {
        my @extra_error = keys %extra_error;
        my $extra_text = join(", ",@extra_error);
        $testname .= ": $extra_text";
      }
      if ($answer) {
        $failed++;
        push @failed_tests, $testname;
      }
    }
  } # Filename matches F90/f90/c
  
  # If any tests failed, print the failure message
  if ($failed > 0) {
    my $failure_text = join("\n  ",@failed_tests);

    # The space before the colon makes the filename easier to cut and paste
    $message .= "File $modified_file :\n  $failure_text\n";
    print $message;
    $exit++;
  }
  $message = '';    
  
} # Loop over files

# Print message for a success if no files have failed
if ($exit == 0) {
  print "No modified files appear to have failed the compliance tests\n";
}

# Exit with number of fails, if it's zero (a UNIX success) it passed
exit $exit;


############################### SUBROUTINES ###################################

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

sub read_file {
  my $file = shift;
  open (my $fh, '<', $file) or die "Cannot read $file: $!\n";
  chomp (my @lines = <$fh>);
  close $fh;
  return @lines;
}

