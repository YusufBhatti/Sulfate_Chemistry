# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************
#
# Script:     make_meta.pm
# Purpose:    Write the rose-meta.conf file
# Code owner: UM System Development team

package make_meta;

use strict;
use warnings;
use fileio;


sub write_meta {

# Generate the rose-meta.conf file

  my $hashref = shift;

# Transfer relevant hash values to variables for ease of use
  my $suitedir = $hashref->{suitedir};  # Suite output directory

  my @meta; # File contents for writing

  push @meta, "[jinja2:suite.rc:UM hosts]\n";
  push @meta, "ns=UM hosts\n";
  push @meta, "\n";
  push @meta, "[jinja2:suite.rc=COMPUTE_HOST]\n";  # Formerly HPC_HOST
  push @meta, "compulsory=true\n";
  push @meta, "description=Host to use for compilation and/or model run\n";
  push @meta, "help=The remote host on which to run the model.\n",
              "    =\n",
              "    =On platforms which require an fcm_make2 task this is\n",
              "    =also the host which will be used to compile executables.\n";
  push @meta, "ns=UM hosts\n";
  push @meta, "sort-key=b\n";
  push @meta, "\n";
  push @meta, "[jinja2:suite.rc=EXTRACT_HOST]\n";  # Formerly DESKTOP_HOST
  push @meta, "compulsory=true\n";
  push @meta, "description=Host to use for code extraction\n";
  push @meta, "help=The local host on which to extract the model's source code.\n",
              "    =\n",
              "    =On platforms which don't require an fcm_make2 task\n",
              "    =this is also the host which will be used to compile executables.\n";
  push @meta, "ns=UM hosts\n";
  push @meta, "sort-key=a\n";

  fileio::write_file("$suitedir/meta/rose-meta.conf",@meta);
}

1;
