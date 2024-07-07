# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************
#
# Script:     macros.pm
# Purpose:    Rose macros for rosebud-generated suites
# Code owner: UM System Development team

package macros;

use strict;
use warnings;
use fileio;


sub fixes {

# Apply the fixes macro to an app
# Arguments:
  my $app       = shift;     # Path to app directory
  my $macrofile = shift;     # File for macro output
  my $meta      = shift;     # Path to user-provided meta-data

  my $header = "=" x 80 . "\n Applying fixer macro to $app:\n";

  fileio::append_file($macrofile, $header);

  system "rose macro --non-interactive -F $meta -C $app >> $macrofile";

  return;
}

sub stash {

# Apply the STASH macro to the runtime app
# Arguments:
  my $app       = shift;     # Path to app directory
  my $meta      = shift;     # Path to user-provided meta-data

  system "rose macro --non-interactive -q $meta -C $app stash_indices.TidyStashTransformPruneDuplicated";

  return;
}

sub profs {

# Apply the STASH profile macro to the runtime app
# Arguments:
  my $app       = shift;     # Path to app directory
  my $meta      = shift;     # Path to user-provided meta-data

  system "rose macro --non-interactive -q $meta -C $app stash_requests.StashProfilesRemoveUnused";

  return;
}

sub upgrade {

# Upgrade an app to a specified version.
# If upgrading across multiple releases, do the upgrade in steps so we
# can run the fixer macro after each individual upgrade.
# Arguments:
  my $app        = shift;     # Path to app directory
  my $macrofile  = shift;     # File for macro output
  my $meta       = shift;     # Path to user-provided meta-data
  my $upgrade    = shift;     # Version to upgrade to

  # Extract the number of the target upgrade:
  my $version = substr($upgrade, 2);

  # Get the list of available upgrades:
  my @releases = `rose app-upgrade $meta -C $app`;
  my @available;    # List of available upgrades (vn9.0, vn9.1, ...)

  foreach my $line (@releases) {
    # Store each available upgrade version:
    if ($line =~ /[^=]\s+(vn\d+\.\d+)/) {
      push @available, $1;
    }
  }

  # Is our upgrade target in the list of available upgrades?
  my $can_upgrade = 0;
  foreach my $version (@available) {
    $can_upgrade = 1 if ($upgrade =~ /$version/);
  }

  # If we can't upgrade, tell the user and go no further.
  unless ($can_upgrade) {
    my $appname = $app;
    $appname =~ s!.*/(\w+)!$1!;  # Relying on greedy behaviour to get basename.
    print "[WARN] Cannot upgrade $appname app to $upgrade - invalid version.\n";
    print "[WARN] Upgrade of $appname app disabled.\n";
    return;
  }

  my $current = 'vn8.6';   # Current app version. Always true.
  foreach my $nextversion (@available) {

    last if ($current eq $upgrade);   # Stop if we've reached the target vn.
    my $header = "=" x 80
      . "\n Applying $current-$nextversion upgrade macro to $app:\n";

    fileio::append_file($macrofile, $header);

    system "rose app-upgrade --non-interactive $meta -C $app $nextversion >>$macrofile 2>&1";

    # Run the fixes macro to ensure the upgrade is valid and re-dump the file:
    fixes($app, $macrofile, $meta);

    $current = $nextversion;
  }

  return;
}

1;
