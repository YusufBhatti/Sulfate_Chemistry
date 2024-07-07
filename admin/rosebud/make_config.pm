# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************
#
# Script:     make_config.pm
# Purpose:    Write the fcm make config file and associated app for the UM
# Code owner: UM System Development team

package make_config;

use strict;
use warnings;
use fileio;
use macros;

sub write_fcm_make {

# Write the fcm_make[_um].cfg config and its rose-app.conf file.

  my $hashref = shift;

# Transfer relevant hash values to variables for ease of use
  my $jobdir  = $hashref->{jobdir};          # Job input directory
  my $suitedir = $hashref->{suitedir};       # Suite output directory
  my $platform = $hashref->{platform};       # IBM or Linux (or ...)?
  my $version = $hashref->{version};         # UM version x.y
  my $build_atmos = $hashref->{build_atmos}; # True = atmos build required
  my $build_recon = $hashref->{build_recon}; # True = recon build required
  my $makeapp = $hashref->{makeapp};         # Name of the fcm_make app
  my $unportable = $hashref->{unportable};   # No. of unportable settings
  my $coupled = $hashref->{coupled};         # True = Coupled model
  my $scm = $hashref->{scm};                 # True = SCM job
  my $fixes = $hashref->{fixes};             # True = apply fixes macro
  my $upgrade = $hashref->{upgrade};         # Version to upgrade to
  my $macrofile = $hashref->{macrofile};     # File containing macro output
  my $meta = $hashref->{meta};               # Path to meta-data
  my @extr_scr = @{$hashref->{extr_scr}};    # Contents of EXTR_SCR
  my $envhashref = $hashref->{envhashref};   # Reference to hash of 
                                             #  environment variables

# Local variables:
  my $openmp;           # OpenMP on or off
  my $openmpflag;       # Compiler flag for OpenMP
  my $model;            # String containing app-specific name (atmos, scm etc.)
  my $line;
  my $item;
  my $platform_config_dir;  # Platform-specific config directory
  my $model_config;         # Model/opt-specific config name 
                            #   (atmos-high, scm-debug, etc.)
  my $steplist = "";
  my @extract_location;
  my @make = "include = \$include_config\n\n"; # fcm_make.cfg contents
  my @app;                                     # rose-app.conf contents

# Some defaults for (potential) printing to file:
  my $um_base = "fcm:um_tr";
  my $um_rev  = "vn$version";
  my $jules_base = "fcm:jules_tr";
  my $jules_rev  = "um$version";
  my $opt = "high";

# Set up job-dependent strings:
  if ($scm) {
    $model = "scm";
  } else {
    $model = "atmos";
  }
  if ($platform eq 'ibm') {
    $openmpflag = '-qsmp=omp';
  } elsif ($platform eq 'linux') {
    $openmpflag = '-openmp';
  }

# Set up the app file header.
# Meta-data is hard-wired to vn8.6 (the earliest available).
  push @app, "meta=um-fcm-make/vn8.6\n\n[env]\n";


# 1. Create list of CPP keys for UM 
#  and pick up compiler optimisation level while we're here
  if ($version eq '8.3') {
    push @app, 'keys_model_app=';
  } else {
    push @app, "keys_${model}_app=";
  }

# Now create list of atmos CPP flags:

  my @atmoscfg = fileio::read_file("$jobdir/FCM_UMATMOS_CFG");

  foreach $line (@atmoscfg) {

    $line =~ s/#.*//;               # Discard comments

    if ($line =~ /(\w+)=\1/i) {   # Pick out jobdefs (FLAG=flag)

      $line =~ s/\s*(\w+)=(\w+).*/ $1=$2/;    # Strip whitespace
      chomp $line;
      push @app, $line;

      $openmp = "true" if ($line =~ /C98_1A/); # Check for compile-time OpenMP
    }
    $opt = $1 if ($line =~ /inc \$UM_SVN_BIND\/bind\d+_mpp_(\w+).cfg/);
  }
  push @app, "\n";

# 2. Create list of CPP keys for reconfiguration
  unless ($scm) {
    push @app, "keys_recon_app=";

    my @reconcfg = fileio::read_file("$jobdir/FCM_UMRECON_CFG");

    foreach $line (@reconcfg) {

      $line =~ s/#.*//;               # Discard comments

      if ($line =~ /(\w+)=\1/i) {     # Pick out jobdefs (FLAG=flag)

        $line =~ s/\s*(\w+)=(\w+).*/ $1=$2/;   # Strip whitespace
        chomp $line;
        push @app, $line;

        $openmp = "true" if ($line =~ /C98_1A/); # Check for compile-time OpenMP
      }
    }
    push @app, "\n";
  }


# 3. Set base & revision
  foreach $line (@extr_scr) {
    $um_base    = $1 if ($line =~ /export UM_SVN_URL=(.*)/);
    $um_rev     = $1 if ($line =~ /export UM_VN=(.*)/);
    $jules_base = $1 if ($line =~ /export JULES_SVN_URL=(.*)/);
    $jules_rev  = $1 if ($line =~ /export JULES_VER=(.*)/);
  }

# Don't need to include anything that matches the default settings
  push @app, "um_base=$um_base\n" unless
          ($um_base =~ m/fcm:um[_-]tr/ or
           $um_base =~ m!fcm:um/trunk! or
           $um_base =~ m!svn://fcm\d/UM_svn/UM/trunk!);
  $um_rev = process_base_rev($um_rev);
  push @app, "um_rev=$um_rev\n";

  push @app, "jules_base=$jules_base\n" unless
          ($jules_base =~ m/fcm:jules[_-]tr/ or
           $jules_base =~ m!fcm:jules/trunk! or
           $jules_base =~ m!svn://fcm\d/JULES_svn/JULES/trunk!);
  $jules_rev = process_base_rev($jules_rev);
  push @app, "jules_rev=$jules_rev\n";


# 4. Disable OpenMP in the platform config file, if necessary.
#    For the SCM, the reverse is true and OpenMP is off by default.
  if ($scm) {
    if ($openmp) {
      push @app, "fcflags_omp=$openmpflag\n";
      push @app, "ldflags_omp=$openmpflag\n";
    }

  } else {   # Not an SCM job

    unless ($openmp) {
      push @app, "fcflags_omp=\n";

#   On the IBM at UM 8.3 the OpenMP library flag is still needed because GCOM is
#   compiled with OpenMP. This was included in the x86 config but missed
#   from the pwr7 config. This was fixed in the config file at UM 8.4.
      unless ($version eq '8.3' and $platform eq 'ibm') {
        push @app, "ldflags_omp=\n";
      }
    }
  }


# 5. Any other settings:

# Set optimisation level 
# From UM 8.4 onwards this is included in the config name
  if ($version eq '8.3') {
    push @app, "opt_level=$opt\n" unless ($opt eq 'high');
  }

# Exec names (if not default or otherwise using $RUNID[.exe]):
  my $recon_exec = get_exec_name($jobdir, 'recon');
  my $atmos_exec = get_exec_name($jobdir, 'atmos');
  if ($recon_exec) {
    push @app, "recon_exec=$recon_exec\n" unless ($recon_exec eq 'qxreconf');
  }
  if ($atmos_exec) {
    push @app, "atmos_exec=$atmos_exec\n" unless ($atmos_exec =~ /\$RUNID/i);
  }

# Build steps required (default = all)
#  We always perform the extract since prebuild settings are not transferable.
#  SCM jobs contain only a single preprocess and build, so no checking required.
  unless ($scm) {
    if ($build_recon and not $build_atmos) {
      $steplist = "steplist=extract preprocess-recon build-recon\n";
    } elsif ($build_atmos and not $build_recon) {
      $steplist = "steplist=extract preprocess-atmos build-atmos\n";
    }
    $steplist =~ s/steplist=extract/mirror_steplist=/ if ($platform eq "ibm");
    push @app, $steplist if ($steplist);
  }


# To modify output text when we have multiple configs in play:
  my $plural = ($coupled) ? 's' : '';
  my $um = ($coupled) ? 'UM ' : '';

# Attempt to add any machine overrides...
  my @machine_overrides = get_machine_overrides($jobdir, $version);
  if (@machine_overrides) {
    push @app, @machine_overrides;
    print "[INFO] A machine override has been applied to the ${um}fcm_make config\n";
  }

# ...and any path overrides:
  my @path_overrides = get_path_overrides($jobdir, $envhashref);
  if (@path_overrides) {
    push @app, @path_overrides;
    print "[INFO] A path override has been applied to the fcm_make config$plural\n";
  }



# 6. Pick a platform config file(name):

  if ($platform eq "ibm") {
    $platform_config_dir = "meto-pwr7-xlf";
  } elsif ($platform eq "linux") {
    $platform_config_dir = "meto-x86-ifort";
  } else {
    die "[FAIL] Unrecognised platform, cannot provide a suitable config file.\n";
  }

  if ($version eq '8.3') {
    $model_config = 'um.cfg';
  } else {
    $model_config = "um-$model-$opt.cfg";
  }

  push @app,
    "include_config=$um_base/fcm-make/$platform_config_dir/$model_config\@$um_rev\n";


# 7. Add branches and working copies:
# All UM branches can be read from the FCM_UMATMOS_CFG file, even in the case
# of reconfiguration-only builds.
  foreach my $project (qw/ um jules /) {

    my @branches = get_branches($project, \@atmoscfg);
    my @working_copies = get_working_copies($project, \@atmoscfg);
    if (@working_copies) {
      print "[INFO] \U$project\E working copy in UM fcm-make app.\n";
      $unportable++;
    }

    my @diffs = ();
    push @diffs, @branches, @working_copies;


# Now add the branches/working copies to the app file:

# Always add the UM and JULES input lists even if empty, to be consistent
# with meta-data.
    push @app, "${project}_sources=@diffs\n";

# Add the extract to the make file; requires FCM 2014-01 if empty.
    push @make, "extract.location\{diff\}\[$project\] = \$${project}_sources\n";

  } # End loop over projects


# 8. Attempt to add any file overrides
  my @file_overrides = get_file_overrides($jobdir, $platform, $version);
  if (@file_overrides) {
    push @make, "\n", @file_overrides;
    print "[INFO] A file override has been applied to the ${um}fcm_make config\n";
  }

# Rosebud can only take a guess at how overrides should be translated for the
#  new config files; the user should check they are valid.
  if (@machine_overrides or @path_overrides or @file_overrides) {
    print "[INFO] Please check overrides are valid and compatible with the central configs.\n";
    if (@file_overrides) {
      print "[WARN] All file overrides affecting files used by both the atmosphere model\n       and the reconfiguration must be duplicated and applied to both builds.\n"
    }
  }

# Write the make and app files:
  fileio::write_file("$suitedir/app/$makeapp/file/fcm-make.cfg",@make);
  fileio::write_file("$suitedir/app/$makeapp/rose-app.conf", @app);

# And sort the app file:
  system "rose config-dump -qf $suitedir/app/$makeapp/rose-app.conf";

# Run the fixes macro if requested:
  macros::fixes("$suitedir/app/$makeapp", $macrofile, $meta) if ($fixes);

# And the upgrade macro, if requested:
  macros::upgrade("$suitedir/app/$makeapp", $macrofile, $meta, $upgrade)
    if ($upgrade);

# We return three items:
#  - A flag to indicate whether OpenMP is in use
#  - The location of the UM source for the ocean build configs
#  - The current tally of unportable settings
  return ($openmp, $um_base, $unportable);
}



sub get_branches {

# Subroutine to identify branches used in the job.
# Returns an array of branches for inclusion.

# Arguments
  my $project = shift;   # Project to examine (UM, JULES,...)
  my $arrayref = shift;
    my @config = @$arrayref;   # Contents of the FCM_<model>_CFG file

# Local variables
  my $line;
  my $branch;
  my $temp;

  my @branches = ();  # return array

# Branch information should always appear in path/version pairs,
#  so this should be safe...

  foreach $line (@config) {

# Get branch paths:
    if ($line =~ /repos::($project)::branch/i) {

# Copy it so that we can use the original line to scan for revisions later
      $temp = $line;
# Only branch name remains:
      $temp =~ s/repos::($project)::branch\d+\s*//i;

# FCM can deduce the path to the branch if it's stored in the project
# associated with that namespace and the FCM configuration file is set up 
# correctly (or an extract.location{primary} setting is provided).
# However, IOIPSL branches can come from either IOIPSL or NEMO repositories
# depending on the version, so we disable this feature for those.
      $temp =~ s!.*?/(dev|pkg|test)/!branches/$1/! 
        unless ($project eq 'ioipsl');

# Strip off any unwanted trailing directories:
      $temp =~ s{/src\s*$}{};# if ($project eq "um");
      $temp =~ s{(/NEMOGCM)?/NEMO\s*$}{};# if ($project eq "nemo");
      $temp =~ s{/NEMOGCM/EXTERNAL\s*$}{};# if ($project eq "ioipsl");# NEMO 3.4
      $temp =~ s{/cice\s*$}{};# if ($project eq "cice");

      chomp $temp;
      $branch = $temp;                # Save branch path
    }

# Get branch revisions:
    if ($line =~ /version::($project)::branch/i) {

      $temp = $line;
# Only revision remains:
      $temp =~ s/version::($project)::branch\w+\s*//i;
      chomp $temp;

      $branch = join "@", $branch, $temp;    # Append revision to branch path
      push @branches, $branch;               # Add branch@path to list
    }
  }
  return @branches;
}


sub get_working_copies {

# Subroutine to identify working copies used in the job.
# Returns an array (of at most one entry!) of working copies for inclusion.

# Arguments
  my $project = shift;    # Project to examine (UM, JULES,...)
  my $arrayref = shift;
    my @config = @$arrayref;  # Contents of the FCM_<model>_CFG file

# Local variables
  my $line;
  my @working_copies = ();  # return array

  foreach $line (@config) {

    if ($line =~ /repos::($project)::user/i) {

# Only branch name remains:
      $line =~ s/repos::($project)::user\s*//i;
# Strip off any unwanted trailing directories:
      $line =~ s{/src\s*$}{};# if ($project eq "um");
      $line =~ s{/NEMO\s*$}{};# if ($project eq "nemo");
      $line =~ s{/cice\s*$}{};# if ($project eq "cice");

      chomp $line;
      push @working_copies, $line;

    }
  }
  return @working_copies;
}


sub get_machine_overrides {

# Subroutine to read the list of machine overrides being used and convert
# them into a Rose config/FCM2-friendly format.
# All machine override flags are assumed to be for fortran compiler options and 
# are appended to the existing compiler flags.
# Consequently there is no support for passing on library flags or CPP/FPP
# keys - these will need manual treatment.

# Arguments
  my $jobdir = shift;    # Job input directory
  my $version = shift;   # UM version

  my @machovr = fileio::read_file("$jobdir/USR_MACH_OVRDS");


# Machine overrides:
# This UMUI file is a list of includes which are read in turn:
  my @arguments;
  foreach my $line (@machovr) {
    my @file = glob $1 if ($line =~ /^inc\s+(.*)/);
    my @overrides = fileio::read_file($file[0]);

#   Look for overrides to apply:
    foreach my $override (@overrides) {
      push @arguments, " $2" if ($override =~ /\s*(%\w+\s*)+(.*)/);
    }

  }

# Prepare the override for the app file:
  if (@arguments) {
    if ($version eq '8.3') {
      unshift @arguments, 'fc_overrides=';
    } else {
      unshift @arguments, 'fcflags_overrides=';
    }
    push @arguments, "\n";
  }

  return @arguments;
}

sub get_path_overrides {

# Subroutine to read the list of path overrides being used and convert
# them into a Rose config/FCM2-friendly format.

  my $jobdir = shift;   # Job input directory
  my $envhashref = shift;
    my %envvars = %{$envhashref}; # Hash of environment variables

  my @pathovr = fileio::read_file("$jobdir/USR_PATHS_OVRDS");


# Path overrides:
# A list in the form "%variable  /path/to/use"
  my @overrides;
  foreach my $line (@pathovr) {
    if ($line =~ /%(\w+)\s+(.*)/) {
      my $var  = $1;
      my @path = glob $2;    # Best done in list context
      my $path = $path[0];

# Look for any environment variables that need substituting:
      foreach my $variable (reverse sort keys %envvars) {
        $path =~ s/\$$variable/$envvars{$variable}/;
      }

# The coupled model always includes prism_home (which needs renaming).
# Although a default exists this is always printed both to aid the user and for
# use with the coupling flags, which are specified in the fcm-make.cfg file.
      $var = 'prism_path' if ($var eq 'prism_home');

      push @overrides, "$var=$path\n";
    }

  }

  return @overrides;
}

sub get_file_overrides {

# Subroutine to read the list of machine overrides being used and convert
# them into a Rose config/FCM2-friendly format.
# This will not be perfect but should at least give a fair approximation of
# what the translated override SHOULD look like, as an aid to the user.
# Note that even when it works it may be possible to structure the final 
# override better using those in the repository as examples, e.g. if appending
# a new optimisation level.

  my $jobdir = shift;
  my $platform = shift;
  my $version = shift;

# Set up a hash containing platform/optimisation specific settings:
  my %level;
  if ($platform eq 'ibm') {
    %level = (high => '-O3 -qstrict', 
              safe => '-O3 -qstrict',);
    $level{debug} = ($version eq '8.5') 
                    ? '-O0 -qfullpath -C -qinitauto=7FBFFFFF -qfloat=nans'
                    : '-O0 -qfullpath';

  } elsif ($platform eq 'linux') {
    %level = (high =>  '-O3 -fp-model precise',
              safe =>  '-O2 -fp-model precise',
             );
    $level{debug} = ($version eq '8.5') 
                    ? '-O2 -fp-model strict -C -ftrapuv -traceback'
                    : '-O0 -fp-model precise -traceback';
  }

# Set the machine overrides variable name:
  my $machine_ovr = ($version eq '8.3') ? '$fc_overrides' 
                                        : '$fcflags_overrides';

  my @fileovr = fileio::read_file("$jobdir/USR_FILE_OVRDS");
  my @override_list;


# File overrides:
# This UMUI file is a list of includes which are read in turn:
  my $debug_counter = 0;
  foreach my $line (@fileovr) {
    $line =~ s/#.*//;  # Not sure why anyone'd need this, but it's possible...
    my @file = glob $1 if ($line =~ /^inc\s+(.*)/);
    my @overrides = fileio::read_file($file[0]);

#   Look for overrides to apply:
    foreach my $override (@overrides) {
      $override =~ s/#.*//;     # Discard comments

# This list should be expandable if more varied conversions are needed
#  - Currently 1 entry: Fortran overrides (fflags64):
      if ($override =~ /^\s*bld::tool::fflags((::\w+)+)\s*(.*)/) {
        my $path = $1;
        my $ovr  = $3;

# Process the path
        $path =~ s/^:://;                    # Strip leading ::
        $path =~ s/(.*)/\L$1/g;              # make lowercase
        $path =~ s/^(um::)(.*)/$1src::$2/g;  # insert 'src' if for UM not JULES

# The new file overrides have to match the file extension (if present).
# We can't easily tell if an override is for a file or a directory, so as a 
# first approximation we assume anything nested 5 or more levels deep is a file.
# This may fail occasionally, e.g. for um/src/atmosphere/UKCA/photolib
        $path =~ s/((\w+::){4,}(\w+))/$1\.F90/;

        $path =~ s/::/\//g;         # switch :: to /

        push my @temp, "{fc.flags}[$path] =";   # construct (most of) LHS

# Process the override settings.

# First replace any special flags.
#
# Swap out opt_level flags, checking we have something to substitute it for! 
# $fcflags_level vars vary too much to be useful here.
        if ($ovr =~ /%fflags_opt_(\w+)/) {
          my $opt_level = $1;
          if (exists $level{$opt_level}) {
            $debug_counter++ if ($opt_level eq 'debug');
            $ovr =~ s/%fflags_opt_(\w+)/$level{$1}/;
          } else {
            print "[WARN] Cannot replace %fflags_opt_$opt_level file override; no substitute provided.\n";
          }
        }

# This should only be used to disable OpenMP in a complete override (not
# appending), so can be removed (since it should already contain an opt level):
        $ovr =~ s/%fflags_basic64//;

# By now hopefully the only '%' flag left is fflags64.
# Are we overriding all settings or appending a new argument?
        if ($ovr =~ /(%\w+)\s+\1/) {        # appending arguments

          $ovr =~ s/(%\w+\s+)+(.*)/$2/;
          push @temp, "\$fcflags_common \$fcflags_level \$fcflags_omp $machine_ovr $ovr";

        } else {                             # complete override

          $ovr =~ s/%\w+\s+(.*)/$1/;
          push @temp, "\$fcflags_common $ovr";

        }

# Attempt to apply each override to the correct build:
# We can safely apply any non-recon override to the atmos model. Guessing
# which atmos overrides need to be applied to the recon model is practically
# impossible; as a first approximation we assume they are mutually exclusive.
        if ($path =~ /um\/src\/utility\/qxreconf/) {
          push @override_list, "build-recon.prop@temp\n";
        } else {
          push @override_list, "build-atmos.prop@temp\n";
        }
        
      }  # End 'is a Fortran override' (fflags*)

    } # End loop over overrides in this file

  } # End loop over files

# Warn a MO user if they're using a flag that can change between versions
# (even if it doesn't change here, so they're prepared for an upgrade).
  print "[WARN] Found a file override using %fflags_opt_debug.\n"
      . "[WARN] Debug-level fcm make settings vary between UM 8.4 and UM 8.6.\n"
    if ($debug_counter);

  return @override_list;
}


sub process_base_rev {

# Process the base source revision number/keyword to ensure it
# matches the forms allowed by the meta-data.

  my $rev = shift;    # revision keyword/number

  $rev =~ s/^HEAD$/HEAD/i;       # Use the FCM style of uppercase built-ins
  $rev =~ s/^(vn|um)/\L$1/i;     # Project-defined keywords are lowercase

  return $rev;
}

sub get_exec_name {

# Returns the name of the executable of the supplied job component.
# If no build target is available (no build), returns an empty string.
# Reads the required file directly so multiple file creation subroutines can
# use it.

# Arguments:
  my $jobdir = shift;     # UMUI output directory
  my $target = shift;     # atmos, recon or ocean

# Return variable
  my $exec = '';          # Executable name

  my $filename = 
    ($target eq 'atmos') ? 'FCM_UMATMOS_CFG'  :
    ($target eq 'recon') ? 'FCM_UMRECON_CFG'  :
    ($target eq 'ocean') ? 'FCM_NEMOCICE_CFG' :
                           '';     # Blank default

  die "[FAIL] Cannot find config file for unknown executable build"
    unless ($filename);

  my @config = fileio::read_file("$jobdir/$filename");

  foreach my $line (@config) {
    $line =~ s/#.*//;     # Discard comments
    $exec = $1 if ($line =~ /bld::exe_name::\w+\s*([\$\w.]+)/);
  }

  return $exec;
}

1;
