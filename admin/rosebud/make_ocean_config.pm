# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************
#
# Script:     make_ocean_config.pm
# Purpose:    Write the fcm_make_ocean config file and associated app
# Code owner: UM System Development team

package make_ocean_config;

use strict;
use warnings;
use fileio;
use make_config;
use macros;

sub write_fcm_make_ocean {

# Write the fcm_make_ocean.cfg config and its rose-app.conf file.

  my $hashref = shift;

# Transfer relevant hash values to variables for ease of use
  my $jobdir  = $hashref->{jobdir};            # Job input directory
  my $suitedir = $hashref->{suitedir};         # Suite output directory
  my $platform = $hashref->{platform};         # IBM or Linux (or ...)?
  my $build_atmos = $hashref->{build_atmos};   # True = atmos build required
  my $build_recon = $hashref->{build_recon};   # True = recon build required
  my $makeapp = $hashref->{makeapp};           # Name of the atmos fcm_make app
  my $oceanmakeapp = $hashref->{oceanmakeapp}; # Name of the ocean fcm_make app
  my $um_base = $hashref->{um_base};           # Location of the UM 'trunk'
  my $unportable = $hashref->{unportable};     # No. of unportable settings
  my $scm = $hashref->{scm};                   # True = SCM job
  my $fixes = $hashref->{fixes};               # True = apply fixes macro
  my $upgrade = $hashref->{upgrade};           # Version to upgrade to
  my $macrofile = $hashref->{macrofile};       # File containing macro output
  my $meta = $hashref->{meta};                 # Path to meta-data
  my @extr_scr = @{$hashref->{extr_scr}};      # Contents of EXTR_SCR
  my $envhashref = $hashref->{envhashref};     # Reference to hash of 
                                               #  environment variables
    my %envvars = %{$envhashref};              # Hash of environment variables

# Local variables:
  my $nemo_openmp;           # OpenMP on or off in NEMO build
  my $cice_openmp;           # OpenMP on or off in CICE build
  my $openmpflag;       # Compiler flag for OpenMP
  my $model = "coupled";# String containing app-specific name (atmos, scm etc.)
  my $ocean_exec = '';  # Name of ocean executable
  my $line;
  my $item;
  my $override;         # An override for the config file
  my $platform_config_dir;  # Platform-specific config directory
  my $model_config;         # Opt-specific config name (nemo-cice-safe, etc.)
  my $um_rev  = "vn8.6";  # UM revision for ocean fcm config - hardwired to 8.6.
  my @extract_location;
  my @make;      # fcm_make_ocean.cfg contents
  my @app;       # rose-app.conf contents

# Some defaults for (potential) printing to file:
  my $nemo_base = "fcm:nemo_tr/NEMO";
  my $nemo_rev  = "vn3.2";
  my $ioipsl_base = "fcm:ioipsl_tr";
  my $ioipsl_rev  = "vn3.0";
  my $cice_base = "fcm:cice_tr/cice";
  my $cice_rev  = "vn4.1";
  my $opt = "safe";             # Only available setting for now

  if ($platform eq 'ibm') {
    $openmpflag = '-qsmp=omp';
  } elsif ($platform eq 'linux') {
    $openmpflag = '-openmp';
  }

# Set up the app file header.
# Meta-data is hard-wired to vn8.6 (the earliest available).
  push @app, "meta=nemo-cice-fcm-make/vn8.6\n\n[env]\n";

# Which base version of NEMO are we using?
  my $ocean_version_file = ($envvars{NEMO_VERSION} eq '302')
                           ? 'nemo3.2-cice4.1-version.cfg'
                           : 'nemo3.4-cice4.1m1-version.cfg';
  push @app, "ocean_version_file=$ocean_version_file\n";
  

# Read the Ocean FCM config
  my @oceancfg = fileio::read_file("$jobdir/FCM_NEMOCICE_CFG");

# 1. Create list of FPP keys for NEMO and CICE:

  my @nemokeys = get_ocean_keys('nemo', \@oceancfg);
  my @cicekeys = get_ocean_keys('cice', \@oceancfg);

  push @app, "keys_nemo_app=@nemokeys\n";
  push @app, "keys_cice_app=@cicekeys\n";


# 2. Set source and revision information:
# The Ocean configs are hardwired to extract the code from given NEMO and CICE
# subdirectories so we don't worry about NEMO_SVN_BASE (="NEMO code location"),
# additional src::nemo::base locations (="Additional NEMO source directories"),
# etc. If those change more fundamental config changes are needed anyway.
  foreach $line (@extr_scr) {

# Already have these in the hash as pure numbers, 
# but need valid keywords/revisions
    $nemo_rev = $1 if ($line =~ /export NEMO_VER=([\w.]+)/);
    $cice_rev = $1 if ($line =~ /export CICE_VER=([\w.]+)/);
    $ioipsl_rev = $1 if ($line =~ /export IOIPSL_VER=([\w.]+)/);

    $nemo_base = $1 if ($line =~ m!export NEMO_SVN_URL=(.*)(/NEMO)?!);
    $cice_base = $1 if ($line =~ m!export CICE_SVN_URL=(.*)(/cice)?!);
    $ioipsl_base = $1 if ($line =~ m!export IOIPSL_SVN_URL=(.*)(/NEMOGCM)?!);
  }

# Source information is only output if not the trunk. Revision information is 
# more flexible and always output to make clear what's being used.

    push @app, "nemo_base=$nemo_base\n" unless
            ($nemo_base =~ m!fcm:nemo[_-]tr! or
             $nemo_base =~ m!fcm:nemo/trunk! or
             $nemo_base =~ m!svn://fcm\d/NEMO_svn/NEMO/trunk!);
    $nemo_rev = make_config::process_base_rev($nemo_rev);
    push @app, "nemo_rev=$nemo_rev\n";

    if ($envvars{NEMO_VERSION} eq '302') {
      push @app, "ioipsl_base=$ioipsl_base\n" unless
              ($ioipsl_base =~ m!fcm:ioipsl[_-]tr! or
               $ioipsl_base =~ m!fcm:ioipsl/trunk! or
               $ioipsl_base =~ m!svn://fcm\d/NEMO_svn/IOIPSL/trunk!);
    } else {
      push @app, "ioipsl_base=$ioipsl_base\n" unless
              ($ioipsl_base =~ m!fcm:nemo[_-]tr! or
               $ioipsl_base =~ m!fcm:nemo/trunk! or
               $ioipsl_base =~ m!svn://fcm\d/NEMO_svn/NEMO/trunk!);
    }
    $ioipsl_rev = make_config::process_base_rev($ioipsl_rev);
    push @app, "ioipsl_rev=$ioipsl_rev\n";

    push @app, "cice_base=$cice_base\n" unless
            ($cice_base =~ m!fcm:cice[_-]tr!      or
             $cice_base =~ m!fcm:cice/trunk!      or
             $cice_base =~ m!svn://fcm\d/CICE_svn/CICE/trunk!);
    $cice_rev = make_config::process_base_rev($cice_rev);
    push @app, "cice_rev=$cice_rev\n";

# 3. Additional miscellaneous settings required before main config inclusion:

# UMUI jobs are effectively hardwired to 'safe' optimisation, so we mimic that
# behaviour here. From UM 8.6 Ocean builds have access to the same range of
# optimisation levels as the UM, albeit some are untested.
# No steplist to update as only one component to preprocess and build
# No machine overrides to update as those in the UMUI are for UM builds only.

# Exec name (if not default):
  $ocean_exec = make_config::get_exec_name($jobdir, 'ocean');
  if ($ocean_exec) {
    my $default_ocean_exec = ($envvars{NEMO_VERSION} eq '302') ? 'model.exe' 
                                                               : 'nemo.exe';
    push @app, "ocean_exec=$ocean_exec\n" 
      unless ($ocean_exec eq $default_ocean_exec);
  }

# 4. Read NEMO and CICE compiler flags files:
#
# NEMO and CICE settings are provided via the compiler flags file, with some
# extra flags hardwired into the UMUI and selected via the OASIS switch.
#
# We analyse the build flags (relative to the trunk configs) to see if any
# overrides are required.
# This is not a full analysis; we make use of the fact that the UMUI processes 
# the NEMO/CICE build flags files in a fixed way. If the user has abused
# hand-edits to rewrite the config files rather than providing a new build 
# flags file then results may not be as expected.
#
# Of necessity, this analysis is currently tied to the IBM compiler flags.

# Grab the two build config files that contain the main compiler settings:
  my $nemo_flags_file = get_ocean_flags_file('NEMO',\@oceancfg);
  my $cice_flags_file = get_ocean_flags_file('CICE',\@oceancfg);

  die "[FAIL] Unable to find NEMO compiler flags file in UMUI job\n"
    unless ($nemo_flags_file);
  die "[FAIL] Unable to find CICE compiler flags file in UMUI job\n"
    unless ($cice_flags_file);

# Create blank variables for each type of flag.
# Set before the main include:
  my $cpp;
  my $cc;
  my $fc;
  my $nemo_fflags;    # Applied as a (nemo) global override to fcflags_nemo_all
  my $cice_fflags;    # Applied as a (cice) global override to fcflags_cice_all
  my $flags_coupling;
  my $ldflags_coupling;

# Set after the main include:
  my $nemo_fppflags;  # Main ocean FPP settings
  my $cice_fppflags;  # Needs a separate namespaced override if different
  my $cppflags;       # Only one set shared between sub-models?
  my $ldflags;        # Constructed from many sources

# Needed if %netcdf/[inc/lib] is referred to in an override:
  my $netcdf_path;


# Read through the NEMO flags file to find the main settings:
  my @nemo_flags = fileio::read_file($nemo_flags_file);

  foreach $line (@nemo_flags) {
    $line =~ s/#.*//;                                # Discard comments
    $netcdf_path = tidy_flags($1) if ($line =~ /^\s*%netcdf\s+(.*)/);
    $cpp = tidy_flags($1) if ($line =~ /bld::tool::cpp\s+(.*)/);
    $cc = tidy_flags($1) if ($line =~ /bld::tool::cc\s+(.*)/);
    $fc = tidy_flags($1) if ($line =~ /bld::tool::fc\s+(.*)/);
    $cppflags = tidy_flags($1) if ($line =~ /bld::tool::cppflags\s+(.*)/);
    $nemo_fppflags = tidy_flags($1) if ($line =~ /bld::tool::fppflags\s+(.*)/);
    $nemo_fflags = tidy_flags($1) if ($line =~ /bld::tool::fflags\s+(.*)/);
    $ldflags = tidy_flags($1) if ($line =~ /bld::tool::ldflags\s+(.*)/);
  }

# Coupling flags come from the UMUI's config file:

  my $search = 0;
  my $append = 0;
  my $oasis = 0;
  foreach $line (@oceancfg) {
    $search = 1 if ($line =~ /# OASIS Coupling/);
    if ($search) {
      $line =~ s/#.*//;                                     # Discard comments
# Need a second search flag because this runs over multiple lines:
      if ($line =~ /(bld::tool::cflags)\s+%\1/) {
        $append = 1;
        $line =~ s/(bld::tool::cflags)\s+%\1\s+//; # Free up flags on this line
      }
      if ($append) {
# If the flags continue past this line, remove the \ and keep looking.
# Else, stop adding things.
        if ($line =~ /\\\n/) {
          $line =~ s/\\\n//;
        } else {
          $append = 0;
          $search = 0;
        }
#       Append the line to the list of coupling flags:
        $flags_coupling .= " $line";
      } # End if (append flags)
    } # End if (search file for key phrase)


# Separate to the above, grab the library coupling flags (OASIS load flags):
# Third search flag required!
    $oasis = 1 if ($line =~ /# OASIS load flags/);
    if ($oasis) {
      $line =~ s/#.*//;                                     # Discard comments
      $ldflags_coupling .= " $2"
        if ($line =~ /(bld::tool::ldflags)\s+%\1\s+(.*)/);
# Don't need to disable this search - it should safely run to the end.
    }


  } # End loop over file

  $flags_coupling   = tidy_flags($flags_coupling);
  $ldflags_coupling = tidy_flags($ldflags_coupling);

# That gives us the main NEMO-centric compiler flags, in some fashion.
# Now do it again for CICE. Settings may override NEMO values.
  my @cice_flags = fileio::read_file($cice_flags_file);

# Ignore most namespaced settings; they're assumed to be for standalone builds.
# The main CICE fortran compiler flags are the notable exception.
  foreach $line (@cice_flags) {
    $line =~ s/#.*//;                                # Discard comments
    $netcdf_path = tidy_flags($1) if ($line =~ /^\s*%netcdf\s+(.*)/);
    $cpp = tidy_flags($1) if ($line =~ /bld::tool::cpp\s+(.*)/);
    $cc = tidy_flags($1) if ($line =~ /bld::tool::cc\s+(.*)/);
    $fc = tidy_flags($1) if ($line =~ /bld::tool::fc\s+(.*)/);
    $cppflags = tidy_flags($1) if ($line =~ /bld::tool::cppflags\s+(.*)/);
    $cice_fppflags = tidy_flags($1) if ($line =~ /bld::tool::fppflags\s+(.*)/);
    $cice_fflags = tidy_flags($1) if ($line =~ /bld::tool::fflags::cice\s+(.*)/);# CICE only
    $ldflags = tidy_flags($1) if ($line =~ /bld::tool::ldflags\s+(.*)/);
  }

# Get OpenMP settings:
  if ($nemo_fflags) {
    $nemo_openmp = 1 if ($nemo_fflags =~ /$openmpflag/);
  }
  if ($cice_fflags) {
    $cice_openmp = 1 if ($cice_fflags =~ /$openmpflag/);
  }
# Do we need to disable OpenMP?
  push @app, "fcflags_omp=\n" unless ($nemo_openmp);
  push @app, "ldflags_omp=\n" unless ($nemo_openmp or $cice_openmp);

# All flags have been read in. Now we compare them one by one against the 
# defaults on the trunk.

  my @pre_include;    # Overrides that go before the include = ...
  my @post_include;   # Overrides that go after the include = ...

# For each option:
#  i) Set up a list containing the default (trunk) flags
#  ii) Compare the new options with the defaults (&compare_flags)
#  iii) Return a string containing the new options if needed, empty if not
#  iv) Add the override to the correct place in the file if non-empty

# Set up defaults:
  my @default_cpp;
  my @default_cc;
  my @default_fc;
  my @default_cppflags;
  my @default_nemo_fppflags;
  my @default_nemo_fflags;
  my @default_cice_fppflags;
  my @default_cice_fflags;
  my @default_ldflags;
  my @default_flags_coupling;
  my @default_ldflags_coupling;
  my $check_flags = 1;

  if ($platform eq 'ibm') {
# Despite the defaults in the configs we print some compiler settings
# to aid developers who then have them visible and available to modify - now
# that most settings are stored in the app file they're much easier to manage
# via the GUI and all app variables have help and/or a description in the 
# meta-data.
#
# The exceptions are genuinely basic settings, such as the compilers
# and pre-processors, for which the defaults are usually adequate, or flags
# for which the coupling flags are the main component (these are overridden
# separately).
#
# We do however update any old or deprecated settings up to the new Rose 
# standards (using the ocean configs in UM vn8.6).

# Lists of possible defaults, where we only expect to match one item ('or'):
    @default_cpp = 'xlc';
    @default_cc  = 'xlc_r';
    @default_fc  = qw/ mpxlf90_r mpxlf95_r /;

# List of all defaults, where we may match every item one-to-one ('and'):
    @default_cppflags = qw/ -E -C /;
    @default_nemo_fppflags = qw/ -E -P -traditional-cpp /;
    @default_nemo_fflags = qw/ -qrealsize=8 -qextname -qsuffix=f=f90
                               -qarch=pwr7 -qtune=pwr7 -NS32768 /;
    @default_cice_fppflags = @default_nemo_fppflags;   # Same for both
    @default_cice_fflags = qw/ -qextname -qsuffix=f=f90 -qarch=pwr7 -qtune=pwr7
                               -qfree=f90 -g /;
    @default_ldflags = qw/ /; # Empty; default values are added later by UMUI
    @default_flags_coupling = qw{ -I$prism_path/build/lib/mpp_io 
                                  -I$prism_path/build/lib/psmile.MPI1 
                                  -I$prism_path/build/mod/oasis3.MPI1 };
    @default_ldflags_coupling = qw{ -L$prism_path/lib 
                                    -lanaisg -lanaism -lpsmile.MPI1
                                    -lfscint -lmpp_io -lscrip };
    $netcdf_path = '/projects/um1/lib/netcdf3.20090102' unless ($netcdf_path);


  } else {
# Here be dragons...
    print "[WARN] No checking of ocean compiler flags performed on platform $platform\n";
    $check_flags = 0;
  }

# Platform-independent additions:
  push @default_nemo_fflags, $openmpflag if ($nemo_openmp);
  push @default_ldflags, $openmpflag if ($nemo_openmp or $cice_openmp);


  if ($check_flags) {

# Go through the compiler settings one-by-one, stripping out or updating
# any old values.

### App file variables ###
# Compiler definitions rarely need changing so remove any defaults:
    if ($cpp) {
      foreach my $opt (@default_cpp) {
        $cpp =~ s/$opt//g;
      }
      push @app, "cpp=$cpp\n" unless ($cpp =~ /\s*/);
    }

    if ($cc) {
      foreach my $opt (@default_cc) {
        $cc =~ s/$opt//g;
      }
      push @app, "cc=$cc\n" unless ($cc =~ /\s*/);
    }
  
    if ($fc) {
      foreach my $opt (@default_fc) {
        $fc =~ s/$opt//g;
      }
      push @app, "fc=$fc\n" unless ($fc =~ /\s*/);
    }

# The two fflags overrides are appended onto existing flags rather than 
# replacing them, so we should only take the extra flags that are needed.
    if ($nemo_fflags) {
      $override = compare_flags($nemo_fflags, \@default_nemo_fflags);
      if ($override) {
        $override = extract_flags($override, \@default_nemo_fflags);
        push @app, "fcflags_nemo_overrides=$override\n";
      }
    }

    if ($cice_fflags) {
      $override = compare_flags($cice_fflags, \@default_cice_fflags);
      if ($override) {
        $override = extract_flags($override, \@default_cice_fflags);
        push @app, "fcflags_cice_overrides=$override\n";
      } # End if (override)
    } # End if (cice_fflags)


### Pre-includes ###
# These may use variables defined in the app file but not the central configs.

# Despite the existence of default values we always print the following two
# settings so they are visible to the user, to aid code development.
# (The first of these appears to be dependent only on the version of OASIS.)
#  - We must therefore print prism_path to the app file too.

    $flags_coupling = join ' ', @default_flags_coupling unless ($flags_coupling);
    push @pre_include, "\$flags_coupling = $flags_coupling\n";


    $ldflags_coupling = join ' ', @default_ldflags_coupling unless ($ldflags_coupling);
# NetCDF flags belong in $ldflags_ocean, not here:
    $ldflags_coupling =~ s!-L\$netcdf_[/\w]+!!g;
    $ldflags_coupling =~ s/-lnetcdf\b//g;
# Tidy whitespace:
    $ldflags_coupling =~ s/\s+/ /g;
    $ldflags_coupling =~ s/^\s+//;
    $ldflags_coupling =~ s/\s+$//;
    push @pre_include, "\$ldflags_coupling = $ldflags_coupling\n";

# This also needs printing to the UM's fcm-make.cfg file.
# We can do this here as any macros that have already been run (e.g. fixes, 
# stash) will only have operated on the app file.
    my @makeatmos = fileio::read_file("$suitedir/app/$makeapp/file/fcm-make.cfg");
    foreach $line (@makeatmos) {
      if ($line =~ /^include = /) {
        $line = "\$ldflags_coupling = $ldflags_coupling\n\n" . $line;
      }
    }
    fileio::write_file("$suitedir/app/$makeapp/file/fcm-make.cfg",@makeatmos);
    

### Post-includes ###
# Only the coupling flags really change here; the rest rarely differs,
# so we probably don't need to supply these unless genuinely necessary.
#
# NEMO- or CICE-only settings must be namespaced (or else CICE settings must be
# (re-)printed after any NEMO overrides).
    if ($cppflags) {
      $override = compare_flags($cppflags, \@default_cppflags);
      push @post_include, 
        "preprocess-ocean.prop{cpp.flags} = $override \$cppflags_coupling\n"
        if ($override);
    }

    my $nemo_override = '';
    if ($nemo_fppflags) {
      $override = compare_flags($nemo_fppflags, \@default_nemo_fppflags);
      if ($override) {
        push @post_include, "preprocess-ocean.prop{fpp.flags} = $override \$fppflags_coupling\n";
        $nemo_override = $override;
      }
    }

# CICE flags are a namespaced subset of the NEMO flags, which are applied 
# globally. Therefore we only need to override CICE flags if 
# a) $nemo_fppflags exists, and is different from $cice_fppflags, or
# b) $nemo_fppflags doesn't exist, and $cice_fppflags differs from the default.
# We compare the nicely processed overrides instead of the raw flags so that we
# can do direct string comparison:
    if ($cice_fppflags) {
      $override = compare_flags($cice_fppflags, \@default_cice_fppflags);
      if (   ($nemo_fppflags and ($nemo_override ne $override)) 
          or ($override and not $nemo_fppflags) ) {
        push @post_include, "preprocess-ocean.prop{fpp.flags}[cice] = $override \$fppflags_coupling\n";
      }
    }

    if ($ldflags) {
      $override = compare_flags($ldflags, \@default_ldflags);
      push @post_include, "build-ocean.prop{fc.flags-ld} = -L\$netcdf_lib_path -lnetcdf \$ldflags_omp $override \$ldflags_coupling\n" if ($override);
    }

  } # End if (check_flags)


# 5. Apply overrides to config file:

# For path overrides, we apply the same values as for the UM.
  my @path_overrides = make_config::get_path_overrides($jobdir, $envhashref);
  push @app, @path_overrides if (@path_overrides);


# It's possible that the compiler flag files have created a pre-include override
# without netcdf_inc_path and netcdf_lib_path being set. If so we need to
# plug that gap.
    my $netcdf_path_ovr = 0;
    my $netcdf_inc_path_ovr = 0;
    my $netcdf_lib_path_ovr = 0;
    foreach $line (@path_overrides) {
      $netcdf_path_ovr = 1 if ($line =~ /\$netcdf_path/);
      $netcdf_inc_path_ovr = 1 if ($line =~ /\$netcdf_inc_path/);
      $netcdf_lib_path_ovr = 1 if ($line =~ /\$netcdf_lib_path/);
    }

# Try to keep the final file as neat as possible, so use the 
# _lib_path / _inc_path overrides if we have them.
# Otherwise, substitute for netcdf_path as is.
    my $use_netcdf_path = 0;
    foreach $line (@pre_include) {
      if ($line =~ m!\$netcdf/inc!) {
        if ($netcdf_inc_path_ovr) {
          $line =~ s!\$netcdf/inc\w*!\$netcdf_inc_path!g;
        } else {
          $line =~ s!\$netcdf/!\$netcdf_path/!g;
          $use_netcdf_path = 1;
        }
      } # End if ($netcdf/inc)
      if ($line =~ m!\$netcdf/lib!) {
       if ($netcdf_lib_path_ovr) {
          $line =~ s!\$netcdf/lib\w+!\$netcdf_lib_path!g;
        } else {
          $line =~ s!\$netcdf/!\$netcdf_path/!g;
          $use_netcdf_path = 1;
        }
      } # End if ($netcdf/lib)
    } # End loop over @pre_include

# Lastly for netcdf paths, add the $netcdf_path override itself if needed:
  push @app, "netcdf_path=$netcdf_path\n" 
    if ($use_netcdf_path and not $netcdf_path_ovr);

# Set up the first part of the fcm-make.cfg file:
  push @make, @pre_include, "\n"  if (@pre_include);  # Add pre-include overrides
  push @make, "include = \$include_config\n\n";
  push @make, @post_include, "\n" if (@post_include); # Add post-include overrides


# 6. Pick a platform config file(name):

  if ($platform eq "ibm") {
    $platform_config_dir = "meto-pwr7-xlf";
  } elsif ($platform eq "linux") {
    $platform_config_dir = "meto-x86-ifort";
  } else {
    die "[FAIL] Unrecognised platform, cannot provide a suitable config file.\n";
  }

  $model_config = "nemo-cice-$opt.cfg";

  push @app,
    "include_config=$um_base/fcm-make/$platform_config_dir/$model_config\@$um_rev\n";


# Despite my best efforts, I'd always encourage the user to check these...
  print "[INFO] Check the $oceanmakeapp app settings are suitable before compiling.\n";

# 7. Add branches and working copies:
  my @modifications = ();
  foreach my $project (qw/ nemo ioipsl cice /) {

    my @branches = make_config::get_branches($project, \@oceancfg);
    my @working_copies = make_config::get_working_copies($project, \@oceancfg);
    if (@working_copies) {
      print "[INFO] \U$project\E working copy in Ocean config file.\n";
      $unportable++;
    }

    my @diffs = ();
    push @diffs, @branches, @working_copies;

# Practically all NEMO 3.2 jobs will need Richard's branch to fix the #include
# of mpif.h, as we've removed the include paths for it:
    if ($project eq 'nemo' and $envvars{NEMO_VERSION} eq '302') {

# Don't need it if it's already there:
      my $add_branch = 1;
      foreach my $diff (@diffs) {
        $add_branch = 0 if ($diff =~ /vn3.2_mpif_for_mpxlf95/);
      }
      push @diffs, 'branches/dev/frrh/vn3.2_mpif_for_mpxlf95@9574' 
        if ($add_branch);
    } # End special NEMO 3.2 treatment


# Now add the branches/working copies to the app file:

# Always add the input lists even if empty, so we can adjust the meta-data
# to require them in future if necessary.
    push @app, "${project}_sources=@diffs\n";

# Add the extract to the make file; requires FCM 2014-01 if empty.
    push @make, "extract.location\{diff\}\[$project\] = \$${project}_sources\n";

  }  # end loop over projects


# No check for file overrides; these are expected to be part of the NEMO/CICE 
# compiler flags files.


# 8. Look for any extra NEMO source directories that are needed.
# CICE is all but fixed, as is IOIPSL, so we leave those as the default.

  my @nemoextract = fileio::read_file("$jobdir/FCM_EXTRACT_NEMO_CFG");
  my @extracts = qw/ OPA_SRC /;     # Always required
  foreach $line (@nemoextract) {
    $line =~ s/#.*//;         # Discard comments
    if ($line =~ /^src::nemo::base\s*(\w+)/) {
      push @extracts, $1
    }
  }
  @extracts = sort(@extracts);

# @extracts now contains a sorted list of source directories.
# We could test against the default values and decide if the app needs this 
# setting, but because
# a) the default NEMO 3.4 $nemo_path_incl list is wrong at UM 8.6, and
# b) to show developers that they now need to include OPA_SRC in the list
# we include it in every app regardless.
  my $tree = ($envvars{NEMO_VERSION} eq '302') ? 'NEMO' : 'NEMOGCM/NEMO';

  foreach my $dir (@extracts) {
    $dir = "$tree/$dir";
  }
  push @app, "nemo_path_incl=@extracts\n";


# 9. Finish up.

# Write the make and app files:
  fileio::write_file("$suitedir/app/$oceanmakeapp/file/fcm-make.cfg",@make);
  fileio::write_file("$suitedir/app/$oceanmakeapp/rose-app.conf", @app);

# And sort the app file:
  system "rose config-dump -qf $suitedir/app/$oceanmakeapp/rose-app.conf";

# Run the fixes macro if requested:
  macros::fixes("$suitedir/app/$oceanmakeapp", $macrofile, $meta) if ($fixes);

# And the upgrade macro, if requested:
  macros::upgrade("$suitedir/app/$oceanmakeapp", $macrofile, $meta, $upgrade)
    if ($upgrade);

# Return the current tally of unportable items:
  return $unportable;
}


#######################################################################


sub get_ocean_keys {

# Read the FPP (etc.) keys of an ocean sub-model from the FCM_NEMOCICE_CFG file

# Arguments
  my $model = shift;               # nemo or cice
  my $arrayref = shift;
    my @oceancfg = @$arrayref;     # Contents of FCM_NEMOCICE_CFG

# Return array
  my @tidykeys = "";              # Final list of FPP keys

# Local variables
  my $line;                        # Loop variable
  my $search = 0;                  # Flag for searching for keys
  my @includefile;                 # An include file containing further keys
  my @include = "";                # Include file contents
  my @keys = "";                   # (Temporary) list of FPP keys
  my @defaults = "";               # List of keys included in the trunk config

  foreach $line (@oceancfg) {
    $line =~ s/#.*//;                                # Discard comments

    $search = 1 if ($line =~ /%${model}_jobdefs/);   # Turn search on
    if ($search) {
      push @keys, $1 if ($line =~ /^\s*([\w=]+)\s*\\?\n/); # Add key

      @includefile = glob $1 if ($line =~ /\s*inc\s*(.*)\s*/);# Get include file
      if (@includefile) {
        @include = fileio::read_file($includefile[0]); # Read include file
        last;    # Done searching
      } # End if (includefile)
    } # End if (search)
  } # End of loop over UMUI-generated config

# Go through the include file for further keys
  $search = 0;
  foreach $line (@include) {

#   If the line only contains comments, skip to the next line
#   (so we don't break the system of checking for continuation lines):
    next if ($line =~ /^\s*#/);

    $line =~ s/#.*//;                                       # Discard comments
    $search = 1 if ($line =~ /bld::tool::fppkeys::$model/); # Search on
    $line =~ s/bld::tool::fppkeys::$model//;     # Free up key on this line

# Can't guarantee contents of this file, 
# so read only as long as there are continuation characters present.
    if ($search) {
      push @keys, $1 if ($line =~ /^\s*([\w=]+)\s*\\?\n/); # Add key
      chomp $line;
      $search = 0 unless ($line =~ /\\$/);   # Stop unless continuation
    }
  } # End loop over include file


# Remove any defaults already in the nemo-cice config:
  if ($model eq 'nemo') {
    @defaults = qw(key_mpp_mpi key_cice);
  } elsif ($model eq 'cice') {
    @defaults = qw(LINUX ORCA_GRID CICE_IN_NEMO coupled ncdf);
  }

  foreach my $key (@keys) {
    foreach my $default (@defaults) {
      $key = "" if ($key eq $default);
    }
  }

# Cleaning up: remove blank entries
  foreach my $key (@keys) {
    push @tidykeys, $key if ($key);
  }

# Return the list of FPP keys
  return @tidykeys;
}


sub get_ocean_flags_file {

# Returns the file containing the compiler flags for a given Ocean sub-model.

# Arguments
  my $model = shift;               # NEMO or CICE (upper case!)
  my $arrayref = shift;
    my @oceancfg = @$arrayref;     # Contents of FCM_NEMOCICE_CFG

# Local variables
  my $search = 0;                  # Flag for searching for include files
  my @path;                        # Raw path to include file

# Assume the first 'inc' command after the relevant comment is for the include
# file we're after:
  foreach my $line (@oceancfg) {

    if ($search) {
      if ($line =~ /^\s*inc\s*(.*)\s*/) {
        @path = glob $1;
        $search = 0;
      }
    }
    $search = 1 if ($line =~ /# Define $model b(ui)?ld declarations/);
  }

  return $path[0];
}


sub tidy_flags {

# Tidy up a string of compiler settings into a form that's ready for
# further processing, updating or removing obsolete, deprecated or
# problematic flags. The tidied string is returned.

  my $flags = shift;     # Flags read in as a single string

# Stuff that needs killing off or changing:
  $flags =~ s/-berok\b//g;               # Harmful!
  $flags =~ s/prism_home/prism_path/g;   # Standardised name
  $flags =~ s/netcdf_home/netcdf_path/g; # Standardised name
  $flags =~ s/%/\$/g;                    # Switch from FCM1 to FCM2 syntax
  $flags =~ s/-O\d\b//g;                 # Set via choice of config file
  $flags =~ s/-qstrict\b//g;             # Included in the optimisation flags
  $flags =~ s/-traditional\b/-traditional-cpp/g; # Old version is deprecated
  $flags =~ s/=pwr6/=pwr7/g;             # Current, better optimised settings
                                         #   (small risk of change in answers)

# We need to stop explicitly setting the include paths for mpif.h files
# and addressing the root problem (using C-style #include statements instead of
# a Fortran INCLUDE).
  $flags =~ s!-I/opt/ibmhpc/pecurrent[/\.\w]+thread64!!g;

# Clean up whitespace.
# (shortens all entries to '' if the above has already removed everything)
  $flags =~ s/\s+/ /g;
  $flags =~ s/^\s+//;
  $flags =~ s/\s+$//;

  return $flags;
}


sub compare_flags {

# Compare the provided compiler options with some default, and return
# a string containing whatever settings are needed.
# Returning a blank string indicates no override is needed.
#
# This analysis assumes:
# i) The input flags are set up sensibly (i.e. no duplication of items)
# ii) Any 'extra' items added at the move to Rose have been appended to the
# input list before entering this routine.

  my $newflags = shift;   # Flags read in as a single string
  my $arrayref = shift;
    my @defaults = @$arrayref;     # List of default flags

# Create an array of flags, and a copy for later.
  my @flags = split / /, $newflags;

# Three (overlapping) possibilities:
#  - option has extra flags (regardless of whether it's missing others)
#  - option has fewer flags (with no extras)
#  - flags are the same as the default

# If they have different numbers of things in them, they must differ;
# return immediately.
  return $newflags if (@flags != @defaults);

# If we're here, there are the same numbers of options in both
# and so either they're the same as the default or there's a new flag present.
# Look for extra options by matching and removing entries.
  my $counter = 0;
  foreach my $flag (@flags) {
    foreach my $default (@defaults) {
      $flag = '' if ($flag eq $default);
    }
    $counter++ if ($flag ne '');
  }

# If anything was left, return the new flags we read in.
  return $newflags if ($counter);

# If we made it to here, presumably we match the default.
# Return an empty string.
    return '';
}


sub extract_flags {

# Extract the new flags in an override that aren't present in the default list.

  my $override = shift;    # String containing processed overrides
  my $arrayref = shift;
    my @defaults = @$arrayref;     # List of default flags

  my $new_flags;           # Return value; string of flags to be appended

# Run through the list, deleting anything present in the default set:
  my @override = split / /, $override;
  foreach my $flag (@override) {
    foreach my $default (@defaults) {
      $flag = '' if ($flag eq $default);
    }
  }

# Turn back into a single string and tidy the formatting:
  $new_flags = join " ", @override;
  $new_flags =~ s/^\s*(.*?)\s*$/$1/;
  $new_flags =~ s/\s{2,}/ /g;

  return $new_flags;
}

1;
