# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************
#
# Script:     make_app.pm
# Purpose:    Create rose-app.conf file and copy any other required input files
#             to the file/ subdirectory
# Code owner: UM System Development team
#
# A temporary space (e.g. in /var/tmp) is created and the namelist files are
# copied there for processing with rose namelist-dump and rose config-dump; 
# symlinks for the coupled model are also passed through a temporary file here.
# This allows rosebud to be run from a directory in which the user doesn't
# have write access.


package make_app;

use strict;
use warnings;
use fileio;
use macros;
use stashc_proc;
use Cwd qw/ cwd /;
use File::Temp qw/ tempdir /;

sub write_app {

# Write the rose-app.conf file and copy other inputs to file/ subdirectory.

  my $hashref = shift;

# Transfer relevant hash values to variables for ease of use
  my $verbose = $hashref->{verbose};           # True = verbose output
  my $jobdir = $hashref->{jobdir};             # job input directory
  my $suitedir = $hashref->{suitedir};         # suite output directory
  my $version = $hashref->{version};           # UM version x.y
  my $runid = $hashref->{runid};               # UMUI job ID
  my $host2 = $hashref->{host2};               # Compute host (ocean namelists)
  my $recontask = $hashref->{recontask};       # recon task name
  my $runapp = $hashref->{runapp};             # atmos/recon app name
  my $keepdefaults = $hashref->{defaults};     # True = print 'default' env vars
  my $build_atmos = $hashref->{build_atmos};   # True = build atmos exec
  my $build_recon = $hashref->{build_recon};   # True = build recon exec
  my $build_ocean = $hashref->{build_ocean};   # True = build ocean exec
  my $run_atmos = $hashref->{run_atmos};       # True = run atmos model
  my $run_recon = $hashref->{run_recon};       # True = run reconfiguration
  my $run_coupled = $hashref->{run_coupled};   # True = run coupled model
  my $oceanmakeapp = $hashref->{oceanmakeapp}; # Name of the ocean fcm_make app
  my $scm = $hashref->{scm};                   # True = SCM job
  my $coupled = $hashref->{coupled};           # True = Coupled model
  my $fixes = $hashref->{fixes};               # True = apply fixes macro
  my $stash = $hashref->{stash};               # True = apply STASH macro
  my $profs = $hashref->{profs};               # True = apply STASH profile macro
  my $upgrade = $hashref->{upgrade};           # Version to upgrade to
  my $macrofile = $hashref->{macrofile};       # File containing macro output
  my $meta = $hashref->{meta};                 # Path to meta-data
  my @script = @{$hashref->{script}};          # Contents of SCRIPT
  my @submit = @{$hashref->{submit}};          # Contents of SUBMIT


  my @app;              # app file contents
  my %envvars;          # environment variables subsection of app file

# flags, loop and other temporary variables
  my $copyscript;
  my $copyenv;
  my $file;
  my $item;          
  my $line;
  my $pair;
  my $variable;
  my $value;
  my $nemo_nl_file = "";  # Location of NEMO namelist file
  my $cice_nl_file = "";  # Location of CICE namelist file
  my @symlinks;           # Save a list of symlinks that the UMUI creates

# Lists of things that need copying into app:
  my @scriptvars = (
#                    'ANCIL_VERSIONS',      # Not accessed by name
                    'LOADMODULE',           # Path to atmos exec
                    'LOADRECON',            # Path to recon exec
                    'OCNMODULE',            # Path to ocean exec
                    'PRINT_STATUS',
                    'RCF_DEL_MPP_OUTPUT',   # Replaced by new var
                    'RCF_PRINTSTATUS',
                    'RCF_TIMER',
#                    'RUNID',               # Currently ignored
                    'STASHMSTR',            # Has a default value in SCRIPT
#                    'UM_ANCIL_FILENAMES',  # Not accessed by name
                    'UM_DEL_MPP_OUTPUT',    # Replaced by new var
                    'VN'
                    );
  my @submitvars = (
                    'FLUME_IOS_NPROC',
#                    'FLUME_STASH_NPROC', # Unsupported
#                    'OMP_NUM_THREADS',   # Task-specific; set in suite.rc
                    'NMPPE',             #  = UM_ATM_NPROCX
                    'NMPPN',             #  = UM_ATM_NPROCY
                    'NPROC_MAX',
                    'RCF_NPROCX',
                    'RCF_NPROCY',
                    'UM_ATM_NPROCX',
                    'UM_ATM_NPROCY',
                    'UM_THREAD_LEVEL',
                   );

  if ($coupled) {
# Additions for the coupled model:
    my @coupled_submitvars = (
                              'NEMO_NPROC', 'NEMO_IPROC', 'NEMO_JPROC',
                              'CICE_NPROC',
                              'PRISM_NPROC', 'CPL_TASK_SPACING',
                             );
    push @submitvars, @coupled_submitvars;
  }

# 0) Copy stuff that can't go into the config file:
  unless ($scm) {
    system "cp $jobdir/PRESM_A $suitedir/app/$runapp/file/user.PRESM_A";
    system "cp $jobdir/INITFILEENV $suitedir/app/$runapp/file/INITFILEENV";
  }
# Because coupled model info goes to work and not share, if we reconfigure we 
# have to edit the location of the start dump so it can be passed between tasks.
  if ($coupled and $run_recon) {
    my $envfile = "$suitedir/app/$runapp/file/INITFILEENV";
    my @initfileenv = fileio::read_file($envfile);
    foreach $line (@initfileenv) {
# Preserve only the filename of the start dump:
      $line =~ s!(export ASTART=)([^/]+/)+(.+)!$1\$ROSE_DATA/$3!;
    }
    fileio::write_file($envfile, @initfileenv);
  }


# 1) Basic header info
  push @app, "meta=um-atmos/vn8.6\n\n";
  push @app, "[command]\n";
  if ($scm) {
    push @app, "default=um-scm\n";
  } elsif ($coupled) {
    push @app, "default=um-coupled\n$recontask=um-recon\n";
  } else {
    push @app, "default=um-atmos\n$recontask=um-recon\n";
  }

# 2) Environment variables
# Start scanning the UMUI for environment variables that need exporting

# Loop through SCRIPT and SUBMIT:
  foreach $line (@script) {

# User-defined script inserts, UMUI-inserted variables, etc.
# Some variables with hard-wired names only needed for coupled models.
    $copyscript = 1 if ($line =~ /# Output Choices - Environment variables/
                    or  ($coupled and (   $line =~ /export OASIS=/
                                       or $line =~ /export NEMO_VERSION=/
                                       or $line =~ /export CICE_NL=/
                                       or $line =~ /export NAMCOUPLE_HOME=/))); 
    if ($copyscript) {

# Skip variables that are set to themselves (because the user picked the same
# name the UMUI uses):
      next if ($line =~ /^\s*export\s*(\w+)=\$\1\s*$/);

      if ($line =~ /^\s*export\s*(\w+)=(.*)/) {   # Strip 'export'
        my $variable = $1;
        $envvars{$variable} = $2;
        $envvars{$variable} =~ s/{|}//g;         # Strip curly braces
        $envvars{$variable} =~ s/\$?(\w+:-)+//g; # Strip tree of var1:-var2:-...
                                                 #  (Rose can't handle these)
      }
    }

# User-defined environment variables are always needed,
# some variables with hard-wired names only needed for coupled models
    $copyscript = 0 if ($line =~ /# User defined output directories/
                    or  ($coupled and (   $line =~ /export NEMO_CPL_TYPE=/
                                       or $line =~ /export NEMO_SEPARATE_MEANS=/
                                       or $line =~ /export CICE_KMT=/     
                                       or $line =~ /export RMP_DIR=/) ) );

# UMUI-generated variables in SCRIPT:
    foreach $item (@scriptvars) {
      if ($line =~ /$item=/) {
        if ($line =~ /.*$item=\$\{$item:-(.*)\}/) {
          $envvars{$item} = $1;
        } elsif ($line =~ /.*$item=(.*)/) {
          $envvars{$item} = $1;
        } else {
          print "[WARN] Cannot evaluate environment variable $item in line\n"
                . "$line\n";
        }
      }
    }

# SCRIPT vars needing special treatment:
# Ancil files:
    $envvars{UM_ANCIL_FILENAMES} = $1 if 
      ($line =~ /echo "ERROR: the Ancil filenames version (.*) not found"/);
# Ancil versions:
    $envvars{ANCIL_VERSIONS} = $1 if
      ($line =~ /echo "ERROR: the Ancil versions file (.*) not found"/);

# Save symlinks for the coupled model later:
    push @symlinks, $line if ($coupled and $line =~ /^\s*ln -s -f/);

  } # End loop through SCRIPT


# Process the vars to pick out anything in still curly braces, 
#  unresolved user variables, etc.
  %envvars = process_vars(\%envvars);


# On to SUBMIT:
  foreach $line (@submit) {
    foreach $item (@submitvars) {
      $envvars{$item} = $2 if ($line =~ /\s*(export)?\s*$item=(\w+)\s*#?/)
    }
  }

# Decide what to do with executable names: filename only or full path?

# Exec names (if not using any of the 'standard' names)
#   (No Rose support for serial reconfiguration yet)
  my $recon_exec = make_config::get_exec_name($jobdir, 'recon');
  my $atmos_exec = make_config::get_exec_name($jobdir, 'atmos');
  $envvars{RECON_EXEC} = $recon_exec if ($recon_exec);
  $envvars{ATMOS_EXEC} = $atmos_exec if ($atmos_exec);

# We may need to provide the full path to either exec:
  $envvars{RECON_EXEC} = locate_exec('recon', $build_recon, $run_recon, 
                                    $envvars{RECON_EXEC}, $envvars{LOADRECON});
  delete $envvars{LOADRECON};    # No more use for this.

# For a coupled model, the atmos task could be triggered two ways:
  my $run_model = 'true' if ($run_atmos or $run_coupled);
  $envvars{ATMOS_EXEC} = locate_exec('atmos', $build_atmos, $run_model, 
                                    $envvars{ATMOS_EXEC}, $envvars{LOADMODULE});
  delete $envvars{LOADMODULE};   # No more use for this.

  if ($coupled) {
    my $ocean_exec = make_config::get_exec_name($jobdir, 'ocean');
    $envvars{OCEAN_EXEC} = $ocean_exec if ($ocean_exec);
    $envvars{OCEAN_EXEC} = locate_exec('ocean', $build_ocean, $run_coupled,
                      $envvars{OCEAN_EXEC}, $envvars{OCNMODULE}, $oceanmakeapp);
    delete $envvars{OCNMODULE};  # No more use for this.
  }

# Any other loose ends:

  if ($coupled) {
# This flag must be inserted by hand (the old scripts use the equivalent 
# command line argument, which is no longer an option):
    $envvars{MP_PGMMODEL} = 'mpmd';

# Set COUPLER variable by checking number of OASIS processors.
# If PRISM_NPROC isn't present for some odd reason, assume an MCT job.
    if (exists $envvars{PRISM_NPROC} and $envvars{PRISM_NPROC} > 0) {
      $envvars{COUPLER} = 'OASIS3';
    } else {
      $envvars{COUPLER} = 'OASIS3-MCT';
    }
 
# DATAW is needed here because coupled models have less control over where
# output goes. DATAM is explicitly set to the same (which is the default
# behaviour in um-atmos).
    $envvars{DATAM} = '$CYLC_TASK_WORK_DIR';
    $envvars{DATAW} = '$CYLC_TASK_WORK_DIR';

# Save the ocean namelist file names for later then delete them so they
# don't break the model.
      $nemo_nl_file = $envvars{NEMO_NL} if (exists $envvars{NEMO_NL});
      $cice_nl_file = $envvars{CICE_NL} if (exists $envvars{CICE_NL});
      delete $envvars{NEMO_NL};
      delete $envvars{CICE_NL};

  } else { # Not a coupled job

# As of UM 9.1 these are near-mandatory due to INITFILEENV retirement
# - many of the items in nlcfiles will likely contain DATAW.
# We use the same default as in um-atmos.
    $envvars{DATAW} = '$ROSE_DATA';
    $envvars{DATAM} = '$ROSE_DATA';

  } # End if (coupled)

# Any fixes needed
  if ($version eq '8.3') {

# Export RECONTMP for some rcf jobs, e.g. AQUM
    my $newrunid = 'atmos';              # The default in um-recon
# Next line not used unless RUNID is read from SCRIPT.
    $newrunid = $envvars{RUNID} if (exists $envvars{RUNID});
    $envvars{RECONTMP} = "$newrunid.recontmp";

  } # End vn8.3 fixes

### ######################################## ###
### Post-processing of environment variables ###
### ######################################## ###

# Remove anything the SCM doesn't need:
  %envvars = remove_nonscm(\%envvars) if ($scm);

# Any conversions between old and new env var names:
  %envvars = convert_vars($version, \%envvars);

# Strip anything that's set to the default value, unless asked not to:
  %envvars = remove_defaults($version, \%envvars) unless ($keepdefaults);

# Variables are not added to the app until after we have created any symlinks 
# for the coupled model, so that any entries that are not needed directly can
# be removed from the list before it is added to the file.


# 3) File creation I: Setup

# Give each conversion a unique temporary space, according to 
# $TMPDIR:-/var/tmp:/tmp.
# We do this here so the any symlinks can also be written here and parsed
# when we run rose config on the combined list of links+namelists.
  my $tmpparent;
  $tmpparent = $ENV{TMPDIR} or $tmpparent = '';
  unless ($tmpparent) {
    print "[WARN] Cannot get value of \$TMPDIR from environment\n" if ($verbose);
  }

  unless (-d $tmpparent) {
    $tmpparent = "/var/tmp";
    unless (-d $tmpparent) {
        $tmpparent = "/tmp";
        unless (-d $tmpparent) {
          die "[FAIL] Unable to find a temporary directory for writing.\n"
            . "[FAIL] Set \$TMPDIR in your environment and try again.\n";
      }
    }
  }

  my $tmpdir = tempdir("$runid-rose.XXXXX", DIR => "$tmpparent");
  unless (-d $tmpdir) {
    mkdir "$tmpdir", 0755 or die "[FAIL] Cannot create directory $tmpdir: $!";
  }
  print "[INFO] Writing temporary suite files to $tmpdir\n" if ($verbose);
  my $tmpfile = 'rose-app.tmp'; # Temporary file for constructing the app


# 4) File creation II: Symlinks (Coupled models only):

# OASIS and NEMO require symlinks to several input files
  if ($coupled) {

    my @link;         # A single symlink that can be added to the app file
    my @links;        # All the symlinks to be added
    my $target;       # Symlink target
    my $name;         # Symlink name
    my %oasis_files;  # Grids files to create symlinks to

# We saved a list of symlinks while looping through SCRIPT earlier.
# The last two entries are hardwired and unnecessary (the UM scripts now create
# those links).
    pop @symlinks;
    pop @symlinks;

# Process the list
# Find the relevant lines, create a list of links to make, then add them:
    foreach $line (@symlinks) {
      if ($line =~ /ln -s -f\s*([^\s]+)\s*([^\s]+)/) {
        $target = $1;
        $name = $2;

# We only want the basename of each link name:
        $name =~ s!(([^/]+/)+)(.+)!$3!;

        @link = create_link($target, $name, $host2);
        push @links, @link;
      }
    } # End foreach line

# Some links were previously created by the scripts but are now in the app.
# A hash of files to insert: key=UMUI env var, value=link name
    %oasis_files = (NC_GRIDS => 'grids.nc',
                    NC_MASKS => 'masks.nc',
                    NC_AREAS => 'areas.nc',
                    NC_ANGLES => 'angles.nc',
                    NAMCOUPLE_DIR => 'cf_name_table.txt',
                   );

    while ( (my $key, $name) = each %oasis_files) {
      if (exists $envvars{$key} and $envvars{$key} ne "") {
# There is one special case where the name was appended to the variable later:
        if ($key eq 'NAMCOUPLE_DIR') {
          $target = "$envvars{$key}/cf_name_table.txt";
        } else {
          $target = $envvars{$key};
        }
        @link = create_link($target, $name, $host2);
        push @links, @link;
      }
    }

# Write the links to the temporary file:
    fileio::write_file("$tmpdir/$tmpfile", @links);

  } # End if (coupled)


# Now that we've finished using the environment variables we can remove 
# any intermediate entries and add what's left to the app file:
  %envvars = remove_redundant(\%envvars) if ($coupled);

  push @app, "\n[env]\n";
  foreach my $variable (sort keys %envvars) {
    push @app, "$variable=$envvars{$variable}\n";
  }
  push @app, "\n";


# 5) File creation III: Namelists

# Add the namelist information to the temporary file.

# Begin by creating a list of namelist files to add.
# Some files (e.g. RECONA, SCM_SET) may be missing depending on the job setup.
  my @nlfiles=('SHARED','SIZES','RECONA','INITHIS', 'SCM_SET');
  push @nlfiles, 'IOSCNTL' unless ($scm);

  my @filestocopy = (@nlfiles, 'CNTLALL', 'CNTLGEN', 'CNTLATM');
  push @filestocopy, 'STASHC' unless ($scm);

  foreach my $file (@filestocopy) {
    system "cp $jobdir/$file $tmpdir/$file" if (-e "$jobdir/$file");
  }
  my $pwd = cwd();
  chdir $tmpdir or die "[FAIL] Cannot chdir to $tmpdir: $!";
 

#  i) Generate a temporary file containing the processed namelists
  my $rosecommand1='rose namelist-dump -l';
  my $exit_status = 0;
  foreach my $file (@nlfiles) {
    if (-e $file) {
      $exit_status = system "$rosecommand1 $file >> $tmpfile";
      die "[FAIL] Cannot convert namelist into Rose application format"
          if ($exit_status);
    }
  }

  if ($coupled) {
# Coupled models require two extra namelist files on the target machine.
# These must be retrieved and may need a little extra work before Rose can 
# process them.

    foreach my $project (qw/ nemo cice /) {

      $file = ($project eq 'nemo') ? $nemo_nl_file :
              ($project eq 'cice') ? $cice_nl_file :
                                     die "[FAIL] Unknown Ocean sub-model.";

      my $ocean_file = process_ocean_namelist($project, $host2, $file, $tmpdir);
      if (-e $ocean_file) {
        $exit_status = system "$rosecommand1 $ocean_file >> $tmpfile";
        die "[FAIL] Cannot convert namelist into Rose application format"
            if ($exit_status);
        unlink $ocean_file;
      }  # End if (file exists)
    }  # End foreach (project)
  }  # End if (coupled)

  if ($scm) {
#   SCM only:
    foreach my $file (qw(CNTLALL CNTLATM)) {
      if (-e $file) {
        $exit_status = system "$rosecommand1 $file >> $tmpfile";
        die "[FAIL] Cannot convert namelist into Rose application format"
            if ($exit_status);
      }
    }
  } else {
#   Non-SCM jobs:
    system "cat SHARED SIZES CNTLALL CNTLGEN CNTLATM > NAMELIST";
    $exit_status = system "$rosecommand1 NAMELIST >> $tmpfile";
    die "[FAIL] Cannot convert namelist into Rose application format"
        if ($exit_status);
    unlink "NAMELIST";


#  ii) Convert STASHC to the new format for Rose
#      The following subroutine is essentially rose_stashc_proc3.pl
#      as taken from the UM @ r45681
    if ( stashc_proc::process_stash('./STASHC') ) {
      $exit_status = system "$rosecommand1 STASHC >> $tmpfile";
      die "[FAIL] Cannot convert namelist into Rose application format"
          if ($exit_status);
    } else {
      print "[WARN] Could not prepare STASHC file for Rose\n";
    }
  }  # End non-SCM jobs

#  iii) Reformat the temporary file and remove duplicated namelists 
   system "rose config-dump -qf $tmpfile";

#  iv) Other tidying up
   my @repeatednl = qw(streq domain time use items r2swclnl r2lwclnl upanca 
                       trans);
   my @namelists = fileio::read_file($tmpfile);
   foreach $line (@namelists) {

     $line =~ s/item=0,0,/item=/;
     foreach $item (@repeatednl) {
       if ($line =~ /[^\[]namelist:$item\(\d+\)[^\]]/) {
         $line =~ s/(namelist:$item\(\d+\) *)+/namelist:$item\(:\) /g;
         $line =~ s/ *$//; # Remove trailing whitespace added by previous line
       }
     }
   }
  push @app, @namelists;

# 6) Finish up: write the file.
  unlink "$tmpfile";
  unlink @nlfiles, 'CNTLALL', 'CNTLGEN', 'CNTLATM', 'STASHC', 'STASHC_preROSE';
  chdir $pwd or die "[FAIL] Cannot chdir to $pwd: $!";
  rmdir $tmpdir or print "[WARN] Cannot remove temporary directory $tmpdir: $!\n";
  fileio::write_file("$suitedir/app/$runapp/rose-app.conf",@app);

# Run the fixes macro if requested:
  macros::fixes("$suitedir/app/$runapp", $macrofile, $meta) if ($fixes);

# And the STASH macro, if requested:
  if ($stash) {
    macros::stash("$suitedir/app/$runapp", $meta);
    print "[INFO] Applied stash_indices.TidyStashTransformPruneDuplicated macro.\n";
    print "[WARN] The STASH macro has caused a loss of bit comparison in some suites.\n";
  }

# And the STASH profile macro, if requested:
  if ($profs) {
    macros::profs("$suitedir/app/$runapp", $meta);
    print "[INFO] Applied stash_requests.StashProfilesRemoveUnused macro.\n";
  }

# And the upgrade macro, if requested:
  macros::upgrade("$suitedir/app/$runapp", $macrofile, $meta, $upgrade)
    if ($upgrade);

# Return the environment variables for use with path overrides later.
  return \%envvars;
}


sub process_vars {

# Tidy up environment variables into a Rose-friendly form

  my $hashref = shift;
    my %envvars = %$hashref;   # Environment variables to be processed


# Clean up each entry.
# Done as a reverse sort so that if we have e.g. two variables VAR and VAR_1
# occurrences of the longer variable will be matched and substituted first.
  foreach my $variable (reverse sort keys %envvars) {
    $envvars{$variable} =~ s/{|}//g;   # Strip curly braces
    $envvars{$variable} =~ s/#.*//;    # Strip trailing comments
    $envvars{$variable} =~ s/\s*$//;   # Strip trailing whitespace

# Check if any of the other entries are using this variable, 
#  and substitute if necessary:
    foreach my $key (keys %envvars) {
      $envvars{$key} =~ s/\$$variable/$envvars{$variable}/;
    }

  }

  return %envvars;
}

sub locate_exec {

# Decide whether to provide the full path to an exec.
# Assume any paths containing $DATAW refer to a local build (e.g. 
# $DATAW/bin/<exec>) and should not be part of the end Rose suite; the user 
# probably forgot the build was deactivated.

  my $target = shift;    # atmos, recon or ocean
  my $build = shift;     # True = build this exec
  my $run = shift;       # True = run this exec
  my $exec = shift;      # exec name (basename)
  my $module = shift;    # Full path to exec
  my $app = shift;       # Optional: Name of the ocean make app

  my @globbed;

# The full path to the Ocean exec is always needed.
  if ($build) {

    $exec = "\$CYLC_SUITE_SHARE_DIR/$app/build-ocean/bin/$exec"
      if ($target eq 'ocean');

  } else {     # No build
# No Rose build.
# This is the only time the atmos and recon tasks need to worry about the path.

    if ($module =~ /\$DATAW/) {

      if ($run) {
        die "[FAIL] UMUI job has no $target build step but runs an executable in \$DATAW.\n"
          . "[FAIL] Enable $target compilation or provide a complete path to the $target exec.\n";
      } else {
#     Neither run nor build; for atmos and recon tasks the exec basename is 
#     returned ready for future builds.
#     For ocean builds we add the 'default' path for future builds
#     (although the UMUI may not allow coupled jobs that do this...?)
        $exec = "\$CYLC_SUITE_SHARE_DIR/$app/build-ocean/bin/$exec"
          if ($target eq 'ocean');

      } # End if (run)

    } else {   # No $DATAW

# A non-standard path, not containing $DATAW. 
# Preserve this path, regardless of run settings.
      @globbed = glob $module;
      $exec = $globbed[0];
    } # End if (contains DATAW)
  } # End if (build)

  return $exec;
}


sub remove_nonscm {

# Remove environment variables that the SCM doesn't use
  my $hashref = shift;
    my %envvars = %$hashref;   # Environment variables to be examined

  my @umvars = qw(UM_DEL_MPP_OUTPUT RCF_DEL_MPP_OUTPUT 
                  NEMO_NPROC
                  NMPPN NMPPE
                  RCF_NPROCX RCF_NPROCY
                  RCF_PRINTSTATUS RCF_TIMER);

  foreach my $umvar (@umvars) {
    delete $envvars{$umvar} if (exists $envvars{$umvar});
  }

  return %envvars;
}

sub convert_vars {

# Convert environment variables and/or their values from UMUI to Rose settings
# For consistency we use: atmos not um; recon not rcf.

  my $version = shift;
  my $hashref = shift;
    my %envvars = %$hashref;   # Environment variables to be converted

# Set up a hash containing the substitutions to be made:
  my %replacements = ( UM_DEL_MPP_OUTPUT  => 'ATMOS_KEEP_MPP_STDOUT', # reverse
                       RCF_DEL_MPP_OUTPUT => 'RECON_KEEP_MPP_STDOUT', # reverse
                       NMPPE              => 'UM_ATM_NPROCX',
                       NMPPN              => 'UM_ATM_NPROCY',
                       ATM_CPL_TYPE       => 'ATMOS_COUPLE_TYPE',    # coupled
                       NEMO_CPL_TYPE      => 'OCEAN_COUPLE_TYPE',    # coupled
                       CPL_TASK_SPACING   => 'COUPLED_TASK_SPACING', # coupled
                       NAMCOUPLE_HOME     => 'NAMCOUPLE_DIR',        # coupled
                       CICE_ATM_DATA      => 'CICE_ATMOS_DATA', # CICE-only?
                       CICE_OCN_DATA      => 'CICE_OCEAN_DATA', # CICE-only
                     );

# The UM executables may need renaming:
  my $atmos_exec;
  my $recon_exec;
  if ($version eq '8.3') {
    $atmos_exec = 'um.exe';
    $recon_exec = 'qxreconf'; # Simplifies logic; name doesn't change at UM 8.3
  } else {
    $atmos_exec = 'um-atmos.exe';
    $recon_exec = 'um-recon.exe';
  }


# Loop over each variable looking for things to do:
  foreach my $variable (keys %envvars) {

# Some vars have their values altered.
# Do this first so we can still easily access the ones that get renamed.
    if ($variable =~ /_DEL_MPP_OUTPUT/) {
      $envvars{$variable} = ($envvars{$variable} eq 'false') ? 'true' : 'false';

    } elsif ($variable eq 'ATMOS_EXEC') {
# - We assume anything containing $RUNID is a standard name (e.g. $RUNID or 
#   $RUNID.exe) and can be overridden; trial/operational renames are usually 
#   hardwired.
#   If not defined at all, give it the default value:
      $envvars{$variable} = $atmos_exec unless ($envvars{$variable});
      $envvars{$variable} = $atmos_exec if ($envvars{$variable} =~ /\$RUNID/i)

    } elsif ($variable eq 'RECON_EXEC') {
#   (No support for serial reconfiguration yet.)
#   If not defined at all, give it the default value:
      $envvars{$variable} = $recon_exec unless ($envvars{$variable});
      $envvars{$variable} = $recon_exec if ($envvars{$variable} eq 'qxreconf');

    } elsif ($variable eq 'OCEAN_EXEC') {
#   Will contain the full path to the exec.
      $envvars{$variable} =~ s/(model|nemo)\.exe$/nemo-cice.exe/; # NEMO 3.2|3.4

    } elsif ($variable eq 'NEMO_RESTART_DATE') {
      $envvars{$variable} = ($envvars{$variable} eq 'Y') ? 'true' : 'false';

    } # End tests for changing variable values

# Loop over list of replacement names:
    while ( (my $oldname, my $newname) = each %replacements) {

      if ($variable eq $oldname) {
        $envvars{$newname} = $envvars{$oldname};
        delete $envvars{$oldname};
      } # End if (replace variable)

    } # End while loop (%replacements)

  } # End foreach (%envvars)

  return %envvars;
}


sub remove_defaults {

# Remove any environment variables that are set to their default value

  my $version = shift;         # UM version x.y
  my $hashref = shift;
    my %envvars = %$hashref;   # Environment variables to be processed

  my %defaults = ();

# Defaults that vary by version:
  if ($version eq '8.3') {
    $defaults{ATMOS_EXEC} = 'um.exe';
    $defaults{RECON_EXEC} = 'qxreconf';
  } else {
    $defaults{ATMOS_EXEC} = 'um-atmos.exe';
    $defaults{RECON_EXEC} = 'um-recon.exe';
  }
  if ($version eq '8.6') {
    $defaults{DR_HOOK} = '0';
    $defaults{DR_HOOK_OPT} = 'noself';
  }
    
  while ( my ($variable, $value) = each %envvars ) {
    my $default = $defaults{$variable};
    delete $envvars{$variable} if 
      ((exists $defaults{$variable}) and ($value eq $default));
  }

  return %envvars;
}


sub remove_redundant {

# Remove intermediate and other now obsolete variables inserted by the UMUI:
# Only required by coupled models.

  my $hashref = shift;
    my %envvars = %$hashref;   # Environment variables to be processed

  my @redundant = qw/ NC_ANGLES NC_AREAS NC_GRIDS NC_MASKS NEMO_NL_ICE 
                      OASIS OASIS_MPI_TYPE USE_GRIDS_DIRECT /;

  foreach my $variable (@redundant) {
    delete $envvars{$variable};
  }
  return %envvars;
}

sub process_ocean_namelist {

# Tidy up a NEMO or CICE namelist file so it can be put through Rose.

  my $project = shift;  # nemo or cice
  my $host = shift;
  my $infile = shift;
  my $outdir = shift;   # Should be cwd, but here for portability

  my $outfile;          # Processed output file
  my @nlfile;           # Contents of namelist file

  if ($project eq 'nemo') {
    $outfile = 'namelist';
  } elsif ($project eq 'cice') {
    $outfile = 'ice_in';
  } else {
     die "[FAIL] Cannot create namelist file for unknown project."
  }

# Copy the file from the target machine:
  system "scp \$(rose host-select $host):$infile $outdir/$outfile 1>/dev/null";
  die "[FAIL] Could not copy $project namelist file from target machine $host\n"
    unless (-e "$outdir/$outfile");

# Clean up any unprotected strings.
# There's no point over-engineering this; the problem items can only be those
# that are modified by the UMUI/scripts and should be labelled as such.
  @nlfile = fileio::read_file("$outdir/$outfile");
  foreach my $line (@nlfile) {

# Add quotes to any strings missing them:
    if ($line =~ /set_by_umui/i) {
      $line =~ s/set_by_umui/set_by_um/gi;  # Since we're here anyway...
      $line =~ s/(set_by_um)/'$1'/i 
        unless ($line =~ /["']\s*set_by_um\s*["']/i);
    }

  }
# Replace the file with a clean version and return its name:
  fileio::write_file("$outdir/$outfile", @nlfile);

  return $outfile;
}


sub create_link {

# Create a symlink for insertion into the app file.

  my $target = shift;    # Link target
  my $name = shift;      # Link name
  my $host = shift;      # Target machine

  my @link;              # Output for the app file

# Rose can't deal with paths beginning with ~, so expand these:
  if ($target =~ /^~/) {
    $target=`ssh \$(rose host-select $host) "echo $target"`;
    chomp $target;
  }

  push @link, "[file:$name]\n";
  push @link, "mode=symlink\n";
  push @link, "source=$target\n\n";

  return @link;
}

1;
