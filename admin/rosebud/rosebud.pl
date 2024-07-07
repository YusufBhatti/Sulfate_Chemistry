#!/usr/bin/env perl
# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************
#
# Script:     rosebud.pl
# Purpose:    Convert UMUI processed output into a Rose suite
# Code owner: UM System Development team
#
# Creates:
#   app/<um_name>/rose-app.conf
#   app/<um_name>/file/<imports>
#   app/<fcm_make_um_name>/rose-app.conf
#   app/<fcm_make_um_name>/file/fcm-make.cfg
#   app/<fcm_make_ocean_name>/rose-app.conf [if coupled]
#   app/<fcm_make_ocean_name>/file/fcm-make.cfg [if coupled]
#   meta/rose-meta.conf
#   rose-suite.conf
#   rose-suite.info
#   suite.rc

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use Cwd;

use lib "$ENV{UMDIR}/rosebud";
use fileio;
use make_app;
use make_meta;
use make_config;
use make_ocean_config;
use make_suite;

# Main data hash to be passed to subroutines:
my %data;

# Generic looper variables
my $item;
my $line;

my $compute;      # Temporary compute host store
my $model = 'atmos'; # Model type (atmos, scm, coupled, etc.)

# Command-line arguments (including hash entries):
$data{access}       = "*";  # Access-list (commit permissions) for rose-suite.info
$data{defaults}     = "";   # Print variables that are set to default values?
$data{host2}        = "";        # Compute host
$data{host1}        = '$ROSE_ORIG_HOST'; # Extract host
$data{owner}        = `whoami`;  # Suite owner for rose-suite.info
$data{issues}       = "";        # Issue-list for rose-suite.info
$data{meta}         = "";   # Meta-data path for applying macros
# $data{suitedir}           # Output Rose suite directory
# $data{subproject}         # Sub-project for rose-suite.info
# $data{ibmhost}            # IBM compute host
# $data{jobdir}             # Input UMUI job directory
# $data{linuxhost}          # Linux compute host
# $data{appendrunid}        # Flag to add $RUNID to suite title
# $data{title};             # Suite title for rose-suite.info
my $name            = "";   # A subscript name/ID for distinguishing apps
my $verbose;                # Increase verbosity for debugging
my $spoiler;                # Enable spoilers
my $help;                   # Print help
my $man;                    # Go to man page

# Output variables:
my $ifupgrade = '';         # Modifier text if upgrading
my $macrologs = 'fixer';    # List of applied macros with logged output

# Other notable hash entries:
# version  - UM version (x.y)
# runid    - UMUI job $RUNID

chomp $data{owner};

# List of possible platforms to run on.
# Used for matching possible platforms in arguments to '-c'.
# This should match valid values of PLATFORM in the UMUI.
my @platforms = qw(ibm linux);

# Set up app and task names:
$data{makeapp}   = "fcm_make"; # Name of the fcm_make app
$data{runapp}    = "um";       # Name of the atmos/recon app
$data{atmostask} = "atmos";    # Name of the atmos task
$data{recontask} = "recon";    # Name of the recon task

# Process arguments
if (@ARGV > 0) {
  GetOptions('a|access=s'      => \$data{access},
             'd|defaults'      => \$data{defaults},
             'c|computehost=s' => \$data{host2},
             'e|extracthost=s' => \$data{host1},
             'i|issues=s'      => \$data{issues},
             'j|jobdir=s'      => \$data{jobdir},
             'n|name=s'        => \$name,
             'o|owner=s'       => \$data{owner},
             'p|subproject=s'  => \$data{subproject},
             'r|runid'         => \$data{appendrunid},
             's|suitedir=s'    => \$data{suitedir},
             't|title=s'       => \$data{title},
             'v|verbose'       => \$verbose,
             'F|fix'           => \$data{fixes},
             'M|meta=s'        => \$data{meta},
             'S|stash'         => \$data{stash},
             'P|profs'         => \$data{profs},
             'U|upgrade=s'     => \$data{upgrade},
             'spoiler'         => \$spoiler,
             'help|?'          => \$help,
             'man'             => \$man) or pod2usage(2);
}

pod2usage(1) if $help;
pod2usage(-exitval => 0, -verbose => 2) if $man;

# If there are any remaining arguments, abort with help
if (@ARGV > 0) {
  print STDERR "[FAIL] Unrecognised arguments: @ARGV\n";
  pod2usage(2);
}

# Add verbose flag to hash for use in subroutines
$data{verbose} = ($verbose) ? 1 : 0;

# Set input directory
my $pwd = cwd();
$data{jobdir} = $pwd unless ($data{jobdir});
print "[INFO] Reading job from: $data{jobdir}\n";
print "[INFO] It was his sled!\n" if ($spoiler);


# Sanity check: do we have some UMUI output that rosebud can process?
# 1: SCRIPT
my @script = fileio::read_file("$data{jobdir}/SCRIPT"); # Will abort here if no
                                                        # UMUI output to read
# Store file contents in data hash:
$data{script} = \@script;

# Pick up the title (unless provided) and UM version:
foreach $line (@script) {
  $data{version} = $1 if ($line =~ /export VN=(\d.\d)/);
  unless ($data{title}) {
    $data{title} = $1 if ($line =~ /export JOB_LINE='(.*)'$/);
  }
}


# Do we continue further?
if ($data{version} le '8.2') {
  die "[FAIL] rosebud is intended for use with jobs at UM 8.3 or later.\n"
    . "[FAIL] Please upgrade or choose a more recent job for Rose suite conversion.\n";
}

# Read miscellaneous (non-science) settings that may be needed often:

# 2: MAIN_SCR and EXTR_SCR
my @main_scr = fileio::read_file("$data{jobdir}/MAIN_SCR");
my @extr_scr = fileio::read_file("$data{jobdir}/EXTR_SCR");
# Store file contents in data hash:
$data{main_scr} = \@main_scr;
$data{extr_scr} = \@extr_scr;

foreach $line (@main_scr) {
  $data{runid} = $1 if ($line =~ /export RUNID=(\w+)/);
}
$data{suitedir} = "$pwd/$data{runid}-rose" unless ($data{suitedir});
print "[INFO] Writing suite to: $data{suitedir}\n";


# Check macro settings:
if ($data{meta}) {
  print "[INFO] Added to meta-data search path: $data{meta}\n";
# Warn against a common mistake:
  print "[WARN] Meta-data path addition does not end in 'rose-meta'\n" 
    unless ($data{meta} =~ m{rose-meta(/)?$});

# We only ever use this with the -M flag, so:
  $data{meta} = "-M $data{meta}";
} 


if ($data{upgrade}) {
  if ($data{version} eq '8.6') {

    # Running an upgrade requires us to fix the job first:
    $data{fixes} = 'true' unless ($data{fixes});

    # Amend output text for later:
    $ifupgrade = ' before upgrade';
    $macrologs .= '/upgrade';

  } else {

    # Cannot upgrade this UM version.
    print "[WARN] Upgrade macro can only be applied to UM 8.6 jobs. Macro disabled.\n";
    $data{upgrade} = '';

  }
} # End if (upgrade)


if ($data{fixes}) {
  if ($data{version} eq '8.6') {

  # Location of a file to store macro output:
    $data{macrofile} = "$data{suitedir}\.macros";
    unlink ($data{macrofile}) if (-e $data{macrofile});

  } else {
    print "[WARN] Fixer macro can only be applied to UM 8.6 jobs. Macro disabled.\n";
    $data{fixes} = '';
  }
}



# Set up app and task names:
if ($name) {
   $data{makeapp}   = $data{makeapp}   . "_" . $name;
   $data{runapp}    = $data{runapp}    . "_" . $name;
   $data{atmostask} = $data{atmostask} . "_" . $name;
   $data{recontask} = $data{recontask} . "_" . $name;
}

# Set some basic information: which platform, which job steps?
# - Basic atmos/recon tasks:
$data{platform}    = "";
$data{build_recon} = "";
$data{build_atmos} = "";
$data{build_ocean} = "";  # Separate atmos and ocean builds...
$data{run_recon}   = "";
$data{run_atmos}   = "";
$data{run_coupled} = "";  # ...but a single, joint run flag.

# Logical switches for other models
$data{scm} = "";
$data{coupled} = "";

# No. of settings that will break generic suite conversion:
$data{unportable} = 0;

# 3: SUBMIT
# If we have got this far we have a UMUI job of some kind. However if it's an 
# SCS component there won't be a SUBMIT file, so we pre-empt that failure
# with a helpful message.
my @submit;
if (-e "$data{jobdir}/SUBMIT") {
  @submit = fileio::read_file("$data{jobdir}/SUBMIT");
} else {
  die "[FAIL] Cannot find file $data{jobdir}/SUBMIT\n"
    . "[FAIL] If this job is an SCS component please disable External control and\n"
    . "       edit any fields set by the SCS to entries valid for standalone use.\n";
}

# Store file contents in data hash:
$data{submit} = \@submit;
foreach $line (@submit) {
  if ($line =~ /TARGET_MC=[^\$]/) {
    $data{platform}  = "ibm"   if ($line =~ /ibm/);
    $data{platform}  = "linux" if ($line =~ /linux/);
  }
    $data{run_atmos} = "true"  if ($line =~ /RUN_ATM=true/);
    $data{run_recon} = "true"  if ($line =~ /RUN_RCF=true/);
    $data{run_coupled} = "true" if ($line =~ /RUN_UM_NEMO_CICE=true/);
    $compute = $1 if ($line =~ /RHOST_NAME=(.*)/);
    $data{scm} = "true" if ($line =~ /SCM_ATMOS=true/);
}
$compute =~ s/(^\s*)(.*?)(\s*$)/$2/;  # strip whitespace

# SCM support was only added at UM 8.4:
if ($data{scm} and $data{version} eq '8.3') {
    die "[FAIL] The UM does not support SCM jobs in Rose at UM 8.3.\n"
      . "[FAIL] Please upgrade your job to UM 8.4 or later.\n";
}  

# Print extract host here (so it goes with the compute host):
print "[INFO] Overriding extract host with: $data{host1}\n"
  unless ($data{host1} eq '$ROSE_ORIG_HOST');


# Now we know the platform, we can parse the compute host argument:
# Requires a comma-separated list
if ($data{host2}) {

# If $hosts2 contains multiple names, we need to parse them:
# This assumes hostnames are 'sensible' (but env vars are allowed):
  my @hosts;
  push @hosts, $1 while ($data{host2} =~ s/(\w+:\$?\w+)//);

# Now split each key:value pair for hosts, if there were any:
  if (@hosts) {

    foreach my $pair (@hosts) {
      $pair =~ /(\w+):(\$?\w+)/;
      my $key   = $1;
      my $value = $2;

      my $known;      # Flag for detecting unknown platforms (possible typos)
      foreach $item (@platforms) {   # List of all possible platforms
        if ($key eq $item) {
          $known = 1;                # Matched a valid platform...

# ...but is it the particular $platform we read from the UMUI?
          $data{host2} = $value if ($key eq $data{platform});
        }
      }
      print "[WARN] Compute host supplied for unknown platform: $key\n"
        unless ($known);
    }

  } else {
# Do nothing; we assume a single value, read verbatim.
  }
}   # endif ($data{host2})


# Now either print out the overriding value, if we still have one, 
# or begin processing the compute host read from the UMUI:
if ($data{host2}) {

  print "[INFO] Overriding compute host with: $data{host2}\n";

} else {

# Evaluate any environment variables:
  if ($compute =~ /(\$\w+)/) {
    if ($compute =~ /^\$HOSTNAME$/) {

#   If the job should be run on the current machine, preserve that behaviour:
      print "[INFO] Replaced $compute with \$ROSE_ORIG_HOST as the compute host\n";
      $compute = '$ROSE_ORIG_HOST';

    } else {

#   Otherwise, evaluate whatever was supplied:
      $compute =~ s/\$(.*)/$1/;
      $compute = $ENV{$compute};
      print "[INFO] Evaluated $1 to $compute as the compute host\n";
    }
  }  # End 'if evaluate/substitute env vars'

  if ($data{platform} eq 'linux') {

#   Special Linux settings:
#     - if the existing machine is a server, set it to a random server.
#     - if the existing machine is a desktop, set it to $ROSE_ORIG_HOST.
#   The latter setting is to avoid putting unnecessary loads on the servers.
    if ($compute =~ /^[cerw]ls\d+$/) {
      $compute = 'random';
      print "[INFO] Compute host is a Linux server, resetting to random\n";

    } elsif ($compute =~ /^[cerw]ld\d+$/) {
      $compute = '$ROSE_ORIG_HOST';
      print "[INFO] Compute host is a Linux desktop, resetting to \$ROSE_ORIG_HOST\n";
    }

  }

# Assign the UMUI-defined compute host:
  $data{host2} = $compute;

}  # End 'unless override compute host'


# Brief sanity check.
# Sufficiently important as to be printed even if it's not being used.
print "[WARN] Failed to determine a compute host from UMUI job settings\n"
  unless ($compute);

# After all that, do we have a compute host?  
die "[FAIL] No suitable compute host specified" unless ($data{host2});


# 4: COMP_SWITCHES
my @switches = fileio::read_file("$data{jobdir}/COMP_SWITCHES");
foreach my $line (@switches) {
  $data{build_atmos} = "true" if ($line =~ /COMP_UMATMOS=true/);
  $data{build_recon} = "true" if ($line =~ /COMP_UMRECON=true/);
  $data{build_ocean} = "true" if ($line =~ /COMP_UM_NEMO_CICE=true/);
}

# Set coupling flag if either build or run steps for coupled model are true.
$data{coupled} = "true" if ($data{build_ocean} or $data{run_coupled});

# Coupled support was only added at UM 8.5:
if ($data{coupled} and $data{version} le '8.4') {
    die "[FAIL] The UM does not support Coupled jobs in Rose at UM 8.4 or earlier.\n"
      . "[FAIL] Please upgrade your job to UM 8.5 or later.\n";
}  


# For non-atmos suites:
#  - manually override any conflicting UMUI settings
#  - rename main task and apps to reflect the change of model
if ($data{scm}) {
  $data{build_recon} = ""; # Often still set to true but ignored by the UMUI
  $data{run_recon}   = ""; # Should already be false, but make certain
  $data{atmostask} =~ s/atmos/scm/;   # Rename main run task
  $data{runapp}    =~ s/um/scm/;      # Rename app
}

if ($data{coupled}) {
  $data{atmostask} =~ s/atmos/coupled/; # Rename main run task
  $data{runapp}    =~ s/um/coupled/;    # Rename app
# Now there are two fcm_make apps:
  $data{oceanmakeapp} = $data{makeapp};
  $data{oceanmakeapp} =~ s/fcm_make/fcm_make_ocean/;
  $data{makeapp}      =~ s/fcm_make/fcm_make_um/;
# Coupled sub-models version info:
# What if these aren't "vnx.y"? Can we infer the version from elsewhere?
  foreach $line (@extr_scr) {
    $data{nemoversion} = $2 if ($line =~ /export NEMO_VER=(vn)?(\w+)/);
    $data{ciceversion} = $2 if ($line =~ /export CICE_VER=(vn)?(\w+)/);
    $data{ioipslversion} = $2 if ($line =~ /export IOIPSL_VER=(vn)?(\w+)/);
  }
} else { 
# Default for non-coupled jobs
  $data{oceanmakeapp} = "";
}
  

# Now the basic platform/job settings should be available for use.

# Create suite directory structure:
fileio::create_dirs(\%data);
print "[INFO] Created Rose directory structure\n" if ($verbose);


print "[INFO] Writing individual suite files:\n" if ($verbose);

# Construct the rose-suite.conf file
$data{host2} = make_suite::write_conf(\%data);
print "[INFO] suite.conf file written\n" if ($verbose);


# Construct the app directory:
$data{envhashref} = make_app::write_app(\%data);
print "[INFO] app.conf file written\n" if ($verbose);


# Construct the meta directory:
make_meta::write_meta(\%data);
print "[INFO] meta.conf file written\n" if ($verbose);


# Contruct the UM fcm_make config:
($data{openmp}, $data{um_base}, $data{unportable}) 
  = make_config::write_fcm_make(\%data);
# Name changes if two configs are present:
my $modifier = ($data{coupled}) ? '_um' : ''; 
print "[INFO] fcm_make$modifier.cfg written\n" if ($verbose);

if ($data{coupled}) {
# Construct the ocean fcm_make config
# (also writes library flags for coupling to the UM fcm_make config):
  $data{unportable} = make_ocean_config::write_fcm_make_ocean(\%data);
  print "[INFO] fcm_make_ocean.cfg written\n" if ($verbose);
}

# Construct the suite.info file:
make_suite::write_info(\%data);
print "[INFO] suite.info written\n" if ($verbose);

# Construct the suite.rc file
$data{unportable} = make_suite::write_rc(\%data);
print "[INFO] suite.rc written\n" if ($verbose);

print "[INFO] Finished writing suite.\n";

if ($data{fixes}) {
  print "[INFO] Fixer macro applied to all apps$ifupgrade.\n";
  print "[WARN] Application of the fixer macro may have changed results.\n";
  print "[INFO] List of $macrologs macro changes: $data{macrofile}\n";
}

if ($data{unportable}) {
  my $plural = ($data{unportable} > 1) ? 's' : '';
  print "[WARN] Suite contains $data{unportable} setting$plural that may break portability between users.\n";
}

unless ($data{version} eq '8.6') {
  my $nextversion = $data{version} + 0.1;
  print "[INFO] This suite cannot be upgraded. Upgrade the UMUI job to get a vn$nextversion suite.\n";
}

print "[INFO] You can add this suite to the Rosie repository by running 'rosie create'\n"
    . "       and copying the suite files into the new directory created in ~/roses.\n";


__END__

=head1 NAME

rosebud

=head1 SYNOPSIS

rosebud.pl [-a access-list] [-c computehost] [-d] [-e extracthost]
           [-i issue-list] [-j jobdir] [-n name] [-o owner] 
           [-p subproject] [-r] [-s suitedir] [-t title] [-v]
           [-F] [-M meta-path] [-P] [-S] [-U version]

rosebud.pl -?|--help

rosebud.pl --man

=head1 DESCRIPTION

Convert a processed UMUI job into a Rose suite. 

B<rosebud> can be used with UM 8.3 -- UM 8.6 jobs.

The output suite can be run immediately, copied into a blank suite created 
with B<rosie create>, or integrated by hand into an existing suite.

Use the -j or -s options to run B<rosebud> on UMUI output in a directory where 
you do not have write permission.

=head1 OPTIONS

=head2 Program information

=over 8

=item B<-?, --help>

Print a message explaining the available command-line options, then exit.

=item B<--man>

Open the full manual page, including command-line options, examples and known 
caveats.

=back

=head2 Suite creation

=over 8

=item B<-a, --access> access-list

The value of "access-list" in rose-suite.info; a space-separated list of users 
with commit access to the trunk that will be generated once this Rose suite is 
lodged in the repository. Defaults to "*".

=item B<-c, --computehost> hostname

=item B<-c, --computehost> platform:hostname[,platform:hostname,...]

Name of the host that will perform the compile and/or run steps (COMPUTE_HOST).
Defaults to the remote host specified in the UMUI job. Note that in Linux 
suites it is the extract host that performs the compilation.

In the first form, the single hostname is used to override the compute host
for any UMUI job.

In the second form, a comma-separated list of platform:hostname pairs
is supplied. For each pair, the hostname will only be used to override the UMUI 
setting if the job is for that platform. At the Met Office, available platforms 
are: ibm, linux.

The special hostname "random" will use rose host-select to choose either a
least-loaded Linux server or a random HPC cluster, as appropriate. This  
facility should only be used when needed; please do not send small, 
desktop-capable Linux jobs to the compute servers by default.

=item B<-d, --defaults>

Supplying this option will cause the resulting rose-app.conf file to contain
any environment variables which are already set to their default values and
are not needed by the suite. By default these variables are not included.

=item B<-e, --extracthost> Linux-hostname

Name of the Linux host that will perform any extracts or mirroring of files
(EXTRACT_HOST). Defaults to $ROSE_ORIG_HOST, which is the name of the host that
launched Rose. Note that in Linux jobs the extract host also performs the 
compilation as a single make task.

Any branches or working copies in the job must be directly visible to this host,
e.g. a branch in /data/local on another machine will cause the extract to fail
when the job is run.

=item B<-i, --issues> issue-list

Space-separated list of tickets/issues associated with the Rose suite. 

=item B<-j, --jobdir> directory

Directory containing the UMUI job to be converted; more specifically, the 
directory containing the UMUI's processed output, typically 
$HOME/umui_jobs/$RUNID. This can be a relative or absolute path. Defaults to 
the current directory.

=item B<-n, --name> name

A name which will be appended to the apps and tasks generated, in the form
<app>_<name> and <task>_<name>. This can be used to give each app and task a 
unique name, which may be helpful if the generated files are intended to be 
merged into a larger suite.

=item B<-o, --owner> userid

The value of "owner" in rose-suite.info; the user ID of the owner of the suite.
The owner has full commit access and only the owner can transfer ownership
or delete the suite. Defaults to the current user.

=item B<-p, --subproject> sub-project

A subdivision of project which can be used to help identify the suite.
The 'sub-project' field is not included in rose-suite.info unless this option 
is provided. B<rosebud> always sets the 'project' field to 'um'.

=item B<-r, --runid>

Supplying this option will cause the RUNID of the UMUI job to be appended to 
the 'title' field in rose-suite.info. This may be useful in quickly tracing the 
pre-Rose history of a suite. The default behaviour is not to append this 
information.

The RUNID of the donor job will be printed in a separate 'description' field in
rose-suite.info whether or not this option is supplied.

=item B<-s, --suitedir> directory

The name of the directory in which the new Rose suite will be created. This can
be a relative or absolute path, and will be created if it does not exist. Any 
previous output from this script, or any other files of the same names, will be
overwritten. The default behaviour is to create the Rose suite in a directory 
"<RUNID>-rose" in the current directory.

=item B<-t, --title> text

The value of "title" in rose-suite.info; a brief description of the suite.
If not supplied this will default to the UMUI job description. This option can 
be used in conjunction with -r.

=item B<-v, --verbose>

Enable verbose output.

=back

=head2 Rose macros

=over 8

=item B<-F, --fix>

Apply the Rose fixer macro to all generated apps. This can only be used with
UM 8.6 jobs and will be applied using vn8.6 meta-data. There is no meta-data
available for fixing UM apps at earlier versions.

Running the fixer macro is a precondition for upgrading a suite to UM 9.0 or 
later, but may change results for a small minority of suites whose namelists 
are inconsistent with the meta-data.

=item B<-M, --meta> path

Prepend this path to the meta-data search path when applying Rose macros.

=item B<-P, --profs>

Apply the Rose macro "stash_requests.StashProfilesRemoveUnused" to the 
runtime app. The macro will be applied using vn8.6 meta-data. This will remove
unused profile namelists from the app and is recommended to reduce unnecessary
complication in the app.

=item B<-S, --stash>

Apply the Rose macro "stash_indices.TidyStashTransformPruneDuplicated" to the 
runtime app. The macro will be applied using vn8.6 meta-data. This macro is one
of two STASH transformer macros and is the one recommended for general use;
application of one of the two macros is necessary to use the Rose STASH GUI.

CAUTION: Application of the STASH transformer macros has been found to cause a
loss of bit comparison between different processor decompositions in some
suites.

=item B<-U, --upgrade> version

Apply the Rose upgrade macro to all generated apps. Apps will be upgraded to
the specified version. This can only be used with UM 8.6 jobs; there is no
meta-data available for upgrading UM apps at earlier versions.

Supplying this argument also implies B<--fix>. The fixer macro will be applied
both before and after each separate single-version upgrade to ensure meta-data
consistency at each stage. If the B<--stash> or B<--profs> arguments are also 
present the STASH or STASH profile transformer macros will be applied first. 

=back

=head1 EXAMPLES

=over 8

=item B<rosebud>

Create a Rose suite $RUNID-rose in the current directory, based on the UMUI
output in the current directory (where $RUNID is the UMUI's 5-letter 
identifier for that job).

=item B<rosebud -j ~fred/umui_jobs/uvxwz -s ../my_rose_suite>

Create a Rose suite in the directory '../my_rose_suite' based on the UMUI 
output in directory '~fred/umui_jobs/uvxwz'.

=item B<rosebud -e els049 -c hpc2f>

Create a Rose suite using els049 to perform the extract and hpc2f to compile 
and/or run the model.

=item B<rosebud -c ibm:hpc2f,linux:els002>

Create a Rose suite that will compile and/or run the model on hpc2f for an IBM 
job, or run the model on els002 for a Linux job.

=item B<rosebud -c linux:\$ROSE_ORIG_HOST>

Create a Rose suite that will run the model on $ROSE_ORIG_HOST (the host that 
launched Rose) for a Linux job, and use the existing UMUI setting for an IBM 
job. Note that the '$' is escaped to prevent the shell evaluating it.

=item B<rosebud -c random>

=item B<rosebud -c ibm:random,linux:random>

Create a Rose suite that will compile and/or run the model on a random HPC 
cluster for an IBM job, or run the model on the least-loaded Linux server for
a Linux job.

=item B<rosebud -t "My first suite" -r>

Create a Rose suite with the title "My first suite - derived from UMUI job 
$RUNID" in the resulting rose-suite.info file.

=item B<rosebud -n omp>

Create a Rose suite with the app names fcm_make_omp & um_omp and the task names
fcm_make_omp, recon_omp, etc.

=item B<rosebud -F -M ~fred/UM_working_copy/rose-meta>

Create a Rose suite and - if the donor job was at UM 8.6 - run the Rose fixer
macro on all apps using the vn8.6 meta-data contained in the indicated 
directory.

=item B<rosebud -S -U vn9.0>

Create a Rose suite and - if the donor job was at UM 8.6 - run the recommended
STASH transformer macro on the runtime app, then run the Rose upgrade macro to
upgrade all apps to vn9.0, using the default trunk meta-data in all cases.

=back

=head1 CAVEATS 

Users may wish to be aware of the following issues.

=over 8

=item B<General use>

Where possible, B<rosebud> produces generic suites that any user 
can run regardless of who owns the donor UMUI job. Users will be notified of
any job settings that may disrupt this process (e.g. working copies of source 
code, some LoadLeveler directives) so that they can find and address these 
settings if necessary.

If any steps (compile, reconfigure or run) are inactive in the UMUI job the 
resulting Rose suite will function but may lack information relating to the
missing tasks. It is easier to turn on all steps in the UMUI and then edit the 
resulting suite.rc graph to remove unwanted tasks, than it is to add in the 
missing tasks later.

When converting coupled jobs an ssh connection to the target machine is
required so that the NEMO and CICE namelist files can be copied and transferred
into the relevant rose-app.conf file.

=item B<Suitable jobs>

The SCM is supported from UM 8.4 onwards.

Limited coupled model support is available at UM 8.5: the only tested 
configurations at UM 8.5 are for jobs using NEMO 3.2 and OASIS3; those using 
NEMO 3.4 and/or OASIS3-MCT may require further UM branches or changes to the UM
configuration files. 

Full coupled model support for both NEMO 3.2 and NEMO 3.4, and for both the
OASIS3 and OASIS3-MCT couplers, is available at UM 8.6.

B<rosebud> does not support:

B<*> SCS component jobs; the UMUI output from such jobs does not contain enough
information to construct a Rose suite. Deactivate the external control switch, 
fill in any missing inputs, and try again.

B<*> CRUNs or jobs using automatic resubmission; these will be converted into a
suite that does not resubmit any tasks. The suite can then be modified to add
automatic resubmission.

B<*> serial or 32-bit builds; the UM's new Rose scripts do not yet support these
features.

B<*> Script inserts; whatever the insert does, the same effect should be 
achievable by scripting a solution within the suite itself.

=item B<ENDGame>

At B<UM 8.3> jobs converting between New Dynamics and ENDGame grids use two 
STASHmaster files, STASHMSTR_IN and STASHMSTR_OUT. To ensure the correct 
STASHmaster is used the rose-app.conf file must be edited to rename
one of these variables:

When converting B<to> ENDGame, rename STASHMSTR_OUT to STASHMSTR. 

When converting B<from> ENDGame, rename STASHMSTR_IN to STASHMSTR.

From UM 8.4 onwards no special treatment of ENDGame jobs is necessary.

=item B<Overrides>

Rose suites run from a new set of machine configuration files written in the
FCM2 syntax. B<rosebud> will attempt to convert any overrides in a UMUI job
into the new syntax, but you are strongly advised to check the resulting
fcm_make.cfg file before running the suite. In particular you should ensure any
file overrides are applied to both the atmosphere B<and> reconfiguration builds
if the affected files are used by both components, and that any overrides 
applied to an ocean build match those provided in the NEMO/CICE compiler
flags files.

At UM 8.5 the "debug" optimisation level settings in the new (FCM2) config 
files used at the Met Office differ from the old (FCM1) settings. File 
overrides using the special "%fflags_opt_<level>" variables will be converted 
automatically into the new format appropriate to that UM version to maintain
consistency with UMUI-based builds. At UM 8.6 the previous "debug" settings 
were restored and a new optimisation level, "rigorous", was added.

=back

=cut
