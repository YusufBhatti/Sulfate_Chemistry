# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************
#
# Script:     make_suite.pm
# Purpose:    Write the top-level suite files: 
#               suite.rc, rose-suite.conf, rose-suite.info
# Code owner: UM System Development team

package make_suite;

use strict;
use warnings;
use fileio;


sub write_conf {

# Generate the rose-suite.conf file

  my $hashref = shift;

# Transfer relevant hash values to variables for ease of use
  my $suitedir = $hashref->{suitedir};  # Suite output directory
  my $platform = $hashref->{platform};  # IBM or Linux (or ...)?
  my $host1   = $hashref->{host1};      # Extract host
  my $host2   = $hashref->{host2};      # Compute host

  my @conf    = "";  # File contents for writing

# Platform specific choices:
  if ($platform eq 'ibm') {

    if ($host2 =~ /^random$/i) {

# Choose a random hpc and login node
      $host2 = 'hpc';

    } else {

# choose a random login node on the desired machine
      $host2 =~ s/(\w+)-l\d+/$1/;   # Remove any explicit choice of login node
    }

  } elsif ($platform eq 'linux') {

# Linux: pick a random server if:
#  - the user specifically requests it, or
#  - the UMUI job was already set to use a server;
    $host2 = 'linux' if ($host2 =~ /^random$/i);

#  Otherwise do nothing, because:
#  - the user has selected a specific machine, or
#  - this has already been set to $ROSE_ORIG_HOST.
#  There is no need to pile work on the servers that would run on a desktop.
  }

# Write the file:
  push @conf, "[jinja2:suite.rc]\n";
  push @conf, "COMPUTE_HOST='$host2'\n";
  push @conf, "EXTRACT_HOST='$host1'\n";

  fileio::write_file("$suitedir/rose-suite.conf", @conf);

# Return the modified value to the parent hash for any future use.
  return $host2;
}



sub write_info {

# Generate the rose-suite.info file

  my $hashref = shift;

# Transfer relevant hash values to variables for ease of use
  my $suitedir = $hashref->{suitedir};  # Suite output directory
#   Job-configurable settings
  my $title = $hashref->{title};        # Suite title
  my $appendrunid = $hashref->{appendrunid};  # True = append $RUNID to suite title
  my $runid = $hashref->{runid};        # Previous value of $RUNID
  my $owner = $hashref->{owner};        # Suite owner
  my $access = $hashref->{access};      # Access-list (commit permissions)
  my $issues = $hashref->{issues};      # Issue-list
  my $subproject = $hashref->{subproject};    # Sub-project (if provided)


  my @info = "";          # Contents of file for writing
  my $project = "um";     # Project; fixed until coupled model support is added


#  Now write the file out, line-by-line
  push @info, "#-------------------------------------------------------------------------------\n";
  push @info, "# SUITE DISCOVERY INFORMATION FILE\n";
  push @info, "#-------------------------------------------------------------------------------\n";
  push @info, "# The \"owner\", \"project\" and \"title\" fields are compulsory. Any KEY=VALUE pairs\n";
  push @info, "# can be added. The following are known to have special meanings:\n";
  push @info, "#\n";
  push @info, "# description: A long description of the suite. E.g.:\n";
  push @info, "# description=This is line 1 of the long description of the suite.\n";
  push @info, "#             This is line 2 of the long description of the suite.\n";
  push @info, "#             This is line 3 of the long description of the suite.\n";
  push @info, "#\n";
  push @info, "# sub-project: A subdivision of \"project\", if relevant.\n";
  push @info, "#-------------------------------------------------------------------------------\n";
  push @info, "\n";
  push @info, "# access-list: A space-separated list of users with commit access to the trunk.\n";
  push @info, "# For example:\n";
  push @info, "# access-list=jane bob fred\n";
  push @info, "access-list=$access\n";
  push @info, "# description: A long description of the suite.\n";
  push @info, "description=Generated from UMUI job $runid using the script:\n";
  push @info, "            $0\n";
  push @info, "# issue-list: A space-separated list of tickets/issues relevant to the suite.\n";
  push @info, "# For example:\n";
  push @info, "# issue-list=UM:#4300 VAR:#1211\n";
  push @info, "issue-list=$issues\n";
  push @info, "# User ID of the owner of the suite.\n";
  push @info, "# The owner has full commit access to the suite.\n";
  push @info, "# Only the owner can pass the suite's ownership to someone else.\n";
  push @info, "# Only the owner can delete the suite.\n";
  push @info, "owner=$owner\n";
  push @info, "# The name of a project associated with the suite.\n";
  push @info, "project=$project\n";
  if ($subproject) {
    push @info, "# sub-project: A subdivision of \"project\".\n";
    push @info, "sub-project=$subproject\n";
  }
  push @info, "# A short title for the suite.\n";
  push @info, "title=$title";

  push @info, " - derived from UMUI job $runid" if ($appendrunid);
  push @info, "\n";


# And write it to file:
  fileio::write_file("$suitedir/rose-suite.info", @info);
}


sub write_rc {

# Generate the suite.rc file.

# The basic structure of each suite.rc file can be more or less the same
#  regardless of what combination of build or run steps it contains, because 
#  only the bits specified by the graph will actually be executed. Therefore 
#  - with the exception of suites for different models, e.g. the SCM - 
#  we do not bother checking each separate step to see if it's needed - and so
#  some files may come out looking rather odd if information is missing, though 
#  they should work fine.

  my $hashref = shift;

# Transfer relevant hash values to variables for ease of use
  my $jobdir  = $hashref->{jobdir};          # Job input directory
  my $suitedir = $hashref->{suitedir};       # Suite output directory
  my $platform = $hashref->{platform};       # IBM or Linux (or ...)?
  my $version = $hashref->{version};         # UM version x.y
  my $build_atmos = $hashref->{build_atmos}; # True = atmos build required
  my $build_recon = $hashref->{build_recon}; # True = recon build required
  my $build_ocean = $hashref->{build_ocean}; # True = ocean build required
  my $run_atmos = $hashref->{run_atmos};     # True = run the atmos model
  my $run_recon = $hashref->{run_recon};     # True = run the reconfiguration
  my $run_coupled = $hashref->{run_coupled}; # True = run the coupled model
  my $recontask = $hashref->{recontask};     # Name of the recon task
  my $atmostask = $hashref->{atmostask};     # Name of the atmos or coupled task
  my $makeapp   = $hashref->{makeapp};       # Name of the fcm_make app
  my $oceanmakeapp = $hashref->{oceanmakeapp};  # Name of the ocean fcm_make app
  my $runapp    = $hashref->{runapp};        # Name of the atmos/recon app
  my $openmp = $hashref->{openmp};           # True = OpenMP on
  my $unportable = $hashref->{unportable};   # No. of unportable settings
  my $scm = $hashref->{scm};                 # True = SCM job
  my $coupled = $hashref->{coupled};         # True = Coupled model
  my @submit = @{$hashref->{submit}};        # Contents of SUBMIT

# Local variables:
  my @rc = "";                # suite.rc file contents
  my @graph = "";             # Suite graph
  my @llheader = "";          # LoadLeveler header for suite.rc insertion
  my $llheader_ref;           # A reference for the LoadLeveler header
  my $warnings;               # No. of warnings of unportable settings 
                              # returned from get_llheader
  my $local_family;           # Family name for inheritance of local tasks
  my $remote_family;          # Fmaily name for inheritance of remote tasks


# These simply have to share a name with the appropriate app:
  my $maketask  = $makeapp;            # Name of the fcm_make task
  my $make2task = $maketask;       
  $make2task =~ s/fcm_make/fcm_make2/; # Name of the fcm_make2 task

  my $oceanmaketask = "";
  my $oceanmake2task = "";
  if ($coupled) {
    $oceanmaketask  = $oceanmakeapp;  # Name of ocean fcm_make task
    $oceanmake2task = $oceanmaketask;
    $oceanmake2task =~ s/fcm_make/fcm_make2/; # Name of ocean fcm_make2 task
  }

  my $defaultclass = 'parallel';    # Serial or parallel
  my $group = "";                   # IBM account group override (optional)
  my $graphsize = 1;                # Number of (active) lines in the graph


# Collect up any information we will need later (setting safe defaults first):
  my $progenv;
  if ($platform eq "ibm") {
    $progenv = ($version eq '8.6') ? "prg_13_1_0_8a" : "prg_12_1_0_10";
  } elsif ($platform eq "linux") {
    $progenv = "prg_ifort-12.0";
  }
  my $compilejobs = 5;
  my $recon_stack = 1000000;
  my $atmos_stack = 4000000;
  my $omp_num_threads = 2;
  my $omp_stacksize   = "1g";

  my @fcm_build = fileio::read_file("$jobdir/FCM_BLD_COMMAND");

# Various UMUI variables we may need:

  foreach my $line (@fcm_build) {
# The easiest way to detect the programming environments is that they both 
#  source something with 'prg_' in it. Apologies to whoever this eventually
#  breaks for.
    $progenv     = $1 if ($line =~ /\. ((.*?)prg_.*)/);
    $compilejobs = $1 if ($line =~ /BLD_CMD="fcm build -v \d -j (\d) "/);
  }
# Add a (temporary?) fix for the module symbol file errors found on the HPC:
  $progenv .= "; module load aix/libc-compile-fix" if ($platform eq "ibm");

  if ($platform eq "linux") {
    foreach my $line (@submit) {
      $recon_stack = $1*1000000 if ($line =~ /^RCF_STACK=(.*)/);
      $atmos_stack = $1*1000000 if ($line =~ /^STACK=(.*)/);
      $omp_num_threads = $1     if ($line =~ /export OMP_NUM_THREADS=(\d+)/);
      $omp_stacksize   = $1     if ($line =~ /export OMP_STACKSIZE=((\d+)(\w?))/);
    } 
  } elsif ($platform eq "ibm") {
    foreach my $line (@submit) {
#     Group override: if selected, the UMUI inserts this in all components 
#     regardless of compile & run settings. (Non-portable; generates a warning.)
      $group = $1 if ($line =~ /#@\s*group\s*=\s*(\w+)/);
    }
  }


# Basic header stuff:
  push @rc, "#!jinja2\n";
  push @rc, "[cylc]\n";
  push @rc, "    UTC mode = True\n";
  push @rc, "\n";
  push @rc, "[scheduling]\n";
  push @rc, "    [[dependencies]]\n";


# Build graph according to compile/run settings.
# Assume a 2-step build is the default for most platforms.
#
# NOTE: if the user elects to build only the reconfiguration then run only the 
# coupled model, there is an unwanted dependency between the build and the 
# model. Given how long the reconfiguration takes to build, and the likelihood 
# of this combination arising at all, removing this dependency is low priority.
#
# Is the graph going to need a second line?
  if ($build_ocean and ($build_atmos or $build_recon or $run_recon)) {
    $graphsize = 2;
  }
  push @graph, "        graph = \"";

  push @graph, "\"\"\n                " if ($graphsize > 1);
  if ($build_atmos or $build_recon) {
    push @graph, "$maketask => ";
    push @graph, "$make2task => " if ($platform ne "linux");
  }
  push @graph, "$recontask => " if ($run_recon);

# The following task can refer to both atmosphere-only and coupled tasks:
  push @graph, "$atmostask" if ($run_atmos);
  push @graph, "$atmostask" 
    if ($run_coupled and ($build_atmos or $run_recon or not $build_ocean));
# The final 'not' accounts for the case of a one-line graph with no ocean build.

  $graph[-1] =~ s/ => $//;   # strip any odd ' => ' left over

# Coupled model additions.
# Note: a second line is only triggered if there's an ocean build, and the 
# run-only case is dealt with above. So this block is protected by
# if (build_ocean) rather than if (coupled) to simplify subsequent logic.
  if ($build_ocean) {
    push @graph, "\n                " if ($graphsize > 1);
    push @graph, "$oceanmaketask => ";
    push @graph, "$oceanmake2task => " if ($platform ne "linux");
    push @graph, "$atmostask" if ($run_coupled);
    $graph[-1] =~ s/ => $//;   # strip any odd ' => ' left over
    push @graph, "\n                \"\"" if ($graphsize > 1);
  } # End coupled model additions

  push @graph, "\"\n";       # add trailing quote, newline

  push @rc, @graph;

# More common-to-all-jobs Rose stuff (<= technical name):
  push @rc, "\n";
  push @rc, "[runtime]\n";
  push @rc, "    [[root]]\n";
  push @rc, "        init-script = \"\"\"\n";
  push @rc, "export CYLC_VERSION={{CYLC_VERSION}}\n";
  push @rc, "export ROSE_VERSION={{ROSE_VERSION}}\n";
  push @rc, "\"\"\"\n";
  push @rc, "        script = \"rose task-run --verbose\"\n";
  push @rc, "        [[[events]]]\n";
  push @rc, "            mail events = submission failed, submission timeout, submission retry, retry, failed, timeout\n";
  push @rc, "            submission timeout = PT12H  # 12 hours\n";
  push @rc, "            execution timeout  = PT3H   # 3 hours\n";

# The rest of the file will be highly platform-dependent.
# UM-specific instructions go in the families for easy merger with other suites.
# Prompted by a user finding that the pre-command scripting in [[root]] was
# incompatible with 'moo' commands in archiving tasks.

  $local_family = 'LINUX_UM';
  $remote_family = ($platform eq 'linux') ? $local_family : 'HPC_UM';

# Only HPC platforms require different local and remote families:
  if ($platform ne 'linux') {
    push @rc, "    [[$local_family]]\n";
    push @rc, "        [[[job]]]\n";
    push @rc, "            batch system = background\n";
    push @rc, "        [[[remote]]]\n";
    push @rc, "            host = \$(rose host-select {{ EXTRACT_HOST }})\n";
  }

  push @rc, "    [[$remote_family]]\n";
  push @rc, "        pre-script = \". $progenv\"\n";

  push @rc, "        [[[job]]]\n";
  if ($platform eq "ibm") {
    push @rc, "            batch system = loadleveler\n";
  } elsif ($platform eq "linux") {
    push @rc, "            batch system = background\n";
  }

  push @rc, "        [[[remote]]]\n";
  push @rc, "            host = \$(rose host-select {{ COMPUTE_HOST }})\n";

# Some defaults
  if ($platform eq "ibm") {
    push @rc, "        [[[directives]]]\n";
    push @rc, "            class        = $defaultclass\n";
    if ($group) {
      push @rc, "            group        = $group\n";
      print "[INFO] LoadLeveler directive 'group' set to $group for all tasks.\n";
      $unportable++;
    }
    push @rc, "            job_type     = $defaultclass\n";
    push @rc, "            notification = never\n";
    
  } elsif ($platform eq "linux") {
# The MPI dir should really be deduced from the job, but in practice it's fixed,
#  at least for the UM versions we are considering.
    push @rc, "        [[[environment]]]\n";
    push @rc, "            PATH = /home/h01/frum/mpi/mpich2-1.4.1/ifort-12.0/bin:\$PATH\n";
  }


##### UM FCM_MAKE STEP #####

  push @rc, "    [[$maketask]]\n";
  push @rc, "        inherit = $local_family\n";

  if ($platform eq "linux") {
    push @rc, "        [[[remote]]]\n";
    push @rc, "            host = \$(rose host-select {{ EXTRACT_HOST }})\n";
    push @rc, "        [[[environment]]]\n";
    push @rc, "            ROSE_TASK_N_JOBS  = $compilejobs\n";
  }


##### UM FCM_MAKE2 STEP #####     ---   No mirror on Linux.

  if ($platform ne "linux") {
    push @rc, "    [[$make2task]]\n";
    push @rc, "        inherit = $remote_family\n";
    push @rc, "        [[[environment]]]\n";
    push @rc, "            ROSE_TASK_N_JOBS  = $compilejobs\n";

    if ($platform eq "ibm") {
      ($llheader_ref, $warnings) 
        = get_llheader(\@submit, $defaultclass, $make2task);
      $unportable += $warnings;
      push @rc, @$llheader_ref;
    }
    
  }   # End make2 step


  if ($coupled) {
##### OCEAN FCM_MAKE STEP ##### --- Only for coupled models

    push @rc, "    [[$oceanmaketask]]\n";
    push @rc, "        inherit = $local_family\n";
  
    if ($platform eq "linux") { # Well, this might work one day...
      push @rc, "        [[[remote]]]\n";
      push @rc, "            host = \$(rose host-select {{ EXTRACT_HOST }})\n";
      push @rc, "        [[[environment]]]\n";
      push @rc, "            ROSE_TASK_N_JOBS  = $compilejobs\n";
    }

##### OCEAN FCM_MAKE2 STEP #####  ---   No mirror on Linux.

# In the UMUI, a single job is used to build both everything, so the make 
# settings read here are the same as for the UM. In Rose they are separate tasks
# and the ocean build is much smaller, so we can use the HPC more efficiently
# by requesting a smaller wallclock time.

    if ($platform ne "linux") {
      push @rc, "    [[$oceanmake2task]]\n";
      push @rc, "        inherit = $remote_family\n";
      push @rc, "        [[[environment]]]\n";
      push @rc, "            ROSE_TASK_N_JOBS  = $compilejobs\n";

      if ($platform eq "ibm") {
        ($llheader_ref, $warnings) 
          = get_llheader(\@submit, $defaultclass, $oceanmake2task);
        $unportable += $warnings;
        @llheader = @$llheader_ref;

# Reduce values due to smaller build requirements:
        foreach my $line (@llheader) {
          if ($line =~ /wall_clock_limit = "(\d+):(\d+):(\d+),(\d+):(\d+):(\d+)"/) {
        # Set soft limit to approx. MAX(10%, 8min) (ignoring seconds):
            my $ten_percent = ($4*3600 + $5*60)*0.1;
            my $eight_minutes = 480;
            my $soft_limit = ($ten_percent >= $eight_minutes) ? $ten_percent
                                                              : $eight_minutes;
          # Convert back to xx:yy:zz
            my $hours   = sprintf("%02d",$soft_limit/3600);
            my $minutes = sprintf("%02d",($soft_limit-$hours*3600) / 60);
            my $seconds = "00";
            $soft_limit = "$hours:$minutes:$seconds";
            my $wall_clock_limit = set_wall($soft_limit,$hours, $minutes);

            $line = "            wall_clock_limit = $wall_clock_limit\n";
          }
        # Reduce memory by 50%
          if ($line =~ /ConsumableMemory\((\d+)/) {
            my $memory = $1/2;
            $line =~ s/ConsumableMemory\(\d+/ConsumableMemory\($memory/;
          }
        } # End llheader modifications
          
        push @rc, @llheader;
      } # End if (ibm)
      
    }   # End ocean make2 step

  } # End coupled-only section

##### RECONFIGURATION STEP #####

  unless ($scm) {
    push @rc, "    [[$recontask]]\n";
    push @rc, "        inherit = $remote_family\n";

    if ($platform eq "linux") {

      push @rc, "        [[[environment]]]\n";
      push @rc, "            decfort_dump_flag = y\n";
      push @rc, "            OMP_NUM_THREADS   = 1\n" if ($openmp);
      push @rc, "            OMP_STACKSIZE     = 1g\n" if ($openmp);
      push @rc, "            ROSE_LAUNCHER_ULIMIT_OPTS = -s $recon_stack -c unlimited\n";
      push @rc, "            ROSE_TASK_APP     = $runapp\n";

    } elsif ($platform eq "ibm") {

      push @rc, "        [[[environment]]]\n";
      push @rc, "            ROSE_TASK_APP    = $runapp\n";

      ($llheader_ref, $warnings) 
        = get_llheader(\@submit, $defaultclass, $recontask, $openmp);
      $unportable += $warnings;
      push @rc, @$llheader_ref;

  }   # End recon step
}   # End test for SCM


##### ATMOS STEP #####

  push @rc, "    [[$atmostask]]\n";
  push @rc, "        inherit = $remote_family\n";

  if ($platform eq "ibm") {

    push @rc, "        [[[environment]]]\n";
    push @rc, "            ROSE_TASK_APP    = $runapp\n"
      unless ($runapp eq $atmostask);  # True in the default coupled setup

    ($llheader_ref, $warnings) 
      = get_llheader(\@submit, $defaultclass, $atmostask, $openmp);
    $unportable += $warnings;
    push @rc, @$llheader_ref;

  } elsif ($platform eq "linux") {

    push @rc, "        [[[environment]]]\n";
    push @rc, "            decfort_dump_flag = y\n";
    push @rc, "            OMP_NUM_THREADS   = $omp_num_threads\n" if ($openmp);
    push @rc, "            OMP_STACKSIZE     = $omp_stacksize\n" if ($openmp);
    push @rc, "            ROSE_LAUNCHER_ULIMIT_OPTS = -s $atmos_stack -c unlimited\n";
    push @rc, "            ROSE_TASK_APP     = $runapp\n";
  }   # End atmos step

  fileio::write_file("$suitedir/suite.rc", @rc);
  return $unportable;
}


sub get_llheader {

# Retrieve/piece together a LoadLeveler header

# Arguments
  my $arrayref = shift;
    my @submit = @$arrayref;  # Contents of SUBMIT
  my $defaultclass = shift;   # The root default: serial or parallel
  my $step = shift;           # Make(2), recon or atmos/coupled/scm task
  my $openmp = shift;         # Optional flag; present = OpenMP on

# Local variables:
  my @llheader = "";    # Output for the suite.rc file
  my $search = 0;       # Flag for searching ll info
  my $warnings = 0;     # No. of unportable settings found locally

  my $minutes = 0;      # For wall_clock_limit
  my $statement = "";   # a LL statement/command

# UMUI variables that need adding/translating manually.
# Safe defaults added where appropriate.
  my $wall_clock_limit = "";
  my $class = "";
  my $memory = "";
  my $node = "";
  my $tasks_per_node = 0;
  my $total_tasks = 0;
  my $nmppe;
  my $nmppn;
  my $nemo_nproc = 0;
  my $prism_nproc = 0;
  my $flume_ios_nproc = 0;
  my $n_smt = 1;
  my $nthreads_per_task = 1;
  my $total_threads;
  my $affin = "cpu";
  my $parallel_threads = "";      # LL directive for 2+ threads
  my $omp_num_threads = "";       # Environment variable for 1 thread

# Flags for searching for LL headers:
  my $makekey  = 'cat >>\$comp_header<<EOF';
  my $reconkey = 'cat >>\$rcf_header<<EOF';
  my $umkey    = 'cat >>\$run_header<<EOF';

# Keywords to match and transfer to Rose:
# The rcf flag tasks_per_node is transferred to total_tasks for consistency
#   with the UM.
  my @keywords = qw( class job_type resources node task_affinity 
                     wall_clock_limit total_tasks node_usage network.MPI );
                    # parallel_threads

  push @llheader, "        [[[directives]]]\n";

  foreach my $line (@submit) {

# Pick up any UMUI variables we might need, and set flags marking the 
#  search limits for each LL region.
#   It's not sensible to loop over UMUI variables in the same way as keywords - 
#   too many UMUI variables have multiple possibilities (e.g. AFFIN) that need
#   special logic anyway.
    if ($step =~ /^fcm_make2/) {

      $wall_clock_limit = set_wall($1, $2, $3) 
             if ($line =~ /COMP_TIME_LIMIT=((\d+):(\d+):\d+)/);
      $search = 1 if ($line =~ /$makekey/);

    } elsif ($step =~ /^recon/) {

      $class = $1  if ($line =~ /LL_CLASS=(.*)/);
      $memory = $1 if ($line =~ /RCF_MEMORY=(.*)/);
      $node = $1   if ($line =~ /RCF_NNODES=(.*)/);
      $wall_clock_limit = set_wall($1, $2, $3) 
             if ($line =~ /RCF_TIME_LIMIT=((\d+):(\d+):\d+)/);
# Used to calculated total_tasks:
      $tasks_per_node = $1 if ($line =~ /RCF_NPES_PER_NODE=(.*)/);
      $search = 1 if ($line =~ /$reconkey/);

    } elsif ($step =~ /^(atmos|scm|coupled)/) {

      $class = $1  if ($line =~ /LL_CLASS=(.*)/);
      $memory = $1 if ($line =~ /MEMORY=(.*)/);
# CRUN support to be added...
      $wall_clock_limit = set_wall($1, $2, $3) 
             if ($line =~ /NRUN_TIME_LIMIT=((\d+):(\d+):\d+)/);

# Special cases begin...
#   Calculation for $num_nnodes, arranged by logical deduction:
      $nmppe = $1 if ($line =~ /^NMPPE=(\d+)/);    # PEs E-W
      $nmppn = $1 if ($line =~ /^NMPPN=(\d+)/);    # PEs N-S
      $nemo_nproc  = $1 if ($line =~ /^NEMO_NPROC=(\d+)/);  # Not yet supported
      $prism_nproc = $1 if ($line =~ /^PRISM_NPROC=(\d+)/); # Not yet supported
      $flume_ios_nproc = $1 if ($line =~ /^FLUME_IOS_NPROC=(\d+)/);
      $total_tasks = ($nmppe * $nmppn)    # = TOTAL_PE_REQ
                     + $nemo_nproc + $prism_nproc + $flume_ios_nproc
                   if ($line =~ /TOTAL_PE_REQ=/);  # So only calculated once
                                                   # all vars are available
      $nthreads_per_task = $1 if ($line =~ /^NTHREADS_PER_TASK=(\d+)/);
      $total_threads = $total_tasks * $nthreads_per_task
                   if ($line =~ /TOTAL_THREADS=/); # Ditto

      $n_smt = $1 if ($line =~ /^N_SMT=(\d)/);
      if ($line =~ /if test \$N_SMT = 1/) {
        if ($n_smt) {
          $node  = int(($total_threads+63) / 64);
          $affin = "cpu";
        } else {
          $node  = int(($total_threads+31) / 32);
          $affin = "core";
        }
#       Can now also calculate this:
        if ($nthreads_per_task > 1) {
          $parallel_threads = "            parallel_threads = $nthreads_per_task";
        } else {
#         Including for the no-OpenMP case:
          $omp_num_threads = "            OMP_NUM_THREADS   = 1";
       }
      }

      $search = 1  if ($line =~ /$umkey/);

    } else {

#  Not make(2), recon or atmos
      die "[FAIL] Cannot copy LoadLeveler header for unknown job step $step";

    }

#   Begin processing header info
#   Loop over each potential keyword
#   Will not process continuation lines!
    if ($search) {
      foreach my $keyword (@keywords) {

#   Potential keyword matches:
        if ($line =~ /#@\s*($keyword\s*=\s*.*)/) {
          $statement = "            $1";

#   Fill in the UMUI environment variables...
          $statement =~ s/\${(RCF_)?MEMORY}/$memory/;
          $statement =~ s/\${(RCF|NUM)_NNODES}/$node/;
          $statement =~ s/\${(TOTAL_PE_REQ)}/$total_tasks/;
          $statement =~ s/\${AFFIN}/$affin/;
          $statement =~ s/\${NTHREADS_PER_TASK}/$nthreads_per_task/;

# ...And other special cases
          $statement =~ s/=.*/= $wall_clock_limit/ 
                      if ($keyword eq "wall_clock_limit");
          $statement =~ s/=.*/= $class/
                      if ($statement =~ /LL_CLASS/);

# Remove things that are the same as a default:
          if ($keyword eq 'class' or $keyword eq 'job_type') {
            next if ($statement =~ /=\s*$defaultclass/);
          }

#    Housekeeping!
          $statement =~ s/=/ =/ if ($keyword eq "node_usage");
          if ($keyword eq "network.MPI") {   # add quotes needed by Rose
            $statement =~ s/= /= "/;
            $statement =~ s/(.$)/$1"/;
          }

# Warn against any LoadLeveler directives that may break for some users:
          if ($keyword eq 'class' and $statement =~ /=\s*(run_\w+)/) {
            print "[INFO] LoadLeveler directive 'class' set to $1 in $step task\n";
            $warnings++;
          }

# Finished this statement.
          push @llheader, "$statement\n";
        }  # End matched keyword
      } # End loop over LL keywords

      $search = 0 if ($line =~ /^EOF/);   # Stop processing lines
    }  # End of search region
  }  # End loop over SUBMIT contents

# For recon, we use total_tasks rather than tasks_per_node.
# Note that the UMUI incorrectly (though possibly deliberately) always sets 
#  tasks_per_node=32 if rcf_npes > 32, and simply scales up a node at a time.
#  (e.g. Asking for 3x16 PEs gives you 2 nodes at 32 tasks_per_node = 64 PEs.)
#  This behaviour is mimicked here to avoid anything suddenly breaking,
#  although it's unlikely such a job would work anyway due to decomposition
#  mismatches.
  if ($step =~ /^recon/) {
    $total_tasks = $tasks_per_node * $node;
    push @llheader, "            total_tasks      = $total_tasks\n";

# Number of OMP threads: the UMUI hardwires this to 1 for the reconfiguration.
#  This replaces "#@ parallel_threads = 1" both for consistency with the UM and
#  slight performance gain. Does no harm if OMP is off.
# OMP_NUM_THREADS goes in the [[[environment]]] section, hence unshift
    unshift @llheader, "            OMP_NUM_THREADS  = 1\n" if ($openmp);

  } elsif ($step =~ /^(atmos|scm|coupled)/) {
    push    @llheader, "$parallel_threads\n" if ($parallel_threads);
    unshift @llheader, "$omp_num_threads\n"  if ($omp_num_threads and $openmp);
  }

  return (\@llheader, $warnings);
}


sub set_wall {

# Process a wall clock limit by adding a hard limit
# Input: xx:yy:zz,xx,yy       (3 strings)
# Output: "uu:vv:ww,xx:yy:zz" (Quoted string)

  my $limit = shift;   # Input/output string
  my $hours = shift;   # Number of hours
  my $minutes = shift; # Number of minutes

# Add 1 min for the hardlimit, ensuring a leading zero where required.
# Softlimits >= 99:59:00 seem unlikely.
  if ($minutes eq '59') {
    $minutes = "00";                           # Reset minutes to 00
    $hours = sprintf("%02d", $hours + 1);      # and increment hours by 1
  } else {
    $minutes = sprintf("%02d", $minutes + 1);  # Increment minutes by 1
  }
  $limit =~ s/((\d+):(\d+):(\d+))/"$hours:$minutes:$4,$1"/; # Add hardlimit

  return $limit;
}

1;
