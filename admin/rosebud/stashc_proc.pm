# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************
#
# This file extracted from the UM project @ r45681
# Script modified to work as a module for rosebud
#
#  Script: rose_stashc_proc3.pl
#
#  Purpose:  Pre process the UMUI STASHC file into a format as 
#            reproposed for ROSE use and as detailed by #4280
#            Backup the current UMUI STASHC
#
#            Implementation in perl for ease of text processing.
#
# Code Owner: See Unified Model Code Owner's HTML page
# This file belongs in section: Runtime
#
# Tested under perl v5.10.1
#
# External documentation: none
#
# Interface and arguments: see echo_options
#
# End of header -------------------------------------------------------

package stashc_proc;

# Modules:
use strict;
use warnings;
###use Getopt::Long;


#-----------------------------------------------------------------
# subroutines
#-----------------------------------------------------------------
# scrape the STASHC file for profile NAMEs
# create array/list of profile name strings.

sub scrape_names {

  my @out=();
  my $profile=shift@_;
  for (@_) {
    if (/$profile/) {
      # search for string
      /NAME=\"(.+)\"/;
      push (@out, $1);
    
    } 
  }
  return \@out;  #reference of @out.
}

# specific scrape for usage so to try and identify satsh macros.
# once identified set a useful package name.

sub scrape_usage {

  our (%usage);

  my @out=();
  
  %usage=();
  my ($profile, @stashc) = @_; 
  
  chomp @stashc;                       
  my $new_line=join " ", @stashc;         
  my @split_lines= split /\//, $new_line; # split on end of namelist character.
 
  for (@split_lines) {
    if (/$profile/) {
      # search for usage string and unit number
       my $pack="";
       if (/IUNT=150/) {
         $pack = "PPVAR LS";
       } elsif (/IUNT=86/ || /IUNT=87/)  {
         $pack = "UARS"
       } elsif (/IUNT=81/ || /IUNT=82/)  {
         $pack = "SURGE"
       } elsif (/IUNT=100/ || /IUNT=101/ || /IUNT=85/)  {
         $pack = "FOAM"
       } elsif (/IUNT=102/)  {
         $pack = "CX columns and BGerr"
       } elsif (/IUNT=98/ && /NAME="UPUKCA"/)  {
         $pack = "UKCA Coupling Macro"
       } elsif (/IUNT=103/)  {
         $pack = "CEH River flow diags"
       } elsif (/IUNT=15/ && /NAME="UPMEAN"/)  {
         $pack = "Nudging code diagnostics"
       } elsif (/IUNT=164/)  {
         $pack = "MakeBC diagnostics"
       }
      /NAME=\"(.+)\"/;
      $usage{$1}=$pack;
      push (@out, $1);
      
    }
  }
  
  #while ( my ($key,$value) = each %usage) {
  #  print "$key => $value\n";
  #}  
  return \@out;  #reference of @out.
}




# rename any profile NAME duplicates 
# use a hash of names to identify any duplicates and rename
# if it has a keypair package name then ensure this is recorded for the new name
# warn if new name exceeds the UM code 8 char limit.
sub dups {

our (%usage);  # hash is set up within the scrape_usage subroutine.

  my @uniq = ();
  my %seen = ();
  my $pack="";
  
  for (@_) {
    unless ($seen{$_}) {
        # if we get here, we have not seen it before
        $seen{$_} = 1;
    } else {
        $seen{$_}++;
        if (length($_) < 8) {   # UM assumes char length max of 8
          #s/^./$seen{$_}/;
          $pack=$usage{$_};
          $_ = $_."$seen{$_}"; 
          $usage{$_}=$pack;   
        } else {
          print "[WARN] char length > 8 chars $_ \n";
          $pack=$usage{$_};
          $_ = $_."$seen{$_}";
          $usage{$_}=$pack;
        }
    }  
     push (@uniq, $_);    
  }
  return \@uniq; #reference of @uniq.
}


# update profile with unique names and write to new format STASHC.

sub profile_write{

  # subroutine input
  # 1) filename of output file
  # 2) profile type eg "TIME", "DOMAIN" or "USE".
  # 3) new name for NAME namelist entry, eg ITIM 
  # 4) reference to the list of unique profile names
  # 5) the STASHC file as an array.
  
 
  my ($stashwrite, $profile_type, $lead, $profile, @stashc) = @_;
  my $count=0;
 
  my @profile = @$profile;
  
  chomp @stashc;                       
  my $new_line=join " ", @stashc;         
  my @split_lines= split /\//, $new_line; # split on end of namelist character.
  
  for (@split_lines) {
    if (/$profile_type/) {
    s/\n/ /g;
    s/\r/ /g;
    s/^[ ]*//;  # remove leading spaces  
    my @splite = split ",", $_;
    $splite[0]=~s/NAME=\".+\"/$lead="$profile[$count]"/;  
    $_ = join "," , @splite;
    $count=$count+1;

    s/ $/\//;   # append / to end of string.      
    print $stashwrite " $_\n";  # write
    }
#   
  }
print $stashwrite "\n";
return 1;
}

#------------------------------------------------
# Subroutine: echo_options
#
# Description: Prints the script's usage syntax.
#------------------------------------------------
sub echo_options {
    print "FATAL ERROR\n";
    print "USAGE:       
             [-filein path to STASHC to be altered (default ./STASHC)]\n";
    exit 1;
}
#------------------------------------------------
# Subroutine: preproc_forwardslash
#
# Description: Replaces forward slashes with underscore to prevent the script
#              mistaking them for namelist end tokens.
#------------------------------------------------
sub preproc_forwardslash {
    my @lines = @_;
    foreach my $line (@lines) {
      if ($line =~/NAME\s*=\s*\".+?\/.+\"/) {
        $line =~ s!(NAME\s*=\s*\".+?)/(.+\")!$1_$2!;
      }
    }
    return @lines;
}

#------------------------------------------------
# Subroutine: stash_rose
#
# Description: checks stashc to process is not already ROSE converted.
#------------------------------------------------
sub stash_rose {
    my @lines = @_;
    for (@lines){
      if (/! ROSE STASHC version/) {
      print "[WARN] STASHC file is already of ROSE type.\n";
      print "[WARN] No STASH processing done.\n";
      print "[FAIL] No suite generated due to unexpected STASH settings.\n";
      exit
      }
    }
    return 1  
}

##### Converted to a subroutine for use with the Rosebud script

sub process_stash {

#------------------------------------------------------------------------
# Inputs to script
# check validty of usage
#------------------------------------------------------------------------

#  my $filein='./STASHC';
   my $filein = shift;

#  my $return= GetOptions ('filein=s' => \$filein);

# check for errors from GetOpts, if return code is not true stop.
#  unless ($return) {
#          echo_options();
#  }

# check input STASHC file exists
  unless (-e $filein) {
          print "[FAIL] STASHC file $filein doesn't exist!\n";
          echo_options();
} 

# the user does not provide any filenames without using options.
#  if ( $#ARGV ne -1) { echo_options(); }


#------------------------------------------------------------------------
# Check stash format, is it already ROSE complient?
#------------------------------------------------------------------------
# open STASHC file and load into array.
  open my $stashin,'< ', $filein or die "[FAIL] Cannot open file $filein $!" ; 
  my @lines=<$stashin>;
  close $stashin;

  my $process_req = stash_rose(@lines);

#------------------------------------------------------------------------
# Set fileout to be filein and rename filein
# process the STASHC file
#------------------------------------------------------------------------

  our (%usage);  # hash is set up within the scrape_usage subroutine.
  my $fileout=$filein;
  my $file_bkup=$filein."_preROSE";

  my @rename=();
  @rename=`mv $filein $file_bkup`;

  if ($?) {
          print "[FAIL] Unable to create backup STASHC file\n";
          exit;
  }else{
#          print "-------------------------------\n";
#          print "Processing STASHC file for ROSE\n";
#          print "-------------------------------\n";
#          print "Backup of STASHC file created! \n";
#          print $file_bkup, "\n"
   }
    

# open STASHC file and load into array.
  open $stashin,'< ', $file_bkup or die "[FAIL] Cannot open file $file_bkup $!" ; 
  @lines=<$stashin>;
  close $stashin;

# open file for reformatted STASHC output.
  open my $stashout,'> ', $fileout or die "[FAIL] Cannot open file $fileout $!" ; 

# Replace forward slashes in names with underscore to prevent errors later
  @lines = preproc_forwardslash(@lines);

  print $stashout "! ROSE STASHC version\n";

#---------------------------------------------------------------
# initialise working arrays for the items to be processed
  my $NUM_DOM_ADD=0;  # domain profile counter
  my $NUM_TIM_ADD=0;  # time profile counter
  my $NUM_USE_ADD=0;  # use profile counter

  my @NUM_DOM=(0);  # array containing stashnum counters
  my @NUM_TIM=(0);
  my @NUM_USE=(0);
  my $i=0;

# scrape the available profile names 
# and create unique set of names for each profile list.

  my $time= scrape_names("&TIME",@lines);   # ref to time
     $time= dups(@$time);
  my $domain= scrape_names("&DOMAIN",@lines); #ref to domain
     $domain= dups(@$domain);
  my $use= scrape_usage("&USE",@lines);# ref to use
     $use= dups(@$use);
#---------------------------------------------------------------

#processing of the stash requests, applying unique name tags 
# as scraped in the above.

  for (@lines) {

# need to know which batch of stash requests being updated
# search for STASHNUM create array of number of profiles 
    if (/STASHNUM/) {
      /NUM_DOM=[ ]*([0-9]+),[ ]*NUM_TIM=[ ]*([0-9]+),[ ]*NUM_USE=[ ]*([0-9]+)/;
      push (@NUM_DOM, $1);
      push (@NUM_TIM, $2);
      push (@NUM_USE, $3);
    }  
  
# update stash requests with the unique name tags.  
    if (/\&STREQ/) { 
        s/IMOD=[ ]*[0-9],//; #remove unecessary IMOD detail
        s/IDOM=([0-9]+)/DOM_NAME=\"$$domain[$1+$NUM_DOM_ADD-1]\"/;
        s/ITIM=([0-9]+)/TIM_NAME=\"$$time[$1+$NUM_TIM_ADD-1]\"/;
        s/IUSE=([0-9]+)/USE_NAME=\"$$use[$1+$NUM_USE_ADD-1]\", PACKAGE=\"$usage{$$use[$1+$NUM_USE_ADD-1]}\",/;
        print $stashout $_; 
    }

# add to total number of profiles for next block of stashnum.
    if (/STASHNUM/) {
    $NUM_DOM_ADD=$NUM_DOM_ADD+$NUM_DOM[$i];
    $NUM_TIM_ADD=$NUM_TIM_ADD+$NUM_TIM[$i];
    $NUM_USE_ADD=$NUM_USE_ADD+$NUM_USE[$i];
    $i=$i+1;
    }  
 
  }  
  print $stashout "\n";


#---------------------------------------------------------------
# append the stash request output with updated profiles using unique names.
  profile_write($stashout,"&DOMAIN","DOM_NAME",$domain,@lines) ;
  profile_write($stashout,"&TIME","TIM_NAME",$time,@lines) ;
  profile_write($stashout,"&USE","USE_NAME",$use,@lines) ;

  if ($?) {
          print "[FAIL] STASHC processing failed\n";
          exit;
#  }else{
#          print "-------------------------------\n";
#          print "Processing STASHC file complete\n";
#          print "-------------------------------\n";
  }


  close ($stashout)

}
#----------------------------------------------------------------
1;
