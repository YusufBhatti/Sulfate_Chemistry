proc main { } {
# Main program

#------------------------------------------------------------
# Control Display Area
#------------------------------------------------------------
frame .view

#------------------------------------------------------------
# LAM area control panel
#------------------------------------------------------------

frame .view.grd -borderwidth 2 -relief groove

frame .view.grd.scr
frame .view.grd.top
frame .view.grd.scr.lc
frame .view.grd.scr.lw
frame .view.grd.scr.le

label .view.grd.top.title -text "Grid" -relief groove -padx 20 -pady 6

labelentry .view.grd.scr.lc.ew "Number of col's" 4
labelentry .view.grd.scr.lc.ns "Number of rows" 4
labelentry .view.grd.scr.lw.dx "Col' spacing" 5
labelentry .view.grd.scr.lw.dy "Row spacing" 5
labelentry .view.grd.scr.le.latlc "First lat" 7
labelentry .view.grd.scr.le.lonlc "First lon" 7
button .view.grd.apply -command "drawarea" -text "Apply" -takefocus 0

# If user presses return in any window then interpet this as apply
bind .view.grd.scr.le.lonlc.entry <KeyPress> \
   {if {[string compare "%A" "\r"] == 0} {drawarea}}
bind .view.grd.scr.le.latlc.entry <KeyPress> \
   {if {[string compare "%A" "\r"] == 0} {drawarea}}
bind .view.grd.scr.lc.ew.entry <KeyPress> \
   {if {[string compare "%A" "\r"] == 0} {drawarea}}
bind .view.grd.scr.lc.ns.entry <KeyPress> \
   {if {[string compare "%A" "\r"] == 0} {drawarea}}
bind .view.grd.scr.lw.dx.entry <KeyPress> \
   {if {[string compare "%A" "\r"] == 0} {drawarea}}
bind .view.grd.scr.lw.dy.entry <KeyPress> \
   {if {[string compare "%A" "\r"] == 0} {drawarea}}

menubutton .view.grd.areas -text "Predefined LAM areas..." -relief raised -menu .view.grd.areas.m
  menu .view.grd.areas.m
 .view.grd.areas.m add command -label "12km NAE" -command "setnae"
 .view.grd.areas.m add command -label "4km UK Mes" -command "setuk4"
 .view.grd.areas.m add command -label "Previous area" -command "setarea"

button .view.grd.drag -command "setm1" -text "Drag New Area" -takefocus 0

pack .view.grd.top .view.grd.top.title
pack .view.grd.scr.le.latlc .view.grd.scr.le.lonlc -side top -pady 2 \
  -padx 2 -anchor e
pack .view.grd.scr.lc.ew .view.grd.scr.lc.ns -side top -pady 2 -padx 2 \
  -anchor e
pack .view.grd.scr.lw.dx .view.grd.scr.lw.dy -side top -pady 2 -padx 2 \
   -anchor e
pack .view.grd.scr.lc .view.grd.scr.lw .view.grd.scr.le -side left
pack  .view.grd.scr
pack .view.grd.apply .view.grd.areas .view.grd.drag -side left -expand yes

#------------------------------------------------------------
# Pole control
#------------------------------------------------------------

frame .view.pl -borderwidth 2 -relief groove

frame .view.pl.top
frame .view.pl.scr
frame .view.pl.scr.le
frame .view.pl.scr.lw

label .view.pl.top.title -text "Coords of Rotated Pole" -relief groove -padx 20 -pady 6

# Pole lat+lon entry boxes
labelentry .view.pl.scr.le.plat "Latitude " 7
labelentry .view.pl.scr.le.plon "Longitude" 7

# And the apply-button
button .view.pl.scr.lw.apply -command "drawmap" -text "Apply" -takefocus 0

# If user presses return in either window then interpet this as apply
bind .view.pl.scr.le.plon.entry <KeyPress> \
   {if {[string compare "%A" "\r"] == 0} {drawmap}}
bind .view.pl.scr.le.plat.entry <KeyPress> \
   {if {[string compare "%A" "\r"] == 0} {drawmap}}

# Menu for selecting pre-defined poles
menubutton .view.pl.poles -text "Predefined LAM poles..." -relief raised \
  -menu .view.pl.poles.m
menu .view.pl.poles.m
 .view.pl.poles.m add command -label "12km NAE" -command "setnaep"
 .view.pl.poles.m add command -label "4km UK Mes" -command "setuk4p"
 .view.pl.poles.m add command -label "Standard lat-lon" -command "setstanp"
 .view.pl.poles.m add command -label "Previous coords" -command "setpole;cross"

# This button repositions the pole so that the point clicked lies on 
# the intersection of the meridian and equator.
button .view.pl.click -command "setm2" -text "Centre of View" -takefocus 0

# Define the relative positions of the components of the pole 
# control section
pack .view.pl.top .view.pl.top.title
pack .view.pl.scr.le.plat .view.pl.scr.le.plon -side top \
  -pady 2 -padx 2 -anchor e
pack .view.pl.scr.lw.apply -side top -pady 2 -padx 2 -anchor e
pack .view.pl.scr.lw .view.pl.scr.le -side left
pack .view.pl.scr
pack .view.pl.poles .view.pl.click -side left -expand yes

#------------------------------------------------------------
# Control for reading a grid from a UMUI job
#------------------------------------------------------------

frame .view.rd -borderwidth 2 -relief groove

frame .view.rd.scr
frame .view.rd.scr.le
frame .view.rd.scr.lw

label .view.rd.title -text "Read in grid from UM job" -relief groove -padx 20 -pady 6

# Username
labelentry .view.rd.scr.le.user "Owner's username" 6
# Job ID
labelentry .view.rd.scr.le.job "Job ID" 6
# If user presses return in either window, then interpret this as Apply
bind .view.rd.scr.le.user.entry <KeyPress> \
   {if {[string compare "%A" "\r"] == 0} {getgridfromjob}}
bind .view.rd.scr.le.job.entry <KeyPress> \
   {if {[string compare "%A" "\r"] == 0} {getgridfromjob}}

# Apply button
button .view.rd.scr.lw.apply -command "getgridfromjob" -text "Read & apply" -takefocus 0

# A comment
label .view.rd.comment -text "(reads ~user/umui_jobs/jobID/SIZES)"

# Define the relative positions of the components
pack .view.rd.scr.le.user .view.rd.scr.le.job -side top -pady 2 -padx 2 \
   -anchor e
pack .view.rd.scr.lw.apply -side top -pady 2 -padx 2 -anchor e
pack .view.rd.scr.le .view.rd.scr.lw -side left
pack .view.rd.title .view.rd.scr .view.rd.comment

#------------------------------------------------------------
# Define the relative positions of the above grid control sections
#------------------------------------------------------------
pack .view.grd -side left -fill both -expand 1
pack .view.pl -side left -fill both -expand 1
pack .view.rd -side left -fill both -expand 1

#------------------------------------------------------------
# Viewing controls
#------------------------------------------------------------

global zoom; set zoom 1.

frame .mou -relief raised -borderwidth 2

# Quit button
button .mou.quit -command "writeOptions ; destroy ." -text "Quit" -takefocus 0

# Options menu
menubutton .mou.options -text "Options..." -menu .mou.options.m \
  -relief raised -borderwidth 2
menu .mou.options.m
.mou.options.m add check -label "Display lat-lon" -variable latlon \
   -onvalue 1 -offvalue 0 -command "drawlatlon"
.mou.options.m add check -label "Display grid" -variable grid -onvalue 1 \
   -offvalue 0 -command "drawarea"
.mou.options.m add check -label "Lowres coasts" -variable coast -onvalue 0 \
   -offvalue 1 -command "drawmap"
.mou.options.m add command -label "Colours & line widths" -command "clrr"
.mou.options.m add command -label "Print" -command "print"
#.mou.options.m add command -label "Save set-up" -command "writeOptions"

# Start in New Dynamics mode
global nd_mode; new_dynamics_mode 1

# Variable to indicate (to Tk Widgets) if we are in new-dynamics mode
# (variable to indicate this to other parts of Lampos is $nd_mode)
global mode_indctr ; set mode_indctr 1

# ND vs OD menu
menubutton .mou.nd_od -text "ND vs. OD..." -menu .mou.nd_od.m  -relief raised -borderwidth 2    
menu .mou.nd_od.m
.mou.nd_od.m add radiobutton -label \
   "ND mode (first lat=lat of bottom of area)" \
   -command "new_dynamics_mode 1" -variable mode_indctr -value 1
.mou.nd_od.m add radiobutton -label "OD mode (first lat=lat of top of area)" \
   -command "new_dynamics_mode 0" -variable mode_indctr -value 0
.mou.nd_od.m add separator
   .mou.nd_od.m add command -label \
   "Add (nrows-1)*row_spacing to first latitude" -command "frstlat_nd2od"
.mou.nd_od.m add command -label \
   "Subtract (nrows-1)*row_spacing from first latitude" \
   -command "frstlat_od2nd"

button .mou.help -command "gv_helpbrowser" -text "Help" -takefocus 0
button .mou.zoomin -command "zoomin"       -text "Zoom In" -takefocus 0
button .mou.zoomou -command "zoomout"      -text "Zoom Out" -takefocus 0

button .mou.flip -command "flipview" -text "Flip View" -takefocus 0

# Grid options menu
#menubutton .mou.grid -text "Grid options..." -menu .mou.grid.m \
#  -relief raised -borderwidth 2
#  menu .mou.grid.m
#  .mou.grid.m add check -label "Display lat-lon" -variable latlon \
#     -onvalue 1 -offvalue 0 -command "drawlatlon"
#  .mou.grid.m add check -label "Display grid" -variable grid -onvalue 1 \
#     -offvalue 0 -command "drawarea"
#  .mou.grid.m add check -label "Lowres coasts" -variable coast -onvalue 0 \
#     -offvalue 1 -command "drawmap"

# Define the relative positions of the buttons + menus
pack .mou.zoomin -pady 2 -padx 0  -side left -expand yes
pack .mou.zoomou -pady 2 -padx 0 -side left -expand yes
pack .mou.flip -pady 2 -padx 0 -side left -expand yes
#pack .mou.view -pady 2 -padx 0 -side left -expand yes
#pack .mou.grid -pady 2 -padx 0 -side left -expand yes
#pack .mou.print -pady 2 -padx 0 -side left -expand yes
pack .mou.options -pady 2 -padx 0 -side left -expand yes
pack .mou.nd_od -pady 2 -padx 0 -side left -expand yes
pack .mou.help -pady 2 -padx 0 -side left -expand yes
pack .mou.quit -pady 2 -padx 0 -side left -expand yes

#------------------------------------------------------------
# Other stuff
#------------------------------------------------------------

global xshift yshift
    set sss [ wm maxsize .]
    scan $sss "%d %d" xshift yshift
    set xshift [expr 0.40*$xshift]
    set yshift [expr 0.35*$yshift]

global scale; set scale [ expr $xshift / 180. ]

global ratio; set ratio 1

global filelat; set filelat "data_latlon"
    setnae
    setnaep
canvas .test -height [expr 2 * $yshift] -width [expr 2 * $xshift] -background white

pack  .view  .mou -fill x  -side bottom
pack .test  -expand 1 -fill both -anchor se

global xpan;set xpan 0
global ypan; set ypan 0
    bind .test <KeyPress-Down>  "pandown"
    bind .test <KeyPress-Up>  "panup"
    bind .test <KeyPress-Left>  "panleft"
    bind .test <KeyPress-Right>  "panright"
    bind . <Any-Enter> "focus .test"
    bind .test <2> "focus .test; grabScreen .test %x %y"
    bind .test <B2-Motion> "moveScreen .test %x %y"

global home 
set home [exec whoami]

# Set colours
global col_back col_lat col_coast col_rect latlon width_coast width_lat width_rect

set col_back #f8f4cc
set col_lat  #f29052
set col_coast #2e0a9a
set col_rect black
set latlon 1
set width_coast 2
set width_lat 2
set width_rect 3

readOptions


# Start up with tk logo

image create photo logo -file  ../../data/pwrdLogo200.gif 
.test create image $xshift $yshift -image logo -tag tl

setdisabled
update idletasks

global xmin ymin
set xt [wm geometry . ] 
regsub -all {[x,-,+]} $xt { } xt
set xmin [lindex $xt 0]
set ymin [lindex $xt 1]
wm minsize . $xmin $ymin
#puts stdout "$xmin $ymin"
wm aspect . [expr $xmin] $ymin [expr $xmin] $ymin
wm title . "LAMPOS - Version 4.2"

zoomin;zoomin;zoomin;zoomin;zoomin;zoomin
.test configure -background $col_back
drawmap

setnormal
if {$latlon == 0} {.test delete tl}

}


#################################################################
proc new_dynamics_mode { mode } {
#################################################################
# Set nd_mode to indicate new or old dynamics mode (i.e. whether a
# grid's position is defined by it's bottom left or top left corner
# respectively.)

# This variable is referred to by the radiobuttons when Tk is figuring
# out if they should be highlighted.
global mode_indctr
set mode_indctr $mode

#if { $mode == 1 } {
#  .mou.nd_od.m entryconfigure 5 -state disabled
#  .mou.nd_od.m entryconfigure 4 -state normal
#} else {
#  .mou.nd_od.m entryconfigure 4 -state disabled
#  .mou.nd_od.m entryconfigure 5 -state normal
#}

global nd_mode
set nd_mode $mode

}

#################################################################
proc frstlat_od2nd { } {
#################################################################
# Change first latitude from latitude of TLC (OD) to that of BLC (ND)
global nd_mode latlc dy ns
#if { $nd_mode == 0 } {
  setarea
  set latlc [ expr $latlc - ( $ns - 1 ) * $dy ]
  .view.grd.scr.le.latlc.entry delete 0 end
  .view.grd.scr.le.latlc.entry insert 0 $latlc
#}
}
#################################################################
proc frstlat_nd2od { } {
#################################################################
# Change first latitude from latitude of BLC (ND) to that of TLC (OD)
global nd_mode latlc dy ns
#if { $nd_mode == 1 } {
  setarea
  set latlc [ expr $latlc + ( $ns - 1 ) * $dy ]
  .view.grd.scr.le.latlc.entry delete 0 end
  .view.grd.scr.le.latlc.entry insert 0 $latlc
#}
}
#################################################################

#################################################################
proc getgridfromjob { } {
#################################################################
# Given a username and a job name, extract the grid parameters and
# reset Lampos accordingly.

global nd_mode

# Record the original value of nd_mode
set original_nd_mode $nd_mode

# Get the job and user from the entry boxes
set job [.view.rd.scr.le.job.entry get]
set user [.view.rd.scr.le.user.entry get]

## Attempt to get the user's home area from the passwd file
#set user_home [exec grep $user /etc/passwd | cut -f 6 -d :]

# Attempt to get the user's home area by tilde-expansion
set user_home [exec /usr/bin/env bash -c "eval echo ~$user"]

# Does the user's home directory exist?
set exist [exec ls -dl $user_home]

# The grid is in ~/umui_jobs/SIZES
set file $user_home/umui_jobs/$job/SIZES


# Does the file exist?
set exist [exec ls -l $file]

# Read the grid parameters from the file with grep and perl...

# First try to extract row_length assuming New Dynamics formatting.
catch { exec grep -i global_ROW_LENGTH= $file | perl -pe s{.*global_ROW_LENGTH=\\s*(\\S+),.*}{\\1}i } row_length
# Did this work?
if { [string is integer $row_length] == 0 } {
  # No. Try the Old Dynamics format
  catch { exec grep -i ROW_LENGTH= $file | perl -pe s{.*ROW_LENGTH=\\s*(\\S+),.*}{\\1}i } row_length 
  if { [string is integer $row_length] == 0 } {
    # Did not find it this way either.
    error "Could not determine grid from file"
    return
  } else {
    # It appears we are in Old Dyanmics mode. 
    new_dynamics_mode 0
  }
} else {
    # It appears we are in New Dyanmics mode. 
    new_dynamics_mode 1
}
# The formatting of the rows variable in the file depends on whether
# this is Old or New Dynamics
if { $nd_mode == 1 } {
  set rows [exec grep -i global_ROWS= $file | perl -pe s{.*global_ROWS=\\s*(\\S+),.*}{\\1}i] 
} else {
  set rows [exec grep -i P_ROWS= $file | perl -pe s{.*P_ROWS=\\s*(\\S+),.*}{\\1}i] 
}
set dx [exec grep -i EWSPACEA= $file | perl -pe s{.*EWSPACEA=\\s*(\\S+),.*}{\\1}i]
set dy [exec grep -i NSSPACEA= $file | perl -pe s{.*NSSPACEA=\\s*(\\S+),.*}{\\1}i]
set latlc [exec grep -i FRSTLATA= $file | perl -pe s{.*FRSTLATA=\\s*(\\S+),.*}{\\1}i]
set lonlc [exec grep -i FRSTLONA= $file | perl -pe s{.*FRSTLONA=\\s*(\\S+),.*}{\\1}i]
set plat [exec grep -i POLELATA= $file | perl -pe s{.*POLELATA=\\s*(\\S+),.*}{\\1}i]
set plon [exec grep -i POLELONA= $file | perl -pe s{.*POLELONA=\\s*(\\S+),.*}{\\1}i]

# Set the parameters
.view.grd.scr.le.lonlc.entry delete 0 end
.view.grd.scr.le.lonlc.entry insert 0 $lonlc
.view.grd.scr.le.latlc.entry delete 0 end
.view.grd.scr.le.latlc.entry insert 0 $latlc
.view.grd.scr.lc.ew.entry delete 0 end
.view.grd.scr.lc.ew.entry insert 0 $row_length
.view.grd.scr.lc.ns.entry delete 0 end
.view.grd.scr.lc.ns.entry insert 0 $rows
.view.grd.scr.lw.dx.entry delete 0 end
.view.grd.scr.lw.dx.entry insert 0 $dx
.view.grd.scr.lw.dy.entry delete 0 end
.view.grd.scr.lw.dy.entry insert 0 $dy
.view.pl.scr.le.plon.entry delete 0 end
.view.pl.scr.le.plon.entry insert 0 $plon
.view.pl.scr.le.plat.entry delete 0 end
.view.pl.scr.le.plat.entry insert 0 $plat

# Redraw
drawmap
drawarea

# Did we change from New to Old Dynamics mode or vice versa?
if { $original_nd_mode == 1 && $nd_mode == 0 } {
  tk_messageBox -message "Switched to Old Dynamics mode." -type ok
}
if { $original_nd_mode == 0 && $nd_mode == 1 } {
  tk_messageBox -message "Switched to New Dynamics mode." -type ok
}

}

