# Postscript printing procs

proc gettext { w } {
global plat plon dx dy ew ns latlc lonlc postprint prmode xshift yshift xmin ymin prcol nd_mode

set xt [wm geometry . ] 

regsub -all {[x,-,+]} $xt { } xt

set xwin [lindex $xt 0]
set ywin [lindex $xt 1]
set ywin [expr $ywin -($ymin - 2*$yshift)]
#puts stdout "$xshift $yshift"

# The meaning of latlc differs between new (UM vn>=5) and old (UM
# vn<5) dynamics.
if { $nd_mode == 0 } {
  set latlc_defn "Latitude of TLC"
  set lonlc_defn "Longitude of TLC"
} else {
  set latlc_defn "Latitude of BLC"
  set lonlc_defn "Longitude of BLC"
}
.test create text [ expr $xwin/2 ]  [ expr $ywin -15]   -text "Latitude of pole $plat; Longitude of pole $plon; Grid length $dx x $dy; Dimensions $ew x $ns; $latlc_defn $latlc; $lonlc_defn $lonlc" -tags te -font -Adobe-Times-Medium-R-Normal--*-180-*

.test create rectangle 2 2 [ expr $xwin] [ expr $ywin ] -tags te -width 2

# to file

if { $prmode == 0} {
set postprint [ $w.frame.e1.entry get]
#puts stdout $postprint
.test postscript -file $postprint -pageheight 15.0c -pagex 10.5c -pagey 14.5c -rotate 1 -colormode $prcol
}

# to printer

if { $prmode == 1 } {
set postprint [ $w.frame.e3.entry get]
# This option works, but can fail with large memory requirements
#
#set post [ .test postscript -pageheight 15.0c -pagex 10.5c -pagey 14.5c -rotate 1 ]
#exec lp -d $postprint <<$post
#
# Use intermediate file instead
set user [ exec whoami]
.test postscript -file /tmp/lampos.$user.ps -pageheight 15.0c -pagex 10.5.0c -pagey 14.5c -rotate 1 -colormode $prcol
exec lp -d $postprint /tmp/lampos.$user.ps
exec rm /tmp/lampos.$user.ps

}

destroy $w


.test delete te
#.test addtag tt all
}


proc print {{w .e1}} {
    catch {destroy $w}
    toplevel $w

    global postprint prmode home prcol
    wm title $w "Postscript file"
    wm iconname $w "File"
    wm transient $w .
    message $w.msg -font -Adobe-times-medium-r-normal--*-180* -aspect 400 \
	    -text "Enter Name of File or Printer"
    frame $w.frame -borderwidth 10
    frame $w.ok
    button $w.ok.ok -text OK -command "gettext $w"
    button $w.ok.dis -text Dismiss -command "destroy $w"
    pack $w.ok.ok $w.ok.dis -side left -fill both -expand 1
    pack $w.msg $w.frame $w.ok -side top -fill both
    frame $w.frame.e2
    frame $w.frame.e4

    labelentry $w.frame.e1 "File   " 40
    labelentry $w.frame.e3 "Printer" 40
    radiobutton $w.frame.e2.f1 -variable prmode -value 0 -text File
    radiobutton $w.frame.e2.f2 -variable prmode -value 1 -text Print
    radiobutton $w.frame.e4.f1 -variable prcol -value "gray" -text Greyscale
    radiobutton $w.frame.e4.f2 -variable prcol -value "color" -text Colour

    $w.frame.e2.f2 invoke
    $w.frame.e4.f2 invoke
    pack $w.frame.e2.f1 -side left -fill x -anchor e
    pack $w.frame.e2.f2 -side left -fill x -anchor w
    pack $w.frame.e4.f1 -side left -fill x -anchor e
    pack $w.frame.e4.f2 -side left -fill x -anchor w

    pack $w.frame.e1  -side top -pady 5 -fill x 
    pack $w.frame.e3  -side top -pady 5 -fill x
    pack $w.frame.e4  -side right -pady 5 -fill x
    pack $w.frame.e2  -side top -pady 5 -fill x
    bind $w <Any-Enter> "focus $w.frame.e1"
    $w.frame.e1.entry delete 0 end
    $w.frame.e1.entry insert 0 ~$home/lampos.ps
    $w.frame.e3.entry delete 0 end
}

