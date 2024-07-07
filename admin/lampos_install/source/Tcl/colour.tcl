
proc clrr { {w .clrpick } } {
catch {destroy $w}
toplevel $w
wm transient $w .
wm title $w "Color Selection Dialog"
wm iconname $w "colors"
set font -Adobe-times-medium-r-normal--*-180*
#label $w.msg -font $font -wraplength 4i -justify left -text "Press the buttons below to choose the colors for the folloing widgets."
#pack $w.msg -side top

frame $w.buttons
pack $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text OK -command "dissWin $w"
pack $w.buttons.dismiss -side left -expand 1
set null " "
frame $w.lat
frame $w.coast
frame $w.rect

button $w.back -text "background" \
    -command \
    "setColor .test $w.back background configure -background z col_back"

button $w.lat.but -text "Colour" \
    -command \
    "setColor .test $w.back foreground itemconfigure tl -fill col_lat"

button $w.coast.but -text "Colour" \
    -command \
    "setColor .test $w.back foreground itemconfigure tg -fill col_coast"

button $w.rect.but -text "Colour" \
    -command \
    "setColor .test $w.back foreground itemconfigure ta -outline col_rect"


label $w.lat.txt -justify left -text "Lat-long" -anchor e
label $w.coast.txt -justify left -text "Coast    " -anchor e
label $w.rect.txt -justify left -text "Area    " -anchor e

labelentry $w.lat.lab "Line Width" 3
labelentry $w.coast.lab "Line Width" 3
labelentry $w.rect.lab "Line Width" 3

pack $w.lat.txt -expand 1 -fill x -side left
pack $w.lat.but -expand 1 -fill x -side left
pack $w.lat.lab -side left -expand 1 -fill x

pack $w.coast.txt -side left -expand 1 -fill x
pack $w.coast.but -side left -expand 1 -fill x
pack $w.coast.lab -side left -expand 1 -fill x

pack $w.rect.txt -side left -expand 1 -fill x
pack $w.rect.but -side left -expand 1 -fill x
pack $w.rect.lab -side left -expand 1 -fill x

pack $w.back $w.coast $w.lat $w.rect -side top -anchor c -pady 2m
global width_coast width_lat width_rect
$w.lat.lab.entry insert 0 $width_lat
$w.coast.lab.entry insert 0 $width_coast
$w.rect.lab.entry insert 0 $width_rect

      set width_coast [ $w.coast.lab.entry get]
      set width_lat [ $w.lat.lab.entry get]
      set width_rect [ $w.rect.lab.entry get]

}

proc dissWin {w} {
global width_coast width_lat width_rect

      set width_coast [ getround $w.coast.lab $width_coast 1 5 3 ]
      set width_lat   [ getround $w.lat.lab $width_lat 1 5 3 ]
      set width_rect   [ getround $w.rect.lab $width_rect 1 5 3 ]

     .test itemconfigure tl -width $width_lat
     .test itemconfigure tg -width $width_coast
     .test itemconfigure ta -width $width_rect

      destroy $w
}


proc setColor {w button name option1 option2 option3 a } {

global col_back col_lat col_coast col_rect

    if { [string compare $a "col_back"] == 0 } {
      set initialColor $col_back
    }
    if {[string compare $a "col_coast"] == 0 } {
      set initialColor $col_coast
    }
    if {[string compare $a "col_lat"] == 0 } {
      set initialColor $col_lat
    }
    if {[string compare $a "col_rect"] == 0 } {
      set initialColor $col_rect
    }

    grab $w
#    set initialColor [$button cget -$name]
    set color [tk_chooseColor -title "Choose a $name color" -parent $w \
        -initialcolor $initialColor]
    if [string compare $color ""] {


    if { [string compare $a "col_back"] == 0 } {
      $w $option1 $option2 $color
      set col_back $color 
    }
    if { [string compare $a "col_coast"] == 0 } {

      $w $option1 $option2 $option3 $color
      set col_coast $color 
    }
    if { [string compare $a "col_lat"] == 0 } {

     $w $option1 $option2 $option3 $color 
      set col_lat $color 
    }
    if { [string compare $a "col_rect"] == 0 } {

     $w $option1 $option2 $option3 $color 
      set col_rect $color
      cross 
    }

    }
    grab release $w


}

