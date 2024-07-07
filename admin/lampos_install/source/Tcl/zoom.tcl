# Zooming and panning routines

proc zoomin { } {

global zoom xshift yshift xpan ypan
set zoom [ expr $zoom*1.2]
.test scale tg [expr $xshift +$xpan] [expr $yshift+$ypan] 1.2 1.2
.test scale ta [expr $xshift +$xpan] [expr $yshift+$ypan] 1.2 1.2
.test scale tl [expr $xshift +$xpan] [expr $yshift+$ypan] 1.2 1.2
.test scale tf [expr $xshift +$xpan] [expr $yshift+$ypan] 1.2 1.2
#puts stdout "$zoom"

}

proc zoomout { } {
global zoom xshift yshift xpan ypan
set zoom [ expr $zoom/1.2]
if { $zoom < 1.0 } {
    set zoom 1.0
    return
}
.test scale tg [expr $xshift +$xpan] [expr $yshift+$ypan] [ expr 1./1.2] [ expr 1./1.2]
.test scale ta [expr $xshift +$xpan] [expr $yshift+$ypan] [ expr 1./1.2] [ expr 1./1.2]
.test scale tl [expr $xshift +$xpan] [expr $yshift+$ypan] [ expr 1./1.2] [ expr 1./1.2]
.test scale tf [expr $xshift +$xpan] [expr $yshift+$ypan] [ expr 1./1.2] [ expr 1./1.2]
#puts stdout "$zoom"

}
proc pandown {} {
    global ypan
     set ypan [expr $ypan+10]
     .test move tt 0 10
}
proc panup {  } {
     global ypan
     set ypan [expr $ypan-10]
     .test move tt 0 -10
}
proc panleft { } {
     global xpan
    set xpan [expr $xpan-10]
    .test move tt -10 0
}
proc panright { } {
     global xpan
     set xpan [expr $xpan+10]
    .test move tt 10 0
}

# Utility procedures for panning screen with mouse.

proc grabScreen {c x y} {
    global areaX1 areaY1
    set areaX1 [$c canvasx $x]
    set areaY1 [$c canvasy $y]

}
proc moveScreen {c x y} {
    global areaX1 areaY1 areaX2 areaY2
    set areaX2 [$c canvasx $x]
    set areaY2 [$c canvasy $y]
     global xpan ypan
     set xpan [expr $xpan+$areaX2-$areaX1]
     set ypan [expr $ypan+$areaY2-$areaY1]
    .test move tt [expr $areaX2-$areaX1] [expr $areaY2-$areaY1]
     set areaX1 $areaX2
     set areaY1 $areaY2
    
}
