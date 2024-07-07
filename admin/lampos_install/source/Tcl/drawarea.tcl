# Draw outline of limited area (or grid), deleting previous one

proc drawarea { } {

global ew ns dx dy lonlc latlc scale xshift yshift zoom grid xpan ypan
global col_rect width_rect
global nd_mode

    .test delete ta
    .test delete tf

# Draw LAM outline

    if {$grid == 0} {
    set ew      [ getround .view.grd.scr.lc.ew $ew 1 9999 4]
    set ns      [ getround .view.grd.scr.lc.ns $ns 1 9999 4]
    set dx	[ getround .view.grd.scr.lw.dx $dx 0. 60. 5]
    set dy	[ getround .view.grd.scr.lw.dy $dy 0. 60. 5]
    set latlc	[ getround .view.grd.scr.le.latlc $latlc -90. 90. 7]
    set lonlc	[ getround .view.grd.scr.le.lonlc $lonlc 0. 360. 7]
    if {$lonlc > 180. } {set lonlc [expr -360. + $lonlc]}

    set x1 [ expr $lonlc * $scale + $xshift+$xpan/$zoom]
    set x2 [ expr $x1 + $dx * ($ew - 1) * $scale]
    if {$nd_mode == 0} {
      # latlc holds the top-edge latitude of the area
      set y1 [ expr -$latlc * $scale + $yshift+$ypan/$zoom]
      set y2 [ expr $y1 + $dy * ($ns - 1) * $scale]
    } else {
      # latlc holds the bottom-edge latitude of the area
      set y2 [ expr -$latlc * $scale + $yshift+$ypan/$zoom]
      set y1 [ expr $y2 - $dy * ($ns - 1) * $scale]
    }

    .test create rectangle $x1 $y1 $x2 $y2 -tags ta -width $width_rect -outline $col_rect
    .test scale ta $xshift $yshift $zoom $zoom

    }

# Draw LAM grid

    if { $grid == 1} {
    set skip 1
    set x1 [ expr $lonlc * $scale + $xshift+$xpan/$zoom]
    set y1 [ expr -$latlc * $scale + $yshift + $ypan/$zoom]
    if {$nd_mode == 1} {
      # latlc holds bottom-edge, not top-edge, latitude. Compensate.
      set y1 [ expr $y1 - $dy * ($ns - 1) * $scale ]
    }
    set zoom2 [ expr .01]
    for {set i 1} {$i <= $ew} {set i [expr $i+$skip]} {
    for {set j 1} {$j <= $ns} {set j [expr $j+$skip]} {
    set x2 [ expr $x1 + $dx * ($ew - $i) * $scale]
    set y2 [ expr $y1 + $dy * ($ns - $j) * $scale]
    .test create rectangle [ expr $x2 -$zoom2] [ expr $y2 -$zoom2] [ expr $x2 + $zoom2] [ expr $y2 + $zoom2] -tags tf

    }
    }
    .test scale tf $xshift $yshift $zoom $zoom
    }

    .test addtag tt all
}

