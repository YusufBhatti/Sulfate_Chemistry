

# Utility procedures for stroking out a rectangle and printing what's
# underneath the rectangle's area.

proc itemMark {c x y} {
    global areaX1 areaY1
    set areaX1 [$c canvasx $x]
    set areaY1 [$c canvasy $y]
    $c delete area
}

proc itemStroke {c x y} {
    global areaX1 areaY1 areaX2 areaY2 nd_mode
    set x [$c canvasx $x]
    set y [$c canvasy $y]
    if {($areaX1 != $x) && ($areaY1 != $y)} {
	$c delete area

	$c addtag area withtag [$c create line $areaX1 $areaY1 $areaX1 $y $x $y $x $areaY1 $areaX1 $areaY1\
		-width 2 -fill red -stipple @../../data/grey.5]
	set areaX2 $x
	set areaY2 $y
    }

}

proc itemsUnderArea {c} {
    global areaX1 areaY1 areaX2 areaY2
    global scale xshift yshift lonlc latlc dx dy ew ns zoom xpan ypan nd_mode

    # Make sure that Y2>Y1 and X2>X1
    if {$areaY2 < $areaY1} {
      set tmp $areaY1
      set areaY1 $areaY2
      set areaY2 $tmp
    }
    if {$areaX2 < $areaX1} {
      set tmp $areaX1
      set areaX1 $areaX2
      set areaX2 $tmp
    }

    $c delete area
    $c delete ta
    $c create rectangle $areaX1 $areaY1 $areaX2 $areaY2 -tags ta -width 2 -outline red
    set lonlcl [ format %.2f [ expr ( $areaX1 -$xshift -$xpan ) / ($scale*$zoom)] ]
    # If in New Dynamics mode then latlcl should be the bottom latitude, 
    # else it shuld be the top
    if {$nd_mode == 1} {
      set latlcl [ format %.2f [expr  ( $yshift + $ypan - $areaY2) / ($scale*$zoom)] ]
    } {
      set latlcl [ format %.2f [expr  ( $yshift + $ypan - $areaY1) / ($scale*$zoom)] ]
    }

    set ewl [ format %.0f [ expr ( $areaX2-$areaX1)/ ( $dx *  $scale*$zoom) +1] ]
    set nsl [ format %.0f [ expr ( $areaY2 -$areaY1) / ( $dy * $scale*$zoom) + 1] ]
    if {$lonlcl < 0. } {set lonlcl [expr 360. + $lonlcl]}

    .view.grd.scr.le.lonlc.entry delete 0 end
    .view.grd.scr.le.lonlc.entry insert 0 $lonlcl
    .view.grd.scr.le.latlc.entry delete 0 end
    .view.grd.scr.le.latlc.entry insert 0 $latlcl
    .view.grd.scr.lc.ew.entry delete 0 end
    .view.grd.scr.lc.ew.entry insert 0 $ewl
    .view.grd.scr.lc.ns.entry delete 0 end
    .view.grd.scr.lc.ns.entry insert 0 $nsl
    bind .test <1> ""
    bind .test <B1-Motion> ""
    bind .test <3> ""
    setnormal
}


# Set up mouse buttons and diable screen buttons for new area identification

proc setm1 { } {
    bind .test <1> "focus .test; itemMark .test %x %y"
    bind .test <B1-Motion> "itemStroke .test %x %y"
    bind .test <3> "itemsUnderArea .test" 
    setdisabled  

}


