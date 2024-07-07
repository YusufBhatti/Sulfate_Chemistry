
proc labelentry {path text length} {
    frame $path
    label $path.label -text $text
    entry $path.entry -width $length -relief sunken
    pack $path.label -side left -expand y
    pack $path.entry -side right -expand y
}


# Input variable, constraining to limits and length

proc getround { w va min max len } {
    set va	[$w.entry get]
    if { $va < $min } {
       set va $min
       $w.entry delete 0 end
       $w.entry insert 0 $va
    }
    if { $va > $max } {
       set va $max
       $w.entry delete 0 end
       $w.entry insert 0 $va
    }
    if { [ string length $va ] > $len } {
       set va [string range $va 0 [expr $len -1 ] ]
       $w.entry delete 0 end
       $w.entry insert 0 $va
    }
    return $va
}


# Reset boxes containing LAM parameters

proc setarea { } {
global dx dy
global lonlc
global latlc
global ew ns dt levs
    set rlonlc $lonlc
    if {$rlonlc < 0. } {set rlonlc [expr 360. + $rlonlc]}    
    .view.grd.scr.lw.dx.entry delete 0 end
    .view.grd.scr.lw.dx.entry insert 0 $dx
    .view.grd.scr.lw.dy.entry delete 0 end
    .view.grd.scr.lw.dy.entry insert 0 $dy
    .view.grd.scr.le.lonlc.entry delete 0 end
    .view.grd.scr.le.lonlc.entry insert 0 $rlonlc
    .view.grd.scr.le.latlc.entry delete 0 end
    .view.grd.scr.le.latlc.entry insert 0 $latlc
    .view.grd.scr.lc.ew.entry delete 0 end
    .view.grd.scr.lc.ew.entry insert 0 $ew
    .view.grd.scr.lc.ns.entry delete 0 end
    .view.grd.scr.lc.ns.entry insert 0 $ns
}

# Reset boxes containing Polar coordinate values

proc setpole { } {
global plon plat 
    set rlon $plon
    set rlat $plat
# allow pole latitude > 90 deg
#    if {$rlat > 90. } {set rlat [expr -180. + $rlat]}
    if {$rlon < 0. } {set rlon [expr 360. + $rlon]}
    .view.pl.scr.le.plon.entry delete 0 end
    .view.pl.scr.le.plon.entry insert 0 $rlon
    .view.pl.scr.le.plat.entry delete 0 end
    .view.pl.scr.le.plat.entry insert 0 $rlat
}

# Disable & reset screen buttons

proc setnormal { } {
#    .mou.area configure -state normal
    .view.grd.areas configure -state normal
#    .mou.view configure -state normal
     .view.pl.poles configure -state normal
    .mou.zoomin configure -state normal
    .mou.zoomou configure -state normal
#    .mou.print configure -state normal
#    .mou.grid configure -state normal
    .view.grd.apply configure -state normal
    .view.pl.scr.lw.apply configure -state normal
}
proc setdisabled { } {
#    .mou.area configure -state disabled
    .view.grd.areas configure -state disabled
#    .mou.view configure -state disabled
    .view.pl.poles configure -state disabled
    .mou.zoomin configure -state disabled
    .mou.zoomou configure -state disabled
#    .mou.print configure -state disabled
#    .mou.grid configure -state disabled
    .view.grd.apply configure -state disabled
    .view.pl.scr.lw.apply configure -state disabled
}


