
# Return current screen coords and convert to new polar values

proc ItemGet { c x y} {
    global scale xshift yshift plon plat zoom xpan ypan
    set areaX1 [$c canvasx $x]
    set areaY1 [$c canvasy $y]

    set rlon [ expr (  ($areaX1 -$xshift-$xpan)  / ($scale*$zoom))] 
    set rlat [ expr  (  ($yshift+$ypan - $areaY1) / ($scale*$zoom))] 
#puts stdout "$plon $plat $rlon $rlat"
#    set fid [open "|eqtoll.exec $plon $plat $rlon $rlat" r ]
    set fid [open "|echo $plon $plat $rlon $rlat | ../../bin/eqtoll.exec" r ]
    set rlon  [read $fid 8]
    set rlat  [read $fid 8]
#puts stdout "$rlon $rlat"
    set rlon [ format %.2f [ expr (  $rlon - 180.)] ]
    set rlat [ format %.2f [expr  (  90. -  $rlat)] ]
#    if {$rlat > 90. } {set rlat [expr -180. + $rlat]}
       
    
    if { $rlon < -180.} {
       set rlon [ expr (360. + $rlon) ]
    }
    
    if {$rlon < 0. } {set rlon [expr 360. + $rlon]}    
    
    .view.pl.scr.le.plon.entry delete 0 end
    .view.pl.scr.le.plon.entry insert 0 $rlon
    .view.pl.scr.le.plat.entry delete 0 end
    .view.pl.scr.le.plat.entry insert 0 $rlat
}


# Set up mouse buttons and disable screen buttons for new centre of view

proc setm2 { } {
    bind .test <1> "focus .test; ItemGet .test %x %y; cross_at_pointer .test %x %y "
    bind .test <3> "bind .test <3> \"\";bind .test <1> \"\"; setnormal "
    setdisabled
}

