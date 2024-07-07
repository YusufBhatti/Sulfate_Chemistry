# Draw map background, deleting previous one and updating area outline 

proc drawmap { } {
global plon plat scale xshift yshift file zoom latlon xpan ypan coast col_coast width_coast
    set xpan 0
    set ypan 0

    set code [checklatlon]
    if {$code == 0} {return}
    .test delete tg

#puts stdout "$scale $plon $plat $xshift $yshift $file"
    
    set plon	[.view.pl.scr.le.plon.entry get]
    set plat	[.view.pl.scr.le.plat.entry get]
    if {$plat < 0. } {set plon [expr 180. + $plon]}
    if {$plon > 180. } {set plon [expr -360. + $plon]}
    

#puts stdout "$scale $plon $plat $xshift $yshift $file"

#    set fid [open "|lltoeq.exec $scale $plon $plat $xshift $yshift $file" r ]
    set fid [open "|echo $scale $plon $plat $xshift $yshift $coast | ../../bin/lltoeq.exec" r ]

    set X1  [read $fid 7]
    set Y1  [read $fid 7]
    set X2  [read $fid 7]
    set Y2  [read $fid 7]
    set X3  [read $fid 7]
    set Y3  [read $fid 7]

    set c  [read $fid 1]

    while { ![eof $fid] } {
    .test create line [expr $X1] [expr $Y1] [expr $X2] [expr $Y2] [expr $X3] [expr $Y3]  -tags tg -fill $col_coast -width $width_coast

    set X1  [read $fid 7]
    set Y1  [read $fid 7]
    set X2  [read $fid 7]
    set Y2  [read $fid 7]
    set X3  [read $fid 7]
    set Y3  [read $fid 7]

    set c  [read $fid 1]

    }

    catch { close $fid } mess

   .test scale tg $xshift $yshift $zoom $zoom
    if { $latlon == 1} {
      drawlatlon
    }
    drawarea
    cross
    .test addtag tt all


}


