proc drawlatlon { } {
global plon plat scale xshift yshift filelat zoom latlon xpan ypan col_lat width_lat
    .test delete tl
if { $latlon == 1} {
#puts stdout "$scale $plon $plat $xshift $yshift $filelat"

    set code [checklatlon]
    if {$code == 0} {return}
    set plon	[.view.pl.scr.le.plon.entry get]
    set plat	[.view.pl.scr.le.plat.entry get]
    if {$plat < 0. } {set plat [expr 180. + $plat]}
    if {$plon > 180. } {set plon [expr -360. + $plon]}    

#puts stdout "$scale $plon $plat $xshift $yshift $filelat"


    set fid [open "|echo $scale $plon $plat $xshift $yshift 2 | ../../bin/lltoeq.exec" r ]
    set X1  [read $fid 7]
    set Y1  [read $fid 7]
    set X2  [read $fid 7]
    set Y2  [read $fid 7]
    set X3  [read $fid 7]
    set Y3  [read $fid 7]


    set c  [read $fid 1]
    while { ![eof $fid] } {
#puts stdout "$X1 $Y1 $X2 $Y2 $xpan $ypan $zoom"
    .test create line [ expr $X1+$xpan/$zoom] [ expr $Y1+$ypan/$zoom] [ expr $X2+$xpan/$zoom] [ expr $Y2+$ypan/$zoom] [ expr $X3+$xpan/$zoom] [ expr $Y3+$ypan/$zoom] -tags tl -fill $col_lat -width $width_lat
    set X1  [read $fid 7]
    set Y1  [read $fid 7]
    set X2  [read $fid 7]
    set Y2  [read $fid 7]
    set X3  [read $fid 7]
    set Y3  [read $fid 7]


    set c  [read $fid 1]

    }
    .test scale tl $xshift $yshift $zoom $zoom
    catch { close $fid } mess

    }
    .test addtag tt all
}



