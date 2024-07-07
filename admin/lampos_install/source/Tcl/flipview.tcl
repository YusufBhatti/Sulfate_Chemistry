
# Draw equivalent domain (alternative grid rotation)
# Author: Hilary Oliver (2000), Stuart Moore (2007)  NIWA

proc flipview { } {

global lonlc latlc plon plat 


set lonlc  [ getround .view.grd.scr.le.lonlc $lonlc 0. 360. 7]
set lonlc  [ expr 180 + $lonlc ]


if { $lonlc > 360 } { set lonlc [ expr $lonlc - 360 ] } 

setarea
drawarea

    set plon   [.view.pl.scr.le.plon.entry get]
    set plat   [.view.pl.scr.le.plat.entry get]

    set plon [ expr $plon - 180 ]
    set plat [ expr 180 - $plat ]

    setpole
    drawmap

}

