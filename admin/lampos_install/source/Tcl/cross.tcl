# Display cross
proc cross { } {
global xshift yshift col_rect
.test delete tc
.test create line  [ expr $xshift - 4 ] $yshift [expr $xshift +4] $yshift -width 2 -tags tc  -fill $col_rect

.test create line  $xshift [ expr $yshift -4 ]  $xshift [ expr $yshift + 4 ] -width 2 -tags tc  -fill $col_rect
.test addtag tt all
}
proc cross_at_pointer {c x y } {
    $c delete tc
    set areaX1 [$c canvasx $x]
    set areaY1 [$c canvasy $y]
$c create line  [ expr $areaX1 - 4 ] $areaY1 [expr $areaX1 +4] $areaY1 -width 2 -tags tc -fill red

$c create line  $areaX1 [ expr $areaY1 -4 ]  $areaX1 [ expr $areaY1 + 4 ] -width 2 -tags tc -fill red

.test addtag tt all
}
