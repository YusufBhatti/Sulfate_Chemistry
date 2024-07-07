# Read and write options file

proc readOptions { } {
global home col_lat col_coast col_back latlon col_rect
global width_lat width_coast width_rect

if { [file exists ~$home/.lampos.options] == 1 } {
set f [open  ~$home/.lampos.options r ]
gets $f width_lat
gets $f width_rect
gets $f width_coast
gets $f col_lat
gets $f col_coast
gets $f col_back
gets $f col_rect
gets $f latlon
close $f
}
}

proc writeOptions { } {
global home col_lat col_coast col_back latlon col_rect
global width_lat width_coast width_rect


set f [open  ~$home/.lampos.options w ]
puts $f $width_lat
puts $f $width_rect
puts $f $width_coast
puts $f $col_lat
puts $f $col_coast
puts $f $col_back
puts $f $col_rect
puts $f $latlon
close $f
#destroy .

}


