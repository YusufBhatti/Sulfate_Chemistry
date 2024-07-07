# Check that correct range of lan or lon has been entered 

proc checklatlon { } {

    set icode 1
    set plon	[.view.pl.scr.le.plon.entry get]
    set plat	[.view.pl.scr.le.plat.entry get]

    if { $plat < -90. || $plat > 180.} {
    set icode 0
    }
    if { $plon < 0. || $plon >360.} {
    set icode 0
    }    
    if { $icode == 0} {   
    tk_messageBox -icon error -message "Latitude or Longitude outside permitted range:\n-90< lat <180\n0< lon <360" -type ok
    }
    return $icode
}