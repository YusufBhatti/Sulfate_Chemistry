# UK LAM parameters


# NAE regional model parameters

proc setnae { } {
new_dynamics_mode 1
global dx;  set dx .11
global dy;  set dy .11
global lonlc; set lonlc 326.22 
global latlc; set latlc -20.070
global ew; set ew  600
global ns; set ns  300
new_dynamics_mode 1
    setarea
}


proc setuk4 { } {
new_dynamics_mode 1
global dx;  set dx .036
global dy;  set dy .036
global lonlc; set lonlc 354.3
global latlc; set latlc -4.452
global ew; set ew  288
global ns; set ns  360
new_dynamics_mode 1
    setarea
}


# Polar rotation NAE regional model

proc setnaep { } {
global plon; set plon 177.5
global plat; set plat 37.5
    setpole
}

# Polar rotation UK4 model

proc setuk4p { } {
global plon; set plon 177.5
global plat; set plat 37.5
    setpole
}



# Polar rotation for standard lat-lon

proc setstanp { } {
global plon; set plon 180.
global plat; set plat 90.
    setpole
}


