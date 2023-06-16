# use with VMD something like this:
# vmd -e view-qmmm.tcl
# customize representations to taste

set material Opaque

# para elegir el archivo
#set fm1 [open [lindex $argv 0]]
# elegimos el archivo automaticamente
set fm1 [open qmmask.dat]
set m1 "serial [read $fm1]"
close $fm1
#set m1 [lindex $argv 0]
#set m2 [lindex $argv 1]
#set m3 [lindex $argv 2]
set m1nh "$m1 and not hydrogen"
#set m4 "$m1 or $m2 or $m3"

#display rendermode GLSL
#axes location off
display projection orthographic
#display depthcue off
#color Display Background white
#display projection perspective

#set fn [lindex $argv 3]
#mol new $fn type parm7
mol delrep 0 top

mol representation DynamicBonds 1.7 0.2 12.0
mol color Name
mol selection $m1nh
mol material $material
mol addrep top

mol representation DynamicBonds 1.3 0.2 12.0
mol color Name
mol selection $m1
mol material $material
mol addrep top

mol representation Lines 2.0
mol color Name
mol selection $m1
mol material $material
mol addrep top

#mol representation DynamicBonds 2.0 0.20 50.0
#mol color Name
#mol selection $m4
#mol material Transparent
#mol addrep top

mol representation VDW 0.2 12.0
mol color Name
mol selection $m1
mol material $material
mol addrep top

#adatom
#mol representation VDW 1.0 50.0
#mol color colorid 4
#mol selection $m2
#mol material $material
#mol addrep top

#surface
#mol representation VDW 1.0 50.0
#mol color colorid 8
#mol selection $m3
#mol material $material
#mol addrep top

#foreach filename [lrange $argv 4 $argc-1] {
#mol addfile $filename type netcdf waitfor all
#}
