# use with VMD something like this:
# vmd -e view-qmmm.tcl
# customize representations to taste

set material Opaque

# pasamos el nombre del archivo como input
#set fm1 [open [lindex $argv 0]]
# elegimos el archivo automaticamente
set fm1 [open qmmask.dat]
set m1 "[read $fm1]"
close $fm1
set m1nh "$m1 and not hydrogen"

#display rendermode GLSL
#axes location off
display projection orthographic
#display depthcue off
#color Display Background white
#display projection perspective

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

mol representation VDW 0.2 12.0
mol color Name
mol selection $m1
mol material $material
mol addrep top

mol representation Licorice 0.15
mol color Name
mol selection "backbone"
mol material $material
mol addrep top

mol representation Lines 0.15
mol color Name
mol selection "water"
mol material $material
mol addrep top

mol representation Lines 0.15
mol color Name
mol selection "resname TFA"
mol material $material
mol addrep top
