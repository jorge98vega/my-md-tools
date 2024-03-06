# use with VMD something like this:
# vmd -e figure-background.tcl -args 'molecules' parm7 netcdf... 
# customize representations to taste

set materialfore AOEdgy
set materialback AOChalky

set selectionbonded [lindex $argv 0]

axes location off
display rendermode GLSL
display projection orthographic
display depthcue off
color Display Background white

set parmfile [lindex $argv 1]
mol new $parmfile type parm7
mol delrep 0 top

mol representation DynamicBonds 1.6 0.2 50.0
mol color Name
mol selection $selectionbonded and (not name H C N O or not protein)
mol material $materialfore
mol addrep top

mol representation VDW 0.2 50.0
mol color Name
mol selection $selectionbonded and (not name H C N O or not protein)
mol material $materialfore
mol addrep top

mol representation HBonds 3.2 16.0 8.0
mol color Name
mol selection $selectionbonded and name NZ HZ1 HZ2 HZ3 O H1 H2 OD1 OD2 OH HH
mol material $materialfore
mol addrep top

mol representation Licorice 0.2 50.0
mol color Name
mol selection (backbone or name H)
mol material $materialback
mol addrep top

mol representation Licorice 0.1 50.0
mol color Name
mol selection protein
mol material $materialback
mol addrep top

color Name C gray
color Name F green
color Name H silver

color change rgb gray 0.450 0.450 0.450
color change rgb silver 0.900 0.900 0.900

material change Outline $materialfore 1.7

material change Ambient $materialback 0.65
material change Diffuse $materialback 0.35

foreach trajfile [lrange $argv 2 $argc-1] {
mol addfile $trajfile type netcdf waitfor all
}
