# use with VMD something like this:
# vmd -e view-waters.tcl -args iWATs_canal.dat 4990 5000
# dibuja los frames desde 4990 hasta 4999
# (recuerda que los frames comienzan en 0)
# customize representations to taste
# info: http://www.theochem.ruhr-uni-bochum.de/~legacy.akohlmey/cpmd-vmd/part4.html

display projection orthographic
set material Opaque
set n [molinfo 0 get numframes]
set all [atomselect 0 {all}]
set Owats [atomselect 0 {water and name O}]

# elegimos el archivo y leemos la info
set first [lindex $argv 1]
set last [lindex $argv 2]
set f [open [lindex $argv 0]]
for {set i $first} {$i < $last} {incr i} {
    set atomsinsel($i) [gets $f]
}
close $f

# creamos user con la info
for {set i 0} {$i < $first} {incr i} {
    $all frame $i
    $all set user 0
}
for {set i $last} {$i < $n} {incr i} {
    $all frame $i
    $all set user 0
}
for {set i $first} {$i < $last} {incr i} {
    set insel {}
    $all frame $i
    $all set user 0
    foreach atom [$Owats get index] {
	# check if atom is in atomsinsel
	if {[lsearch -exact [split $atomsinsel($i) " "] $atom] >= 0} {
	    lappend insel 1
	} else {
	    lappend insel 0
	}
    }
    $Owats frame $i
    $Owats set user $insel
    unset insel
    #$atomsinsel($i) delete
}

# clean
$Owats delete
$all delete
unset Owats all atomsinsel i n atom first last

# selections
mol delrep 0 top

mol representation Lines 2.0
mol color Name
mol selection "backbone"
mol material $material
mol addrep top

mol representation Lines 2.0
mol color Name
mol selection "protein"
mol material $material
mol addrep top

mol representation Lines 2.0
mol color Name
mol selection {same residue as user > 0}
mol material $material
mol addrep top
mol selupdate 2 0 on
