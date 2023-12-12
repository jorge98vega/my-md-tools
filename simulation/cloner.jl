#!/bin/sh
#=
exec julia -O3 "$0" -- $@
=#


using Printf


include("cloner.in")
Ntotal = Nrings*Natoms # Number of atoms in the system


function rotate(pos, center, alpha)
    rotmat = [cos(alpha) -sin(alpha) 0
              sin(alpha) cos(alpha)  0
              0          0           1]
    return rotmat*(pos - center) + center
end


function clon(file, names, resnames, resids, pos, com, t, alpha, n)
    for ring in 0:Nrings-1
        for atom in 1:Natoms
	    index = ring*Natoms + atom
	    atid = (n-1)*Ntotal + index
	    resid = (n-1)*Nrings*Nres + resids[index];
	    print(file, "ATOM", repeat(" ", 7-length(string(atid))))
	    print(file, atid, "  ", names[index], repeat(" ", 4-length(names[index])))
	    print(file, resnames[index], repeat(" ", 6-length(string(resid))))
	    print(file, resid, "    ")
	    rotpos = rotate(pos[3*(index-1) .+ (1:3)], com, alpha)
	    newpos = rotpos + t
	    # newpos = pos[3*(index-1) .+ (1:3)] + t
	    for dim in 1:3
	        print(file, repeat(" ", 8-length(@sprintf("%.3f", newpos[dim]))))
	        @printf(file, "%.3f", newpos[dim])
	    end
	    println(file, "  0.00  0.00")
	end
	println(file, "TER   ")
    end
end


function main()
    # Initialize
    names = Vector{String}(undef, Ntotal)
    resnames = Vector{String}(undef, Ntotal)
    resids = zeros(Int64, Ntotal)
    pos = zeros(Float64, 3*Ntotal)
    com = zeros(Float64, 3)

    # Read the pdb information
    open(pdbfile) do file
        regex = r"\h([^\W\d]+\d*-?)\h+(\w+-?)\h+[^\W\d]*\h*(\w+)\h+((?:\h+-?\d*\.?\d*){3})"
	i = 1
	for line in eachline(file)
	    m = match(regex, line)
	    if !isnothing(m)
	        names[i] = m.captures[1]
	        resnames[i] = m.captures[2]
	        resids[i] = parse(Int64, m.captures[3])
	        pos[3*(i-1) .+ (1:3)] = parse.(Float64, split(m.captures[4]))
		com[1:3] += pos[3*(i-1) .+ (1:3)]
		i += 1
		i > Ntotal && break
	    end
	end
    end
    com[1:3] = com[1:3]/Ntotal

    # Write the new pdb with the cloned system (traslated)
    open(newfile, "w") do file
        print(file, "\n")
	for n in 1:size(ts)[1]
	    clon(file, names, resnames, resids, pos, com,
		 ts[n, :], as[n]*pi/180, n)
	end
	println(file, "END   ")
    end
end


main()
