#!/bin/sh
#=
exec julia -O3 "$0" -- $@
=#


# Sustituye un residuo por otro


using Printf


include("cloner.in")
Ntotal = Nrings*Natoms # Number of atoms in the system (in rings)
# Nrings es el número total de anillos en el sistema a limpiar
# Nres es el número de residuos en un anillo, tras ellos ponemos un TER
# Natoms es el número de átomos en un anillo

function clon(file, names, resnames, resids, pos, oldresname, subspdb)
    subsnames = Vector{String}()
    subsresnames = Vector{String}()
    subsresids = Vector{Int64}()
    subspos = Vector{Float64}()
    open(subspdb) do subsfile 
        regex = r"\h([^\W\d]+\d*-?)\h+(\w+-?)\h+[^\W\d]*\h*(\w+)\h+((?:\h+-?\d*\.?\d*){3})"
        resid = 0
        previous_resid = ""
        for line in eachline(subsfile)
            m = match(regex, line)
            if !isnothing(m)
                push!(subsnames, m.captures[1])
                push!(subsresnames, m.captures[2])
                current_resid = m.captures[3]
                if current_resid != previous_resid
                    previous_resid = current_resid
                    resid += 1
                end
                push!(subsresids, resid)
                subspos = vcat(subspos, parse.(Float64, split(m.captures[4])))
            end
        end
    end
    
    subsatoms = 0
    lastsubs = 0
    push!(resids, 0) # dummy para escribir los TER
    for index in 1:length(names)
        atid = index + subsatoms
        resid = resids[index];
	if resnames[index] == oldresname
	    if resid == lastsubs
		subsatoms -= 1
		continue
	    end
	    lastsubs = resid
	    subsatoms += length(subsnames)-1
	    for subsindex in 1:length(subsnames)
		resid += subsresids[subsindex]-1
		print(file, "ATOM", repeat(" ", 7-length(string(atid))))
		print(file, atid, "  ", subsnames[subsindex], repeat(" ", 4-length(subsnames[subsindex])))
		print(file, subsresnames[subsindex], repeat(" ", 6-length(string(resid))))
		print(file, resid, "    ")
		newpos = pos[3*(index-1) .+ (1:3)] + subspos[3*(subsindex-1) .+ (1:3)]
		for dim in 1:3
		    print(file, repeat(" ", 8-length(@sprintf("%.3f", newpos[dim]))))
                    @printf(file, "%.3f", newpos[dim])
		end
		println(file, "  0.00  0.00")
		atid += 1
	    end
	else
	    print(file, "ATOM", repeat(" ", 7-length(string(atid))))
            print(file, atid, "  ", names[index], repeat(" ", 4-length(names[index])))
            print(file, resnames[index], repeat(" ", 6-length(string(resid))))
            print(file, resid, "    ")
            newpos = pos[3*(index-1) .+ (1:3)]
            for dim in 1:3
                print(file, repeat(" ", 8-length(@sprintf("%.3f", newpos[dim]))))
		@printf(file, "%.3f", newpos[dim])
	    end
	    println(file, "  0.00  0.00")
	end
	
	if index <= Ntotal
	    if index%Natoms == 0
		println(file, "TER   ")
	    end
	else
	    if resid != resids[index+1]
		println(file, "TER   ")
	    end
	end
    end
end


function main(oldresname, subspdb)
    # Initialize
    names = Vector{String}()
    resnames = Vector{String}()
    resids = Vector{Int64}()
    pos = Vector{Float64}()

    # Read the pdb information
    open(pdbfile) do file
        regex = r"\h([^\W\d]+\d*-?)\h+(\w+-?)\h+[^\W\d]*\h*(\w+)\h+((?:\h+-?\d*\.?\d*){3})"
	resid = 0
	previous_resid = ""
	for line in eachline(file)
	    m = match(regex, line)
	    if !isnothing(m)
	        push!(names, m.captures[1])
	        push!(resnames, m.captures[2])
		current_resid = m.captures[3]
		if current_resid != previous_resid
		    previous_resid = current_resid
		    resid += 1
		end
	        push!(resids, resid)
	        pos = vcat(pos, parse.(Float64, split(m.captures[4])))
	    end
	end
    end

    # Write the new pdb with the cloned system (traslated)
    open(newfile, "w") do file
        print(file, "\n")
	clon(file, names, resnames, resids, pos, oldresname, subspdb)
	println(file, "END   ")
    end
end


main(ARGS[1], ARGS[2])
