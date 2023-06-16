.PHONY: null min heat equil prod mmpbsa
.DEFAULT: null

null:
	@:

min:
ifeq ($(INPUT),)
	$(error Error: La variable INPUT no está declarada.\
	Uso: "make min INPUT=job")
endif
	cd min/; sbatch -J ${INPUT}_min run_min.sh ${INPUT}

heat:
ifeq ($(INPUT),)
	$(error Error: La variable INPUT no está declarada.\
	Uso: "make heat INPUT=job")
endif
	cd heat/; sbatch -J ${INPUT}_heat run_heat.sh ${INPUT}

equil:
ifeq ($(INPUT),)
	$(error Error: La variable INPUT no está declarada.\
	Uso: "make equil INPUT=job GPU=index")
endif
ifeq ($(GPU),)
	$(error Error: La variable GPU no está declarada.\
	Uso: "make equil INPUT=job GPU=index")
endif
	cd equil/; sbatch -J ${INPUT}_equil run_equil.sh ${INPUT} ${GPU}

prod:
ifeq ($(INPUT),)
	$(error Error: La variable INPUT no está declarada.\
	Uso: "make prod INPUT=job GPU=index")
endif
ifeq ($(GPU),)
	$(error Error: La variable GPU no está declarada.\
	Uso: "make prod INPUT=job GPU=index")
endif
	cd prod/; sbatch -J ${INPUT}_prod run_prod.sh ${INPUT} ${GPU}

mmpbsa:
ifeq ($(INPUT),)
	$(error Error: La variable INPUT no está declarada.\
	Uso: "make mmpbsa INPUT=job")
endif
	cd mmpbsa; mkdir output
	cd mmpbsa; cp run_mmpbsa_mpi.sh MMPBSA_MPI.sh output/
	cd mmpbsa; cp complex.top receptor.top ligand.top output/
	cd mmpbsa; cp ${INPUT}_mdestab_extract.coord output/
	cd mmpbsa/output/; sbatch -J ${INPUT}_mmpbsa run_mmpbsa_mpi.sh ${INPUT}
