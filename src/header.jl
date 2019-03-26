using ArgParse
include("main.jl")
	
	s = ArgParseSettings()
	@add_arg_table s begin 
		"--np"
			help = "np"
			nargs = '?'
			arg_type = Int
		"--nx"              
			help = "nx"
                        nargs = '?'
                        arg_type = Int
		"--ny" 
			help = "ny"
                        nargs = '?'
                        arg_type = Int
		"--nz" 
			help = "nz"
                        nargs = '?'
                        arg_type = Int
		"--DEFAULT" 
			help = "default"
                        action = :store_true
		"--rt" 
			help = "rt"
                        nargs = '?'
                        arg_type = Int
		"--USE_MPI"
			help = "mpi"
                        action = :store_true
	end
	parsed_args =parse_args(s)
	println(parsed_args)
	main(parsed_args)
