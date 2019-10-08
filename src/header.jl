using ArgParse
using Logging
using Distributed
Distributed.@everywhere include("main.jl")

#TODO: USE_MPI is redundant for np
s = ArgParseSettings()
@add_arg_table s begin 
    "--np"
        help     = "np"
        nargs    = '?'
        arg_type = Int
        default  = 1
    "--nx"              
        help     = "nx"
        nargs    = '?'
        arg_type = Int
    "--ny" 
        help     = "ny"
        nargs    = '?'
        arg_type = Int
    "--nz" 
        help     = "nz"
        nargs    = '?'
        arg_type = Int
    "--npx"              
        help     = "npx"
        nargs    = '?'
        arg_type = Int
    "--npy" 
        help     = "npy"
        nargs    = '?'
        arg_type = Int
    "--npz" 
        help     = "npz"
        nargs    = '?'
        arg_type = Int
    "--zl"              
        help     = "zl"
        nargs    = '?'
        arg_type = Int
    "--zu" 
        help     = "zu"
        nargs    = '?'
        arg_type = Int
    "--pz" 
        help     = "pz"
        nargs    = '?'
        arg_type = Int
    "--DEFAULT" 
        help   = "default"
        action = :store_true
    "--rt" 
        help     = "rt"
        nargs    = '?'
        arg_type = Int
        default  = 10
end

parsed_args = parse_args(s)

@info "Arguments" parsed_args

ENV["JULIA_DEBUG"] = "all"

@everywhere main(parsed_args)
