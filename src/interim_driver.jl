using Distributed
Distributed.@everywhere include("header.jl")
Distributed.@everywhere include("main.jl")

hpcg_args = header_calling_hpcg()
ENV["JULIA_DEBUG"] = "all"
main(hpcg_args)
