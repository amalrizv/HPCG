node_names      =["172.16.2.102"]
nprocs          = 2
addprocs(1)
for i=1:length(node_names)
        addprocs([(nodenames[i], nprocs)])
end










#= Generalized function for adding n procs on each hostname_list which should be an array
function add_n_procs_to_hostnames(hostname_list, procs_on_each_node)

for i=1:length(hostname_list)
        addprocs([(hostname_list[i], procs_on_each_node)])
end

=#
~                                                                                                                      

