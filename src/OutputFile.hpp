#=
 @file Output_File.hpp

 HPCG output file classes
=#

using DataStructures

#The OutputFile class for the uniform collecting and reporting of performance data for HPCG

#=

  The OutputFile class facilitates easy collecting and reporting of
  key-value-formatted data that can be then registered with the HPCG results
  collection website. The keys may have hierarchy key1::key2::key3=val with
  double colon :: as a separator. A sample output may look like this (note how
  "major" and "micro" keys repeat with different ancestor keys):

\code

version=3.2.1alpha
version::major=3
version::minor=2
version::micro=1
version::release=alpha
axis=xyz
axis::major=x
axis::minor=y

\endcode
=#
mutable struct OutputFile
  descendants::Any  	# This is a linked list of descendant elements of the type OutputFile
  name::String 		# name of the benchmark
  version::String 	# version of the benchmark
  key::String 		# the key under which the element is stored
  value::String 	# the value of the stored element
  eol::String 		# end-of-line character sequence in the output file
  keySeparator::String 	# character sequence to separate keys in the output file
  			# Recursively generate output string from descendant list, and their descendants and so on
end
