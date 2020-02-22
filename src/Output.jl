include("OutputFile_header.jl")
using DataStructures
using Dates
function OutputFile(of, name_arg, version_arg)
	of.name = name_arg
	of.version = version_arg
	of.edol = "\n"
	of.keySeperator = "::"
end


function OutputFile(of) 
	of.eol = "\n"
	of.keySeperator = "::"
end

function del_OutputFile(of) 
  for it in of.descendents
    it = nothing
  end
end

 # descendents is an array of strings
function OutputFile_add(of, key_arg, value_arg::String) 
    of = OutputFile 
    of.key = key_arg
    of.value = value_arg
end


function OutputFile_add(of, key_arg, value_arg::Float64) 
    of = OutputFile
    of.key = key_arg
    of.value = string(value_arg)
end


function OutputFile_add(of,key_arg, value_arg::Int64) 
    of = OutputFile
    of.key = key_arg
    of.value = string(value_arg)
end


function OutputFile_setKeyValue(of,key_arg, value_arg) 
  of.key = key_arg
  of.value = value_arg
end


function OutputFile_get(of, key_arg) 
  for kv in of.descendents
    if k == key_arg
      return v
    end
  end
  return 0
end


function generateRecursive(of,prefix) 
  result = ""

  result = result *prefix * of.key * "=" * of.value * of.eol

  for kv in of.descendents
    result = result + generateRecursive(of, prefix * of.key * of.keySeparator)
  end

  return result
end


function generate(of) 
  result = of.name * "\nversion=" * version * eol

  result = result * "LANG=julia" * eol

  for kv in of.descendents 
    result = result * generateRecursive(of, "")
  end
  #time this using mytimer()
  #Check ReadHPCGDat.jl
  date = now()
  sdate =string(date)

  filename = name *"_" *version * "_"
  filename += string(sdate) + ".exp"

  ofstream myfile(filename.c_str())
  myfile = fopen(filename, "w")
  write(myfile, result)
  close(myfile)

  return result
end

function allocKeyVal(key_arg, value_arg) 
  of = OutputFile()
  setKeyValue(of, key_arg, value_arg)
  return of
end
