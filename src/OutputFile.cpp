include("OutputFile_header.jl")
using DataStructures
function OutputFile(of, name_arg, version_arg)

  #: name(name_arg), version(version_arg), eol("\n"), keySeparator("::") {}
end

function OutputFile(of, name_arg, version_arg)
end

function OutputFile(void) 
	#: eol("\n"), keySeparator("::") {}
end

#does this delete all the OutputFile structs listed in descendents
function del_OutputFile(of) 
  for it in of.descendents
    it = nothing
  end
end


function OutputFile_add(of, key_arg, value_arg::String) 
#  descendants.push_back(allocKeyVal(key_arg, value_arg));
end


function OutputFile_add(of, key_arg, value_arg::Float64) 
#  stringstream ss;
#  ss << value_arg;
#  descendants.push_back(allocKeyVal(key_arg, ss.str()));
end


function OutputFile_add(of,key_arg, value_arg::Int64) 
   stringstream ss;
#  ss << value_arg;
#  descendants.push_back(allocKeyVal(key_arg, ss.str()));
end


function OutputFile_setKeyValue(of,key_arg, value_arg) 
  of.key = key_arg
  of.value = value_arg
end


function OutputFile_get(of, key_arg) 
  for it in of.descendents
    if it.key == key_arg
      return it
    end
  end
  return 0
end


function OutputFile_generateRecursive(of,prefix) 
  string result = ""

  result += prefix + key + "=" + value + eol

 for it in of.descendents
    result += generateRecursive(of, prefix + key + keySeparator);
  end

  return result;
end


function OutputFile_generate(of) 
  string result = name + "\nversion=" + version + eol;

  for it in of.descendent 
    result += generateRecursive(it, "")
  end
  #time this using mytimer()
  time_t rawtime;
  time(&rawtime);
  tm * ptm = localtime(&rawtime);
  #Check ReadHPCGDat.jl
  char sdate[25];
  #use tm_mon+1 because tm_mon is 0 .. 11 instead of 1 .. 12
  # use Date Time package
  write(sdate,ptm.tm_year + 1900, ptm.tm_mon+1, ptm.tm_mday, ptm.tm_hour, ptm.tm_min,ptm.tm_sec);

  filename = name *"_" *version * "_";
  filename += string(sdate) + ".txt";

  ofstream myfile(filename.c_str());
  myfile << result;
  myfile.close();

  return result;
end

function OutputFile_allocKeyVal(key_arg, value_arg) 
  of = OutputFile()
  setKeyValue(of, key_arg, value_arg)
  return of;
}
