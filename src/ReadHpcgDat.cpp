function ReadHpcgDat(localDimensions, secondsPerRun, localProcDimensions) 
  hpcgStream = open("hpcg.dat", "r")
  readuntil(hpcgStream, '\n')
  readuntil(hpcgStream, '\n')


  for i=1:3
    localDImensions[i] =  strip(readuntil(hpcgStream, ' '))
    if localDimensions[i] < 16
      localDimensions[i] = 16
    end
  end
  readuntil(hpcgStream, '\n') # skip the rest of the second line

  secondsPerRun  = strip(readuntil(hpcgStream,' ' ))

  if secondsPerRun[1] < 0
    secondsPerRun[1] = 30 * 60 # 30 minutes
  end

  readuntil(hpcgStream, '\n') #skip the rest of the third line

  for i= 1:3
    # the user didn't specify (or values are invalid) process dimensions
    localProcDimensions[i] = strip(readuntil(hpcgStream, ' '))
    if localProcDimensions[i] < 1
      localProcDimensions[i] = 0 # value 0 means: "not specified" and it will be fixed later
    end
  close(hpcgStream)

  return 0
end
