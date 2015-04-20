from numpy import *
filepath = 'C:\\Users\\Nath\\Nathaniel\\Auburn\\Research\\Jamming and Spoofing Detection\\Code\\Data\\USRP Converter\\Converted Data\\'
name     = 'April_9_5M_Sim03042015_int16_Simulator_Ant2_NoInterference_sat109_GSPDOLink.dat'
filename = filepath + name
blksize = int(5000000*.1)
processLength = int(100/.1)

cosv4 = matrix([ 1, 0, -1, 0])
sinv4 = matrix([ 0, 1, 0, -1])
repmat = matrix([1] * int(blksize/(size(cosv4))))

si = transpose(transpose(sinv4) * (repmat))
co = transpose(transpose(cosv4) * (repmat))
modsine   = array(si.reshape(1,size(si)))
modcosine = array(co.reshape(1,size(co)))

fi = open(filename,"r")
totmax = zeros((processLength,1),int16)
myarray = zeros((1,2*blksize))
print("Determining Maximum Value")
for counter in range(0,processLength):
	print(int(counter/processLength*100)," percent complete         \r",end = '')
	myarray = fromfile(fi,dtype=int16,count = 2*blksize)
	rawSigI = myarray[0:2*blksize:2]
	rawSigQ = myarray[1:2*blksize:2]
	
	outputWaveform = ((rawSigI) * (modcosine) - (rawSigQ) * modsine)
	totmax[counter] = amax(outputWaveform)
print(100)
maximum = amax(totmax)

fi.close()	

fi = open(filename,"r")
cleared = open('outFile.dat','w')
cleared.close()
fo = open('outFile.dat',"a")
print("Processing Data...")
for counter in range(0,processLength):	
	print(int(counter/processLength*100)," percent complete         \r",end = '')
	myarray = fromfile(fi,dtype=int16,count = 2*blksize)
	rawSigI = myarray[0:2*blksize:2]
	rawSigQ = myarray[1:2*blksize:2]

	
	outputWaveform = int8(127 * ((rawSigI) * (modcosine) - (rawSigQ) * modsine) / maximum)
	outputWaveform.tofile(fo)
		
	
	
print(100)

fi.close()	
fo.close()

