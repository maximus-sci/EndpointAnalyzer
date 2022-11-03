#as of version 1.0, this script now works properly with EndpointAnalyzer,v1
#this takes normalized bedgraph files and then outputs that read depth relative to control
import os
import math
import sys


bgFiles = [f for f in os.listdir(os.getcwd()) if ((f.endswith(".bg")) & ("norm" in f) & ("rel" not in f))]
if (len(bgFiles) > 0):
	print("List of normalized files:")
	for i in range(len(bgFiles)):
		print(str(i+1)+") "+bgFiles[i])
else:
	print("No normalized files found -- make sure you run the normalizer first")

missingCntrlName = True
while missingCntrlName:
	try:
		cntrl_name = input("Enter the name or number of the control file (or type 'quit'):")
		try:
			confirm = False
			while (confirm == False):
				print("\nYou selected:",bgFiles[int(cntrl_name)-1])
				correctFile = input("if this is correct, type yes:")
				if (correctFile.lower() in ["yes","y"]):
					confirm = True
				else:
					print("please restart")
					sys.exit(0)
			cntrl_name = bgFiles[int(cntrl_name)-1]
			cntrl      = open(cntrl_name,"r")
		except:
			continue
		print("Control file found!")
		print("\nUsing",cntrl_name,"to normalize remaining files\n\n")
		missingCntrlName = False
		cntDic     = {}
	except:
		if cntrl_name == 'quit':
			sys.exit(0)
		print("File",cntrl_name,"not found, please try again")
		continue
	

for line in cntrl:
    lineParts        = line.strip().split("\t")
    cntDic[lineParts[1]] = float(lineParts[3])

bgFiles.remove(cntrl_name)
print("List of files corrected relative to " +cntrl_name+":\n"+ "\n".join(bgFiles))
print("****************************************")
for f in bgFiles:
    print("Processing "+f+" relative to control")
    inFile  = open(f,"r")
    name    = f.split(".")[0]
    #L2File  = open("Relative_L2FC_"+name,"w+", newline = "\n")
    subFile = open(name+".relCntrlSubtract.norm.filt.bg","w+", newline = "\n")
    for line in inFile:
        line = line.rstrip()
        lineParts = line.split("\t")
        newLine = lineParts[0:3]
        #print(float(lineParts[3])," divided by ",float(cntDic[lineParts[1]]))
        #log2    = newLine + [str(math.log(float(lineParts[3])/float(cntDic[lineParts[1]]),2))]
        subtract= newLine + [str(float(lineParts[3])-float(cntDic[lineParts[1]]))]
        #L2File.write("\t".join(log2)+"\n")
        subFile.write("\t".join(subtract)+"\n")
    inFile.close()
    #L2File.close()
    subFile.close()  