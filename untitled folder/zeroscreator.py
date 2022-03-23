import sys
newfile=open(sys.argv[1], 'w')

for k in range(0,(int)(sys.argv[2])):
	for j in range(0,(int)(sys.argv[3])):
		newfile.write("0\t")
	newfile.write("\n")
