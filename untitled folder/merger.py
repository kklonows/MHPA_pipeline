import sys
filenames = sys.argv[1:]
output = open("mergez.txt", "w")
files = []
for arg in sys.argv[1:]:
    files.append(open(arg, "r"))
num_files = len(files)
num_empty = 0
while True:
    num_empty = 0
    line = []
    for file in files:
        item = file.readline() # returns empty string after EOF
        if not item: 
            item = "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" # or other marker value
            num_empty+=1
        line.append(item)
    if num_empty == num_files:
        break
    output.write("\t".join([x.strip() for x in line]))
    output.write("\n")
for file in files:
    file.close()
output.close()
