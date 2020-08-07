

SP="test"

import sys
import numpy
import os

loc= sys.argv[-13]
batches = int(sys.argv[-12])
stringent = sys.argv[-11]
path = sys.argv[-10]
out_path = sys.argv[-9]
REF  = sys.argv[-8]
TYPE = sys.argv[-7]
ext = sys.argv[-6]
FILE = sys.argv[-5]
freq = sys.argv[-4]
SCORE = sys.argv[-3]
length = sys.argv[-2]
IDENTIFIANTS =  sys.argv[-1]

SAVE="yes"
if stringent=="yes":
	SAVE="no"

species=[SP]



print(sys.argv)


sub=[]
if FILE != "none":
	f=open(FILE ,"r")
	for l in f:
		st = l.strip("\n").strip("\r") 
		if st.endswith(ext):
			pass
		else:
			st = st + ext
		sub.append(st)
	f.close()
else:
	tmp = os.listdir(path)
	for file in tmp: 
		if file.endswith(ext) :
			sub.append(file)

sub.remove(REF)
divide = int(round(len(sub) / float(batches),0)-1)





if loc == "divide.py":
	loc = " "

partitions={}
j=0
i=0
while i < batches:
	part =  "sub"  + str(i+1)
	partitions[part] = [REF]
	partitions[part].extend(sub[j:j+divide])
	h=open(out_path + part + ".txt","w")
	for id in partitions[part]:
		h.write(id + "\n")
	h.close()
	os.system("python " + loc + "small_cruncher.py  "   + stringent + " "  + path + " " + out_path + " " + REF + " " + TYPE + " " + ext + " " + part + ".txt" + " " + freq + " " + SCORE + " " + length + " " + IDENTIFIANTS)
	i+=1
	j+=divide


















