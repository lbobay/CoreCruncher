
import sys
import os


REF = sys.argv[-3]
path = sys.argv[-2]
out_path = sys.argv[-1]


files=os.listdir(out_path)

partitions=[]
for stuff in files:
	if stuff.startswith("sub") and stuff.endswith(".txt"):
		partitions.append(stuff)

genomes={}
for sub in partitions:
	genomes[sub] =[]
	f=open(out_path + sub ,"r")
	for l in f:
		genomes[sub].append(l.strip("\n"))
	f.close()

compteur={}
build={}
parent={}
sub1=partitions[0]
f=open(out_path + "families_" + sub1 ,"r")
for l in f:
	a=l.strip("\n").split("\t")
	fam=a[0]
	build[fam] = a
	compteur[fam]=1
	for id in a[1:]:
		st = id.split("&")[0]
		if st == REF:
			parent[id] = fam

f.close()



for sub in  partitions[1:]:
	f=open(out_path + "families_" + sub ,"r")
	for l in f:
		a=l.strip("\n").split("\t")
		ID="-"
		thing=a[0]
		for id in a[1:]:
			st = id.split("&")[0]
			#print(st,REF)
			if st == REF:
				ID = id
		#print(ID)
		a.remove(thing)
		if ID in a:
			a.remove(ID)
		if ID in parent:
			#print("TAG1")
			fam = parent[ID]
			build[fam].extend(a)
			compteur[fam]+=1
	f.close()


h=open(out_path + "families_core.txt","w")
for fam in build:
	#print("BUILD: ",fam,compteur[fam] )
	if compteur[fam] == len(partitions):
		h.write(fam + "\t" + "\t".join(build[fam]) + "\n")

h.close()



sp="test"

seq={}
if 1==1:
	seq[sp]={}
	for fam in build:
		for id in build[fam]:
			seq[sp][id]=""


ext = "." +  REF.split(".")[-1]

if 1==1:
	tmp = os.listdir(path)
	for file in tmp:
		if file.endswith(ext) :
			if file in sub:
				toto={}
				f=open(path + file ,"r")
				for l in f:
					if l[0]==">":
						tag=0
						id = l.strip(">").strip("\n").split(" ")[0].split("\t")[0]
						if IDENTIFIANTS == "unique":
							id=id
						else:
							id = file + "&" + id
						if id in seq[sp]:
							tag=1
							toto[id]=[]
					elif tag==1:
						toto[id].append(l.strip("\n"))
				f.close()
				for id in toto:
					seq[sp][id] = "".join(toto[id]).upper()

nb = 0
for fam in build:
	if compteur[fam] == len(partitions):
		nb+=1
		#print(fam,compteur[fam],len(partitions))
		h=open(out_path + "core/" + fam + ext ,"w") 
		for id in build[fam][1:]:
			h.write(">" + id + "\n" + seq[sp][id] + "\n")
		h.close()


print("Fincal Core Genome= ",nb,"genes")
print("This core genome was computed into ",len(partitions)," subsets to deal with memry issues")








	





