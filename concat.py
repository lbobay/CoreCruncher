
import os
import sys
sp ="SP"
species=[sp]



folder = sys.argv[-2]
if folder[-1]=="/":
	pass
else:
	folder=folder + "/"
	
ext =  sys.argv[-1]

sp = "SP"
species=[sp]
strains={}
strains[sp]=[]


bad={}
bad[sp]={}



parent={}
orthologues={}
for sp in species:
	orthologues[sp]={}
	f=open(folder + 'families_core.txt',"r")
	for l in f:
		a=l.strip("\n").split("\t")
		if a[0] in bad[sp]:
			pass
		else:
			orthologues[sp][a[0]]=a[1:]
			for id in a[1:]:
				linker = ext + "&"
				st = id.split("&")[0]
				parent[id] = st
				if st not in strains[sp]:
					strains[sp].append(st)
				
	f.close()

print(strains[sp])

genes={}
for sp in species:
	genes[sp]=[]
	files = os.listdir(folder + "core")
	for stuff in files:
		if stuff.endswith(".align"):
			genes[sp].append(stuff)

#print(genes[sp])

numbers={}
for sp in species:
	numbers[sp]=[]
	for stuff in genes[sp]:
		nb = int(stuff.lstrip("fam").split(ext)[0])
		numbers[sp].append(nb)
	numbers[sp].sort()

#print(numbers[sp])

exclusion=[]


tirets='------------------------------------------------------------------'
tmp={}
for sp in species:
	print(sp)
	tmp[sp]={}
	for nb in numbers[sp]:
		ortho="fam" + str(nb) + ext + ".align"
		#print(ortho)
		f=open(folder + 'core/'  + ortho ,"r")
		memo=[]
		for l in f:
			if l[0] == ">":
				st = l.strip("\n").strip(">").rstrip('\r').split("&")[0]
				memo.append(st)
				if st in tmp[sp]:
					pass
				else:
					tmp[sp][st]=''
			else:
				tmp[sp][st]+=l.strip("\n").upper()
				longueur = len(tmp[sp][st])
		f.close()
		for st in strains[sp]:
			if st in tmp[sp]:
				pass
			else:
				tmp[sp][st]=''
			if len(tmp[sp][st]) < longueur:
				while len(tmp[sp][st]) < longueur:
					#print st, ' ',len(tmp[sp][st]),' ',longueur
					tmp[sp][st]+='-'



concat={}
for sp in species:
	concat[sp]={}
	for id in tmp[sp]:
		concat[sp][id] = tmp[sp][id]


tmp={}






for sp in species:
	h=open(folder + 'concat' + ext,"w")
	for st in strains[sp]:
		h.write(">" + st + "\n")
		i=0
		while i < len(concat[sp][st]):		# MODIF 
			h.write(concat[sp][st][i:i+60] + "\n")
			i+=60
	h.close()


