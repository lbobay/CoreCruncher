

import os

########################################################################
import math

 
def mean( echantillon ) :
    size = len( echantillon )
    moyenne = float(sum( echantillon )) / float(size)
    return moyenne


def stat_variance( echantillon ) :
    n = float(len( echantillon )) # taille
    mq = mean( echantillon )**2
    s = sum( [ x**2 for x in echantillon ] )
    variance = s / n - mq
    return variance


def stat_ecart_type( echantillon ) :
    variance = stat_variance( echantillon )
    ecart_type = math.sqrt( variance )
    return ecart_type


def median_indice(sorted_list):
	sorted_list.sort()
	indices = []
	list_size = len(sorted_list)
	median = 0
	if list_size % 2 == 0:
		indices.append(int(list_size / 2) - 1)  # -1 because index starts from 0
		indices.append(int(list_size / 2))
		median = (sorted_list[indices[0]] + sorted_list[indices[1]]) / 2
	else:
		indices.append(int(list_size / 2))
		median = sorted_list[indices[0]]
	return indices


def quartile1( samples ):
	ind = median_indice(samples)
	Q1 = numpy.median(samples[:ind[0]])
	return Q1

def quartile3( samples ):
	ind = median_indice(samples)
	Q3 = numpy.median(samples[ind[-1] + 1:])
	return Q3

########################################################################


SP="test"

import sys
import numpy

path = sys.argv[-10]
out_path = sys.argv[-9]
REF  = sys.argv[-8]
TYPE = sys.argv[-7]
ext = sys.argv[-6]
FILE = sys.argv[-5]
freq = int(sys.argv[-4])
score = int(sys.argv[-3])
length = int(sys.argv[-2])
IDENTIFIANTS =  sys.argv[-1]


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


print("load sequences")
check=[]
genomes={}
size={}
parent={}
seq={}
for sp in species:
	nb=0
	seq[sp]={}
	genomes[sp]=[]
	size[sp]={}
	tmp = os.listdir(path)
	for file in tmp:
		if file.endswith(ext) :
			if file in sub:
				nb+=1
				genomes[sp].append(file)
				toto={}
				f=open(path + file ,"r")
				for l in f:
					if l[0]==">":
						id = l.strip(">").strip("\n").split(" ")[0]
						if IDENTIFIANTS == "unique":
							id=id
						else:
							id = file + "&" + id
						toto[id]=[]
						parent[id]=file
					else:
						toto[id].append(l.strip("\n"))
				f.close()
				for id in toto:
					seq[sp][id] = "".join(toto[id]).upper()
					size[sp][id] = len(seq[sp][id])
					del seq[sp][id]


seq={}
test={}
toto={}
parent={}


print("IDENTIFIANTS= ",IDENTIFIANTS)
if IDENTIFIANTS == 	"Redundant":
	print("WARNING: Identical gene IDs have been found in different genomes. The gene IDs used in the output files will be a combination of genome IDs & gene IDs")

south = float(length)/100
north = 2 - south
print("length overlap= ",south," ",north)
print("build dico")
dico={}
for sp in species:
	ref = REF
	#genomes[sp].append(ref)
	dico[sp]={}
	for st1 in genomes[sp]:
		dico[sp][st1]={}
		for st2 in genomes[sp]:
			dico[sp][st1][st2]={}
			if st1 == ref or st2==ref:
				if st1!=st2:
					nb=0
					f=open(out_path +  "BBH/" + st1 + "-" + st2 ,"r")
					for l in f:
						a=l.strip("\n").split("\t")
						if len(a) > 3:
							nb+=1
							if IDENTIFIANTS == "unique":
								id1 =  a[0].split(" ")[0]
								id2 =  a[1].split(" ")[0]				
							else:
								id1 = st1 + "&" + a[0].split(" ")[0]
								id2 = st2 + "&" + a[1].split(" ")[0]
							#print(id1," ",id2," ",a[2:]
							if id1 in dico[sp][st1][st2]:
								pass
							else:
								dico[sp][st1][st2][id1] = {}
							s1,s2 = float(size[sp][id1]),float(size[sp][id2])
							if s1/s2 >= south and s1/s2 <=north:
								dico[sp][st1][st2][id1][id2] = float(a[2])
						else:
							"Usearch issue: File BBH/" + st1 + "-" + st2," did not finish running"
					f.close()
					if nb == 0:
						"Usearch issue: File BBH/" + st1 + "-" + st2," is empty. You should remove empty genome or relaunch Usearch"

size={}

print("Build best")
best={}
for sp in species:
	best[sp]={}
	for st1 in genomes[sp]:
		best[sp][st1]={}
		for st2 in genomes[sp]:
			#print(st1," ",st2
			best[sp][st1][st2]={}
			if st1 != st2:
				if st1 == ref or st2==ref:
					for id1 in dico[sp][st1][st2]:
						#print("dico ",id1
						memo =""
						top=0
						for id2 in dico[sp][st1][st2][id1]:
							if dico[sp][st1][st2][id1][id2]  > top:
								top =  dico[sp][st1][st2][id1][id2]
								memo = id2
						if memo != "":
							best[sp][st1][st2][id1] = [memo,top]

# Check best
#for sp in best:
#	for st1 in best[sp]:
#		for st2 in best[sp][st1]:
#			#print(st1," ",st2
#			if st1 != st2:
#				for id1 in best[sp][st1][st2]:
#					#print("best ", id1
#					pass

print("Build bbh")
bbh={}
for sp in species:
	bbh[sp]={}
	for st1 in genomes[sp]:
		bbh[sp][st1]={}
		for st2 in genomes[sp]:
			bbh[sp][st1][st2]={}
			nb=0
			if st1 != st2:
				if st1 == ref or st2==ref:				
					for id1 in best[sp][st1][st2]:
						#print(id1
						id2 = best[sp][st1][st2][id1][0]
						#print(id1 ," ",best[sp][st1][st2][id1]
						try:
							rev = best[sp][st2][st1][id2][0]
							if rev == id1:
								bbh[sp][st1][st2][id1] = best[sp][st1][st2][id1]
								nb+=1
						except KeyError:
							pass
					#print(st1," ",st2," ",nb
	
best={}

# Check bbh
#for sp in species:
#	for st1 in genomes[sp]:
#		#print("Tag1 ",st1)
#		for st2 in genomes[sp]:
#			#print(st1," ",st2)
#			if st1 != st2:
#				if st1 == ref or st2==ref:
#					nb=0
#					tmp=[]
#					for id1 in bbh[sp][st1][st2]:
#						nb+=1
#						tmp.append(bbh[sp][st1][st2][id1][1])
#					#print("Check bbh: ",st1," ",st2," ",nb," ",sum(tmp) / float(len(tmp)))

# initiate
link={}
families={}
for sp in species:
	link[sp],families[sp]={},{}
	f=open(path + REF,"r")
	nb = 0
	for l in f:
		if l[0] == ">":
			nb+=1
			if IDENTIFIANTS == "unique":
				resu = l.strip("\n").strip(">").split(" ")[0]
			else:
				resu =   REF + "&" + l.strip("\n").strip(">").split(" ")[0]
			fam = "fam" + str(nb)
			link[sp][resu] = fam
			families[sp][fam] = [resu]
			#print(fam," ",resu
	f.close()


problematic=[]
for sp in species:
	for st1 in genomes[sp]:
		for st2 in genomes[sp]:
			tag=0
			#print(st1," ",st2
			if st1 != st2 and REF in [st1,st2] :
				nb=0
				#print("Tag2 ",st1," ",st2	
				for id1 in bbh[sp][st1][st2]:
					id2 =bbh[sp][st1][st2][id1][0]
					resu1,resu2=id1,id2
					#resu1 = st1 + "&" + id1
					#resu2 = st2 + "&" + id2
					score = bbh[sp][st1][st2][id1][1]
					if resu1 in link[sp]:
						fam1 = link[sp][resu1]
						if resu2 in link[sp]:
							fam2 = link[sp][resu2]
							if fam1==fam2:
								pass
							else:
								print("PROBLEM ",fam1," ",fam2)
								if fam1 not in problematic:
									problematic.append(fam1)
								if fam2 not in problematic:
									problematic.append(fam2)
						else:
							link[sp][resu2] = fam1
							if 1==1:
								nb+=1
								#print("second round ",id2," ",nb
							if resu2 not in families[sp][fam1]:
								families[sp][fam1].append(resu2)
								
				if 1==1:
					#print(nb," pairs"
					pass


link={}

vertical={}
for sp in families:
	vertical[sp]={}
	for fam in families[sp]:
		vertical[sp][fam]=[]
		id_ref = families[sp][fam][0]
		for resu in families[sp][fam]:
			if IDENTIFIANTS == "unique":
				st = parent[resu]
			else:
				st = resu.split("&")[0]
			if st != REF:
				try:
					score = bbh[sp][REF][st][id_ref][1]
					vertical[sp][fam].append(score)
				except KeyError:
					pass
		vertical[sp][fam].sort()
		
		#print("Family= ",fam," ",len(families[sp][fam])


selected={}
total=0
core={}
for sp in species:
	core[sp]=[]
	cutoff = len(genomes[sp]) * float(freq) / 100
	if cutoff - int(cutoff)> 0:
		cutoff = int(cutoff) + 1.0
	print("cutoff of " ,freq , "% > core gene if present in at least ",int(cutoff)," genomes over ",len(genomes[sp]))
	for fam in families[sp]:
		total+=1
		if fam not in problematic:
			#print(fam," ",families[sp][fam]
			if len(families[sp][fam]) >= cutoff:
				core[sp].append(fam)
				for resu in families[sp][fam]:
					selected[resu]= "y"



							
print(len(core[sp])," core genes / ",total," ")

# Now check for paralogs based on horizontal and vertical distributions
# Detect double outliers
MEDIAN,SD,QUARTILES={},{},{}
horizontal={}
for sp in species:
	horizontal[sp]={}
	MEDIAN[sp],SD[sp],QUARTILES[sp]={},{},{}
	for st1 in genomes[sp]:
		horizontal[sp][st1]={}
		MEDIAN[sp][st1],SD[sp][st1],QUARTILES[sp][st1]={},{},{}
		for st2 in genomes[sp]:
			if st1 != st2 and REF in [st1,st2] :
				horizontal[sp][st1][st2]=[]
				for id1 in bbh[sp][st1][st2]:
					id2 = bbh[sp][st1][st2][id1][0]
					resu1,resu2=id1,id2
					#resu1 = st1 + "&" + id1
					#resu2 = st2 + "&" + id2
					if resu1 in selected and resu2 in selected:
						score = bbh[sp][st1][st2][id1][1]
						horizontal[sp][st1][st2].append(score)
				horizontal[sp][st1][st2].sort()
				if len(horizontal[sp][st1][st2]) >= 4:
					MEDIAN[sp][st1][st2],SD[sp][st1][st2]=numpy.median(horizontal[sp][st1][st2]),stat_ecart_type(horizontal[sp][st1][st2])
					QUARTILES[sp][st1][st2] = [quartile1(horizontal[sp][st1][st2]),quartile3(horizontal[sp][st1][st2])]
					#print("horizontal: ",st1," ",st2," ",sum(horizontal[sp][st1][st2])/len(horizontal[sp][st1][st2])," ",horizontal[sp][st1][st2][:10]
				else:
					print("BBH file ",st1," ",st2," seems empty")
					MEDIAN[sp][st1][st2],SD[sp][st1][st2]=0,0
					QUARTILES[sp][st1][st2]=[0,0]


horizontal={}

# Two steps:
# Step 1 identify families with potential paralogs
# Step 2 check if score is in the range of orthologs (in case of divergent genome)
filtered_core={}
check=[]
for sp in species:
	compteur=0
	filtered_core[sp]=[]
	#h=open("../test/reference.fa","w")
	#g=open("../test/core_ref.fa","w")
	total=0
	nb=0
	ref = REF
	tag=0
	for fam in core[sp]:
		total+=1
		for resu in families[sp][fam]:
			if IDENTIFIANTS == "unique":
				id = resu
				st=parent[id]
			else:
				a= resu.split("&")
				st,id= a[0],resu
			if st == ref:
				ref_id = id
		# Step 1
		paralogs,orthologs=[],[]
		for resu in families[sp][fam]:
			if IDENTIFIANTS == "unique":
				st,id= parent[resu],resu
			else:
				a= resu.split("&")
				st,id= a[0],resu
			if st != ref:
				#print(ref," ",st," ",ref_id," ",dico[sp][ref][st][ref_id]
				if ref_id in dico[sp][ref][st]:
					for para in dico[sp][ref][st][ref_id]:
						if para != id:
							#if fam == "fam1582":
								#print("paralog: ",fam," ",para," ",dico[sp][ref][st][ref_id][para]
							paralogs.append(dico[sp][ref][st][ref_id][para])
						else:
							#if fam == "fam1582":
								#print("ortholog: ",fam," ",para," ",dico[sp][ref][st][ref_id][para]," ",mean(horizontal[sp][ref][st])
							orthologs.append(dico[sp][ref][st][ref_id][para])
			#print(ref_id," ",tmp
		if len(paralogs) == 0:
			print(fam," no paralogs")
			tag=1
			#h.write(">" + fam + "&" + ref_id + "\n" + seq[sp][ref_id] + "\n")
			if ref_id in selected:
				pass
				#g.write(">" + fam + "&" + ref_id + "\n" + seq[sp][ref_id] + "\n")
				if fam not in check:
					check.append(fam)
				else:
					print(fam," already there")
		elif len(paralogs) > 0:
			#print(fam," ",max(paralogs)," ",min(orthologs)
			if max(paralogs)<=min(orthologs):
				print(fam," max(paralogs)<=min(orthologs)")
				tag=1
				#h.write(">" + fam + "&" + ref_id + "\n" + seq[sp][ref_id] + "\n")
				#g.write(">" + fam + "&" + ref_id + "\n" + seq[sp][ref_id] + "\n")
			else:
				#print("\n Case: ", fam," ",paralogs," ",orthologs)
				nb+=1
		if tag==1:
			#print("Check ",fam," ",len(families[sp][fam])," ",families[sp][fam]
			outlier="no"
			memo=[]
			vert = numpy.median(vertical[sp][fam])
			try:
				vert_sd = stat_ecart_type(vertical[sp][fam])
			except ValueError:
				#print(vertical[sp][fam]
				vert_sd = 3.0
			try:
				Q1_vert,Q3_vert= quartile1(vertical[sp][fam]),quartile3(vertical[sp][fam])
				down_vert = Q1_vert - 1.5 * (Q3_vert - Q1_vert)
				up_vert = Q3_vert + 1.5 * (Q3_vert - Q1_vert)
			except ValueError:
				down_vert = vert - 3*vert_sd
				up_vert = vert + 3*vert_sd
			id_ref = families[sp][fam][0]
			for st in genomes[sp]:
				if st != REF:
					try:
						m = MEDIAN[sp][REF][st]
						sd = SD[sp][REF][st]
						#down,up = m - 3*sd , m + 3*sd
						Q1,Q3 = QUARTILES[sp][REF][st][0],QUARTILES[sp][REF][st][1]
						down,up = Q1 - 1.5 * (Q3 - Q1) ,  Q3 + 1.5 * (Q3 - Q1)
						score = bbh[sp][REF][st][id_ref][1]
						if score >= down and score <= up:
							pass
						else:
							#print(score," ",down,"-",up
							if score >= down_vert and score <= up_vert:
								#print("Horizontal outlier"
								pass
							else:
								memo.append([score,down,up])
								outlier ="y"
								print(fam," double outlier ",score," ",down,"-",up,"    and   ",down_vert," ",up_vert," SD= ",vert_sd)
					except KeyError:
						pass
			if outlier=="y":
				compteur+=1
			else:
				filtered_core[sp].append(fam)
			
						
	#h.close()
	#g.close()

vertical={}
dico={}
core={}
bbh={}

print("Excluded: ", nb," genes")
			
print("There are ",compteur," families with double outliers\n")

seq={}
for sp in species:
	seq[sp]={}
	for fam in filtered_core[sp]:
		for id in families[sp][fam]:
			seq[sp][id]=""


for sp in species:
	tmp = os.listdir(path)
	for file in tmp:
		if file.endswith(ext) :
			if file in sub:
				toto={}
				f=open(path + file ,"r")
				for l in f:
					if l[0]==">":
						tag=0
						id = l.strip(">").strip("\n").split(" ")[0]
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


print("Writing")
for sp in species:
	g=open(out_path + "families_core.txt","w")
	for fam in filtered_core[sp]:
		resu = fam
		h=open(out_path + "core/" + fam + ext ,"w") 
		for id in families[sp][fam]:
			#if IDENTIFIANTS == "unique":
			#	h.write(">" + parent[id] + "&" + id + "\n" + seq[sp][id]  + "\n")
			#else:
			h.write(">" + id + "\n" + seq[sp][id]  + "\n")
			#if IDENTIFIANTS == "unique":
			#	resu += "\t" + parent[id] + "&" + id
			#else:
			resu += "\t" + id
		h.close()
		g.write(resu + "\n")

print("\n######################################")
print("Final core genome= ",len(filtered_core[sp])," genes")
print("ouput written in ",out_path)
print("The file 'families_core.txt' contains the list of orthologous genes")
print("The directory 'core' contains the sequences each orthologous gene")
print("Reference genome used for the analysis= ",REF)
print("######################################")










