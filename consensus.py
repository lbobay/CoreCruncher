

import sys

args = sys.argv


if "-folders" in args:
	i = args.index("-folders")
else:
	print("\nIncorrect options. Exiting...")
	print("USAGE:     python consensus.py   -folders  FOLDER1  FOLDER2  -in genome_folder  -out  output_path\n" )
	print("-folders   The two output folders of CoreCruncher built on different pivot genomes to build the consensus")
	print("-in        Folder with the original gene or protein sequences of the genomes")
	print("-out       PATH to write the output files")
	exit()


folders=[]
j=i+1
while j < len(args):
	fofo = args[j]
	if fofo=="-out" or fofo == "-in":
		break
	if fofo.endswith("/"):
		pass
	else:
		fofo = fofo + "/"
	folders.append(fofo)
	j+=1

if "-out" in args:
	i=args.index("-out")
	outpath= args[i+1]
	if outpath[-1] == "/":
		pass
	else:
		outpath=outpath + "/"
else:
	outpath = ""



if "-in" in args:
	i=args.index("-in")
	path= args[i+1]
	if path[-1] == "/":
		pass
	else:
		path=outpath + "/"
else:
	print("Pleasse privide the input folder where the genome sequences are located with -in ")
	print("Exiting...")
	exit()


print(len(folders)," folders will be analyzed:",folders)

reference={}
for fofo in folders:
	f=open(fofo + "summary.txt","r")
	for l in f:
		a=l.strip("\n").split(" ")
		if a[0] == "Pivot":
			reference[fofo]=a[-1]
	f.close()



parent={}
dico={}
for fofo in folders:
	parent[fofo]={}
	dico[fofo]={}
	f=open(fofo + "families_core.txt","r")
	for l in f:
		a=l.strip("\n").split("\t")
		dico[fofo][a[1]] = a[1:]
		for id in a[1:]:
			parent[fofo][id] = a[1]
	f.close()

for fofo in dico:
	print("Core genome of ",fofo,len(dico[fofo].keys()))

final={}
maybe={}
new={}

nomatch=0
match=0
exact=0
inexact=0
link={}
done=[]
MEMO=[]

for stuff1 in folders:
	for stuff2 in folders:
		if stuff1 != stuff2:
			tmp=[stuff1,stuff2]
			tmp.sort()
			if tmp not in done:
				done.append(tmp)
				fofo1 = tmp[0]
				fofo2 = tmp[1]
				MEMO.append(fofo1)
				ref1 = reference[fofo1]
				ref2 = reference[fofo2]
				for fam1 in dico[fofo1]:
					for id1 in dico[fofo1][fam1]:
						st1 = id1.split("&")[0]
						if st1 == ref1:
							pass
						if st1== ref2:
							link[dico[fofo1][fam1][0]] = id1
							link[id1] = dico[fofo1][fam1][0]
					if dico[fofo1][fam1][0] in link:
						fam2 = link[dico[fofo1][fam1][0]]
						if fam2 in dico[fofo2]:
							tmp=[]
							for id2 in dico[fofo2][fam2]:
								if id2 in dico[fofo1][fam1]:
									tmp.append(id2)
							#print(len(tmp),len(dico[fofo1][fam1]),len( dico[fofo2][fam2]) )
							if len(tmp) == min( len(dico[fofo1][fam1]),len( dico[fofo2][fam2]) ) :
								exact+=1
								final[fam1] = dico[fofo1][fam1]
							else:
								inexact +=1
								tmp1,tmp2=[],[]
								for id in dico[fofo1][fam1]:
									resu = parent[fofo1][id]
									if resu not in tmp1:
										tmp1.append(resu)
								for id in dico[fofo2][fam2]:
									resu = parent[fofo2][id]
									if resu not in tmp2:
										tmp2.append(resu)			
								#print(len(tmp1),len(tmp2))	
								if len(tmp1) == 1 and len(tmp2) == 1:
									if len(dico[fofo1][fam1]) <= len(dico[fofo2][fam2])	:
										final[fam1] = dico[fofo1][fam1]
									else:
										final[fam1] = dico[fofo2][fam2]	
								#maybe[fam1] = [dico[fofo1][fam1],dico[fofo2][fam2]]
							match+=1
						else:
							#print("No match:",dico[fofo1][fam1][0])
							nomatch+=1
							tmp1,tmp2=[],[]
							for id in dico[fofo1][fam1]:
								resu1 = parent[fofo1][id]
								if resu1 not in tmp1:
									tmp1.append(resu1)
								if id in parent[fofo2]:
									resu2 = parent[fofo2][id]
									if resu2 not in tmp2:
										tmp2.append(resu2)
							#print("TAG2",len(tmp1),len(tmp2))
							if len(tmp2)==0:
								final[fam1] = dico[fofo1][fam1]
							elif tmp1[0] == tmp2[0]:
								final[fam1] = dico[fofo1][fam1]
								#print(tmp1,tmp2)
							else:
								print("TAG3")
										
					else:
						#print("No match:",dico[fofo1][fam1][0])
						tmp1,tmp2=[],[]
						for id in dico[fofo1][fam1]:
							resu1 = parent[fofo1][id]
							if resu1 not in tmp1:
								tmp1.append(resu1)
							if id in parent[fofo2]:
								resu2 = parent[fofo2][id]
								if resu2 not in tmp2:
									tmp2.append(resu2)
						if len(tmp2)==0:
							final[fam1] = dico[fofo1][fam1]
						else:
							if tmp1[0] == tmp2[0]:
								print("TAG4")
								final[fam1] = dico[fofo1][fam1]
						#print("TAG1",len(tmp1),len(tmp2))
						nomatch+=1

print(nomatch,match)
print(exact,inexact)



memory={}
for fam in final:
	for id in final[fam]:
		if id not in memory:
			memory[id]="y"
		else:
			print("BIG PROBLEM ",id)



print("partial=",len(final.keys()))

new_gene=0
new={}
done=[]
print(MEMO)
for stuff1 in folders:
	for stuff2 in folders:
		if stuff1 != stuff2:
			if stuff1 in MEMO:
				fofo1 = stuff2
			else:
				fofo1 = stuff1
			if fofo1 not in done:
				print(fofo1)
				done.append(fofo1)
				ref1 = reference[fofo1]
				ref2 = reference[fofo2]
				for fam1 in dico[fofo1]:
					for id1 in dico[fofo1][fam1]:
						st1 = id1.split("&")[0]
						if st1 == ref1:
							pass
						if st1== ref2:
							link[dico[fofo1][fam1][0]] = id1
							link[id1] = dico[fofo1][fam1][0]
				for fam1 in dico[fofo1]:
					if fam1 not in memory:
						tag=0
						for fam in final:
							if tag ==0:
								for id in final[fam]:
									if id == fam1:
										tag=1
										break
						if tag==0:
							new[fam1] = list(dico[fofo1][fam1])
							#print("NEW GENE")
							new_gene+=1

print("new genes",new_gene)




for fam in new:
	tag=0
	for id in new[fam]:
		if id in memory:
			tag=1
	if tag==1:
		print(fam," already exists")
	else:
		final[fam] = list(new[fam])
		for id in new[fam]:
			memory[id]="y"




sp="test"
seq={}
if 1==1:
	nb=0
	seq[sp]={}
	tmp = os.listdir(path)
	for file in tmp:
		if 1==1: #if file.endswith(ext) :
			if 1==1:
				nb+=1
				genomes[sp].append(file)
				toto={}
				f=open(path + file ,"r")
				for l in f:
					if l[0]==">":
						id = l.strip(">").strip("\n").split(" ")[0].split("\t")[0]
						id = file + "&" + id
						toto[id]=[]
					else:
						toto[id].append(l.strip("\n"))
				f.close()
				for id in toto:
					if id in memory:
						seq[sp][id] = "".join(toto[id]).upper()


tmp = list(seq[sp].keys())
id = tmp[0]

SEQ = seq[id]
A = SEQ.count("A")
C = SEQ.count("C")
G = SEQ.count("G")
T = SEQ.count("T")

tot = A+C+G+T

if float(tot) / len(SEQ) > 0.9:
	ext = ".fa"
else:
	ext = ".prot"


print("CONSENSUS CORE GENOME:",len(final.keys()))

print("Writing output file consensus_core.txt in ",outpath)
h=open(outpath + "consensus_core.txt","w")

try:
	os.mkdir(oupath + "core")
except OSError:
	pass

nb=0
for fam in final:
	nb+=1
	name = "fam" + str(nb)
	h.write(name + "\t" + "\t".join(final[fam]) + "\n")
	g=open(oupath + "core/" + fam + ext)
	for id in final[fam]:
		g.write(">" + id + "\n" + seq[sp][id] + "\n")
	g.close()

h.close()






