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

def median( echantillon) :
	echantillon.sort()
	size = len( echantillon )
	if len( echantillon ) % 2 == 0:
		M= float(echantillon[size / 2 - 1] + echantillon[size / 2]) / 2
	else:
		M= echantillon[size / 2]
	return M

def ninetyfive( echantillon) :
	echantillon.sort()
	size = len( echantillon )
	i95 = int( float(size) * 95/100 ) - 1
	return echantillon[i95]


########################################################################


def translate(seq):
	d={}
	d["TTT"],d["TTC"],d["TTA"],d["TTG"]="F","F","L","L"
	d["CTT"],d["CTC"],d["CTA"],d["CTG"]="L","L","L","L"
	d["ATT"],d["ATC"],d["ATA"],d["ATG"]="I","I","I","M"
	d["GTT"],d["GTC"],d["GTA"],d["GTG"]="V","V","V","V"
	d["TCT"],d["TCC"],d["TCA"],d["TCG"]="S","S","S","S"
	d["CCT"],d["CCC"],d["CCA"],d["CCG"]="P","P","P","P"
	d["ACT"],d["ACC"],d["ACA"],d["ACG"]="T","T","T","T"
	d["GCT"],d["GCC"],d["GCA"],d["GCG"]="A","A","A","A"
	d["TAT"],d["TAC"],d["TAA"],d["TAG"]="Y","Y","*","*"
	d["CAT"],d["CAC"],d["CAA"],d["CAG"]="H","H","Q","Q"
	d["AAT"],d["AAC"],d["AAA"],d["AAG"]="N","N","K","K"
	d["GAT"],d["GAC"],d["GAA"],d["GAG"]="D","D","E","E"
	d["TGT"],d["TGC"],d["TGA"],d["TGG"]="C","C","*","W"
	d["CGT"],d["CGC"],d["CGA"],d["CGG"]="R","R","R","R"
	d["AGT"],d["AGC"],d["AGA"],d["AGG"]="S","S","R","R"
	d["GGT"],d["GGC"],d["GGA"],d["GGG"]="G","G","G","G"
	i=0
	tmp=[]
	while i in range(len(seq)-3):
		codon= seq[i:i+3]
		if codon not in d.keys():
			tmp.append("X")
		else:
			tmp.append(d[codon])
			if d[codon]=="*":
				break
		i += 3
	prot="".join(tmp)
	return	prot



import os

import sys

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

score = float(score)/100

species=['Analysis']


print(species)


strains={}
for sp in species:
	strains[sp]=[]



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

liste={}
for sp in species:
	files = os.listdir(path)
	for truc in files:
		if truc.endswith(ext):
			if FILE != "none":
				if truc in sub:
					strains[sp].append(truc)
			else:
				strains[sp].append(truc)
	liste[sp] = strains[sp]



bbh={}
for sp in species:
	bbh[sp]={}
	BBH = os.listdir(out_path + 'BBH/')
	for file in BBH:
		#bbh[sp][file]='y'
		if 1==1:
			nb=0
			tag=0
			f=open(out_path + 'BBH/' + file , "r")
			for l in f:
				nb+=1
				a=l.strip("\n").split("\t")
				if len(a)<12:
					tag=1
			f.close()
			if nb >= 100 and tag != 1:
				bbh[sp][file]='y'
			else:
				print("Relaunch Usearch for ",file)

if TYPE=="DNA":
	option=" -strand plus "
elif TYPE =="protein":
	option=""

parent,dico={},{}
for sp in species:
	if 1==1:
		nb = len(liste[sp])
		if nb >= 1:
			print(sp,' ',nb,' strains')
			if 1==1:
				if 1==1:
					prot1 = REF
					for prot2 in liste[sp]:
						resu1= prot1 + '-' + prot2
						resu2= prot2 + '-' + prot1
						print(sp, ' ',prot1,' ',prot2)
						if 1==1:
							if resu1 in bbh[sp] and  resu2 in bbh[sp] :
								pass
							elif resu1 in bbh[sp]:
								print(sp,'/',nb)
								os.system('usearch61  -usearch_global  '   + path + prot2 + ' -db  ' + path  + prot1  + option +  ' -id ' + str(score) + ' -maxaccepts 3  -blast6out ' + out_path + 'BBH/' + prot2 + '-' + prot1 )
							elif resu2 in bbh[sp]:
								os.system('usearch61  -usearch_global  '  + path + prot1 + ' -db  ' + path  + prot2  + option + ' -id ' + str(score) + ' -maxaccepts 3  -blast6out ' + out_path + 'BBH/' + prot1 + '-' + prot2 )
							else:
								os.system('usearch61  -usearch_global  '  + path  + prot1 + ' -db  ' + path   + prot2  + option + ' -id ' + str(score) + ' -maxaccepts 3  -blast6out ' + out_path + 'BBH/' + prot1 + '-' + prot2 )
								os.system('usearch61  -usearch_global  '  + path + prot2 + ' -db  ' + path  + prot1 +  option + ' -id ' + str(score) + ' -maxaccepts 3  -blast6out ' + out_path + 'BBH/' + prot2 + '-' + prot1 )

		else:
			print(sp, ' empty')






















	
