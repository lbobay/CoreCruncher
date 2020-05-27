

print("USAGE: python master.py  -in input_folder   -out  output_folder   OPTIONAL: -freq  frequency_across_genomes  -ref ref_genome -prog usearch/blast -id unique/combined   -ext .fa/.fasta/.prot/.faa   -list path_to_genome_list   -score identify_score   -length sequence_length_conservation  -restart yes/no    -align muscle/mafft  -stringent yes/no") 
print("Use -h to see the different options\n")

import os
import sys

from datetime import datetime
startTime = datetime.now()


loc=""
for stuff in sys.argv:
	if "master.py" in stuff:
		loc = stuff.split("master.py")[0]
		print("Location=",loc)

path =  "/Users/lbobay/Desktop/dossiers/pan/test/database/"

out_path = "/Users/lbobay/Desktop/dossiers/pan/test/output/"

REF = "NA"

arguments =   sys.argv


if "-h" in arguments:
	print("\nUSAGE:")
	print("-in      input folder")
	print("-out     output folder")
	print("\nOPTIONS: ")
	print("-freq     	Minimum frequency of the gene across genomes to be considered core (default= 90%, an ortholog is considered a core gene even if it is missing in 10% of the set of genomes)")
	print("-score    	Identity score used by usearch or blast to define orthologs in % (Default= 70)")
	print("-length   	Minimum sequence length conservation used by to define orthologs (default= 80%)")
	print("-prog     	Program to use to compare sequences: usearch or blast (default= usearch)")
	print("-ref      	Reference genome (default: first genome in folder will be used as reference). If you want to specify the reference genome to use, speficy the name of the file in the folder (e.g. -ref genome1.prot)")
	print("-id       	Type of gene IDs in output files. Choose 'unique' if the same gene IDs are not found in different genomes or 'combined' to combine genome ID & gene ID (default= 'combined').  ")
	print("-ext	     	File extensions .fa/.fasta/.prot/.faa (default: will try to find it automatically)")
	print("-list     	Path to a file containing the list of genomes to analyze (default: none, all the genomes in the folder will be analyzed by default)")
	print("-restart  	Restart analysis from scratch: yes or no (default= no). If yes is chosen, the program will erase the usearch output files and relaunch usearch or blast")
	print("-align    	Align core gene sequences with specified program (muscle or mafft) and merge all the core genes into a single concatenate. Example=  -align musclev0.0.0 or -align mafft)")
	print("-stringent   Define a stringent core genome: yes or no (default= no). By default, core genes with paralogs will be conserved in the core genome and the paralogous sequences will be removed. If stringent is chosen, the core gene will be entirely removed from the core genome)")

	exit()

if "-in" in arguments:
	i = arguments.index("-in")
	path = arguments[i+1]
	if path[-1] == "/":
		pass
	else:
		path += "/"

if "-out" in arguments:
	i = arguments.index("-out")
	out_path = arguments[i+1]
	if out_path[-1] == "/":
		pass
	else:
		out_path += "/"

if "-freq" in arguments:
	i = arguments.index("-freq")
	freq = arguments[i+1]
else:
	freq="90"

if "-ref" in arguments:
	i = arguments.index("-ref")
	REF = arguments[i+1]

if "-ext" in arguments:
	i = arguments.index("-ext")
	ext = arguments[i+1]
else:
	ext = "NA"
	
if "-list" in arguments:
	i = arguments.index("-list")
	FILE = arguments[i+1]
else:
	FILE = "none"

PROG = "usearch"
if "-prog" in arguments:
	i = arguments.index("-prog")
	PROG = arguments[i+1].lower()

if "usearch" in PROG:
	PROG="usearch"
elif "blast" in PROG:
	PROG="blast"

if "-score" in arguments:
	i = arguments.index("-score")
	score= arguments[i+1]
else:
	score="70"

if "-length" in arguments:
	i = arguments.index("-length")
	length = arguments[i+1]
else:
	length = "80"

if "-id" in arguments:
	i = arguments.index("-id")
	IDENTIFIANTS = arguments[i+1]
else:
	IDENTIFIANTS= "combined"

if "-restart" in arguments:
	i = arguments.index("-restart")
	restart = arguments[i+1]
else:
	restart= "no"

ALIGN="muscle"
if "-align" in arguments:
	i = arguments.index("-align")
	ALIGN = arguments[i+1]


if "-stringent" in arguments:
	i = arguments.index("-stringent")
	stringent = arguments[i+1]
else:
	stringent= "no"

print("input = ", path)
print("output = ", path)


print("PARAMETERS: " + path + " " + out_path + " " + REF  + " " + ext + " " + FILE + " " + freq + " " + score + " " + length + " " + stringent)

try:
	os.mkdir(out_path)
except OSError:
	pass

try:
	os.mkdir(out_path + "BBH")
except OSError:
	pass

if restart == "yes":
	os.system("rm " + out_path + "BBH/*" )	


try:
	os.mkdir(out_path + "core")
except OSError:
	os.system("rm -r " + out_path + "core" )	
	os.mkdir(out_path + "core")

tmp = os.listdir(path)


extensions=[".fa",".fasta",".prot",".faa",".fna",".fnn",""]
if ext in extensions:
	extensions.remove(ext)


ext_found=[]
for stuff in tmp:
	if stuff[0] != ".":
		alternative = "." + stuff.split(".")[-1]
		if alternative not in ext_found:
			ext_found.append(alternative)

if ext == "NA":
	if len(ext_found)==1:
		ext = ext_found[0]
		print("Files with extension",ext, "have been found. These files will be used for analysis. You can specify a different extension with -ext")
	else:
		tmp=[]
		for stuff in ext_found:
			if stuff in extensions:
				if stuff not in tmp:
					tmp.append(stuff)
		if len(tmp)==1:
			ext = tmp[0]
			print("Files with extension",ext, "have been found. These files will be used for analysis. You can specify a different extension with -ext")
		else:
			print("Please specify the extension of the files that you want to analyze, the following have been found:",ext_found)

memo=[]
tag=0
files = []
for stuff in tmp:
	if stuff.endswith(ext):
		files.append(stuff)
		tag=1
	else:
		for point in extensions:
			if stuff.endswith(point):
				if point not in memo:
					memo.append(point)

if len(memo) > 0:
	print("\n######################################")
	print("WARNING: other files have been detected with the following extensions: "," ".join(memo))
	print("The analysis will be conducted on the ",ext, " files")
	print("If you want to analyze other files, specify the file extensions to use with the option -ext  (e.g.  -ext .fasta)")
	print("######################################\n")

if tag==0:
	print("\n######################################")
	print("No files have been found with extension " + ext)
	print("Please specify the file extension with option -ext (e.g. -ext .prot)")
	print("Exiting...")
	print("######################################\n")
	exit()

toto,seq,size={},{},{}
check=[]
nb=0
for stuff in files:
	if 1==1:
		nb+=1
		if nb <= 3 or nb >= len(files) - 1:
			f=open(path + stuff ,"r")
			for l in f:
				if l[0]==">":
					id = l.strip(">").strip("\n").split(" ")[0]
					id = stuff + "&" + id
					toto[id]=[]
				else:
					toto[id].append(l.strip("\n"))
			f.close()
			for id in toto:
				seq[id] = "".join(toto[id]).upper()
				size[id] = len(seq[id])
				break
			SEQ = seq[id]
			A = SEQ.count("A")
			C = SEQ.count("C")
			G = SEQ.count("G")
			T = SEQ.count("T")
			dash = SEQ.count("-")
			tot = A + C + G + T + dash
			ratio = float(tot) / size[id]
			check.append(ratio)

#print("Check= ",check

if min(check) < 0.5 and max(check) < 0.5 :
	print("TYPE= Protein sequences")
	TYPE="protein"
elif min(check) > 0.5 and max(check) > 0.5 :
	print("TYPE= DNA sequences")
	TYPE="DNA"
elif min(check) < 0.5 and max(check) > 0.5 :
	print("PROBLEM: BOTH DNA and protein sequences were given. Exiting...")
	exit()

print("Building core genome with ",len(files)," genomes")

if REF =="NA":
	REF= files[0]

print("Reference genome= ", REF)

print(datetime.now() - startTime)

print("TAG")
if PROG=="usearch":
	os.system("python " + loc + "usearch_core.py  " + path + " " + out_path + " " + REF + " " + TYPE + " " + ext + " " + FILE + " " + freq + " " + score + " " + length + " " + IDENTIFIANTS)
elif PROG=="blast":
	os.system("python " + loc + "blast_core.py  " + path + " " + out_path + " " + REF + " " + TYPE + " " + ext + " " + FILE + " " + freq + " " + score + " " + length + " " + IDENTIFIANTS)
else:
	print("UNKNOWN PROGRAM",PROG)

print("TAG2")
print("Sequence search complete")
print(datetime.now() - startTime)

os.system("python " + loc + "big_cruncher.py  " + stringent + " "  + path + " " + out_path + " " + REF + " " + TYPE + " " + ext + " " + FILE + " " + freq + " " + score + " " + length + " " + IDENTIFIANTS)


if "-align" in arguments:
	print("Launch alignments with ",ALIGN)
	os.system("python " + loc + "align.py " + ALIGN + " " + out_path + " " + ext)
	os.system("python " + loc + "concat.py " + out_path + " " + ext )




print(datetime.now() - startTime)














