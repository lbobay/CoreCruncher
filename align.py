
import os
import sys
import subprocess

ALIGN = sys.argv[-3]
out_path = sys.argv[-2]
ext = sys.argv[-1]

FILES = os.listdir(out_path + "core")




stuff = ALIGN.lower()
if "muscle" in stuff:
	version = subprocess.check_output("muscle -version", shell=True)
	version = version.lower()
	sub = int(version.split(".")[0][-1])
	print("Calling muscle with '",ALIGN,"' with default parameters")
	for seq in FILES:
		if seq.endswith(ext):
			if sub <= 4:
				os.system(ALIGN + " -in " + out_path + "core/" + seq + "  -out " + out_path + "core/" + seq + ".align")
			else:
				os.system(ALIGN + " -align " + out_path + "core/" + seq + "  -output " + out_path + "core/" + seq + ".align")
elif "mafft" in stuff:
	print("Calling mafft with '",ALIGN,"' with default parameters")
	for seq in FILES:
		if seq.endswith(ext):
			os.system(ALIGN  + "  --quiet  " + out_path + "core/" + seq + "  > " + out_path + "core/" + seq + ".align")
else:
	print("Alignment program '",ALIGN,"' is unknown. Exiting...")
	exit()


