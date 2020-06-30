# CoreCruncher
Requires Usearch (32 bits) or BLAST (tested successfully with BLAST 2.7.1+)
Please download USEARCH: https://www.drive5.com/usearch/ or BLAST

OPTIONAL: can align core genomes and provide sequence concatenate if alignment program is available
Alignment program muscle or mafft must be executable and in /usr/local/bin

################################################################
IMPORTANT: Rename the usearch executable 'usearch61'
copy or move usearch61 into /usr/local/bin/
Or change the lines 176 to 181 in the script usearch_core.py
################################################################

Requires the python library numpy

Simple example:
python corecruncher_master.py -in example -out out_folder
