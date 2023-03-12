0.README

This is for the 340abs data set, H-chains and L-chains separate runs

PREPARE:

1. Create a working directory, i.e. "run-011623/"
2. cd "run-011623/"
3. create sub-directory "cntct/", "results/"
4. make a link to the PDB directory, "ln -s ../../pdb-340abs pdb", or copy it over
5. cp "340abs-pdbid.txt" here, and generate a "abs-namelist.txt" by "awk '{print $1}' 340abs-pdbid.txt"
6. cp "340abs-sequence-H-aligned-merge.txt" here
7. make sure the script folder is in the top of this directory, i.e. "../Scripts/"

RUN:

A. perl ../Scripts/evaluate-dis-v2.pl 340abs-pdbid.txt HorL >& run-abs-date.log &  (may take 45-60 minutes)

B. csh ../Scripts/run-post-process-dis.csh absH|absL 340absH|340absL
	(make sure change the name of "cntact" to "cntctH or cntctL")

C. csh ../Scripts/run-v7-examples.csh abs 340absH-contact-dis+fcs+cdr-fix-$date.txt 340abs-$date.
	(340absH or 340absL or 340absA)

BSA:

In a directory "bsa/", (a) cp over "pdbid.txt"; (b) link ../../pdb to here; (c) make a temp directory "tmp/":
perl ../../Scripts/compute-bsa.pl nbs|abs pdbid.txt H|L >& run-bsa-date.txt &


RESULTS:

The default directory is "results/" save to another name for different chain.
check the results in the directory, "results/"

CLUSTERS:

1. create subdirectory: "Clusters/", Copy over "340absH-0210-v7-6u-es.txt"
2. Run "./Scripts/es-clustering-v2.pl u 340absH-0210-v7-6u-es.txt 0.85 2 > es-clustering-85.txt"
3. (a) "ln -s ../pdb-340abs pdb"; (b) creatieve subdirectory "pdb1/"; (c) "cp ../340abs-pdbid.txt ." 
4. Run "./Scripts/pick-extract-nameedpdb-v2.pl pdbid.list es-clustering-99.txt cluster=#"

5. Selective clustering by specify classes or name or pdb, or ES set (4,5,6)
   "./Scripts/es-clustering-sel.pl 340abs-pdbid.list 340absH-0210-v7-6u-es.txt ES=class1" 
   # "selection" ES=name|pdbid|classid|4,5,6,..

For specific selections:

1. ./Scripts/es-clustering-sel.pl 340abs-pdbid.txt 340absH-0214-v7-6u-es.txt sim=0.90 sel=8,9,12,13,16,18,19 w=18 > es-cluster-90-S2.txt
2. ./Scripts/pick-extract-namedpdb-v2.pl 340abs-pdbid.txt es-cluster-90-S2.txt cluster=S2-90
3. awk '{print $2}' es-cluster-90-S2.txt > pickup-list.txt
4. ./Scripts/epitope-statis+fcs+cdr-v7d4.pl abs ../../340absH-contact-dis+fcs+cdr-fix-0207.txt 5a pickup-list.txt > foo.txt
5. egrep -v "RBD|CDR" foo.txt > es-cluster-90-S2-for-weblogo.txt
6. https://weblogo.berkeley.edu/logo.cgi
7. Loading up "es-cluster-90-S2-for-weblogo.txt",
8. Setup resolution = 392 bits; the plot range 332-529, Logo size: 36x18 cm
9. Individual ES range by the residue numbers.

