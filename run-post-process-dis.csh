#!/bin/csh
# file: run-post-process-dis.csh $option data-name
# processes after "evaluate-dis.pl"
#
# $option: "nbs", "absH" or "absL". "abs" means "absH"
# $da: data name i.e. "217abs" or "42nbs"
#
# Jiansheng Jiang, 06-11-22
# Jiansheng Jiang, 02-07-23, v2
#
#

# option
set opi=$1
# data name, i.e. "217abs" or "42nbs"
set da=$2
set daf=${da}"-contact-dis"
# data directory name
set cnt="cntct"
# date stamp
set pid=`date "+%m%d"`


# CDR tags file
set cdr=${da}"-define-cdr-tagged-v3.txt"
# FCS definition file
set fcs="../Scripts/Epitope-ES-CL-define5.txt"

# output data files
set dfile1=${daf}"-"${pid}".txt"
set dfile2=${daf}"+cdr-"${pid}".txt"
set dfile3=${daf}"+fcs+cdr-"${pid}".txt"
set dfile4=${daf}"+fcs+cdr-fix-"${pid}".txt"



# after "run-evaluate-dis.csh"
# The folder "cntct/" containts all individual contact distances for each antibodies


# 1. gather together for post process
cat ${cnt}/*contact-dis-*.txt > ${dfile1}
# This data file format is
# PDB_NAME	CHAINS	RBD	Rid	AB	Rid	DIST	(7 data terms)


# 2. add, read & write the CDR tags
#
# this needs a file of CDRloop defined tags that was generated with
# an sequence alignment file by "Cluster Omega" (merged)
# i.e. "217abs-sequence-aligned-merge.txt"
#set seqf=$2
#set seqf="217abs-sequence-H-aligned-merge.txt"
perl ../Scripts/define-align-cdr-v3.pl ${opi} ${da}-sequence-aligned-merge.txt > ${cdr}
# "define-align-cdr-v2.pl" is generating the CDR tags to antibody's sequence
#
# Now add this cdr-tags to the contact-dis data
perl ../Scripts/readnwrite-cdr-tags-v2.pl ${opi} ${cdr} ${dfile1} > ${dfile2}
#
# This data file format is
#PDB_NAME	CHAINS	RBD	Rid	AB	Rid	DIST	CDR	(8 data terms)


# 3. add, read & write the ES numbers
#
perl ../Scripts/readnwrite-es-v2.pl ${fcs} ${dfile2} > ${dfile3}
# This data file format is
#PDB_NAME	CHAINS	RBD	Rid	FCS	AB	Rid	DIST	CDR 	(9 data iterms)


# 4. fix some problems in RBD residue numbering - only for abs
# i.e. Beta-26,Beta-32,Beta-49,Beta-50
perl ../Scripts/fix-rbd-numbering.pl ${dfile3} > ${dfile4}

#


exit


