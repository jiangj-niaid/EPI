#!/bin/csh
# file: run-cns-dis.csh pdb name A B
# compute the distance between chain A and chain B within < 5.00 A
# select any two chains: A and B, distance cutoff C<5.0
#
# Jiansheng Jiang, 10-05-21
# Jiansheng Jiang, 03-30-22, v2
#
# Assume a "./pdb/" folder is in the current working directory stores all pdb files
#1. For each $pdb;
#2. cns < contact.inp > contact-$pdb.out
#3. mv contact.list contact-$pdb.list
#4. parse-contact.pl contact-$pdb.list > contact-dis-$pdb.txt


set p=$1
set n=$2
set pn=${p}"_"${n}
set A=$3
set B=$4
set C=5.0

set pdb=${p}".pdb"
set inp="contact-"${A}${B}".inp"
set out="contact-"${pn}".out"
set lis="contact-"${pn}"_"${A}${B}".list"
set txt="contact-dis-"${pn}"_"${A}${B}".txt"

# Assume a default CNS template "contact.inp" in "../Script/"
# Select-chain.pl modify the chain selection in "contact.inp" 
perl ../Scripts/select-chains.pl contact ${A} ${B} ${C}

if (-e "model.pdb") then
  rm -rf model.pdb
endif

# Assume a "./pdb/" folder in the current working directory
egrep -v "HETATM" ./pdb/${pdb} > model.pdb
# some pdb file missing "END" - fixing
#set ed=`tail -1 model.pdb`
#if (${ed} !~ "END") then
  echo "END" >> model.pdb
  echo "" >> model.pdb
#endif

# run CNS calculation of contact distance at the interaction surface
# CNS 1.3 should be in the path, i.e.
# /Users/jiangji/programs/osx/cns/cns_solve_1.3/mac-intel-darwin/bin/cns
# The distance can be modified in "contact.inp"
cns < ${inp} > ${out}

# "parse-contact.pl" parses the output list of distances
# results in a subfolder "./cntct/"
mv contact.list ${lis}
perl ../Scripts/parse-contact.pl ${lis} ${pn} ${A} ${B} > cntct/${txt}

sleep 2

rm -rf model.pdb
rm -rf ${inp} ${out} ${lis}

exit


