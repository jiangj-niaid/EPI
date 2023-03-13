# EPI ![image](https://user-images.githubusercontent.com/75588673/224588219-5080cb8f-83c0-4d87-86ac-88bada48967d.png)

Description

To analysis the binding interface of antibody/nanobody in complex with Receptor Binding Domain (RBD) of Spike (SARS-CoV-2), we developed a software package of EPI (Epitope-Paratope Interaction).  The package is a mixed scripts of C-shell, perl and python language and utilize CNS 1.3 (http://cns-online.org/v1.3/) (Brunger, Adams et al. 1998) and PISA (Proteins, Interfaces, Structures and Assemblies (Krissinel and Henrick 2007)).  The package is able to analyze large number of structures of complexes of antibody and antigen from PDB (Protein Data Bank).
![image](https://user-images.githubusercontent.com/75588673/224588248-6a721908-05f3-4666-94a2-7e06f913df78.png)

Features

1.	For a giving list of antibody and spike structures (pdb ids) that are downloaded from PDB,  we compute the contacted interface (the distance is cut-off 5.0 angstroms), build up a sub-database as the bases for further analysis
2.	Identify the epitope sites (common or frequently contact sites on RBD surface)
3.	Analysis the contacted amino acids of antibody in binding to Epitope Sites (ES) on the RBD
4.	Statistical analysis the CDR loops contributions
5.	Clustering analysis and identify the binding motifs of antibody
6.	Write out Pymol scripts for graphical display
7.	The user can make query the epitope-paratope interaction data for a specific antibody/pdb-id/class/ES-sites, and ask for the similar antibody those bind to the same ES sites for comparison

![image](https://user-images.githubusercontent.com/75588673/224588369-bad6d523-0867-471a-9967-b2e40d6c88d1.png)

Documentation

Further Reading

Acknowledgements

![image](https://user-images.githubusercontent.com/75588673/224588575-0feb7e83-38e0-4d45-ac55-ccfcfec8612d.png)
