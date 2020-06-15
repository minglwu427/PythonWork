# Python_Work_Script
Script for work

1. The NGS Analysis Data - This project is terminated.  The development of this code was kept internal.  
2. PeptideRegisCode - This code is used to turn a commericial dataframe of raw peptide we rodered and convert it into a csv file for our automation robot.  The code to do this function is in Automated_peptide_Regis.py.  In addition, Table_To_Plate.py are codes that convert data reading in table format and convert to plate format.
3. Sequence Separation - this code is used to split a string of petitde (Or a protein more precisely) into sequence of 9 and 10.  This is particularly useful for  immunology study.
4. Titration Line Graph.py was a script that takes raw data of cell base assay and return the result of graph via pyplot.  This script was retired after our company bought the license for prism.
5. Primer Designing tool - This project compare two genetic sequences and compare them with each other to define their mutation region.  Once the mutation region is found, it automatically generate primers by running a Tm calculator to get the best primer.   The inputs are the two csv files named Parent and Prodigy.  while the main codes (main.py), contains preprocessing (Join Tables/Cleaning etc), Primer selection (Compare the DNA sequence and generate primer), Tm_Calculator (Contains the Tm calculator base on previous publications).  The final output are in the output primer folder, the DNA sequence are compared with each other their primers information are shown in the data frame.    
