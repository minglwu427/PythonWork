# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from Bio import AlignIO
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
from Bio import Align
from Bio.pairwise2 import format_alignment
from collections import defaultdict
import pandas as pd
from Bio.Alphabet import generic_dna
import numpy as np
from Tm_Calculation import *

"""
Open a fastq file and covert it Translation and turn it into output
"""

def Extract(lst): 
    return [item[0] for item in lst] 

table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', } 


inv_map = defaultdict(list)
{inv_map[v].append(k) for k, v in table.items()}


mother = list(SeqIO.parse("Deletion Example\Deletion Example DNA Sequences.fasta", "fasta"))
mother = mother[0]

#oh = open("out.txt", "w")

def translation_to_AA(fasta_file):
    translate = []
    OG_SEQ = []
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences = record.translate()
            sequences.id = record.id
            sequences.name = record.name
            translate.append(sequences)
            OG_SEQ.append(str(record.seq))
    
    
    SeqIO.write(translate, "AA_translate.fasta", "fasta")
    return translate,OG_SEQ


translate, OG_SEQ = translation_to_AA("Delete_Parent.fasta")

def align_AA(translate):
    OG_AA = []
    PP_AA = []
    mother = translate[0]
    for AA_seq in translate:
        pair_result = pairwise2.align.globalxx(str(mother.seq), str(AA_seq.seq))
        #alignments.append([pair_result[0][0], pair_result[0][1]])
        OG_AA.append(pair_result[0][0])
        PP_AA.append(pair_result[0][1])
        
    return OG_AA, PP_AA


OG_AA, PP_AA = align_AA(translate)


final_dataframe = pd.DataFrame()
final_dataframe['OG SEQ'] = OG_SEQ
final_dataframe['OG AA'] = OG_AA
final_dataframe['PP AA'] = PP_AA

def define_mutation(OG_Seq, PP_Seq):
    if OG_Seq.count('-')  == 0 and PP_Seq.count('-')  != 0 :
        return 'Deletion'
    elif OG_Seq.count('-') !=0 and PP_Seq.count('-')==0:
        return 'Insertion'
    elif OG_Seq.count('-') !=0 and PP_Seq.count('-') !=0 and OG_Seq.count('-')<PP_Seq.count('-'):
        return "Deletion and Substitution"
    elif OG_Seq.count('-') !=0 and PP_Seq.count('-') !=0 and OG_Seq.count('-')>PP_Seq.count('-'):
        return 'Insertion and Substitution'
    elif OG_Seq.count('-') !=0 and PP_Seq.count('-') !=0 and OG_Seq.count('-')==PP_Seq.count('-'):
        return  "Substitution"
    elif OG_Seq.count('-')==0 and PP_Seq.count('-')==0 :
        return "No Mutations"
   
def assign_mutation(final_dataframe,OG_SEQ,OG_AA,PP_AA):
    definition=[]
    for i in range(len(OG_SEQ)):
        definition.append(define_mutation(OG_AA[i], PP_AA[i]))
    final_dataframe['Mutation Type'] = definition
    return final_dataframe

final_dataframe = assign_mutation(final_dataframe, OG_SEQ, OG_AA,PP_AA)
mother = str(mother.seq)



def reverse_translation(final_dataframe, mother_sequence):
    OG_translated_DNA = []
    PP_translated_DNA= []
    for i in range(len(final_dataframe)): 
        index = list(final_dataframe.iloc[i])
        OG_DNA_RT = ""
        PP_DNA_RT = ""
        OG_nuc = mother_sequence
        PP_nuc = index[0]
        OG_AA_aligned = index[1]
        PP_AA_aligned = index[2]
        Mutation = index[3]
        if Mutation == 'Deletion':
            for j in range(len(OG_AA_aligned)):
                if OG_nuc[3*j: 3*j+3 ] == PP_nuc[3*j: 3*j+3 ]: # handle no mistake
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3 ]
                    PP_DNA_RT = PP_DNA_RT + OG_nuc[3*j: 3*j+3 ]
                elif (OG_nuc[3*j: 3*j+3 ] != PP_nuc[3*j: 3*j+3 ]) and OG_AA_aligned[j] == PP_AA_aligned : #handle silent Mutation
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3]
                    PP_DNA_RT = PP_DNA_RT + PP_nuc[3*j: 3*j+3]
                elif OG_nuc[3*j: 3*j+3 ] != PP_nuc[3*j: 3*j+3 ] and PP_AA_aligned[j] == '-' :
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3 ]
                    PP_DNA_RT = PP_DNA_RT + "---"
                    PP_nuc = PP_nuc[:j] + '---' + PP_nuc[j:]
            OG_translated_DNA.append(OG_DNA_RT)
            PP_translated_DNA.append(PP_DNA_RT)
            
        if Mutation == 'Insertion':
            for j in range(len(PP_AA_aligned)):
                if OG_nuc[3*j: 3*j+3 ] == PP_nuc[3*j: 3*j+3 ]:
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3 ]
                    PP_DNA_RT = PP_DNA_RT + OG_nuc[3*j: 3*j+3 ]
                elif (OG_nuc[3*j: 3*j+3 ] != PP_nuc[3*j: 3*j+3 ]) and OG_AA_aligned[j] == PP_AA_aligned[j] : 
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3]
                    PP_DNA_RT = PP_DNA_RT + PP_nuc[3*j: 3*j+3]
                elif OG_nuc[3*j: 3*j+3 ] != PP_nuc[3*j: 3*j+3 ] and OG_AA_aligned[j] == '-' : #handle insertion
                    OG_DNA_RT = OG_DNA_RT + "---"#OG_nuc[3*j: 3*j+3 ]
                    PP_DNA_RT = PP_DNA_RT + PP_nuc[3*j: 3*j+3 ]
                    OG_nuc = OG_nuc[:j] + '---' + OG_nuc[j:]
            OG_translated_DNA.append(OG_DNA_RT)
            PP_translated_DNA.append(PP_DNA_RT)
        if Mutation == 'Substitution':
            OG_AA_aligned = OG_AA_aligned.replace('-','')
            PP_AA_aligned = PP_AA_aligned.replace('-','')
            for j in range(len(PP_AA_aligned)):
                if OG_nuc[3*j: 3*j+3 ] == PP_nuc[3*j: 3*j+3 ]:
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3 ]
                    PP_DNA_RT = PP_DNA_RT + OG_nuc[3*j: 3*j+3 ]
                elif (OG_nuc[3*j: 3*j+3 ] != PP_nuc[3*j: 3*j+3 ]) and OG_AA_aligned[j] == PP_AA_aligned[j] : 
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3]
                    PP_DNA_RT = PP_DNA_RT + PP_nuc[3*j: 3*j+3]
                else: #as long as they are reverse translating its own line everything should be fine
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3]
                    PP_DNA_RT = PP_DNA_RT + PP_nuc[3*j: 3*j+3]
            OG_translated_DNA.append(OG_DNA_RT)
            PP_translated_DNA.append(PP_DNA_RT)
        if Mutation == 'Insertion and Substitution':
            OG_AA_aligned.count('-')
            OG_AA_aligned = OG_AA_aligned.replace('-','',PP_AA_aligned.count('-'))
            PP_AA_aligned = PP_AA_aligned.replace('-','',PP_AA_aligned.count('-'))
            for j in range(len(PP_AA_aligned)):
                if OG_nuc[3*j: 3*j+3 ] == PP_nuc[3*j: 3*j+3 ]:
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3 ]
                    PP_DNA_RT = PP_DNA_RT + OG_nuc[3*j: 3*j+3 ]
                elif (OG_nuc[3*j: 3*j+3 ] != PP_nuc[3*j: 3*j+3 ]) and OG_AA_aligned[j] == PP_AA_aligned[j] : 
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3]
                    PP_DNA_RT = PP_DNA_RT + PP_nuc[3*j: 3*j+3]
                elif OG_nuc[3*j: 3*j+3 ] != PP_nuc[3*j: 3*j+3 ] and OG_AA_aligned[j] == '-' : #handle insertion
                    OG_DNA_RT = OG_DNA_RT + "---"#OG_nuc[3*j: 3*j+3 ]
                    PP_DNA_RT = PP_DNA_RT + PP_nuc[3*j: 3*j+3 ]
                    OG_nuc = OG_nuc[:j] + '---' + OG_nuc[j:]
                else: #as long as they are reverse translating its own line everything should be fine
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3]
                    PP_DNA_RT = PP_DNA_RT + PP_nuc[3*j: 3*j+3]
            OG_translated_DNA.append(OG_DNA_RT)
            PP_translated_DNA.append(PP_DNA_RT)
        if Mutation == 'No Mutations':
            OG_DNA_RT = OG_nuc
            PP_DNA_RT = OG_nuc
            OG_translated_DNA.append(OG_DNA_RT)
            PP_translated_DNA.append(PP_DNA_RT)
    return OG_translated_DNA, PP_translated_DNA 

OG_translated_DNA, PP_translated_DNA = reverse_translation(final_dataframe, mother)

final_dataframe['OG_translated_DNA'] = OG_translated_DNA
final_dataframe['PP_translated_DNA'] = PP_translated_DNA

del(translate, OG_translated_DNA, PP_translated_DNA) 

def generate_raw_primers(final_dataframe):
    forward_primer_list = []
    backward_primer_list = []
    mutation_region_list = []
    for i in range(len(final_dataframe)): 
        index = list(final_dataframe.iloc[i])
        first_error_location = 0
        find_first_error = False
        find_last_error = False
        last_error_location = 0
        OG_DNA_RT = index[4]
        PP_DNA_RT = index[5]
        for j in range(len(OG_DNA_RT)):
            if (OG_DNA_RT[j] != PP_DNA_RT[j]) and find_first_error == False: 
                find_first_error = True
                first_error_location = j
            elif  (OG_DNA_RT[j:j+20] == PP_DNA_RT[j:j+20]) and find_first_error == True and find_last_error == False:
                find_last_error = True
                last_error_location = j
            else:
                continue
        if find_first_error == True and find_last_error == True :
            MM_region = PP_DNA_RT[first_error_location:last_error_location].replace('-',"")
            backward_primer = PP_DNA_RT[first_error_location-30: first_error_location]
            backward_primer = str(Seq(backward_primer).reverse_complement()) #reverse compliment
            forward_primer  = PP_DNA_RT[last_error_location: last_error_location+30]
            mutation_region_list.append(MM_region)
            forward_primer_list.append(forward_primer)
            backward_primer_list.append(backward_primer)
        else:
            mutation_region_list.append('No Error')
            forward_primer_list.append('No Error')
            backward_primer_list.append('No Error')
    return forward_primer_list,  backward_primer_list, mutation_region_list
    

forward_primer_list,backward_primer_list, mutation_region_list=generate_raw_primers(final_dataframe)
final_dataframe['forward_primer'] = forward_primer_list
final_dataframe['backward_primer'] = backward_primer_list
final_dataframe['mutation_region'] = mutation_region_list

del(forward_primer_list,backward_primer_list, mutation_region_list) 




"""
    OG_codon = OG_SEQ[0][3*i: 3*i+3 ]
    if OG_codon in inv_map[OG_AA[1][i]]:
        print('True')
    else:
        print('False')


primer1 = 'GTGGGCCGGTTCACCATCAGCCGGGACAAC'
(primer1.count('C') + primer1.count('G'))/len(primer1)
C = primer1.count('C')
A = primer1.count('A')
G = primer1.count('G')
T = primer1.count('T')
Temperature = 4*(G+C) + 2*(A+T)
Temperature = Temperature+ np.log(500E-9)
#Temperature = 64.9 + (4*(G+C-16.4)/(A+T+G+C))
"""

#final_dataframe.to_excel('Datafame.xlsx')