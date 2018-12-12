# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 16:17:59 2018

@author: Ming Lun Wu
"""
from Bio import SeqIO # module needed for sequence input
import re
import pandas as pd
import matplotlib.pyplot as plt

def main():
    fastq = open("Test.fastq", "rU")
    Bcode_list = [
            "AACACTTGGGCA",
            "AACAACCGACAC",
            "AAAATGTCTGCC",
            "AACGGATAACTG",
            "AACCCTAACCTA",
            'AAGGCTCTTACA',
            'AATAGGTGCAAC',
            'AAATACACGCCA',
            'AATTGAGAGCGA',
            'AAGAAACCTCCG']
    Pools = [[] for i in range(len(Bcode_list))]
    
    Sep_Seq_Data = SortAndReturn(Bcode_list,fastq,Pools)
    #print(len(Data.loc["AACACTTGGGCA"][100].seq))
    SeqLen = len(Sep_Seq_Data.loc['AACACTTGGGCA'][0].seq) # this is a bad way to do these code 
    SeqAndCA = SeqBpAna(Sep_Seq_Data, SeqLen)
    
    BpBarGraphs = []
    for eachseq in SeqAndCA:
        Graph = SequenceGraph(eachseq,250)
        BpBarGraphs.append(Graph)
    print(BpBarGraphs)
    
def SequenceGraph(BpCountData,NumOfMaxBP):
    """
    This function turns a dataframe of ATCG in row and BP location in column and turns it
    into a stacked barchart to provide a visualization of the sequence data.
    """
    transData = BpCountData.transpose()
    transData.plot(kind= 'bar', stacked = True)
    plt.xticks([])
    return plt.show()
    #plt.savefig('First sequence', dpi = 1000)

    """
    ResultToCSV(SeqAndCA)
    """
    

    
def SortAndReturn(Bcode_list,fastq,Pools):
    """
    this set of code gather a fastq data and a list of barcode and generate a panda 
    dataframe with index of barcode and each sequence read
    """
    for record in SeqIO.parse(fastq, "fastq"):
        for Bcode in Bcode_list:
            Match = re.search(r'%s' %Bcode, str(record.seq)) 
            if Match:
                Pools[Bcode_list.index(Bcode)].append(record)

    dataframe = pd.DataFrame(Pools)
    dataframe['Barcode'] = Bcode_list
    dataframe = dataframe.set_index('Barcode')
    return dataframe

def SeqBpAna(SeqData,SeqLen):
    """
    This function turns analysis your separated sequence data and see at each position the occurance of each
    pair

    SeqData = an array of data that separated the each sequence according to its Barcode
    SeqLen = How many bp are in your sequence
    """
    bpCountA = []
    for bcode,sequences in SeqData.iterrows(): # Turn 
        bpCount = pd.DataFrame(0,index= {"A","C","T","G"}, columns = range(0,SeqLen),dtype = int)
        for sequ in sequences:
            if (sequ is not None):
                for pos in range(len(sequ.seq)):
                    bp = sequ.seq[pos]
                    bpCount.loc[bp][pos] =    bpCount.loc[bp][pos] + 1
        bpCountA.append(bpCount)
    return bpCountA

def ResultToCSV(SeqAndCA):
    """
    Saved Data into CSV
    """
    i = 0
    for dataframe in SeqAndCA:
        i = i + 1
        #plt.savefig("Well" + fset['Well Position'][0]+'.png')
        dataframe.to_csv('C:/Users/Ming Lun Wu/Desktop/Python Stuff/NGS Analysis Data/Seq%s.csv' %i)

if __name__ == "__main__":
    main()

