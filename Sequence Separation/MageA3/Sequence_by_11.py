# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 13:57:24 2018

@author: Ming Lun Wu
"""

import re
import pandas as pd

def fragmentation(seq_file,frag_size):
    seq_file = open(seq_file)
    sequence = seq_file.readline()
    sequence = re.sub('[*\n-]+',"",sequence)
    fragments = []
    for i in range(len(sequence)):
        fragments.append(sequence[i:i+frag_size])
    fragments = [ frag for frag in fragments if len(frag) == frag_size]
    return fragments    

def main():
    fragments = fragmentation('C:/Users/Ming Lun Wu/Desktop/sequence.txt', 11)
    fragments = pd.DataFrame(fragments, columns = ['sequence'])

    fragments.to_csv('C:/Users/Ming Lun Wu/Desktop/sequence_separated.csv')
    
if __name__ == '__main__':
    main()