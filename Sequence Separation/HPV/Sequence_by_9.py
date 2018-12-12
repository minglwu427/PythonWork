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
    filesN = ['HPV16E6','HPV16E7','HPV18E6','HPV18E7']
    fragmentsfiles = pd.DataFrame()
    for files in filesN:
        fragments = fragmentation('%s.txt' %files, 9)
        temp_fragments = pd.DataFrame(fragments,columns = [files+' by 9'])
        fragmentsfiles = pd.concat([fragmentsfiles,temp_fragments], ignore_index = False, axis = 1)
        
        fragments = fragmentation('%s.txt' %files, 10)
        temp_fragments = pd.DataFrame(fragments,columns = [files+' by 10'])
        fragmentsfiles = pd.concat([fragmentsfiles,temp_fragments], ignore_index = False, axis = 1)
        
    fragmentsfiles.to_csv('C:/Users/Ming Lun Wu/Desktop/HPV_Seq_Separated.csv')

if __name__ == '__main__':
    main()