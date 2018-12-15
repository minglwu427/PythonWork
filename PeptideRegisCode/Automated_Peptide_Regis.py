# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 14:56:24 2018

@author: Ming Lun Wu
"""

import pandas as pd
import re

def Load_Clean_Order_Form(FileName = 'Order.xlsx'):
    Order = pd.read_excel(FileName) 
    Order['Sequence'] = Order['Sequence'].str.replace("\t",'')
    return Order


def check_duplicate_sequence(Order):
    Seq = Order['Sequence'].unique()
    if len(Seq) == len(Order):
        print('There is no repeating')
    else:
        print('There are repeating sequence')
    
def match_plate_with_order(Order, plate_start_num):
    plate=[]
    row =[]
    column = []
    temp_row = 1
    temp_column = 1
    temp_plate = plate_start_num
    index = len(Order)
    index_counter = 0
    while index_counter < index:
        index_counter = index_counter + 1 
        row.append(temp_row)
        column.append(temp_column)
        plate.append(temp_plate)
        temp_row = temp_row + 1 
        if temp_row == 13:
            temp_row = 1
            temp_column = temp_column +1
            if temp_column ==  9:
                temp_column = 1
                temp_plate = temp_plate + 1
    Order['Plate'] = plate
    Order['Column'] = column
    Order['Row'] = row
    column_dict = {1:'A',2:'B',3:'C',4:'D',5:'E',6:'F',7:'G',8:'H'}
    Order['Column'] = Order['Column'].map(column_dict)
    return Order


def check_repeating(master_list, order):
    master_Seq = master_list['Sequence'].unique()
    order['Repeating'] = order['Sequence'].isin(master_Seq)
    return order

def assign_a2pep_num(Order, sample_start_num, ali_start_num,master_list): 
    ali_counting_num = ali_start_num
    sample_counting_num = sample_start_num
    for i in Order.index:
        if Order.loc[i,'Repeating']  == False :
            Order.at[i, 'A2 Pep Number'] = 'PEP' + str(sample_counting_num).zfill(8) + '-0'
            Order.at[i, 'lot Number'] = 0
            Order.at[i, 'Pep_Regis_Num'] = 'PEP' + str(sample_counting_num).zfill(8)
            sample_counting_num+= 1
            """
            Still need to add Unique Plate Number, 
            Positioning depend on how the Freezer is configurated
            """
        else:
            temp = master_list[master_list['Sequence'] == Order.loc[i,'Sequence']]
            temp = temp.sort_values(by=['lot Number'])
            Order.at[i, 'lot Number'] = temp['lot Number'].iloc[-1] + 1
            Order.at[i, 'Pep_Regis_Num'] =  temp['Pep_Regis_Num'].iloc[-1]
            Order.at[i, 'A2 Pep Number'] = Order.at[i, 'Pep_Regis_Num'] +'-'+  str(Order.at[i, 'lot Number'])
            # find the repeating number and add lot plus one
        Order.loc[i,'Pep_Ali_ID'] = 'A' + str(ali_counting_num).zfill(7) 
        Order.loc[i,'Vender Name'] = ''
        Order.loc[i,'Vender ID'] = ''
        Order.loc[i,'Source Protein'] = ''
        Order.loc[i,'UniProt ID'] = ''

        ali_counting_num += 1
    return Order

def main():
    master = pd.read_excel('TestMaster.xlsx') 
    Order = Load_Clean_Order_Form()
    check_duplicate_sequence(Order)
    Order = match_plate_with_order(Order, plate_start_num = 1 )
    sample_start_num = 0 
    aliqout_start_num = 0 
    Order = check_repeating(master, Order)
    Order = assign_a2pep_num(Order, sample_start_num, aliqout_start_num, master)
    
    '''
    Note this code can't currently handle if you were to order two peptide with
    identical sequence at once

    '''
if __name__ == "__main__":
    main()
