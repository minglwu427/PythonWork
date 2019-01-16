# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 11:54:37 2019

@author: Ming Lun Wu
"""
import pandas as pd
import re

def Load_Clean_Table_Form(FileName):
    Table = pd.read_excel(FileName) 
    return Table

def match_plate_with_order(Order, plate_start_num, vertical):
    """
    Adopted code from Automate Pep Regis Order is Basically Table here
    """
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
        if vertical == True:
            temp_row = temp_row + 1 
            if temp_row == 9:
                temp_row = 1
                temp_column = temp_column +1
                if temp_column ==  13:
                    temp_column = 1
                    temp_plate = temp_plate + 1
        else:
            temp_column = temp_column + 1 
            if temp_column ==  13:
                temp_column = 1
                temp_row = temp_row +1
                if temp_row == 9:
                    temp_row = 1
                    temp_plate = temp_plate + 1
    Order['Plate'] = plate
    Order['Column'] = column
    Order['Row'] = row
    column_dict = {1:'A',2:'B',3:'C',4:'D',5:'E',6:'F',7:'G',8:'H'}
    Order['Row'] = Order['Row'].map(column_dict)
    return Order

def clean_slice_dataframe(Table):
    plates = []
    TempTable = Table[['Order ID', 'Column','Row', 'Volume (uL)','Plate']]
    Number_of_Plate = Table['Plate'].max()
    for i in range(1,Number_of_Plate+1):
        plates.append(TempTable[TempTable['Plate'] == i])
    return plates

def pivot_dataframe_to_map(plates_dataframe):
    volume_maps = []
    orderID_maps = []
    for plates in plates_dataframe:
        temp_volume = plates.pivot(index = 'Row', columns ='Column', values = 'Volume (uL)')
        temp_orderID = plates.pivot(index = 'Row', columns ='Column', values = 'Order ID')
        volume_maps.append(temp_volume)
        orderID_maps.append(temp_orderID)
    return volume_maps, orderID_maps
    
def dfs_tabs(df_list, sheet_list, file_name ):
    """
    This is a script copy from stackoverflow with minor modification
    to fix it into our purposes.  It help you write in multiple dataframe 
    across separate tabs/sheets
    
    shout out of TomDobbs for posting this answer online
    https://stackoverflow.com/questions/32957441/putting-many-python-pandas-dataframes-to-one-excel-worksheet
    """

    writer = pd.ExcelWriter(file_name,engine='xlsxwriter')   
    for dataframe, sheet in zip(df_list, sheet_list):
        dataframe.to_excel(writer, sheet_name=sheet, startrow=0 , startcol=0)   
    writer.save()
    
def multiple_dfs(maps, sheet_list, file_name, spaces):
    writer = pd.ExcelWriter(file_name,engine='xlsxwriter')   
    
    for df_list,sheet in zip(maps,sheet_list):
        row = 0
        for dataframe in df_list:
            dataframe.to_excel(writer,sheet_name=sheet,startrow=row , startcol=0)   
            row = row + len(dataframe.index) + spaces + 1
    writer.save()

def main():
    Table = Load_Clean_Table_Form('U9724DL040.xlsx')
    Table = match_plate_with_order(Table, 1, vertical=False)
    plates_dataframe= clean_slice_dataframe(Table)
    volume_maps, orderID_maps = pivot_dataframe_to_map(plates_dataframe)
    maps = [volume_maps, orderID_maps]
    sheets_name = ['volume_maps', 'orderID_maps']
    multiple_dfs(maps, sheets_name, 'multi-test.xlsx',spaces=2)
    '''
    Table = Table[Table['Plate']==1]
    TempTable = Table[['Order ID', 'Column','Row']]
    New_Table = TempTable.pivot(index = 'Row', columns ='Column', values = 'Order ID')
    '''
    
    #pd.DataFrame.pivot_table(TempTable,values = 'Item number',index = 'Row', columns = 'Column')
if __name__ == "__main__":
    main()
