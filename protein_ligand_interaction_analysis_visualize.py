"""
This script read count csv
generate visualization script of chimerax
"""

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

def generate_rgb_values(num_bins:int)->list:
    """
    input number of bins output a list of hex values of color
    """
    # Create a colormap that transitions from blue to white to red
    cmap = plt.get_cmap('bwr')

    # Generate an array of equally spaced values between 0 and 1
    values = np.linspace(0, 1, num_bins)

    # Get the RGB values for each value in the array
    rgb_values = [cmap(value)[:3] for value in values]

    # Convert RGB values to hexadecimal format
    hex_values = ['#%02x%02x%02x' % (int(r*255), int(g*255), int(b*255)) for r, g, b in rgb_values]

    return hex_values


def generate_count_bin(csv_input:str, csv_output:str, bin_num:int,group_by:list)->pd.DataFrame:
    
    df = pd.read_csv(csv_input)
    # generate count
    new_df = df.groupby(group_by).size().reset_index(name='count')
    # classifer bins' id by count
    new_df['count_bins'] = pd.cut(new_df['count'], bins=bin_num, labels=False)
    new_df.to_csv(csv_output,index=False)
    return new_df


def generate_chimerax_script(df:pd.DataFrame,cxc_output:str, hex_values:list):
    
    new_df = df
    with open(cxc_output, "w") as f:
        f.write("color #1 blue tar c \n")
        for index, row in new_df.iterrows():
            print(f"Row {index+1}: count_bins = {row['count_bins']}")
            row_bin = int(row['count_bins'])
            #cmd =f"color #1:{row['global_resindex']} {hex_values[row_bin]} tar c \n"
            cmd =f"color #1:{row['resnum']} {hex_values[row_bin]} tar c \n"
            print(cmd)
            f.write(cmd)
        print(f"write to {f.name}")


def main(csv_input:str, 
         csv_output:str, 
         cxc_output:str, 
         bin_num:int=10, 
         group_by_list:list=['global_resindex', 'resname', 'resnum']):
    
    """
    input:
    csv_input: str, path to input csv
    csv_output: str, path to output csv
    cxc_output: str, path to output cxc
    bin_num: int, number of bins
    group_by_list: list, list of columns to group by
    
    output:
    None"""
    hex_values = generate_rgb_values(bin_num)
    new_df = generate_count_bin(csv_input, csv_output, bin_num, group_by_list)
    generate_chimerax_script(new_df,cxc_output, hex_values)

    print(f"input csv file \n> {csv_input}")
    print(f"output csv file \n> {csv_output}")
    print(f"output cxc file \n> {cxc_output}")
    print(f"number of bins \n> {bin_num}")
    print(f"group by list \n> {group_by_list}")
    print("Done")


if __name__ == "__main__":
#    # input
#    bin_num = 10
#    #csv_data_dir = "data/ligand_contact/distance_4_5_eachligand/"
#    #case_name= "contacteachlig_GJZX_ACLY_0001_100w_GX1_1100"
#    #group_by_list = ['global_resindex', 'resname', 'resnum']
#    
#    csv_data_dir = 'data/ligand_contact/distance_4_5_allligand/'
#    case_name = 'contacteachlig_GJZX_ACLY_0001_100w_GX1'
#    
#    group_by_list = ['resname', 'resnum']
#
#    csv_input = os.path.join(csv_data_dir, f"{case_name}.csv")
#    csv_output = os.path.join(csv_data_dir,f'{case_name}_count_bin{bin_num}.csv')
#    cxc_output = os.path.join(csv_data_dir,f'{case_name}_count_bin{bin_num}.cxc')
#
#
#    # processing
#    main(csv_input, csv_output, cxc_output, bin_num, group_by_list)


    import glob
    from logging_util import logging
    from functools import partial

    logging = partial(logging, logfile=__file__+'.log')


    bin_num = 10
    group_by_list = ['resname', 'resnum']
    
#    # round1
#    for i in glob.glob('data/ligand_contact/distance_4_5_eachligand/*.csv'):
#        csv_input = i
#        csv_output = i.replace('.csv', f'_count_bin{bin_num}.csv').replace('data/ligand_contact/distance_4_5_eachligand','data/ligand_contact/distance_4_5_eachligand_analysis')
#        cxc_output = i.replace('.csv', f'_count_bin{bin_num}.cxc').replace('data/ligand_contact/distance_4_5_eachligand','data/ligand_contact/distance_4_5_eachligand_analysis')
#        try:
#            main(csv_input, csv_output, cxc_output, bin_num, group_by_list)
#            print('done')
#            logging(f"done {i}")
#        except:
#            print(f"error in {i}")
#            logging(f"error in {i}")
#    
#    for i in glob.glob('data/ligand_contact/distance_4_5_allligand/*.csv'):
#        csv_input = i
#        csv_output = i.replace('.csv', f'_count_bin{bin_num}.csv').replace('data/ligand_contact/distance_4_5_allligand','data/ligand_contact/distance_4_5_allligand_analysis')
#        cxc_output = i.replace('.csv', f'_count_bin{bin_num}.cxc').replace('data/ligand_contact/distance_4_5_allligand','data/ligand_contact/distance_4_5_allligand_analysis')
#
#        try:
#            main(csv_input, csv_output, cxc_output, bin_num, group_by_list)
#            print('done')
#            logging(f"done {i}")
#        except:
#            print(f"error in {i}")
#            logging(f"error in {i}")
#
    # round2
    for i in glob.glob('data/ligand_contact/distance_4_5_allligand_combine_groupby_highIC50/*.csv'):
        csv_input = i
        csv_output = i.replace('.csv', f'_count_bin{bin_num}.csv').replace('data/ligand_contact/distance_4_5_allligand_combine_groupby_highIC50','data/ligand_contact/distance_4_5_allligand_combine_groupby_highIC50_analysis')
        cxc_output = i.replace('.csv', f'_count_bin{bin_num}.cxc').replace('data/ligand_contact/distance_4_5_allligand_combine_groupby_highIC50','data/ligand_contact/distance_4_5_allligand_combine_groupby_highIC50_analysis')
        try:
            main(csv_input, csv_output, cxc_output, bin_num, group_by_list)
            print('done')
            logging(f"done {i}")
        except:
            print(f"error in {i}")
            logging(f"error in {i}")

    for i in glob.glob('data/ligand_contact/distance_4_5_allligand_combine_groupby_lowIC50/*.csv'):
        csv_input = i
        csv_output = i.replace('.csv', f'_count_bin{bin_num}.csv').replace('data/ligand_contact/distance_4_5_allligand_combine_groupby_lowIC50','data/ligand_contact/distance_4_5_allligand_combine_groupby_lowIC50_analysis')
        cxc_output = i.replace('.csv', f'_count_bin{bin_num}.cxc').replace('data/ligand_contact/distance_4_5_allligand_combine_groupby_lowIC50','data/ligand_contact/distance_4_5_allligand_combine_groupby_lowIC50_analysis')
        try:
            main(csv_input, csv_output, cxc_output, bin_num, group_by_list)
            print('done')
            logging(f"done {i}")
        except:
            print(f"error in {i}")
            logging(f"error in {i}")

    data_in_dir = 'data/ligand_contact/distance_4_5_allligand_combine'
    for i in glob.glob(f'{data_in_dir}/*.csv'):
        csv_input = i
        csv_output = i.replace('.csv', f'_count_bin{bin_num}.csv').replace(data_in_dir,f'{data_in_dir}_analysis')
        cxc_output = i.replace('.csv', f'_count_bin{bin_num}.cxc').replace(data_in_dir,f'{data_in_dir}_analysis')
        try:
            main(csv_input, csv_output, cxc_output, bin_num, group_by_list)
            print('done')
            logging(f"done {i}")
        except:
            print(f"error in {i}")
            logging(f"error in {i}")