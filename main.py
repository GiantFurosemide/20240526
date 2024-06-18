import json
import glob
import os
import pandas as pd
import MDAnalysis as mda
from csv_util import *
from json_util import *
from protein_ligand_interaction_analysis import contact_count,contact_count_for_pool
from multiprocessing import Pool

raw_data_list = read_json('data_all_ligand_name.json')


case_id2ligand_name = {}
for rec in raw_data_list:
    case_id2ligand_name[rec['case_id']] = rec['ligand_name']


input_list = [(rec['case_id'],[rec['gro']]+rec['traj']) for rec in raw_data_list]
input_dict = {}
for data in input_list:
    input_dict[data[0]] = data[1]


#for record in input_dict.items():
#    case_id, input_list = record
#    u = mda.Universe(*input_list)
#    detect_distance = 4
#    selection = f'protein and around {detect_distance} resname {case_id2ligand_name[case_id]}'
#    out_prefix = '/media/muwang/新加卷/muwang/work/analysis/GJZX_ACLY/20240526/data/ligand_contact'
#    out_csv = os.path.join(out_prefix,f"contact_{case_id}_{case_id2ligand_name[case_id]}.csv")
#    contact_count(u,selection,out_csv=out_csv)

with Pool(6) as pool:
    pool.map(contact_count_for_pool, input_dict.items())
