import json
import glob
import os
import pandas as pd
import MDAnalysis as mda
from csv_util import *
from json_util import *
def make_residue_df(selected_atom_group:mda.AtomGroup,
                    trajectory_index:int,
                    df:pd.DataFrame=None,
                    df_columns:tuple=None):
    """
    record contact information of selected atom group FOR ONE frame(define by trajectory_index)
    :param selected_atom_group: mda.AtomGroup
    :param trajectory_index: fram index of selected atom group
    :param df: input dataframe for further concatenate
    :param df_columns: keys ('global_resindex', 'resname', 'resnum','frame_index') for dataframe
    :return: df with contact information of selected atom group
    """

    if df_columns is None:
        df_columns = ('global_resindex', 'resname', 'resnum','frame_index')
    nr_residues = len(selected_atom_group.residues)
    for i in range(nr_residues):

        resname = selected_atom_group.residues[i].resname
        resindex = selected_atom_group.residues[i].resindex
        resnum = selected_atom_group.residues[i].resnum
        #print(f" {resindex} {resname} {resnum} {trajectory_index}")
        new_data = pd.DataFrame([[resindex,resname,resnum,trajectory_index]],columns=df_columns)
        df = pd.concat([df,new_data])
    return df


def contact_count(u:mda.Universe,
                 selected_atom_group:str="protein and around 3.5 resname OXY",
                 df_columns:tuple=('global_resindex', 'resname', 'resnum','frame_index'),
                 out_csv='ligand_contact_output.csv',
                 )->pd.DataFrame:
    """
    main function for contact count of selected atom group for trajectory and save to csv
    :param u: mda.Universe object
    :param selected_atom_group: selection prompt for selected atom group
    :param df_columns: keys ('global_resindex', 'resname', 'resnum','frame_index') for dataframe
    :param out_csv: output csv file name
    :return: df with contact information of selected atom group
    """
    df = pd.DataFrame(columns=df_columns)
    for ts in u.trajectory:
        protein_contact_ligand = u.select_atoms(selected_atom_group)
        df = make_residue_df(protein_contact_ligand,ts.frame,df=df,df_columns=df_columns)
        print(f'extracting ... {protein_contact_ligand.n_residues}records of frame {ts.frame}')

    df.to_csv(out_csv,index=False)
    print(f"contact data saved to {out_csv}")
    return df

def contact_count_for_pool(record):
    case_id, input_list = record
    print(f'loading university {case_id}(case_id)')

    raw_data_list = read_json('data_all_ligand_name.json')
    case_id2ligand_name = {}
    for rec in raw_data_list:
        case_id2ligand_name[rec['case_id']] = rec['ligand_name']


    u = mda.Universe(*input_list)
    detect_distance = 4
    selection = f'protein and around {detect_distance} resname {case_id2ligand_name[case_id]}'
    out_prefix = '/media/muwang/新加卷/muwang/work/analysis/GJZX_ACLY/20240526/data/ligand_contact'
    out_csv = os.path.join(out_prefix,f"contact_{case_id}_{case_id2ligand_name[case_id]}.csv")
    print(f"> case id: {case_id}")
    contact_count(u,selection,out_csv=out_csv)

if __name__ == '__main__':

##### round1: contact analysis from traj of 'all'

#    raw_data_list = read_json('data_all_ligand_name.json')
#
#
#    case_id2ligand_name = {}
#    for rec in raw_data_list:
#        case_id2ligand_name[rec['case_id']] = rec['ligand_name']
#
#
#    input_list = [(rec['case_id'],[rec['gro']]+rec['traj']) for rec in raw_data_list]
#    input_dict = {}
#    for data in input_list:
#        input_dict[data[0]] = data[1]
#
#
#    for record in input_dict.items():
#        case_id, input_list = record
#        u = mda.Universe(*input_list)
#        detect_distance = 4
#        selection = f'protein and around {detect_distance} resname {case_id2ligand_name[case_id]}'
#        out_prefix = '/media/muwang/新加卷/muwang/work/analysis/GJZX_ACLY/20240526/data/ligand_contact'
#        out_csv = os.path.join(out_prefix,f"contact_{case_id}_{case_id2ligand_name[case_id]}.csv")
#        contact_count(u,selection,out_csv=out_csv)


##### round2: contact analysis from traj of 'protein o ligand'
    raw_data_list = read_json('data_all_ligand_name.json')
    # raw_data_list : [{'top':, 'gro':, 'traj':, 'case_id','ligand_name':,'case_id': ,'out_gro':,'out_traj': },...]
    for record in raw_data_list:
        gro_in = record['gro']
        traj_list_in = record['traj']
        selection_str = f"protein or resname {record['ligand_name']}"
        out_gro = record['out_gro']
        out_traj = record['out_traj']
        ligand_name = record['ligand_name']
        case_id = record['case_id']

        done_list = [('GJZX_ACLY_0001_100w','GX1'),
                     ('GJZX_ACLY_0002_100w','GX1'),
                     ('GJZX_ACLY_0021_100w','AAF'),
                     ('GJZX_ACLY_0023_100','AAH'),
                     ('GJZX_ACLY_0028_100w','AAM'),
                     ('GJZX_ACLY_0030_100w','AAO'),
                     ('GJZX_ACLY_0033_100w','AAR'),
                     ('GJZX_ACLY_0062_100w','AAR')
                     ]


        if (case_id, ligand_name) in done_list:
            continue


        # for each ligand (count contact)
        for res_num in range(1100,1100+24):
            u = mda.Universe(out_gro,out_traj)
            detect_distance = 4.5
            selection = f'protein and around {detect_distance} (resname {ligand_name } and resnum {res_num})'
            out_prefix = '/media/muwang/新加卷/muwang/work/analysis/GJZX_ACLY/20240526/data/ligand_contact/distance_4_5_eachligand'
            out_csv = os.path.join(out_prefix,f"contacteachlig_{case_id}_{ligand_name}_{res_num}.csv")
            contact_count(u,selection,out_csv=out_csv)


    for record in raw_data_list:
        gro_in = record['gro']
        traj_list_in = record['traj']
        selection_str = f"protein or resname {record['ligand_name']}"
        out_gro = record['out_gro']
        out_traj = record['out_traj']
        ligand_name = record['ligand_name']
        case_id = record['case_id']

        # for all ligand (count contact)
        u = mda.Universe(out_gro, out_traj)
        detect_distance = 4.5
        selection = f'protein and around {detect_distance} resname {ligand_name}'
        out_prefix = '/media/muwang/新加卷/muwang/work/analysis/GJZX_ACLY/20240526/data/ligand_contact/distance_4_5_allligand'
        out_csv = os.path.join(out_prefix, f"contacteachlig_{case_id}_{ligand_name}.csv")
        contact_count(u, selection, out_csv=out_csv)
