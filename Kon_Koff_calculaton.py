import os
import logging
import glob

import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import numpy as np
import pandas as pd

from json_util import read_json
from csv_util import csv2df,records2csv



def get_logger():
    """
    get logger
    """
    logging.basicConfig(filename=f'{__file__}.log',
                        level=logging.DEBUG,
                        format='[%(asctime)s] - %(name)s - %(levelname)s : %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filemode='w')
    logger = logging.getLogger()
    return logger


def frame_index_list(df: pd.DataFrame, nr_total_frames: int) -> list:
    # Initialize a list of zeros with length nr_total_frames
    output_list = [0]*nr_total_frames

    # Set the value to 1 for indices that exist in df['frame_index']
    for index in df['frame_index']:
        if index < nr_total_frames:
            output_list[index] = 1

    return np.array(output_list)


def calculate_KonKoffKD(u: mda.Universe, 
                        ligand_name: str, 
                        n_ligands: int, 
                        temperature: float, 
                        csv_in: str, 
                        case_id:str,
                        bound_durations: list,
                        unbound_durations: list):
    # 计算Kon和Koff
    logger.info(f"Calculating Kon and Koff for {ligand_name}")
    total_time = len(u.trajectory) * u.trajectory.dt * 1e-12  # 总模拟时间 (s)
    NA = 6.022e23  # 阿伏伽德罗常数 (mol^-1)
    concentration_ligands = (n_ligands/NA)/ (u.dimensions[0] * u.dimensions[1] * u.dimensions[2] * 1e-27)  # 配体浓度 (M=mol/L) (1 Å^3 = 1e-27 L) (M)
    kon = len(bound_durations) / (total_time * concentration_ligands)  # 结合事件数除以总时间和配体数 (M^-1s^-1)
    koff = 1 / (np.mean(unbound_durations) * u.trajectory.dt * 1e-12)  # 解离速率常数 (s^-1)
    kD = koff / kon  # 解离常数 (M)
    delta_G = np.log(kD) * 8.3144598 * temperature # 自由能 delta_G =RTln(KD) (J/mol) (R=8.314 J/(mol K)) 
    
    # write out kon koff KD to csv
    analysis_out =  "_".join(csv_in.replace('.csv','').split('_')[:-1]) + '_ALL_KonKoffKD.csv'
    records = [{'ligand_name': ligand_name,'case_id':case_id,'raw_data_csv':csv_in,'kon':kon,'koff':koff,'kD':kD ,'delta_G':delta_G}]
    records2csv(records,analysis_out)
    logger.info(f"wrote out {analysis_out}")
    
    # write out kon koff KD to txt
    analysis_out = "_".join(csv_in.replace('.csv','').split('_')[:-1]) + '_ALL_KonKoffKD.txt'
    record = records[0]
    with open(analysis_out, 'w') as f:
        f.write(f"case_id:{record['case_id']}\n")
        f.write(f"ligand_name:{record['ligand_name']}\n")
        f.write(f"raw_data_csv:{record['raw_data_csv']}\n")
        f.write(f"kon(M^-1s^-1):{record['kon']:.6e}\n")
        f.write(f"koff(s^-1):{record['koff']:.6e}\n")
        f.write(f"kD(M):{record['kD']:6e}\n")
        f.write(f"delta_G(J):{record['delta_G']:.6e}\n")
    logger.info(f"wrote out {analysis_out}")


###################### main #######################

logger = get_logger()
data = read_json('data_all_ligand_name.json')

for data_record in data:
    
    top = data_record['out_gro']
    traj_in = data_record['out_traj']
    ligand_name = data_record['ligand_name']
    data_for_u = [top,traj_in]
    case_id = data_record['case_id']
    data_prefix = 'data/ligand_contact/distance_4_5_eachligand'
    temperature = 303.15 # 温度 (K)

    # 加载模拟轨迹
    logger.info(f"Loading trajectory {traj_in}")
    u = mda.Universe(*data_for_u)
    nr_frames = len(u.trajectory)
    logger.info(f"Done loading trajectory {traj_in}")

    # 定义受体和配体的选择
    #receptor = u.select_atoms('protein')
    ligands = u.select_atoms(f'resname {ligand_name}')  # 选择所有配体

    n_ligands = len(ligands.residues)


    # 距离阈值（单位：Å）
    # cutoff = 4.5

    # 初始化变量以记录结合和解离事件的持续时间
    all_bound_durations = []
    all_unbound_durations = []

    # 分别计算每个配体的接触情况
    csv_files = glob.glob(f'{data_prefix}/*{case_id}_{ligand_name}_????.csv')
    
    if csv_files:
        pass
    else:
        continue

    for csv in csv_files:
        #ligand = ligands.residues[i].atoms
        
        # !!!! contact.Contacts() 无法使用，需要修改 
        ## 计算接触次数
        #logging(f"Calculating contacts for ligand {i + 1} of {n_ligands}")
        #contact_analysis = contacts.Contacts(u, select=(receptor, ligand), refgroup=(receptor, ligand), radius=cutoff)
        #contact_analysis.run()
        #logging(f"Done calculating contacts for ligand {i + 1} of {n_ligands}")

        # 提取接触时间序列
        # timeseries = contact_analysis.timeseries[:, 1]
        logger.info(f"Loading csv {csv}")
        df = csv2df(csv)
        time_series = frame_index_list(df,nr_frames)
        logger.info(f"Done loading csv {csv}")
        # 分析结合和解离事件
        bound_states = time_series > 0
        # 结合状态变化点
        bound_changes = np.where(np.concatenate(([bound_states[0]], bound_states[:-1] != bound_states[1:], [True])))[0]
        # 非结合状态变化点
        unbound_changes = np.where(np.concatenate(([~bound_states[0]], ~bound_states[:-1] != ~bound_states[1:], [True])))[0]

        # 计算结合和解离的持续时间
        bound_duration = np.diff(bound_changes)[::2]
        unbound_duration = np.diff(unbound_changes)[::2]
        
        ###################### for each ligands #######################

        calculate_KonKoffKD(u, ligand_name, n_ligands, temperature, csv, case_id, bound_duration, unbound_duration)
        ###################### for each ligands #######################
        
        # 记录所有配体的持续时间
        all_bound_durations.extend(bound_duration)
        all_unbound_durations.extend(unbound_duration)
        logger.info(f"Done calculating durations for ligand {csv}")
    

    ###################### for all ligands #######################

    calculate_KonKoffKD(u, ligand_name, n_ligands, temperature, csv, case_id, all_bound_durations, all_unbound_durations)
    ###################### for all ligands #######################


