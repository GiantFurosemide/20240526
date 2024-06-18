"""
This script extracts protein ligand from trajectory file,
and saves it in a new gro file and new xtc/trr file.
"""

import json
import glob
import os
import pandas as pd
import MDAnalysis as mda
from csv_util import *
from json_util import *


def extract_selection_to_gro_traj(gro_in:str, traj_list_in:list, selection_str:str ,out_gro:str, out_traj:str=None):
    """
    extract selected atom from trajectory file and save it in a new gro file and trajectory file.
    :param gro_in: input gro file
    :param traj_list_in: a list of input trajectory files
    :param selection_str: a selection string
    :param out_gro: output gro file
    :param out_traj: output trajectory file
    :return: none
    """

    print(f"input gro file:\n> {gro_in}\n input trajectory list:\n> {traj_list_in}")
    print(f"selection prompt is:\n> {selection_str}")
    print(f"out new gro file:\n> {out_gro}")
    print(f"out new trajectory:\n> {out_traj}")

    u = mda.Universe(*[gro_in, *traj_list_in])
    selection = u.select_atoms(selection_str)
    selection.write(out_gro)

    if out_traj is not None:
        with mda.Writer(out_traj, selection.n_atoms) as W:
            for ts in u.trajectory:
                W.write(selection)
                if int(ts.frame) % 200 == 0:
                    print(f"writing trajectory {ts.frame} (1 ns/frame)...")

    print(f"input gro file:\n> {gro_in}\n input trajectory list:\n> {traj_list_in}")
    print(f"selection prompt is:\n> {selection_str}")
    print(f"out new gro file:\n> {out_gro}")
    print(f"out new trajectory:\n> {out_traj}")
    print("done")


def extract_selection_to_gro_of_frame(gro_in:str, traj_list_in:list, selection_str:str ,out_gro:str, frame:int=0):
    """
    extract selected atom from trajectory file of FRAME N and save to a new gro file.
    frame 0 is the first frame of the trajectory file.
    :param gro_in: nput gro file
    :param traj_list_in: a list of input trajectory files
    :param selection_str: a selection string
    :param out_gro: output gro file
    :param frame: frame to extract
    :return:
    """
    print(f"input gro file:\n> {gro_in}\n input trajectory list:\n> {traj_list_in}")
    print(f"selection prompt is:\n> {selection_str}")
    print(f"out new gro file:\n> {out_gro}")
    print(f"frame to extract:\n> {frame}")

    u = mda.Universe(*[gro_in, *traj_list_in])
    selection = u.select_atoms(selection_str)
    #selection.write(out_gro)
    number_traj = len(u.trajectory)
    if frame < 0:
        out_gro = out_gro.replace('.gro', f"_{frame+number_traj}.gro")
    else:
        out_gro = out_gro.replace('.gro', f"_frame{frame}.gro")
    print(f"frame to extract:\n> {frame}")

    u.trajectory[frame]
    selection.write(out_gro)

    print(f"input gro file:\n> {gro_in}\n input trajectory list:\n> {traj_list_in}")
    print(f"selection prompt is:\n> {selection_str}")
    print(f"out new gro file:\n> {out_gro}")
    print(f"frame to extract:\n> {frame}")
    print("done")


if __name__ == "__main__":
    raw_data_list = read_json('data_all_ligand_name.json')
    # raw_data_list : [{'top':, 'gro':, 'traj':, 'case_id','ligand_name':,'case_id': ,'out_gro':,'out_traj': },...]

    # round 1 extract new protein and ligand
    for record in raw_data_list:
        gro_in = record['gro']
        traj_list_in = record['traj']
        selection_str = f"protein or resname {record['ligand_name']}"
        out_gro = record['out_gro']
        #out_traj = record['out_traj']
        out_traj = None
        extract_selection_to_gro_traj(gro_in, traj_list_in, selection_str, out_gro, out_traj)

#    # round 2 extraxt last frame
#    for record in raw_data_list:
#        gro_in = record['gro']
#        traj_list_in = record['traj']
#        selection_str = "protein"
#        out_gro = record['out_gro']
#        out_gro= out_gro.replace('.gro', "_SEL_protein.gro")
#        frame=-1
#
#        extract_selection_to_gro_of_frame(gro_in, traj_list_in, selection_str, out_gro, frame=frame)
