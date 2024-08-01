import MDAnalysis as mda
from sklearn.cluster import DBSCAN
import numpy as np

# 任务一
def get_coordinates(u, resid):
    protein = u.select_atoms('protein')
    ca = protein.select_atoms('name CA and resnum {}'.format(resid))
    side_chain = protein.select_atoms('not name CA C N O and resnum {}'.format(resid))
    ca_coord = ca[0].position # first atom
    side_chain_coord = side_chain[0].position # first atom
    return ca_coord, side_chain_coord

# 任务二
def cluster_coordinates(coords, eps=0.5, min_samples=5):
    clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(coords)
    labels = clustering.labels_
    return labels

# 任务三
def calculate_weight_sum(labels, weights: np.ndarray):
    unique_labels = np.unique(labels)
    group_weights = []
    for label in unique_labels:
        group_weight = np.sum(weights[labels == label])
        group_weights.append((label, group_weight))
    group_weights.sort(key=lambda x: x[1], reverse=True)
    return group_weights



#使用示例
u = mda.Universe('data/structure/step4.1_equilibration.gro')  # 你的PDB文件
from csv_util import csv2df
df = csv2df('data/ligand_contact/distance_4_5_allligand_combine_groupby_lowIC50_analysis/combined_AA1_AA2_AAM_count_bin10.csv')



coords = [] # in Angstrom
weights = []   # 你的weight
resids = df['resnum'].values  # 你的resid列表

for i in df['count'].values :
    weights.append(i)
    #weights.append(i)

for resid in resids:
    try:
        ca_coord, side_chain_coord = get_coordinates(u, resid)
        coords.append(ca_coord)
        #coords.append(side_chain_coord)
    except ValueError as e :
        print('Residue {} not found'.format(resid))
        continue


labels = cluster_coordinates(coords, eps=4.5, min_samples=5)

group_weights = calculate_weight_sum(labels, np.array(weights))
for group, weight in group_weights:
    print('Group: {}, Weight: {}'.format(group, weight))

