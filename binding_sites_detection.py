import MDAnalysis as mda
from sklearn.cluster import DBSCAN
import numpy as np

# 任务一
def get_coordinates(u, resid):
    protein = u.select_atoms('protein')
    ca = protein.select_atoms('name CA and resid {}'.format(resid))
    side_chain = protein.select_atoms('not name CA C N O and resid {}'.format(resid))
    ca_coord = ca.center_of_mass()
    side_chain_coord = side_chain.center_of_mass()
    return ca_coord, side_chain_coord

# 任务二
def cluster_coordinates(coords, eps=1.0, min_samples=5):
    clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(coords)
    labels = clustering.labels_
    return labels

# 任务三
def calculate_weight_sum(labels, weights):
    unique_labels = np.unique(labels)
    group_weights = []
    for label in unique_labels:
        group_weight = np.sum(weights[labels == label])
        group_weights.append((label, group_weight))
    group_weights.sort(key=lambda x: x[1], reverse=True)
    return group_weights

# 使用示例
u = mda.Universe('protein.pdb')  # 你的PDB文件
coords = []
weights = []
resids = []  # 你的resid列表
for resid in resids:
    ca_coord, side_chain_coord = get_coordinates(u, resid)
    coords.append(ca_coord)
    coords.append(side_chain_coord)
    weights.append(1.0)  # 你的weight
    weights.append(1.0)  # 你的weight
labels = cluster_coordinates(coords)
group_weights = calculate_weight_sum(labels, weights)
for group, weight in group_weights:
    print('Group: {}, Weight: {}'.format(group, weight))
