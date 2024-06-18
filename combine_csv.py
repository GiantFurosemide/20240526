import glob
from logging_util import logging
from csv_util import combine_csv

data_in_dir = 'data/ligand_contact/distance_4_5_allligand'
data_out_dir = 'data/ligand_contact/distance_4_5_allligand_combine'

ligand_name_list = ['AA1',
                    'AA2',
                    'AAM',
                    'AAT',
                    'AAA',
                    'AAU',
                    'AAK',
                    'AAP',
                    'AAE',
                    'AAL',
                    'AAY',
                    'AAJ',
                    'AAW',
                    'AAI',
                    'AAN',
                    'AAX',
                    'AAG',
                    'AAV',
                    'AA0',
                    'AAQ',
                    'AAS',
                    'AAC',
                    'AAB',
                    'AAD',
                    'AAZ',
                    'AAR',
                    'AAO',
                    'AAH',
                    'AAF',
                    'GX1',
                    ]

for i in ligand_name_list:
    csv_files = glob.glob(f'{data_in_dir}/*{i}.csv')   
    if len(csv_files) == 0:
        continue
    combine_csv(csv_files, f'{data_out_dir}/combined_{i}.csv')
    logging(f'combine {csv_files} --> {data_out_dir}/combined_{i}.csv', logfile=__file__+'.log')


# for AA1 AA2 AAM
data_in_dir = 'data/ligand_contact/distance_4_5_allligand_combine'
data_out_dir = 'data/ligand_contact/distance_4_5_allligand_combine_groupby_lowIC50'   
ligand_name_list = 'AA1 AA2 AAM'

csv_files_in = []
for i in ligand_name_list.split():
    csv_files += glob.glob(f'{data_in_dir}/combined_{i}.csv') 
tmp_name = '_'.join(ligand_name_list.split())
combine_csv(csv_files, f"{data_out_dir}/combined_{tmp_name}.csv")
logging(f'combine {csv_files} --> {data_out_dir}/combined_{tmp_name}.csv', logfile=__file__+'.log')

# for AAR AAO
data_in_dir = 'data/ligand_contact/distance_4_5_allligand_combine'
data_out_dir = 'data/ligand_contact/distance_4_5_allligand_combine_groupby_highIC50' 
ligand_name_list = 'AAR AAO'

csv_files_in = []
for i in ligand_name_list.split():
    csv_files += glob.glob(f'{data_in_dir}/combined_{i}.csv') 
tmp_name = '_'.join(ligand_name_list.split())
combine_csv(csv_files, f"{data_out_dir}/combined_{tmp_name}.csv")
logging(f'combine {csv_files} --> {data_out_dir}/combined_{tmp_name}.csv', logfile=__file__+'.log')
