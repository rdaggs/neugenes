import os

# directories
root = 'machine/root/directory'
root_directory = os.path.join(root,'model')
root_directory_new = root
output_directory = os.path.join(root_directory,'model')

# AllenSDK routing
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
manifest_file_path = os.path.join(root_directory,'mcc/mouse_connectivity_manifest.json')
annotation_file_path = os.path.join(root_directory,'mcc/annotation/ccf_2017/structure_masks/resolution_25')
mask_weights = os.path.join(root_directory,'mcc/mask_weights.json')
mcc = MouseConnectivityCache(manifest_file=manifest_file_path, resolution=25)
reference_space =  mcc.get_reference_space()
tree = mcc.get_structure_tree()
acronym_map = tree.value_map(lambda x: x['id'], lambda y: y['acronym'])
name_map = tree.get_name_map() 

# Brain atlas structure query 
names_file_path = os.path.join(root,'model','mcc/names.txt')

def id_to_name(id):
    return tree.get_structures_by_name([f'{name_map[id]}'])

def id_to_acronymn(id):
    return tree.get_structures_by_name([f'{name_map[id]}'])[0]['acronym']



# test dataset path 

LIGHT_BRAIN = ['SSp-ll', 'SSp-tr', 'VISC', 'AUDpo', 'VISal', 'VISam', 'VISl', 'VISp', 'VISpl', 
               'VISpm', 'VISli', 'VISpor', 'ACAd', 'ACAv', 'PL', 'ILA', 'ORB', 'VISrl', 'TEa', 'AId', 
               'AIp', 'AIv', 'RSPagl', 'RSPd', 'RSPv', 'PTLp', 'PERI', 'ECT', 'TT', 'DP', 'PIR', 
               'NLOT', 'COA', 'PAA', 'TR', 'HPF', 'FC', 'ENTl', 'ENTm', 'PAR', 'POST', 'PRE', 'SUB', 
               'ProS', 'CLA', 'EP', 'LA', 'BMA', 'ACB', 'FS', 'LSc', 'LSr', 'LSv', 'SH', 'sAMY', 'MEA', 
               'GPe', 'GPi', 'SI', 'MA', 'NDB', 'SPF', 'PP', 'MG', 'SGN', 'AM', 'AD', 'MD', 'PT', 'RE', 
               'RH', 'PCN', 'SubG', 'SO', 'ASO', 'PVH', 'ARH', 'ADP', 'AHA', 'MEPO', 'PS', 'SCH', 'VMH', 
               'TU', 'ZI', 'VTA', 'RN', 'SOC', 'TRN', 'LDT', 'CN', 'PARN', 'RPA', 'DN', 'PVT', 'BLA']

FULL_BRAIN = [  'SSp-ll', 'SSp-tr', 'VISC', 'AUDpo', 'VISal', 'VISam', 'VISl', 'VISp', 'VISpl', 'VISpm', 
                'VISli', 'VISpor', 'ACAd', 'ACAv', 'PL', 'ILA', 'ORB', 'VISrl', 'TEa', 'AId', 'AIp', 'AIv', 
                'RSPagl', 'RSPd', 'RSPv', 'PTLp', 'PERI', 'ECT', 'TT', 'DP', 'PIR', 'NLOT', 'COA', 'PAA', 'TR', 
                'HPF', 'FC', 'ENTl', 'ENTm',  'PAR', 'POST', 'PRE', 'SUB', 'ProS', 'CLA', 'EP', 'LA', 'BMA', 
                'ACB', 'FS', 'LSc', 'LSr', 'LSv', 'SH', 'sAMY', 'MEA', 'GPe', 'GPi', 'SI', 'MA', 'NDB',  'SPF', 
                'PP', 'MG', 'SGN', 'AM', 'AD', 'MD', 'PT', 'RE', 'RH', 'PCN', 'SubG', 'SO', 'ASO', 'PVH', 'ARH', 
                'ADP', 'AHA', 'MEPO', 'PS', 'SCH', 'VMH', 'TU', 'ZI', 'VTA', 'RN', 'SOC', 'TRN', 'LDT', 'CN', 
                'PARN', 'RPA', 'DN', 'PVT', 'BLA', 'ARH', 'DMH', 'ENT', 'FRP', 'HPF', 'HY', 'ILA', 'LHA', 'LS', 
                'MPO', 'MPN', 'MD', 'MBmot', 'ACB', 'NTS', 'PB', 'PVH', 'PeF', 'PIR', 'PH', 'PL', 'SSp', 'VISp', 
                'RSP', 'MOs', 'SFO', 'SCH', 'VMH', 'VTA', 'TH', 'VAL', 'VP', 'MG', 'LGd', 'LP', 'PO', 'POL', 'VM', 
                'VPL', 'VPM', 'PVT', 'AV', 'AM', 'AD', 'IAM', 'IAD', 'LD', 'IMD', 'MTN', 'PCN', 'CL', 'PF', 'RE', 
                'Xi', 'ILM', 'RH', 'CM', 'SGN', 'GENd', 'GENv', 'VPLpc', 'VPMpc', 'LGd-sh', 'LGd-co', 'LGd-ip', 
                'LGv', 'IntG', 'SPFm', 'SPFp', 'SPA', 'PP', 'ATN', 'APN', 'NOT', 'OP', 'PPT', 'RPF', 'CUN', 'III', 
                'IV', 'VI', 'XII', 'I5', 'IF', 'IPR', 'PRC','VTN']






