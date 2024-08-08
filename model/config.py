import os

# directories
root_directory = '/Your/Root/Directory/Here'
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
names_file_path = os.path.join(root_directory,'mcc/names.txt')

def id_to_name(id):
    return tree.get_structures_by_name([f'{name_map[id]}'])

def id_to_acronymn(id):
    return tree.get_structures_by_name([f'{name_map[id]}'])[0]['acronym']




