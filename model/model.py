import pandas as pd
import os
import cv2
import config 
from data_preprocessing import spaces_in_filenames, nd2_to_tiff
from signal_postprocessing import (to_csv,
                                   update_structure_weights,
                                   calibrate_expression)
from utils import (acronymn_to_id, 
                   atlas_registration,
                   generate_alignment_matrix,
                   process_image,
                   process_mask,
                   cells_in_ROI,
                   download_overlay)
os.chdir(config.root_directory)


def process(base_dir,structure_acronymns):
    # preprocess dataset 
    base_dir = os.path.join(os.path.dirname(config.root_directory),base_dir)
    spaces_in_filenames(base_dir)
    R  = nd2_to_tiff(base_dir)
    data_path = os.path.join(base_dir,R)
    structures = acronymn_to_id(structure_acronymns)
    
    # brain registration
    registration = atlas_registration(data_path,ensemble=True,index_order=False,index_spacing=False,section_thickness = None)

    # iterative base 
    result = []
    im_n = 1
    num = [f for f in os.listdir(data_path) if f.lower().endswith(('.jpeg', '.tif', '.jpg', '.png'))]   # count
    
    # output paths
    save = os.path.join(data_path,f'output')


    os.makedirs(save, exist_ok=True)

    for filename in os.listdir(data_path):
        path = os.path.join(data_path,filename)

        if os.path.isfile(path) and path.lower().endswith(('.jpeg', '.tif', '.jpg', '.png')):
            print(f'processing image {im_n}/{len(num)} --> {filename}')
            image = cv2.imread(path)
            cell_yx = process_image(path)
            alignment = generate_alignment_matrix(registration,filename)
            
            # iteratively process structures 
            image_result = {'Filename': filename}
            for i in range(len(structures)):
                id,ac = structures[i],structure_acronymns[i]
                count = 0
                _, mask_yx,_ = process_mask(id, alignment,image)
                local_yx, count = cells_in_ROI(mask_yx,cell_yx,threshold=9)
                image_result[f'mask_{id}'] = count
                # download_overlay(save,filename,mask_yx,local_yx,ac,image)
            
            result.append(image_result)
            im_n+=1
    
    to_csv(result,os.path.join(config.output_directory,f'{data_path}/result_raw.csv'),acronym_map = config.acronym_map)
    structure_weights = update_structure_weights(result)
    calibrated_result = calibrate_expression(result,structure_weights)
    to_csv(calibrated_result,os.path.join(config.output_directory,f'{data_path}/result_norm.csv'),acronym_map = config.acronym_map)
