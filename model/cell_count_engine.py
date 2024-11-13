import os
import cv2
import argparse
import model.config as config
from model.data_preprocessing import spaces_in_filenames, nd2_to_tiff
from model.signal_postprocessing import (to_csv,update_structure_weights,calibrate_expression,plot_brain_region_histogram)
from utils import (acronymn_to_id, 
                   atlas_registration,
                   generate_alignment_matrix,
                   process_image,
                   process_mask,
                   cells_in_ROI,
                   download_overlay)

os.chdir(config.root_directory)



def main():
    parser = argparse.ArgumentParser(description='Regional Cell Count of Brain Dataset')
    parser.add_argument('--input_dir', type=str, required=True, help='Input directory containing the png brain scans')
    parser.add_argument('--structures', nargs='+', required=True, help='List of structure acronyms (e.g., VPL, BLA, PVT) or "FULL_BRAIN" to include all currently supported structures')
    args = parser.parse_args()
    if args.structures == ['FULL_BRAIN']: args.structures = config.FULL_BRAIN
    if args.structures == ['LIGHT_BRAIN']: args.structures = config.LIGHT_BRAIN
    process(args.input_dir, args.structures)


def process(base_dir,structure_acronymns,overlay=False):
    # preprocess dataset 
    if(structure_acronymns == "FULL_BRAIN"): structure_acronymns = FULL_BRAIN
    base_dir = os.path.join(os.path.dirname(config.root_directory),base_dir)
    spaces_in_filenames(base_dir)
    R  = nd2_to_tiff(base_dir)
    data_path = os.path.join(base_dir,R)

    # brain registration
    structures = acronymn_to_id(structure_acronymns)
    registration = atlas_registration(data_path,ensemble=True,index_order=False,index_spacing=False,section_thickness = None)

    # iterative base, output paths
    result = []
    im_n = 1
    num = [f for f in os.listdir(data_path) if f.lower().endswith(('.jpeg', '.tif', '.jpg', '.png'))]   # count
    save = os.path.join(data_path,f'output')

    for filename in os.listdir(data_path):
        path = os.path.join(data_path,filename)

        if os.path.isfile(path) and path.lower().endswith(('.jpeg', '.tif', '.jpg', '.png')):
            print(f'processing image {im_n}/{len(num)} --> {filename}')
            image = cv2.imread(path)
            cell_yx = process_image(path)
            alignment = generate_alignment_matrix(registration,filename)
            
            image_result = {'Filename': filename}
            for i in range(len(structures)): # iteratively process structures 
                id,ac,count = structures[i],structure_acronymns[i],0
                _, mask_yx,_ = process_mask(id, alignment,image)
                local_yx, count = cells_in_ROI(mask_yx,cell_yx,threshold=9)
                image_result[f'mask_{id}'] = count
            
            result.append(image_result)
            im_n+=1
    
    
    # exporting results     
    fn_1 = os.path.join(config.root_directory,f'{data_path}/result_raw.csv')
    fn_2 = os.path.join(config.root_directory, f'{data_path}/result_norm.csv')
    to_csv(result,fn_1,acronym_map = config.acronym_map)
    
    structure_weights = update_structure_weights(result)
    calibrated_result = calibrate_expression(result,structure_weights,stabilizing_parameter=0.03)
    to_csv(calibrated_result,fn_2,acronym_map = config.acronym_map)
    plot_brain_region_histogram(fn_2,structure_acronymns)
    return result

if __name__ == "__main__":
    main()