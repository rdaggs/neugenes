import os 
import pandas as pd
import nrrd
import json
import numpy as np
import matplotlib.pyplot as plt
import model.config as config
import model.utils as utils
os.chdir(config.root_directory)


def to_csv(data, output_path, acronym_map):
    """
        Find cells within the region of interest (ROI) defined by mask coordinates.

        Args:
            data: Coordinates of the mask.
            output_path: Coordinates of the cells.
            acronym_map: Distance threshold for considering cells within the ROI.

        Returns:
            CSV containing each images quantification within the ROI 
    """
    rows = []

    for entry in data:
        filename = entry['Filename']
        areas = {key: value for key, value in entry.items() if key != 'Filename'}
        structure_ids = [int(key.split('_')[1]) for key in areas.keys()]
        
        try:
            acronyms = [acronym_map[structure_id] for structure_id in structure_ids]
        
        except KeyError as e:
            print(f"Structure ID {e.args[0]} not found in acronym_map. Skipping entry.")
            continue

        row = {'Filename': filename}
        for acronym, activation in zip(acronyms, areas.values()):
            row[acronym] = activation

        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)

def plot_brain_region_histogram(filepath, structures):
    data = pd.read_csv(filepath)
    brain_regions = data.columns[1:]
    region_sums = data[structures].sum()
    
    scaled_region_sums = (region_sums / region_sums.max()) * 100
    sorted_regions = scaled_region_sums.sort_values(ascending=False)

    plt.figure(figsize=(18, 8))
    plt.bar(sorted_regions.index, sorted_regions.values, color='blue')
    plt.title('Regional Expression Levels Calibrated', fontsize=14)
    plt.ylabel('Calibrated Expression Intensity (Scaled to 100)', fontsize=12)
    plt.xlabel('Brain Region (as defined in Allen CCFv3)', fontsize=12)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(filepath.replace(".csv", "") + "_histogram.png")

def update_structure_weights(data):
    """
        Create a weight file representative of each masks voxels. This will
        be used for calibrating intensity relative to size

        Args:
            data: resulting 

        Returns:
            mask_sums: weight file 
    """

    
    # check if the weight file exists; if not, initialize it with an empty dictionary
    if not os.path.exists(config.mask_weights):
        with open(config.mask_weights, 'w') as json_file:
            json.dump({}, json_file)
        print(f'Initialized an empty JSON file at {config.mask_weights}')

    # existing data from the JSON file
    with open(config.mask_weights, 'r') as json_file:
        try:
            mask_sums = json.load(json_file)
        except json.JSONDecodeError:
            mask_sums = {}
    
    # gather ids
    ids = set()
    for record in data:
        for key in record.keys():
            if key.startswith('mask_'):
                numeric_id = ''.join([char for char in key if char.isdigit()])
                ids.add(numeric_id)
    ids = [int(id) for id in ids]

    # check if exists
    for id in ids:
        if str(id) in mask_sums:
            continue
        
        try:
            path = os.path.join(config.annotation_file_path, f'structure_{id}.nrrd')
            mask, header = nrrd.read(path)
            sum_val = np.sum(mask == 1)
            mask_sums[str(id)] = int(sum_val)
        
        except FileNotFoundError:
            print(f'downloading unregistered mask {id}. might take a moment')
            utils.gen_mask(id)
            path = os.path.join(config.annotation_file_path, f'structure_{id}.nrrd')
            mask, _ = nrrd.read(path)
            sum_val = np.sum(mask == 1)
            mask_sums[str(id)] = int(sum_val)
        
        except Exception as e:
            print(f'external error occurred: {e}...consider re-establishing connection with AllenSDK API')

    with open(config.mask_weights, 'w') as json_file:
        json.dump(mask_sums, json_file)

    return mask_sums


def calibrate_expression(data,structure_weights):
    """
        To calibrate expression density, it is important to consider the
        size of the region when reviewing its cell count. This function takes
        expression data and returns an adjusted version relative to the size 
        of each structure

        Args:
            data: output of inference model 
            structure_weights: size of each structure in voxels

        Returns:
            [calibrat3d] output from inference 
    """
    # normalize the weights
    min_val = min(structure_weights.values())
    max_val = max(structure_weights.values())
    normalized_weights = {key: (value - min_val) / (max_val - min_val) * 1 for key, value in structure_weights.items()}

    # update expression values
    for x in data:
        for key in x.keys():
            if key.startswith('mask_'):
                id = ''.join([char for char in key if char.isdigit()])
                w = normalized_weights[id]
                # this transformation is meant to accentuate smaller regions
                if w != 0:
                    # x[key] = x[key]/(w + 0.27)
                    x[key] /= w + 0.05
                else:
                    # break if division error is posed 
                    continue
    
    # feedback updated data
    return data