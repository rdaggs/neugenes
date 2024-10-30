import pandas as pd
import os
import cv2
import numpy as np
import nrrd
import ipywidgets as widgets
import math
import model.config as config
from PIL import Image
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from DeepSlice import DSModel
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from allensdk.core.reference_space import ReferenceSpace
import functools
from model.data_preprocessing import match_aspect_ratio

from requests.exceptions import HTTPError
os.chdir(config.root_directory)

def get_ids():
    """
        Gathering acronyms and converting corresponding AllenSDK 
        structure ids for mask generation
    """
    ac_to_id = {}
    with open(config.names_file_path,'r') as file:
        for line in file:
            parts = line.split()
            if len(parts) >= 2:

                acronym = parts[0]
                structure_number = parts[1]
                ac_to_id[acronym] = structure_number
    
    return ac_to_id

def validate_id(id):
    """
        In cooperation with neuroanatomists, some structures were deemed incorrectly labeled 
        or generated and this function aims to catch some of these errors. If any of these 
        decisions plauge your dataset, feel free to reach out. I have mutated nrrd files 
        locally so please reach out if you'd like me to email these to you. (stripping Hippocampal 
        formation of its most posterior region to resemble hippocampus)

        Args:
            id: Raw id number  

        Returns:
            id: Corrected or Validated id for AllenSDK
    """

    # adjustment to hippocampal formation
    corrections = {
        # adjustment to hippocampal formation
        1089: 375,
        1080: 375
    }
    
    return corrections.get(id, id)
    
    
def acronymn_to_id(structures_ac):
    ac_to_id = get_ids()

    # list of structures
    if isinstance(structures_ac, list):
        id_list = []

        for acronym in structures_ac:
            if acronym in ac_to_id:
                id_list.append(validate_id(ac_to_id[acronym]))
            else:
                print(f"Functionality for structure corresponding to acronym '{acronym}' is either not supported or not found in the dictionary. Skipping.")
        
        return id_list
    
    # singular structure 
    elif isinstance(structures_ac, str):
        if structures_ac in ac_to_id:
            return validate_id(ac_to_id[structures_ac])
        else:
            print(f"Functionality for structure corresponding to acronym '{structures_ac}' is either not supported or not found in the dictionary. Skipping.")
            return None

    else:
        raise TypeError("Input should be a list of acronyms or a single acronym, not an AllenSDK structure id.")


def id_to_acronym(id):
    """
    Converts a given ID back into its corresponding acronym using the reverse mapping.
    
    Args:
        id (int): The AllenSDK structure ID to be converted.
    
    Returns:
        acronym (str): The corresponding acronym for the given ID.
    """
    ac_to_id = get_ids()
    id_to_ac = {v: k for k, v in ac_to_id.items()}

    if isinstance(id,list):
        r = []
        for x in id:
            r.append(id_to_ac[x])
        return r

    elif isinstance(id,int):
        return id_to_ac[str(id)]
    
    elif isinstance(id,str):
        return id_to_ac[id]

def atlas_registration(dir,ensemble,index_order,index_spacing,section_thickness):
    """

        Args:
            dir: directory where dataset is located 

            index_order: To reorder your sections according to their number. If your section numbers 
                         are the precise index which the sections were cut (ie; 1, 2, 4, indicates 
                         that section 3 has been left out of the series)

            index_spacing: 
            
            micron_spacing: if you know the exact Thickness of your sections in microns you can include 
                            this here as well, but don't worry, if you dont include it we will estimate 
                            and tell you our guess (known-->) enforce_index_spacing(section_thickness = 25)

        Returns:
            anchoring: predctions of 9 anchored Allen CCF coordinates for each image in directory
    """
    model = DSModel('mouse')

    # section thickness 
    if section_thickness: 
        model.predict(dir, ensemble, section_numbers=True)
    
    if not section_thickness: 
        model.predict(dir, ensemble, section_numbers=False)

    model.propagate_angles()
    
    if index_order:
        model.enforce_index_order()
        
    if index_spacing and section_thickness:
        model.enforce_index_spacing(section_thickness)

    model.save_predictions(dir + '/results')
    anchoring = pd.read_csv(os.path.join(dir, 'results.csv')) 
    return anchoring


def process_image(image):
    """
        Gathering the (y,x) location of each expressed cell body in an image 
        for processing against masks 

        Args:
            image: brain scan as a cv2 (numpy.ndarray) object 
                   or a previously 

        Returns:
            activation_locations: a coordinate-wise mapping of each cell 
                                  activation in brain scan
    """
    if isinstance(image,str):
        image = cv2.imread(image)

    image = np.array(image)

    # functionality to pass in cv2 object
    if not isinstance(image, str):
        image = image
        img_boxing = image.copy()
    
    # string location of image
    if not isinstance(image,np.ndarray):
        image = np.array(image)
        image = cv2.imread(image)
        img_boxing = image.copy()

    # thresholding 
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    threshold_value = dynamic_threshold_value(gray)
    _, threshold_for_contour = cv2.threshold(gray, threshold_value, 255, cv2.THRESH_TOZERO)
    contours, _ = cv2.findContours(threshold_for_contour, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    
    activations = []
    activation_locations = []
    count = 0

    # evaluating each activation
    
    for contour in contours:
        x, y, w, h = cv2.boundingRect(contour)
        
        if (w > 5 and h > 5) and (w < 60 and h < 60) and len(contour) > 4:
            count += 1
            
            try:
                ellipse = cv2.fitEllipse(contour)
                cv2.ellipse(img_boxing, ellipse, (0, 255, 0), 2)
            
            except cv2.error as e:
                print(f"Error fitting ellipse: {e}")
                continue            

            M = cv2.moments(contour)
            
            # center of activation = activation value 
            if M['m00'] != 0:  # avoid division by zero
                centroid_x = int(M['m10'] / M['m00'])
                centroid_y = int(M['m01'] / M['m00'])
                activation_value = gray[centroid_y, centroid_x]
                activations.append(activation_value)

                # store coordinates 
                loc = [centroid_y, centroid_x]
                activation_locations.append(loc)

            activations.append(activation_value)

            # store coordinates 
            loc = [centroid_y,centroid_x]
            activation_locations.append(loc)
    
    activation_locations = np.array(activation_locations)

    return activation_locations 


def dynamic_threshold_value(image):
    min,max = np.min(image),np.max(image)
    # this regression earns 0.91 classification error --> starting place
    threshold = int((0.94 * max) - 30.13)
    return threshold


def generate_alignment_matrix(data, filename):
    """
        Args:
            data: Data is organized where each row represents filename-->alignment
                  so that the alignment can be fetched for a given matrix given the result

            filename: filename in dataframe that you want alignment from
        
        Returns:
            alignment: alignment matrix representative of anchoring coordinates in CCF space 
                       in format [ox,oy......vy,vz]
    """

    row = data[data['Filenames'] == filename]

    alignment = row[
                       ['ox', 
                        'oy', 
                        'oz', 
                        'ux', 
                        'uy', 
                        'uz', 
                        'vx', 
                        'vy', 
                        'vz']
                   ].values.flatten().tolist()
    
    return alignment

    
def generate_target_slice(ouv, atlas):
    """
        Args:
            alignment: CCF alignment matrix for anchoring to atlas
            volume: a 3d structure to be partitioned according to anchored alignment

        Returns:
            regions: slice of given 3d volume according to anchored alignment 

        Thank you to Harry Carey for this function
    """
    width = None
    height = None
    ox, oy, oz, ux, uy, uz, vx, vy, vz = ouv
    width = np.floor(math.hypot(ux,uy,uz)).astype(int) + 1
    height = np.floor(math.hypot(vx,vy,vz)).astype(int) + 1
    data = np.zeros((width, height), dtype=np.uint32).flatten()
    xdim, ydim, zdim = atlas.shape
    y_values = np.arange(height)
    x_values = np.arange(width)
    hx = ox + vx * (y_values / height)
    hy = oy + vy * (y_values / height)
    hz = oz + vz * (y_values / height)
    wx = ux * (x_values / width)
    wy = uy * (x_values / width)
    wz = uz * (x_values / width)
    lx = np.floor(hx[:, None] + wx).astype(int) 
    ly = np.floor(hy[:, None] + wy).astype(int) 
    lz = np.floor(hz[:, None] + wz).astype(int) 
    valid_indices = (0 <= lx) & (lx < xdim) & (0 <= ly) & (ly < ydim) & (0 <= lz) & (lz < zdim)
    valid_indices = valid_indices.flatten()
    lxf = lx.flatten()
    lyf = ly.flatten()
    lzf = lz.flatten() 
    valid_lx = lxf[valid_indices]
    valid_ly = lyf[valid_indices]
    valid_lz = lzf[valid_indices]
    atlas_slice = atlas[valid_lx,valid_ly,valid_lz]
    data[valid_indices] = atlas_slice
    data_im = data.reshape((height, width))
    return data_im


def volume_to_registered_slice(alignment, volume):
    
    """
        Args:
            alignment: CCF alignment matrix for anchoring to atlas
            volume: a 3d structure to be partitioned according to anchored alignment

        Returns:
            regions: slice of given 3d volume according to anchored alignment 
    """
    
    # anchoring
    Ox, Oy, Oz, Ux, Uy, Uz, Vx, Vy, Vz = alignment

    # just for mouse for now should switch to fetch from volume shape
    bounds = [455, 527, 319]
    X_size = np.sqrt(np.sum(np.square((Ux, Uy, Uz))))
    Z_size = np.sqrt(np.sum(np.square((Vx, Vy, Vz))))

    X_size = np.round(X_size).astype(int)
    Z_size = np.round(Z_size).astype(int)

    Uarange = np.arange(0, 1, 1 / X_size)
    Varange = np.arange(0, 1, 1 / Z_size)

    Ugrid, Vgrid = np.meshgrid(Uarange, Varange)

    Ugrid_x = Ugrid * Ux
    Ugrid_y = Ugrid * Uy
    Ugrid_z = Ugrid * Uz
    Vgrid_x = Vgrid * Vx
    Vgrid_y = Vgrid * Vy
    Vgrid_z = Vgrid * Vz

    X_Coords = (Ugrid_x + Vgrid_x).flatten() + Ox
    Y_Coords = (Ugrid_y + Vgrid_y).flatten() + Oy
    Z_Coords = (Ugrid_z + Vgrid_z).flatten() + Oz

    X_Coords = np.round(X_Coords).astype(int)
    Y_Coords = np.round(Y_Coords).astype(int)
    Z_Coords = np.round(Z_Coords).astype(int)

    out_bounds_Coords = (
        (X_Coords > bounds[0])
        | (Y_Coords > bounds[1])
        | (Z_Coords > bounds[2])
        | (X_Coords < 0)
        | (Y_Coords < 0)
        | (Z_Coords < 0)
    )
    X_pad = X_Coords.copy()
    Y_pad = Y_Coords.copy()
    Z_pad = Z_Coords.copy()

    X_pad[out_bounds_Coords] = 0
    Y_pad[out_bounds_Coords] = 0
    Z_pad[out_bounds_Coords] = 0

    regions = volume[X_pad, Y_pad, Z_pad]

    # fixing rounding errors
    C = len(regions)

    compare = C - X_size * Z_size

    # x is off by one
    if abs(compare) == X_size:
        if compare > 0:
            Z_size += 1
        if compare < 0:
            Z_size -= 1

    # z is off by one
    elif abs(compare) == Z_size:
        if compare > 0:
            X_size += 1
        if compare < 0:
            X_size -= 1

    # both are off by one in the same direction
    elif abs(compare) == X_size + Z_size + 1:
        if compare > 0:
            Z_size += 1
            X_size += 1
        if compare < 0:
            Z_size -= 1
            X_size -= 1

    regions = regions.reshape((abs(Z_size), abs(X_size)))

    return regions

def mask_aspect_ratio(base_image, slice_to_change):
    """
        Args:
            base_image:       image of aspect ratio that the mask needs to match
            slice_to_change:  slice of 3d alignment 

        Returns:
            slice_padded:     slice with adjusted padding and matched aspect ratio so that
                              when slice is resized we do not lose shape of masked ROI 
    """
    
    # aspect ratio of the base image
    target_height, target_width = base_image.shape[:2]
    target_aspect_ratio = target_width / target_height

    # dimensions of slice   
    current_height, current_width = slice_to_change.shape[:2]
    current_aspect_ratio = current_width / current_height
    
    # slice is wider than target aspect ratio --> adjust height
    if current_aspect_ratio > target_aspect_ratio:
        
        new_height = int(current_width / target_aspect_ratio)
        padding_top = (new_height - current_height) // 2
        padding_bottom = new_height - current_height - padding_top
        padding_left = 0
        padding_right = 0
    
    # slice is taller than target aspect ratio --> adjust width
    else:
        
        new_width = int(current_height * target_aspect_ratio)
        padding_left = (new_width - current_width) // 2
        padding_right = new_width - current_width - padding_left
        padding_top = 0
        padding_bottom = 0
    
    # applying specified padding
    slice_padding = cv2.copyMakeBorder(
        slice_to_change,
        padding_top, 
        padding_bottom,
        padding_left, 
        padding_right,
        cv2.BORDER_CONSTANT, value=0
    )
    
    
    # convert to multi channel image for layering in visualization/debugging if single channel
    if len(slice_to_change.shape) == 2:
        slice_dim_add = np.stack((slice_padding,) * 3, axis=-1)

    # already multi-channel image
    else: 
        slice_dim_add = slice_padding

    slice_padded = slice_dim_add
    return slice_padded

def process_mask(structure_id,alignment,image):
    """
        Args:
            structure_id: int representing position of desired mask
            alignment: CCF alignment matrix for anchoring to atlas
            image: a cv2 numpy.ndarray object to be crossed with image

        Returns:
            regions: slice of given 3d volume according to anchored alignment 
    """
    if not isinstance(image,np.ndarray):
        image = np.asarray(image)

    try:
        # generate and refine 3d volume 
        structure_3d, _ = config.mcc.get_structure_mask(structure_id)

        if np.sum(structure_3d) == 0:
            #print(f'Mask for {structure_id} is empty...if this happens please remove from')
            z = 4

        transposed = np.transpose(structure_3d, (2,0,1))
        mirrored = transposed[:, ::-1, ::-1]
        volume = np.array(mirrored)
        slice = generate_target_slice(alignment, volume).astype(np.uint8)
        slice_resize = cv2.resize(slice, (image.shape[1], image.shape[0]), interpolation=cv2.INTER_NEAREST)
            
        cords = np.where(slice_resize)

        """ 
                [[y0,y1,y2...yn-1,yn],[x0,x1,x2...xn-1,xn]] 
                        ↓ ↓ transformation ↓ ↓ 
                [[y0,x0],[y1,x1]...[yn-1,xn-1],[yn,xn]] 
        """
        restructured_coords = np.column_stack((cords[0], cords[1]))

    except HTTPError as http_err:
        print(f"HTTP error for structure {structure_id}: {http_err}")
        restructured_coords = np.empty((0, 2))
        slice_resize = np.zeros_like(image, dtype=np.uint8)
        slice = np.zeros_like(image, dtype=np.uint8)
    
    except ValueError as val_err:
        print(f"ValueError: {val_err}")
        restructured_coords = np.empty((0, 2))
        slice_resize = np.zeros_like(image, dtype=np.uint8)
        slice = np.zeros_like(image, dtype=np.uint8)
    
    except Exception as err:
        print(f"error for structure {structure_id}: {err}")
        restructured_coords = np.empty((0, 2))
        slice_resize = np.zeros_like(image, dtype=np.uint8)
        slice = np.zeros_like(image, dtype=np.uint8)

    return slice_resize, restructured_coords, slice


def cells_in_ROI(mask_coords, cell_coords, threshold=9):
    """
    Find cells within the region of interest (ROI) defined by mask coordinates.

    Args:
        mask_coords: Coordinates of the mask.
        cell_coords: Coordinates of the cells.
        threshold: Distance threshold for considering cells within the ROI.

    Returns:
        np.ndarray: Coordinates of cells within the ROI.
        int: Number of cells within the ROI.
    """
    
    mask_tree = cKDTree(mask_coords)
    cells_in_mask = []
    for cell in cell_coords:
        distance, _ = mask_tree.query(cell)
        if distance <= threshold:
            cells_in_mask.append(cell)
    
    cells_in_mask = np.array(cells_in_mask)
    
    return np.array(cells_in_mask), len(cells_in_mask)


def matrix_old(data, filename):
    
    """
        Args:
            data: Data is organized where each row represents filename-->alignment
                  so that the alignment can be fetched for a given matrix given the result

            filename: filename in dataframe that you want alignment from

        Returns:
            alignment: alignment matrix representative of anchoring coordinates in CCF space 
    """

    ox = data['ox'].values
    oy = data['oy'].values
    oz = data['oz'].values
    ux = data['ux'].values
    uy = data['uy'].values
    uz = data['uz'].values
    vx = data['vx'].values
    vy = data['vy'].values
    vz = data['vz'].values

    row = data.loc[filename]

    ox = row['ox']
    oy = row['oy']
    oz = row['oz']
    ux = row['ux']
    uy = row['uy']
    uz = row['uz']
    vx = row['vx']
    vy = row['vy']
    vz = row['vz']

    alignment = [ox,oy,oz,ux,uy,uz,vx,vy,vz]

    return alignment


def save_figure(image, filename):
    plt.figure(figsize=(10, 10))
    plt.imshow(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))
    plt.axis('off')
    plt.savefig(filename, dpi=600, bbox_inches='tight', pad_inches=0)
    plt.close()


def save_image(data_dir, structures, filename, overlay_activations: bool, overlay_mask: bool, save_dir): 
    image = cv2.imread(os.path.join(data_dir, filename))
    
    clean_image_path = os.path.join(save_dir, f"{os.path.splitext(filename)[0]}_clean.png")
    plt.figure(figsize=(10, 10))
    plt.imshow(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))  
    plt.axis('off')
    plt.title(filename)
    plt.savefig(clean_image_path, dpi=600, bbox_inches='tight', pad_inches=0)
    plt.close()
    
    for acronymn in structures:
        download_structure(data_dir, acronymn, filename, overlay_activations, overlay_mask, save_dir, image)


def download_structure(data_dir, 
                   structure, 
                   filename, 
                   overlay_activations: bool, 
                   overlay_mask: bool, 
                   save_dir, 
                   image):
    structure_id = acronymn_to_id(structure)
    data = pd.read_csv(os.path.join(data_dir, 'results.csv')) 
    alignment = generate_alignment_matrix(data, filename)
    mask, mask_coords = process_mask(structure_id, alignment, image)
    cell_yx = process_image(image)
    local_activations, count = cells_in_ROI(mask_coords, cell_yx, threshold=9)    

    # enforce local_activations as a 2d array
    if local_activations.ndim == 1 and len(local_activations) == 0:
        print(f'no cells detected in roi <{structure}>')
        local_activations = np.empty((0, 2))
    
    # enforce mask_coords as a 2d array
    if mask_coords.ndim == 1 and len(mask_coords) == 0:
        print('No mask detected on this coronal axis')
        mask_coords = np.empty((0, 2))
    
    # image with mask overlay
    if overlay_mask and mask_coords.size > 0:
        mask_image_path = os.path.join(save_dir, f"{os.path.splitext(filename)[0]}_{structure}_mask.png")
        plt.figure(figsize=(10, 10))
        plt.imshow(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))
        plt.scatter(mask_coords[:, 1], mask_coords[:, 0], s=1, c='red')
        plt.title(f'Predicted mask of {structure} in {filename}')
        plt.axis('off')
        plt.savefig(mask_image_path, dpi=600, bbox_inches='tight', pad_inches=0)
        plt.close()

    # image with activations overlay
    if overlay_activations and local_activations.size > 0:
        activations_image_path = os.path.join(save_dir, f"{os.path.splitext(filename)[0]}_{structure}_mask_activations.png")
        plt.figure(figsize=(10, 10))
        plt.imshow(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))
        plt.scatter(local_activations[:, 1], local_activations[:, 0], s=3, edgecolor='green', facecolor='none', marker='s')
        plt.title(f'Cells in {structure} in {filename}')
        plt.axis('off')
        plt.savefig(activations_image_path, dpi=600, bbox_inches='tight', pad_inches=0)
        plt.close()


def download_directory(data_dir, structures, overlay_activations: bool, overlay_mask: bool, save_dir):
    for filename in os.listdir(data_dir):
        if filename.endswith(('.png', '.jpg', '.jpeg', '.tif', '.bmp')):
            print(f"visualizing results for {filename}")
            save_image(data_dir, structures, filename, overlay_activations, overlay_mask, save_dir)


def download_overlay(data_dir,filename,mask_yx,local_yx,ac,image):
    """
        Given mask predictions and corresponding local cellular activity, display a mask

        Args:
            save_path: Coordinates of the mask.
            mask_yx: Coordinates of the cells.
            local_yx: Distance threshold for considering cells within the ROI.

        Returns:
            N/A 
    """
    file_base,_ = os.path.splitext(filename)

    # enforce local_activations as a 2d array
    if local_yx.ndim == 1 and len(local_yx) == 0:
        local_yx = np.empty((0, 2))        # no cells detected in ROI

    # enforce mask_coords as a 2d array
    if mask_yx.ndim == 1 and len(mask_yx) == 0:
        mask_yx = np.empty((0, 2))         # no mask detected on this coronal axis 


    # save mask overlay
    if mask_yx.size > 0:
        save_path = os.path.join(data_dir,f"{file_base}_mask_{ac}.png")        
        plt.figure(figsize=(10, 10))
        plt.imshow(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))
        plt.scatter(mask_yx[:, 1], mask_yx[:, 0], s=1, c='red')
        plt.title(f'Predicted mask of {ac} in {filename}')
        plt.axis('off')
        plt.savefig(save_path, dpi=600, bbox_inches='tight', pad_inches=0)
        plt.close()

    # save local activations overlay
    if local_yx.size > 0:
        save_path = os.path.join(data_dir, f"{file_base}_cells_{ac}.png")
        plt.figure(figsize=(10, 10))
        plt.imshow(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))
        plt.scatter(local_yx[:, 1], local_yx[:, 0], s=3, edgecolor='green', facecolor='none', marker='s')
        plt.title(f'Cells in {ac} in {filename}')
        plt.axis('off')
        plt.savefig(save_path, dpi=600, bbox_inches='tight', pad_inches=0)
        plt.close()


def gen_mask(structure_id):
    # alternate routing for masking
    mask_writer = functools.partial(ReferenceSpace.check_and_write, '/Users/riley/Desktop/NeuGenes/Models/model/mcc/annotation/ccf_2017/structure_masks/resolution_25')
    mask_generator = config.reference_space.many_structure_masks([structure_id], mask_writer)
        
    for m in mask_generator:
        print('downloaded')

def display_structure_3d(id):
    """ 
        Given an Allen structure id, return the 
    
    """
    structure_id = id
    gen_mask(id)
    path = f'/Users/riley/Desktop/NeuGenes/Models/model/mcc/annotation/ccf_2017/structure_masks/resolution_25/structure_{structure_id}.nrrd'
    data,_ = nrrd.read(path)
    
    def display_slice(slice_index):
        plt.imshow(data[slice_index], cmap='gray')
        plt.show()

    widgets.interact(display_slice, slice_index=(0, data.shape[0] - 1))


# function to query last image of tiff stack, save in new directory
def registration_directory(input_dir,channels):
    
    output_dir = f'{input_dir}_register'
    os.makedirs(output_dir, exist_ok=True)
    for fn in os.listdir(input_dir):
        if fn.lower().endswith('.tiff') or fn.lower().endswith('.tif'):
            
            fn_path = os.path.join(input_dir, fn)
            tiff_stack = Image.open(fn_path) # gathering last image in tiff stack, saving
            path = os.path.join(output_dir,fn.rsplit('.', 1)[0] + '.png')
            
            if(os.path.exists(path)):
                continue
            
            image = tiff_stack.seek(len(channels)-1)
            resized_image = match_aspect_ratio(image, (480, 358))
            resized_image.save(path, format='PNG')
    
    os.chdir(input_dir)
    return output_dir


def process_tiff_stack(filename,channels,alignment,structures):
    
    image_result = {}
    for c, channel in enumerate(channels):
        channel_result = {}
        
        # opening tiff stack
        tiff_stack = Image.open(filename)
        tiff_stack.verify()
        tiff_stack = Image.open(filename)
        tiff_stack.seek(c)

        # resize to match registration
        resized_image = match_aspect_ratio(tiff_stack, (480, 358))
        cell_yx = process_image(resized_image)

        for id in structures:
            ac = id_to_acronym(id)
            _, mask_yx, _ = process_mask(id, alignment, resized_image)
            _, count = cells_in_ROI(mask_yx, cell_yx, threshold=9)
            channel_result[ac] = count
    
        image_result[channel] = channel_result

    return {filename.rsplit('.', 1)[0] + '.tif': image_result}
        


   
