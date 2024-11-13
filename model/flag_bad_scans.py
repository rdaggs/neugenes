import os
import cv2 
import argparse
import model.config as config
from model.utils import dynamic_threshold_value

os.chdir(config.root_directory)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir',type=str,required=True,help='Input directory containing the png brain scans')
    args = parser.parse_args()
    flag_scans(args.input_dir)

def flag_scans(input_dir):
    data_path = os.path.join(os.path.dirname(config.root_directory),input_dir)

    for fn in os.listdir(data_path):
        path = os.path.join(data_path,fn)

        if os.path.isfile(path) and path.lower().endswith(('.jpeg', '.tif', '.jpg', '.png')):
            
            # multiple pieces of tissue
            image = cv2.imread(path)
            threshold = dynamic_threshold_value(image)
            _, thresh = cv2.threshold(image, threshold, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
            contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            num_bodies = len(contours)
            if num_bodies > 2:
                print(f"revisit scan titled {fn} and crop extraneous tissue to prevent fault in registration ")

            if threshold > 220:
                print(f"revisit scan titled {fn} and adjust exposure so that expression intensity can be recorded more accurately")
                


if __name__ == "__main__":
    main()
