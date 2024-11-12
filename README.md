# NeuGenes
Colocalization analysis, Brain Registration to Allen CCFv3, Protein signaling for axonal tracing

## Installation

Clone the repo
```bash
git clone https://github.com/rdaggs/neugenes.git
cd neugenes
```

Then run (we recommend using a python virtual environment)

```bash
pip install --upgrade pip
set -e
python3.11 -m venv venv
source venv/bin/activate
```

Update the setuptools
```bash
# if errors arise with pip, try this
pip install --upgrade setuptools
```

Run the following command to integrate the neccessary build distribution of tensorflow 2.15 (dependent machine architecture):

```bash
pip install tensorflow-2.15.0-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
pip install tensorflow-2.15.0-cp311-cp311-manylinux_2_17_aarch64.manylinux2014_aarch64.whl 
pip install tensorflow-2.15.0-cp310-cp310-macosx_10_15_x86_64.whl
pip install tensorflow-2.15.0-cp310-cp310-macosx_12_0_arm64.whl 
```

Download preliminary structures
(or wait and let processor interface directly with API for downloading the structures)
```bash
unzip resolution_25.zip -d model.config.mask_zip_file

```

## Data Preprocessing

Some of the image cleaning and preprocessing is handled internally. 
The output of standard nikon microscopy equipment (xx.nd2) is not supported yet so you'll need
to process and save as tiff, jpeg, png etc.
```bash
python model.flag_bad_scans.py
```
Given you have multiple datasets, condition your dataset to only compare the axially corresponding scans.
i.e. if coronal axis of brain scans in control = [456,228,....111] and stress = [461,219,....41], this
alternate preprocessing step can eliminate the final scan for a more precise comparison (given a threshold score
of your choice)
```bash
python model.condition_coronal_correlation.py
```




## Inference

Run
```bash
python -u main.py [--dir DATASET_DIRECTORY] [--structures BLA,PVT,SCN.....]
```

for instance:

```bash
python -u main.py -u --input_dir Npsr1-cre_Ai14      --structures BLA,PVT,HPF,LS
python -u main.py -u --input_dir DREADD_m2.7_CONTROL --structures FULL_BRAIN

```

