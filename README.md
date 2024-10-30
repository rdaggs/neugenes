# NeuGenes
Colocalization analysis, Brain Registration to Allen CCFv3, Protein signaling for axonal tracing

## Installation

Clone the repo
```bash
git clone git@github.com:rdaggs/neugenes.git
cd neugenes
```

Then run (we recommend using a python virtual environment)

```bash
pip install --upgrade pip
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
```bash
unzip resolution_25.zip -d model.config.mask_zip_file
# alternatively you can wait and let processing interface directly with the API for downloading the structures
```

## Data Preprocessing

Some of the image cleaning and preprocessing is handled internally. 
The output of standard nikon microscopy equipment (xx.nd2) is not supported yet so you'll need
to process and save as tiff, jpeg, png etc.
```bash
python model.flag_bad_scans.py
```
Watch for broken tissue and pay attention to expression brightness. Some of the contour in the cell expression can be lost if you don't adjust contrast/brightness. 




## Inference

Run
```bash
python -u main.py [--dir DATASET_DIRECTORY] [--structures BLA,PVT,SCN.....]
```

for instance:

```bash
python -u main.py fit --data Npsr1-cre_Ai14 --structures BLA,PVT,HPF,LS
python -u main.py fit --data DREADD_m2.7 --structures FULL_BRAIN

```

