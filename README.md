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
pip install --upgrade setuptools
```

Run the following command to install the neugenes library:

```bash
pip install -e .
```

Download preliminary structures
```bash
unzip resolution_25.zip -d model.config.mask_zip_file
```

## Data Preprocessing

Some of the image cleaning and preprocessing is handled internally. 
The output of standard nikon microscopy equipment (xx.nd2) is not supported yet so you'll need
to process and save as tiff, jpeg, png etc.
```bash
python -u flag_bad_scans.py
```
Watch for bad scans 
Some of the contour in the cell expression can be lost 
if you don't adjust contrast/brightness. 




## Inference

Run
```bash
python -u main.py [--dir DATASET_DIRECTORY] [--structures BLA,PVT,SCN.....]
```

for instance:

```bash
python -u main.py fit --data Npsr1-cre_Ai14 --structures BLA,PVT,SCN,HPF,NTm,ENTmv,PAR,POST,PRE,SUB.....
```

