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

## Inference

Run
```bash
python -u main.py [--dir DATASET_DIRECTORY] [--structures BLA,PVT,SCN.....]
```

for instance:

```bash
python -u main.py fit --data Npsr1-cre_Ai14 --structures BLA,PVT,SCN,HPF,NTm,ENTmv,PAR,POST,PRE,SUB.....
```

