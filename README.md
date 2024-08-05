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

Run the following command to install the domain-exploration library:

```bash
pip install -e .
```

## Training

Run
```bash
python main.py fit [--data DATA_MODULE] [--data.domain DOMAIN] [--model MODEL]
```

for instance:

```bash
python main.py fit --data MNISTDataModule --data.domain pixel --model MLP_MNIST
```
