#!/bin/bash
set -e

python3.11 -m venv venv
source venv/bin/activate
pushd model
pip install -r requirements.txt
python -u main.py --dir images --structures 'BLA,PVT,MH,CA,'
popd
t
