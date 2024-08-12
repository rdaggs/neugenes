#!/bin/bash
set -e

python3.11 -m venv venv
source venv/bin/activate

pushd model
pip install -r requirements.txt
python -u main.py --dir [YOUR_DIRECTORY_HERE] --structures MOp,MOs,SSp-ll,SSp-tr,SSs,VISC,AUDpo,VISal,VISam,VISl,VISp,VISpl,VISpm,VISli,VISpor,ACAd,ACAv,PL,ILA,ORB,VISrl,TEa,AId,AIp,AIv,RSPagl,RSPd,RSPv,PTLp,PERI,ECT,TT,DP,PIR,NLOT,COA,PAA,TR,HPF,FC,ENTl,ENTm,ENTmv,PAR,POST,PRE,SUB,ProS,CLA,EP,LA,BMA,ACB,FS,LSc,LSr,LSv,SH,sAMY,MEA,PAL,GPe,GPi,PALv,SI,MA,PALm,MA,NDB,BSTa,BSTp,SPF,PP,MG,SGN,AM,AD,MD,PT,RE,RH,PCN,SubG,SO,ASO,NC,PVH,ARH,PVR,ADP,AHA,MEPO,PS,SCH,MEZ,VMH,TU,ZI,VTA,RN,SOC,TRN,LDT,CN,PARN,RPA,DN,PVT,BLA
popd

# processing time ≈ 1 min per image 


