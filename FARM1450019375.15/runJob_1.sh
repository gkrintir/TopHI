#!/bin/bash
export X509_USER_PROXY=/afs/cern.ch/user/g/gkrintir/github/HI/CMSSW_7_5_7_patch2/src/UserCode/diall/proxyforprod
export OUTDIR=/afs/cern.ch/user/g/gkrintir/github/HI/CMSSW_7_5_7_patch2/src/UserCode/diall/out
cd /afs/cern.ch/user/g/gkrintir/github/HI/CMSSW_7_5_7_patch2/src
eval `scram r -sh`
cd -
cp /afs/cern.ch/user/g/gkrintir/github/HI/CMSSW_7_5_7_patch2/src/UserCode/diall/FARM1450019375.15/ExampleAnalysisParameters_cfg.py .
runExample ExampleAnalysisParameters_cfg.py 0 1
export OUTPUT=AnaResults_1.root
cp AnaResults.root $OUTDIR/$OUTPUT
rm AnaResults.root
