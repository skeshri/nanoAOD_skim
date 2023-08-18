#!/usr/bin/env python
import os,sys

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from H4Lmodule import *
from H4LCppModule import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *

# 4l ggH Sample
# testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/120000/3ED05633-EBB7-4A44-8F9D-CD956490BCFD.root"

# 2l2q ggH Sample 1000 GeV
## UL2018
testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/GluGluHToZZTo2L2Q_M1000_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/70000/F920B631-39E7-EE4C-9D8D-B242A2103115.root"

## UL2017
# testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluHToZZTo2L2Q_M1000_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/40000/4E6C744B-023A-E14A-9ADD-663C7101E98E.root"

# 2l2nu ggH Sample 1000 GeV
# testfile = "root://cms-xrd-global.cern.ch/"

entriesToRun = 0 # 0 if need to run over all entries else put number of entries to run
isMCTrueFalse = True
# Keep DownloadFileToLocalThenRun=True this should reduce the file read error from eos.
DownloadFileToLocalThenRun=True
if testfile.find("UL18") != -1:
    year = 2018
if testfile.find("UL17") != -1:
    year = 2017
H4LCppModule = lambda: HZZAnalysisCppProducer(year)
#p=PostProcessor(".",[testfile],"",None,[H4LCppModule()],provenance=True,fwkJobReport=False,haddFileName="nano_M125.root",maxEntries=entriesToRun,prefetch=DownloadFileToLocalThenRun,outputbranchsel="keep_and_drop.txt")
p=PostProcessor(".",[testfile],"",None,[H4LCppModule()],provenance=True,fwkJobReport=False,haddFileName="nano_M125_cpp.root",maxEntries=entriesToRun,prefetch=DownloadFileToLocalThenRun,outputbranchsel="keep_and_drop.txt")

p.run()
print "DONE"
