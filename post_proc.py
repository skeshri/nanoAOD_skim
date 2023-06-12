#!/usr/bin/env python
import os,sys

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from HZZAnalysisModule import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
from JetSFMaker import *

testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv2/GluGluHToZZTo2L2Q_M3000_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/50000/11F21FB5-7535-B249-8F59-1173DC384F33.root"


entriesToRun = 1000  # 0 if need to run over all entries else put number of entries to run
isMCTrueFalse = True
# Keep DownloadFileToLocalThenRun=True this should reduce the file read error from eos.
DownloadFileToLocalThenRun=True
p=PostProcessor(".",[testfile],"",None,[HZZAnalysisModule()],provenance=True,fwkJobReport=False,haddFileName="nano.root",maxEntries=entriesToRun,prefetch=DownloadFileToLocalThenRun)


p.run()
print "DONE"