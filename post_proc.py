#!/usr/bin/env python
import os,sys

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from wvAnalysisModule import *

testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18NanoAODv5/WplusTo2JWminusToLNuJJ_EWK_LO_SM_MJJ100PTJ10_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/Nano1June2019_102X_upgrade2018_realistic_v19-v1/250000/F97E2CDF-3050-DD40-8B35-8CB1BCF48695.root"

if testfile.find("SingleMuon") != -1 or testfile.find("EGamma") != -1:
  print "==> Processing a data file..."
  p=PostProcessor(".",[testfile],None,None,[wvAnalysisModule()],provenance=False,fwkJobReport=False,jsonInput="Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt",)
  #p=PostProcessor(".",[testfile],None,None,[wvAnalysisModule()],provenance=False,fwkJobReport=False,jsonInput="Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt",maxEntries=21000,)

else:
  print "==> Processing a MC file..."
  from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
  jetmetCorrector = createJMECorrector(isMC=True, dataYear=2018, jesUncert="All", redojec=True, jetType = "AK4PFchs")
  fatJetCorrector = createJMECorrector(isMC=True, dataYear=2018, jesUncert="All", redojec=True, jetType = "AK8PFPuppi")
  p=PostProcessor(".",[testfile],"","keep_and_drop.txt",[jetmetCorrector(),fatJetCorrector(),wvAnalysisModule()],provenance=True,)
  #p=PostProcessor(".",[testfile],"","keep_and_drop.txt",[jetmetCorrector(),fatJetCorrector(),wvAnalysisModule()],provenance=True,maxEntries=21000,)

p.run()

print "DONE"
#os.system("ls -lR")
