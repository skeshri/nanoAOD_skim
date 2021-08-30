#!/usr/bin/env python
import os,sys

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from HHAnalysisModule import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *

testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAODv7/GluGluToHHTo2G4Q_node_SM_TuneCUETP8M1_PSWeights_13TeV-madgraph-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/100000/04C6ACF7-FA1B-4147-BC43-49AB7273C745.root"


entriesToRun = 0  # 0 if need to run over all entries else put number of entries to run
isMCTrueFalse = True
Year = 2017
# Keep DownloadFileToLocalThenRun=True this should reduce the file read error from eos.
DownloadFileToLocalThenRun=True

if testfile.find("SinglePhoton") != -1 or testfile.find("EGamma") != -1 or testfile.find("SingleElectron") != -1 or testfile.find("DoubleMuon") != -1 or testfile.find("MuonEG") != -1 or testfile.find("DoubleEG") != -1:
  isMCTrueFalse=False
  if testfile.find("Run2016") != -1:
    Year=2016
    jsonFileName="Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt"
  if testfile.find("Run2017") != -1:
    Year=2017
    jsonFileName="Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt"
  if testfile.find("Run2018") != -1:
    Year=2018
    jsonFileName="Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt"
  print "\n===> Running over ",Year," data...\n"
  print "===> JSON File: ",jsonFileName
  # jetmetCorrector = createJMECorrector(isMC=isMCTrueFalse, dataYear=Year, jesUncert="All", redojec=True, jetType = "AK4PFchs")
  jetmetCorrector = createJMECorrector(isMC=isMCTrueFalse, dataYear=Year, jesUncert="Merged", jetType = "AK4PFchs")
  # fatJetCorrector = createJMECorrector(isMC=isMCTrueFalse, dataYear=Year, jesUncert="All", redojec=True, jetType = "AK8PFPuppi")
  fatJetCorrector = createJMECorrector(isMC=isMCTrueFalse, dataYear=Year, jesUncert="Merged", jetType = "AK8PFPuppi")
  p=PostProcessor(".",[testfile],None,None,[HHAnalysisModule(),fatJetCorrector(),jetmetCorrector()],provenance=False,fwkJobReport=False,jsonInput=jsonFileName,maxEntries=entriesToRun,haddFileName="nano_allEvt.root",prefetch=DownloadFileToLocalThenRun)
else:
  print "==> Processing a MC file..."
  isMCTrueFalse=True
  if testfile.find("RunIIAutumn18NanoAODv") != -1: year = 2018
  if testfile.find("RunIIFall17NanoAODv") != -1: year = 2017
  if testfile.find("RunIISummer16NanoAODv") != -1: year = 2016
  # jetmetCorrector = createJMECorrector(isMC=isMCTrueFalse, dataYear=year, jesUncert="All", redojec=True, jetType = "AK4PFchs")
  # fatJetCorrector = createJMECorrector(isMC=isMCTrueFalse, dataYear=year, jesUncert="All",  jetType = "AK8PFPuppi")
  fatJetCorrector = createJMECorrector(isMC=isMCTrueFalse, dataYear=year, jesUncert="Merged", jetType = "AK8PFPuppi")
  jetmetCorrector = createJMECorrector(isMC=isMCTrueFalse, dataYear=year, jesUncert="Merged", jetType = "AK4PFchs")

  p=PostProcessor(".",[testfile],"","keep_and_drop_inclusive.txt",[HHAnalysisModule(),fatJetCorrector(),jetmetCorrector()],provenance=True,fwkJobReport=False,maxEntries=entriesToRun,haddFileName="nano_allEvt.root",prefetch=DownloadFileToLocalThenRun)

p.run()
print "DONE"
#os.system("ls -lR")
