#!/usr/bin/env python
import os,sys

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from wvInclusiveAnalysisModule import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *

testfile = "root://cms-xrd-global.cern.ch//store/data/Run2018A/SingleMuon/NANOAOD/Nano25Oct2019-v1/20000/D03C6AE0-73AD-A940-B8CA-779A621D4853.root"


entriesToRun = 500  # 0 if need to run over all entries else put number of entries to run
isMCTrueFalse = True
Year = 2016
# Keep DownloadFileToLocalThenRun=True this should reduce the file read error from eos.
DownloadFileToLocalThenRun=True

if testfile.find("SingleMuon") != -1 or testfile.find("EGamma") != -1 or testfile.find("SingleElectron") != -1 or testfile.find("DoubleMuon") != -1 or testfile.find("MuonEG") != -1 or testfile.find("DoubleEG") != -1:
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
  fatJetCorrector = createJMECorrector(isMC=isMCTrueFalse, dataYear=Year, jesUncert="All", redojec=True, jetType = "AK8PFPuppi")
  p=PostProcessor(".",[testfile],None,None,[wvInclusiveAnalysisModule(),fatJetCorrector()],provenance=False,fwkJobReport=False,jsonInput=jsonFileName,maxEntries=entriesToRun,haddFileName="nano.root",prefetch=DownloadFileToLocalThenRun)
else:
  print "==> Processing a MC file..."
  isMCTrueFalse=True
  if testfile.find("RunIIAutumn18NanoAODv") != -1: year = 2018
  if testfile.find("RunIIFall17NanoAODv") != -1: year = 2017
  if testfile.find("RunIISummer16NanoAODv") != -1: year = 2016
  # jetmetCorrector = createJMECorrector(isMC=isMCTrueFalse, dataYear=year, jesUncert="All", redojec=True, jetType = "AK4PFchs")
  fatJetCorrector = createJMECorrector(isMC=isMCTrueFalse, dataYear=year, jesUncert="All", redojec=True, jetType = "AK8PFPuppi")
  p=PostProcessor(".",[testfile],"","keep_and_drop_inclusive.txt",[wvInclusiveAnalysisModule(),fatJetCorrector()],provenance=True,fwkJobReport=False,maxEntries=entriesToRun,haddFileName="nano.root",prefetch=DownloadFileToLocalThenRun)

p.run()
print "DONE"
#os.system("ls -lR")
