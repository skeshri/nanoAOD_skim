#!/usr/bin/env python
import os,sys

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from HZZAnalysisModule import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
from JetSFMaker import *

testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer19UL17NanoAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_mc2017_realistic_v6-v2/270000/07E11789-6AC3-9943-A770-53E366C5A7A7.root"


entriesToRun = 100  # 0 if need to run over all entries else put number of entries to run
isMCTrueFalse = True
Year = 2016
# Keep DownloadFileToLocalThenRun=True this should reduce the file read error from eos.
DownloadFileToLocalThenRun=False

if testfile.find("SingleMuon") != -1 or testfile.find("EGamma") != -1 or testfile.find("SingleElectron") != -1 or testfile.find("DoubleMuon") != -1 or testfile.find("MuonEG") != -1 or testfile.find("DoubleEG") != -1:
  isMCTrueFalse=False
  if testfile.find("Run2016") != -1:
    Year=2016
    jsonFileName="Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt"
  if testfile.find("Run2017") != -1:
    Year=2017
    jsonFileName="Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt"
  if testfile.find("Run2018") != -1:
    Year=2018
    jsonFileName="Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt"
  jsonFileName = "data/jsonFiles/"+jsonFileName
  print "\n===> Running over ",Year," data...\n"
  print "===> JSON File: ",jsonFileName
  # jetmetCorrector = createJMECorrector(isMC=isMCTrueFalse, dataYear=Year, jesUncert="Merged", jetType = "AK4PFchs")
  fatJetCorrector = createJMECorrector(isMC=isMCTrueFalse, dataYear=Year, jesUncert="Merged", jetType = "AK8PFPuppi")
  p=PostProcessor(".",[testfile],None,None,[HZZAnalysisModule(),fatJetCorrector()],provenance=False,fwkJobReport=False,jsonInput=jsonFileName,maxEntries=entriesToRun,haddFileName="nano.root",prefetch=DownloadFileToLocalThenRun)
else:
  print "==> Processing a MC file..."
  isMCTrueFalse=True
  if testfile.find("RunIIAutumn18NanoAODv") != -1: year = 2018
  if testfile.find("VVjj_2018v") != -1: year = 2018
  if testfile.find("RunIIFall17NanoAODv") != -1: year = 2017
  if testfile.find("RunIISummer19UL17NanoAOD") != -1: year = 2017
  if testfile.find("VVjj_2017v") != -1: year = 2017
  if testfile.find("RunIISummer16NanoAODv") != -1: year = 2016
  if testfile.find("VVjj_2016v") != -1: year = 2016
  # jetmetCorrector = createJMECorrector(isMC=isMCTrueFalse, dataYear=year, jesUncert="Merged", jetType = "AK4PFchs")
  fatJetCorrector = createJMECorrector(isMC=isMCTrueFalse, dataYear=year, jesUncert="Merged", jetType = "AK8PFPuppi")
  if year == 2016:
    era="Legacy2016"
    sfFileName="DeepCSV_2016LegacySF_V1.csv"
  if year == 2017:
    era="2017"
    sfFileName="DeepCSV_94XSF_V5_B_F.csv"
  if year == 2018:
    era="2018"
    sfFileName="DeepCSV_102XSF_V2.csv"
  btagSF = lambda: btagSFProducer(era,algo="deepcsv",selectedWPs=['L','M','T','shape_corr'],sfFileName=sfFileName)
  puidSF = lambda: JetSFMaker("%s" % year)
  p=PostProcessor(".",[testfile],"",None,[HZZAnalysisModule(),fatJetCorrector(),btagSF(),puidSF()],provenance=True,fwkJobReport=False,maxEntries=entriesToRun,haddFileName="nano.root",prefetch=DownloadFileToLocalThenRun)

p.run()
print "DONE"
#os.system("ls -lR")
