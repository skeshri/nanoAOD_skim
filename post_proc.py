#!/usr/bin/env python
import os,sys

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from HZZAnalysisModule import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
from JetSFMaker import *

testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv2/GluGluHToZZTo2L2Q_M3000_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/50000/11F21FB5-7535-B249-8F59-1173DC384F33.root"
# testfile = "file:/tmp/rasharma/11F21FB5-7535-B249-8F59-1173DC384F33.root"

entriesToRun = 100  # 0 if need to run over all entries else put number of entries to run
isMCTrueFalse = True
# Keep DownloadFileToLocalThenRun=True this should reduce the file read error from eos.
DownloadFileToLocalThenRun=True

if testfile.find("SingleMuon") != -1 or testfile.find("EGamma") != -1 or testfile.find("SingleElectron") != -1 or testfile.find("DoubleMuon") != -1 or testfile.find("MuonEG") != -1 or testfile.find("DoubleEG") != -1:
  isMCTrueFalse=False
  Year=2018
  if testfile.find("Run2016") != -1:
    Year=2016
    jsonFileName="Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
  if testfile.find("Run2017") != -1:
    Year=2017
    jsonFileName="Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
  if testfile.find("Run2018") != -1:
    Year=2018
    jsonFileName="Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"
  jsonFileName = "data/jsonFiles/"+jsonFileName
  print "\n===> Running over ",Year," data...\n"
  print "===> JSON File: ",jsonFileName
  # jetmetCorrector = createJMECorrector(isMC=isMCTrueFalse, dataYear=Year, jesUncert="Merged", jetType = "AK4PFchs")
  # fatJetCorrector = createJMECorrector(isMC=isMCTrueFalse, dataYear=Year, jesUncert="Merged", jetType = "AK8PFPuppi")
  # p=PostProcessor(".",[testfile],None,None,[HZZAnalysisModule(),fatJetCorrector()],provenance=False,fwkJobReport=False,jsonInput=jsonFileName,maxEntries=entriesToRun,prefetch=DownloadFileToLocalThenRun)
  p=PostProcessor(".",[testfile],None,None,[HZZAnalysisModule()],provenance=False,fwkJobReport=False,jsonInput=jsonFileName,maxEntries=entriesToRun,prefetch=DownloadFileToLocalThenRun)
else:
  print "==> Processing a MC file..."
  isMCTrueFalse=True
  year = 2018
  if testfile.find("VVjj_2018v") != -1: year = 2018
  if testfile.find("RunIIAutumn18NanoAODv") != -1: year = 2018
  if testfile.find("RunIISummer20UL18") != -1: year = 2018

  if testfile.find("RunIISummer19UL17NanoAOD") != -1: year = 2017
  if testfile.find("RunIIFall17NanoAODv") != -1: year = 2017
  if testfile.find("RunIISummer20UL17") != -1: year = 2017
  if testfile.find("VVjj_2017v") != -1: year = 2017

  if testfile.find("RunIISummer16NanoAODv") != -1: year = 2016
  if testfile.find("RunIISummer20UL16") != -1: year = 2016
  if testfile.find("VVjj_2016v") != -1: year = 2016

  # jetmetCorrector = createJMECorrector(isMC=isMCTrueFalse, dataYear=year, jesUncert="Merged", jetType = "AK4PFchs")
  # fatJetCorrector = createJMECorrector(isMC=isMCTrueFalse, dataYear=year, jesUncert="Total", jetType = "AK8PFPuppi")
  # fatJetCorrector = createJMECorrector(isMC=isMCTrueFalse, dataYear=year, jesUncert="Merged", jetType = "AK4PFchs")

  if year == 2016:
    era="Legacy2016"
    sfFileName="DeepCSV_2016LegacySF_V1.csv"
  if year == 2017 or year == 'UL2017':
    era="2017"
    sfFileName="DeepCSV_94XSF_V5_B_F.csv"
  if year == 2018:
    era="2018"
    sfFileName="DeepCSV_102XSF_V2.csv"
  btagSF = lambda: btagSFProducer(era,algo="deepcsv",selectedWPs=['L','M','T','shape_corr'],sfFileName=sfFileName)
  puidSF = lambda: JetSFMaker("%s" % year)
  # p=PostProcessor(".",[testfile],"",None,[HZZAnalysisModule(),fatJetCorrector(),btagSF(),puidSF()],provenance=True,fwkJobReport=False,maxEntries=entriesToRun,prefetch=DownloadFileToLocalThenRun)
  p=PostProcessor(".",[testfile],"",None,[HZZAnalysisModule(),btagSF(),puidSF()],provenance=True,fwkJobReport=False,maxEntries=entriesToRun,prefetch=DownloadFileToLocalThenRun)
  # p=PostProcessor(".",[testfile],"",None,[HZZAnalysisModule()],provenance=True,fwkJobReport=False,maxEntries=entriesToRun,prefetch=DownloadFileToLocalThenRun)

p.run()
print "DONE"
# os.system("ls -lR")
