#!/usr/bin/env python
import os,sys

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from H4Lmodule import *
from H4LCppModule import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
from JetSFMaker import *

ifRunningOnCondor = False

testfilelist = []

if ifRunningOnCondor:
    testfile = "root://cms-xrd-global.cern.ch//store/data/Run2018A/SingleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9_GT36-v1/2820000/06EFD7A6-A6F1-E145-B1F8-C6794FE4EDDD.root"
    testfilelist.append(testfile)

else:
    if len(sys.argv) > 1:
       InputFileList = sys.argv[1]
    else:
       InputFileList = "ExampleInputFileList.txt"
    with open(InputFileList, 'r') as file:
      for line in file:
        # Remove newline characters
        line = line.strip()
        # Append the line to the list with the "root://cms-xrd-global.cern.ch//" prefix
        testfilelist.append("root://cms-xrd-global.cern.ch/" + line)

# Set  entriesToRun = 0 if need to run over all entries else put number of entries to run
entriesToRun = 0 if ifRunningOnCondor else 100

isMC = True
isFSR = False
jsonFileName = ""
# Keep DownloadFileToLocalThenRun=True this should reduce the file read error from eos.
DownloadFileToLocalThenRun=True

if testfilelist[0].find("/data/") != -1:
    isMC = False

if testfilelist[0].find("UL18") != -1 or testfilelist[0].find("UL2018") != -1: # UL2018 for identification of 2018 UL data and UL18 for identification of 2018 UL MC
    year = 2018
    cfgFile = 'Input_2018.yml'
    jsonFileName="golden_Json/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"
    sfFileName="DeepCSV_102XSF_V2.csv"

if testfilelist[0].find("UL17") != -1 or testfilelist[0].find("UL2017") != -1:
    year = 2017
    cfgFile = 'Input_2017.yml'
    jsonFileName="golden_Json/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
    sfFileName="DeepCSV_94XSF_V5_B_F.csv"

if testfilelist[0].find("UL16") != -1 or testfilelist[0].find("UL2016") != -1:
    sfFileName="DeepCSV_2016LegacySF_V1.csv"

H4LCppModule = lambda: HZZAnalysisCppProducer(year,cfgFile, isMC, isFSR)
print("Input json file: {}".format(jsonFileName))
print("Input cfg file: {}".format(cfgFile))
print("isMC: {}".format(isMC))
print("isFSR: {}".format(isFSR))

if isMC:
    jetmetCorrector = createJMECorrector(isMC=isMC, dataYear=year, jesUncert="All", jetType = "AK4PFchs")
    fatJetCorrector = createJMECorrector(isMC=isMC, dataYear=year, jesUncert="All", jetType = "AK8PFPuppi")
    btagSF = lambda: btagSFProducer("UL"+str(year), algo="deepjet",selectedWPs=['L','M','T','shape_corr'], sfFileName=sfFileName)
    puidSF = lambda: JetSFMaker("%s" % year)
    # p=PostProcessor(".",testfilelist, None, None,[H4LCppModule(), jetmetCorrector(), fatJetCorrector(), btagSF(), puidSF()], provenance=True,fwkJobReport=False,haddFileName="nano_M125_cpp.root", maxEntries=entriesToRun, prefetch=DownloadFileToLocalThenRun, outputbranchsel="keep_and_drop.txt")
    p=PostProcessor(".",testfilelist, None, None,[H4LCppModule(), jetmetCorrector(), fatJetCorrector(), puidSF()], provenance=True,fwkJobReport=False,haddFileName="nano_M125_cpp.root", maxEntries=entriesToRun, prefetch=DownloadFileToLocalThenRun, outputbranchsel="keep_and_drop.txt")
else:
    jetmetCorrector = createJMECorrector(isMC=isMC, dataYear=year, jesUncert="All", jetType = "AK4PFchs")
    fatJetCorrector = createJMECorrector(isMC=isMC, dataYear=year, jesUncert="All", jetType = "AK8PFPuppi")
    p=PostProcessor(".",testfilelist, None, None,[H4LCppModule(), jetmetCorrector(), fatJetCorrector()], provenance=False, fwkJobReport=False,haddFileName="nano_M125_cpp.root", jsonInput=jsonFileName, maxEntries=entriesToRun, prefetch=DownloadFileToLocalThenRun, outputbranchsel="keep_and_drop_data.txt")

p.run()
print "DONE"
