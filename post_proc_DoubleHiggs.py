#!/usr/bin/env python
"""This requires few user inputs.

Attributes:
    DownloadFileToLocalThenRun (bool): Set true if you want to download the file to local first then run.
    entriesToRun (int): Number of events to run. If you want to run over full statistics then set this to 0.
    testfile (str): Name of root file to run.
"""

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from HHAnalysisModule import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
from CheckYear_Samples import CheckYear

testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAODv7/GluGluToRadionToHHTo2G4Q_M-3000_narrow_TuneCP5_PSWeights_13TeV-madgraph-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/120000/DB64186F-0685-E445-ACF2-0C039BFD6301.root"


entriesToRun = 0  # 0 if need to run over all entries else put number of entries to run
# Keep DownloadFileToLocalThenRun=True this should reduce the file read error from eos.
DownloadFileToLocalThenRun = True

if testfile.find("SinglePhoton") != -1 or testfile.find("EGamma") != -1 or testfile.find("SingleElectron") != -1 or testfile.find("DoubleMuon") != -1 or testfile.find("MuonEG") != -1 or testfile.find("DoubleEG") != -1:
    isMC = False
    Year, jsonFileName = CheckYear(testfile, isMC)
    print "\n===> Running over ", Year, " data...\n"
    print "===> JSON File: ", jsonFileName
    # jetmetCorrector = createJMECorrector(isMC=isMC, dataYear=Year, jesUncert="All", redojec=True, jetType = "AK4PFchs")
    # fatJetCorrector = createJMECorrector(isMC=isMC, dataYear=Year, jesUncert="All", redojec=True, jetType = "AK8PFPuppi")
    jetmetCorrector = createJMECorrector(isMC=isMC,
                                         dataYear=Year,
                                         jesUncert="Merged",
                                         jetType="AK4PFchs")
    fatJetCorrector = createJMECorrector(isMC=isMC,
                                         dataYear=Year,
                                         jesUncert="Merged",
                                         jetType="AK8PFPuppi")
    p = PostProcessor(".", [testfile], None, None,
                     [HHAnalysisModule(), fatJetCorrector(), jetmetCorrector()],
                      provenance=False,
                      fwkJobReport=False,
                      jsonInput=jsonFileName,
                      maxEntries=entriesToRun,
                      haddFileName="nano_allEvt.root",
                      prefetch=DownloadFileToLocalThenRun)
else:
    print "==> Processing a MC file..."
    isMC = True
    year = CheckYear(testfile, isMC)
    print "\n===> Running over ", year, " MC...\n"
    # jetmetCorrector = createJMECorrector(isMC=isMC, dataYear=year, jesUncert="All", redojec=True, jetType = "AK4PFchs")
    # fatJetCorrector = createJMECorrector(isMC=isMC, dataYear=year, jesUncert="All",  jetType = "AK8PFPuppi")
    fatJetCorrector = createJMECorrector(isMC=isMC,
                                         dataYear=year,
                                         jesUncert="Merged",
                                         jetType="AK8PFPuppi")
    jetmetCorrector = createJMECorrector(isMC=isMC,
                                         dataYear=year,
                                         jesUncert="Merged",
                                         jetType="AK4PFchs")

    p = PostProcessor(".", [testfile], "",
                      "keep_and_drop_inclusive.txt",
                      [HHAnalysisModule(), fatJetCorrector(), jetmetCorrector()],
                      provenance=True,
                      fwkJobReport=False,
                      maxEntries=entriesToRun,
                      haddFileName="nano_allEvt.root",
                      prefetch=DownloadFileToLocalThenRun)

p.run()
print "DONE"
# os.system("ls -lR")
