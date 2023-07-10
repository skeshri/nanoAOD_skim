from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
import os
from Helper import *
ROOT.PyConfig.IgnoreCommandLineOptions = True


# MHT producer, unclean jets only (no lepton overlap cleaning, no jet selection)
class HZZAnalysisCppProducer(Module):
    def __init__(self):
        base = "$CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim"
        ROOT.gSystem.Load("%s/JHUGenMELA/MELA/data/slc7_amd64_gcc700/libJHUGenMELAMELA.so" % base)
        ROOT.gSystem.Load("%s/JHUGenMELA/MELA/data/slc7_amd64_gcc700/libjhugenmela.so" % base)
        ROOT.gSystem.Load("%s/JHUGenMELA/MELA/data/slc7_amd64_gcc700/libmcfm_707.so" % base)
        ROOT.gSystem.Load("%s/JHUGenMELA/MELA/data/slc7_amd64_gcc700/libcollier.so" % base)
        if "/H4LTools_cc.so" not in ROOT.gSystem.GetLibraries():
            print("Load C++ module")
            base = "$CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim"
            if base:
                ROOT.gROOT.ProcessLine(
                    ".L %s/src/H4LTools.cc+O" % base)
            else:
                base = "$CMSSW_BASE//src/PhysicsTools/NanoAODTools"
                ROOT.gSystem.Load("libPhysicsToolsNanoAODTools.so")
                ROOT.gROOT.ProcessLine(
                    ".L %s/interface/H4LTools.h" % base)
        if "/RoccoR_cc.so" not in ROOT.gSystem.GetLibraries():
            base = "$CMSSW_BASE//src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim"
            if base:
                ROOT.gROOT.ProcessLine(
                    ".L %s/src/RoccoR.cc+O" % base)
            else:
                base = "$CMSSW_BASE/src/PhysicsTools/NanoAODTools"
                ROOT.gSystem.Load("libPhysicsToolsNanoAODTools.so")
                ROOT.gROOT.ProcessLine(
                    ".L %s/interface/RoccoR.h" % base)
        #ROOT.gROOT.ProcessLine(
        #    ".L %s/JHUGenMELA/MELA/src/Mela.cc+O" % base)
        self.worker = ROOT.H4LTools()
        pass

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.initReaders(inputTree)  # initReaders must be called in beginFile
        self.out = wrappedOutputTree
        self.out.branch("mass4l",  "F")
        self.out.branch("pT4l",  "F")
        self.out.branch("eta4l",  "F")
        self.out.branch("phi4l",  "F")
        self.out.branch("massZ1",  "F")
        self.out.branch("pTZ1",  "F")
        self.out.branch("etaZ1",  "F")
        self.out.branch("phiZ1",  "F")
        self.out.branch("massZ2",  "F")
        self.out.branch("pTZ2",  "F")
        self.out.branch("etaZ2",  "F")
        self.out.branch("phiZ2",  "F")
        self.out.branch("D_CP",  "F")
        self.out.branch("D_0m",  "F")
        self.out.branch("D_0hp",  "F")
        self.out.branch("D_int",  "F")
        self.out.branch("D_L1",  "F")
        self.out.branch("D_L1Zg",  "F")


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    # this function gets the pointers to Value and ArrayReaders and sets
    # them in the C++ worker class
    def initReaders(self, tree):
        self.nElectron = tree.valueReader("nElectron")
        self.Electron_pt = tree.arrayReader("Electron_pt")
        self.Electron_eta = tree.arrayReader("Electron_eta")
        self.Electron_phi = tree.arrayReader("Electron_phi")
        self.Electron_mass = tree.arrayReader("Electron_mass")
        self.Electron_dxy = tree.arrayReader("Electron_dxy")
        self.Electron_dz = tree.arrayReader("Electron_dz")
        self.Electron_sip3d = tree.arrayReader("Electron_sip3d")
        self.Electron_mvaFall17V2Iso_WP90 = tree.arrayReader("Electron_mvaFall17V2Iso_WP90")
        self.Electron_pdgId = tree.arrayReader("Electron_pdgId")
        self.worker.SetElectrons(self.nElectron, self.Electron_pt, self.Electron_eta, self.Electron_phi, self.Electron_mass,self.Electron_dxy, self.Electron_dz, self.Electron_sip3d,
                                 self.Electron_mvaFall17V2Iso_WP90,self.Electron_pdgId)

        self.nMuon = tree.valueReader("nMuon")
        self.Muon_pt = tree.arrayReader("Muon_pt")
        self.Muon_eta = tree.arrayReader("Muon_eta")
        self.Muon_phi = tree.arrayReader("Muon_phi")
        self.Muon_mass = tree.arrayReader("Muon_mass")
        self.Muon_isGlobal = tree.arrayReader("Muon_isGlobal")
        self.Muon_isTracker = tree.arrayReader("Muon_isTracker")
        self.Muon_dxy = tree.arrayReader("Muon_dxy")
        self.Muon_dz = tree.arrayReader("Muon_dz")
        self.Muon_sip3d = tree.arrayReader("Muon_sip3d")
        self.Muon_ptErr = tree.arrayReader("Muon_ptErr")
        self.Muon_nTrackerLayers = tree.arrayReader("Muon_nTrackerLayers")
        self.Muon_isPFcand = tree.arrayReader("Muon_isPFcand")
        self.Muon_pdgId = tree.arrayReader("Muon_pdgId")
        self.Muon_charge = tree.arrayReader("Muon_charge")
        self.Muon_pfRelIso03_all = tree.arrayReader("Muon_pfRelIso03_all")
        self.Muon_genPartIdx = tree.arrayReader("Muon_genPartIdx")
        self.worker.SetMuons(self.nMuon, self.Muon_pt, self.Muon_eta, self.Muon_phi, self.Muon_mass, self.Muon_isGlobal, self.Muon_isTracker,
                              self.Muon_dxy, self.Muon_dz, self.Muon_sip3d, self.Muon_ptErr, self.Muon_nTrackerLayers, self.Muon_isPFcand, self.Muon_pdgId,self.Muon_charge, self.Muon_pfRelIso03_all, self.Muon_genPartIdx)

        self.nFsrPhoton = tree.valueReader("nFsrPhoton")
        self.FsrPhoton_pt = tree.arrayReader("FsrPhoton_pt")
        self.FsrPhoton_eta = tree.arrayReader("FsrPhoton_eta")
        self.FsrPhoton_phi = tree.arrayReader("FsrPhoton_phi")
        self.FsrPhoton_dROverEt2 = tree.arrayReader("FsrPhoton_dROverEt2")
        self.FsrPhoton_relIso03 = tree.arrayReader("FsrPhoton_relIso03")
        self.FsrPhoton_muonIdx = tree.arrayReader("FsrPhoton_muonIdx")
        self.worker.SetFsrPhotons(self.nFsrPhoton,self.FsrPhoton_dROverEt2,self.FsrPhoton_eta,self.FsrPhoton_phi,self.FsrPhoton_pt,
                                  self.FsrPhoton_relIso03)

        self.nGenPart = tree.valueReader("nGenPart")
        self.GenPart_pt = tree.arrayReader("GenPart_pt")
        self.worker.SetGenParts(self.nGenPart,self.GenPart_pt)

        self.nJet = tree.valueReader("nJet")
        self.Jet_pt = tree.arrayReader("Jet_pt")
        self.Jet_eta = tree.arrayReader("Jet_eta")
        self.Jet_phi = tree.arrayReader("Jet_phi")
        self.Jet_mass = tree.arrayReader("Jet_mass")
        self.Jet_btagDeepC = tree.arrayReader("Jet_btagDeepB")
        self.Jet_jetId = tree.arrayReader("Jet_jetId")
        self.Jet_puId = tree.arrayReader("Jet_puId")
        self.worker.SetJets(self.nJet,self.Jet_pt,self.Jet_eta,self.Jet_phi,self.Jet_mass,self.Jet_btagDeepC,self.Jet_jetId,self.Jet_puId)
        # self._ttreereaderversion must be set AFTER all calls to
        # tree.valueReader or tree.arrayReader
        self._ttreereaderversion = tree._ttreereaderversion

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail,
        go to next event)"""
        # do this check at every event, as other modules might have read
        # further branches
        if event._tree._ttreereaderversion > self._ttreereaderversion:
            self.initReaders(event._tree)
        # do NOT access other branches in python between the check/call to
        # initReaders and the call to C++ worker code



        keepIt = False

        passedTrig=False
        passedFullSelection=False
        passedZ4lSelection=False
        passedQCDcut=False
        passedZ1LSelection=False
        passedZ4lZ1LSelection=False
        passedZ4lZXCRSelection=False
        passedZXCRSelection=False
        passedFiducialSelection=False
        nZXCRFailedLeptons=0
        isMC = True
        passedTrig = PassTrig(event, 2018)
        if (passedTrig==False): return keepIt
        self.worker.MuonPtCorrection(isMC)
        foundZZCandidate = self.worker.ZZSelection_2l2q()
        if (foundZZCandidate):
            keepIt = True
            pTZ1 = self.worker.Z1.Pt()
            etaZ1 = self.worker.Z1.Eta()
            phiZ1 = self.worker.Z1.Phi()
            massZ1 = self.worker.Z1.M()
            pTZ2 = self.worker.Z2.Pt()
            etaZ2 = self.worker.Z2.Eta()
            phiZ2 = self.worker.Z2.Phi()
            massZ2 = self.worker.Z2.M()
            D_CP = self.worker.D_CP
            D_0m = self.worker.D_0m
            D_0hp = self.worker.D_0hp
            D_int = self.worker.D_int
            D_L1 = self.worker.D_L1
            D_L1Zg = self.worker.D_L1Zg

            ZZsystem = ROOT.TLorentzVector()
            ZZsystem = self.worker.Z1 + self.worker.Z2
            pT4l = ZZsystem.Pt()
            eta4l = ZZsystem.Eta()
            phi4l = ZZsystem.Phi()
            mass4l = ZZsystem.M()
            self.out.fillBranch("mass4l",mass4l)
            self.out.fillBranch("pT4l",pT4l)
            self.out.fillBranch("eta4l",eta4l)
            self.out.fillBranch("phi4l",phi4l)
            self.out.fillBranch("massZ1",massZ1)
            self.out.fillBranch("pTZ1",pTZ1)
            self.out.fillBranch("etaZ1",etaZ1)
            self.out.fillBranch("phiZ1",phiZ1)
            self.out.fillBranch("massZ2",massZ2)
            self.out.fillBranch("pTZ2",pTZ2)
            self.out.fillBranch("etaZ2",etaZ2)
            self.out.fillBranch("phiZ2",phiZ2)
            self.out.fillBranch("D_CP",D_CP)
            self.out.fillBranch("D_0m",D_0m)
            self.out.fillBranch("D_0hp",D_0hp)
            self.out.fillBranch("D_int",D_int)
            self.out.fillBranch("D_L1",D_L1)
            self.out.fillBranch("D_L1Zg",D_L1Zg)



        return keepIt


# define modules using the syntax 'name = lambda : constructor' to avoid
# having them loaded when not needed

H4LCppModule = lambda: HZZAnalysisCppProducer()
