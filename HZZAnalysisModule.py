import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class HZZAnalysisProducer(Module):
    def __init__(self):
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("ZLep_type",  "F");
        self.out.branch("ZLep_mass",  "F");
        self.out.branch("ZLep_pt",  "F");
        self.out.branch("ZLep_eta",  "F");
        self.out.branch("ZLep_phi",  "F");
        self.out.branch("ZHad_mass",  "F");
        self.out.branch("ZHad_pt",  "F");
        self.out.branch("ZHad_eta",  "F");
        self.out.branch("ZHad_phi",  "F");
        self.out.branch("BoostedZHad",  "I");
        self.out.branch("ZZ_mass",  "F");
        self.out.branch("ZZ_pt",  "F");
        self.out.branch("ZZ_eta",  "F");
        self.out.branch("ZZ_phi",  "F");
        """process event, return True (go to next module) or False (fail, go to next event)"""
        #pass
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        """nanoAOD skimming is done considering the final events selection
        for the vv semileptonic final state.

        For the vv semileptonic final state we should be either one or two
        tight leptons, either one Fat jet and two small radius jet or four
        small radius jets.

        Arguments:
            event {instance of event} -- instance of event

        Returns:
            boolean -- if the event passes skimming then it returns true and
                       go to the next module else returns false and go to
                       the next event.
        """
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        jets = Collection(event, "Jet")
        fatJets = Collection(event, "FatJet")
        keepIt = True

        SelectedElectrons = [x for x in electrons if x.pt > 15 and x.cutBased>3]
        SelectedMuons = [x for x in muons if x.pt > 15 and x.tightId>=1]
        SelectedJets = [x for x in jets if x.pt > 20]
        SelectedFatJets = [x for x in fatJets if x.pt > 200]

        nLeptons = len(SelectedElectrons)+len(SelectedMuons)

        if not(nLeptons>=2 and (len(SelectedFatJets)>=1 or len(SelectedJets)>=2)):
            keepIt = False
        else:
            keepIt = True

        ZLep_type = -1
        VLeptons = []
        BoostedZHad = -1
        LV_ZLep = ROOT.TLorentzVector(0,0,0,0)
        LV_ZHad = ROOT.TLorentzVector(0,0,0,0)

        if nLeptons == 2 and len(SelectedFatJets) >=1:
            BoostedZHad = 1
            if len(SelectedElectrons) ==2:
                if SelectedElectrons[0].charge * SelectedElectrons[1].charge < 0:
                    ZLep_type = 0
                    VLeptons = [SelectedElectrons[0],SelectedElectrons[1]]
            if len(SelectedMuons) ==2:
                if SelectedMuons[0].charge * SelectedMuons[1].charge < 0:
                    ZLep_type = 1
                    VLeptons = [SelectedMuons[0],SelectedMuons[1]]


            for Lepton in VLeptons:
                LV_lepton = ROOT.TLorentzVector()
                LV_lepton.SetPtEtaPhiM(Lepton.pt,Lepton.eta,Lepton.phi,Lepton.mass)
                LV_ZLep = LV_ZLep + LV_lepton

            # Use only leading pT fat jet
            LV_ZHad.SetPtEtaPhiM(SelectedFatJets[0].pt,SelectedFatJets[0].eta,SelectedFatJets[0].phi,SelectedFatJets[0].msoftdrop)

        elif nLeptons == 2 and len(SelectedJets) >=2:
            BoostedZHad = 0
            if len(SelectedElectrons) ==2:
                if SelectedElectrons[0].charge * SelectedElectrons[1].charge < 0:
                    ZLep_type = 0
                    VLeptons = [SelectedElectrons[0],SelectedElectrons[1]]
            if len(SelectedMuons) ==2:
                if SelectedMuons[0].charge * SelectedMuons[1].charge < 0:
                    ZLep_type = 1
                    VLeptons = [SelectedMuons[0],SelectedMuons[1]]

            LV_ZLep = ROOT.TLorentzVector()
            for Lepton in VLeptons:
                LV_lepton = ROOT.TLorentzVector()
                LV_lepton.SetPtEtaPhiM(Lepton.pt,Lepton.eta,Lepton.phi,Lepton.mass)
                LV_ZLep = LV_ZLep + LV_lepton

            # Choose two leading pT jets
            LV_jets = ROOT.TLorentzVector()
            LV_jets.SetPtEtaPhiM(SelectedJets[0].pt,SelectedJets[0].eta,SelectedJets[0].phi,SelectedJets[0].mass)
            LV_ZHad = LV_jets
            LV_jets.SetPtEtaPhiM(SelectedJets[1].pt,SelectedJets[1].eta,SelectedJets[1].phi,SelectedJets[1].mass)
            LV_ZHad = LV_ZHad + LV_jets

        self.out.fillBranch("ZLep_type",ZLep_type)
        self.out.fillBranch("ZLep_mass",LV_ZLep.M())
        self.out.fillBranch("ZLep_pt",LV_ZLep.Pt())
        self.out.fillBranch("ZLep_eta",LV_ZLep.Eta())
        self.out.fillBranch("ZLep_phi",LV_ZLep.Phi())
        self.out.fillBranch("ZHad_mass",LV_ZHad.M())
        self.out.fillBranch("ZHad_pt",LV_ZHad.Pt())
        self.out.fillBranch("ZHad_eta",LV_ZHad.Eta())
        self.out.fillBranch("ZHad_phi",LV_ZHad.Phi())
        self.out.fillBranch("BoostedZHad",BoostedZHad)
        self.out.fillBranch("ZZ_mass",(LV_ZLep+LV_ZHad).M())
        self.out.fillBranch("ZZ_pt",(LV_ZLep+LV_ZHad).Pt())
        self.out.fillBranch("ZZ_eta",(LV_ZLep+LV_ZHad).Eta())
        self.out.fillBranch("ZZ_phi",(LV_ZLep+LV_ZHad).Phi())

        return keepIt


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
HZZAnalysisModule = lambda : HZZAnalysisProducer() #(jetSelection= lambda j : j.pt > 30)