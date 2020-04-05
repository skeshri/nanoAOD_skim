import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class wvAnalysisProducer(Module):
    def __init__(self):
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        #self.out.branch("EventMass",  "F");
	#pass
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        jets = Collection(event, "Jet")
        fatJets = Collection(event, "FatJet")
        keepIt = True
        eventElectrons = 0
        eventMuons = 0
        eventJets = 0
        eventFatJets = 0
        
        for lep in muons :
            if lep.softId and abs(lep.dz)<0.1 and abs(lep.dxy < 0.02) and lep.pfRelIso04_all < 0.4 and lep.pt > 20:
                eventMuons += 1
        for lep in electrons :
            if lep.cutBased == 2 and lep.pt > 20 :
                eventElectrons += 1
        for jet in jets :
            if jet.pt > 20:
               eventJets += 1
        for fatjet in fatJets :
            if fatjet.pt > 20:
               eventFatJets += 1
        
        if not ( ((eventElectrons >= 1 or eventMuons >=1) and (eventMuons+eventElectrons < 3)) and ( (eventJets >= 2 and eventFatJets >= 1) or (eventJets >= 4 and eventFatJets == 0) )  ):
            keepIt = False
        
        return keepIt


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
wvAnalysisModule = lambda : wvAnalysisProducer() #(jetSelection= lambda j : j.pt > 30) 
 
