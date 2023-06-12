#ifndef H4LTools_h
#define H4LTools_h

#include <utility>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <vector>

class H4LTools {
    public:
      H4LTools(){}
      int elePtcut = 7;
      int MuPtcut = 5;
      int sip3dCut = 4;
      float Zmass = 91.1876;
      void SetElectrons(TTreeReaderValue<unsigned> *nElectron_, TTreeReaderArray<float> *Electron_pt_, TTreeReaderArray<float> *Electron_eta_,
                        TTreeReaderArray<float> *Electron_phi_, TTreeReaderArray<float> *Electron_mass_, TTreeReaderArray<float> *Electron_dxy_,TTreeReaderArray<float> *Electron_dz_,
                        TTreeReaderArray<float> *Electron_sip3d_, TTreeReaderArray<bool> *Electron_mvaFall17V2Iso_WP90_, TTreeReaderArray<int> *Electron_pdgId_){
        nElectron = nElectron_; 
        Electron_pt = Electron_pt_; 
        Electron_phi = Electron_phi_;
        Electron_eta = Electron_eta_;
        Electron_mass = Electron_mass_;
        Electron_dxy = Electron_dxy_;
        Electron_dz = Electron_dz_;
        Electron_sip3d = Electron_sip3d_;
        Electron_mvaFall17V2Iso_WP90 = Electron_mvaFall17V2Iso_WP90_;
        Electron_pdgId = Electron_pdgId_;
      }

      void SetMuons(TTreeReaderValue<unsigned> *nMuon_, TTreeReaderArray<float> *Muon_pt_, TTreeReaderArray<float> *Muon_eta_,
                        TTreeReaderArray<float> *Muon_phi_, TTreeReaderArray<float> *Muon_mass_, TTreeReaderArray<bool> *Muon_isGlobal_, TTreeReaderArray<bool> *Muon_isTracker_,
                        TTreeReaderArray<float> *Muon_dxy_, TTreeReaderArray<float> *Muon_dz_,TTreeReaderArray<float> *Muon_sip3d_, TTreeReaderArray<float> *Muon_ptErr_,
                        TTreeReaderArray<int> *Muon_nTrackerLayers_, TTreeReaderArray<bool> *Muon_isPFcand_, TTreeReaderArray<int> *Muon_pdgId_, TTreeReaderArray<float> *Muon_pfRelIso03_all_){
        nMuon = nMuon_; 
        Muon_pt = Muon_pt_; 
        Muon_phi = Muon_phi_;
        Muon_eta = Muon_eta_;
        Muon_mass = Muon_mass_;
        Muon_isGlobal = Muon_isGlobal_;
        Muon_isTracker = Muon_isTracker_;
        Muon_dxy = Muon_dxy_;
        Muon_dz = Muon_dz_;
        Muon_sip3d = Muon_sip3d_;
        Muon_ptErr = Muon_ptErr_;
        Muon_nTrackerLayers = Muon_nTrackerLayers_;
        Muon_isPFcand = Muon_isPFcand_;
        Muon_pdgId = Muon_pdgId_;
        Muon_pfRelIso03_all = Muon_pfRelIso03_all_;
      }

      std::vector<unsigned int> goodLooseElectrons2012();
      std::vector<unsigned int> goodLooseMuons2012();
      std::vector<unsigned int> goodMuons2015_noIso_noPf(std::vector<unsigned int> Muonindex);
      std::vector<unsigned int> goodElectrons2015_noIso_noBdt(std::vector<unsigned int> Electronindex);
      std::vector<bool> passTight_BDT_Id();
      std::vector<bool> passTight_Id();

      bool ZZSelection();
      TLorentzVector Z1;
      TLorentzVector Z2;


    private:

      TTreeReaderValue<unsigned> *nElectron = nullptr;
      TTreeReaderArray<float> *Electron_pt = nullptr;
      TTreeReaderArray<float> *Electron_eta = nullptr;
      TTreeReaderArray<float> *Electron_phi = nullptr;
      TTreeReaderArray<float> *Electron_mass = nullptr;
      TTreeReaderArray<float> *Electron_dxy = nullptr;
      TTreeReaderArray<float> *Electron_dz = nullptr;
      TTreeReaderArray<float> *Electron_sip3d = nullptr;
      TTreeReaderArray<bool> *Electron_mvaFall17V2Iso_WP90 = nullptr;
      TTreeReaderArray<int> *Electron_pdgId = nullptr;

      TTreeReaderValue<unsigned> *nMuon = nullptr;
      TTreeReaderArray<float> *Muon_pt = nullptr;
      TTreeReaderArray<float> *Muon_eta = nullptr;
      TTreeReaderArray<float> *Muon_phi = nullptr;
      TTreeReaderArray<float> *Muon_mass = nullptr;
      TTreeReaderArray<bool> *Muon_isGlobal = nullptr;
      TTreeReaderArray<bool> *Muon_isTracker = nullptr;
      TTreeReaderArray<bool> *Muon_isPFcand = nullptr;
      TTreeReaderArray<int> *Muon_nTrackerLayers = nullptr;
      TTreeReaderArray<float> *Muon_sip3d = nullptr;
      TTreeReaderArray<float> *Muon_dxy = nullptr;
      TTreeReaderArray<float> *Muon_dz = nullptr;
      TTreeReaderArray<float> *Muon_ptErr = nullptr;
      TTreeReaderArray<int> *Muon_pdgId = nullptr;
      TTreeReaderArray<float> *Muon_pfRelIso03_all = nullptr;


      

      

};

#endif

