#include "../interface/H4LTools.h"
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <vector>

std::vector<unsigned int> H4LTools::goodLooseElectrons2012(){
    std::vector<unsigned int> LooseElectronindex;
    unsigned nE = (*nElectron).Get()[0];
    for (unsigned int i=0; i<nE; i++){
        if (((*Electron_pt)[i]>elePtcut)&&(fabs((*Electron_eta)[i])<2.5)){
            LooseElectronindex.push_back(i);
        }
    }

    return LooseElectronindex;
}

std::vector<unsigned int> H4LTools::goodLooseMuons2012(){
    std::vector<unsigned int> LooseMuonindex;
    unsigned nMu = (*nMuon).Get()[0];
    for (unsigned int i=0; i<nMu; i++){
        if ((Muon_Pt_Corrected[i]>MuPtcut)&&(fabs((*Muon_eta)[i])<2.4)&&(((*Muon_isGlobal)[i]||(*Muon_isTracker)[i]||(*Muon_isPFcand)[i]))){
            LooseMuonindex.push_back(i);
        }
    }

    return LooseMuonindex;
}
std::vector<unsigned int> H4LTools::goodMuons2015_noIso_noPf(std::vector<unsigned int> Muonindex){
    std::vector<unsigned int> bestMuonindex;
    for (unsigned int i=0; i<Muonindex.size(); i++){
        if ((Muon_Pt_Corrected[Muonindex[i]]>MuPtcut)&&(fabs((*Muon_eta)[Muonindex[i]])<2.4)&&((*Muon_isGlobal)[Muonindex[i]]||(*Muon_isTracker)[Muonindex[i]])){
            if ((*Muon_sip3d)[Muonindex[i]]<sip3dCut){
                if((fabs((*Muon_dxy)[Muonindex[i]])<0.5)&&(fabs((*Muon_dz)[Muonindex[i]])<1)){
                    bestMuonindex.push_back(Muonindex[i]);
                }
            }
        }
    }
    
    return bestMuonindex;
}
std::vector<unsigned int> H4LTools::goodElectrons2015_noIso_noBdt(std::vector<unsigned int> Electronindex){
    std::vector<unsigned int> bestElectronindex;
    for (unsigned int i=0; i<Electronindex.size(); i++){
        if (((*Electron_pt)[Electronindex[i]])>elePtcut){
            if((*Electron_sip3d)[Electronindex[i]]<sip3dCut){
                if((fabs((*Electron_dxy)[Electronindex[i]])<0.5)&&(fabs((*Electron_dz)[Electronindex[i]])<1)){
                    bestElectronindex.push_back(Electronindex[i]);
                }
            }
        }
    }

    return bestElectronindex;
}
std::vector<bool> H4LTools::passTight_BDT_Id(){
    std::vector<bool> tightid;
    float cutVal,mvaVal;
    cutVal = 1000;
    mvaVal = -1;
    unsigned nE = (*nElectron).Get()[0];
    for (unsigned int i=0; i<nE; i++){
        if((*Electron_pt)[i]<10){
            if(fabs((*Electron_eta)[i])<0.8) cutVal = 0.9044286167;
            if((fabs((*Electron_eta)[i])>=0.8)&&(fabs((*Electron_eta)[i])<1.479)) cutVal = 0.9094166886;
            if(fabs((*Electron_eta)[i])>=1.479) cutVal = 0.9443653660;
        }
        else{
            if(fabs((*Electron_eta)[i])<0.8) cutVal = 0.1968600840;
            if((fabs((*Electron_eta)[i])>=0.8)&&(fabs((*Electron_eta)[i])<1.479)) cutVal = 0.0759172100;
            if(fabs((*Electron_eta)[i])>=1.479) cutVal = -0.5169136775;
        }

        mvaVal = (*Electron_mvaFall17V2Iso_WP90)[i];
        if(mvaVal > cutVal){
            tightid.push_back(true);
        }
        else{
            tightid.push_back(false);
        }
    
    }
    
    return tightid;
    
}
std::vector<bool> H4LTools::passTight_Id(){
    std::vector<bool> tightid;
    unsigned nMu = (*nMuon).Get()[0];
    for (unsigned int i=0; i<nMu; i++){
        if (Muon_Pt_Corrected[i]<200){
            tightid.push_back((*Muon_isPFcand)[i]);
        }
        else{
            tightid.push_back((*Muon_isPFcand)[i]||((((*Muon_ptErr)[i]/Muon_Pt_Corrected[i])<0.3)&&(fabs((*Muon_dxy)[i])<0.2)&&(fabs((*Muon_dz)[i])<0.5)&&((*Muon_nTrackerLayers)[i]>5)));
        }

    }

    return tightid;
}

std::vector<unsigned int> H4LTools::goodFsrPhotons(){
    std::vector<unsigned int> goodFsrPhoton;
    unsigned nFsr = (*nFsrPhoton).Get()[0];
    for (unsigned int i=0; i<nFsr; i++){
        if(((*FsrPhoton_pt)[i]>2)&&(fabs((*FsrPhoton_eta)[i])<2.4)&&((*FsrPhoton_relIso03)[i]<0.8)){
            goodFsrPhoton.push_back(i);
        }
    }
    return goodFsrPhoton;
}

std::vector<unsigned int> H4LTools::SelectedJets(std::vector<unsigned int> ele, std::vector<unsigned int> mu){
    std::vector<unsigned int> goodJets;
    unsigned nJ = (*nJet).Get()[0];
    for(unsigned int i=0;i<nJ;i++){
        int overlaptag=0;
        TLorentzVector jettest;
        jettest.SetPtEtaPhiM((*Jet_pt)[i],(*Jet_eta)[i],(*Jet_phi)[i],(*Jet_mass)[i]);
        for(unsigned int ie=0;ie<ele.size();ie++){
            TLorentzVector eletest;
            eletest.SetPtEtaPhiM((*Electron_pt)[ele[ie]],(*Electron_eta)[ele[ie]],(*Electron_phi)[ele[ie]],(*Electron_mass)[ele[ie]]);
            if(eletest.DeltaR(jettest)<0.4) overlaptag++;
        }
        for(unsigned int im=0;im<mu.size();im++){
            TLorentzVector mutest;
            mutest.SetPtEtaPhiM((*Muon_pt)[mu[im]],(*Muon_eta)[mu[im]],(*Muon_phi)[mu[im]],(*Muon_mass)[mu[im]]);
            if(mutest.DeltaR(jettest)<0.4) overlaptag++;
        }
    
      if(overlaptag==0){
        if(((*Jet_pt)[i]>30)&&(fabs((*Jet_eta)[i])<4.7)){
            //std::cout<<"jetindex: "<<i<<"jetID "<<(*Jet_jetId)[i]<<" puID "<<(*Jet_puId)[i]<<std::endl;
            if(((*Jet_jetId)[i]>0)&&((*Jet_puId)[i]==7))
            {
                goodJets.push_back(i);
            }
        }
      }
    }
    return goodJets;
}

unsigned H4LTools::doFsrRecovery(TLorentzVector Lep){
    // This Function returns the index for the possible FsrPhoton
    unsigned int FsrIdx = 999; //only Idx>0 works, pay attention!
    std::vector<unsigned int> BestFsrPhotons;
    BestFsrPhotons = goodFsrPhotons();
    float dRl,dRlOverPt;
    dRl = 999;
    dRlOverPt = 999;
    for(unsigned int i=0;i<BestFsrPhotons.size();i++){
        TLorentzVector fsrcand;
        fsrcand.SetPtEtaPhiM((*FsrPhoton_pt)[BestFsrPhotons[i]],(*FsrPhoton_eta)[BestFsrPhotons[i]],(*FsrPhoton_phi)[BestFsrPhotons[i]],0);
        float dRlC,dRlOverPtC;
        dRlC = fsrcand.DeltaR(Lep);
        dRlOverPtC = dRl/(fsrcand.Pt()*fsrcand.Pt());
        if ((dRlC<0.5)&&(dRlOverPtC<0.012)){
            if(dRlOverPtC<dRlOverPt){
                dRl = dRlC;
                dRlOverPt = dRlOverPtC;
                FsrIdx = BestFsrPhotons[i];
            }
        }
    }

    return FsrIdx;
    
}

std::vector<TLorentzVector> H4LTools::BatchFsrRecovery(std::vector<TLorentzVector> LepList){
    
    std::vector<TLorentzVector> LepFsrList;
    
    for(unsigned int i=0;i<LepList.size();i++){
        int Fsrtag;
        Fsrtag = doFsrRecovery(LepList[i]);
        if (Fsrtag>900){
            LepFsrList.push_back(LepList[i]);
            continue;
        }
        TLorentzVector FsrPhoton;
        FsrPhoton.SetPtEtaPhiM((*FsrPhoton_pt)[Fsrtag],(*FsrPhoton_eta)[Fsrtag],(*FsrPhoton_phi)[Fsrtag],0);
        TLorentzVector LepFsrRecovery;
        LepFsrRecovery = FsrPhoton + LepList[i];
        LepFsrList.push_back(LepFsrRecovery);
    }
    return LepFsrList;
}

float H4LTools::ApplyRoccoR(bool isMC, int charge, float pt, float eta, float phi, float genPt, float nLayers)
{

    float scale_factor;
    if(isMC && nLayers > 5)
    {
        if(genPt > 0)
            scale_factor = calibrator->kSpreadMC(charge, pt, eta, phi, genPt);
        else{
            TRandom3 rand;
            rand.SetSeed(abs(static_cast<int>(sin(phi)*100000)));
            
            double u1;
            u1 = rand.Uniform(1.);
            scale_factor = calibrator->kSmearMC(charge, pt, eta, phi, nLayers, u1);
        }
    }
    else
        scale_factor = calibrator->kScaleDT(charge, pt, eta, phi);
    
    return scale_factor;
    
}
void H4LTools::MuonPtCorrection(bool isMC){
    unsigned nMu = (*nMuon).Get()[0];
    Muon_Pt_Corrected.clear();
    for (unsigned int i=0; i<nMu; i++){
        float scalefactor;
        scalefactor = ApplyRoccoR(isMC, (*Muon_charge)[i],(*Muon_pt)[i],(*Muon_eta)[i],(*Muon_phi)[i],(*GenPart_pt)[(*Muon_genPartIdx)[i]],(*Muon_nTrackerLayers)[i]);
        Muon_Pt_Corrected.push_back(((*Muon_pt)[i])*scalefactor);
    }
    return;
}

bool H4LTools::ZZSelection(){
    bool foundZZCandidate = false;
 
    std::vector<unsigned int> looseEle,looseMu,bestEle,bestMu, tighteleforjetidx, tightmuforjetidx;
    looseEle = goodLooseElectrons2012();
    looseMu = goodLooseMuons2012();
    bestEle = goodElectrons2015_noIso_noBdt(looseEle);
    bestMu = goodMuons2015_noIso_noPf(looseMu);
    std::vector<unsigned int> Electronindex;
    std::vector<unsigned int> Muonindex;
    Electronindex = bestEle;
    Muonindex = bestMu;
    std::vector<bool> AllEid;
    std::vector<bool> AllMuid;
    AllEid = passTight_BDT_Id();
    AllMuid = passTight_Id();
    std::vector<unsigned int> jetidx;
    for (unsigned int iuj=0;iuj<looseEle.size();iuj++){
        if(AllEid[looseEle[iuj]]) tighteleforjetidx.push_back(looseEle[iuj]);
    }
    for (unsigned int juj=0;juj<looseMu.size();juj++){
        if(AllMuid[looseMu[juj]]) tightmuforjetidx.push_back(looseMu[juj]);
    }
    jetidx = SelectedJets(tighteleforjetidx,tightmuforjetidx);
    

    TLorentzVector z1,z2;
    std::vector<TLorentzVector> Elelist;
    std::vector<TLorentzVector> Mulist;
    std::vector<TLorentzVector> ElelistFsr;
    std::vector<TLorentzVector> MulistFsr;
    std::vector<int> Elechg;
    std::vector<int> Muchg;
    std::vector<float> Muiso;
    std::vector<bool> Eid;
    std::vector<bool> muid;
    
    for(unsigned int ie=0; ie<Electronindex.size();ie++){
        if((*Electron_pdgId)[Electronindex[ie]]>0){
            Elechg.push_back(-1);
        }
        else{
            Elechg.push_back(1);
        }
        TLorentzVector Ele;
        Ele.SetPtEtaPhiM((*Electron_pt)[Electronindex[ie]],(*Electron_eta)[Electronindex[ie]],(*Electron_phi)[Electronindex[ie]],(*Electron_mass)[Electronindex[ie]]);
        Elelist.push_back(Ele);
        Eid.push_back(AllEid[Electronindex[ie]]);
    }

    for(unsigned int imu=0; imu<Muonindex.size();imu++){
        if((*Muon_pdgId)[Muonindex[imu]]>0){
            Muchg.push_back(-1);
        }
        else{
            Muchg.push_back(1);
        }
        TLorentzVector Mu;
        Mu.SetPtEtaPhiM(Muon_Pt_Corrected[Muonindex[imu]],(*Muon_eta)[Muonindex[imu]],(*Muon_phi)[Muonindex[imu]],(*Muon_mass)[Muonindex[imu]]);
        Mulist.push_back(Mu);
        muid.push_back(AllMuid[Muonindex[imu]]);
        Muiso.push_back((*Muon_pfRelIso03_all)[Muonindex[imu]]);
    }
    
    ElelistFsr = BatchFsrRecovery(Elelist);
    MulistFsr = BatchFsrRecovery(Mulist);
    int nTightEle = 0;
    int nTightMu = 0;
    int nTightEleChgSum = 0;
    int nTightMuChgSum = 0;
    
    std::vector<int> TightEleindex;
    std::vector<int> TightMuindex;

    for(unsigned int ae=0; ae<Eid.size();ae++){
        if(Eid[ae]==true){
            nTightEle++;
            TightEleindex.push_back(ae);
            nTightEleChgSum += Elechg[ae];
        }
    }

    for(unsigned int amu=0; amu<muid.size();amu++){
        float RelIsoNoFsr;
        RelIsoNoFsr = Muiso[amu];
        unsigned int FsrMuonidx;
        FsrMuonidx = doFsrRecovery(Mulist[amu]);
        if(FsrMuonidx<900){
            TLorentzVector fsrmuon;
            fsrmuon.SetPtEtaPhiM((*FsrPhoton_pt)[FsrMuonidx],(*FsrPhoton_eta)[FsrMuonidx],(*FsrPhoton_phi)[FsrMuonidx],0);
            if(Mulist[amu].DeltaR(fsrmuon)>0.01){
              RelIsoNoFsr = RelIsoNoFsr - (*FsrPhoton_pt)[FsrMuonidx]/Mulist[amu].Pt();  
            }
        }
        if((muid[amu]==true)&&(RelIsoNoFsr<0.35)){
            nTightMu++;
            TightMuindex.push_back(amu);
            nTightMuChgSum += Muchg[amu];
        }
    }

   
    if((nTightMu+nTightEle)<4){
        return foundZZCandidate;
    } 
    
    if((abs(nTightEleChgSum)+abs(nTightMuChgSum))>(nTightMu+nTightEle-4)){
        return foundZZCandidate;
    } 
    

    //Find Z candidates
    std::vector<TLorentzVector> Zlist;
    std::vector<int> Zflavor; //mu->13, e->11
    std::vector<int> Zlep1index;
    std::vector<int> Zlep2index;
    std::vector<float> Zlep1pt;
    std::vector<float> Zlep1eta;
    std::vector<float> Zlep1phi;
    std::vector<float> Zlep1mass;
    std::vector<float> Zlep1chg;
    std::vector<float> Zlep2pt;
    std::vector<float> Zlep2eta;
    std::vector<float> Zlep2phi;
    std::vector<float> Zlep2mass;
    std::vector<float> Zlep2chg;
    std::vector<float> Zlep1ptNoFsr;
    std::vector<float> Zlep1etaNoFsr;
    std::vector<float> Zlep1phiNoFsr;
    std::vector<float> Zlep1massNoFsr;
    std::vector<float> Zlep2ptNoFsr;
    std::vector<float> Zlep2etaNoFsr;
    std::vector<float> Zlep2phiNoFsr;
    std::vector<float> Zlep2massNoFsr;

    if(TightEleindex.size()>1){
        for(unsigned int ke=0; ke<(TightEleindex.size()-1);ke++){
            for(unsigned int je=ke+1;je<TightEleindex.size();je++){
                if ((Elechg[TightEleindex[ke]]+Elechg[TightEleindex[je]])==0){
                    TLorentzVector Zcan;
                    Zcan = ElelistFsr[TightEleindex[ke]] + ElelistFsr[TightEleindex[je]];
                    if((Zcan.M()>12)&&(Zcan.M()<120)){
                        Zlist.push_back(Zcan);
                        Zlep1index.push_back(TightEleindex[ke]);
                        Zlep2index.push_back(TightEleindex[je]);
                        Zflavor.push_back(11);
                        Zlep1pt.push_back(ElelistFsr[TightEleindex[ke]].Pt());
                        Zlep2pt.push_back(ElelistFsr[TightEleindex[je]].Pt());
                        Zlep1eta.push_back(ElelistFsr[TightEleindex[ke]].Eta());
                        Zlep2eta.push_back(ElelistFsr[TightEleindex[je]].Eta());
                        Zlep1phi.push_back(ElelistFsr[TightEleindex[ke]].Phi());
                        Zlep2phi.push_back(ElelistFsr[TightEleindex[je]].Phi());
                        Zlep1mass.push_back(ElelistFsr[TightEleindex[ke]].M());
                        Zlep2mass.push_back(ElelistFsr[TightEleindex[je]].M());
                        Zlep1ptNoFsr.push_back(Elelist[TightEleindex[ke]].Pt());
                        Zlep2ptNoFsr.push_back(Elelist[TightEleindex[je]].Pt());
                        Zlep1etaNoFsr.push_back(Elelist[TightEleindex[ke]].Eta());
                        Zlep2etaNoFsr.push_back(Elelist[TightEleindex[je]].Eta());
                        Zlep1phiNoFsr.push_back(Elelist[TightEleindex[ke]].Phi());
                        Zlep2phiNoFsr.push_back(Elelist[TightEleindex[je]].Phi());
                        Zlep1massNoFsr.push_back(Elelist[TightEleindex[ke]].M());
                        Zlep2massNoFsr.push_back(Elelist[TightEleindex[je]].M());
                        Zlep1chg.push_back(Elechg[TightEleindex[ke]]);
                        Zlep2chg.push_back(Elechg[TightEleindex[je]]);
                    }
                }
            }
        }
    }

    if(TightMuindex.size()>1){
        for(unsigned int kmu=0; kmu<(TightMuindex.size()-1);kmu++){
            for(unsigned int jmu=kmu+1;jmu<TightMuindex.size();jmu++){
                if ((Muchg[TightMuindex[kmu]]+Muchg[TightMuindex[jmu]])==0){
                    TLorentzVector Zcan;
                    Zcan = MulistFsr[TightMuindex[kmu]] + MulistFsr[TightMuindex[jmu]];
                    if((Zcan.M()>12)&&(Zcan.M()<120)){
                        Zlist.push_back(Zcan);
                        Zlep1index.push_back(TightMuindex[kmu]);
                        Zlep2index.push_back(TightMuindex[jmu]);
                        Zflavor.push_back(13);
                        Zlep1pt.push_back(MulistFsr[TightMuindex[kmu]].Pt());
                        Zlep2pt.push_back(MulistFsr[TightMuindex[jmu]].Pt());
                        Zlep1eta.push_back(MulistFsr[TightMuindex[kmu]].Eta());
                        Zlep2eta.push_back(MulistFsr[TightMuindex[jmu]].Eta());
                        Zlep1phi.push_back(MulistFsr[TightMuindex[kmu]].Phi());
                        Zlep2phi.push_back(MulistFsr[TightMuindex[jmu]].Phi());
                        Zlep1mass.push_back(MulistFsr[TightMuindex[kmu]].M());
                        Zlep2mass.push_back(MulistFsr[TightMuindex[jmu]].M());
                        Zlep1ptNoFsr.push_back(Mulist[TightMuindex[kmu]].Pt());
                        Zlep2ptNoFsr.push_back(Mulist[TightMuindex[jmu]].Pt());
                        Zlep1etaNoFsr.push_back(Mulist[TightMuindex[kmu]].Eta());
                        Zlep2etaNoFsr.push_back(Mulist[TightMuindex[jmu]].Eta());
                        Zlep1phiNoFsr.push_back(Mulist[TightMuindex[kmu]].Phi());
                        Zlep2phiNoFsr.push_back(Mulist[TightMuindex[jmu]].Phi());
                        Zlep1massNoFsr.push_back(Mulist[TightMuindex[kmu]].M());
                        Zlep2massNoFsr.push_back(Mulist[TightMuindex[jmu]].M());
                        Zlep1chg.push_back(Muchg[TightMuindex[kmu]]);
                        Zlep2chg.push_back(Muchg[TightMuindex[jmu]]);
                    }
                }
            }
        }
    }
    
    if(Zlist.size()<2){
        return foundZZCandidate;
    }
    
    //Find ZZ candidate
    std::vector<int> Z1CanIndex;
    std::vector<int> Z2CanIndex;

    for (unsigned int m=0; m<(Zlist.size()-1); m++){
        for (unsigned int n=m+1; n<Zlist.size(); n++){
            if (Zflavor[m]==Zflavor[n]){
               if ((Zlep1index[m] == Zlep1index[n])||(Zlep2index[m] == Zlep1index[n])) continue;  //non-overlapping
               if ((Zlep1index[m] == Zlep2index[n])||(Zlep2index[m] == Zlep2index[n])) continue;
            }
            if (Zlist[m].DeltaR(Zlist[n])<0.02) continue; //ghost removal
            bool nPassPt20;
            int nPassPt10;
            nPassPt20 = (Zlep1pt[m]>20) || (Zlep2pt[m]>20) || (Zlep1pt[n]>20) || (Zlep2pt[n]>20);
            nPassPt10 = 0;
            if (Zlep1pt[m]>10) nPassPt10 += 1; 
            if (Zlep2pt[m]>10) nPassPt10 += 1; 
            if (Zlep1pt[n]>10) nPassPt10 += 1; 
            if (Zlep2pt[n]>10) nPassPt10 += 1; 
            if (nPassPt10 < 2) continue;
            if (nPassPt20 == false) continue; //lep Pt requirements

            if ((Zlep1chg[m]+Zlep1chg[n])==0){
                TLorentzVector lepA,lepB,lepAB;
                lepA.SetPtEtaPhiM(Zlep1ptNoFsr[m],Zlep1etaNoFsr[m],Zlep1phiNoFsr[m],Zlep1massNoFsr[m]);
                lepB.SetPtEtaPhiM(Zlep1ptNoFsr[n],Zlep1etaNoFsr[n],Zlep1phiNoFsr[n],Zlep1massNoFsr[n]);
                lepAB = lepA + lepB;
                if(lepAB.M()<4) continue;  //QCD Supressionas
            }
            if ((Zlep1chg[m]+Zlep2chg[n])==0){
                TLorentzVector lepA,lepB,lepAB;
                lepA.SetPtEtaPhiM(Zlep1ptNoFsr[m],Zlep1etaNoFsr[m],Zlep1phiNoFsr[m],Zlep1massNoFsr[m]);
                lepB.SetPtEtaPhiM(Zlep2ptNoFsr[n],Zlep2etaNoFsr[n],Zlep2phiNoFsr[n],Zlep2massNoFsr[n]);
                lepAB = lepA + lepB;
                if(lepAB.M()<4) continue;
            }
            if ((Zlep2chg[m]+Zlep1chg[n])==0){
                TLorentzVector lepA,lepB,lepAB;
                lepA.SetPtEtaPhiM(Zlep2ptNoFsr[m],Zlep2etaNoFsr[m],Zlep2phiNoFsr[m],Zlep2massNoFsr[m]);
                lepB.SetPtEtaPhiM(Zlep1ptNoFsr[n],Zlep1etaNoFsr[n],Zlep1phiNoFsr[n],Zlep1massNoFsr[n]);
                lepAB = lepA + lepB;
                if(lepAB.M()<4) continue;
            }
            if ((Zlep2chg[m]+Zlep2chg[n])==0){
                TLorentzVector lepA,lepB,lepAB;
                lepA.SetPtEtaPhiM(Zlep2ptNoFsr[m],Zlep2etaNoFsr[m],Zlep2phiNoFsr[m],Zlep2massNoFsr[m]);
                lepB.SetPtEtaPhiM(Zlep1ptNoFsr[n],Zlep1etaNoFsr[n],Zlep1phiNoFsr[n],Zlep1massNoFsr[n]);
                lepAB = lepA + lepB;
                if(lepAB.M()<4) continue;
            }

            if ((Zlist[m].M()<40) && (Zlist[n].M()<40))  continue; //Z1 mass

            TLorentzVector zZ1,zZ2;
            if (fabs(Zlist[m].M()-Zmass)<fabs(Zlist[n].M()-Zmass)){
                zZ1 = Zlist[m];
                zZ2 = Zlist[n];
            }
            else{
                zZ1 = Zlist[n];
                zZ2 = Zlist[m];
            }    
            
            bool passSmartCut = true;
            if (Zflavor[m]==Zflavor[n]){
                TLorentzVector Za,Zb,lepM1,lepM2,lepN1,lepN2;
                int lepM1chg,lepM2chg,lepN1chg,lepN2chg;
                lepM1.SetPtEtaPhiM(Zlep1pt[m],Zlep1eta[m],Zlep1phi[m],Zlep1mass[m]);
                lepM2.SetPtEtaPhiM(Zlep2pt[m],Zlep2eta[m],Zlep2phi[m],Zlep2mass[m]);
                lepN1.SetPtEtaPhiM(Zlep1pt[n],Zlep1eta[n],Zlep1phi[n],Zlep1mass[n]);
                lepN2.SetPtEtaPhiM(Zlep2pt[n],Zlep2eta[n],Zlep2phi[n],Zlep2mass[n]);
                lepM1chg = Zlep1chg[m];
                lepM2chg = Zlep2chg[m];
                lepN1chg = Zlep1chg[n];
                lepN2chg = Zlep2chg[n];
                if(lepM1chg == lepN1chg){
                    Za = lepM1 + lepN2;
                    Zb = lepN1 + lepM2;
                }
                    
                else{
                    Za = lepM1 + lepN1;
                    Zb = lepN2 + lepM2;
                }
                if (fabs(Za.M()-Zmass)<fabs(Zb.M()-Zmass)){
                    if ( (fabs(Za.M()-Zmass)<abs(zZ1.M()-Zmass)) && (Zb.M()<12) ) passSmartCut=false;
                }        
                                    
                else{
                    if ( (fabs(Zb.M()-Zmass)<fabs(zZ1.M()-Zmass)) && (Za.M()<12) ) passSmartCut=false;
                }
            }
            if (passSmartCut==false) continue ;
            if (zZ1.M()+zZ2.M()<70) continue;
            foundZZCandidate = true;
            if(Zlist[m].M()>Zlist[n].M()){
                Z1CanIndex.push_back(m);
                Z2CanIndex.push_back(n);
            }
                
            else{
                Z1CanIndex.push_back(n);
                Z2CanIndex.push_back(m);
            }
                    
           
        }
    }
    if(foundZZCandidate == false){
        return foundZZCandidate;
    }
    int Z1index,Z2index; 
    Z1index = Z1CanIndex[0];
    Z2index = Z2CanIndex[0];
    float Z2Ptsum;
    Z2Ptsum = Zlep1pt[Z2index] + Zlep2pt[Z2index];
    if(Z1CanIndex.size()>1){
        for(unsigned int iz=0;iz<Z1CanIndex.size();iz++){
            if (Z1index==Z1CanIndex[iz]){
                if((Zlep1pt[Z2CanIndex[iz]] + Zlep2pt[Z2CanIndex[iz]])>Z2Ptsum){
                    Z1index = Z1CanIndex[iz];
                    Z2index = Z2CanIndex[iz];
                    Z2Ptsum = Zlep1pt[Z2index] + Zlep2pt[Z2index];
                }
            }
            if(fabs(Zlist[Z1CanIndex[iz]].M()-Zmass)<fabs(Zlist[Z1index].M()-Zmass)){
                Z1index = Z1CanIndex[iz];
                Z2index = Z2CanIndex[iz];
                Z2Ptsum = Zlep1pt[Z2index] + Zlep2pt[Z2index];
            }
        }
    }
       
    
    Z1 = Zlist[Z1index];
    Z2 = Zlist[Z2index];
    
    TLorentzVector ZZsystem;
    ZZsystem = Z1+Z2;
    float massZZ;
    massZZ = ZZsystem.M();
         

    unsigned int jet1index, jet2index;
    if(jetidx.size()>0)
    {
        if(jetidx.size()==1)
        {
            jet1index = jetidx[0];
        }
        if(jetidx.size()==2)
        {
            jet1index = jetidx[0];
            jet2index = jetidx[1];
            if((*Jet_pt)[jetidx[1]]>(*Jet_pt)[jetidx[0]])
            {
                jet1index = jetidx[1];
                jet2index = jetidx[0];
            }
        }
        if(jetidx.size()>2)
        {
            jet1index = jetidx[0];
            jet2index = jetidx[1];
            if((*Jet_pt)[jetidx[1]]>(*Jet_pt)[jetidx[0]])
            {
                jet1index = jetidx[1];
                jet2index = jetidx[0];
            }
            for (unsigned int pj=2;pj<jetidx.size();pj++){
                if(((*Jet_pt)[jetidx[pj]]>jet1index)&&((*Jet_pt)[jetidx[pj]]>jet2index)){
                    jet1index = jetidx[pj];
                }
                if((*Jet_pt)[jetidx[pj]]>jet2index){
                    jet2index = jetidx[pj];
                }
            }
        }
    }
    TLorentzVector Jet1,Jet2;
    SimpleParticleCollection_t associated;
    if(jetidx.size()>0){
        Jet1.SetPtEtaPhiM((*Jet_pt)[jet1index],(*Jet_eta)[jet1index],(*Jet_phi)[jet1index],(*Jet_mass)[jet1index]);
        associated.push_back(SimpleParticle_t(0, Jet1));
        if(jetidx.size()>1){
           Jet2.SetPtEtaPhiM((*Jet_pt)[jet2index],(*Jet_eta)[jet2index],(*Jet_phi)[jet2index],(*Jet_mass)[jet2index]);
           associated.push_back(SimpleParticle_t(0, Jet2));
        }
    }
    
    
        
    SimpleParticleCollection_t daughters;
    TLorentzVector Lep1,Lep2,Lep3,Lep4;
    Lep1.SetPtEtaPhiM(Zlep1pt[Z1index],Zlep1eta[Z1index],Zlep1phi[Z1index],Zlep1mass[Z1index]);
    Lep2.SetPtEtaPhiM(Zlep2pt[Z1index],Zlep2eta[Z1index],Zlep2phi[Z1index],Zlep2mass[Z1index]);
    Lep3.SetPtEtaPhiM(Zlep1pt[Z2index],Zlep1eta[Z2index],Zlep1phi[Z2index],Zlep1mass[Z2index]);
    Lep4.SetPtEtaPhiM(Zlep2pt[Z2index],Zlep2eta[Z2index],Zlep2phi[Z2index],Zlep2mass[Z2index]);
    
    daughters.push_back(SimpleParticle_t((-1)*Zflavor[Z1index]*Zlep1chg[Z1index], Lep1));
    daughters.push_back(SimpleParticle_t((-1)*Zflavor[Z1index]*Zlep2chg[Z1index], Lep2));
    daughters.push_back(SimpleParticle_t((-1)*Zflavor[Z2index]*Zlep1chg[Z2index], Lep3));
    daughters.push_back(SimpleParticle_t((-1)*Zflavor[Z2index]*Zlep2chg[Z2index], Lep4));
    me_0plus_JHU=999.0; me_qqZZ_MCFM=999.0; p0plus_m4l=999.0; bkg_m4l=999.0; D_bkg_kin=999.0; D_bkg=999.0;
    D_bkg_kin_vtx_BS=999.0;
    
    p0minus_VAJHU=999.0; pg1g4_VAJHU=999.0; Dgg10_VAMCFM=999.0, D_g4=999.0; D_g1g4=999.0; D_0m=999.0; D_CP=999.0; D_0hp=999; D_int=999.0;D_L1=999.0; D_L1_int=999.0; D_L1Zg=999.0; D_L1Zgint=999.0;
    p0plus_VAJHU=9999.0; p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen=999.0; pDL1_VAJHU=999.0; pD_L1Zgint=999.0; p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen=999.0; p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen=999.0, p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen=999.0, p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen=999.0, p0plus_VAJHU=999.0; 
    
    mela->setInputEvent(&daughters, &associated, 0, 0);
    mela->setCurrentCandidateFromIndex(0);
    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(me_0plus_JHU, true);
    
    mela->setProcess(TVar::H0minus, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p0minus_VAJHU, true);
    // additional probabilities   GG_SIG_ghg2_1_ghz2_1_JHUGen
    mela->setProcess(TVar::H0hplus, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p0plus_VAJHU, true); 
    
    // p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen, Couplings:ghg2=1,0;ghz1=1,0;ghz2=1,0 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghz2_1_JHUGen
    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    (mela->selfDHggcoupl)[0][gHIGGS_GG_2][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_VV_1][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_VV_2][0]=1.;
    mela->computeP(p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen, true);    //FIXME 

    p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen -= p0plus_VAJHU+me_0plus_JHU;  


    // p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen
    // Couplings:ghg2=1,0;ghz1_prime2=10000,0
    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    (mela->selfDHggcoupl)[0][gHIGGS_GG_2][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_VV_1_PRIME2][0]=10000.;
    // (mela->selfDHzzcoupl)[0][3][0]=1.;
    mela->computeP(p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen, true);    //FIXME

	// p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen/1e8, ghg2=1,0;ghz1=1,0;ghz1_prime2=10000,0, Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen
	mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    (mela->selfDHggcoupl)[0][gHIGGS_GG_2][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_VV_1][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_VV_1_PRIME2][0]=10000.;
    mela->computeP(p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen, true);    //FIXME 
    p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen -= p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen+me_0plus_JHU;

	// p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen, ghg2=1,0;ghzgs1_prime2=10000,0	
	mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    (mela->selfDHggcoupl)[0][gHIGGS_GG_2][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_ZA_1_PRIME2][0]=10000.;
    // (mela->selfDHzzcoupl)[0][3][0]=1.;
    mela->computeP(p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen, true);    //FIXME

	// p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen, ghg2=1,0;ghzgs1_prime2=10000,0
	mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    (mela->selfDHggcoupl)[0][gHIGGS_GG_2][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_VV_1][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_ZA_1_PRIME2][0]=10000.;
    mela->computeP(p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen, true);    //FIXME

	p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen -= me_0plus_JHU+p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen;

	// pg1g4_VAJHU=0.0;
    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    (mela->selfDHggcoupl)[0][0][0]=1.;
    (mela->selfDHzzcoupl)[0][0][0]=1.;
    (mela->selfDHzzcoupl)[0][3][0]=1.;
    mela->computeP(pg1g4_VAJHU, true);

    pg1g4_VAJHU -= me_0plus_JHU+p0minus_VAJHU;

    mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
    mela->computeP(me_qqZZ_MCFM, true);
                    
    mela->computeD_gg(TVar::MCFM, TVar::D_gg10, Dgg10_VAMCFM);
                    
    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
    mela->computePM4l(TVar::SMSyst_None, p0plus_m4l);
                    
    mela->setProcess(TVar::bkgZZ, TVar::JHUGen, TVar::ZZGG);
    mela->computePM4l(TVar::SMSyst_None, bkg_m4l);
    D_0m = me_0plus_JHU / (me_0plus_JHU + (p0minus_VAJHU * pow(getDg4Constant(massZZ),2)));
	D_CP = pg1g4_VAJHU / (2 * sqrt(me_0plus_JHU * p0minus_VAJHU ));
    D_0hp = me_0plus_JHU / (me_0plus_JHU + (p0plus_VAJHU * pow(getDg2Constant(massZZ),2)));
	D_int = p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen / (2 * sqrt(me_0plus_JHU * p0plus_VAJHU));
    D_L1 = me_0plus_JHU / (me_0plus_JHU + ((p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen/1e8) * pow(getDL1Constant(massZZ),2)));
	D_L1Zg = me_0plus_JHU / (me_0plus_JHU + ((p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen/1e8) * pow(getDL1ZgsConstant(massZZ),2)));
    mela->resetInputEvent();
    return foundZZCandidate;
    

}

float H4LTools::getDg4Constant(float ZZMass){
    return spline_g4->Eval(ZZMass);
}

float H4LTools::getDg2Constant(float ZZMass){
    return spline_g2->Eval(ZZMass);
}

float H4LTools::getDL1Constant(float ZZMass){
    return spline_L1->Eval(ZZMass);
}

float H4LTools::getDL1ZgsConstant(float ZZMass){
    return spline_L1Zgs->Eval(ZZMass);
}