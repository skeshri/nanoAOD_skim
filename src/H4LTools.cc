#include "../interface/H4LTools.h"
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <vector>

std::vector<unsigned int> H4LTools::goodLooseElectrons2012(){
    std::vector<unsigned int> LooseElectronindex;
    for (unsigned int i=0; i<Electron_pt.size(); i++){
        if ((Electron_pt[i]>elePtcut)&&(fabs(Electron_eta[i])<eleEtacut)){
            LooseElectronindex.push_back(i);
        }
    }

    return LooseElectronindex;
}

std::vector<unsigned int> H4LTools::goodLooseMuons2012(){
    std::vector<unsigned int> LooseMuonindex;
    for (unsigned int i=0; i<Muon_eta.size(); i++){
        if ((Muon_pt[i]>MuPtcut)&&(fabs(Muon_eta[i])<MuEtacut)&&((Muon_isGlobal[i]||Muon_isTracker[i]||Muon_isPFcand[i]))){
            LooseMuonindex.push_back(i);
        }
    }

    return LooseMuonindex;
}
std::vector<unsigned int> H4LTools::goodMuons2015_noIso_noPf(std::vector<unsigned int> Muonindex){
    std::vector<unsigned int> bestMuonindex;
    for (unsigned int i=0; i<Muonindex.size(); i++){
        if ((Muon_pt[Muonindex[i]]>MuPtcut)&&(fabs(Muon_eta[Muonindex[i]])<MuEtacut)&&(Muon_isGlobal[Muonindex[i]]||Muon_isTracker[Muonindex[i]])){
            if (Muon_sip3d[Muonindex[i]]<Musip3dCut){
                if((fabs(Muon_dxy[Muonindex[i]])<MuLoosedxycut)&&(fabs(Muon_dz[Muonindex[i]])<MuLoosedzcut)){
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
        if ((Electron_pt[Electronindex[i]])>elePtcut){
            if(Electron_sip3d[Electronindex[i]]<elesip3dCut){
                if((fabs(Electron_dxy[Electronindex[i]])<eleLoosedxycut)&&(fabs(Electron_dz[Electronindex[i]])<eleLoosedzcut)){
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
    for (unsigned int i=0; i<Electron_pt.size(); i++){
        if(Electron_pt[i]<10){
            if(fabs(Electron_eta[i])<0.8) cutVal = eleBDTWPLELP;
            if((fabs(Electron_eta[i])>=0.8)&&(fabs(Electron_eta[i])<1.479)) cutVal = eleBDTWPMELP;
            if(fabs(Electron_eta[i])>=1.479) cutVal = eleBDTWPHELP;
        }
        else{
            if(fabs(Electron_eta[i])<0.8) cutVal = eleBDTWPLEHP;
            if((fabs(Electron_eta[i])>=0.8)&&(fabs(Electron_eta[i])<1.479)) cutVal = eleBDTWPMEHP;
            if(fabs(Electron_eta[i])>=1.479) cutVal = eleBDTWPHEHP;
        }

        mvaVal = Electron_mvaFall17V2Iso[i];
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
    //unsigned nMu = (*nMuon).Get()[0];
    for (unsigned int i=0; i<Muon_pt.size(); i++){
        if (Muon_pt[i]<MuHighPtBound){
            tightid.push_back(Muon_isPFcand[i]);
        }
        else{
            tightid.push_back(Muon_isPFcand[i]||(((Muon_ptErr[i]/Muon_pt[i])<MuTightpTErrorcut)&&(fabs(Muon_dxy[i])<MuTightdxycut)&&(fabs(Muon_dz[i])<MuTightdzcut)&&(Muon_nTrackerLayers[i]>MuTightTrackerLayercut)));
        }

    }

    return tightid;
}

std::vector<unsigned int> H4LTools::goodFsrPhotons(){
    std::vector<unsigned int> goodFsrPhoton;
    for (unsigned int i=0; i<FsrPhoton_pt.size(); i++){
        if((FsrPhoton_pt[i]>fsrphotonPtcut)&&(fabs(FsrPhoton_eta[i])<fsrphotonEtacut)&&(FsrPhoton_relIso03[i]<fsrphotonIsocut)){
            goodFsrPhoton.push_back(i);
        }
    }
    return goodFsrPhoton;
}

std::vector<unsigned int> H4LTools::SelectedJets(std::vector<unsigned int> ele, std::vector<unsigned int> mu){
    std::vector<unsigned int> goodJets;
    for(unsigned int i=0;i<Jet_pt.size();i++){
        if((Jet_pt[i]>JetPtcut)&&(fabs(Jet_eta[i])<JetEtacut)){
            if((Jet_jetId[i]>0)&&(Jet_puId[i]==7)){
                int overlaptag=0;
                TLorentzVector jettest;
                jettest.SetPtEtaPhiM(Jet_pt[i],Jet_eta[i],Jet_phi[i],Jet_mass[i]);
                for(unsigned int ie=0;ie<ele.size();ie++){
                    TLorentzVector eletest;
                    eletest.SetPtEtaPhiM(Electron_pt[ele[ie]],Electron_eta[ele[ie]],Electron_phi[ele[ie]],Electron_mass[ele[ie]]);
                    if(eletest.DeltaR(jettest)<0.4) overlaptag++;
                }
                for(unsigned int im=0;im<mu.size();im++){
                    TLorentzVector mutest;
                    mutest.SetPtEtaPhiM(Muon_pt[mu[im]],Muon_eta[mu[im]],Muon_phi[mu[im]],Muon_mass[mu[im]]);
                    if(mutest.DeltaR(jettest)<0.4) overlaptag++;
                }
                if(overlaptag==0) goodJets.push_back(i);
            }
        }
    } 
    njets_pt30_eta4p7 = goodJets.size();
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
        fsrcand.SetPtEtaPhiM(FsrPhoton_pt[BestFsrPhotons[i]],FsrPhoton_eta[BestFsrPhotons[i]],FsrPhoton_phi[BestFsrPhotons[i]],0);
        float dRlC,dRlOverPtC;
        dRlC = fsrcand.DeltaR(Lep);
        if ((dRlC<fsrphotondRlcut)&&(FsrPhoton_dROverEt2[BestFsrPhotons[i]]<fsrphotondRlOverPtcut)){
            if(FsrPhoton_dROverEt2[BestFsrPhotons[i]]<dRlOverPt){
                dRl = dRlC;
                dRlOverPt = FsrPhoton_dROverEt2[BestFsrPhotons[i]];
                FsrIdx = BestFsrPhotons[i];
            }
        }
    }

    return FsrIdx;
    
}
unsigned H4LTools::doFsrRecovery_Run3(std::vector<unsigned int> goodfsridx, unsigned lepidx, int lepflavor){//lepflavor 11 or 13
    
    unsigned matchedfsridx = 999;
    if(lepflavor == 11){
        for(unsigned fsridx=0; fsridx<goodfsridx.size(); fsridx++){
            if(FsrPhoton_electronIdx[goodfsridx[fsridx]] == lepidx){
                matchedfsridx = fsridx;
                break;
            }
        }
    }
    if(lepflavor == 13){
        for(unsigned fsridx=0; fsridx<goodfsridx.size(); fsridx++){
            if(FsrPhoton_muonIdx[goodfsridx[fsridx]] == lepidx){
                matchedfsridx = fsridx;
                break;
            }
        }
    }
    return matchedfsridx;
}
void H4LTools::BatchFsrRecovery_Run3(){
    unsigned fsridx;
    std::vector<unsigned> fsrlist;
    fsrlist = goodFsrPhotons();
    for(unsigned int i=0; i<Electron_pt.size(); i++){
        TLorentzVector fsr,lep;
        lep.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],Electron_mass[i]);
        fsridx = doFsrRecovery_Run3(fsrlist,i,11);
        if(fsridx<900){
            fsr.SetPtEtaPhiM(FsrPhoton_pt[fsrlist[fsridx]], FsrPhoton_eta[fsrlist[fsridx]], FsrPhoton_phi[fsrlist[fsridx]], 0);
            lep = lep + fsr;
            Electrondressed_Run3.push_back(lep);
        }
        else{
            Electrondressed_Run3.push_back(lep);
        }
    }
    for(unsigned int j=0; j<Muon_pt.size(); j++){
        TLorentzVector fsr,lep;
        lep.SetPtEtaPhiM(Muon_pt[j],Muon_eta[j],Muon_phi[j],Muon_mass[j]);
        fsridx = doFsrRecovery_Run3(fsrlist,j,13);
        if(fsridx<900){
            fsr.SetPtEtaPhiM(FsrPhoton_pt[fsrlist[fsridx]], FsrPhoton_eta[fsrlist[fsridx]], FsrPhoton_phi[fsrlist[fsridx]], 0);
            lep = lep + fsr;
            Muondressed_Run3.push_back(lep);
        }
        else{
            Muondressed_Run3.push_back(lep);
        }
    }
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
        FsrPhoton.SetPtEtaPhiM(FsrPhoton_pt[Fsrtag],FsrPhoton_eta[Fsrtag],FsrPhoton_phi[Fsrtag],0);
        TLorentzVector LepFsrRecovery;
        LepFsrRecovery = FsrPhoton + LepList[i];
        LepFsrList.push_back(LepFsrRecovery);
    }
    return LepFsrList;
}

std::vector<float> H4LTools::ElectronFsrPt(){
    std::vector<float> lepPt;
    for (unsigned int i=0;i<Electrondressed_Run3.size();i++){
        lepPt.push_back(Electrondressed_Run3[i].Pt());
    }
    return lepPt;
}

std::vector<float> H4LTools::ElectronFsrEta(){
    std::vector<float> lepEta;
    for (unsigned int i=0;i<Electrondressed_Run3.size();i++){
        lepEta.push_back(Electrondressed_Run3[i].Eta());
    }
    return lepEta;
}

std::vector<float> H4LTools::ElectronFsrPhi(){
    std::vector<float> lepPhi;
    for (unsigned int i=0;i<Electrondressed_Run3.size();i++){
        lepPhi.push_back(Electrondressed_Run3[i].Phi());
    }
    return lepPhi;
}

std::vector<float> H4LTools::MuonFsrPt(){
    std::vector<float> lepPt;
    for (unsigned int i=0;i<Muondressed_Run3.size();i++){
        lepPt.push_back(Muondressed_Run3[i].Pt());
    }
    return lepPt;
}

std::vector<float> H4LTools::MuonFsrEta(){
    std::vector<float> lepEta;
    for (unsigned int i=0;i<Muondressed_Run3.size();i++){
        lepEta.push_back(Muondressed_Run3[i].Eta());
    }
    return lepEta;
}

std::vector<float> H4LTools::MuonFsrPhi(){
    std::vector<float> lepPhi;
    for (unsigned int i=0;i<Muondressed_Run3.size();i++){
        lepPhi.push_back(Muondressed_Run3[i].Phi());
    }
    return lepPhi;
}

/*std::vector<TLorentzVector> H4LTools::ElectronFsr(){
    std::vector<TLorentzVector> leplist;
    std::vector<TLorentzVector> leplistfsr;
    for(unsigned int i=0;i<Electron_pt.size();i++){
        TLorentzVector Lep;
        Lep.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],Electron_mass[i]);
        leplist.push_back(Lep);
    }
    leplistfsr = BatchFsrRecovery(leplist);
    return leplistfsr;
}

std::vector<TLorentzVector> H4LTools::MuonFsr(){
    std::vector<TLorentzVector> leplist;
    std::vector<TLorentzVector> leplistfsr;
    for(unsigned int i=0;i<Muon_pt.size();i++){
        TLorentzVector Lep;
        Lep.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],Muon_mass[i]);
        leplist.push_back(Lep);
    }
    leplistfsr = BatchFsrRecovery(leplist);
    return leplistfsr;
}
std::vector<float> H4LTools::ElectronFsrPt(){
    std::vector<float> lepPt;
    std::vector<TLorentzVector> leplistfsr;
    leplistfsr = ElectronFsr();
    for (unsigned int i=0;i<leplistfsr.size();i++){
        lepPt.push_back(leplistfsr[i].Pt());
    }
    return lepPt;
}

std::vector<float> H4LTools::ElectronFsrEta(){
    std::vector<float> lepEta;
    std::vector<TLorentzVector> leplistfsr;
    leplistfsr = ElectronFsr();
    for (unsigned int i=0;i<leplistfsr.size();i++){
        lepEta.push_back(leplistfsr[i].Eta());
    }
    return lepEta;
}

std::vector<float> H4LTools::ElectronFsrPhi(){
    std::vector<float> lepPhi;
    std::vector<TLorentzVector> leplistfsr;
    leplistfsr = ElectronFsr();
    for (unsigned int i=0;i<leplistfsr.size();i++){
        lepPhi.push_back(leplistfsr[i].Phi());
    }
    return lepPhi;
}

std::vector<float> H4LTools::MuonFsrPt(){
    std::vector<float> lepPt;
    std::vector<TLorentzVector> leplistfsr;
    leplistfsr = MuonFsr();
    for (unsigned int i=0;i<leplistfsr.size();i++){
        lepPt.push_back(leplistfsr[i].Pt());
    }
    return lepPt;
}

std::vector<float> H4LTools::MuonFsrEta(){
    std::vector<float> lepEta;
    std::vector<TLorentzVector> leplistfsr;
    leplistfsr = MuonFsr();
    for (unsigned int i=0;i<leplistfsr.size();i++){
        lepEta.push_back(leplistfsr[i].Eta());
    }
    return lepEta;
}

std::vector<float> H4LTools::MuonFsrPhi(){
    std::vector<float> lepPhi;
    std::vector<TLorentzVector> leplistfsr;
    leplistfsr = MuonFsr();
    for (unsigned int i=0;i<leplistfsr.size();i++){
        lepPhi.push_back(leplistfsr[i].Phi());
    }
    return lepPhi;
}*/

void H4LTools::LeptonSelection(){
    
    looseEle = goodLooseElectrons2012();
    looseMu = goodLooseMuons2012();
    bestEle = goodElectrons2015_noIso_noBdt(looseEle);
    bestMu = goodMuons2015_noIso_noPf(looseMu);
    Electronindex = bestEle;
    Muonindex = bestMu;
    AllEid = passTight_BDT_Id();
    AllMuid = passTight_Id();
    for (unsigned int iuj=0;iuj<looseEle.size();iuj++){
        if(AllEid[looseEle[iuj]]) tighteleforjetidx.push_back(looseEle[iuj]);
    }
    for (unsigned int juj=0;juj<looseMu.size();juj++){
        if(AllMuid[looseMu[juj]]) tightmuforjetidx.push_back(looseMu[juj]);
    }
    jetidx = SelectedJets(tighteleforjetidx,tightmuforjetidx);
    
    for(unsigned int ie=0; ie<Electronindex.size();ie++){
        if(Electron_pdgId[Electronindex[ie]]>0){
            Elechg.push_back(-1);
        }
        else{
            Elechg.push_back(1);
        }
        TLorentzVector Ele;
        Ele.SetPtEtaPhiM(Electron_pt[Electronindex[ie]],Electron_eta[Electronindex[ie]],Electron_phi[Electronindex[ie]],Electron_mass[Electronindex[ie]]);
        Elelist.push_back(Ele);
        ElelistFsr.push_back(Electrondressed_Run3[Electronindex[ie]]);
        Eiso.push_back(Electron_pfRelIso03_all[Electronindex[ie]]);
        Eid.push_back(AllEid[Electronindex[ie]]);
    }

    for(unsigned int imu=0; imu<Muonindex.size();imu++){
        if(Muon_pdgId[Muonindex[imu]]>0){
            Muchg.push_back(-1);
        }
        else{
            Muchg.push_back(1);
        }
        TLorentzVector Mu;
        Mu.SetPtEtaPhiM(Muon_pt[Muonindex[imu]],Muon_eta[Muonindex[imu]],Muon_phi[Muonindex[imu]],Muon_mass[Muonindex[imu]]);
        Mulist.push_back(Mu);
        MulistFsr.push_back(Muondressed_Run3[Muonindex[imu]]);
        muid.push_back(AllMuid[Muonindex[imu]]);
        Muiso.push_back(Muon_pfRelIso03_all[Muonindex[imu]]);
    }
    
    //ElelistFsr = BatchFsrRecovery(Elelist); 
    //MulistFsr = BatchFsrRecovery(Mulist);
    
    for(unsigned int ae=0; ae<Eid.size();ae++){
        float RelEleIsoNoFsr;
        RelEleIsoNoFsr = Eiso[ae];
        unsigned FsrEleidx;
        FsrEleidx = doFsrRecovery_Run3(goodFsrPhotons(), Electronindex[ae], 11);
        Elelep_index.push_back(Lepointer);
        Lepointer++;
        lep_pt.push_back(Elelist[ae].Pt());
        lep_eta.push_back(Elelist[ae].Eta());
        lep_phi.push_back(Elelist[ae].Phi());
        lep_mass.push_back(Elelist[ae].M());
        if (Electron_charge[Electronindex[ae]]<0) lep_id.push_back(11);
        else lep_id.push_back(-11);
        if (isFSR && (FsrEleidx < 900)){
            TLorentzVector fsrele;
            fsrele.SetPtEtaPhiM(FsrPhoton_pt[FsrEleidx],FsrPhoton_eta[FsrEleidx],FsrPhoton_phi[FsrEleidx],0);
            if(Elelist[ae].DeltaR(fsrele)>0.01){
                RelEleIsoNoFsr = RelEleIsoNoFsr - FsrPhoton_pt[FsrEleidx]/Elelist[ae].Pt(); 
            }
          //FsrEleidx = doFsrRecovery(Elelist[ae]);
          /*if(FsrEleidx<900){
              TLorentzVector fsrele;
              fsrele.SetPtEtaPhiM(FsrPhoton_pt[FsrEleidx],FsrPhoton_eta[FsrEleidx],FsrPhoton_phi[FsrEleidx],0);
              std::cout<<"Ele correction: "<< std::endl;
              if(Elelist[ae].DeltaR(fsrele)>0.01){
                RelEleIsoNoFsr = RelEleIsoNoFsr - FsrPhoton_pt[FsrEleidx]/Elelist[ae].Pt(); 
              }
          }*/
        }
        if((Eid[ae]==true)&&(RelEleIsoNoFsr<9999)){
            nTightEle++;
            TightEleindex.push_back(ae);
            nTightEleChgSum += Elechg[ae];
            istightele.push_back(true);
            lep_tightId.push_back(1);
        }
        else{
            istightele.push_back(false);
            lep_tightId.push_back(0);
        }
        if (isMC){
            lep_genindex.push_back(Electron_genPartIdx[Electronindex[ae]]);
            lep_matchedR03_PdgId.push_back(GenPart_pdgId[Electron_genPartIdx[Electronindex[ae]]]);
            lep_matchedR03_MomId.push_back(motherID(Electron_genPartIdx[Electronindex[ae]]));
            lep_matchedR03_MomMomId.push_back(motherID(motheridx(Electron_genPartIdx[Electronindex[ae]])));
        } 
        else {
            lep_genindex.push_back(-1);
            lep_matchedR03_PdgId.push_back(-1);
            lep_matchedR03_MomId.push_back(-1);
            lep_matchedR03_MomMomId.push_back(-1);
        }
    }

    for(unsigned int amu=0; amu<muid.size();amu++){
        float RelIsoNoFsr;
        RelIsoNoFsr = Muiso[amu];
        unsigned int FsrMuonidx;
        FsrMuonidx = doFsrRecovery_Run3(goodFsrPhotons(), Muonindex[amu], 13);
        Mulep_index.push_back(Lepointer);
        Lepointer++;
        lep_pt.push_back(Mulist[amu].Pt());
        lep_eta.push_back(Mulist[amu].Eta());
        lep_phi.push_back(Mulist[amu].Phi());
        lep_mass.push_back(Mulist[amu].M());
        if (Muon_charge[Muonindex[amu]]<0) lep_id.push_back(13);
        else lep_id.push_back(-13);
        if (isFSR && (FsrMuonidx < 900)){
            TLorentzVector fsrmuon;
            fsrmuon.SetPtEtaPhiM(FsrPhoton_pt[FsrMuonidx],FsrPhoton_eta[FsrMuonidx],FsrPhoton_phi[FsrMuonidx],0);
            if(Mulist[amu].DeltaR(fsrmuon)>0.01){
                RelIsoNoFsr = RelIsoNoFsr - FsrPhoton_pt[FsrMuonidx]/Mulist[amu].Pt();  
            }              
          /*if(FsrMuonidx<900){
              TLorentzVector fsrmuon;
              fsrmuon.SetPtEtaPhiM(FsrPhoton_pt[FsrMuonidx],FsrPhoton_eta[FsrMuonidx],FsrPhoton_phi[FsrMuonidx],0);
              std::cout<<"muon FSR recovered"<<endl;
              if(Mulist[amu].DeltaR(fsrmuon)>0.01){
                RelIsoNoFsr = RelIsoNoFsr - FsrPhoton_pt[FsrMuonidx]/Mulist[amu].Pt();  
              }
          }*/
        }
        if((muid[amu]==true)&&(RelIsoNoFsr<0.35)){
            nTightMu++;
            TightMuindex.push_back(amu);
            nTightMuChgSum += Muchg[amu];
            istightmu.push_back(true);
            lep_tightId.push_back(1);
        }
        else{
            istightmu.push_back(false);
            lep_tightId.push_back(0);
        }
        if (isMC) {
            lep_genindex.push_back(Muon_genPartIdx[Muonindex[amu]]);
            lep_matchedR03_PdgId.push_back(GenPart_pdgId[Muon_genPartIdx[Muonindex[amu]]]);
            lep_matchedR03_MomId.push_back(motherID(Muon_genPartIdx[Muonindex[amu]]));
            lep_matchedR03_MomMomId.push_back(motherID(motheridx(Muon_genPartIdx[Muonindex[amu]])));
        }
        else {
            lep_genindex.push_back(-1);
            lep_matchedR03_PdgId.push_back(-1);
            lep_matchedR03_MomId.push_back(-1);
            lep_matchedR03_MomMomId.push_back(-1);
        }
        
    }
    
    
}
bool H4LTools::findZCandidate(){
    
    TLorentzVector z1,z2;
    
    if (nTightEle>=4) {
        cut4e++;
        flag4e = true;
    }
    else if (nTightMu>=4){
        cut4mu++;
        flag4mu = true;
    } 
    else if ((nTightMu>=2)&&(nTightEle>=2)){
        cut2e2mu++;
        flag2e2mu = true;
    } 
    
    
    if(Electronindex.size()>1){
        for(unsigned int ke=0; ke<(Electronindex.size()-1);ke++){
            for(unsigned int je=ke+1;je<Electronindex.size();je++){
                if ((Elechg[Electronindex[ke]]+Elechg[Electronindex[je]])==0){
                    TLorentzVector Zcan;
                    Zcan = ElelistFsr[Electronindex[ke]] + ElelistFsr[Electronindex[je]];
                    if((Zcan.M()>MZcutdown)&&(Zcan.M()<MZcutup)){
                        Zlist.push_back(Zcan);
                        Zlep1index.push_back(Electronindex[ke]);
                        Zlep2index.push_back(Electronindex[je]);
                        Zlep1lepindex.push_back(Elelep_index[ke]);
                        Zlep2lepindex.push_back(Elelep_index[je]);
                        Zflavor.push_back(11);
                        Zlep1pt.push_back(ElelistFsr[Electronindex[ke]].Pt());
                        Zlep2pt.push_back(ElelistFsr[Electronindex[je]].Pt());
                        Zlep1eta.push_back(ElelistFsr[Electronindex[ke]].Eta());
                        Zlep2eta.push_back(ElelistFsr[Electronindex[je]].Eta());
                        Zlep1phi.push_back(ElelistFsr[Electronindex[ke]].Phi());
                        Zlep2phi.push_back(ElelistFsr[Electronindex[je]].Phi());
                        Zlep1mass.push_back(ElelistFsr[Electronindex[ke]].M());
                        Zlep2mass.push_back(ElelistFsr[Electronindex[je]].M());
                        Zlep1ptNoFsr.push_back(Elelist[Electronindex[ke]].Pt());
                        Zlep2ptNoFsr.push_back(Elelist[Electronindex[je]].Pt());
                        Zlep1etaNoFsr.push_back(Elelist[Electronindex[ke]].Eta());
                        Zlep2etaNoFsr.push_back(Elelist[Electronindex[je]].Eta());
                        Zlep1phiNoFsr.push_back(Elelist[Electronindex[ke]].Phi());
                        Zlep2phiNoFsr.push_back(Elelist[Electronindex[je]].Phi());
                        Zlep1massNoFsr.push_back(Elelist[Electronindex[ke]].M());
                        Zlep2massNoFsr.push_back(Elelist[Electronindex[je]].M());
                        Zlep1chg.push_back(Elechg[Electronindex[ke]]);
                        Zlep2chg.push_back(Elechg[Electronindex[je]]);
                        Zlep1istight.push_back(istightele[ke]);
                        Zlep2istight.push_back(istightele[je]);
                    }
                }
            }
        }
    }
    

    if(Muonindex.size()>1){
        for(unsigned int kmu=0; kmu<(Muonindex.size()-1);kmu++){
            for(unsigned int jmu=kmu+1;jmu<Muonindex.size();jmu++){
                if ((Muchg[Muonindex[kmu]]+Muchg[Muonindex[jmu]])==0){
                    TLorentzVector Zcan;
                    Zcan = MulistFsr[Muonindex[kmu]] + MulistFsr[Muonindex[jmu]];
                    if((Zcan.M()>MZcutdown)&&(Zcan.M()<MZcutup)){
                        Zlist.push_back(Zcan);
                        Zlep1index.push_back(Muonindex[kmu]);
                        Zlep2index.push_back(Muonindex[jmu]);
                        Zlep1lepindex.push_back(Mulep_index[kmu]);
                        Zlep2lepindex.push_back(Mulep_index[jmu]);
                        Zflavor.push_back(13);
                        Zlep1pt.push_back(MulistFsr[Muonindex[kmu]].Pt());
                        Zlep2pt.push_back(MulistFsr[Muonindex[jmu]].Pt());
                        Zlep1eta.push_back(MulistFsr[Muonindex[kmu]].Eta());
                        Zlep2eta.push_back(MulistFsr[Muonindex[jmu]].Eta());
                        Zlep1phi.push_back(MulistFsr[Muonindex[kmu]].Phi());
                        Zlep2phi.push_back(MulistFsr[Muonindex[jmu]].Phi());
                        Zlep1mass.push_back(MulistFsr[Muonindex[kmu]].M());
                        Zlep2mass.push_back(MulistFsr[Muonindex[jmu]].M());
                        Zlep1ptNoFsr.push_back(Mulist[Muonindex[kmu]].Pt());
                        Zlep2ptNoFsr.push_back(Mulist[Muonindex[jmu]].Pt());
                        Zlep1etaNoFsr.push_back(Mulist[Muonindex[kmu]].Eta());
                        Zlep2etaNoFsr.push_back(Mulist[Muonindex[jmu]].Eta());
                        Zlep1phiNoFsr.push_back(Mulist[Muonindex[kmu]].Phi());
                        Zlep2phiNoFsr.push_back(Mulist[Muonindex[jmu]].Phi());
                        Zlep1massNoFsr.push_back(Mulist[Muonindex[kmu]].M());
                        Zlep2massNoFsr.push_back(Mulist[Muonindex[jmu]].M());
                        Zlep1chg.push_back(Muchg[Muonindex[kmu]]);
                        Zlep2chg.push_back(Muchg[Muonindex[jmu]]);
                        Zlep1istight.push_back(istightmu[kmu]);
                        Zlep2istight.push_back(istightmu[jmu]);      
                    }
                }
            }
        }
    }
    Zsize = Zlist.size(); 
    for (unsigned int znofsr = 0; znofsr<Zsize; znofsr++){
        TLorentzVector Zlep1nofsr,Zlep2nofsr,Zcannofsr;
        Zlep1nofsr.SetPtEtaPhiM(Zlep1ptNoFsr[znofsr],Zlep1etaNoFsr[znofsr],Zlep1phiNoFsr[znofsr],Zlep1massNoFsr[znofsr]);
        Zlep2nofsr.SetPtEtaPhiM(Zlep2ptNoFsr[znofsr],Zlep2etaNoFsr[znofsr],Zlep2phiNoFsr[znofsr],Zlep2massNoFsr[znofsr]);
        Zcannofsr = Zlep1nofsr + Zlep2nofsr;
        Zlistnofsr.push_back(Zcannofsr);
    }
    
    for (unsigned int tightZidx = 0; tightZidx<Zsize; tightZidx++){
        if (Zlep1istight[tightZidx]&&Zlep2istight[tightZidx]){
            Zistight.push_back(true);
            nTightZ++; 
        } 
        else{
            Zistight.push_back(false);
       
        }
    }
    if (Zsize>0){
        return true;
    }
    else{
        return false;
    }
    


}

void H4LTools::findZ1LCandidate(){
    if(!(Electronindex.size()+Muonindex.size()==3)) return;
    if(Zsize==0) return;
    if(nTightZ<1) return;
    double minZ1DeltaM=9999.9;
    for (unsigned int i=0;i<Zsize;i++){
        TLorentzVector lep_i1, lep_i2, lep_j1;
        TLorentzVector lep_i1_nofsr, lep_i2_nofsr, lep_j1_nofsr;
        float lep_i1_chg = 0, lep_i2_chg = 0, lep_j1_chg = 0;
        lep_i1.SetPtEtaPhiM(Zlep1pt[i],Zlep1eta[i], Zlep1phi[i], Zlep1mass[i]);
        lep_i2.SetPtEtaPhiM(Zlep2pt[i],Zlep2eta[i], Zlep2phi[i], Zlep2mass[i]);
        lep_i1_nofsr.SetPtEtaPhiM(Zlep1ptNoFsr[i],Zlep1etaNoFsr[i], Zlep1phiNoFsr[i], Zlep1mass[i]);
        lep_i2_nofsr.SetPtEtaPhiM(Zlep2ptNoFsr[i],Zlep2etaNoFsr[i], Zlep2phiNoFsr[i], Zlep2mass[i]);
        lep_i1_chg = Zlep1chg[i]; lep_i2_chg = Zlep2chg[i];
        bool foundj1 = false;
        for (unsigned int findj1e=0;findj1e<Electronindex.size();findj1e++){
            if (abs(ElelistFsr[Electronindex[findj1e]].Eta()-lep_i1.Eta())<0.001) continue;
            if (abs(ElelistFsr[Electronindex[findj1e]].Phi()-lep_i1.Phi())<0.001) continue;
            if (abs(ElelistFsr[Electronindex[findj1e]].Eta()-lep_i2.Eta())<0.001) continue;
            if (abs(ElelistFsr[Electronindex[findj1e]].Phi()-lep_i2.Phi())<0.001) continue;
            foundj1 = true;
            lep_j1_chg = Elechg[Electronindex[findj1e]];
            lep_j1.SetPtEtaPhiM(ElelistFsr[Electronindex[findj1e]].Pt(),ElelistFsr[Electronindex[findj1e]].Eta(), ElelistFsr[Electronindex[findj1e]].Phi(), ElelistFsr[Electronindex[findj1e]].M());
            lep_j1_nofsr.SetPtEtaPhiM(Elelist[Electronindex[findj1e]].Pt(),Elelist[Electronindex[findj1e]].Eta(), Elelist[Electronindex[findj1e]].Phi(), Elelist[Electronindex[findj1e]].M());
        }
        for (unsigned int findj1mu=0;findj1mu<Muonindex.size();findj1mu++){
            if (abs(MulistFsr[Muonindex[findj1mu]].Eta()-lep_i1.Eta())<0.001) continue;
            if (abs(MulistFsr[Muonindex[findj1mu]].Phi()-lep_i1.Phi())<0.001) continue;
            if (abs(MulistFsr[Muonindex[findj1mu]].Eta()-lep_i2.Eta())<0.001) continue;
            if (abs(MulistFsr[Muonindex[findj1mu]].Phi()-lep_i2.Phi())<0.001) continue;
            foundj1 = true;
            lep_j1_chg = Muchg[Muonindex[findj1mu]];
            lep_j1.SetPtEtaPhiM(MulistFsr[Muonindex[findj1mu]].Pt(),MulistFsr[Muonindex[findj1mu]].Eta(),MulistFsr[Muonindex[findj1mu]].Phi(),MulistFsr[Muonindex[findj1mu]].M());
            lep_j1_nofsr.SetPtEtaPhiM(Mulist[Muonindex[findj1mu]].Pt(),Mulist[Muonindex[findj1mu]].Eta(),Mulist[Muonindex[findj1mu]].Phi(),Mulist[Muonindex[findj1mu]].M());
        }
        if(!foundj1) continue;
        //Check PtCut
        if((lep_i1.Pt()<20)&&(lep_i2.Pt()<20)) continue;
        if((lep_i1.Pt()<10) || (lep_i2.Pt()<10)) continue;
        //check dRCut
        if(lep_i1.DeltaR(lep_i2)<0.02) continue;
        if(lep_j1.DeltaR(lep_i2)<0.02) continue;
        if(lep_i1.DeltaR(lep_j1)<0.02) continue;

        // Check M(l+,l-)>4.0 GeV for any OS pair
        // Do not include FSR photons
        if ((lep_i1_chg*lep_j1_chg) < 0){
            TLorentzVector lepi1j1nofsr;
            lepi1j1nofsr = lep_i1_nofsr + lep_j1_nofsr;
            if(lepi1j1nofsr.M()<4.0) continue;
        }

        if ((lep_i2_chg*lep_j1_chg) < 0){
            TLorentzVector lepi2j1nofsr;
            lepi2j1nofsr = lep_i2_nofsr + lep_j1_nofsr;
            if(lepi2j1nofsr.M()<4.0) continue;
        }
        
        if(!Zistight[i]) continue;
        if((Zlist[i].M()<40)||(Zlist[i].M()>120)) continue;
        double Z1DeltaM = abs(Zlist[i].M()-Zmass);
        if ( Z1DeltaM<=minZ1DeltaM ) {

            minZ1DeltaM = Z1DeltaM;
            TLorentzVector Z1L;
            Z1L = Zlist[i]+lep_j1;
            mass3l = Z1L.M();
            Z1LZ1index = i;
            passedZ1LSelection = true;
        }
    }
    return;
}

bool H4LTools::ZZSelection(){
    
    bool foundZZCandidate = false;
    bool proper_Leps = false;
    int Nmm,Nmp,Nep,Nem;
    Nmm=0;Nmp=0;Nep=0;Nem=0;
    for (unsigned int iec=0;iec<Electronindex.size();iec++){
        if(Elechg[Electronindex[iec]]>0) Nep++;
        else Nem++;
    }
    for (unsigned int imc=0;imc<Muonindex.size();imc++){
        if(Muchg[Muonindex[imc]]>0) Nmp++;
        else Nmm++;
    }
    if (Nmm>=2 && Nmp>=2) proper_Leps=true;
    if (Nem>=2 && Nep>=2) proper_Leps=true;
    if (Nem>0 && Nep>0 && Nmm>0 && Nmp>0) proper_Leps=true;
    if (!proper_Leps) return foundZZCandidate;
    if(Zlist.size()==0) return foundZZCandidate;
    if(Zsize<2) return foundZZCandidate;

    bool foundSRCandidate=false;    
    //Find ZZ candidate
    std::vector<int> Z1CanIndex;
    std::vector<int> Z2CanIndex;
    int ghosttag = 0, QCDtag=0, lepPtTag = 0;
    for (unsigned int m=0; m<(Zlist.size()-1); m++){
        for (unsigned int n=m+1; n<Zlist.size(); n++){
            if(!(Zistight[m] || Zistight[n])) continue;
            if (Zflavor[m]==Zflavor[n]){
               if ((Zlep1index[m] == Zlep1index[n])||(Zlep2index[m] == Zlep1index[n])) continue;  //non-overlapping
               if ((Zlep1index[m] == Zlep2index[n])||(Zlep2index[m] == Zlep2index[n])) continue;
            }
            if (Zlist[m].DeltaR(Zlist[n])<0.02) continue; //ghost removal
            ghosttag++;
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
            lepPtTag++;
            if ((Zlep1chg[m]+Zlep1chg[n])==0){
                TLorentzVector lepA,lepB,lepAB;
                lepA.SetPtEtaPhiM(Zlep1ptNoFsr[m],Zlep1etaNoFsr[m],Zlep1phiNoFsr[m],Zlep1massNoFsr[m]);
                lepB.SetPtEtaPhiM(Zlep2ptNoFsr[n],Zlep2etaNoFsr[n],Zlep2phiNoFsr[n],Zlep2massNoFsr[n]);
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
            QCDtag++;
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
            if (zZ1.M()+zZ2.M()<MZZcut) continue;
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
    if(ghosttag){
        if (flag2e2mu) cutghost2e2mu++;
        if (flag4e) cutghost4e++;
        if (flag4mu) cutghost4mu++;
    }
    if(lepPtTag){
        if (flag2e2mu) cutLepPt2e2mu++;
        if (flag4e) cutLepPt4e++;
        if (flag4mu) cutLepPt4mu++;
    }
    if(QCDtag){
        if (flag2e2mu) cutQCD2e2mu++;
        if (flag4e) cutQCD4e++;
        if (flag4mu) cutQCD4mu++;
    }
    if(foundZZCandidate == false){
        return foundZZCandidate;
    }
    if (flag2e2mu) cutZZ2e2mu++;
    if (flag4e) cutZZ4e++;
    if (flag4mu) cutZZ4mu++;
    int Z1index,Z2index; 
    Z1index = Z1CanIndex[0];
    Z2index = Z2CanIndex[0];
    float Z2PtsumSR,Z2PtsumCR;
    double minZ1DeltaM_CR=99999.9;double minZ1DeltaM_SR=99999.9;double max_D_bkg_kin_CR=0.0;double max_D_bkg_kin_SR=0.0;
    Z2PtsumSR = -1;Z2PtsumCR = -1;
    if(Z1CanIndex.size()>1){
        for(unsigned int iz=0;iz<Z1CanIndex.size();iz++){
            if(isFiducialAna){
                if (Zistight[Z1CanIndex[iz]] && Zistight[Z2CanIndex[iz]]) foundSRCandidate=true; //it is a global switch.
                if (foundSRCandidate){ //SR has the priority
                    if (!(Zistight[Z1CanIndex[iz]] && Zistight[Z2CanIndex[iz]])) continue;
                    if (Z1index==Z1CanIndex[iz]){
                        if((Zlep1pt[Z2CanIndex[iz]] + Zlep2pt[Z2CanIndex[iz]])>Z2PtsumSR){
                            Z1index = Z1CanIndex[iz];
                            Z2index = Z2CanIndex[iz];
                            Z2PtsumSR = Zlep1pt[Z2index] + Zlep2pt[Z2index];
                        }
                    }
                    if(fabs(Zlist[Z1CanIndex[iz]].M()-Zmass)<minZ1DeltaM_SR){
                        Z1index = Z1CanIndex[iz];
                        Z2index = Z2CanIndex[iz];
                        Z2PtsumSR = Zlep1pt[Z2index] + Zlep2pt[Z2index];
                        minZ1DeltaM_SR = fabs(Zlist[Z1CanIndex[iz]].M()-Zmass);
                    }    
                    passedFullSelection = true;
		    passedZXCRSelection = false;
                }
                else {
                    if (Z1index==Z1CanIndex[iz]){
                        if((Zlep1pt[Z2CanIndex[iz]] + Zlep2pt[Z2CanIndex[iz]])>Z2PtsumCR){
                            Z1index = Z1CanIndex[iz];
                            Z2index = Z2CanIndex[iz];
                            Z2PtsumCR = Zlep1pt[Z2index] + Zlep2pt[Z2index];
                        }
                    }
                    if(fabs(Zlist[Z1CanIndex[iz]].M()-Zmass)<minZ1DeltaM_CR){
                        Z1index = Z1CanIndex[iz];
                        Z2index = Z2CanIndex[iz];
                        Z2PtsumCR = Zlep1pt[Z2index] + Zlep2pt[Z2index];
                        minZ1DeltaM_CR = fabs(Zlist[Z1CanIndex[iz]].M()-Zmass);
                    }    
                    passedZXCRSelection = true;
                }
                
            }
            else{
                if (Zistight[Z1CanIndex[iz]] && Zistight[Z2CanIndex[iz]]) foundSRCandidate=true; //it is a global switch.
                bool same4l = false;
                bool foundZ11=false; bool foundZ12=false; bool foundZ21=false; bool foundZ22=false;
                vector<int> lep_Hindex,lep_Hflavor;
                lep_Hindex.clear();lep_Hflavor.clear();;
                lep_Hindex.push_back(Zlep1index[Z1index]);
                lep_Hindex.push_back(Zlep2index[Z1index]);
                lep_Hindex.push_back(Zlep1index[Z2index]);
                lep_Hindex.push_back(Zlep2index[Z2index]);
                lep_Hflavor.push_back(Zflavor[Z1index]);
                lep_Hflavor.push_back(Zflavor[Z1index]);
                lep_Hflavor.push_back(Zflavor[Z2index]);
                lep_Hflavor.push_back(Zflavor[Z2index]);
                for(int l = 0; l < 4; l++){
                    if ((lep_Hindex[l]==Zlep1index[Z1CanIndex[iz]])&&(lep_Hflavor[l]==Zflavor[Z1CanIndex[iz]])) foundZ11 = true;
                    if ((lep_Hindex[l]==Zlep2index[Z1CanIndex[iz]])&&(lep_Hflavor[l]==Zflavor[Z1CanIndex[iz]])) foundZ12 = true;
                    if ((lep_Hindex[l]==Zlep1index[Z2CanIndex[iz]])&&(lep_Hflavor[l]==Zflavor[Z2CanIndex[iz]])) foundZ21 = true;
                    if ((lep_Hindex[l]==Zlep2index[Z2CanIndex[iz]])&&(lep_Hflavor[l]==Zflavor[Z2CanIndex[iz]])) foundZ22 = true;
                }
                same4l = (foundZ11 && foundZ12 && foundZ21 && foundZ22);
                if ((!same4l) && foundSRCandidate){ //SR has the priority
                    SimpleParticleCollection_t daughtersSR;
                    SimpleParticleCollection_t associatedSR;
                    TLorentzVector Lep1SR,Lep2SR,Lep3SR,Lep4SR;
                    
                    Lep1SR.SetPtEtaPhiM(Zlep1pt[Z1CanIndex[iz]],Zlep1eta[Z1CanIndex[iz]],Zlep1phi[Z1CanIndex[iz]],Zlep1mass[Z1CanIndex[iz]]);
                    Lep2SR.SetPtEtaPhiM(Zlep2pt[Z1CanIndex[iz]],Zlep2eta[Z1CanIndex[iz]],Zlep2phi[Z1CanIndex[iz]],Zlep2mass[Z1CanIndex[iz]]);
                    Lep3SR.SetPtEtaPhiM(Zlep1pt[Z2CanIndex[iz]],Zlep1eta[Z2CanIndex[iz]],Zlep1phi[Z2CanIndex[iz]],Zlep1mass[Z2CanIndex[iz]]);
                    Lep4SR.SetPtEtaPhiM(Zlep2pt[Z2CanIndex[iz]],Zlep2eta[Z2CanIndex[iz]],Zlep2phi[Z2CanIndex[iz]],Zlep2mass[Z2CanIndex[iz]]);

                    daughtersSR.push_back(SimpleParticle_t((-1)*Zflavor[Z1CanIndex[iz]]*Zlep1chg[Z1CanIndex[iz]], Lep1SR));
                    daughtersSR.push_back(SimpleParticle_t((-1)*Zflavor[Z1CanIndex[iz]]*Zlep2chg[Z1CanIndex[iz]], Lep2SR));
                    daughtersSR.push_back(SimpleParticle_t((-1)*Zflavor[Z2CanIndex[iz]]*Zlep1chg[Z2CanIndex[iz]], Lep3SR));
                    daughtersSR.push_back(SimpleParticle_t((-1)*Zflavor[Z2CanIndex[iz]]*Zlep2chg[Z2CanIndex[iz]], Lep4SR));

                    float D_bkg_kin_tmp=-999; 
                    mela->setInputEvent(&daughtersSR, &associatedSR, 0, 0);
                    mela->setCurrentCandidateFromIndex(0);

                    float me_0plus_JHU_tmp, me_qqZZ_MCFM_tmp;
                    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
                    mela->computeP(me_0plus_JHU_tmp, true);            
                    mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
                    mela->computeP(me_qqZZ_MCFM_tmp, true);
                    D_bkg_kin_tmp = me_0plus_JHU_tmp / (me_0plus_JHU_tmp + me_qqZZ_MCFM_tmp);

                    mela->resetInputEvent(); 
                    if (D_bkg_kin_tmp>max_D_bkg_kin_SR){
                        max_D_bkg_kin_SR = D_bkg_kin_tmp;
                        Z1index = Z1CanIndex[iz];
                        Z2index = Z2CanIndex[iz];
                    }
                    passedFullSelection = true;
                }
                else if(same4l && (foundSRCandidate)){
                    if (Z1index==Z1CanIndex[iz]){
                        if((Zlep1pt[Z2CanIndex[iz]] + Zlep2pt[Z2CanIndex[iz]])>Z2PtsumSR){
                            Z1index = Z1CanIndex[iz];
                            Z2index = Z2CanIndex[iz];
                            Z2PtsumSR = Zlep1pt[Z2index] + Zlep2pt[Z2index];
                        }
                    }
                    if(fabs(Zlist[Z1CanIndex[iz]].M()-Zmass)<minZ1DeltaM_SR){
                        Z1index = Z1CanIndex[iz];
                        Z2index = Z2CanIndex[iz];
                        Z2PtsumSR = Zlep1pt[Z2index] + Zlep2pt[Z2index];
                        minZ1DeltaM_SR = fabs(Zlist[Z1CanIndex[iz]].M()-Zmass);
                    }    
                    passedFullSelection = true;
                }
                else if(((!same4l) && (!foundSRCandidate))){
                    SimpleParticleCollection_t daughtersCR;
                    SimpleParticleCollection_t associatedCR;
                    TLorentzVector Lep1CR,Lep2CR,Lep3CR,Lep4CR;
                    
                    Lep1CR.SetPtEtaPhiM(Zlep1pt[Z1CanIndex[iz]],Zlep1eta[Z1CanIndex[iz]],Zlep1phi[Z1CanIndex[iz]],Zlep1mass[Z1CanIndex[iz]]);
                    Lep2CR.SetPtEtaPhiM(Zlep2pt[Z1CanIndex[iz]],Zlep2eta[Z1CanIndex[iz]],Zlep2phi[Z1CanIndex[iz]],Zlep2mass[Z1CanIndex[iz]]);
                    Lep3CR.SetPtEtaPhiM(Zlep1pt[Z2CanIndex[iz]],Zlep1eta[Z2CanIndex[iz]],Zlep1phi[Z2CanIndex[iz]],Zlep1mass[Z2CanIndex[iz]]);
                    Lep4CR.SetPtEtaPhiM(Zlep2pt[Z2CanIndex[iz]],Zlep2eta[Z2CanIndex[iz]],Zlep2phi[Z2CanIndex[iz]],Zlep2mass[Z2CanIndex[iz]]);

                    daughtersCR.push_back(SimpleParticle_t((-1)*Zflavor[Z1CanIndex[iz]]*Zlep1chg[Z1CanIndex[iz]], Lep1CR));
                    daughtersCR.push_back(SimpleParticle_t((-1)*Zflavor[Z1CanIndex[iz]]*Zlep2chg[Z1CanIndex[iz]], Lep2CR));
                    daughtersCR.push_back(SimpleParticle_t((-1)*Zflavor[Z2CanIndex[iz]]*Zlep1chg[Z2CanIndex[iz]], Lep3CR));
                    daughtersCR.push_back(SimpleParticle_t((-1)*Zflavor[Z2CanIndex[iz]]*Zlep2chg[Z2CanIndex[iz]], Lep4CR));

                    float D_bkg_kin_tmp=-999; 
                    mela->setInputEvent(&daughtersCR, &associatedCR, 0, 0);
                    mela->setCurrentCandidateFromIndex(0);

                    float me_0plus_JHU_tmp, me_qqZZ_MCFM_tmp;
                    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
                    mela->computeP(me_0plus_JHU_tmp, true);            
                    mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
                    mela->computeP(me_qqZZ_MCFM_tmp, true);
                    D_bkg_kin_tmp = me_0plus_JHU_tmp / (me_0plus_JHU_tmp + me_qqZZ_MCFM_tmp);

                    mela->resetInputEvent(); 
                    if (D_bkg_kin_tmp>max_D_bkg_kin_CR){
                        max_D_bkg_kin_CR = D_bkg_kin_tmp;
                        Z1index = Z1CanIndex[iz];
                        Z2index = Z2CanIndex[iz];
                    }
                    passedZXCRSelection = true;
                }
                else if(same4l && (!foundSRCandidate)){
                    if (Z1index==Z1CanIndex[iz]){
                        if((Zlep1pt[Z2CanIndex[iz]] + Zlep2pt[Z2CanIndex[iz]])>Z2PtsumCR){
                            Z1index = Z1CanIndex[iz];
                            Z2index = Z2CanIndex[iz];
                            Z2PtsumCR = Zlep1pt[Z2index] + Zlep2pt[Z2index];
                        }
                    }
                    if(fabs(Zlist[Z1CanIndex[iz]].M()-Zmass)<minZ1DeltaM_CR){
                        Z1index = Z1CanIndex[iz];
                        Z2index = Z2CanIndex[iz];
                        Z2PtsumCR = Zlep1pt[Z2index] + Zlep2pt[Z2index];
                        minZ1DeltaM_CR = fabs(Zlist[Z1CanIndex[iz]].M()-Zmass);
                    }    
                    passedFullSelection = true;
                }

            }
                     
    	    
        }
            
    }
    
       
    
    Z1 = Zlist[Z1index];
    Z2 = Zlist[Z2index];

    Z1nofsr = Zlistnofsr[Z1index];
    Z2nofsr = Zlistnofsr[Z2index];
    if(!Zlep1istight[Z1index]) nfailedleptons++;
    if(!Zlep2istight[Z1index]) nfailedleptons++;
    if(!Zlep1istight[Z2index]) nfailedleptons++;
    if(!Zlep2istight[Z2index]) nfailedleptons++;
    ZZsystem = Z1+Z2;
    ZZsystemnofsr = Z1nofsr+Z2nofsr;

    /*if(abs(ZZsystemnofsr.M()-ZZsystem.M())>0.000001){
        std::cout<<"FSR works "<<abs(ZZsystemnofsr.M()-ZZsystem.M())<<std::endl;
        std::cout<<"FSR: "<<ZZsystem.M()<<" noFSR:"<<ZZsystemnofsr.M()<<std::endl;
    }*/
    
    float massZZ;
    if (isFSR) massZZ = ZZsystem.M();
    else massZZ = ZZsystemnofsr.M();
    if ((massZZ>HiggscutDown)&&(massZZ<HiggscutUp)){
        if (flag2e2mu) cutm4l2e2mu++;
        if (flag4e) cutm4l4e++;
        if (flag4mu) cutm4l4mu++;
    }

    unsigned int jet1index, jet2index;
    jet1index = 99;
    jet2index = 99;
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
            if(Jet_pt[jetidx[1]]>Jet_pt[jetidx[0]])
            {
                jet1index = jetidx[1];
                jet2index = jetidx[0];
            }
        }
        if(jetidx.size()>2)
        {
            jet1index = jetidx[0];
            jet2index = jetidx[1];
            if(Jet_pt[jetidx[1]]>Jet_pt[jetidx[0]])
            {
                jet1index = jetidx[1];
                jet2index = jetidx[0];
            }
            for (unsigned int pj=2;pj<jetidx.size();pj++){
                if((Jet_pt[jetidx[pj]]>jet1index)&&(Jet_pt[jetidx[pj]]>jet2index)){
                    jet1index = jetidx[pj];
                }
                if(Jet_pt[jetidx[pj]]>jet2index){
                    jet2index = jetidx[pj];
                }
            }
        }
    }
    TLorentzVector Jet1,Jet2;
    SimpleParticleCollection_t associated;
    if(jetidx.size()>0){
        Jet1.SetPtEtaPhiM(Jet_pt[jet1index],Jet_eta[jet1index],Jet_phi[jet1index],Jet_mass[jet1index]);
        associated.push_back(SimpleParticle_t(0, Jet1));
        pTj1 = Jet1.Pt();
        etaj1 = Jet1.Eta();
        phij1 = Jet1.Phi();
        mj1 = Jet1.M();
        if(jetidx.size()>1){
           Jet2.SetPtEtaPhiM(Jet_pt[jet2index],Jet_eta[jet2index],Jet_phi[jet2index],Jet_mass[jet2index]);
           associated.push_back(SimpleParticle_t(0, Jet2));
           pTj2 = Jet2.Pt();
           etaj2 = Jet2.Eta();
           phij2 = Jet2.Phi();
           mj2 = Jet2.M();
        }
    }
    
    
        
    SimpleParticleCollection_t daughters;
    TLorentzVector Lep1,Lep2,Lep3,Lep4;
    
    Lep1.SetPtEtaPhiM(Zlep1pt[Z1index],Zlep1eta[Z1index],Zlep1phi[Z1index],Zlep1mass[Z1index]);
    Lep2.SetPtEtaPhiM(Zlep2pt[Z1index],Zlep2eta[Z1index],Zlep2phi[Z1index],Zlep2mass[Z1index]);
    Lep3.SetPtEtaPhiM(Zlep1pt[Z2index],Zlep1eta[Z2index],Zlep1phi[Z2index],Zlep1mass[Z2index]);
    Lep4.SetPtEtaPhiM(Zlep2pt[Z2index],Zlep2eta[Z2index],Zlep2phi[Z2index],Zlep2mass[Z2index]);
   
    lep_Hindex[0] = Zlep1lepindex[Z1index];
    lep_Hindex[1] = Zlep2lepindex[Z1index];
    lep_Hindex[2] = Zlep1lepindex[Z2index];
    lep_Hindex[3] = Zlep2lepindex[Z2index];
    pTL1 = Lep1.Pt();
    etaL1 = Lep1.Eta();
    phiL1 = Lep1.Phi();
    massL1 = Lep1.M();
    pTL2 = Lep2.Pt();
    etaL2 = Lep2.Eta();
    phiL2 = Lep2.Phi();
    massL2 = Lep2.M();
    pTL3 = Lep3.Pt();
    etaL3 = Lep3.Eta();
    phiL3 = Lep3.Phi();
    massL3 = Lep3.M();
    pTL4 = Lep4.Pt();
    etaL4 = Lep4.Eta();
    phiL4 = Lep4.Phi();
    massL4 = Lep4.M();
    
    
    daughters.push_back(SimpleParticle_t((-1)*Zflavor[Z1index]*Zlep1chg[Z1index], Lep1));
    daughters.push_back(SimpleParticle_t((-1)*Zflavor[Z1index]*Zlep2chg[Z1index], Lep2));
    daughters.push_back(SimpleParticle_t((-1)*Zflavor[Z2index]*Zlep1chg[Z2index], Lep3));
    daughters.push_back(SimpleParticle_t((-1)*Zflavor[Z2index]*Zlep2chg[Z2index], Lep4));
    me_0plus_JHU=999.0; me_qqZZ_MCFM=999.0; p0plus_m4l=999.0; bkg_m4l=999.0; D_bkg_kin=999.0; D_bkg=999.0;
    D_bkg_kin_vtx_BS=999.0;
    
    p0minus_VAJHU=999.0; pg1g4_VAJHU=999.0; Dgg10_VAMCFM=999.0; D_g4=999.0; D_g1g4=999.0; D_0m=999.0; D_CP=999.0; D_0hp=999; D_int=999.0;D_L1=999.0; D_L1_int=999.0; D_L1Zg=999.0; D_L1Zgint=999.0;
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
