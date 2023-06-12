#include "../interface/H4LTools.h"
#include <TLorentzVector.h>
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
        if (((*Muon_pt)[i]>MuPtcut)&&(fabs((*Muon_eta)[i])<2.4)&&(((*Muon_isGlobal)[i]||(*Muon_isTracker)[i]||(*Muon_isPFcand)[i]))){
            LooseMuonindex.push_back(i);
        }
    }

    return LooseMuonindex;
}
std::vector<unsigned int> H4LTools::goodMuons2015_noIso_noPf(std::vector<unsigned int> Muonindex){
    std::vector<unsigned int> bestMuonindex;
    //std::cout<<Muonindex.size()<<std::endl;
    for (unsigned int i=0; i<Muonindex.size(); i++){
        if (((*Muon_pt)[Muonindex[i]]>MuPtcut)&&(fabs((*Muon_eta)[Muonindex[i]])<2.4)&&((*Muon_isGlobal)[Muonindex[i]]||(*Muon_isTracker)[Muonindex[i]])){
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
        if ((*Muon_pt)[i]<200){
            tightid.push_back((*Muon_isPFcand)[i]);
        }
        else{
            tightid.push_back((*Muon_isPFcand)[i]||((((*Muon_ptErr)[i]/(*Muon_pt)[i])<0.3)&&(fabs((*Muon_dxy)[i])<0.2)&&(fabs((*Muon_dz)[i])<0.5)&&((*Muon_nTrackerLayers)[i]>5)));
        }

    }

    return tightid;
}

bool H4LTools::ZZSelection(){
    bool foundZZCandidate = false;
 
    std::vector<unsigned int> looseEle,looseMu,bestEle,bestMu;
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
    

    TLorentzVector z1,z2;
    std::vector<TLorentzVector> Elelist;
    std::vector<TLorentzVector> Mulist;
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
        Mu.SetPtEtaPhiM((*Muon_pt)[Muonindex[imu]],(*Muon_eta)[Muonindex[imu]],(*Muon_phi)[Muonindex[imu]],(*Muon_mass)[Muonindex[imu]]);
        Mulist.push_back(Mu);
        muid.push_back(AllMuid[Muonindex[imu]]);
        Muiso.push_back((*Muon_pfRelIso03_all)[Muonindex[imu]]);
    }

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
        if((muid[amu]==true)&&(Muiso[amu]<0.35)){
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

    if(TightEleindex.size()>1){
        for(unsigned int ke=0; ke<(TightEleindex.size()-1);ke++){
            for(unsigned int je=ke+1;je<TightEleindex.size();je++){
                if ((Elechg[TightEleindex[ke]]+Elechg[TightEleindex[je]])==0){
                    TLorentzVector Zcan;
                    Zcan = Elelist[TightEleindex[ke]] + Elelist[TightEleindex[je]];
                    if((Zcan.M()>12)&&(Zcan.M()<120)){
                        Zlist.push_back(Zcan);
                        Zlep1index.push_back(TightEleindex[ke]);
                        Zlep2index.push_back(TightEleindex[je]);
                        Zflavor.push_back(11);
                        Zlep1pt.push_back(Elelist[TightEleindex[ke]].Pt());
                        Zlep2pt.push_back(Elelist[TightEleindex[je]].Pt());
                        Zlep1eta.push_back(Elelist[TightEleindex[ke]].Eta());
                        Zlep2eta.push_back(Elelist[TightEleindex[je]].Eta());
                        Zlep1phi.push_back(Elelist[TightEleindex[ke]].Phi());
                        Zlep2phi.push_back(Elelist[TightEleindex[je]].Phi());
                        Zlep1mass.push_back(Elelist[TightEleindex[ke]].M());
                        Zlep2mass.push_back(Elelist[TightEleindex[je]].M());
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
                    Zcan = Mulist[TightMuindex[kmu]] + Mulist[TightMuindex[jmu]];
                    if((Zcan.M()>12)&&(Zcan.M()<120)){
                        Zlist.push_back(Zcan);
                        Zlep1index.push_back(TightMuindex[kmu]);
                        Zlep2index.push_back(TightMuindex[jmu]);
                        Zflavor.push_back(13);
                        Zlep1pt.push_back(Mulist[TightMuindex[kmu]].Pt());
                        Zlep2pt.push_back(Mulist[TightMuindex[jmu]].Pt());
                        Zlep1eta.push_back(Mulist[TightMuindex[kmu]].Eta());
                        Zlep2eta.push_back(Mulist[TightMuindex[jmu]].Eta());
                        Zlep1phi.push_back(Mulist[TightMuindex[kmu]].Phi());
                        Zlep2phi.push_back(Mulist[TightMuindex[jmu]].Phi());
                        Zlep1mass.push_back(Mulist[TightMuindex[kmu]].M());
                        Zlep2mass.push_back(Mulist[TightMuindex[jmu]].M());
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
                lepA.SetPtEtaPhiM(Zlep1pt[m],Zlep1eta[m],Zlep1phi[m],Zlep1mass[m]);
                lepB.SetPtEtaPhiM(Zlep1pt[n],Zlep1eta[n],Zlep1phi[n],Zlep1mass[n]);
                lepAB = lepA + lepB;
                if(lepAB.M()<4) continue;  //QCD Supressionas
            }
            if ((Zlep1chg[m]+Zlep2chg[n])==0){
                TLorentzVector lepA,lepB,lepAB;
                lepA.SetPtEtaPhiM(Zlep1pt[m],Zlep1eta[m],Zlep1phi[m],Zlep1mass[m]);
                lepB.SetPtEtaPhiM(Zlep2pt[n],Zlep2eta[n],Zlep2phi[n],Zlep2mass[n]);
                lepAB = lepA + lepB;
                if(lepAB.M()<4) continue;
            }
            if ((Zlep2chg[m]+Zlep1chg[n])==0){
                TLorentzVector lepA,lepB,lepAB;
                lepA.SetPtEtaPhiM(Zlep2pt[m],Zlep2eta[m],Zlep2phi[m],Zlep2mass[m]);
                lepB.SetPtEtaPhiM(Zlep1pt[n],Zlep1eta[n],Zlep1phi[n],Zlep1mass[n]);
                lepAB = lepA + lepB;
                if(lepAB.M()<4) continue;
            }
            if ((Zlep2chg[m]+Zlep2chg[n])==0){
                TLorentzVector lepA,lepB,lepAB;
                lepA.SetPtEtaPhiM(Zlep2pt[m],Zlep2eta[m],Zlep2phi[m],Zlep2mass[m]);
                lepB.SetPtEtaPhiM(Zlep1pt[n],Zlep1eta[n],Zlep1phi[n],Zlep1mass[n]);
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
        


    return foundZZCandidate;


}