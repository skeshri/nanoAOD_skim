#include <iostream>
#include <set>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>

void removeDuplicates_new() {
    //TString prefix = "/raid/raid9/qguo/Run2/after/Run2_2/new/CMSSW_10_2_18/src/data_2018_NotBestMelaCand";
    //TString prefix = "/publicfs/cms/data/hzz/guoqy/newNTuple_UL/2018/Data/DataUL2018_all";
    //TString prefix = "/publicfs/cms/user/qyguo/lacked_Data1718/18/DoubleMuon_Run18A_1";
    TString prefix = "/eos/user/y/yujil/SkimNanoAOD_2022_ZXCRv4_Data/DataZXCR_2022";
    TString filename = prefix+".root";

    std::cout<<filename<<std::endl;

    TFile *oldfile = new TFile(filename);
    //oldfile->cd("Ana");
    //TTree *oldtree = (TTree*)oldfile->Get("Ana/passedEvents");
    
    TTree *oldtree = (TTree*)gDirectory->Get("Events");
    //TTree *oldtree = (TTree*)oldfile->Get("passedEvents");

    Long64_t nentries = oldtree->GetEntries();
    std::cout<<nentries<<" total entries."<<std::endl;
    UInt_t Run, LumiSect;
    ULong64_t Event;
    bool passedZ4lSelection;
    oldtree->SetBranchAddress("run",&Run);
    oldtree->SetBranchAddress("luminosityBlock",&LumiSect);
    oldtree->SetBranchAddress("event",&Event);

    //Create a new file + a clone of old tree in new file
    

    std::set<TString> runlumieventSet;
    std::set<UInt_t> runSet;
    std::set<UInt_t> lumiSet;
    std::set<ULong64_t> eventSet;
    int nremoved = 0;
    int filetag = 0;
    for (int j=0; j<(nentries/100000.);j++){
        TFile *newfile = new TFile(
                prefix+"_noDuplicates_"+filetag+".root"
                ,"recreate");
        TTree *newtree = oldtree->CloneTree(0);
            
        for (Long64_t i=100000*j;i<100000*(j+1); i++) {
            if (i==nentries) break;
            if (i%100000==0) std::cout<<i<<"/"<<nentries<<std::endl;
            oldtree->GetEntry(i);

            /*TString s_Run  = std::to_string(Run);
            TString s_Lumi = std::to_string(LumiSect);
            TString s_Event = std::to_string(Event);
            TString runlumievent = s_Run+":"+s_Lumi+":"+s_Event;*/
            
            if ((runSet.find(Run)==runSet.end()) || (lumiSet.find(LumiSect)==lumiSet.end()) || (eventSet.find(Event)==eventSet.end())) {
                runSet.insert(Run);
                lumiSet.insert(LumiSect);
                eventSet.insert(Event);
                newtree->Fill();
            } else {
                nremoved++;
            }

        }
        filetag++;  
        newtree->AutoSave();
        delete newtree;
        delete newfile;
        //cout<<"Saved "<<filetag<<" files"<<endl;
    }
    std::cout<<nremoved<<" duplicates."<<std::endl;
    //newtree->Print();
    //newtree->AutoSave();
    //delete oldfile;
}

