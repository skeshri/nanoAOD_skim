#include <iostream>
#include <unordered_set>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

void removeDuplicates()
{
    TString outputFile = "/eos/user/r/rasharma/nanoAOD_ntuples/SkimNanoAOD_2022_ZXCRv2/Data_noDuplicates.root";
    TString infiles[] = {
        "/eos/user/r/rasharma/nanoAOD_ntuples/SkimNanoAOD_2022_ZXCRv2/DoubleMuon.root",
        "/eos/user/r/rasharma/nanoAOD_ntuples/SkimNanoAOD_2022_ZXCRv2/SingleMuon.root",
        "/eos/user/r/rasharma/nanoAOD_ntuples/SkimNanoAOD_2022_ZXCRv2/MuonEG.root",
        "/eos/user/r/rasharma/nanoAOD_ntuples/SkimNanoAOD_2022_ZXCRv2/Muon.root",
        "/eos/user/r/rasharma/nanoAOD_ntuples/SkimNanoAOD_2022_ZXCRv2/EGamma.root"
        };

    TChain *oldtree = new TChain("Events");
    for (TString infile : infiles)
    {
        oldtree->Add(infile);
    }

    Long64_t nentries = oldtree->GetEntries();
    std::cout << nentries << " total entries." << std::endl;

    ULong64_t Run, LumiSect, Event;
    oldtree->SetBranchAddress("run", &Run);
    oldtree->SetBranchAddress("luminosityBlock", &LumiSect);
    oldtree->SetBranchAddress("event", &Event);

    TFile *newfile = new TFile(outputFile, "recreate");
    TTree *newtree = oldtree->CloneTree(0);

    std::unordered_set<std::string> runlumieventSet; // Changed to unordered_set for performance
    int nremoved = 0;
    for (Long64_t i = 0; i < nentries; i++)
    {
        if (i % 100000 == 0)
            std::cout << i << "/" << nentries << std::endl;
        oldtree->GetEntry(i);

        TString s_Run = std::to_string(Run);
        TString s_Lumi = std::to_string(LumiSect);
        TString s_Event = std::to_string(Event);
        TString runlumievent = s_Run + ":" + s_Lumi + ":" + s_Event;

        if (runlumieventSet.insert(runlumievent).second)
        {
            newtree->Fill();
        }
        else
        {
            nremoved++;
        }
    }

    std::cout << nremoved << " duplicates removed." << std::endl;
    newtree->Print();
    newtree->AutoSave();
    delete newfile;
}
