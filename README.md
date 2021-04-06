# NanoAOD Skim
nanoAOD skiming code for H->ZZ->2l2Q studies.

## Code setup

1. Step: 1: Get CMSSW release

   ```bash
   cmsrel CMSSW_10_2_22
   cd CMSSW_10_2_22/src
   cmsenv
   ```

2. Step: 2: Get  official nanoAODTools

   ```bash
   git clone git@github.com:cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
   cd PhysicsTools/NanoAODTools
   git checkout 079c9e18c14c9d71ffe6d0cc4b42f15d97c29efc
   ```

3. Step: 3: Get our analysis repository

   ```bash
   cd $CMSSW_BASE/src
   git clone git@github.com:osWW-VBS/nanoAOD_skim.git PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
   cd PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
   git submodule init
   git submodule update
   cd -
   cmsenv
   # patch PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/nanoAOD_tools.patch
   cp PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/data/btag/*.csv PhysicsTools/NanoAODTools/data/btagSF/.
   scram b
   voms-proxy-init -voms cms
   ```

   (Optional: Fix git repo)

   ```bash
   find PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/.git/ -name "*.py*" -delete
   ```

4. Step: 4: interactive running

   ```bash
   cd PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
   python post_proc.py
   ```

5. batch job submission.
   1. Crab-job submission (Not tested)
      ```bash
      cd crab/
      voms-proxy-init -voms cms --valid 200:00
      source /cvmfs/cms.cern.ch/crab3/crab.sh
      crab submit -c crab_cfg.py
      ```

   2. Step: 5 (b): Condor-job submission
      1. In file `condor_setup.py`, specify correct input text file from which you need to take input NanoAOD DAS names. Also, updated the output EOS path. Then do following:

         ```bash
         cd $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
         # Edit condor_setup.py, then
         python condor_setup_lxplus.py
         # Set proxy before submitting the condor jobs.
         voms-proxy-init -voms cms --valid 200:00
         condor_submit <Files-created-from-above-command>.jdl
         ```

         for condor jobs submission on lxplus one need to export the proxy as shown below before submitted the condor jobs:
         ```bash
         voms-proxy-init -voms cms --valid 200:00
         ```
         then one get something like
         ```
         Created proxy in /tmp/x509up_u95168
         ```
         Now, need to set this proxy to your `X509_USER_PROXY` environment variable like:
         ```bash
         cp /tmp/x509up_u117617 /afs/cern.ch/user/<NameInitial>/<UserName>/
         export X509_USER_PROXY=/afs/cern.ch/user/<NameInitial>/<UserName>/x509up_u95168
         ```
         where x509up_u95168 would be replaced by whatever your proxy name is.


