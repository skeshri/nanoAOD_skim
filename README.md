# nanoAOD_vvVBS
nanoAOD skiming code for vv semi-leptonic VBS studies


## Code setup

1. Step: 1: Get CMSSW release

   ```bash
   cmsrel CMSSW_10_2_15
   cd CMSSW_10_2_15/src
   cmsenv
   ```
   
2. Step: 2: Get  official nanoAODTools

   ```bash
   git clone git@github.com:cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
   ```
   
3. Step: 3: Get our analysis repository

   ```bash
   git clone git@github.com:ram1123/nanoAOD_vvVBS.git PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_vvVBS
   cmsenv
   scram b
   voms-proxy-init -voms cms
   ```
   
4. Step: 4: interactive running

   ```bash
   cd PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_vvVBS   
   python post_proc.py
   ```
   
5. Step: 5 (a): Crab-job submission   

   ```bash
   cd crab/
   voms-proxy-init -voms cms --valid 200:00
   source /cvmfs/cms.cern.ch/crab3/crab.sh
   crab submit -c crab_cfg.py
   ```
   
   or,
   
   Step: 5 (b): Condor-job submission

   ```bash
   python new_condor_setup.py
   voms-proxy-init -voms cms --valid 200:00
   condor_submit new_condor.jdl
   ```

# To-Do List

