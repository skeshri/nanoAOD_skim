# NanoAOD Skim
nanoAOD skiming code for H->ZZ->2l2Q studies.

## Code setup

1. Step: 1: Get CMSSW release

   ```bash
   cmsrel CMSSW_10_6_20
   cd CMSSW_10_6_20/src
   cmsenv
   ```

2. Step: 2: Get  official nanoAODTools

   ```bash
   git clone git@github.com:cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
   cd PhysicsTools/NanoAODTools
   git checkout e963c7080606607777b083f59d752fe67b766265
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
   voms-proxy-init --voms cms --valid 168:00
   ```

   (Optional: Fix git repo)

   ```bash
   find PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/.git/ -name "*.py*" -delete
   ```

4. Step: 4: interactive running

   ```bash
   git clone -b v2.3.5 https://github.com/JHUGen/JHUGenMELA
   sh JHUGenMELA/MELA/setup.sh -j 8
   cd JHUGenMELA/MELA
   make
   cd $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
   python post_proc.py
   ```

