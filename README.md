# nanoAOD_vvVBS
nanoAOD skiming code for vv semi-leptonic VBS studies

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
   git clone git@github.com:ram1123/nanoAOD_vvVBS.git PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_vvVBS
   cd PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_vvVBS
   git checkout Feature-HH
   git submodule init
   git submodule update
   cd -
   cmsenv
   scram b
   voms-proxy-init -voms cms
   ```

4. Step: 4: interactive running

   ```bash
   cd PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_vvVBS
   python post_proc_DoubleHiggs.py
   ```

5. batch job submission.
   1. Step: 5 (b): Condor-job submission

      1. Setting up a voms Proxy: To access grid files to run the macro on, we must run the following commands:

         ```bash
         cmsenv
         voms-proxy-init --voms cms --valid 168:00
         ```

         after the voms command, you should receive an output similar to:

         ```
         Created proxy in /tmp/x509up_u95168
         ```

         to set this proxy to your `X509_USER_PROXY` environment variable for the example above, simply use the command:

         ```bash
         cp /tmp/x509up_u95168 ~/
         export X509_USER_PROXY=~/x509up_u95168
         ```

         where `x590up_u95168` would be replaced by whatever your proxy name is.

      1. In file `condor_setup.py`, specify correct input text file from which you need to take input NanoAOD DAS names. Also, updated the output EOS path. Then do following:
         ```bash
         cd $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_vvVBS
         # Edit condor_setup.py, then
         python condor_setup.py  --DASNames=HH_sample_list_v7_2017_campaign.txt
         # Set proxy before submitting the condor jobs.
         condor_submit <Files-created-from-above-command>.jdl
         ```


