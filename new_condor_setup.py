import subprocess
import tarfile
import os
from datetime import datetime

current_datetime = datetime.now()

#Initial_path = '/eos/uscms/store/user/rasharma/nanoAOD_skim/2018/Test'
Initial_path = '/eos/uscms/store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2018'

condor_file_name = 'submit_condor_jobs_lnujj'

# Function to create a tar file
def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
            tar.add(source_dir, arcname=os.path.basename(source_dir))


dirName =(str(current_datetime.year)[-2:]
         +str(format(current_datetime.month,'02d'))
         +str(format(current_datetime.day,'02d'))
         +"_"
         +str(format(current_datetime.hour,'02d'))
         +str(format(current_datetime.minute,'02d'))
         +str(format(current_datetime.second,'02d'))
         )
print dirName
output_log_path = 'condor_logs'+os.sep+dirName
os.system('mkdir -p '+output_log_path)
    
post_proc_to_run = "post_proc.py"
command = "python "+post_proc_to_run
# Get CMSSW directory path and name
cmsswDirPath = os.environ['CMSSW_BASE']
CMSSWRel = cmsswDirPath.split("/")[-1]

# create tarball of present working CMSSW base directory
os.system('rm -f CMSSW*.tgz')
make_tarfile(CMSSWRel+".tgz", cmsswDirPath)

print("copying the ",CMSSWRel+".tgz  file to eos...\n")
os.system('xrdcp -f ' + CMSSWRel+".tgz" + ' root://cmseos.fnal.gov/'+Initial_path+'/' + CMSSWRel+".tgz")

with open('input_data_Files/samples_list_das.dat') as in_file:
  count = 0
  outjdl_file = open(condor_file_name+".jdl","w")
  outjdl_file.write("Executable = "+condor_file_name+".sh\n")
  outjdl_file.write("Universe = vanilla\n")
  outjdl_file.write("Notification = ERROR\n")
  outjdl_file.write("Should_Transfer_Files = YES\n")
  outjdl_file.write("WhenToTransferOutput = ON_EXIT\n")
  outjdl_file.write("Transfer_Input_Files = Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt, "+post_proc_to_run+"\n")
  outjdl_file.write("x509userproxy = $ENV(X509_USER_PROXY)\n")
  for lines in in_file:
     if lines[0] == "#": continue
     count = count +1
     #if count > 1: break
     print "="*51,"\n"
     print "==>  Sample : ",count
     sample_name = lines.split('/')[1]
     campaign = lines.split('/')[2].split('-')[0]
     print "==> sample_name = ",sample_name
     print "==> campaign = ",campaign
     ########################################
     #
     #      Create output directory
     #
     ########################################
     if sample_name.find("SingleMuon") != -1 or sample_name.find("EGamma") != -1:
       output_string = sample_name + os.sep + campaign + os.sep + dirName
       output_path = Initial_path + os.sep + output_string
       os.system("xrdfs root://cmseos.fnal.gov/ mkdir "+Initial_path + os.sep + sample_name)
       os.system("xrdfs root://cmseos.fnal.gov/ mkdir "+Initial_path + os.sep + sample_name + os.sep + campaign)
       os.system("xrdfs root://cmseos.fnal.gov/ mkdir "+ Initial_path + os.sep + sample_name + os.sep + campaign + os.sep + dirName)
     else:
       output_string = sample_name+os.sep+dirName
       output_path = Initial_path+output_string
       os.system("xrdfs root://cmseos.fnal.gov/ mkdir "+Initial_path + os.sep + sample_name)
       os.system("xrdfs root://cmseos.fnal.gov/ mkdir "+Initial_path + os.sep + sample_name+os.sep+dirName)
     print "output_path = ",output_path

     ########################################
     #print 'dasgoclient --query="file dataset='+lines.strip()+'"'
     #print "..."
     output = os.popen('dasgoclient --query="file dataset='+lines.strip()+'"').read()

     for root_file in output.split():
       #print "=> ",root_file
       outjdl_file.write("Output = "+output_log_path+"/"+sample_name+"_$(Process).stdout\n")
       outjdl_file.write("Error  = "+output_log_path+"/"+sample_name+"_$(Process).stdout\n")
       outjdl_file.write("Log  = "+output_log_path+"/"+sample_name+"_$(Process).log\n")
       outjdl_file.write("Arguments = "+("root://cms-xrd-global.cern.ch/"+root_file).replace('/','\/')+" "+output_path+"  "+Initial_path+"\n")
       outjdl_file.write("Queue \n")
  outjdl_file.close();


outScript = open(condor_file_name+".sh","w");
outScript.write('#!/bin/bash');
outScript.write("\n"+'echo "Starting job on " `date`');
outScript.write("\n"+'echo "Running on: `uname -a`"');
outScript.write("\n"+'echo "System software: `cat /etc/redhat-release`"');
outScript.write("\n"+'source /cvmfs/cms.cern.ch/cmsset_default.sh');
outScript.write("\n"+'echo "copy cmssw tar file from store area"');
outScript.write("\n"+'xrdcp -s root://cmseos.fnal.gov/${3}/'+CMSSWRel +'.tgz  .');
outScript.write("\n"+'tar -xf '+ CMSSWRel +'.tgz' );
outScript.write("\n"+'rm '+ CMSSWRel +'.tgz' );
outScript.write("\n"+'cd ' + CMSSWRel + '/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_vvVBS/' );
#outScript.write("\n"+'echo "====> List files : " ');
#outScript.write("\n"+'ls -alh');
outScript.write("\n"+'rm *.root');
outScript.write("\n"+'scramv1 b ProjectRename');
outScript.write("\n"+'eval `scram runtime -sh`');
outScript.write("\n"+'sed -i "s/testfile = .*/testfile = \\"${1}\\"/g" '+post_proc_to_run);
#sed 's/testfile = .*/testfile = "root:\/\/cms-xrd-global.cern.ch\/\/store\/data\/Run2018A\/EGamma\/NANOAOD\/Nano1June2019-v1\/40000\/7F8E5EFA-DF69-3442-8741-B04F7313B3B1.root"/g' post_proc.py
#outScript.write("\n"+'sed  -i "s/\/store\/mc\/RunIIFall17NanoAODv5\/ZZTo4L_13TeV_powheg_pythia8\/NANOAODSIM\/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1\/110000\/00BB722F-1A9E-C64B-B2D0-A5818F1661C3.root/${1}/g" post_proc_zprime.py');
outScript.write("\n"+command);
outScript.write("\n"+'echo "====> List root files : " ');
outScript.write("\n"+'ls *.root');
outScript.write("\n"+'echo "====> copying *.root file to stores area..." ');
outScript.write("\n"+'xrdcp -f *.root root://cmseos.fnal.gov/${2}');
outScript.write("\n"+'rm *.root');
outScript.write("\n"+'cd ${_CONDOR_SCRATCH_DIR}');
outScript.write("\n"+'rm -rf ' + CMSSWRel);
outScript.write("\n");
outScript.close();
os.system("chmod 777 "+condor_file_name+".sh");

print "===> Set Proxy Using:";
print "\tvoms-proxy-init --voms cms --valid 168:00";
print "\"condor_submit "+condor_file_name+".jdl\" to submit";
#os.system("condor_submit "+condor_file_name+".jdl")
