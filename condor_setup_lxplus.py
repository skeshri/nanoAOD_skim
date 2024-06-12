"""
# How to run:
python3 condor_setup_lxplus.py
"""
import argparse
import os
import sys
import datetime

sys.path.append("Utils/.")

from color_style import style

def main(args):

    # Variables from argparse
    submission_name = args.submission_name
    use_custom_eos = args.use_custom_eos
    use_custom_eos_cmd = args.use_custom_eos_cmd
    InputFileFromWhereReadDASNames = args.input_file
    EOS_Output_path = args.eos_output_path
    if EOS_Output_path == "":
        # Get the username and its initial and set the path as /eos/user/<UserInitials>/<UserName>/nanoAOD_ntuples
        username = os.environ['USER']
        user_initials = username[0:1]
        EOS_Output_path = '/eos/user/'+user_initials+'/'+username+'/nanoAOD_ntuples'
    if submission_name != "":
        EOS_Output_path = EOS_Output_path + '/' + submission_name
    condor_log_path = args.condor_log_path

    # Get top-level directory name from PWD
    TOP_LEVEL_DIR_NAME = os.path.basename(os.getcwd())
    condor_file_name = args.condor_file_name
    # Add time stamp to the condor_file_name
    now = datetime.datetime.now()
    condor_file_name = condor_file_name + "_" + now.strftime("%Y%m%d_%H%M%S")
    condor_file_name = condor_file_name + "_" + submission_name

    condor_queue = args.condor_queue
    DontCreateTarFile = args.DontCreateTarFile

    # Create log files
    import infoCreaterGit
    SumamryOfCurrentSubmission = raw_input("\n\nWrite summary for current job submission: ")
    infoLogFiles = infoCreaterGit.BasicInfoCreater('summary.dat',SumamryOfCurrentSubmission)
    infoLogFiles.generate_git_patch_and_log()

    # Get CMSSW directory path and name
    cmsswDirPath = os.environ['CMSSW_BASE']
    CMSSWRel = cmsswDirPath.split("/")[-1]

    # Create directories for storing log files and output files at EOS.
    import fileshelper
    dirsToCreate = fileshelper.FileHelper( (condor_log_path + '/condor_logs/'+submission_name).replace("//","/"), EOS_Output_path)
    output_log_path = dirsToCreate.create_log_dir_with_date()
    storeDir = dirsToCreate.create_store_area(EOS_Output_path)
    dirName = dirsToCreate.dir_name

    # create tarball of present working CMSSW base directory
    if not DontCreateTarFile: os.system('rm -f CMSSW*.tgz')
    import makeTarFile
    if not DontCreateTarFile: makeTarFile.make_tarfile(cmsswDirPath, CMSSWRel+".tgz")
    print("copying the "+CMSSWRel+".tgz  file to eos path: "+storeDir+"\n")
    os.system('xrdcp ' + CMSSWRel+".tgz" + '  root://eosuser.cern.ch/'+storeDir+'/' + CMSSWRel+".tgz")

    post_proc_to_run = "post_proc.py"
    command = "python "+post_proc_to_run
    condor_arguments_list = []  # A list that contains all the arguments to be passed for each job

    with open('input_data_Files/'+InputFileFromWhereReadDASNames) as in_file:
        outjdl_file = open(condor_file_name+".jdl","w")
        outjdl_file.write("+JobFlavour   = \""+condor_queue+"\"\n")
        outjdl_file.write("Executable = "+condor_file_name+".sh\n")
        outjdl_file.write("Universe = vanilla\n")
        outjdl_file.write("Notification = ERROR\n")
        outjdl_file.write("Should_Transfer_Files = NO\n")
        outjdl_file.write("x509userproxy = $ENV(X509_USER_PROXY)\n")
        outjdl_file.write("MY.SingularityImage = \"/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-cat/cmssw-lxplus/cmssw-el7-lxplus:latest/\"\n")
        outjdl_file.write("Output = {output_log_path}/$(logtxt)_$(Process).stdout\n".format(output_log_path=output_log_path))
        outjdl_file.write("Error  = {output_log_path}/$(logtxt)_$(Process).err\n".format(output_log_path=output_log_path))
        outjdl_file.write("Log  = {output_log_path}/$(logtxt)_$(Process).log\n".format(output_log_path=output_log_path))
        outjdl_file.write('Arguments = "$(infile) $(outfile) $(eospath) $(outfilename)"\n')
        count = 0
        count_jobs = 0
        output_string_list = []
        for SampleDASName in in_file:
            if SampleDASName[0] == "#": continue
            count = count +1
            # if count > 1: break
            print(style.RED +"="*51+style.RESET+"\n")
            print ("==> Sample : ",count)
            sample_name = SampleDASName.split('/')[1]
            # if SampleDASName contains `ext` then add it
            if 'ext' in SampleDASName:
                sample_name = sample_name + '_ext' + SampleDASName.split('ext')[-1].split('/')[0]
            print("==> sample_name = ",sample_name)
            campaign = SampleDASName.split('/')[2].split('-')[0]
            campaign_ifNeeded = SampleDASName.split('/')[2].split('-')[1]
            print("==> campaign = ",campaign)
            ########################################
            #
            #      Create output directory
            #
            ########################################
            if (SampleDASName.strip()).endswith("/NANOAOD"): # if the sample name ends with /NANOAOD, then it is a data if it ends with /NANOAODSIM then it is a MC. As the line contain the "\n" at the end, so we need to use the strip() function.
                output_string = sample_name + os.sep + campaign + os.sep + dirName
                # if output_string is already present in the output_string_list, then append the campaign_ifNeeded
                if output_string in output_string_list:
                    output_string = sample_name + os.sep + campaign + os.sep + dirName + "_" + campaign_ifNeeded
                output_path = EOS_Output_path + os.sep + output_string
                print("==> output_path = ",output_path)
                os.system("mkdir -p "+ output_path)
                infoLogFiles.send_git_log_and_patch_to_eos(output_path)
            else:
                output_string = campaign + os.sep + sample_name + os.sep + dirName
                # if output_string is already present in the output_string_list, then append the campaign_ifNeeded
                if output_string in output_string_list:
                    output_string = campaign + os.sep + sample_name + os.sep + dirName + "_" + campaign_ifNeeded
                    if output_string in output_string_list:
                        output_string = campaign + os.sep + sample_name + os.sep + dirName + "_" + campaign_ifNeeded + "_ext" + SampleDASName.split('ext')[-1].split('/')[0]
                output_path = EOS_Output_path+ os.sep + output_string
                print("==> output_path = ",output_path)
                os.system("mkdir -p "+output_path)
                infoLogFiles.send_git_log_and_patch_to_eos(output_path)

            output_string_list.append(output_string)

            #  print "==> output_path = ",output_path

            ########################################
            # print 'dasgoclient --query="file dataset='+SampleDASName.strip()+'"'
            # print "..."
            if use_custom_eos:
                xrd_redirector = 'root://cms-xrd-global.cern.ch/'
                output = os.popen(use_custom_eos_cmd + SampleDASName.strip()).read()
            else:
                xrd_redirector = 'root://cms-xrd-global.cern.ch/'
                output = os.popen('dasgoclient --query="file dataset='+SampleDASName.strip()+'"').read()

            count_root_files = 0
            for root_file in output.split():
                # print "=> ",root_file
                count_root_files+=1
                count_jobs += 1
                condor_arguments_list.append(
                    (
                        xrd_redirector + root_file,
                        output_path,
                        EOS_Output_path,
                        (root_file.split("/")[-1]).split(".")[0],
                        output_path.split("/")[-2], # This argument is used for the log file name
                    )
                )
                if args.debug:
                    # break the for loop after 1 iteration to submit only 1 job
                    break
            if args.debug:
                # break the for loop after 1 iteration to submit only 1 job
                break
            print("Number of files: ",count_root_files)
            print("Number of jobs (till now): ",count_jobs)
        outjdl_file.write("queue infile, outfile, eospath, outfilename, logtxt from {}.txt\n".format(condor_file_name))
        outjdl_file.close();

    # Write all condor jobs arguments from list to a file with same name as condor_file_name but with .txt extension
    with open(condor_file_name+".txt", "w") as f:
        for item in condor_arguments_list:
            f.write("{}\n".format(",".join(item)))

    outScript = open(condor_file_name+".sh","w");
    outScript.write('#!/bin/bash');
    outScript.write("\n"+'echo "Starting job on " `date`');
    outScript.write("\n"+'echo "Running on: `uname -a`"');
    outScript.write("\n"+'echo "System software: `cat /etc/redhat-release`"');
    outScript.write("\n"+'source /cvmfs/cms.cern.ch/cmsset_default.sh');
    outScript.write("\n"+'echo "====> List input arguments : " ');
    outScript.write("\n"+'echo "1. nanoAOD ROOT file: ${1}"');
    outScript.write("\n"+'echo "2. EOS path to store output root file: ${2}"');
    outScript.write("\n"+'echo "3. EOS path from where we copy CMSSW: ${3}"');
    outScript.write("\n"+'echo "4. Output root file name: ${4}"');
    outScript.write("\n"+'echo "========================================="');
    outScript.write("\n"+'echo "copy cmssw tar file from store area"');
    outScript.write("\n"+'xrdcp -f  root://eosuser.cern.ch/${3}/'+CMSSWRel +'.tgz  .');
    outScript.write("\n"+'tar -xf '+ CMSSWRel +'.tgz' );
    outScript.write("\n"+'rm '+ CMSSWRel +'.tgz' );
    outScript.write("\n"+'cd ' + CMSSWRel + '/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/'+TOP_LEVEL_DIR_NAME+'/' );
    outScript.write("\n"+'rm *.root');
    outScript.write("\n"+'scramv1 b ProjectRename');
    outScript.write("\n"+'eval `scram runtime -sh`');
    outScript.write("\n"+'echo "========================================="');
    outScript.write("\n"+'echo "cat post_proc.py"');
    outScript.write("\n"+'echo "..."');
    outScript.write("\n"+'cat post_proc.py');
    outScript.write("\n"+'echo "..."');
    outScript.write("\n"+'echo "========================================="');
    if args.NOsyst:
        outScript.write(
            "\n"
            + command
            + " --entriesToRun 0  --inputFile ${1} --outputFile ${4}_hadd.root --cutFlowFile ${4}.json --DownloadFileToLocalThenRun True  --NOsyst"
        )
    else:
        outScript.write(
            "\n"
            + command
            + " --entriesToRun 0  --inputFile ${1} --outputFile ${4}_hadd.root --cutFlowFile ${4}.json --DownloadFileToLocalThenRun True"
        )
    outScript.write("\n"+'echo "====> List root files : " ');
    outScript.write("\n"+'ls -ltrh *.root');
    outScript.write("\n"+'ls -ltrh *.json');
    outScript.write("\n"+'echo "====> copying *.root file to stores area..." ');
    outScript.write("\n"+'if ls ${4}_hadd.root 1> /dev/null 2>&1; then');
    outScript.write("\n"+'    echo "File ${4}_hadd.root exists. Copy this."');
    outScript.write("\n"+'    echo "xrdcp -f ${4}_hadd.root  root://eosuser.cern.ch/${2}/${4}_Skim.root"');
    outScript.write("\n"+'    xrdcp -f ${4}_hadd.root  root://eosuser.cern.ch/${2}/${4}_Skim.root');
    outScript.write("\n"+'    echo "xrdcp -f ${4}.json  root://eosuser.cern.ch/${2}/cutFlow_${4}.json"');
    outScript.write("\n"+'    xrdcp -f ${4}.json  root://eosuser.cern.ch/${2}/cutFlow_${4}.json');
    outScript.write("\n"+'else');
    outScript.write("\n"+'    echo "Something wrong: file ${4}_hadd.root does not exists, please check the post_proc.py script."');
    outScript.write("\n"+'fi');
    outScript.write("\n"+'rm *.root');
    outScript.write("\n"+'cd ${_CONDOR_SCRATCH_DIR}');
    outScript.write("\n"+'rm -rf ' + CMSSWRel);
    outScript.write("\n");
    outScript.close();
    os.system("chmod 777 "+condor_file_name+".sh");

    print("\n#===> Set Proxy Using:")
    print("voms-proxy-init --voms cms --valid 168:00")
    print("\n# It is assumed that the proxy is created in file: /tmp/x509up_u48539. Update this in below two lines:")
    print("cp /tmp/x509up_u48539 ~/")
    print("export X509_USER_PROXY=~/x509up_u48539")
    print("\n#Submit jobs:")
    print("condor_submit "+condor_file_name+".jdl")
    # os.system("condor_submit "+condor_file_name+".jdl")

# Below patch is to format the help command as it is
class PreserveWhitespaceFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Condor Job Submission", formatter_class=PreserveWhitespaceFormatter)
    parser.add_argument("--submission_name", default="SkimNanoAOD", help="String to be changed by user.")
    parser.add_argument("--use_custom_eos", default=False, action='store_true', help="Use custom EOS.")
    parser.add_argument("--DontCreateTarFile", default=False, action='store_true', help="Create tar file of CMSSW directory.")
    parser.add_argument("--use_custom_eos_cmd", default='eos root://cmseos.fnal.gov find -name "*.root" /store/group/lnujj/VVjj_aQGC/custom_nanoAOD', help="Custom EOS command.")
    parser.add_argument("--input_file", default='', required=True,  help="Input file from where to read DAS names.")
    parser.add_argument("--eos_output_path", default='', help="EOS path for output files. By default it is `/eos/user/<UserInitials>/<UserName>/nanoAOD_ntuples`")
    parser.add_argument("--condor_log_path", default='./', help="Path where condor log should be saved. By default is the current working directory")
    parser.add_argument("--condor_file_name", default='submit_condor_jobs', help="Name for the condor file.")
    parser.add_argument("--condor_queue", default="tomorrow", help="""
                        Condor queue options: (Reference: https://twiki.cern.ch/twiki/bin/view/ABPComputing/LxbatchHTCondor#Queue_Flavours)

                        name            Duration
                        ------------------------
                        espresso            20min
                        microcentury     1h
                        longlunch           2h
                        workday 8h        1nd
                        tomorrow           1d
                        testmatch          3d
                        nextweek           1w
                        """)

    parser.add_argument("--post_proc", default="post_proc.py", help="Post process script to run.")
    parser.add_argument("--transfer_input_files", default="keep_and_drop.txt", help="Files to be transferred as input.")
    parser.add_argument("--NOsyst", default=False, action='store_true', help="Run without systematics.")
    parser.add_argument("--debug", default=False, action='store_true', help="Debug mode.")

    args = parser.parse_args()
    main(args)
# condor_setup_lxplus.py
