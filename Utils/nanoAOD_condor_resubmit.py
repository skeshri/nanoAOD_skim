# Example Usage:
#
# python Utils/nanoAOD_condor_resubmit.py -d condor_logs/SkimNanoAOD_2022_ZXCR/240312_135155/ -s /eos/user/r/rasharma/nanoAOD_ntuples/SkimNanoAOD_2022_ZXCR/ -i submit_condor_jobs_lnujj_SkimNanoAOD_2022_ZXCR.jdl -n 1
#

import os
from optparse import OptionParser
import ROOT
from ROOT import TFile
import re

DEBUG = True

def files_to_remove(files,dir):
    filelist_to_remove = []
    for file in files:
        print("Checking file: {}".format(file))
        try:
            tfile = TFile.Open(file);
        except:
            pass
        if tfile:
            if (tfile.IsZombie()):
                filelist_to_remove.append(file)
        else:
            print('File could not be opened, adding it to missing files')
            filelist_to_remove.append(file)

    if DEBUG: print("Files to remove: {}".format(filelist_to_remove))
    return filelist_to_remove

def list_files(file_name):
    with open(file_name) as f:
        datafile = f.readlines()
        file_list = []
        for line in datafile:
            if "USER_" in line:
                start_pt = line.find("output")
                end_pt = line.find(".root", start_pt + 1)
                res = line[start_pt:end_pt+5]
                if "/" in res :
                    res = res[res.find("/")+1:]
                file_list.append(res)
    return file_list

def get_files_from_jdl(path_jdl):
    with open(path_jdl) as myfile:
        content = myfile.read()

    flist = re.findall("[a-z0-9A-Z-]+.root", content)
    return flist

def list_root(directory):
    flist = []
    flistWithPath = []
    for root, directories, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith("_Skim.root"):
                flist.append(filename.replace("_Skim",""))
                fileWithPath = os.path.join(root,filename)  # Get file name with path
                flistWithPath.append(fileWithPath)
    return flist,flistWithPath

def submit_missing(InputJdlFile,resubmit=True):
    bashCommand = "condor_submit {}".format(InputJdlFile)
    if resubmit :
        print('Resubmitting now!')
        os.system(bashCommand)
    else :
        print('Ready to resubmit, please set resubmit to True if you are ready : ')
        print(bashCommand)

def prepare_runJobs_missing(FailedJobRootFile,InputJdlFile,CondorLogDir,EOSDir,Resubmit_no):
    if DEBUG: print("FailedJobRootFile: {}".format(FailedJobRootFile))
    if DEBUG: print("InputJdlFile: {}".format(InputJdlFile))
    if DEBUG: print("CondorLogDir: {}".format(CondorLogDir))
    if DEBUG: print("EOSDir: {}".format(EOSDir))

    bashCommand = "cp {}  original_{}".format(InputJdlFile,InputJdlFile)
    os.system(bashCommand)

    outjdl_fileName = InputJdlFile.replace(".jdl","")+'_resubmit_'+Resubmit_no+".jdl"

    outjdl_file = open(outjdl_fileName,"w")

    with open(InputJdlFile, 'r') as myfile:
        """Copy the main part of original jdl file to new jdl file.
        All the lines before "Output = " should be copied to new jdl file.
        """
        for line in myfile:
            # Check if line starts with "Output = "
            if line.startswith("Output = "):
                break
            outjdl_file.write(line)

    for RootFiles in FailedJobRootFile:
        if DEBUG: print(RootFiles)
        bashCommand = "grep {} {}/*.stdout".format(RootFiles.replace(".root",""), CondorLogDir)
        if DEBUG: print(bashCommand)
        grep_stdout_files = os.popen(bashCommand).read()
        if DEBUG: print("~~"*51)
        if DEBUG: print(grep_stdout_files.strip())
        if DEBUG: print(len(grep_stdout_files))
        if DEBUG: print("~~"*51)
        OldRefFile = ""
        if grep_stdout_files.strip() != "":
            if DEBUG: print("==> ",grep_stdout_files.strip().split(':')[0].replace('.stdout',''))
            if grep_stdout_files.strip().split(':')[0].replace('.stdout','').split('_')[-2] == "resubmit":
                OldRefFile = grep_stdout_files.strip().split(':')[0].replace('.stdout','').split('_')[-4]
            else:
                OldRefFile = grep_stdout_files.strip().split(':')[0].replace('.stdout','').split('_')[-1]
        grepCommand_GetJdlInfo = 'grep -A1 -B3 "{}" {}'.format(RootFiles, InputJdlFile)
        if DEBUG: print(grepCommand_GetJdlInfo)
        grep_condor_jdl_part = os.popen(grepCommand_GetJdlInfo).read()
        if DEBUG: print("=="*51)
        if DEBUG: print(grep_condor_jdl_part)
        updateString = grep_condor_jdl_part.replace('$(Process)',OldRefFile+'_$(Process)'+ '_resubmit_' +Resubmit_no)
        if DEBUG: print("=="*51)
        if DEBUG: print(updateString)
        if DEBUG: print("=="*51)
        outjdl_file.write(updateString)
    outjdl_file.close()
    return outjdl_fileName

def main():
    parser = OptionParser()
    parser.add_option("-d", "--dir", dest="dir",default="",help="log directory path")
    parser.add_option("-s", "--stage-dest", dest="stage_dest",help="directory output files were staged to")
    parser.add_option("-i", "--input", dest="input",default="all_root.jdl",help="input jdl file with all root files present")
    parser.add_option("-r", "--resubmit", action="store_true",  dest="resubmit",default=False,help="resubmit")
    parser.add_option("-n", "--resubmit_no", dest="resubmit_no",default=1,help="resubmit counter")

    (options, args) = parser.parse_args()

    if options.stage_dest is not None:
        stageDir = os.path.abspath(options.stage_dest)
    else:
        stageDir = dir

    full_output =  get_files_from_jdl(options.input)
    print(full_output)
    present_output, present_output_WithPath =  list_root(stageDir)
    print(present_output)
    print("length(jdl file): {}".format(len(full_output)))
    print("Length(output root file): {}".format(len(present_output)))
    not_finished = list(set(full_output) - set(present_output))
    print(not_finished)
    corrupted_files = files_to_remove(present_output_WithPath,stageDir)
    not_finished += corrupted_files
    print(not_finished)
    print('Number of missing files : {}'.format(len(not_finished)))
    print('Missing the following files : {}'.format(not_finished))
    jdlfile = prepare_runJobs_missing(not_finished,options.input,options.dir,stageDir,str(options.resubmit_no))
    print(jdlfile)
    print('Submitting missing jobs : ')
    submit_missing(jdlfile,options.resubmit)

if __name__ == "__main__":
    main()
