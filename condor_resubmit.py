"""Condor resubmit script.

# Three steps to submit condor script

1. Set input information
2. `grep` the **search string** ("preselected entries from root:") from the log files.
  1. If it is not present it means jobs failed.
  2. \todo check any other signature of the failed jobs and add it.
  3. Make list of failed jobs log file.
3. Get the initial 7 lines from the main jdl file on which original jobs was submitted.
  1. \todo number 7 in the previous statement is hard-coded. Make it flexible.
4. Add the above 7 lines into new jdl file.
5. To add the failed jobs information in the new jdl file, do following:
  1. From the log files check the name of root files that failed.
  2. Search the same root file from the main jdl file and grep the corresponding
    necessary lines from the main jdl file and add it to the new jdl file.
"""

import os

parser = argparse.ArgumentParser()
parser.add_argument('--LogFilePath', type=str,
                    default='condor_logs/Run2017_v7_5Aug20200/200805_032620/',
                    help='String to be added in the output file name')
parser.add_argument('--CondorFileName', type=str,
                    default='submit_condor_jobs_lnujj_Run2017_v7_5Aug20200',
                    help='String to be added in the output file name')
parser.add_argument('--ResubmitCount', type=int,
                    default=1,
                    help='String to be added in the output file name')
parser.add_argument('--StringSearch', type=str,
                    default='preselected entries from root:',
                    help='String to be added in the output file name')
args = parser.parse_args()


"""path of log file directory"""
path = args.LogFilePath

"""Name of main condor jdl/sh file name"""
condor_file_name = args.CondorFileName

"""This variable `Resubmit_no` is going to append in the new jdl file.
New jdl file name is the main jdl file + _resubmit_ + Resubmit_no

New JDL FILE name : condor_file_name + "_resubmit_" + Resubmit_no
"""
Resubmit_no = args.ResubmitCount

"""This string is going to be searched in the log files.

If this string is not present in the log file this means the jobs are failed.
"""
string_to_search = args.StringSearch

grepCommand = 'grep -L  "'+string_to_search+'" '+ path + os.sep +'*_1.stdout'
grepCommand = grepCommand.replace('//','/')
print 'grep command: ',grepCommand
output = os.popen(grepCommand).read()
print('output:')

outjdl_file = open(condor_file_name+'_resubmit_'+Resubmit_no+".jdl","w")
with open(condor_file_name+".jdl") as myfile:
    head = [next(myfile) for x in xrange(7)]

for lines in head:
  outjdl_file.write(lines)

for lines in output.split():
  print "==> ",lines.strip()
  # print "==> ",lines.strip().split('/')[-1].replace('.stdout','')
  if lines.strip().split('/')[-1].replace('.stdout','').split('_')[-2] == "resubmit":
    OldRefFile = lines.strip().split('/')[-1].replace('.stdout','').split('_')[-4]
  else:
    OldRefFile = lines.strip().split('/')[-1].replace('.stdout','').split('_')[-1]
  # print "====> ",OldRefFile
  grep_output = os.popen('grep -E "Running.*root" '+lines.strip()).read()
  # print grep_output
  root_file = grep_output.split()[5].split('/')[-1]
  grepCommand_GetJdlInfo = 'grep -A1 -B3 "'+root_file+'" '+condor_file_name+'.jdl'
  # print grepCommand_GetJdlInfo
  grep_condor_jdl_part = os.popen(grepCommand_GetJdlInfo).read()
  updateString = grep_condor_jdl_part.replace('$(Process)',OldRefFile+'_$(Process)'+ '_resubmit_' +Resubmit_no)
  # print updateString
  outjdl_file.write(updateString)
outjdl_file.close()

print "===> ",outjdl_file
