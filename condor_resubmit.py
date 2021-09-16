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
import argparse

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
parser.add_argument("--ErrorStrings", nargs="+",
                    default=["Disk quota exceeded", "error writing"],
                    help='''Possible error messages using which one can
                    identify the failed jobs''')
parser.add_argument('--StringSearch', type=str,
                    default='Selected entries from root:',
                    help='''If this string is there in the log file then
                    we interpret that this job is success.
                    ''')

# subparser = parser.add_subparsers(dest='command')
# ErrorSearch   = subparser.add_parser('ErrorSearch')
# SuccessSearch = subparser.add_parser('SuccessSearch')

# ErrorSearch.add_argument("--ErrorStrings", nargs="+",
#                          default=["Disk quota exceeded", "error writing"],
#                          help='''Possible error messages using which one can
#                               identify the failed jobs
#                               ''')
# SuccessSearch.add_argument('--StringSearch', type=str,
#                            default='Sected entries from root:',
#                            help='''If this string is there in the log file then
#                            we interpret that this job is success.
#                            ''')

args = parser.parse_args()


"""Name of main condor jdl/sh file name"""
condor_file_name = args.CondorFileName
if condor_file_name.endswith('.jdl'):
    condor_file_name = condor_file_name.replace('.jdl', '')
if condor_file_name.endswith('.sh'):
    condor_file_name = condor_file_name.replace('.sh', '')

"""This variable `Resubmit_no` is going to append in the new jdl file.
New jdl file name is the main jdl file + _resubmit_ + Resubmit_no

New JDL FILE name : condor_file_name + "_resubmit_" + Resubmit_no
"""
Resubmit_no = args.ResubmitCount

"""This string is going to be searched in the log files.

If this string is not present in the log file this means the jobs are failed.
"""
# string_to_search = args.StringSearch
string_to_search = '\|'.join(args.ErrorStrings)

# grepCommand = 'grep -L  "' + string_to_search + '" ' +  args.LogFilePath + os.sep + '*.stdout'  # -L for not listing files that contatins string
grepCommand = 'grep  -l "' + string_to_search + '" ' +  args.LogFilePath + os.sep + '*.stdout' # -l to just show filenames
grepCommand = grepCommand.replace('//', '/')
print('grep command: {}'.format(grepCommand))
output = os.popen(grepCommand).read()
# print('output: \n{}\n{}{}'.format("="*121, output, "="*121))

outjdl_file = open(condor_file_name + '_resubmit_' + str(Resubmit_no) + ".jdl", "w")
with open(condor_file_name + ".jdl") as myfile:
    head = [next(myfile) for x in xrange(7)]

for lines in head:
    outjdl_file.write(lines)

countLines = 0
for lines in output.split():
    print("==="*108)
    print("==> Reading File: {}".format(lines.strip()))
    FileName = lines.strip()
    print("Length: {}".format(len(FileName)))
    print("==> File name: {}".format(FileName))
    SplitArray = FileName.split('/')
    print("==> Splitted array based on '/': {}".format(SplitArray))
    FileNameWithoutExt = SplitArray[-1].replace('.stdout', '')
    print("==> File name: {}".format(FileNameWithoutExt))
    if FileNameWithoutExt.split('_')[-2] == "resubmit":
        OldRefFile = FileNameWithoutExt.split('_')[-4]
    else:
        OldRefFile = FileNameWithoutExt.split('_')[-1]
    print "====> OldRefFile: ", OldRefFile
    grep_output = os.popen('grep -E "Filename.*root" ' + FileName).read()
    # print "====> grep_output: ", grep_output
    # print "====> grep_output: ", grep_output.split()
    # print "====> grep_output: ", grep_output.split()
    # root_file = grep_output.split()[5].split('/')[-1]
    rootFileNameWithoutPath = grep_output.split()[1].split('/')[-1]
    grepCommand_GetJdlInfo = 'grep -A1 -B3 "' + str(rootFileNameWithoutPath) + '" ' + str(condor_file_name) + '.jdl'
    # print "===***> ", grepCommand_GetJdlInfo
    grep_condor_jdl_part = os.popen(grepCommand_GetJdlInfo).read()
    # print "===***> ", grep_condor_jdl_part
    updateString = grep_condor_jdl_part.replace('$(Process)', str(OldRefFile) + '_$(Process)' + '_resubmit_' + str(Resubmit_no))
    print updateString
    outjdl_file.write(updateString)
    countLines += 1
outjdl_file.close()

print "===> ", outjdl_file
print "===> total files read: ", countLines
