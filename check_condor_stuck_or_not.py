import os

output = os.popen('condor_q -submitter rasharma').read()

error_check_string = 'Server responded with an error'

#Oprint output.split("\t")
lpcschedd = ""
print type(output)
for outputs in output.split('\n'):
  if outputs.find('Submitter') != -1:
    lpcschedd = outputs.split()[2].split('.')[0]
  if outputs.find('rasharma') != -1:
    condor_tail = "condor_tail "+outputs.split()[0]+" -name "+lpcschedd
    print "\n","-"*51,"\n\n"
    print outputs,"\n\n"
    print condor_tail
    print "\n"
    os.system(condor_tail)
print "\n\n"
