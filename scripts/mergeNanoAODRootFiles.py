#!/usr/bin/env python
import os
import glob
import ROOT
import subprocess
import sys
import logging
import argparse

logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')

def system(command):
    logging.debug("Executing command#12: {command}".format(command=command))
    #return subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)

def system_with_terminal_display(command):
    logging.info("Executing command#16: {command}".format(command=command))
    return subprocess.call(command, shell=True)

def isValidRootFile(fname):
    if not os.path.exists(fname):
        return False
    f = ROOT.TFile.Open(fname)
    if not f:
        return False
    isValid = not (f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered) or f.GetListOfKeys().IsEmpty())
    f.Close()
    if not isValid:
        logging.warning("Zombie found: {}".format(fname))
    return isValid

def checkfaulty(fname, ref = None):
    if not ref:
        ref = ROOT.TFile.Open(fname)
    faultyfiles = []  # Initialize an empty list to hold the names of faulty files
    probe = ROOT.TFile.Open(fname)

    if not probe:
        logging.error("Could not open file {fname}".format(fname=fname))
        return False

    for e in ref.GetListOfKeys():
        name = e.GetName()
        try:
            otherObj = probe.GetListOfKeys().FindObject(name).ReadObj()
        except:
            faultyfiles.append(probe.GetName())

    probe.Close()

    if faultyfiles:
        logging.warning("Faulty files found: {faultyfiles}".format(faultyfiles=faultyfiles))
        return False

    return True

def isValidAndFaultFree(fname,ref = None):
    # First check if it's a valid ROOT file
    if not isValidRootFile(fname):
        return False

    # Then check for faulty keys
    return checkfaulty(fname, ref)

def merge_files(targetFile, filesToMerge):
    # If there are more than 100 files, split them into groups of 100
    if len(filesToMerge) > 100:
        logging.info('A lot of files to merge (>100); this might take some time...')
        tempTargets = []
        tempFilesToMerge = [filesToMerge[x:x+100] for x in range(0, len(filesToMerge), 100)]

        temp_directory = os.path.join("/tmp", os.environ['USER'])
        logging.info("temp_directory: {temp_directory}".format(temp_directory=temp_directory))

        for i, batch in enumerate(tempFilesToMerge):
            tempTargetFile = os.path.join(temp_directory, targetFile.replace('.root', '-temp{}.root'.format(i)))
            tempTargets.append(tempTargetFile)
            logging.info("tempTargetFile: {tempTargetFile}".format(tempTargetFile=tempTargetFile))
            logging.debug("batch: {batch}".format(batch=batch))

            # Check if temporary target file already exists and is valid
            if os.path.exists(tempTargetFile):
                if isValidAndFaultFree(tempTargetFile):
                    continue
                else:
                    logging.info("Removing temp hadd file {tempTargetFile}".format(tempTargetFile=tempTargetFile))
                    system_with_terminal_display('rm {tempTargetFile}'.format(tempTargetFile=tempTargetFile))

            system_with_terminal_display('haddnano.py {0} {1}'.format(tempTargetFile, ' '.join(batch)))


        # Merge temporary files into the final target file
        targetFile = os.path.join(temp_directory, targetFile)
        system_with_terminal_display('haddnano.py {0} {1}'.format(targetFile, " ".join(tempTargets)))

        # Remove temporary files
        for tempTarget in tempTargets:
            system_with_terminal_display('rm {tempTarget}'.format(tempTarget=tempTarget))

    # If there are 100 or fewer files, directly merge them
    else:
        logging.info("Merging {0} files into {1}".format(len(filesToMerge), targetFile))
        logging.debug("haddnano.py {targetFile} {filesToMerge}".format(targetFile=targetFile, filesToMerge=' '.join(filesToMerge)))
        system_with_terminal_display('haddnano.py {targetFile} {filesToMerge}'.format(targetFile=targetFile, filesToMerge=' '.join(filesToMerge)))


def main():
    parser = argparse.ArgumentParser(description="Merge ROOT files using haddnano.py.")
    parser.add_argument('-i', '--inputDir', type=str, help="Path of the input directory that contains ROOT files to be merged.")
    parser.add_argument('-o', '--outputDir', type=str, help="Path of the output directory where the merged ROOT file will be saved.")
    parser.add_argument('-f', '--outputFile', type=str, help="Name of the hadd-ed output ROOT file.")

    args = parser.parse_args()

    # Validate and use the argument
    inputDir = args.inputDir
    targetFile = os.path.join(args.outputDir, args.outputFile)
    if not os.path.isdir(inputDir):
        logging.error("The specified input directory does not exist: {inputDir}".format(inputDir=inputDir))
        exit(1)

    logging.info("Input directory: {inputDir}".format(inputDir=inputDir))

    #  List all ROOT files in the specified directory
    filesToMerge = glob.glob(os.path.join(inputDir, '*.root'))
    if not filesToMerge:
        logging.error("No ROOT files found in the specified directory.")
        exit(1)

    # Here you would continue with your logic to merge the files as needed
    logging.info("Found {filesToMerge} ROOT files to merge.".format(filesToMerge=len(filesToMerge)))
    logging.debug("filesToMerge: {0}".format(filesToMerge))
    logging.info("targetFile: {0}".format(targetFile))

    # Check if the target file already exists and is valid
    if os.path.exists(targetFile):
        if isValidAndFaultFree(targetFile):
            logging.info("Seems hadd is already performed. Skipping.")
            exit(0)
        else:
            logging.info("Removing invalid target file {targetFile}".format(targetFile=targetFile))
            system_with_terminal_display('rm {targetFile}'.format(targetFile=targetFile))

    # Merge the files
    filesToMerge = [f for f in filesToMerge if isValidAndFaultFree(f, ROOT.TFile.Open(filesToMerge[0]))]
    logging.info("After zombie: filesToMerge: {0}".format(len(filesToMerge)))
    merge_files(targetFile, filesToMerge)

if __name__ == "__main__":
    main()
