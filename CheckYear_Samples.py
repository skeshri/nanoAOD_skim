"""This file contains list of functions that helps
us to get the year using the dataset name.
"""


def CheckMCYear(MCSampleName):
    """
    Check Year of MC based on campaign string search
    in the path of root file (or DAS full name)

    Args:
        MCSampleName (str): file name or DAS name

    Returns:
        int: year with which the sample belongs
    """
    if MCSampleName.find("RunIIAutumn18NanoAODv") != -1:
        year = 2018
    if MCSampleName.find("RunIIFall17NanoAODv") != -1:
        year = 2017
    if MCSampleName.find("RunIISummer16NanoAODv") != -1:
        year = 2016
    return year


def CheckDataYear(DataSampleName):
    """
    Check Year of Data based on campaign string search
    in the path of root file (or DAS full name)

    Args:
        DataSampleName (str): file name or DAS name

    Returns:
        int, str: returns year with which the data belongs along with the corresponding
                  json file name with relative path
    """
    if DataSampleName.find("Run2016") != -1:
        Year = 2016
        jsonFileName = "data/jsonFiles/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt"
    if DataSampleName.find("Run2017") != -1:
        Year = 2017
        jsonFileName = "data/jsonFiles/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt"
    if DataSampleName.find("Run2018") != -1:
        Year = 2018
        jsonFileName = "data/jsonFiles/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt"
    return Year, jsonFileName


def CheckYear(MCDataSampleName="", ifMC=False):
    """
    Check Year of Data based on campaign string search
    in the path of root file (or DAS full name)

    Args:
        MCDataSampleName (str, optional): file name or DAS name
        ifMC (bool, optional): True if `MCDataSampleName` belongs to MC else False

    Returns:
        TYPE: returns year with which the data belongs
    """
    if ifMC:
        return CheckMCYear(MCDataSampleName)
    else:
        return CheckDataYear(MCDataSampleName)
