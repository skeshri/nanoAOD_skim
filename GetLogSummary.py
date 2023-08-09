import os
import re

DEBUG = False

# Define the path to the directory containing the *.stdout files
log_files_dir = 'condor_logs/Run2018_v9/230808_225816'

# Regular expression patterns to match the desired lines
patterns = {
    "Total": re.compile(r'Total:\s+:(\d+)\s+Events'),
    "PassTrig": re.compile(r'PassTrig:\s+:(\d+)\s+Events'),
    "Pass2eCut": re.compile(r'Pass2eCut:\s+:(\d+)\s+Events'),
    "Pass2muCut": re.compile(r'Pass2muCut:\s+:(\d+)\s+Events'),
    "Pass2lCut": re.compile(r'Pass2lCut:\s+:(\d+)\s+Events'),
    "Pass2l1JCut": re.compile(r'Pass2l1JCut:\s+:(\d+)\s+Events'),
    "Pass2l2jCut": re.compile(r'Pass2l2jCut:\s+:(\d+)\s+Events'),
    "Pass2l1Jor2jCut": re.compile(r'Pass2l1Jor2jCut:\s+:(\d+)\s+Events')
}

# Dictionary to store the sum of events for each category
category_sums = {}

# Loop through each *.stdout file in the directory
for filename in os.listdir(log_files_dir):
    if filename.endswith('.stdout'):
        file_path = os.path.join(log_files_dir, filename)
        if DEBUG: print("\n\nFilename: "+log_files_dir+'/'+filename)

        with open(file_path, 'r') as f:
            content = f.read()

            # Find events for each category based on the specified patterns
            for category, pattern in patterns.items():
                match = pattern.search(content)
                if DEBUG: print("{} {}".format(category, pattern))
                if DEBUG: print("Match: {}".format(match))
                if match:
                    events = int(match.group(1))
                    if DEBUG: print("{:25} {}".format(category, events))
                    if category in category_sums:
                        category_sums[category] += events
                    else:
                        category_sums[category] = events

# Print the sum of events for each category
if DEBUG: print("\n\n")
for category, events in category_sums.items():
    print(f"{category:17}: {events:>7} Events")
