#!/bin/sh

# the argument to this script will be the file created by proxy command, i.e.
# voms-proxy-init --voms cms --valid 168:00
# Created proxy in /tmp/x509up_u48539.
# then the first argument will be x509up_u48539 only.

proxy=${1}
cp /tmp/${proxy} ~/
export X509_USER_PROXY=~/${proxy}
