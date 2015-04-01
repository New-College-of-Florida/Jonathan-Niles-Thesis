#!/bin/sh

# make sure the system is up to date
apt-get update
apt-get install git libssl-dev openssl mysql-client-5.1 mysql-client-core-5.1

# get kentutils
git clone https://github.com/ENCODE-DCC/kentUtils.git
cd kentutils
make

# install all dependencies using python-pip
sudo apt-get install python-pip
pip install https://bitbucket.org/mirnylab/hiclib/get/tip.tar.gz
pip install https://bitbucket.org/mirnylab/mirnylib/get/tip.tar.gz
pip install bx-python
pip install pybedtools

# NOTE:
# You must run the following script root privileges (SUDO).
