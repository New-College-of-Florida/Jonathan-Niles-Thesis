#!/bin/sh

# install all dependencies using python-pip
sudo apt-get install python-pip
pip install https://bitbucket.org/mirnylab/hiclib/get/tip.tar.gz
pip install https://bitbucket.org/mirnylab/mirnylib/get/tip.tar.gz
pip install bx-python
pip install pybedtools
