#!/bin/bash

conda activate crispro

wget https://gitlab.com/bauerlab/crispro/-/archive/main/crispro-main.zip
unzip crispro-main.zip
cd crispro-main
python setup.py install
