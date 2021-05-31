#!/bin/bash

# Usage: create_venv.sh <virtualenv> <python_path>

set -e

if [[  -z "$2" ]]; then
	    PYTHON=`which python3`
    else
	        PYTHON=$2
	fi

	VENV=`readlink -f ${1-$(basename $(pwd))-env}`

	if [[ ! -d $VENV ]]; then
		    $PYTHON -m venv $VENV
	    else
		        echo "$VENV already exists"
			    exit 0
		    fi

		    source $VENV/bin/activate
		    pip install --upgrade pip
		    pip install -r requirements.txt

#install art
wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
tar -xzvf artbinmountrainier2016.06.05linux64.tgz
cp art_bin_MountRainier/art_illumina $VENV/bin
cp art_bin_MountRainier/ART_profiler_illumina/* $VENV/bin
chmod u+x $VENV/bin
rm -fR art_bin_MountRainier artbinmountrainier2016.06.05linux64.tgz
