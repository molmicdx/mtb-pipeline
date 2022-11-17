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

