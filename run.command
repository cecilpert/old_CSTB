#!/bin/bash
export ROOTFOLDER=""
export FLASK_APP=scripts/webapp.py
export FLASK_DEBUG=1
#flask run --host=crispr.arwen-dev.ibcp.fr:3002
/data/software/mobi/flask/0.12.2/bin/flask run -p 3002