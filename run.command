#!/bin/bash
cd ~/CSTB
export FLASK_APP=scripts/webapp.py
export FLASK_DEBUG=1
flask run
