#!/bin/sh

python3 -m pip install --user --upgrade pip

python3 -m pip install --user --upgrade virtualenv

python3 -m virtualenv -p `which python3` env

source ./env/bin/activate

python3 -m pip install --upgrade -r requirements.txt

