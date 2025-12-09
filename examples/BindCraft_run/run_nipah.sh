#!/bin/bash

# this should be run in the BindCraft dir!
python -u ./bindcraft.py \
	--settings 'BindCraft_example/nipah.json' \
	--filters  './settings_filters/default_filters.json' \
	--advanced './settings_advanced/default_4stage_multimer.json'
