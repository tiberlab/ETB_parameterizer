#!/bin/bash
# snapshot.sh - a script to save a snapshot of the current state of
#               the ETB fit process
#
# Usage:
#   ./snapshot.sh <profile_name> <stage_number> <series_number> <save_number>
#
#   e.g. ./snapshot.sh GaN 2 1 01

# check the number of arguments
if [ $# -ne 4 ]
then
    echo "Usage: ./snapshot.sh <profile_name> <stage_number> <series_number> <save_number>"
    exit 1
fi

# check the validity of the arguments
profile_name=$1
stage_number=$2
series_number=$3
save_number=$4

if [[ ! $profile_name =~ ^[A-Za-z0-9_]+$ ]]
then
    echo "Warning: profile_name contains invalid characters. Only alphanumeric and underscore are allowed."
fi

if [[ ! $stage_number =~ ^[0-9]+$ ]] || [ $stage_number -lt 1 ] || [ $stage_number -gt 4 ]
then
    echo "Warning: stage_number should be a positive integer between 1 and 4."
fi

if [[ ! $series_number =~ ^[0-9]+$ ]] || [ $series_number -lt 0 ]
then
    echo "Warning: series_number should be a non-negative integer."
fi

if [[ ! $save_number =~ ^[0-9]+$ ]] || [ $save_number -lt 0 ]
then
    echo "Warning: save_number should be a non-negative integer."
fi

# get the current date
date=$(date +%d%b%Y)

# check if the folder exists or not, if not then create it
if [ ! -d "../archives/$1/Stage$2/series$3_save$4_$date" ]
then
    mkdir ../archives/$1/Stage$2/series$3_save$4_$date
    # if the folder exists already then ask the user if they want to overwrite it
else
    echo "The folder ../archives/$1/Stage$2/series$3_save$4_$date already exists. Do you want to overwrite it? (y/n)"
    read overwrite
    if [ $overwrite != "y" ]
    then
        exit 1
    fi
fi

# copy Stage2.m, log.csv, rep-point.csv and rep-cost.csv to the archive directory
# get the lines in the `if true ... end` of Stage2.m, and write them to a new file
start_line=$(grep -n "if true" Stage$2.m | cut -d: -f1)
end_line=$(grep -n "end%" Stage$2.m | cut -d: -f1)
sed -n "${start_line},${end_line}p" Stage$2.m > ../archives/$1/Stage$2/series$3_save$4_$date/profiles.m
cp ../aux/log_stage$2.csv ../archives/$1/Stage$2/series$3_save$4_$date
cp ../aux/rep-point_stage$2.csv ../archives/$1/Stage$2/series$3_save$4_$date
cp ../aux/rep-costs_stage$2.csv ../archives/$1/Stage$2/series$3_save$4_$date
