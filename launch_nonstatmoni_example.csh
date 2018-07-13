#!/bin/csh -f

echo "will use config file "$1.ini

set gpsb = 1195594218
set gpse = 1195624218
set config_file_fldr = /path/to/config_file_folder
set nonstatmoni_fldr = /path/to/Nonstatmoni

setenv PYTHONPATH $nonstatmoni_fldr:$PYTHONPATH


echo "starting the generation of the BRMSs"
python /path/to/Nonstatmoni/Brms_generator/brms_main.py -i $config_file_fldr/$1.ini -d $config_file_fldr/$1 -b $gpsb -e $gpse -scn equal -sc V1:DQ_META_ITF_Mode -st -7


echo "Looking for the brms results folder..."
set resfolder = "`find $config_file_fldr/$1 -maxdepth 1 -type d -name '*'$gpsb'*' -print -quit`"
echo "Brms results found in "
echo $resfolder
echo "Launching the slow Correlation and Coherence computation"
python /path/to/Nonstatmoni_reader/nonstatmoni_main.py -i $config_file_fldr/$1.ini -d $config_file_folder/$1/report/ -b "$resfolder"/brms.hdf5 -p 6 --ntop 15

echo "done!"
