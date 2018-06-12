# NonStatMoni
insert better description here

## Basic user guide:
* Configuration file: in order to run, a configuration file containing most of the options required for the execution is required. The file config_example.cfg is an example configuration file, with an explanation of most of the possible options.
* Running the first step: the first step consists in the acquisition of the main channel data and in the computation of its brms. Its results are saved in a .hdf5 file that can then be used by the second step. 
In order to launch  the first step one has to execute:
  * 
  ```sh
  python /path/to/NonStatMoni/Brms_generator/brms_main.py -i /path/to/config_file.ini -d /path/to/results/folder -b gpsb -e gpse -scn equal -sc V1:DQ_META_ITF_Mode -st -7 
```

** Notes: a results subfolder named "GPSB - GPSE - resolution xxxx - npoints xx, overlap xxxxxx" will be automatically created in the folder ```/path/to/results/folder``` . Where the gpsb, gpse and resolution, npoints and overlap parameters will be taken from the configuration file. This folder will contain a "brms.hdf5" hdf5 file containing the results of the computation. 

arguments of brms_main.py:
  -h, --help            show this help message and exit
  -i cfgFILE, --init cfgFILE
                        set initialization FILE
  -d OutDir, --dir OutDir
                        output directory
  -b gpsb, --gpsb gpsb  start gps (used if not using segment files)
  -e gpse, --gpse gpse  end gps (used if not using segment files)
  -sc State_Channel, --statechannel State_Channel
                        state channel that can be used to generate time
                        segments that will be analyzed
  -st State_Threshold, --threshold State_Threshold
                        Threshold that should be applied to the state channel
                        in order to generate segments
  -scn State_Condition, --statecondition State_Condition
                        How the state channel should be compared with the
                        threshold i.e. equal or lesser_equal
  -ss science_segment_file, --sciencesegments science_segment_file
                        Json formatted science segments file.If provided,
                        gpsb, gpse, statechannel threshold and state condition
                        arguments will be ignored
  -dqfl [DQ-Flag_files [DQ-Flag_files ...]], --dqflags [DQ-Flag_files [DQ-Flag_files ...]]
                        Json formatted data quality flag files. Periods where
                        any of these flags are active will be excluded from
                        analysis


set resfolder = "`find /users/valentin/TestConfigs/$1 -maxdepth 1 -type d -name '*'$gpsb'*' -print -quit`"
echo $resfolder
python /users/valentin/PycharmProjects/Nonstatmoni_reader/nonstatmoni_main.py -i /users/valentin/TestConfigs/$1.ini -d /users/valentin/TestConfigs/$1/report/ -b "$resfolder"/brms.hdf5 -p 5 --ntop 15 
  ```
* Second step: after the first step, NonStatMoni can use the processed BRMSs to attempt the estimation of low frequency coherences and correlation cohefficients between the BRMSs and slow auxiliary channels data. This step can use the same configuration file as the first step or, optionally, a different one, as long as the configuration file does not contain additional or different BRMS bands or channels with respect to the ones processed in the first step. Additionally, at the end of the computation of coherences and correlation cohefficients, the generation of the relevant plots and results page is also performed in this step. in order to execute this step one shall run
```sh
example_two
```
where x

usage: nonstatmoni_main.py [-h] [-i cfgFILE] [-b brms_file] [-d OutDir]
                           [-p NumberOfProcesses] [-n nTop]

optional arguments:
  -h, --help            show this help message and exit
  -i cfgFILE, --init cfgFILE
                        set initialization FILE
  -b brms_file, --brms brms_file
                        set brms data file
  -d OutDir, --dir OutDir
                        output directory
  -p NumberOfProcesses, --proc NumberOfProcesses
                        number of parallel processes to run
  -n nTop, --ntop nTop  number of top coherences and correlation cohefficients
                        to be shown in the result table




#!/bin/csh -f

echo "will use config "$1.ini

set gpsb = 1195594218
set gpse = 1195624218




echo "working..."
setenv PYTHONPATH /users/valentin/PycharmProjects/Nonstatmoni_reader #:$PYTHONPATH
python /users/valentin/PycharmProjects/Nonstatmoni_reader/Brms_generator/brms_main.py -i /users/valentin/TestConfigs/$1.ini -d /users/valentin/TestConfigs/$1 -b $gpsb -e $gpse -scn equal -sc V1:DQ_META_ITF_Mode -st -7 

set resfolder = "`find /users/valentin/TestConfigs/$1 -maxdepth 1 -type d -name '*'$gpsb'*' -print -quit`"
echo $resfolder
python /users/valentin/PycharmProjects/Nonstatmoni_reader/nonstatmoni_main.py -i /users/valentin/TestConfigs/$1.ini -d /users/valentin/TestConfigs/$1/report3/ -b "$resfolder"/brms.hdf5 -p 6 --ntop 15 

echo "done!"
 
