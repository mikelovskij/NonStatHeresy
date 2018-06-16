# NonStatMoni
insert better description here

## Basic user guide:
* **Configuration file**: in order to run, a configuration file containing most of the options required for the execution is required. The file config_example.ini is an example configuration file, with an explanation of most of the possible options.

* **Running the first step**: the first step consists in the acquisition of the main channel data and in the computation of its brms. Its results are saved in a .hdf5 file that can then be used in the second step. The results of the first step will be placed in a results subfolder named "GPSB - GPSE - resolution xxxx - npoints xx, overlap xxxxxx" that will be automatically created in a folder set by the user.

   * **How to launch the first step** The file brms_main.py has to be executed and requires some arguments, as in the example below:
  ```sh
  python /path/to/NonStatMoni/Brms_generator/brms_main.py -i /path/to/config_file.ini -d /path/to/results/folder -b gpsb -e gpse -scn equal -sc V1:DQ_META_ITF_Mode -st -7 
  ```

  * The full list of arguments that can be set for the first step during the call of brms_main.py is listed below. Wile the -i , -d, -b, -e arguments are all mandatory, one has to either chose a state channel (and optionally a state threshold and state condition) **OR** a science segment json file and a list of DQ_flags files. 
 In the example above, ```-scn equal -sc V1:DQ_META_ITF_Mode -st -7 ``` mean that only periods when the value of the V1:DQ_META_ITF_Mode channel is equal to -7 will be used. 

  ```
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

  ```
```
set resfolder = "`find /users/valentin/TestConfigs/$1 -maxdepth 1 -type d -name '*'$gpsb'*' -print -quit`"
echo $resfolder
 
```
  
* **Second step**: once the results of the first step have been produced, NonStatMoni uses the processed BRMSs to attempt the estimation of low frequency coherences and correlation cohefficients between the BRMSs and slow auxiliary channels data. This step can use the same configuration file as the first step or, optionally, a different one, as long as the configuration file does not contain additional or different BRMS bands or channels with respect to the ones processed in the first step.
In order to speed up the operations it is possible to use multithreading, by choosing the number of processes with the "-p" argument.
Additionally, at the end of the computation of coherences and correlation cohefficients, the generation of the relevant plots and results page is also performed in this step. in order to execute this step one shall run

Example command used to launch the second step:
```sh
 python /path/to/Nonstatmoni/nonstatmoni_main.py -i /path/to/Config.ini -d /path/to/output/results-b /path/to/brms.hdf5 -p 5 --ntop 15
 ```

 ```
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
 ```



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
 
