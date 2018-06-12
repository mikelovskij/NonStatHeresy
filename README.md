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
Notes: a results folder named INSERISCONOME will be created in the folder ```/path/to/results/folder``` . 
parameters 


set resfolder = "`find /users/valentin/TestConfigs/$1 -maxdepth 1 -type d -name '*'$gpsb'*' -print -quit`"
echo $resfolder
python /users/valentin/PycharmProjects/Nonstatmoni_reader/nonstatmoni_main.py -i /users/valentin/TestConfigs/$1.ini -d /users/valentin/TestConfigs/$1/report/ -b "$resfolder"/brms.hdf5 -p 5 --ntop 15 
  ```
* Second step: after the first step, NonStatMoni can use the processed BRMSs to attempt the estimation of low frequency coherences and correlation cohefficients between the BRMSs and slow auxiliary channels data. This step can use the same configuration file as the first step or, optionally, a different one, as long as the configuration file does not contain additional or different BRMS bands or channels with respect to the ones processed in the first step. Additionally, at the end of the computation of coherences and correlation cohefficients, the generation of the relevant plots and results page is also performed in this step. in order to execute this step one shall run
```sh
example_two
```
where x
