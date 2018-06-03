# NonStatMoni
insert better description here

##Basic user guide:
* Configuration file: in order to run, a configuration file containing most of the options required for the execution is required. The file config_example.cfg is an example configuration file, with an explanation of most of the possible options.
* Running the first step: the first step consists in the acquisition of the main channel data and in the computation of its brms. Its results are saved in a .hdf5 file that can then be used by the second step. 
In order to launch  the first step one has to execute:
  * 
  ```sh
  example code .py
  ```
* Second step: after the first step, NonStatMoni can use the processed BRMSs to attempt the estimation of low frequency coherences and correlation cohefficients between the BRMSs and slow auxiliary channels data. This step can use the same configuration file as the first step or, optionally, a different one, as long as the configuration file does not contain additional or different BRMS bands or channels with respect to the ones processed in the first step. Additionally, at the end of the computation of coherences and correlation cohefficients, the generation of the relevant plots and results page is also performed in this step. in order to execute this step one shall run
```sh
example_two
```
where x
