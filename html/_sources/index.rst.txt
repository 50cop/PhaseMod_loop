.. Doc_231214_Soliton_manipulation documentation master file, created by
   sphinx-quickstart on Wed Jan 10 10:07:12 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to the doc *Processing C2--00002 from 14/12/2023*
=========================================================
  
Preamble
--------

This documentation describes some codes used to process data recorded in the PhaseMod recirculating loop the 14/12/2023.
The file structure is as follows:
::

    231214_Soliton_manipulation
    ├── Recordings_browser.mlapp
    ├── Processing of C2_00002_2_Solitons_oscillations          
    │   ├── Generated figures
    │   │   └── ...
    │   ├── Interactive_Hamiltonian_trajectory_calculator
    │   │   └── ...
    │   ├── Processing_two_solitons_oscillations.m
    │   ├── two_solitons_oscillations.mat
    │   └── ...
    ├── subscripts
        ├── pltfcn_EXP
        │   └── ...
        └── ...         
    ├── gif
    └── Doc_231214_Soliton_manipulation          
    

The first stages of the processing of recorded data are not detailed here. There are as follows:

#. Recorded file is openned a first time to perform calibration via the app `Recordings_browser.mlapp` which calls scripts contained in the `subscripts` folder;
#. `Recordings_browser.mlapp` is called again with the correct calibration parameters;
#. Spatio-temporal diagrams of the relevant experimental runs are exported to matlab structures saved in the `two_solitons_oscillations.mat` file;
    
Of this series of experiments, the *C2--00002.trc* recording is the most interesting and has undergone the most thorough processing.



Codes description
-----------------

Below are links to the description of some of the codes used for processing the spatio-temporal diagrams extracted from the recording.


.. toctree::
    :maxdepth: 1
       
    Processing_two_solitons_oscillations
    Interactive_Hamiltonian_trajectory_calculator








