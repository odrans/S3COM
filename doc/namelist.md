## S3COM namelist configuration

The default configuration namelist file is `config_default.nml`. Fields with `NA` as default values must be edited.

#### General 

| Parameter Name | Description | Default |
| --- | --- | --- | 
| fname_in | Input file name | NA |
| path_out | Path for output files | "./" |
| suffix_out | suffix to be used in output files | "" |
| month | Month of the simulations | 1 |
| flag_retrievals | Turn on cloud retrievals (otherwise only simulate measurements)  | FALSE |
| flag_output_atm | Output atmospheric variables (S3COM_atm.nc) | FALSE |
| npoints_it | Number of simultaneous data points sent to RTTOV | 100 |
| nchannels | Number of satellite channels | 2 |

#### RTTOV initialization

| Parameter Name | Description | Default |
| --- | --- | --- | 
| path_rttov | Path to RTTOV | NA | 
| addrefrac | Atmospheric refraction | TRUE |
| ir_scatt_model | Thermal infrared scattering model | 1 |
| vis_scatt_model | Visible scattering model | 1 |
| dom_nstreams | Number of streams used in DOM | 8 |
| dom_rayleigh | Rayleigh multiple scattering | TRUE |
| rttov_nthreads | Number of threads in RTTOV (MPI required) | 1 |

#### RTTOV instruments

Check the RTTOV website for the [list of selected instruments](https://nwp-saf.eumetsat.int/site/software/rttov/documentation/platforms-supported/) and  [channel specifications](https://nwp-saf.eumetsat.int/site/software/rttov/download/coefficients/spectral-response-functions/). Default is MODIS/Aqua.

| Parameter Name | Description | Default |
| --- | --- | --- | 
| platform | Platform ID | 9 |
| satellite | Satellite ID | 2 |
| instrument | Instrument ID | 13 |
| channel_list | List of channels to simulate | 1, 2 |
