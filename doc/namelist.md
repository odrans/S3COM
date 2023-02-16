## S3COM namelist configuration

The default configuration namelist file is `config_default.nml`. Fields in bold (associated with `NA` as default values) must be edited before running S3COM.

#### General 

| Parameter Name | Description | Default |
| --- | --- | --- | 
| **fname_in** | **Input file name** | NA |
| path_out | Path for output files | "./" |
| suffix_out | suffix to be used in output files | "" |
| month | Month of the simulations | 1 |
| flag_retrievals | Turn on cloud retrievals | FALSE |
| flag_output_atm | Output atmospheric variables (S3COM_atm.nc) | FALSE |
| flag_output_jac | Output Jacobian variables (S3COM_jac.nc) | FALSE |
| npoints_it | Number of simultaneous data points sent to RTTOV | 100 |
| nchannels | Number of satellite channels | 2 |
| model_name | Name of the physical model to be used as input | "ICON" |

#### RTTOV instruments

Check the RTTOV website for the [list of selected instruments](https://nwp-saf.eumetsat.int/site/software/rttov/documentation/platforms-supported/) and  [channel specifications](https://nwp-saf.eumetsat.int/site/software/rttov/download/coefficients/spectral-response-functions/). Default is MODIS/Aqua.

| Parameter Name | Description | Default |
| --- | --- | --- | 
| platform | Platform ID | 9 |
| satellite | Satellite ID | 2 |
| instrument | Instrument ID | 13 |
| channel_list | List of channels to simulate (disabled if channel_seq is not 0) | 1, 2 |
| channel_seq | Sequence of channels to simulate (first and last values) | 0, 0 |

#### RTTOV initialization

See the [RTTOV general documentation](https://nwp-saf.eumetsat.int/site/download/documentation/rtm/docs_rttov13/users_guide_rttov13_v1.2.pdf) for more details on some of these parameters.

| Parameter Name | Description | Default |
| --- | --- | --- | 
| **path_rttov** | **Path to RTTOV** | NA |
| rttov_nthreads | Number of threads in RTTOV (MPI required) | 1 |
| do_jacobian_calc | | FALSE |
| do_opdep_calc | Include gaseous absorption | TRUE |
| ir_scatt_model | Infrared scattering model (1: DOM, 2: Chou-scaling) | 1 |
| vis_scatt_model | Visible scattering model (1: DOM, 2: Single-scattering, 3: MFASIS) | 1 |
| dom_nstreams | Number of streams used in DOM | 8 |
| dom_rayleigh | Activate Rayleigh multiple scattering | TRUE |
| ice_scheme | 1: Baum, 2: Baran 2014, 3: Baran 2018 | 3 |
| clw_scheme | 1: OPAC, 2: Deff | 2 |
| mmr_cldaer | cloud / aerosol units. true: kg/kg (cld+aer); false: g.m-3 (cld), cm-3 (aer) | FALSE |
| ozone_data | Origin of ozone data. True: user input, False: RTTOV database | FALSE |
| add_aerosols | Include aerosols | FALSE |
| add_clouds | Include clouds | TRUE |
| add_refrac | Include atmospheric refraction | TRUE |


