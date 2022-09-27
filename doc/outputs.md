## S3COM outputs

Up to three files are created by S3COM:
- `S3COM_[_prefix_]rad.nc` for radiation outputs
- `S3COM_[_prefix_]atm.nc` for atmospheric outputs (optional)
- `S3COM_[_prefix_]ret.nc` for retrieval outputs (optional)

The prefix can optionally be added from the namelist file. Default is empty.


### Dimensions

| Name | Description | Units |
| --- | --- | --- | 
| Latitude | Latitude| Degrees North | 
| Longitude | Longitude | Degrees East |
| height | model level index | - |
| Channel | RTTOV channel index | - |


### Radiation parameters

| Name | Description | Units |
| --- | --- | --- | 
| BRF_total | All-sky upwelling reflectance at TOA* | - | 
| BRF_clear | Clear-sky upwelling reflectance at TOA | - |
| BT_total | All-sky upwelling brightness temperature at TOA | K |
| BT_clear | Clear-sky upwelling brightness temperature at TOA | K |
| Radiance_total | All-sky upwelling radiance at TOA | $\mathrm{W}$ $\mathrm{m}^{-2}$ $\mathrm{sr}^{-1}$ $\mathrm{um}^{-1}$  | 
| Radiance_clear | Clear-sky upwelling radiance at TOA | $\mathrm{W}$ $\mathrm{m}^{-2}$ $\mathrm{sr}^{-1}$ $\mathrm{um}^{-1}$ |
| BRDF | Surface bidirectionnal reflectance density function (prescribed by RTTOV) | - |

*TOA: Top Of Atmosphere

### Atmospheric parameters

| Name | Description | Units |
| --- | --- | --- | 
| ta | Atmospheric Temperature | K | 
| z | Height above sea level | m |
| clc | Cloud fraction | - |
| cdnc | Cloud droplet number concentration | $\mathrm{m}^{-3}$ |
| lwc | Liquid water content | $\mathrm{kg}$ $\mathrm{m}^{-3}$ | 
| reff | Droplet effective radius| $\mu \mathrm{m}$ |


### Atmospheric parameters (To be Completed)

| Name | Description | Units |
| --- | --- | --- | 
| iwp_ret | Retrieved ice water path| kg $\mathrm{m}^{-2}$ | 
