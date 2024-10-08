## Inputs 
## Overview of the thermodynamic variables required from the ICON-LES simulations for RTTOV

|   | Name       | Description                                                         | Units            |
|---|------------|---------------------------------------------------------------------|------------------|
| **1D Variables** |            |                                                             |                  |
| 1 | time       | Sample time (proleptic gregorian) UTC                                | %Y%m%d.%f        |
|   |            |                                                                     |                  |
| **2D (lat x lon)** |            |                                                             |                  |
| 1 | FR_LAND    | Land use class fraction                                              | 0 or 1           |
| 2 | topography_c | Geometric height of the Earth's surface above sea level              | m                |
| 3 | ps         | Surface pressure                                                     | Pa               |
| 4 | t_s        | Skin temperature                                                     | K                |
| 5 | tas        | Temperature at 2m                                                    | K                |
| 6 | huss       | Specific water vapour content at 2m                                  | kg kg⁻¹          |
| 7 | u_10m      | Zonal wind at 10m                                                    | m s⁻¹            |
| 8 | v_10m      | Meridional wind at 10m                                               | m s⁻¹            |
|   |            |                                                                     |                  |
| **3D (150 full level center: height x lat x lon / 151 half level center: height_2 x lat x lon)** |   |              |                  |
| 1 | z_mc       | Geometric height at full level center                                | m                |
| 2 | z_ifc      | Geometric height at half level center                                | m                |
| 3 | ta         | Air temperature at full level center                                 | K                |
| 4 | pres       | Air pressure at full level center                                    | Pa               |
| 5 | hus        | Specific humidity at full level center                               | kg kg⁻¹          |
| 6 | clc        | Cloud cover at full level center (0 - 1)                             | fraction         |
| 7 | clw        | Specific cloud water content at full level center                    | kg kg⁻¹          |
| 8 | qnc        | Cloud droplet number concentration at full level center              | kg⁻¹             |
| 9 | cli        | Specific cloud ice content at full level center                      | kg kg⁻¹          |
| 10 | qr        | Rain mixing ratio at full level center                               | kg kg⁻¹          |
| 11 | qs        | Snow mixing ratio at full level center                               | kg kg⁻¹          |
|   |            |                                                                     |                  |
| **3D variables calculated (151 half level center: height_2 x lat x lon)** |  |              |                  |
| 1 | ta_ifc     | Air temperature at half level center                                 | K                |
| 2 | pres_ifc   | Pressures at half level center                                       | Pa               |
| 3 | hus_ifc    | Specific humidity at half level center                               | kg kg⁻¹          |


*TOA: Top Of Atmosphere correspond to the Height 0

