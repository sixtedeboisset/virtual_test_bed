The results of this simulation are close to those of the first 1D simulation. Consequently, most of the plots are not reported here but the steady state values are compared together:

!table id=results_comparison_1D_3D caption= Compared results of the 1D and 3D simulations
| Loop                     | Parameter                        | 3D simulation                | 1D simulation           | Relative gap       |
| :----------------------- | :------------------------------- | :--------------------------- | :---------------------- | :----------------- |
| Primary loop             | T fluid - core min (K)           |                              | 888.697558              |                    |
|                          | T fluid - core max (K)           |                              | 1196.39737              |    |
|                          | mass flow rate (kg/s)            |                              | 9.35722261              |    |
|                          | p fluid - core out (MPa)         |                              | 9.004203                |    |
|                          | Hw core (avg) W/(K.m^2^)         |                              |                |    |
|                          | Power core (extracted) (MWth)    |                              | 15.000132               |    |
| Secondary loop           | T fluid - comp in  (K)           |                              | 267.897738              |    |
|                          | T fluid - comp out (K)           |                              | 419.690253              |    |
|                          | T fluid - turb in  (K)           |                              |                         |    |
|                          | T fluid - turb out (K)           |                              |                         |    |
|                          | p fluid - comp in  (MPa)         |                              |                         |    |
|                          | p fluid - comp out (MPa)         |                              |                         |    |
|                          | p fluid - turb in  (MPa)         |                              |                         |    |
|                          | p fluid - turb out (MPa)         |                              |                         |    |
|                          | mass flow rate (kg/s)            |                              |                         |    |
|                          | Power hx (extracted) (MWth)      |                              |                         |    |
|                          | Power generator (extracted) (MWe)|                              |                         |    |

The fluid temperature is presented in [3D_core_T_coolant] and [plot_core_T_coolant]. It increases progressively along the core, from 890 K to 1190 K at the core exit. The colors on the 3D representation are put on the solid materials but represent the temperatures in the surrounded coolant channels.

!media media/gcmr/3D_core_balance_of_plant/3D_core_T_coolant.png
      style=display: block;margin-left:auto;margin-right:auto;width:70%;
      id=3D_core_T_coolant
      caption= Fluid temperature profile in the core

!media media/gcmr/3D_core_balance_of_plant/plot_core_T_coolant.png
      style=display: block;margin-left:auto;margin-right:auto;width:70%;
      id=plot_core_T_coolant
      caption= Plot of the fluid temperature in the core


The moderator and fuel temperatures are presented in [3D_core_T] and [plot_core_T].

The maximum solid temperature stays under 1400 K, which respects the requirements. The temperature increases in the first half of the core and a maximum value is reached between the middle and the top of the core. This effect is due to the rise of the fluid temperature, which continues to increase after the middle of the core, where the power density is the highest: the fluid continues to receive power in the second half of the core.

!media media/gcmr/3D_core_balance_of_plant/3D_core_T.png
      style=display: block;margin-left:auto;margin-right:auto;width:70%;
      id=3D_core_T
      caption= Moderator and fuel temperature profile in the core

!media media/gcmr/3D_core_balance_of_plant/plot_core_T.png
      style=display: block;margin-left:auto;margin-right:auto;width:70%;
      id=plot_core_T
      caption= Plot of the moderator and fuel temperature in the core



The heat transfer coefficient is almost the same everywhere in the core. It is provided in [3D_core_htc] and [plot_core_htc]. The colors on the 3D representation are put on the solid materials but represent the heat transfer coefficient with the surrounded coolant channels.

!media media/gcmr/3D_core_balance_of_plant/3D_core_htc.png
      style=display: block;margin-left:auto;margin-right:auto;width:70%;
      id=3D_core_htc
      caption= Heat transfer coefficient profile in the core

!media media/gcmr/3D_core_balance_of_plant/plot_core_htc.png
      style=display: block;margin-left:auto;margin-right:auto;width:70%;
      id=plot_core_htc
      caption= Plot of the heat transfer coefficient in the core
