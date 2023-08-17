# Model description

## Mesh

A complete `Mesh` representing each assembly which constitutes the core pattern is built using the [Mesh](https://mooseframework.inl.gov/syntax/Mesh/index.html) and [MeshGenerator](https://mooseframework.inl.gov/syntax/Mesh/index.html)  MOOSE tools. The input file is `core_mesh.i` and the mesh is in the file `core_mesh_in.e`.

## Transfer - general principle

The general idea is to transfer data between one main app, which solves the heat conduction 3D problem of the core, and a sub app, which solves the 1D thermal hydraulic problem of the loops. These data are about the wall and the fluid temperatures and the heat transfer coefficient (the wall temperature is the temperature of the boundary between the fluid and the graphite in the coolant channels).

This principle is presented in [transfer_principle]:

!media media/gcmr/3D_core_balance_of_plant/transfer_principle.png
      style=display: block;margin-left:auto;margin-right:auto;width:25%;
      id=transfer_principle
      caption= Principle of the transfer

## THM sub app

The input file built for the thermal hydraulic simulation of the two loops is reused for the following simulation. It is nevertheless simplified on the secondary side. In this loop, only the heat exchanger is considered with imposed inlet (temperature, mass flow rate) and outlet (pressure) operating conditions. These ones are copied on the results of the 1D simulation (steady state conditions).

Some adjustments of the pump have also been done to reach better operating values in the core, but the operating parameters are still relatively close to those chosen in the 1D thermalhydraulic model.

The core is also different to be able to couple it with the heat conduction model.

A smaller number of coolant channels than the real one is used in the core. Only one coolant channel is considered per assembly, instead of 18. This simplification is done to have a lighter numerical problem. Consequently, the coolant channels properties must be adapted:

- One coolant channel conducts the same mass flow rate than the 18 real ones per assembly: the area is the sum of the area of each small coolant channel.

- It is the same for the heating perimeter because the heat transfer is proportional to this value.

- The pressure drop must be the same and is linked to the hydraulic diameter: the coolant diameter of the single coolant channel is the same than for a small real coolant channel.

The heat transfer is simulated using a [HeatTransferFromExternalAppTemperature1Phase](https://mooseframework.inl.gov/source/components/HeatTransferFromExternalAppTemperature1Phase.html) heat transfer component.

!listing microreactors/gcmr/3D_core_balance_of_plant/htgr_thm_primary.i block=Components\coolant_ht_r0_c1


Transferred data in the sub app file are defined in the [Postprocessor](https://mooseframework.inl.gov/syntax/Postprocessors/index.html)’ blocks: the wall temperature is received from the main app and the fluid temperature and heat transfer coefficient are sent to the main app.

Around the core, some pipes are added to model the inlet and outlet plena and the flow of fluid in the gap (the `upcomer`) between the vessel and the outer shield (see [thm_geometry]). This pipe is connected as those in the core to the heat conduction problem. The principle is consequently the same.

!listing microreactors/gcmr/3D_core_balance_of_plant/htgr_thm_primary.i block=Components\upcomer_ht_outer_shield

!listing microreactors/gcmr/3D_core_balance_of_plant/htgr_thm_primary.i block=Components\upcomer_ht_vessel

Nevertheless, this step requires to be cautious with the definition of the wall and fluid temperatures and of the heat transfer coefficient:

- The fluid temperatures in the core and in the upcomer are not the same: a `T_fluid_core_var` and a `T_fluid_upcomer_var` are declared and received the values of the `T` variable (fluid temperature) in the corresponding blocks using the [CopyValueAux](https://mooseframework.inl.gov/source/auxkernels/CopyValueAux.html) [AuxKernel](https://mooseframework.inl.gov/syntax/AuxKernels/index.html)

!listing microreactors/gcmr/3D_core_balance_of_plant/htgr_thm_primary.i block=AuxKernels\T_fluid_core_ak

!listing microreactors/gcmr/3D_core_balance_of_plant/htgr_thm_primary.i block=AuxKernels\T_fluid_upcomer_ak

- The wall temperatures differ in the core, on the reflector and on the vessel boundaries. In the case of the core, the elemental variable `T_wall` is kept, but others are defined for the two upcomers boundaries: `T_wall_reflector_var` and `T_wall_vessel_var`. They received their values through the transfer from the mainapp, and are associated to the [Materials](https://mooseframework.inl.gov/syntax/Materials/index.html) properties using [ADCoupledVariableValueMaterial](https://mooseframework.inl.gov/source/materials/CoupledVariableValueMaterial.html) blocks.

!listing microreactors/gcmr/3D_core_balance_of_plant/htgr_thm_primary.i block=Materials\T_wall_vessel

!listing microreactors/gcmr/3D_core_balance_of_plant/htgr_thm_primary.i block=Materials\T_wall_reflector

- The heat transfer coefficients are once again different on the core coolant, reflector and vessel boundaries. New blocks are added in the [Materials](https://mooseframework.inl.gov/syntax/Materials/index.html) to define them at the vessel and reflector boundaries. They are associated with the `Hw_vessel` and `Hw_reflector` [AuxVariables](https://mooseframework.inl.gov/source/variables/AuxVariable.html) using [ADMaterialRealAux](https://mooseframework.inl.gov/source/auxkernels/MaterialRealAux.html) `AuxKernels`.

!listing microreactors/gcmr/3D_core_balance_of_plant/htgr_thm_primary.i block=Materials\Hw_vessel

!listing microreactors/gcmr/3D_core_balance_of_plant/htgr_thm_primary.i block=AuxKernels\Hw_vessel_ak

!listing microreactors/gcmr/3D_core_balance_of_plant/htgr_thm_primary.i block=Materials\Hw_reflector

!listing microreactors/gcmr/3D_core_balance_of_plant/htgr_thm_primary.i block=AuxKernels\Hw_reflector_ak

All these elements are presented in [thm_geometry]. The choice of the appropriate junctions is also important. The [JunctionParallelChannels1Phase](https://mooseframework.inl.gov/source/components/JunctionParallelChannels1Phase.html) can connect more than two [FlowChannels1Phase](https://mooseframework.inl.gov/source/components/FlowChannel1Phase.html) with different sections but these channels must be parallel and with the same orientation. By contrast, the [JunctionOneToOne1Phase](https://mooseframework.inl.gov/source/components/JunctionOneToOne1Phase.html) can connect only two [FlowChannels1Phase](https://mooseframework.inl.gov/source/components/FlowChannel1Phase.html) with the same section but they do not need to be parallel. Finally, the [VolumeJunction1Phase] can connect more than two [FlowChannels1Phase](https://mooseframework.inl.gov/source/components/FlowChannel1Phase.html) with different orientations. It explains the choice of junctions in the model, and the fact that the `upcomer_out` pipe is added to be able to connect the `upcomer` and the `plenum_inlet` which have different sections and orientations.

!media media/gcmr/3D_core_balance_of_plant/thm_geometry.png
      style=display: block;margin-left:auto;margin-right:auto;width:25%;
      id=thm_geometry
      caption= Geometry of the thermalhydraulic system

## Core model using the Heat conduction module

### Definition of the heat conduction problem

The heat conduction problem is implemented using a [Kernels](https://mooseframework.inl.gov/syntax/Kernels/index.html) block which calls itself the `T` [Variable](https://mooseframework.inl.gov/syntax/Variables/index.html):

!listing microreactors/gcmr/3D_core_balance_of_plant/mainapp_core.i block=Kernels

where the power density objects are presented below.

### Power modulation

The power density is not uniform in the fuel channels and among them. Four rings of fuel channels are defined in each assembly and release different amount of power density, as presented in [assembly_pattern]:

!media media/gcmr/3D_core_balance_of_plant/description_assembly_pattern.png
      style=display: block;margin-left:auto;margin-right:auto;width:40%;
      id=assembly_pattern
      caption=Assembly pattern

The average power density for each fuel channel (depending on the ring) is computed and then used to get the power density.


Moreover, the power density also varies along each fuel channel. A cosine power shape is adopted. The power density starts at the half of the average power density at the bottom of the fuel channel, increases and reaches a maximum in the middle of the fuel channel and finally decreases to reach one half of its average value at its top.

This variation in the length of the assembly is defined using the [ParsedFunction](https://mooseframework.inl.gov/source/functions/MooseParsedFunction.html) `density_length_assemb_fn_aux`.

!listing microreactors/gcmr/3D_core_balance_of_plant/mainapp_core.i block=Functions/density_length_assemb_fn_aux

A radial power modulation in the core is also added. It reaches a maximum value in the middle of the core and minimum ones at the boundaries. A polynomial [ParsedFunction](https://mooseframework.inl.gov/source/functions/MooseParsedFunction.html) called `density_rad_core_fn` is used to represent it. It is normalized using the `dens_rad_core_normalization` [Postprocessor] (https://mooseframework.inl.gov/syntax/Postprocessors/index.html).

!listing microreactors/gcmr/3D_core_balance_of_plant/mainapp_core.i block=Functions/density_length_assemb_fn_auxdensity_rad_core_fn

!listing microreactors/gcmr/3D_core_balance_of_plant/mainapp_core.i block=Functions/density_length_assemb_fn_auxdensity_rad_core_fn_aux

Finally, a time dependent power evolution is added to take into account the fact that the core does not deliver directly a 15 MWth power. A [PiecewiseLinear](https://mooseframework.inl.gov/source/functions/PiecewiseLinear.html) function is provided to define it. The thermal power increases from 0 MWth to 15 MWth during 100 s and stays then at its operating value.

!listing microreactors/gcmr/3D_core_balance_of_plant/mainapp_core.i block=Functions/density_time_shape_fn

Ultimately, all these functions and average power densities are multiplied using [ParsedFunctions](https://mooseframework.inl.gov/source/functions/MooseParsedFunction.html) and these ones are associated to [AuxVariables](https://mooseframework.inl.gov/source/variables/AuxVariable.html) using [AuxKernels](https://mooseframework.inl.gov/syntax/AuxKernels/index.html). These power density [AuxVariables] (https://mooseframework.inl.gov/source/variables/AuxVariable.html) are used to solve the problem in the [Kernels] (https://mooseframework.inl.gov/syntax/Kernels/index.html) block.

!listing microreactors/gcmr/3D_core_balance_of_plant/mainapp_core.i block=AuxKernels


### Materials properties

The fuel and graphite properties (specific heat and thermal conductivity) are described using temperature dependent models provided in the [MOOSE resources](https://mooseframework.inl.gov/bison/source/materials/GraphiteMatrixThermal.html) but that are not directly available with the Heat conduction model. They are consequently defined with functions that are called by the [Materials](https://mooseframework.inl.gov/syntax/Materials/index.html) sub blocks, of [HeaConductionMaterial](https://mooseframework.inl.gov/source/materials/HeatConductionMaterial.html) type. Approximative values depending on the temperature [donner la source ?] are provided for the SA508 steel properties using the same method.

!listing microreactors/gcmr/3D_core_balance_of_plant/mainapp_core.i block=Materials

The solid materials densities are defined as constant values, using [Density](https://mooseframework.inl.gov/source/materials/Density.html) sub blocks of the [Materials](https://mooseframework.inl.gov/syntax/Materials/index.html). There is namely no way to simulate the volume change of these solid materials with the temperature change. A temperature-dependent model would induce a virtual loss of materials, which is nonsense.


### Transfers with the sub app

The link between the sub app and the main app is defined in the main app input file. A `thm` sub block in the [MultiApps](https://mooseframework.inl.gov/syntax/MultiApps/index.html) block is written. Its type is [TransientMultiApp](https://mooseframework.inl.gov/source/multiapps/TransientMultiApp.html).

!listing microreactors/gcmr/3D_core_balance_of_plant/mainapp_core.i block=MultiApps

Some postprocessors and variables are defined to measure data that should be shared with the sub app. The wall temperature is sent to the sub app while the heat transfer coefficient and fluid temperature are received from the thermal hydraulic app.

Transferred data are defined in the [Transfers](https://mooseframework.inl.gov/syntax/Transfers/index.html) block. The `from_multi_app` or `to_multi_app` parameter indicates that data are received or sent to the ‘thm’ multiapp defined above.

Firstly, it contains blocks about data that are imported in the main app to be contained in the results file, but that are not used for the main app problem resolution: the fluid minimum and maximum temperatures, the total mass flow rate, the average pressure, and the power received by the fluid among others. Their type is defined as [MultiAppPostProcessorTransfer](https://mooseframework.inl.gov/source/transfers/MultiAppPostprocessorTransfer.html). It also indicates the origin and the destination of data among the main and sub apps postprocessors.

Some sub blocks are added to define received data from the sub app that are used to solve the heat conduction problem: `T_fluid_core_from_thm`, `T_fluid_upcomer_from_thm`, `Hw_core_from_thm`, `Hw_reflector_from_thm` and `Hw_vessel_from_thm`, whose type is [MultiAppGeneralFieldNearestNodeTransfer](https://mooseframework.inl.gov/source/transfers/MultiAppGeneralFieldNearestNodeTransfer.html). It indicates that data used in the main app coming from the sub app are those computed at the nearest nodes from those of the main app.

Finally, `T_wall_to_thm`, `T_wall_vessel_to_thm` and `T_wall_reflector_to_thm` sub blocks are added, and their type is [MultiAppGeneralFieldUserObjectTransfer](https://mooseframework.inl.gov/source/transfers/MultiAppGeneralFieldUserObjectTransfer.html). They receives data computed in a [UserObjects](https://mooseframework.inl.gov/syntax/UserObjects/index.html) block and sends it to the `thm` multiapp.

!listing microreactors/gcmr/3D_core_balance_of_plant/mainapp_core.i block=Transfers

The [UserObjects](https://mooseframework.inl.gov/syntax/UserObjects/index.html) block contains a `T_wall_uo`, a `T_wall_vessel_uo` and a `T_wall_reflector_uo`sub blocks, whose type is [NearestPointLayeredSideAverage](https://mooseframework.inl.gov/source/userobjects/NearestPointLayeredSideAverage.html). It indicates the way to compute along the provided `direction` data transmitted to the sub app on the `boundary`. The discretization is namely different: it is given in the `num_layers` parameter.  In the case of the coolant channels of the core, it receives their positions through a `positions.txt` file provided in the `points_file` parameter.

!listing microreactors/gcmr/3D_core_balance_of_plant/mainapp_core.i block=UserObjects

Ultimately, the fluid temperature and heat transfer coefficient are applied to the boundaries between the fluid and the solid materials. It is done using a [BCs](https://mooseframework.inl.gov/syntax/BCs/index.html) block, which specifies the boundary conditions at the coolant, reflector and vessel boundaries (and also the radiative loss from the vessel, presented below):

!listing microreactors/gcmr/3D_core_balance_of_plant/mainapp_core.i block=BCs

### Radiative loss of power

Some radiative loss is also considered. The reflector of the core transmits heat by radiative transfer to the vessel of the reactor. This is taken into account using a [GapHeatTransfer](https://mooseframework.inl.gov/source/bcs/GapHeatTransfer.html) block which is in the [ThermalContact](https://mooseframework.inl.gov/syntax/ThermalContact/index.html) block.

!listing microreactors/gcmr/3D_core_balance_of_plant/mainapp_core.i block=ThermalContact

The vessel also loses some power by radiation to the environment. To represent it, a `vessel_radiative_loss_BC` of [RadiativeHeatFluxBC](https://mooseframework.inl.gov/source/bcs/RadiativeHeatFluxBC.html) boundary condition type is added in the [BCs](https://mooseframework.inl.gov/syntax/BCs/index.html) block. The external temperature is defined at 300 K.
