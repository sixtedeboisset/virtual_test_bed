# Model description

## Mesh

A complete mesh representing each assembly constituting the core pattern is built using the [Mesh](https://mooseframework.inl.gov/syntax/Mesh/index.html) and `MeshGenerator` MOOSE tools. The input file is `core_mesh.i` and the mesh is in the file `core_mesh_in.e`.

## Transfer - general principle

The general idea of the simulation is to transfer data between one main app, which solves the heat conduction 3D problem of the core, and a sub app, which solves the 1D thermal hydraulic simulation of the loops. These data are about the wall and the fluid temperatures and the heat transfer coefficient.

This principle is presented in [transfer_principle]:

!media media/gcmr/3D_core_balance_of_plant/transfer_principle.png
      style=display: block;margin-left:auto;margin-right:auto;width:25%;
      id=transfer_principle
      caption= Principle of the transfer

## THM sub app

The input file built for the thermal hydraulic simulation of the two loops is almost the same than the one presented in the 1D simulation of the balance of plant. The only changes are about the core.

A smaller number of coolant channels than the real one is used in the core. Only one coolant channel is considered per assembly, instead of 18. This simplification is done to have a lighter numerical problem. Consequently, the coolant channels properties must be adapted:

- One coolant channel conducts the same mass flow rate than the 18 real ones per assembly: the section `A` is the sum of the sections of each small coolant channel.

- It is the same for the heating perimeter `P_hf` because the heat transfer is proportional to this value.

- The pressure drop must be the same and is linked with the hydraulic diameter `D_h`: the diameter of the single coolant channel is the same than for a small real coolant channel.

The validity of this simplification has been tested. Models of one assembly with 18 coolant channels or with only one have been built and the results of the simulations have been compared. They are equivalent, so the simplified model of the assembly coupled with only one coolant channel as been used to build the core.

[BLOCK–channelproperties]

The heat transfer is simulated using a [HeatTransferFromExternalAppTemperature1Phase](https://mooseframework.inl.gov/source/components/HeatTransferFromExternalAppTemperature1Phase.html) heat transfer component.

[BLOCK–heattransfer]

The transferred data in the sub app file are defined in the [Postprocessors](https://mooseframework.inl.gov/syntax/Postprocessors/index.html)’ blocks: the wall temperature is received from the main app and the fluid temperature and heat transfer coefficient are sent to the main app.

[BLOCK–postprocessors]

## Core model

### Definition of the heat conduction problem

The heat conduction problem obeys to the equation presented in the previous section. It is implemented using a [Kernels](https://mooseframework.inl.gov/syntax/Kernels/index.html) block which calls the `T` [Variable](https://mooseframework.inl.gov/syntax/Variables/index.html):

[BLOCK-Variable]

[BLOCK-Kernels]

### Power modulation

The power density is not uniform in the fuel channels and among them. Four rings of fuel channels are defined in each assembly and release different amount of power density, as presented in [assembly_pattern]:

!media media/gcmr/3D_core_balance_of_plant/description_assembly_pattern.png
      style=display: block;margin-left:auto;margin-right:auto;width:40%;
      id=assembly_pattern
      caption=Assembly pattern

This radial variation is defined using [AuxVariables](https://mooseframework.inl.gov/source/variables/AuxVariable.html) that are called by [AuxKernels](https://mooseframework.inl.gov/syntax/AuxKernels/index.html):

[BLOCKS-aux]

Moreover, the power density also varies along each fuel channel. A cosine power shape is adopted. The power density starts at the half of the average power density at the bottom of the fuel channel, increases and reaches a maximum in the middle of the fuel channel and finally decreases to reach one half of its average value at the top. The ratio between the power density along the assembly and the average value for each fuel channel (ratio which is the same for all the fuel channels) is presented in [power_var_length]:

!media media/gcmr/3D_core_balance_of_plant/power_var_length.png
      style=display: block;margin-left:auto;margin-right:auto;width:50%;
      id=power_var_length
      caption=Ratio between the power density and the average density of each fuel channel


This variation in the length of the assembly is defined using the functions `density_ring1_fn`, `density_ring2_fn`, `density_ring3_fn` and `density_ring4_fn` called in the [AuxKernels](https://mooseframework.inl.gov/syntax/AuxKernels/index.html) block.

[BLOCK-functions]

### Materials properties

The fuel and graphite properties (specific heat and heat transfer coefficient) are described using temperature dependent models. They are defined in the [Functions](https://mooseframework.inl.gov/syntax/Functions/index.html) and called by the [Materials](https://mooseframework.inl.gov/syntax/Materials/index.html) sub blocks.

[BLOCK-functions]
[BLOCK-materials]

The graphite and fuel densities are defined as constant values, because there is no way to model the volume change of these solid materials with the temperature change. A temperature-dependent model would induce a virtual loss of materials, which is nonsense.

### Transfers with the sub app

The link between the sub app and the main app is defined in the main app input file. A `thm` sub block in the [MultiApps](https://mooseframework.inl.gov/syntax/MultiApps/index.html) block is written and contains the information about the position of the sub app. Its type is [TransientMultiApp](https://mooseframework.inl.gov/source/multiapps/TransientMultiApp.html).

[BLOCK-multiapps]

The transferred data are defined in the [Transfers](https://mooseframework.inl.gov/syntax/Transfers/index.html) block. The `from_multi_app` or `to_multi_app` parameters indicate that the data are received or sent to the `thm` multiapp defined above.

Firstly, it contains blocks about data that are imported in the main app to be contained in the results file, but that are not used for the main app problem resolution: the fluid minimum and maximum temperatures, the total mass flow rate, the average pressure, and the power received by the fluid.

*****************
LISTE A COMPLETER
*****************

 Their type is defined as [MultiAppPostProcessorTransfer](https://mooseframework.inl.gov/source/transfers/MultiAppPostprocessorTransfer.html). It also indicates the origin and the destination of the data among the main and sub apps postprocessors.

Two sub blocks are added to define the received data from the sub app that are used to solve the heat conduction problem: `T_fluid_from_thm` and `Hw_fluid_from_thm`, whose type is [MultiAppGeneralFieldNearestNodeTransfer](https://mooseframework.inl.gov/source/transfers/MultiAppGeneralFieldNearestNodeTransfer.html). It indicates that the data used in the main app coming from the sub app are those computed at the nearest nodes from those of the main app.

Finally, a `T_wall_to_thm` sub block is added, and its type is [MultiAppGeneralFieldUserObjectTransfer](https://mooseframework.inl.gov/source/transfers/MultiAppGeneralFieldUserObjectTransfer.html). It receives data computed in a [UserObjects](https://mooseframework.inl.gov/syntax/UserObjects/index.html) block and sent it to the `thm` multiapp.

[BLOCK-Transfers]

The [UserObjects](https://mooseframework.inl.gov/syntax/UserObjects/index.html) contains a `T_wall_uo` sub block, whose type is [NearestPointLayeredSideAverage](https://mooseframework.inl.gov/source/userobjects/NearestPointLayeredSideAverage.html). It indicates the way to compute along the provided `direction` the data transmitted to the sub app on the `boundary`. This calculation is done because the discretization is different: it is given in the `num_layers` parameter.  It receives the coolant channels positions through a `positions.txt` file provided in the `points_file` parameter.

[Block-UserObjects]

Ultimately, the fluid temperature and heat transfer coefficient are applied to the boundaries between each coolant channel and the assembly in the core mesh. It is done using a [BCs](https://mooseframework.inl.gov/syntax/BCs/index.html) block, whose type is [CoupledConvectiveHeatFluxBC](https://mooseframework.inl.gov/source/bcs/CoupledConvectiveHeatFluxBC.html).

[BLOCK-BCs]


