###################################################
###########   GEOMETRICAL PARAMETERS   ############
###################################################

radius_fuel = 0.00794 # m
section_fuel_channel = '${fparse pi * radius_fuel * radius_fuel}'
nb_fuel_per_assembly = 42
nb_assembly = 55
assembly_height = 2.
channel_n_elems = 10
upcomer_heated_n_elems = ${channel_n_elems}

# calculation of the core origin used to indicate the correct location of the core

L_plenum_inlet = 0.2
L_upcomer_out = 0.5
x_core_origin = '${fparse L_plenum_inlet + L_upcomer_out}'
# x_upcomer_origin = ${fparse assembly_height + x_core_origin}

fuel_channels_area = '${fparse section_fuel_channel * nb_assembly * nb_fuel_per_assembly}'

radius_core = 0.9

###################################################
###########      POWER PARAMETERS      ############
###################################################

tot_power = 15e6 # W
tot_power_assembly = '${fparse tot_power / (nb_assembly)}' # W

# The goal of that variation is to have fuel channels in the center of the assembly that emit more heat (30% more) than the fuel channels on the boundary of the assembly.
# There are 4 rings of fuel channels. The following calculation steps are :
#   1 - definition of the number of fuel channels per ring
#   2 - definition of the factor applied to the minimal power density: for instance, for ring_4, the emitted power is the minimal one, so the factor is 1 (or 1 + 0%). For ring_1, the fuel channels emit 30% more than the minimal value, so the factor is 1.3 = 1 + 30%
#   3 - the minimal power emitted by a fuel channel (so those of ring_4) is calculated, and then the corresponding power density,
#   4 - using the factors, the power density for each ring is calculated.
#   5 - in the [Kernels] block, these values are associated to the corresponding fuel channels depending on their ring

nb_fuel_ring1 = 6
nb_fuel_ring2 = 6
nb_fuel_ring3 = 12
nb_fuel_ring4 = 18

factor_density_ring1 = 1.3
factor_density_ring2 = 1.25
factor_density_ring3 = 1.2
factor_density_ring4 = 1

power_channel_min_avg = '${fparse tot_power_assembly / (nb_fuel_ring1 * factor_density_ring1 + nb_fuel_ring2 * factor_density_ring2 + nb_fuel_ring3 * factor_density_ring3 + nb_fuel_ring4 * factor_density_ring4)}'
density_channel_min_avg = '${fparse power_channel_min_avg / (section_fuel_channel * assembly_height)}'

density_ring1_avg = '${fparse factor_density_ring1 * density_channel_min_avg}'
density_ring2_avg = '${fparse factor_density_ring2 * density_channel_min_avg}'
density_ring3_avg = '${fparse factor_density_ring3 * density_channel_min_avg}'
density_ring4_avg = '${fparse factor_density_ring4 * density_channel_min_avg}'

# The goal of that variation is to have a sinusoidal variation along the length of the channels of the power density in the fuel channels.
# The power density has to be smaller at the inlet and outlet of the channels, and bigger in the middle.
# The following repartition is applied:
#
#    f(z)
#     ^                       .
#     |         .             ^              .                      f(z) = B + A * sin( pi * z / L )
#     |   .                   |  A                 .                density_ring_i = density_ringi_avg * f(z)
#     |.                      |                       .
#     |   ^                                                         A and B must be chosen so that we have:
#     |   | B                                                       < density_ringi > = density_ringi_avg
#     X-------------------------------------------------> x
#   x = 0                                            x = L = assembly_height
#
# Here we chose B and we compute the implied value of A :

B = 0.5
A = '${fparse (1 - B) * pi / 2}'

# radial variation in the core:

#                    g(r)
#                     ^
#                     |                                              g(r) = (-r^2 + R^2)/(2 * R^2) + C        => g_max = 1
#                    .|.                                                                                      => g_min = 0.5
#           .         |          .
#       .             |              .
#     .               |                .
#                     |                   ^
#                     |                   | C
#   --+---------------X----------------+----------------> r = sqrt(z^2 + y^2)
#    -R                                R
#

C = 0.5

###################################################
###########        TEMPERATURES        ############
###################################################

T_ini = 490 # K

# This tests an action used to exchange T_wall, T_fluid and HTC between
# a heat conduction simulation and a THM simulation

###################################################
###########    GRAPHITE PARAMETERS     ############
###################################################

# A_cp_g = 2253.74089
# B_cp_g = 0.03812164
# C_cp_g = -377700.14
# D_cp_g = -181791871
# E_cp_g = 6.6655e10
# F_cp_g = -6.012e12

###################################################
###########    MATERIALS PROPERTIES    ############
###################################################

density_TRISO_U = 14300
density_TRISO_buffer_layer = 1900
density_TRISO_PyC = 1900
density_TRISO_SiC = 3200

radius_TRISO_U = 4.25e-4
radius_TRISO_buffer_layer = 4.75e-4
radius_TRISO_PyC_inner = 5.1e-4
radius_TRISO_SiC = 5.45e-4
radius_TRISO_PyC_outer = 5.8e-4

mass_TRISO_U = '${fparse (4/3) * pi * density_TRISO_U * radius_TRISO_U * radius_TRISO_U * radius_TRISO_U}'
mass_TRISO_buffer_layer = '${fparse (4/3) * pi * density_TRISO_buffer_layer * (radius_TRISO_buffer_layer * radius_TRISO_buffer_layer * radius_TRISO_buffer_layer - radius_TRISO_U * radius_TRISO_U * radius_TRISO_U)}'
mass_TRISO_PyC_inner = '${fparse (4/3) * pi * density_TRISO_PyC * (radius_TRISO_PyC_inner * radius_TRISO_PyC_inner * radius_TRISO_PyC_inner - radius_TRISO_buffer_layer * radius_TRISO_buffer_layer * radius_TRISO_buffer_layer)}'
mass_TRISO_SiC = '${fparse (4/3) * pi * density_TRISO_SiC * (radius_TRISO_SiC * radius_TRISO_SiC * radius_TRISO_SiC - radius_TRISO_PyC_inner * radius_TRISO_PyC_inner * radius_TRISO_PyC_inner)}'
mass_TRISO_PyC_outer = '${fparse (4/3) * pi * density_TRISO_PyC * (radius_TRISO_PyC_outer * radius_TRISO_PyC_outer * radius_TRISO_PyC_outer - radius_TRISO_SiC * radius_TRISO_SiC * radius_TRISO_SiC)}'

mass_TRISO = '${fparse mass_TRISO_U + mass_TRISO_buffer_layer + mass_TRISO_PyC_inner + mass_TRISO_SiC + mass_TRISO_PyC_outer}'

density_TRISO = '${fparse mass_TRISO / ( (4/3) * pi * radius_TRISO_PyC_outer * radius_TRISO_PyC_outer * radius_TRISO_PyC_outer)}' # kg/m3, homogeneous density over an entire UN TRISO particle

density_graphite = 2260 # kg/m3
fuel_packing_fraction = 0.5 # proportion of TRISO pqrticule in the graphite matrix in the fuel channels

density_fuel = '${fparse density_TRISO * fuel_packing_fraction + ( 1 - fuel_packing_fraction) * density_graphite}' # kg/m3

density_steel = 8000 # kg/m3

[Mesh]
  [file_mesh]
    type = FileMeshGenerator
    file = core_mesh_with_vessel_in.e
  []
[]

[Functions]

  #####
  # 0 #   Solid materials properties
  #####

  [graphite_specific_heat_fn]
    type = PiecewiseLinear
    x = '300 350 400 450  500  550  600  650  700  750  800  850  900  950  1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 1500 1550 1600'
    y = '713 856 990 1111 1218 1310 1390 1459 1520 1573 1619 1661 1697 1730 1759 1786 1809 1831 1851 1869 1885 1900 1914 1927 1939 1950 1960'
  []

  [graphite_th_conductivity_fn]
    type = PiecewiseLinear
    x = '300 350 400 450 500 550 600 650 700  750  800  850  900  950  1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 1500 1550 1600'
    y = '135 130 125 120 115 110 106 102 97.9 94.0 90.3 86.8 83.4 80.2 77.1 74.3 71.6 69.0 66.6 64.4 62.3 60.4 58.7 57.1 55.7 54.5 53.4'
  []

  [fuel_specific_heat_fn]
    type = PiecewiseLinear
    x = '300 350 400 450 500 550 600  650  700  750  800  850  900  950  1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 1500 1550 1600'
    y = '661 731 802 870 930 983 1029 1069 1104 1134 1161 1184 1205 1223 1239 1254 1267 1279 1290 1300 1309 1317 1324 1331 1338 1344 1349'
  []

  [fuel_th_conductivity_fn]
    type = PiecewiseLinear
    x = '300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 1500 1550 1600'
    y = '55  53  51  49  47  45  44  42  41  39  38  36  35  34  33   32   30   29   29   28   27   26   26   25   24   24   23'
  []

  #####
  # 1 #   Constant values used to define the power shape functions
  #####

  # core parameters
  [fuel_channels_area_fn]
    type = ConstantFunction
    value = ${fuel_channels_area}
  []
  [assembly_height_fn]
    type = ConstantFunction
    value = ${assembly_height}
  []
  [radius_core_fn]
    type = ConstantFunction
    value = ${radius_core}
  []

  # average density in the fuel channels depending on their ring position in the assembly
  [density_ring1_avg_fn]
    type = ConstantFunction
    value = ${density_ring1_avg}
  []
  [density_ring2_avg_fn]
    type = ConstantFunction
    value = ${density_ring2_avg}
  []
  [density_ring3_avg_fn]
    type = ConstantFunction
    value = ${density_ring3_avg}
  []
  [density_ring4_avg_fn]
    type = ConstantFunction
    value = ${density_ring4_avg}
  []

  # coefficients of the power shape function for the variation in the length of the assembly
  [A_fn]
    type = ConstantFunction
    value = ${A}
  []
  [B_fn]
    type = ConstantFunction
    value = ${B}
  []

  # coefficient of the power shape function for the radial variation in core
  [C_fn]
    type = ConstantFunction
    value = ${C}
  []

  # position of the core
  [x_core_origin_fn]
    type = ConstantFunction
    value = ${x_core_origin}
  []

  #####
  # 2 #   Power shape functions (variation in the length) depending on the fuel channel ring in the assembly
  #####

  [density_length_assemb_fn_aux]
    type = ParsedFunction
    expression = 'B_fn + A_fn * sin(pi * (x - x_core_origin_fn) / assembly_height_fn)'
    symbol_names = 'A_fn B_fn assembly_height_fn x_core_origin_fn'
    symbol_values = 'A_fn B_fn assembly_height_fn x_core_origin_fn'
  []

  [density_ring1_fn_aux]
    type = ParsedFunction
    symbol_names = 'density_length_assemb_fn_aux density_ring1_avg_fn'
    symbol_values = 'density_length_assemb_fn_aux density_ring1_avg_fn'
    expression = 'density_ring1_avg_fn * density_length_assemb_fn_aux'
  []
  [density_ring2_fn_aux]
    type = ParsedFunction
    symbol_names = 'density_length_assemb_fn_aux density_ring2_avg_fn'
    symbol_values = 'density_length_assemb_fn_aux density_ring2_avg_fn'
    expression = 'density_ring2_avg_fn * density_length_assemb_fn_aux'
  []
  [density_ring3_fn_aux]
    type = ParsedFunction
    symbol_names = 'density_length_assemb_fn_aux density_ring3_avg_fn'
    symbol_values = 'density_length_assemb_fn_aux density_ring3_avg_fn'
    expression = 'density_ring3_avg_fn * density_length_assemb_fn_aux'
  []
  [density_ring4_fn_aux]
    type = ParsedFunction
    symbol_names = 'density_length_assemb_fn_aux density_ring4_avg_fn'
    symbol_values = 'density_length_assemb_fn_aux density_ring4_avg_fn'
    expression = 'density_ring4_avg_fn * density_length_assemb_fn_aux'
  []

  #####
  # 3 #    Radial power shape function in the core: higher in the middle of the core, lower at the boundary
  #####

  # non normalized power shape function
  [density_rad_core_fn_aux]
    type = ParsedFunction
    expression = '(- y * y - z * z + radius_core_fn * radius_core_fn) / (2 * radius_core_fn * radius_core_fn) + C_fn'
    symbol_names = 'radius_core_fn C_fn'
    symbol_values = 'radius_core_fn C_fn'
  []

  # normalized power shape function
  [density_rad_core_fn]
    type = ParsedFunction
    expression = '(fuel_channels_area_fn * assembly_height_fn / dens_rad_core_normalization) * density_rad_core_fn_aux'
    symbol_names = 'fuel_channels_area_fn dens_rad_core_normalization density_rad_core_fn_aux assembly_height_fn'
    symbol_values = 'fuel_channels_area_fn dens_rad_core_normalization density_rad_core_fn_aux assembly_height_fn'
  []

  #####
  # 4 #   power variation with time
  #####

  [density_time_shape_fn]
    type = PiecewiseLinear
    x = '0 100'
    y = '0 1'
  []

  #####
  # 5 #   Final power shape functions
  #####

  [density_ring1_fn]
    type = ParsedFunction
    expression = 'density_ring1_fn_aux * density_rad_core_fn * density_time_shape_fn'
    symbol_names = 'density_ring1_fn_aux density_rad_core_fn density_time_shape_fn'
    symbol_values = 'density_ring1_fn_aux density_rad_core_fn density_time_shape_fn'
  []
  [density_ring2_fn]
    type = ParsedFunction
    expression = 'density_ring2_fn_aux * density_rad_core_fn * density_time_shape_fn'
    symbol_names = 'density_ring2_fn_aux density_rad_core_fn density_time_shape_fn'
    symbol_values = 'density_ring2_fn_aux density_rad_core_fn density_time_shape_fn'
  []
  [density_ring3_fn]
    type = ParsedFunction
    expression = 'density_ring3_fn_aux * density_rad_core_fn * density_time_shape_fn'
    symbol_names = 'density_ring3_fn_aux density_rad_core_fn density_time_shape_fn'
    symbol_values = 'density_ring3_fn_aux density_rad_core_fn density_time_shape_fn'
  []
  [density_ring4_fn]
    type = ParsedFunction
    expression = 'density_ring4_fn_aux * density_rad_core_fn * density_time_shape_fn'
    symbol_names = 'density_ring4_fn_aux density_rad_core_fn density_time_shape_fn'
    symbol_values = 'density_ring4_fn_aux density_rad_core_fn density_time_shape_fn'
  []

[]

[Materials]
  [graphite_thermal]
    type = HeatConductionMaterial
    #       matrix channels   reflector   moderator_pincell   outer_shield
    block = ' 1 2 3 4 5 6        10           101 103             250'
    temp = T
    thermal_conductivity_temperature_function = graphite_th_conductivity_fn # W/(m.K)
    specific_heat_temperature_function = graphite_specific_heat_fn # J/(kg.K)
  []

  [fuel_thermal]
    type = HeatConductionMaterial
    #         fuel_pincell
    block = '301 303 401 403 501 503 601 603'
    temp = T
    thermal_conductivity_temperature_function = fuel_th_conductivity_fn # W/(m.K)
    specific_heat_temperature_function = fuel_specific_heat_fn # J/(kg.K)
  []

  [steel_thermal]
    type = HeatConductionMaterial
    block = '450'
    temp = T
    thermal_conductivity = 45
    specific_heat = 466
  []

  [graphite_density]
    type = Density
    #       matrix channels   reflector   moderator_pincell   outer_shield
    block = ' 1 2 3 4 5 6        10           101 103             250'
    density = ${density_graphite} # kg/m3
  []

  [fuel_density]
    type = Density
    #         fuel_pincell
    block = '301 303 401 403 501 503 601 603'
    density = ${density_fuel} # kg/m3
  []

  [steel_density]
    type = Density
    block = '450'
    density = ${density_steel} # kg/m3
  []
[]

[Variables]
  [T]
  []
[]

[ICs]
  [T_ic]
    type = ConstantIC
    variable = T
    value = ${T_ini}
  []
[]

[AuxVariables]

  [T_fluid_upcomer]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = ${T_ini}
  []
  [T_fluid_core]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = ${T_ini}
  []
  [htc]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = 1000
  []
  [htc_reflector]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = 1000
  []
  [htc_vessel]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = 1000
  []
  [density_ring1]
  []
  [density_ring2]
  []
  [density_ring3]
  []
  [density_ring4]
  []
[]

[AuxKernels]
  [density_ring1_aux]
    type = FunctionAux
    function = density_ring1_fn
    variable = density_ring1
    execute_on = 'TIMESTEP_END'
  []
  [density_ring2_aux]
    type = FunctionAux
    function = density_ring2_fn
    variable = density_ring2
    execute_on = 'TIMESTEP_END'
  []
  [density_ring3_aux]
    type = FunctionAux
    function = density_ring3_fn
    variable = density_ring3
    execute_on = 'TIMESTEP_END'
  []
  [density_ring4_aux]
    type = FunctionAux
    function = density_ring4_fn
    variable = density_ring4
    execute_on = 'TIMESTEP_END'
  []
[]

[Kernels]
  [heat_conduction]
    type = HeatConduction
    variable = T
  []
  [heat_td]
    type = TimeDerivative
    variable = T
  []
  [heat_source_fuel_ring1]
    type = CoupledForce
    variable = T
    block = '301 303'
    v = density_ring1
  []
  [heat_source_fuel_ring2]
    type = CoupledForce
    variable = T
    block = '401 403'
    v = density_ring2
  []
  [heat_source_fuel_ring3]
    type = CoupledForce
    variable = T
    block = '501 503'
    v = density_ring3
  []
  [heat_source_fuel_ring4]
    type = CoupledForce
    variable = T
    block = '601 603'
    v = density_ring4
  []
[]

[BCs]
  [vessel_radiative_loss_BC]
    type = RadiativeHeatFluxBC
    variable = T
    boundary = vessel_boundary_out
    Tinfinity = 300
    boundary_emissivity = 0.3
    view_factor = 1
  []

  [cooling_channels]
    type = CoupledConvectiveHeatFluxBC
    boundary = coolant_boundary
    T_infinity = T_fluid_core
    htc = htc
    variable = T
  []
  [upcomer_reflector_BC]
    type = CoupledConvectiveHeatFluxBC
    boundary = 'reflector_boundary'
    T_infinity = T_fluid_upcomer
    htc = htc_reflector
    variable = T
  []
  [upcomer_vessel_BC]
    type = CoupledConvectiveHeatFluxBC
    boundary = 'vessel_boundary_in'
    T_infinity = T_fluid_upcomer
    htc = htc_vessel
    variable = T
  []
[]

[ThermalContact]
  [gap_ht]
    type = GapHeatTransfer
    variable = T
    primary = reflector_boundary
    secondary = vessel_boundary_in
    emissivity_primary = 1
    emissivity_secondary = 0.3
    gap_conductivity = 0.025
    gap_geometry_type = CYLINDER
    quadrature = true
    cylinder_axis_point_1 = '${x_core_origin} 0 0'
    cylinder_axis_point_2 = '${fparse x_core_origin + assembly_height} 0 0'
  []
[]

[UserObjects]
  [T_wall_uo]
    type = NearestPointLayeredSideAverage
    points_file = 'positions.txt'
    direction = x
    num_layers = ${channel_n_elems}
    boundary = coolant_boundary
    variable = T
  []

  [T_wall_reflector_uo]
    type = NearestPointLayeredSideAverage
    points = '${x_core_origin} 0. 0.'
    direction = x
    num_layers = ${upcomer_heated_n_elems}
    boundary = reflector_boundary
    variable = T
  []

  [T_wall_vessel_uo]
    type = NearestPointLayeredSideAverage
    points = '${x_core_origin} 0. 0.'
    direction = x
    num_layers = ${upcomer_heated_n_elems}
    boundary = vessel_boundary_in
    variable = T
  []

[]

[Executioner]
  type = Transient
  end_time = 4000

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    growth_factor = 1.2
    cutback_factor = 0.8
  []
  dtmin = 0.1
  dtmax = 100

  steady_state_detection = true
  steady_state_start_time = 1900

  nl_abs_tol = 1e-10
  abort_on_solve_fail = true

  solve_type = 'PJFNK'

  automatic_scaling = true
  compute_scaling_once = false
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart '
  petsc_options_value = 'hypre boomeramg 100'

  #fixed_point_abs_tol = 1e-6
  #fixed_point_max_its = 5
  #fixed_point_rel_tol = 1e-3
  #accept_on_max_fixed_point_iteration = true
[]

[MultiApps]
  [thm]
    type = TransientMultiApp
    input_files = htgr_thm_primary.i
    execute_on = 'TIMESTEP_END'
    bounding_box_padding = '0.1 0.1 0.1'
    sub_cycling = true
    # tolerate_failure = true
    max_failures = 400
    max_procs_per_app = 5
    output_in_position = true
  []
[]

[Transfers]
  [T_wall_to_thm]
    type = MultiAppUserObjectTransfer
    to_multi_app = thm
    user_object = T_wall_uo
    variable = T_wall
  []
  [T_wall_reflector_to_thm]
    type = MultiAppUserObjectTransfer
    to_multi_app = thm
    user_object = T_wall_reflector_uo
    variable = T_wall_reflector
  []
  [T_wall_vessel_to_thm]
    type = MultiAppUserObjectTransfer
    to_multi_app = thm
    user_object = T_wall_vessel_uo
    variable = T_wall_vessel
  []

  [T_fluid_core_from_thm]
    type = MultiAppNearestNodeTransfer
    from_multi_app = thm
    source_variable = T_fluid_core_var
    variable = T_fluid_core
  []
  [T_fluid_upcomer_from_thm]
    type = MultiAppNearestNodeTransfer
    from_multi_app = thm
    source_variable = T_fluid_upcomer_var
    variable = T_fluid_upcomer
  []

  [Hw_core_from_thm]
    type = MultiAppNearestNodeTransfer
    from_multi_app = thm
    source_variable = Hw_core_and_hx
    variable = htc
  []
  [Hw_reflector_from_thm]
    type = MultiAppNearestNodeTransfer
    from_multi_app = thm
    source_variable = Hw_reflector
    variable = htc_reflector
  []
  [Hw_vessel_from_thm]
    type = MultiAppNearestNodeTransfer
    from_multi_app = thm
    source_variable = Hw_vessel
    variable = htc_vessel
  []

  # ### avg or sum coolant channels Postprocessors - transfers

  # [T_fluid_max_from_thm]
  #   type = MultiAppPostprocessorTransfer
  #   from_multi_app = thm
  #   from_postprocessor = core_T_out
  #   to_postprocessor = T_fluid_max
  #   reduction_type = maximum
  # []
  # [T_fluid_min_from_thm]
  #   type = MultiAppPostprocessorTransfer
  #   from_multi_app = thm
  #   from_postprocessor = core_T_in
  #   to_postprocessor = T_fluid_min
  #   reduction_type = minimum
  # []
  # [m_dot_in_from_thm]
  #   type = MultiAppPostprocessorTransfer
  #   from_multi_app = thm
  #   from_postprocessor = hx_pri_m_dot_in
  #   to_postprocessor = m_dot_in_fluid
  #   reduction_type = sum
  # []
  # [p_in_avg_from_thm]
  #   type = MultiAppPostprocessorTransfer
  #   from_multi_app = thm
  #   from_postprocessor = core_p_in
  #   to_postprocessor = p_in_fluid_avg
  #   reduction_type = average
  # []
  # [power_fluid_from_thm]
  #   type = MultiAppPostprocessorTransfer
  #   from_multi_app = thm
  #   from_postprocessor = core_power_flux
  #   to_postprocessor = power_fluid
  #   reduction_type = sum
  # []

  # #####################################
  # ########## heat exchanger ###########
  # #####################################

  # [hx_sec_m_dot_in_from_thm]
  #   type = MultiAppPostprocessorTransfer
  #   from_multi_app = thm
  #   from_postprocessor = hx_sec_m_dot_in
  #   to_postprocessor = hx_sec_m_dot_in
  #   reduction_type = maximum
  # []
  # [hx_sec_m_dot_out_from_thm]
  #   type = MultiAppPostprocessorTransfer
  #   from_multi_app = thm
  #   from_postprocessor = hx_sec_m_dot_out
  #   to_postprocessor = hx_sec_m_dot_out
  #   reduction_type = maximum
  # []

  # [hx_sec_T_in_from_thm]
  #   type = MultiAppPostprocessorTransfer
  #   from_multi_app = thm
  #   from_postprocessor = hx_sec_T_in
  #   to_postprocessor = hx_sec_T_in
  #   reduction_type = maximum
  # []
  # [hx_sec_T_out_from_thm]
  #   type = MultiAppPostprocessorTransfer
  #   from_multi_app = thm
  #   from_postprocessor = hx_sec_T_out
  #   to_postprocessor = hx_sec_T_out
  #   reduction_type = maximum
  # []

  # [hx_sec_p_in_from_thm]
  #   type = MultiAppPostprocessorTransfer
  #   from_multi_app = thm
  #   from_postprocessor = hx_sec_p_in
  #   to_postprocessor = hx_sec_p_in
  #   reduction_type = maximum
  # []
  # [hx_sec_p_out_from_thm]
  #   type = MultiAppPostprocessorTransfer
  #   from_multi_app = thm
  #   from_postprocessor = hx_sec_p_out
  #   to_postprocessor = hx_sec_p_out
  #   reduction_type = maximum
  # []

  # [hx_sec_power_in_from_thm]
  #   type = MultiAppPostprocessorTransfer
  #   from_multi_app = thm
  #   from_postprocessor = hx_sec_power_in
  #   to_postprocessor = hx_sec_power_in
  #   reduction_type = maximum
  # []
  # [hx_sec_power_out_from_thm]
  #   type = MultiAppPostprocessorTransfer
  #   from_multi_app = thm
  #   from_postprocessor = hx_sec_power_out
  #   to_postprocessor = hx_sec_power_out
  #   reduction_type = maximum
  # []

  # [hx_sec_power_flux_from_thm]
  #   type = MultiAppPostprocessorTransfer
  #   from_multi_app = thm
  #   from_postprocessor = hx_sec_power_flux
  #   to_postprocessor = hx_sec_power_flux
  #   reduction_type = maximum
  # []

[]

[Outputs]
  exodus = true
  csv = true
[]

[Postprocessors]

  ##### Postprocessors used by the multiapp

  [T_wall_avg]
    type = SideAverageValue
    variable = T
    boundary = coolant_boundary
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [T_fluid_core_avg]
    type = ElementAverageValue
    variable = T_fluid_core
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [htc_avg]
    type = ElementAverageValue
    variable = htc
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # ##### Fluid postprocessors received using the multiapp

  # [T_fluid_core_max]
  #   type = Receiver
  # []
  # [T_fluid_core_min]
  #   type = Receiver
  # []
  # [p_in_fluid_avg]
  #   type = Receiver
  # []
  # [power_fluid]
  #   type = Receiver
  # []
  # [m_dot_in_fluid]
  #   type = Receiver
  # []

  ##### assembly postprocessors

  [fuel_temp_avg]
    type = ElementAverageValue
    variable = T
    block = '301 303 401 403 501 503 601 603'
  []
  [fuel_temp_max]
    type = ElementExtremeValue
    variable = T
    block = '301 303 401 403 501 503 601 603'
  []
  [fuel_temp_min]
    type = ElementExtremeValue
    variable = T
    block = '301 303 401 403 501 503 601 603'
    value_type = min
  []
  [mod_temp_avg]
    type = ElementAverageValue
    variable = T
    #       matrix channels   reflector   moderator_pincell
    block = '1 2 3 4 5 6           10           101 103'
  []
  [mod_temp_max]
    type = ElementExtremeValue
    variable = T
    #       matrix channels   reflector   moderator_pincell
    block = '1 2 3 4 5 6           10           101 103'
    value_type = max
  []
  [mod_temp_min]
    type = ElementExtremeValue
    variable = T
    #       matrix channels   reflector   moderator_pincell
    block = '1 2 3 4 5 6           10           101 103'
    value_type = min
  []
  [power_ring1]
    type = ElementIntegralVariablePostprocessor
    block = '301 303'
    variable = density_ring1
  []
  [power_ring2]
    type = ElementIntegralVariablePostprocessor
    block = '401 403'
    variable = density_ring2
  []
  [power_ring3]
    type = ElementIntegralVariablePostprocessor
    block = '501 503'
    variable = density_ring3
  []
  [power_ring4]
    type = ElementIntegralVariablePostprocessor
    block = '601 603'
    variable = density_ring4
  []
  [power_through_coolant_boundaries]
    type = ConvectiveHeatTransferSideIntegral
    T_fluid_var = T_fluid_core
    htc_var = htc
    T_solid = T
    boundary = coolant_boundary
  []
  #[power_through_reflector_boundary]
  #  type = ConvectiveHeatTransferSideIntegral
  #  T_fluid_var = T_fluid_upcomer
  #  htc_var = htc
  #  T_solid = T
  #  boundary = reflector_boundary
  #[]
  #[power_through_vessel_boundary_in]
  #  type = ConvectiveHeatTransferSideIntegral
  #  T_fluid_var = T_fluid_upcomer
  #  htc_var = htc
  #  T_solid = T
  #  boundary = vessel_boundary_in
  #[]
  # [heat_balance]
  #   type = ParsedPostprocessor
  #   pp_names = 'power_ring1 power_ring2 power_ring3 power_ring4 fuel_power'
  #   function = '(power_ring1+ power_ring2+ power_ring3 + power_ring4 - fuel_power )/fuel_power'
  # []

  [dens_rad_core_normalization]
    type = FunctionElementIntegral
    function = density_rad_core_fn_aux
    block = '301 303 401 403 501 503 601 603'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
