#
#   1 DESCRIPTION
#   2 PARAMETERS
#   3 GLOBAL PARAMS
#   4 PROPERTIES - FLUID PROPERTIES
#                - CLOSURES
#                - MATERIALS
#                - HEAT STRUCTURE MATERIALS
#   5 FUNCTIONS
#   6 COMPONENTS
#   7 POSTPROCESSORS
#   8 CONTROLS - CONTROLS
#              - CONTROL LOGIC
#              - AUX VARIABLES
#              - AUX KERNELS VARIABLES
#   9 EXECUTION - PRECONDITIONING
#               - EXECUTIONER
#               - OUTPUTS
#

##########################
##########################
####                  ####
####        1         ####
####                  ####
####   DESCRIPTION    ####
####                  ####
####                  ####
##########################
##########################

################################################################################################################
################  ********************      GOAL OF THIS SIMULATION     ********************  ##################
################################################################################################################

# The following file simulates how a small HTGR coupled with a Brayton cycle reaches a steady state from an intial power equal to zero to a production of P_core = 15 MWth and P_generator = 2 MWe
# The steady state reached during this simulation (stored using a checkpoint) can be used as a starting point for a load follow transient simulation, simulating a power modulation.
#
# Such a power modulation can occur to adapt the power needed on the grid to the production
#

################################################################################################################
################  ********************      PRIMARY LOOP DESCRIPTION    ********************  ##################
################################################################################################################

# This input file models the primary loop of a High Temperature - Graphite - Gas cooled (HTGR) reactor.
#
# This input models the different components of the primary loop:
#
# 1/ Core: Heat is provided by a cylindrical heat structure, which models a coolant channel. This channel is replicated using the num_rods tool, giving the complete core.
#    Through the core, the fluid receives 15 MWth.
#    The temperature of the fluid increases from approximately 900 K to 1200 K.
#
# 2/ Pressurizer: this component is added to maintain during the simulation a pressure which stays at 90 bar approximately.
#
# 3/ Heat exchanger: the primary fluid transfers through a heat structure its power to a secondary fluid.
#    The temperature of the primary fluid decreases again to approximately 900 K
#    in the heat exchanger, the temperature of the secondary fluid increases from approximately 500 K to 1175 K
#
# 4/ Pump: a pump is used because of the pressure lost in the core and in the heat exchanger
#

################################################################################################################
################  ********************     SECONDARY LOOP DESCRIPTION   ********************  ##################
################################################################################################################

# This input file models an open, recuperated Brayton cycle with a PID
# controlled start up using a coupled motor.
#
# Heat is supplied to the system by a volumetric heat source, and a second heat
# source is used to model a recuperator. The recuperator transfers heat from the
# turbine exhaust gas to the compressor outlet gas.
#
# Initially the fluid and heat structures are at rest at ambient conditions,
# and the shaft speed is zero.
# The transient is controlled as follows:
#   * 0 - few seconds:                     Motor increases shaft speed to approx. 90,000 RPM by PID control
#   * few seconds - end of the simulation: Torque supplied by turbine increases to steady state level
#                                          as working fluid temperature increases. Torque supplied by
#                                          the motor is ramped down to 0 N-m transitioning shaft control
#                                          to the turbine at its rated speed of 95,000 RPM.

##########################
##########################
####                  ####
####        2         ####
####                  ####
####    PARAMETERS    ####
####                  ####
####                  ####
##########################
##########################

##################################################################
################      INITIAL PARAMETERS      ####################
##################################################################

T_ini = 490 #K

vel_ini_pri = 10 # m/s
vel_ini_sec = 0 # m/s

# vel_ini_default = 0. # m/s
# vel_x_ini_default = 0. # m/s
vel_y_ini_default = 0. # m/s
vel_z_ini_default = 0. # m/s

p_sec = 5.5e5

################################################################################################################
################  ********************      PRIMARY LOOP PARAMETERS     ********************  ##################
################################################################################################################

##################################################################
################      INITIAL PARAMETERS      ####################
##################################################################

pri_press = 9e6 # Pa

##################################################################
################     JUNCTIONS PARAMETERS     ####################
##################################################################

volume_jct = 1e-3

##################################################################
###################     CORE PARAMETERS    #######################
##################################################################

# parameters of the coolant channels

core_radius_coolant = 0.00635 # m
core_length_channel = 2 # m
core_channel_n_elems = 50

# numbers of channels and assemblies

core_nb_assembly = 55
core_nb_coolant_per_assembly = 18
core_nb_fuel_per_assembly = 42
core_nb_coolant_tot = '${fparse core_nb_assembly * core_nb_coolant_per_assembly}'

# other parameters of the assembly

core_lattice_pitch = 0.022 # m        optimal size of the lattice pitch to have a good burnup and a long life of the assemnly for UN TRISO
core_radius_fuel = 0.00794 # m

# calculus of the equivalent parameters of a cylindrical heat structure around the coolant channel

core_section_assembly = '${fparse 3 * ( sqrt(3) / 2) * 5 * core_lattice_pitch * 5 * core_lattice_pitch}' # calculus of the assembly section (hexagonal) using the lattice pitch (1/5 of a side of the hexagon)
core_section_fuel_channel = '${fparse pi * core_radius_fuel * core_radius_fuel}'
core_section_mod_and_coolant = '${fparse core_section_assembly - ( core_nb_fuel_per_assembly * core_section_fuel_channel )}'

core_radius_equiv_mod = '${fparse sqrt( core_section_mod_and_coolant / ( pi * core_nb_coolant_per_assembly))}' # calculus of the radius of the moderator in an equivalent cylindrical heat structure around each coolant channel
core_radius_equiv_complete_hs = '${fparse sqrt( core_section_assembly / (pi * core_nb_coolant_per_assembly))}' # calculus of the radius of the complete equivalent cylindrical heat structure around each coolant channel

# power

tot_power = 15000000 # Wth

##################################################################
###############  CORE GEOMETRICAL PARAMETERS  ####################
##################################################################

rank14_y = 0.6668
rank13_y = 0.5716
rank12_y = 0.4763
rank11_y = 0.38105
rank10_y = 0.2859
rank9_y = 0.1905
rank8_y = 0.09526
rank7_y = 0.
rank6_y = -0.09526
rank5_y = -0.1905
rank4_y = -0.2859
rank3_y = -0.38105
rank2_y = -0.4763
rank1_y = -0.5716
rank0_y = -0.6668

odd_rank_column0_z = -0.66
even_rank_column0_z = -0.495
odd_rank_column1_z = -0.33
even_rank_column1_z = -0.165
odd_rank_column2_z = 0.
even_rank_column3_z = 0.165
odd_rank_column3_z = 0.33
even_rank_column4_z = 0.495
odd_rank_column4_z = 0.66

#
#                               column 1       column 3
#                                __/\__         __/\__
#                    column 0   |      |       |      |    column 4
#                     __/\__          *         *          __/\__           rank 14
#                    |      |    *         *         *    |      |          rank 13
#                           *         *         *         *                 rank 12
#                                *         *         *                      rank 11
#                           *         *         *         *                 rank 10
#                     *          *         *         *         *            rank 9
#                           *         *         *         *                 rank 8
#                     *          *         *         *         *            rank 7
#                           *         *         *         *                 rank 6
#                     *          *         *         *         *            rank 5
#                           *         *         *         *                 rank 4
#                                *         *         *                      rank 3
#                           *         *         *         *                 rank 2
#                                *         *         *                      rank 1
#                                     *    |    *                           rank 0
#                                          |
#  y  ^                                    V
#     |                                 column 2
#     |
#     |
#     X--------->
#               z

##################################################################
##################      GAP PARAMETERS      ######################
##################################################################

core_radius = 0.9
upcomer_distance = 0.05
vessel_internal_radius = '${fparse core_radius + upcomer_distance}'

upcomer_out_n_elems = 3
upcomer_n_elems = '${core_channel_n_elems}'

upcomer_area = '${fparse pi * vessel_internal_radius * vessel_internal_radius - pi * core_radius * core_radius}'

upcomer_P_hf_in = '${fparse 2 * pi * core_radius}'
upcomer_P_hf_out = '${fparse 2 * pi * vessel_internal_radius}'
upcomer_perimeter = '${fparse upcomer_P_hf_in + upcomer_P_hf_out}'
upcomer_D_h = '${fparse 4 * upcomer_area / upcomer_perimeter}'

##################################################################
##################     PLENUM PARAMETERS    ######################
##################################################################

plenum_radius = ${vessel_internal_radius}
plenum_D_h = '${fparse 2 * plenum_radius}'
plenum_area = '${fparse pi * plenum_radius * plenum_radius}'

plenum_n_elems = 5

PRI_L_plenum_inlet = 0.2
PRI_L_plenum_outlet = 0.2

##################################################################
###################      HX PARAMETERS     #######################
##################################################################

# hx parameters

hx_dia_pri = 0.003
hx_dia_sec = 0.005
hx_wall_thickness = 0.001
hx_length = 2
hx_n_elems_axial = 10
hx_nb_channels = 20000

##################################################################
###################    PIPES PARAMETERS    #######################
##################################################################

pri_pipes_radius = 0.12
pri_pipes_area = '${fparse pi * pri_pipes_radius * pri_pipes_radius}'
pri_pipes_D_h = '${fparse 2 * pi * pri_pipes_radius}'

pri_pipe1_n_elems = 20
pri_pipe2_n_elems = 10
pri_pipe3_n_elems = 10
pri_pipe4_n_elems = 10
pri_pipe5_n_elems = 3
pri_pipe_prz_n_elems = 10

PRI_L_upcomer_out = 0.5
PRI_L_upcomer = ${core_length_channel}
PRI_L2 = 2.
PRI_L3 = 2.
PRI_L4 = 2.
PRI_L5 = 2.
PRI_L_prz = 2.
PRI_L1 = '${fparse PRI_L_upcomer_out + PRI_L_plenum_inlet + core_length_channel + PRI_L_plenum_outlet + PRI_L2 + PRI_L3 + hx_length + PRI_L4 + PRI_L5 - PRI_L_upcomer}'

##################################################################
###################          PUMP          #######################
##################################################################

pump_area = '${fparse 1 * pri_pipes_area}'
pump_volume = 0.5
pump_inertia_coeff = '1 1 1 1'
pump_inertia_const = 1.61397
pump_omega_rated = 5
pump_speed_cr_I = 1e12
pump_speed_cr_fr = 0
pump_torque_rated = 50
pump_volumetric_rated = 1.9
pump_head_rated = 400
pump_tau_fr_coeff = '0 0 9.084 0'
pump_tau_fr_const = 0
pump_density_rated = 4.5

pri_motor_inertia = 2
pri_motor_torque = 50

shaft_initial_speed = 2

##################################################################
##################  GEOMETRICAL PARAMETERS  ######################
##################################################################

#
#
#                                                                                         ^
#                                                                                         |
#                                                                                         |
#                                                                                         | PRESSURIZER
#                                                                                         |                   HX [SEC]
#                                                                                         |           [...]<=============[...]
#               UPCOMER(out)  PLENUM(inlet)       CORE      PLENUM(outlet)     PIPE2      |    PIPE3          HX [PRI]            PIPE4                 PIPE5
#   (x0,y0)    X------------> -------------> =============> -------------> --------------> --------------> =============> --------------> [PUMP] -------------->
#                                            <------------- <---------------------------------------------------------------------------------------------------
#                                                UPCOMER                                                     PIPE1
#
# Remark: the upcomer is linked by a junction to the "upcomer_out" pipe, even if they are not exactly geometrically joined (it occurs because the upcomer must have the ame length than the core).

pri_x0 = '${fparse PRI_L1 + PRI_L_upcomer + PRI_L_upcomer_out + PRI_L_plenum_inlet}'
pri_y0 = 0.

pri_x_pipe1 = ${pri_x0}
pri_x_upcomer = '${fparse pri_x_pipe1 - PRI_L1}'
pri_x_upcomer_out = '${fparse pri_x_upcomer - PRI_L_upcomer - PRI_L_upcomer_out - PRI_L_plenum_inlet}'
pri_x_plenum_inlet = '${fparse pri_x_upcomer_out + PRI_L_upcomer_out}'
pri_x_core = '${fparse pri_x_plenum_inlet + PRI_L_plenum_inlet}'
pri_x_plenum_outlet = '${fparse pri_x_core + core_length_channel}'
pri_x_pipe2 = '${fparse pri_x_plenum_outlet + PRI_L_plenum_outlet}'
pri_x_pipe3 = '${fparse pri_x_pipe2 + PRI_L2}'
pri_x_hx = '${fparse pri_x_pipe3 + PRI_L3}'
pri_x_pipe4 = '${fparse pri_x_hx + hx_length}'
pri_x_pipe5 = '${fparse pri_x_pipe4 + PRI_L4}'

pri_y_pipe1 = ${pri_y0}
pri_y_upcomer = ${pri_y_pipe1}
pri_y_upcomer_out = '${fparse pri_y_pipe1}'
pri_y_plenum_inlet = ${pri_y_upcomer_out}
pri_y_core = ${pri_y_upcomer_out}
pri_y_plenum_outlet = ${pri_y_upcomer_out}
pri_y_pipe2 = ${pri_y_upcomer_out}
pri_y_pipe3 = ${pri_y_upcomer_out}
pri_y_hx = ${pri_y_upcomer_out}
pri_y_pipe4 = ${pri_y_upcomer_out}
pri_y_pipe5 = ${pri_y_upcomer_out}

################################################################################################################
################  ********************    SECONDARY LOOP PARAMETERS     ********************  ##################
################################################################################################################

sec_x_hx = '${fparse pri_x_hx + hx_length}'
sec_y_hx = '${fparse pri_y_hx + hx_wall_thickness}'

##########################
##########################
####                  ####
####        3         ####
####                  ####
####  GLOBAL PARAMS   ####
####                  ####
####                  ####
##########################
##########################

[GlobalParams]

  # initial_vel = ${vel_ini_default}
  # initial_vel_x = ${vel_x_ini_default}
  initial_vel_y = ${vel_y_ini_default}
  initial_vel_z = ${vel_z_ini_default}

  initial_T = ${T_ini}

  closures = no_closures
  scaling_factor_1phase = '1 1e-3 1e-5'

  gravity_vector = '0 0 0'

  scaling_factor_rhoV = 1
  scaling_factor_rhouV = 1e-2
  scaling_factor_rhovV = 1e-2
  scaling_factor_rhowV = 1e-3
  scaling_factor_rhoEV = 1e-5
  scaling_factor_temperature = 1e-2
  rdg_slope_reconstruction = full
[]

##########################
##########################
####                  ####
####        4         ####
####                  ####
####    PROPERTIES    ####
####                  ####
####                  ####
##########################
##########################

[FluidProperties]
  [he]
    type = IdealGasFluidProperties
    molar_mass = 4e-3
    gamma = 1.67
    k = 0.2556
    mu = 3.22639e-5
  []
  [air]
    type = IdealGasFluidProperties
    molar_mass = 29e-3
    gamma = 1.4
    k = 0.025 # W/(m.K)
    mu = 1.8e-5 # Pa.s
  []
[]

[Closures]
  [no_closures]
    type = Closures1PhaseNone
  []
[]

[Materials]
  [f]
    type = ADWallFrictionChurchillMaterial
    block = 'upcomer upcomer_out plenum_inlet plenum_outlet pri_pipe1 pri_pipe2 pri_pipe3 pri_pipe4 pri_pipe5 pressu/pipe_prz hx/pri hx/sec
    core/coolant_channel_r14_c1
    core/coolant_channel_r14_c3
    core/coolant_channel_r13_c1
    core/coolant_channel_r13_c2
    core/coolant_channel_r13_c3
    core/coolant_channel_r12_c0
    core/coolant_channel_r12_c1
    core/coolant_channel_r12_c3
    core/coolant_channel_r12_c4
    core/coolant_channel_r11_c1
    core/coolant_channel_r11_c2
    core/coolant_channel_r11_c3
    core/coolant_channel_r10_c0
    core/coolant_channel_r10_c1
    core/coolant_channel_r10_c3
    core/coolant_channel_r10_c4
    core/coolant_channel_r9_c0
    core/coolant_channel_r9_c1
    core/coolant_channel_r9_c2
    core/coolant_channel_r9_c3
    core/coolant_channel_r9_c4
    core/coolant_channel_r8_c0
    core/coolant_channel_r8_c1
    core/coolant_channel_r8_c3
    core/coolant_channel_r8_c4
    core/coolant_channel_r7_c0
    core/coolant_channel_r7_c1
    core/coolant_channel_r7_c2
    core/coolant_channel_r7_c3
    core/coolant_channel_r7_c4
    core/coolant_channel_r6_c0
    core/coolant_channel_r6_c1
    core/coolant_channel_r6_c3
    core/coolant_channel_r6_c4
    core/coolant_channel_r5_c0
    core/coolant_channel_r5_c1
    core/coolant_channel_r5_c2
    core/coolant_channel_r5_c3
    core/coolant_channel_r5_c4
    core/coolant_channel_r4_c0
    core/coolant_channel_r4_c1
    core/coolant_channel_r4_c3
    core/coolant_channel_r4_c4
    core/coolant_channel_r3_c1
    core/coolant_channel_r3_c2
    core/coolant_channel_r3_c3
    core/coolant_channel_r2_c0
    core/coolant_channel_r2_c1
    core/coolant_channel_r2_c3
    core/coolant_channel_r2_c4
    core/coolant_channel_r1_c1
    core/coolant_channel_r1_c2
    core/coolant_channel_r1_c3
    core/coolant_channel_r0_c1
    core/coolant_channel_r0_c3'
    D_h = D_h
    f_D = f_D
    mu = mu
    rho = rho
    vel = vel
    roughness = 3e-6
  []

  [Hw_core_and_hx]
    type = ADWallHeatTransferCoefficient3EqnDittusBoelterMaterial
    block = 'hx/pri hx/sec
    core/coolant_channel_r14_c1
    core/coolant_channel_r14_c3
    core/coolant_channel_r13_c1
    core/coolant_channel_r13_c2
    core/coolant_channel_r13_c3
    core/coolant_channel_r12_c0
    core/coolant_channel_r12_c1
    core/coolant_channel_r12_c3
    core/coolant_channel_r12_c4
    core/coolant_channel_r11_c1
    core/coolant_channel_r11_c2
    core/coolant_channel_r11_c3
    core/coolant_channel_r10_c0
    core/coolant_channel_r10_c1
    core/coolant_channel_r10_c3
    core/coolant_channel_r10_c4
    core/coolant_channel_r9_c0
    core/coolant_channel_r9_c1
    core/coolant_channel_r9_c2
    core/coolant_channel_r9_c3
    core/coolant_channel_r9_c4
    core/coolant_channel_r8_c0
    core/coolant_channel_r8_c1
    core/coolant_channel_r8_c3
    core/coolant_channel_r8_c4
    core/coolant_channel_r7_c0
    core/coolant_channel_r7_c1
    core/coolant_channel_r7_c2
    core/coolant_channel_r7_c3
    core/coolant_channel_r7_c4
    core/coolant_channel_r6_c0
    core/coolant_channel_r6_c1
    core/coolant_channel_r6_c3
    core/coolant_channel_r6_c4
    core/coolant_channel_r5_c0
    core/coolant_channel_r5_c1
    core/coolant_channel_r5_c2
    core/coolant_channel_r5_c3
    core/coolant_channel_r5_c4
    core/coolant_channel_r4_c0
    core/coolant_channel_r4_c1
    core/coolant_channel_r4_c3
    core/coolant_channel_r4_c4
    core/coolant_channel_r3_c1
    core/coolant_channel_r3_c2
    core/coolant_channel_r3_c3
    core/coolant_channel_r2_c0
    core/coolant_channel_r2_c1
    core/coolant_channel_r2_c3
    core/coolant_channel_r2_c4
    core/coolant_channel_r1_c1
    core/coolant_channel_r1_c2
    core/coolant_channel_r1_c3
    core/coolant_channel_r0_c1
    core/coolant_channel_r0_c3'
    D_h = D_h
    rho = rho
    vel = vel
    T = T
    T_wall = T_wall
    cp = cp
    mu = mu
    k = k
    Hw = Hw
  []

  [Hw_reflector]
    type = ADWallHeatTransferCoefficient3EqnDittusBoelterMaterial
    block = 'upcomer'
    D_h = D_h
    rho = rho
    vel = vel
    T = T
    T_wall = T_wall_reflector
    cp = cp
    mu = mu
    k = k
    Hw = Hw:1
  []
  [Hw_vessel]
    type = ADWallHeatTransferCoefficient3EqnDittusBoelterMaterial
    block = 'upcomer'
    D_h = D_h
    rho = rho
    vel = vel
    T = T
    T_wall = T_wall_vessel
    cp = cp
    mu = mu
    k = k
    Hw = Hw:2
  []
  [Hw_average]
    type = ADWallHeatTransferCoefficient3EqnDittusBoelterMaterial
    block = 'upcomer'
    D_h = D_h
    rho = rho
    vel = vel
    T = T
    T_wall = T_wall_vessel
    cp = cp
    mu = mu
    k = k
    Hw = Hw
  []

  [T_wall_reflector]
    type = ADCoupledVariableValueMaterial
    coupled_variable = T_wall_reflector_var
    prop_name = T_wall_reflector
    block = 'upcomer'
  []

  [T_wall_vessel]
    type = ADCoupledVariableValueMaterial
    coupled_variable = T_wall_vessel_var
    prop_name = T_wall_vessel
    block = 'upcomer'
  []
[]

[HeatStructureMaterials]
  [graphite]
    type = SolidMaterialProperties
    rho = 2160 # kg/m3
    k = 40 # W/(m.K)
    cp = 2100 # J/(kg.K)         approximate mean specific heat of graphite between 800 K (coolant) and 1400 K (fuel)
  []
  [fuel]
    type = SolidMaterialProperties
    rho = 10970 # kg/m3
    k = 5 # W/(m.K)
    cp = 300 # J/(kg.K)
  []
  [steel]
    type = SolidMaterialProperties
    rho = 8050
    k = 45
    cp = 466
  []
[]

##########################
##########################
####                  ####
####        5         ####
####                  ####
####    FUNCTIONS     ####
####                  ####
####                  ####
##########################
##########################

[Functions]

  ################################################################################################################
  ################  ********************      PRIMARY LOOP FUNCTIONS      ********************  ##################
  ################################################################################################################

  [head_fcn]
    type = PiecewiseLinear
    data_file = bingham_head_data.csv
    format = columns
  []
  [torque_fcn]
    type = PiecewiseLinear
    data_file = bingham_torque_data.csv
    format = columns
  []
[]

##########################
##########################
####                  ####
####        6         ####
####                  ####
####    COMPONENTS    ####
####                  ####
####                  ####
##########################
##########################

[Components]

  ################################################################################################################
  ################  ********************      PRIMARY LOOP COMPONENTS     ********************  ##################
  ################################################################################################################

  [pri_pipe1]
    type = FlowChannel1Phase
    position = '${pri_x_pipe1} ${pri_y_pipe1} 0.'
    orientation = '-1 0 0'
    length = ${PRI_L1}
    n_elems = ${pri_pipe1_n_elems}
    A = '${pri_pipes_area}'
    D_h = '${pri_pipes_D_h}'

    fp = he
    initial_p = ${pri_press}
    initial_vel = ${vel_ini_pri}
  []

  [pri_jct_1_upcomer]
    type = JunctionParallelChannels1Phase
    position = '${pri_x_upcomer} ${pri_y_upcomer} 0.'
    volume = ${volume_jct}
    initial_p = ${pri_press}
    initial_vel_x = ${vel_ini_pri}
    connections = 'pri_pipe1:out upcomer:in'
  []

  [upcomer]
    type = FlowChannel1Phase
    position = '${pri_x_upcomer} ${pri_y_upcomer} 0.'
    orientation = '-1 0 0'
    length = ${PRI_L_upcomer}
    n_elems = ${upcomer_n_elems}
    A = '${upcomer_area}'
    D_h = '${upcomer_D_h}'

    fp = he
    initial_p = ${pri_press}
    initial_vel = ${vel_ini_pri}
  []

  [upcomer_ht_outer_shield]
    type = HeatTransferFromExternalAppTemperature1Phase
    flow_channel = upcomer
    initial_T_wall = ${T_ini}
    P_hf = ${upcomer_P_hf_in}
    T_ext = T_wall_reflector
    #   Hw = Hw_reflector
  []
  [upcomer_ht_vessel]
    type = HeatTransferFromExternalAppTemperature1Phase
    flow_channel = upcomer
    initial_T_wall = ${T_ini}
    P_hf = ${upcomer_P_hf_out}
    T_ext = T_wall_vessel
    # Hw = Hw_vessel
  []

  [pri_jct_upc_out]
    type = JunctionOneToOne1Phase
    connections = 'upcomer:out upcomer_out:in'
  []

  [upcomer_out]
    type = FlowChannel1Phase
    position = '${pri_x_upcomer_out} ${pri_y_upcomer_out} 0.'
    orientation = '1 0 0'
    length = ${PRI_L_upcomer_out}
    n_elems = ${upcomer_out_n_elems}
    A = '${upcomer_area}'
    D_h = '${upcomer_D_h}'

    fp = he
    initial_p = ${pri_press}
    initial_vel = ${vel_ini_pri}
  []

  [pri_jct_upc_out_plm_in]
    type = JunctionParallelChannels1Phase
    position = '${pri_x_upcomer_out} ${pri_y_upcomer_out} 0.'
    volume = ${volume_jct}
    initial_p = ${pri_press}
    initial_vel_x = ${vel_ini_pri}
    connections = 'upcomer_out:out plenum_inlet:in'
  []

  [plenum_inlet]
    type = FlowChannel1Phase
    position = '${pri_x_plenum_inlet} ${pri_y_plenum_inlet} 0.'
    orientation = '1 0 0'
    length = ${PRI_L_plenum_inlet}
    n_elems = ${plenum_n_elems}
    A = '${plenum_area}'
    D_h = '${plenum_D_h}'

    fp = he
    initial_p = ${pri_press}
    initial_vel = ${vel_ini_pri}
  []

  [pri_jct_plenum_inlet_core]
    type = JunctionParallelChannels1Phase
    position = '${pri_x_core} ${pri_y_core} 0.'
    volume = ${volume_jct}
    initial_p = ${pri_press}
    initial_vel_x = ${vel_ini_pri}
    connections = 'plenum_inlet:out
    core/coolant_channel_r14_c1:in
    core/coolant_channel_r14_c3:in
    core/coolant_channel_r13_c1:in
    core/coolant_channel_r13_c2:in
    core/coolant_channel_r13_c3:in
    core/coolant_channel_r12_c0:in
    core/coolant_channel_r12_c1:in
    core/coolant_channel_r12_c3:in
    core/coolant_channel_r12_c4:in
    core/coolant_channel_r11_c1:in
    core/coolant_channel_r11_c2:in
    core/coolant_channel_r11_c3:in
    core/coolant_channel_r10_c0:in
    core/coolant_channel_r10_c1:in
    core/coolant_channel_r10_c3:in
    core/coolant_channel_r10_c4:in
    core/coolant_channel_r9_c0:in
    core/coolant_channel_r9_c1:in
    core/coolant_channel_r9_c2:in
    core/coolant_channel_r9_c3:in
    core/coolant_channel_r9_c4:in
    core/coolant_channel_r8_c0:in
    core/coolant_channel_r8_c1:in
    core/coolant_channel_r8_c3:in
    core/coolant_channel_r8_c4:in
    core/coolant_channel_r7_c0:in
    core/coolant_channel_r7_c1:in
    core/coolant_channel_r7_c2:in
    core/coolant_channel_r7_c3:in
    core/coolant_channel_r7_c4:in
    core/coolant_channel_r6_c0:in
    core/coolant_channel_r6_c1:in
    core/coolant_channel_r6_c3:in
    core/coolant_channel_r6_c4:in
    core/coolant_channel_r5_c0:in
    core/coolant_channel_r5_c1:in
    core/coolant_channel_r5_c2:in
    core/coolant_channel_r5_c3:in
    core/coolant_channel_r5_c4:in
    core/coolant_channel_r4_c0:in
    core/coolant_channel_r4_c1:in
    core/coolant_channel_r4_c3:in
    core/coolant_channel_r4_c4:in
    core/coolant_channel_r3_c1:in
    core/coolant_channel_r3_c2:in
    core/coolant_channel_r3_c3:in
    core/coolant_channel_r2_c0:in
    core/coolant_channel_r2_c1:in
    core/coolant_channel_r2_c3:in
    core/coolant_channel_r2_c4:in
    core/coolant_channel_r1_c1:in
    core/coolant_channel_r1_c2:in
    core/coolant_channel_r1_c3:in
    core/coolant_channel_r0_c1:in
    core/coolant_channel_r0_c3:in'
  []

  [core]

    #####################################
    #######        RANK 14        #######
    #####################################

    [coolant_channel_r14_c1]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank14_y} ${even_rank_column1_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r14_c1]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r14_c1
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r14_c3]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank14_y} ${even_rank_column3_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r14_c3]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r14_c3
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    #####################################
    #######        RANK 13        #######
    #####################################

    [coolant_channel_r13_c1]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank13_y} ${odd_rank_column1_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r13_c1]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r13_c1
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r13_c2]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank13_y} ${odd_rank_column2_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r13_c2]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r13_c2
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r13_c3]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank13_y} ${odd_rank_column3_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r13_c3]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r13_c3
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    #####################################
    #######        RANK 12        #######
    #####################################

    [coolant_channel_r12_c0]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank12_y} ${even_rank_column0_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r12_c0]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r12_c0
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r12_c1]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank12_y} ${even_rank_column1_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r12_c1]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r12_c1
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r12_c3]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank12_y} ${even_rank_column3_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r12_c3]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r12_c3
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r12_c4]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank12_y} ${even_rank_column4_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r12_c4]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r12_c4
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    #####################################
    #######        RANK 11         #######
    #####################################

    [coolant_channel_r11_c1]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank11_y} ${odd_rank_column1_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r11_c1]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r11_c1
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r11_c2]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank11_y} ${odd_rank_column2_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r11_c2]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r11_c2
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r11_c3]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank11_y} ${odd_rank_column3_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r11_c3]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r11_c3
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    #####################################
    #######        RANK 10        #######
    #####################################

    [coolant_channel_r10_c0]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank10_y} ${even_rank_column0_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r10_c0]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r10_c0
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r10_c1]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank10_y} ${even_rank_column1_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r10_c1]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r10_c1
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r10_c3]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank10_y} ${even_rank_column3_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r10_c3]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r10_c3
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r10_c4]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank10_y} ${even_rank_column4_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r10_c4]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r10_c4
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    #####################################
    #######        RANK 9         #######
    #####################################

    [coolant_channel_r9_c0]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank9_y} ${odd_rank_column0_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r9_c0]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r9_c0
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r9_c1]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank9_y} ${odd_rank_column1_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r9_c1]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r9_c1
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r9_c2]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank9_y} ${odd_rank_column2_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r9_c2]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r9_c2
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r9_c3]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank9_y} ${odd_rank_column3_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r9_c3]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r9_c3
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r9_c4]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank9_y} ${odd_rank_column4_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r9_c4]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r9_c4
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    #####################################
    #######        RANK 8         #######
    #####################################

    [coolant_channel_r8_c0]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank8_y} ${even_rank_column0_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r8_c0]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r8_c0
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r8_c1]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank8_y} ${even_rank_column1_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r8_c1]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r8_c1
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r8_c3]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank8_y} ${even_rank_column3_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r8_c3]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r8_c3
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r8_c4]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank8_y} ${even_rank_column4_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r8_c4]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r8_c4
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    #####################################
    #######        RANK 7         #######
    #####################################

    [coolant_channel_r7_c0]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank7_y} ${odd_rank_column0_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r7_c0]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r7_c0
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r7_c1]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank7_y} ${odd_rank_column1_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r7_c1]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r7_c1
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r7_c2]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank7_y} ${odd_rank_column2_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r7_c2]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r7_c2
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r7_c3]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank7_y} ${odd_rank_column3_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r7_c3]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r7_c3
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r7_c4]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank7_y} ${odd_rank_column4_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r7_c4]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r7_c4
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    #####################################
    #######        RANK 6         #######
    #####################################

    [coolant_channel_r6_c0]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank6_y} ${even_rank_column0_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r6_c0]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r6_c0
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r6_c1]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank6_y} ${even_rank_column1_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r6_c1]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r6_c1
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r6_c3]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank6_y} ${even_rank_column3_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r6_c3]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r6_c3
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r6_c4]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank6_y} ${even_rank_column4_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r6_c4]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r6_c4
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    #####################################
    #######        RANK 5         #######
    #####################################

    [coolant_channel_r5_c0]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank5_y} ${odd_rank_column0_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r5_c0]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r5_c0
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r5_c1]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank5_y} ${odd_rank_column1_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r5_c1]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r5_c1
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r5_c2]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank5_y} ${odd_rank_column2_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r5_c2]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r5_c2
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r5_c3]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank5_y} ${odd_rank_column3_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r5_c3]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r5_c3
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r5_c4]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank5_y} ${odd_rank_column4_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r5_c4]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r5_c4
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    #####################################
    #######        RANK 4         #######
    #####################################

    [coolant_channel_r4_c0]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank4_y} ${even_rank_column0_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r4_c0]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r4_c0
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r4_c1]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank4_y} ${even_rank_column1_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r4_c1]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r4_c1
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r4_c3]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank4_y} ${even_rank_column3_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r4_c3]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r4_c3
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r4_c4]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank4_y} ${even_rank_column4_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r4_c4]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r4_c4
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    #####################################
    #######        RANK 3         #######
    #####################################

    [coolant_channel_r3_c1]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank3_y} ${odd_rank_column1_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r3_c1]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r3_c1
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r3_c2]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank3_y} ${odd_rank_column2_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r3_c2]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r3_c2
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r3_c3]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank3_y} ${odd_rank_column3_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r3_c3]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r3_c3
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    #####################################
    #######        RANK 2         #######
    #####################################

    [coolant_channel_r2_c0]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank2_y} ${even_rank_column0_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r2_c0]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r2_c0
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r2_c1]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank2_y} ${even_rank_column1_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r2_c1]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r2_c1
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r2_c3]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank2_y} ${even_rank_column3_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r2_c3]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r2_c3
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r2_c4]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank2_y} ${even_rank_column4_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r2_c4]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r2_c4
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    #####################################
    #######        RANK 1         #######
    #####################################

    [coolant_channel_r1_c1]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank1_y} ${odd_rank_column1_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r1_c1]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r1_c1
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r1_c2]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank1_y} ${odd_rank_column2_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r1_c2]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r1_c2
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r1_c3]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank1_y} ${odd_rank_column3_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r1_c3]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r1_c3
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    #####################################
    #######        RANK 0         #######
    #####################################

    [coolant_channel_r0_c1]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank0_y} ${even_rank_column1_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r0_c1]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r0_c1
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

    [coolant_channel_r0_c3]
      type = FlowChannel1Phase
      position = '${pri_x_core} ${rank0_y} ${even_rank_column3_z}'
      orientation = '1 0 0'
      length = ${core_length_channel}
      n_elems = ${core_channel_n_elems}
      A = '${fparse pi * core_nb_coolant_per_assembly * core_radius_coolant * core_radius_coolant}'
      D_h = '${fparse 2 * core_radius_coolant}'
      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [core_ht_r0_c3]
      type = HeatTransferFromExternalAppTemperature1Phase
      flow_channel = core/coolant_channel_r0_c3
      P_hf = '${fparse core_nb_coolant_per_assembly * 2 * pi * core_radius_coolant}'
      initial_T_wall = ${T_ini}
      # var_type = elemental

    []

  []

  [pri_jct_core_plenum_out]
    type = JunctionParallelChannels1Phase
    position = '${pri_x_pipe2} ${pri_y_pipe2} 0.'
    volume = ${volume_jct}
    initial_p = ${pri_press}
    initial_vel_x = ${vel_ini_pri}
    connections = 'core/coolant_channel_r14_c1:out
    core/coolant_channel_r14_c3:out
    core/coolant_channel_r13_c1:out
    core/coolant_channel_r13_c2:out
    core/coolant_channel_r13_c3:out
    core/coolant_channel_r12_c0:out
    core/coolant_channel_r12_c1:out
    core/coolant_channel_r12_c3:out
    core/coolant_channel_r12_c4:out
    core/coolant_channel_r11_c1:out
    core/coolant_channel_r11_c2:out
    core/coolant_channel_r11_c3:out
    core/coolant_channel_r10_c0:out
    core/coolant_channel_r10_c1:out
    core/coolant_channel_r10_c3:out
    core/coolant_channel_r10_c4:out
    core/coolant_channel_r9_c0:out
    core/coolant_channel_r9_c1:out
    core/coolant_channel_r9_c2:out
    core/coolant_channel_r9_c3:out
    core/coolant_channel_r9_c4:out
    core/coolant_channel_r8_c0:out
    core/coolant_channel_r8_c1:out
    core/coolant_channel_r8_c3:out
    core/coolant_channel_r8_c4:out
    core/coolant_channel_r7_c0:out
    core/coolant_channel_r7_c1:out
    core/coolant_channel_r7_c2:out
    core/coolant_channel_r7_c3:out
    core/coolant_channel_r7_c4:out
    core/coolant_channel_r6_c0:out
    core/coolant_channel_r6_c1:out
    core/coolant_channel_r6_c3:out
    core/coolant_channel_r6_c4:out
    core/coolant_channel_r5_c0:out
    core/coolant_channel_r5_c1:out
    core/coolant_channel_r5_c2:out
    core/coolant_channel_r5_c3:out
    core/coolant_channel_r5_c4:out
    core/coolant_channel_r4_c0:out
    core/coolant_channel_r4_c1:out
    core/coolant_channel_r4_c3:out
    core/coolant_channel_r4_c4:out
    core/coolant_channel_r3_c1:out
    core/coolant_channel_r3_c2:out
    core/coolant_channel_r3_c3:out
    core/coolant_channel_r2_c0:out
    core/coolant_channel_r2_c1:out
    core/coolant_channel_r2_c3:out
    core/coolant_channel_r2_c4:out
    core/coolant_channel_r1_c1:out
    core/coolant_channel_r1_c2:out
    core/coolant_channel_r1_c3:out
    core/coolant_channel_r0_c1:out
    core/coolant_channel_r0_c3:out
    plenum_outlet:in'
  []

  [plenum_outlet]
    type = FlowChannel1Phase
    position = '${pri_x_plenum_outlet} ${pri_y_plenum_outlet} 0.'
    orientation = '1 0 0'
    length = ${PRI_L_plenum_outlet}
    n_elems = ${plenum_n_elems}
    A = '${plenum_area}'
    D_h = '${plenum_D_h}'

    fp = he
    initial_p = ${pri_press}
    initial_vel = ${vel_ini_pri}
  []

  [pri_jct_plenum_out_2]
    type = JunctionParallelChannels1Phase
    position = '${pri_x_pipe2} ${pri_y_pipe2} 0.'
    volume = ${volume_jct}
    initial_p = ${pri_press}
    initial_vel_x = ${vel_ini_pri}
    connections = 'plenum_outlet:out pri_pipe2:in'
  []

  [pri_pipe2]
    type = FlowChannel1Phase
    position = '${pri_x_pipe2} ${pri_y_pipe2} 0.'
    orientation = '1 0 0'
    length = ${PRI_L2}
    n_elems = ${pri_pipe2_n_elems}
    A = '${pri_pipes_area}'
    D_h = '${pri_pipes_D_h}'

    fp = he
    initial_p = ${pri_press}
    initial_vel = ${vel_ini_pri}
  []

  [pressu]
    [pri_jct_2_3_prz]
      type = VolumeJunction1Phase
      connections = 'pri_pipe2:out pri_pipe3:in pressu/pipe_prz:in'
      position = '${pri_x_pipe3} ${pri_y_pipe3} 0.'
      volume = ${volume_jct}
      initial_p = ${pri_press}
      initial_vel_x = ${vel_ini_pri}
    []
    [pipe_prz]
      type = FlowChannel1Phase
      position = '${pri_x_pipe3} ${pri_y_pipe3} 0.'
      orientation = '0 1 0'
      length = ${PRI_L_prz}
      n_elems = ${pri_pipe_prz_n_elems}
      A = '${pri_pipes_area}'
      D_h = '${pri_pipes_D_h}'

      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []
    [prz]
      type = InletStagnationPressureTemperature1Phase
      p0 = ${pri_press}
      T0 = ${T_ini}
      input = 'pressu/pipe_prz:out'
    []
  []

  [pri_pipe3]
    type = FlowChannel1Phase
    position = '${pri_x_pipe3} ${pri_y_pipe3} 0.'
    orientation = '1 0 0'
    length = ${PRI_L3}
    n_elems = ${pri_pipe3_n_elems}
    A = '${pri_pipes_area}'
    D_h = '${pri_pipes_D_h}'

    fp = he
    initial_p = ${pri_press}
    initial_vel = ${vel_ini_pri}
  []

  [pri_jct_3_hx_pri]
    type = JunctionParallelChannels1Phase
    position = '${pri_x_hx} ${pri_y_hx} 0.'
    connections = 'pri_pipe3:out hx/pri:in'
    volume = ${volume_jct}
    initial_p = ${pri_press}
    initial_vel_x = ${vel_ini_pri}
  []

  [hx]
    [pri]
      type = FlowChannel1Phase
      position = '${pri_x_hx} ${pri_y_hx} 0.'
      orientation = '1 0 0'
      length = ${hx_length}
      n_elems = ${hx_n_elems_axial}
      A = '${fparse pi * hx_nb_channels * hx_dia_pri * hx_dia_pri * 0.25}'
      D_h = ${hx_dia_pri}

      fp = he
      initial_p = ${pri_press}
      initial_vel = ${vel_ini_pri}
    []

    [ht_pri]
      type = HeatTransferFromHeatStructure1Phase
      hs = hx/wall
      hs_side = inner
      flow_channel = hx/pri
      P_hf = '${fparse pi * hx_nb_channels * hx_dia_pri}'
    []

    [wall]
      type = HeatStructureCylindrical
      length = ${hx_length}
      n_elems = ${hx_n_elems_axial}
      n_part_elems = 1
      names = 'hx_wall'
      orientation = '1 0 0'
      position = '${pri_x_hx} ${pri_y_hx} 0.'
      widths = '${hx_wall_thickness}'
      materials = 'steel'
      inner_radius = '${fparse hx_dia_pri / 2}'
      num_rods = ${hx_nb_channels}
    []

    [ht_sec]
      type = HeatTransferFromHeatStructure1Phase
      hs = hx/wall
      hs_side = outer
      flow_channel = hx/sec
      P_hf = '${fparse pi * hx_nb_channels * hx_dia_sec}'
    []

    [sec]
      type = FlowChannel1Phase
      position = '${sec_x_hx} ${sec_y_hx} 0.'
      orientation = '-1 0 0'
      length = ${hx_length}
      n_elems = ${hx_n_elems_axial}
      A = '${fparse pi * hx_nb_channels * hx_dia_sec * hx_dia_sec * 0.25}'
      D_h = '${hx_dia_sec}'

      fp = air
      initial_p = ${p_sec}
      initial_vel = ${vel_ini_sec}
    []
  []

  [hx_sec_inlet]
    type = InletMassFlowRateTemperature1Phase
    input = hx/sec:in
    m_dot = 20.5
    T = 470
  []

  [hx_sec_outlet]
    type = Outlet1Phase
    input = hx/sec:out
    p = 5.5e5
  []

  [pri_jct_hx_pri_4]
    type = JunctionParallelChannels1Phase
    position = '${pri_x_pipe4} ${pri_y_pipe4} 0.'
    connections = 'hx/pri:out pri_pipe4:in'
    volume = ${volume_jct}
    initial_p = ${pri_press}
    initial_vel_x = ${vel_ini_pri}
  []

  [pri_pipe4]
    type = FlowChannel1Phase
    position = '${pri_x_pipe4} ${pri_y_pipe4} 0.'
    orientation = '1 0 0'
    length = ${PRI_L4}
    n_elems = ${pri_pipe4_n_elems}
    A = '${pri_pipes_area}'
    D_h = '${pri_pipes_D_h}'

    fp = he
    initial_p = ${pri_press}
    initial_vel = ${vel_ini_pri}
  []

  # [outlet]
  #   type = Outlet1Phase
  #   input = pri_pipe4:out
  #   p = 90e5
  # []

  # [inlet]
  #   type = InletMassFlowRateTemperature1Phase
  #   input = pri_pipe5:in
  #   m_dot = 9.4
  #   T = 890
  # []

  [circ]
    [pump]
      type = ShaftConnectedPump1Phase
      inlet = 'pri_pipe4:out'
      outlet = 'pri_pipe5:in'
      position = '${pri_x_pipe5} ${pri_y_pipe5} 0.'
      A_ref = ${pump_area}
      scaling_factor_rhoEV = 1e-5
      volume = ${pump_volume}
      inertia_coeff = ${pump_inertia_coeff}
      inertia_const = ${pump_inertia_const}
      omega_rated = ${pump_omega_rated}
      speed_cr_I = ${pump_speed_cr_I}
      speed_cr_fr = ${pump_speed_cr_fr}
      torque_rated = ${pump_torque_rated}
      volumetric_rated = ${pump_volumetric_rated}
      head_rated = ${pump_head_rated}
      tau_fr_coeff = ${pump_tau_fr_coeff}
      tau_fr_const = ${pump_tau_fr_const}
      head = head_fcn
      torque_hydraulic = torque_fcn
      density_rated = ${pump_density_rated}
      initial_p = ${pri_press}
      initial_vel_x = ${vel_ini_pri}
    []

    [motor]
      type = ShaftConnectedMotor
      inertia = ${pri_motor_inertia}
      torque = ${pri_motor_torque}
    []

    [shaft]
      type = Shaft
      connected_components = 'circ/motor circ/pump'
      initial_speed = ${shaft_initial_speed}
    []
  []

  [pri_pipe5]
    type = FlowChannel1Phase
    position = '${pri_x_pipe5} ${pri_y_pipe5} 0.'
    orientation = '1 0 0'
    length = ${PRI_L5}
    n_elems = ${pri_pipe5_n_elems}
    A = '${pri_pipes_area}'
    D_h = '${pri_pipes_D_h}'

    fp = he
    initial_p = ${pri_press}
    initial_vel = ${vel_ini_pri}
  []

  [pri_jct_5_1]
    type = JunctionOneToOne1Phase
    connections = 'pri_pipe5:out pri_pipe1:in'
  []

[]

[AuxVariables]
  [Hw_core_and_hx]
    family = monomial
    order = constant
    block = 'core/coolant_channel_r14_c1
    core/coolant_channel_r14_c3
    core/coolant_channel_r13_c1
    core/coolant_channel_r13_c2
    core/coolant_channel_r13_c3
    core/coolant_channel_r12_c0
    core/coolant_channel_r12_c1
    core/coolant_channel_r12_c3
    core/coolant_channel_r12_c4
    core/coolant_channel_r11_c1
    core/coolant_channel_r11_c2
    core/coolant_channel_r11_c3
    core/coolant_channel_r10_c0
    core/coolant_channel_r10_c1
    core/coolant_channel_r10_c3
    core/coolant_channel_r10_c4
    core/coolant_channel_r9_c0
    core/coolant_channel_r9_c1
    core/coolant_channel_r9_c2
    core/coolant_channel_r9_c3
    core/coolant_channel_r9_c4
    core/coolant_channel_r8_c0
    core/coolant_channel_r8_c1
    core/coolant_channel_r8_c3
    core/coolant_channel_r8_c4
    core/coolant_channel_r7_c0
    core/coolant_channel_r7_c1
    core/coolant_channel_r7_c2
    core/coolant_channel_r7_c3
    core/coolant_channel_r7_c4
    core/coolant_channel_r6_c0
    core/coolant_channel_r6_c1
    core/coolant_channel_r6_c3
    core/coolant_channel_r6_c4
    core/coolant_channel_r5_c0
    core/coolant_channel_r5_c1
    core/coolant_channel_r5_c2
    core/coolant_channel_r5_c3
    core/coolant_channel_r5_c4
    core/coolant_channel_r4_c0
    core/coolant_channel_r4_c1
    core/coolant_channel_r4_c3
    core/coolant_channel_r4_c4
    core/coolant_channel_r3_c1
    core/coolant_channel_r3_c2
    core/coolant_channel_r3_c3
    core/coolant_channel_r2_c0
    core/coolant_channel_r2_c1
    core/coolant_channel_r2_c3
    core/coolant_channel_r2_c4
    core/coolant_channel_r1_c1
    core/coolant_channel_r1_c2
    core/coolant_channel_r1_c3
    core/coolant_channel_r0_c1
    core/coolant_channel_r0_c3'
  []

  [Hw_reflector]
    family = monomial
    order = constant
    block = 'upcomer'
  []

  [Hw_vessel]
    family = monomial
    order = constant
    block = 'upcomer'
  []

  [T_wall_reflector_var]
    family = monomial
    order = constant
    block = 'upcomer'
  []

  [T_wall_vessel_var]
    family = monomial
    order = constant
    block = 'upcomer'
  []

  [T_fluid_core_var]
    family = monomial
    order = constant
    block = 'core/coolant_channel_r14_c1
    core/coolant_channel_r14_c3
    core/coolant_channel_r13_c1
    core/coolant_channel_r13_c2
    core/coolant_channel_r13_c3
    core/coolant_channel_r12_c0
    core/coolant_channel_r12_c1
    core/coolant_channel_r12_c3
    core/coolant_channel_r12_c4
    core/coolant_channel_r11_c1
    core/coolant_channel_r11_c2
    core/coolant_channel_r11_c3
    core/coolant_channel_r10_c0
    core/coolant_channel_r10_c1
    core/coolant_channel_r10_c3
    core/coolant_channel_r10_c4
    core/coolant_channel_r9_c0
    core/coolant_channel_r9_c1
    core/coolant_channel_r9_c2
    core/coolant_channel_r9_c3
    core/coolant_channel_r9_c4
    core/coolant_channel_r8_c0
    core/coolant_channel_r8_c1
    core/coolant_channel_r8_c3
    core/coolant_channel_r8_c4
    core/coolant_channel_r7_c0
    core/coolant_channel_r7_c1
    core/coolant_channel_r7_c2
    core/coolant_channel_r7_c3
    core/coolant_channel_r7_c4
    core/coolant_channel_r6_c0
    core/coolant_channel_r6_c1
    core/coolant_channel_r6_c3
    core/coolant_channel_r6_c4
    core/coolant_channel_r5_c0
    core/coolant_channel_r5_c1
    core/coolant_channel_r5_c2
    core/coolant_channel_r5_c3
    core/coolant_channel_r5_c4
    core/coolant_channel_r4_c0
    core/coolant_channel_r4_c1
    core/coolant_channel_r4_c3
    core/coolant_channel_r4_c4
    core/coolant_channel_r3_c1
    core/coolant_channel_r3_c2
    core/coolant_channel_r3_c3
    core/coolant_channel_r2_c0
    core/coolant_channel_r2_c1
    core/coolant_channel_r2_c3
    core/coolant_channel_r2_c4
    core/coolant_channel_r1_c1
    core/coolant_channel_r1_c2
    core/coolant_channel_r1_c3
    core/coolant_channel_r0_c1
    core/coolant_channel_r0_c3'
  []

  [T_fluid_upcomer_var]
    family = monomial
    order = constant
    block = 'upcomer'
  []

[]

[AuxKernels]
  [Hw_core_ak]
    type = ADMaterialRealAux
    variable = Hw_core_and_hx
    property = 'Hw'
    block = 'core/coolant_channel_r14_c1
    core/coolant_channel_r14_c3
    core/coolant_channel_r13_c1
    core/coolant_channel_r13_c2
    core/coolant_channel_r13_c3
    core/coolant_channel_r12_c0
    core/coolant_channel_r12_c1
    core/coolant_channel_r12_c3
    core/coolant_channel_r12_c4
    core/coolant_channel_r11_c1
    core/coolant_channel_r11_c2
    core/coolant_channel_r11_c3
    core/coolant_channel_r10_c0
    core/coolant_channel_r10_c1
    core/coolant_channel_r10_c3
    core/coolant_channel_r10_c4
    core/coolant_channel_r9_c0
    core/coolant_channel_r9_c1
    core/coolant_channel_r9_c2
    core/coolant_channel_r9_c3
    core/coolant_channel_r9_c4
    core/coolant_channel_r8_c0
    core/coolant_channel_r8_c1
    core/coolant_channel_r8_c3
    core/coolant_channel_r8_c4
    core/coolant_channel_r7_c0
    core/coolant_channel_r7_c1
    core/coolant_channel_r7_c2
    core/coolant_channel_r7_c3
    core/coolant_channel_r7_c4
    core/coolant_channel_r6_c0
    core/coolant_channel_r6_c1
    core/coolant_channel_r6_c3
    core/coolant_channel_r6_c4
    core/coolant_channel_r5_c0
    core/coolant_channel_r5_c1
    core/coolant_channel_r5_c2
    core/coolant_channel_r5_c3
    core/coolant_channel_r5_c4
    core/coolant_channel_r4_c0
    core/coolant_channel_r4_c1
    core/coolant_channel_r4_c3
    core/coolant_channel_r4_c4
    core/coolant_channel_r3_c1
    core/coolant_channel_r3_c2
    core/coolant_channel_r3_c3
    core/coolant_channel_r2_c0
    core/coolant_channel_r2_c1
    core/coolant_channel_r2_c3
    core/coolant_channel_r2_c4
    core/coolant_channel_r1_c1
    core/coolant_channel_r1_c2
    core/coolant_channel_r1_c3
    core/coolant_channel_r0_c1
    core/coolant_channel_r0_c3'
  []

  [Hw_reflector_ak]
    type = ADMaterialRealAux
    variable = Hw_reflector
    property = 'Hw:1'
    block = 'upcomer'
  []

  [Hw_vessel_ak]
    type = ADMaterialRealAux
    variable = Hw_vessel
    property = 'Hw:2'
    block = 'upcomer'
  []

  [T_fluid_core_ak]
    type = CopyValueAux
    source = T
    variable = T_fluid_core_var
    block = 'core/coolant_channel_r14_c1
    core/coolant_channel_r14_c3
    core/coolant_channel_r13_c1
    core/coolant_channel_r13_c2
    core/coolant_channel_r13_c3
    core/coolant_channel_r12_c0
    core/coolant_channel_r12_c1
    core/coolant_channel_r12_c3
    core/coolant_channel_r12_c4
    core/coolant_channel_r11_c1
    core/coolant_channel_r11_c2
    core/coolant_channel_r11_c3
    core/coolant_channel_r10_c0
    core/coolant_channel_r10_c1
    core/coolant_channel_r10_c3
    core/coolant_channel_r10_c4
    core/coolant_channel_r9_c0
    core/coolant_channel_r9_c1
    core/coolant_channel_r9_c2
    core/coolant_channel_r9_c3
    core/coolant_channel_r9_c4
    core/coolant_channel_r8_c0
    core/coolant_channel_r8_c1
    core/coolant_channel_r8_c3
    core/coolant_channel_r8_c4
    core/coolant_channel_r7_c0
    core/coolant_channel_r7_c1
    core/coolant_channel_r7_c2
    core/coolant_channel_r7_c3
    core/coolant_channel_r7_c4
    core/coolant_channel_r6_c0
    core/coolant_channel_r6_c1
    core/coolant_channel_r6_c3
    core/coolant_channel_r6_c4
    core/coolant_channel_r5_c0
    core/coolant_channel_r5_c1
    core/coolant_channel_r5_c2
    core/coolant_channel_r5_c3
    core/coolant_channel_r5_c4
    core/coolant_channel_r4_c0
    core/coolant_channel_r4_c1
    core/coolant_channel_r4_c3
    core/coolant_channel_r4_c4
    core/coolant_channel_r3_c1
    core/coolant_channel_r3_c2
    core/coolant_channel_r3_c3
    core/coolant_channel_r2_c0
    core/coolant_channel_r2_c1
    core/coolant_channel_r2_c3
    core/coolant_channel_r2_c4
    core/coolant_channel_r1_c1
    core/coolant_channel_r1_c2
    core/coolant_channel_r1_c3
    core/coolant_channel_r0_c1
    core/coolant_channel_r0_c3'
  []

  [T_fluid_upcomer_ak]
    type = CopyValueAux
    source = T
    variable = T_fluid_upcomer_var
    block = 'upcomer'
  []

[]

##########################
##########################
####                  ####
####        7         ####
####                  ####
####  POSTPROCESSORS  ####
####                  ####
####                  ####
##########################
##########################

[Postprocessors]

  ################################################################################################################
  ################  ********************    POSTPROCESSORS FROM MAIN APP  ********************  ##################
  ################################################################################################################

  [T_wall_avg]
    type = ElementAverageValue
    block = 'core/coolant_channel_r14_c1
    core/coolant_channel_r14_c3
    core/coolant_channel_r13_c1
    core/coolant_channel_r13_c2
    core/coolant_channel_r13_c3
    core/coolant_channel_r12_c0
    core/coolant_channel_r12_c1
    core/coolant_channel_r12_c3
    core/coolant_channel_r12_c4
    core/coolant_channel_r11_c1
    core/coolant_channel_r11_c2
    core/coolant_channel_r11_c3
    core/coolant_channel_r10_c0
    core/coolant_channel_r10_c1
    core/coolant_channel_r10_c3
    core/coolant_channel_r10_c4
    core/coolant_channel_r9_c0
    core/coolant_channel_r9_c1
    core/coolant_channel_r9_c2
    core/coolant_channel_r9_c3
    core/coolant_channel_r9_c4
    core/coolant_channel_r8_c0
    core/coolant_channel_r8_c1
    core/coolant_channel_r8_c3
    core/coolant_channel_r8_c4
    core/coolant_channel_r7_c0
    core/coolant_channel_r7_c1
    core/coolant_channel_r7_c2
    core/coolant_channel_r7_c3
    core/coolant_channel_r7_c4
    core/coolant_channel_r6_c0
    core/coolant_channel_r6_c1
    core/coolant_channel_r6_c3
    core/coolant_channel_r6_c4
    core/coolant_channel_r5_c0
    core/coolant_channel_r5_c1
    core/coolant_channel_r5_c2
    core/coolant_channel_r5_c3
    core/coolant_channel_r5_c4
    core/coolant_channel_r4_c0
    core/coolant_channel_r4_c1
    core/coolant_channel_r4_c3
    core/coolant_channel_r4_c4
    core/coolant_channel_r3_c1
    core/coolant_channel_r3_c2
    core/coolant_channel_r3_c3
    core/coolant_channel_r2_c0
    core/coolant_channel_r2_c1
    core/coolant_channel_r2_c3
    core/coolant_channel_r2_c4
    core/coolant_channel_r1_c1
    core/coolant_channel_r1_c2
    core/coolant_channel_r1_c3
    core/coolant_channel_r0_c1
    core/coolant_channel_r0_c3'
    variable = T_wall
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [htc_avg]
    type = ElementAverageValue
    block = 'core/coolant_channel_r14_c1
    core/coolant_channel_r14_c3
    core/coolant_channel_r13_c1
    core/coolant_channel_r13_c2
    core/coolant_channel_r13_c3
    core/coolant_channel_r12_c0
    core/coolant_channel_r12_c1
    core/coolant_channel_r12_c3
    core/coolant_channel_r12_c4
    core/coolant_channel_r11_c1
    core/coolant_channel_r11_c2
    core/coolant_channel_r11_c3
    core/coolant_channel_r10_c0
    core/coolant_channel_r10_c1
    core/coolant_channel_r10_c3
    core/coolant_channel_r10_c4
    core/coolant_channel_r9_c0
    core/coolant_channel_r9_c1
    core/coolant_channel_r9_c2
    core/coolant_channel_r9_c3
    core/coolant_channel_r9_c4
    core/coolant_channel_r8_c0
    core/coolant_channel_r8_c1
    core/coolant_channel_r8_c3
    core/coolant_channel_r8_c4
    core/coolant_channel_r7_c0
    core/coolant_channel_r7_c1
    core/coolant_channel_r7_c2
    core/coolant_channel_r7_c3
    core/coolant_channel_r7_c4
    core/coolant_channel_r6_c0
    core/coolant_channel_r6_c1
    core/coolant_channel_r6_c3
    core/coolant_channel_r6_c4
    core/coolant_channel_r5_c0
    core/coolant_channel_r5_c1
    core/coolant_channel_r5_c2
    core/coolant_channel_r5_c3
    core/coolant_channel_r5_c4
    core/coolant_channel_r4_c0
    core/coolant_channel_r4_c1
    core/coolant_channel_r4_c3
    core/coolant_channel_r4_c4
    core/coolant_channel_r3_c1
    core/coolant_channel_r3_c2
    core/coolant_channel_r3_c3
    core/coolant_channel_r2_c0
    core/coolant_channel_r2_c1
    core/coolant_channel_r2_c3
    core/coolant_channel_r2_c4
    core/coolant_channel_r1_c1
    core/coolant_channel_r1_c2
    core/coolant_channel_r1_c3
    core/coolant_channel_r0_c1
    core/coolant_channel_r0_c3'
    variable = Hw_core_and_hx
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ################################################################################################################
  ################  ********************    PRIMARY LOOP POSTPROCESSORS   ********************  ##################
  ################################################################################################################

  #####################################################
  ######                  core                   ######
  #####################################################

  ###### pressure

  [core_p_in]
    type = SideAverageValue
    boundary = pri_jct_1_upcomer
    variable = p
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [core_p_out]
    type = SideAverageValue
    boundary = pri_jct_plenum_out_2
    variable = p
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [core_p_difference]
    type = DifferencePostprocessor
    value1 = core_p_out
    value2 = core_p_in
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ###### m_dot

  # [core_m_dot_in]
  #   type = ADFlowJunctionFlux1Phase
  #   boundary = pri_pipe1:in
  #   equation = mass
  #   junction = pri_jct_6_1
  #   connection_index = 1
  #   execute_on = 'INITIAL TIMESTEP_END'
  # []
  # [core_m_dot_out]
  #   type = ADFlowJunctionFlux1Phase
  #   boundary = pri_pipe3:out
  #   equation = mass
  #   junction = pri_jct_3_hx_pri
  #   connection_index = 0
  #   execute_on = 'INITIAL TIMESTEP_END'
  # []
  # [core_m_dot_difference]
  #   type = DifferencePostprocessor
  #   value1 = core_m_dot_out
  #   value2 = core_m_dot_in
  #   execute_on = 'INITIAL TIMESTEP_END'
  # []

  ###### power

  #####################################
  #######        RANK 14        #######
  #####################################

  ### column 1

  [core_r14_c1_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r14_c1
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 3

  [core_r14_c3_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r14_c3
    execute_on = 'INITIAL TIMESTEP_END'
  []

  #####################################
  #######        RANK 13        #######
  #####################################

  ### column 1

  [core_r13_c1_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r13_c1
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 2

  [core_r13_c2_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r13_c2
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 3

  [core_r13_c3_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r13_c3
    execute_on = 'INITIAL TIMESTEP_END'
  []

  #####################################
  #######        RANK 12        #######
  #####################################

  ### column 0

  [core_r12_c0_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r12_c0
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 1

  [core_r12_c1_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r12_c1
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 3

  [core_r12_c3_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r12_c3
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 4

  [core_r12_c4_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r12_c4
    execute_on = 'INITIAL TIMESTEP_END'
  []

  #####################################
  #######        RANK 11        #######
  #####################################

  ### column 1

  [core_r11_c1_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r11_c1
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 2

  [core_r11_c2_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r11_c2
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 3

  [core_r11_c3_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r11_c3
    execute_on = 'INITIAL TIMESTEP_END'
  []

  #####################################
  #######        RANK 10        #######
  #####################################

  ### column 0

  [core_r10_c0_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r10_c0
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 1

  [core_r10_c1_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r10_c1
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 3

  [core_r10_c3_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r10_c3
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 4

  [core_r10_c4_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r10_c4
    execute_on = 'INITIAL TIMESTEP_END'
  []

  #####################################
  #######        RANK 9         #######
  #####################################

  ### column 0

  [core_r9_c0_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r9_c0
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 1

  [core_r9_c1_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r9_c1
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 2

  [core_r9_c2_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r9_c2
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 3

  [core_r9_c3_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r9_c3
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 4

  [core_r9_c4_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r9_c4
    execute_on = 'INITIAL TIMESTEP_END'
  []

  #####################################
  #######        RANK 8         #######
  #####################################

  ### column 0

  [core_r8_c0_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r8_c0
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 1

  [core_r8_c1_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r8_c1
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 3

  [core_r8_c3_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r8_c3
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 4

  [core_r8_c4_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r8_c4
    execute_on = 'INITIAL TIMESTEP_END'
  []

  #####################################
  #######        RANK 7         #######
  #####################################

  ### column 0

  [core_r7_c0_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r7_c0
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 1

  [core_r7_c1_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r7_c1
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 2

  [core_r7_c2_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r7_c2
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 3

  [core_r7_c3_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r7_c3
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 4

  [core_r7_c4_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r7_c4
    execute_on = 'INITIAL TIMESTEP_END'
  []

  #####################################
  #######        RANK 6         #######
  #####################################

  ### column 0

  [core_r6_c0_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r6_c0
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 1

  [core_r6_c1_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r6_c1
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 3

  [core_r6_c3_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r6_c3
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 4

  [core_r6_c4_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r6_c4
    execute_on = 'INITIAL TIMESTEP_END'
  []

  #####################################
  #######        RANK 5         #######
  #####################################

  ### column 0

  [core_r5_c0_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r5_c0
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 1

  [core_r5_c1_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r5_c1
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 2

  [core_r5_c2_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r5_c2
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 3

  [core_r5_c3_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r5_c3
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 4

  [core_r5_c4_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r5_c4
    execute_on = 'INITIAL TIMESTEP_END'
  []

  #####################################
  #######        RANK 4         #######
  #####################################

  ### column 0

  [core_r4_c0_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r4_c0
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 1

  [core_r4_c1_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r4_c1
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 3

  [core_r4_c3_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r4_c3
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 4

  [core_r4_c4_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r4_c4
    execute_on = 'INITIAL TIMESTEP_END'
  []

  #####################################
  #######        RANK 3         #######
  #####################################

  ### column 1

  [core_r3_c1_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r3_c1
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 2

  [core_r3_c2_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r3_c2
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 3

  [core_r3_c3_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r3_c3
    execute_on = 'INITIAL TIMESTEP_END'
  []

  #####################################
  #######        RANK 2         #######
  #####################################

  ### column 0

  [core_r2_c0_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r2_c0
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 1

  [core_r2_c1_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r2_c1
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 3

  [core_r2_c3_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r2_c3
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 4

  [core_r2_c4_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r2_c4
    execute_on = 'INITIAL TIMESTEP_END'
  []

  #####################################
  #######        RANK 1         #######
  #####################################

  ### column 1

  [core_r1_c1_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r1_c1
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 2

  [core_r1_c2_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r1_c2
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 3

  [core_r1_c3_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r1_c3
    execute_on = 'INITIAL TIMESTEP_END'
  []

  #####################################
  #######        RANK 0         #######
  #####################################

  ### column 1

  [core_r0_c1_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r0_c1
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ### column 3

  [core_r0_c3_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse 2 * pi * core_nb_coolant_per_assembly * core_radius_coolant}'
    block = core/coolant_channel_r0_c3
    execute_on = 'INITIAL TIMESTEP_END'
  []

  #####################################
  #######       CORE TOT        #######
  #####################################

  [core_power_in]
    type = ADFlowJunctionFlux1Phase
    boundary = plenum_inlet:in
    equation = energy
    junction = pri_jct_upc_out_plm_in
    connection_index = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [core_power_out]
    type = ADFlowJunctionFlux1Phase
    boundary = plenum_outlet:out
    equation = energy
    junction = pri_jct_plenum_out_2
    connection_index = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [core_power_difference]
    type = DifferencePostprocessor
    value1 = core_power_out
    value2 = core_power_in
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [core_power_flux]
    type = SumPostprocessor
    execute_on = TIMESTEP_END
    values = 'core_r14_c1_power_flux core_r14_c3_power_flux
    core_r13_c1_power_flux core_r13_c2_power_flux  core_r13_c3_power_flux
    core_r12_c0_power_flux core_r12_c1_power_flux core_r12_c3_power_flux core_r12_c4_power_flux
    core_r11_c1_power_flux core_r11_c2_power_flux  core_r11_c3_power_flux
    core_r10_c0_power_flux core_r10_c1_power_flux core_r10_c3_power_flux core_r10_c4_power_flux
    core_r9_c0_power_flux core_r9_c1_power_flux core_r9_c2_power_flux  core_r9_c3_power_flux  core_r9_c4_power_flux
    core_r8_c0_power_flux core_r8_c1_power_flux core_r8_c3_power_flux core_r8_c4_power_flux
    core_r7_c0_power_flux  core_r7_c1_power_flux  core_r7_c2_power_flux   core_r7_c3_power_flux  core_r7_c4_power_flux
    core_r6_c0_power_flux core_r6_c1_power_flux core_r6_c3_power_flux core_r6_c4_power_flux
    core_r5_c0_power_flux core_r5_c1_power_flux core_r5_c2_power_flux  core_r5_c3_power_flux core_r5_c4_power_flux
    core_r4_c0_power_flux core_r4_c1_power_flux core_r4_c3_power_flux core_r4_c4_power_flux
    core_r3_c1_power_flux core_r3_c2_power_flux  core_r3_c3_power_flux
    core_r2_c0_power_flux core_r2_c1_power_flux core_r2_c3_power_flux core_r2_c4_power_flux
    core_r1_c1_power_flux core_r1_c2_power_flux  core_r1_c3_power_flux
    core_r0_c1_power_flux core_r0_c3_power_flux'
  []

  ###### Temperature

  [core_T_in]
    type = SideAverageValue
    boundary = pri_jct_upc_out_plm_in
    variable = T
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [core_T_out]
    type = SideAverageValue
    boundary = pri_jct_plenum_out_2
    variable = T
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [core_T_difference]
    type = DifferencePostprocessor
    value1 = core_T_out
    value2 = core_T_in
    execute_on = 'INITIAL TIMESTEP_END'
  []

  #####################################
  #######        UPCOMER        #######
  #####################################

  [upc_power_in]
    type = ADFlowJunctionFlux1Phase
    boundary = upcomer:in
    equation = energy
    junction = pri_jct_1_upcomer
    connection_index = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [upc_power_out]
    type = ADFlowJunctionFlux1Phase
    boundary = upcomer_out:out
    equation = energy
    junction = pri_jct_upc_out_plm_in
    connection_index = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [upc_power_difference]
    type = DifferencePostprocessor
    value1 = upc_power_out
    value2 = upc_power_in
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ###### T

  [upc_T_in]
    type = SideAverageValue
    boundary = pri_jct_1_upcomer
    variable = T
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [upc_T_out]
    type = SideAverageValue
    boundary = pri_jct_upc_out_plm_in
    variable = T
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [upc_T_difference]
    type = DifferencePostprocessor
    value1 = upc_T_out
    value2 = upc_T_in
    execute_on = 'INITIAL TIMESTEP_END'
  []

  #####################################################
  ######             heat_exchanger              ######
  #####################################################

  ################### primary loop ####################

  ###### pressure

  [hx_pri_p_in]
    type = SideAverageValue
    boundary = pri_jct_3_hx_pri
    variable = p
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [hx_pri_p_out]
    type = SideAverageValue
    boundary = pri_jct_hx_pri_4
    variable = p
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [hx_pri_p_difference]
    type = DifferencePostprocessor
    value1 = hx_pri_p_out
    value2 = hx_pri_p_in
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ###### m_dot

  [hx_pri_m_dot_in]
    type = ADFlowJunctionFlux1Phase
    boundary = hx/pri:in
    equation = mass
    junction = pri_jct_3_hx_pri
    connection_index = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [hx_pri_m_dot_out]
    type = ADFlowJunctionFlux1Phase
    boundary = hx/pri:out
    equation = mass
    junction = pri_jct_hx_pri_4
    connection_index = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [hx_pri_m_dot_difference]
    type = DifferencePostprocessor
    value1 = hx_pri_m_dot_out
    value2 = hx_pri_m_dot_in
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ###### temperature

  [hx_pri_T_in]
    type = SideAverageValue
    boundary = pri_jct_3_hx_pri
    variable = T
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [hx_pri_T_out]
    type = SideAverageValue
    boundary = pri_jct_hx_pri_4
    variable = T
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [hx_pri_T_difference]
    type = DifferencePostprocessor
    value1 = hx_pri_T_out
    value2 = hx_pri_T_in
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ###### energy

  [hx_pri_power_in]
    type = ADFlowJunctionFlux1Phase
    boundary = hx/pri:in
    equation = energy
    junction = pri_jct_3_hx_pri
    connection_index = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [hx_pri_power_out]
    type = ADFlowJunctionFlux1Phase
    boundary = hx/pri:out
    equation = energy
    junction = pri_jct_hx_pri_4
    connection_index = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [hx_pri_power_difference]
    type = DifferencePostprocessor
    value1 = hx_pri_power_out
    value2 = hx_pri_power_in
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [hx_pri_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse hx_nb_channels * pi * hx_dia_pri}'
    block = hx/pri
    execute_on = 'INITIAL TIMESTEP_END'
  []

  #####################################################
  ######                   pump                  ######
  #####################################################

  ###### pressure

  [pump_p_in]
    type = SideAverageValue
    boundary = pri_pipe4:out
    variable = p
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [pump_p_out]
    type = SideAverageValue
    boundary = pri_pipe5:in
    variable = p
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [pump_p_difference]
    type = DifferencePostprocessor
    value1 = pump_p_out
    value2 = pump_p_in
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ##### temperature

  [pump_T_in]
    type = SideAverageValue
    boundary = pri_pipe4:out
    variable = T
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [pump_T_out]
    type = SideAverageValue
    boundary = pri_pipe5:in
    variable = T
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [pump_T_difference]
    type = DifferencePostprocessor
    value1 = pump_T_out
    value2 = pump_T_in
    execute_on = 'INITIAL TIMESTEP_END'
  []

  ###### Core properties

  [Hw_core_max]
    type = ADElementExtremeMaterialProperty
    block = 'core/coolant_channel_r14_c1
    core/coolant_channel_r14_c3
    core/coolant_channel_r13_c1
    core/coolant_channel_r13_c2
    core/coolant_channel_r13_c3
    core/coolant_channel_r12_c0
    core/coolant_channel_r12_c1
    core/coolant_channel_r12_c3
    core/coolant_channel_r12_c4
    core/coolant_channel_r11_c1
    core/coolant_channel_r11_c2
    core/coolant_channel_r11_c3
    core/coolant_channel_r10_c0
    core/coolant_channel_r10_c1
    core/coolant_channel_r10_c3
    core/coolant_channel_r10_c4
    core/coolant_channel_r9_c0
    core/coolant_channel_r9_c1
    core/coolant_channel_r9_c2
    core/coolant_channel_r9_c3
    core/coolant_channel_r9_c4
    core/coolant_channel_r8_c0
    core/coolant_channel_r8_c1
    core/coolant_channel_r8_c3
    core/coolant_channel_r8_c4
    core/coolant_channel_r7_c0
    core/coolant_channel_r7_c1
    core/coolant_channel_r7_c2
    core/coolant_channel_r7_c3
    core/coolant_channel_r7_c4
    core/coolant_channel_r6_c0
    core/coolant_channel_r6_c1
    core/coolant_channel_r6_c3
    core/coolant_channel_r6_c4
    core/coolant_channel_r5_c0
    core/coolant_channel_r5_c1
    core/coolant_channel_r5_c2
    core/coolant_channel_r5_c3
    core/coolant_channel_r5_c4
    core/coolant_channel_r4_c0
    core/coolant_channel_r4_c1
    core/coolant_channel_r4_c3
    core/coolant_channel_r4_c4
    core/coolant_channel_r3_c1
    core/coolant_channel_r3_c2
    core/coolant_channel_r3_c3
    core/coolant_channel_r2_c0
    core/coolant_channel_r2_c1
    core/coolant_channel_r2_c3
    core/coolant_channel_r2_c4
    core/coolant_channel_r1_c1
    core/coolant_channel_r1_c2
    core/coolant_channel_r1_c3
    core/coolant_channel_r0_c1
    core/coolant_channel_r0_c3'
    mat_prop = Hw
    value_type = max
  []
  [Hw_core_min]
    type = ADElementExtremeMaterialProperty
    block = 'core/coolant_channel_r14_c1
    core/coolant_channel_r14_c3
    core/coolant_channel_r13_c1
    core/coolant_channel_r13_c2
    core/coolant_channel_r13_c3
    core/coolant_channel_r12_c0
    core/coolant_channel_r12_c1
    core/coolant_channel_r12_c3
    core/coolant_channel_r12_c4
    core/coolant_channel_r11_c1
    core/coolant_channel_r11_c2
    core/coolant_channel_r11_c3
    core/coolant_channel_r10_c0
    core/coolant_channel_r10_c1
    core/coolant_channel_r10_c3
    core/coolant_channel_r10_c4
    core/coolant_channel_r9_c0
    core/coolant_channel_r9_c1
    core/coolant_channel_r9_c2
    core/coolant_channel_r9_c3
    core/coolant_channel_r9_c4
    core/coolant_channel_r8_c0
    core/coolant_channel_r8_c1
    core/coolant_channel_r8_c3
    core/coolant_channel_r8_c4
    core/coolant_channel_r7_c0
    core/coolant_channel_r7_c1
    core/coolant_channel_r7_c2
    core/coolant_channel_r7_c3
    core/coolant_channel_r7_c4
    core/coolant_channel_r6_c0
    core/coolant_channel_r6_c1
    core/coolant_channel_r6_c3
    core/coolant_channel_r6_c4
    core/coolant_channel_r5_c0
    core/coolant_channel_r5_c1
    core/coolant_channel_r5_c2
    core/coolant_channel_r5_c3
    core/coolant_channel_r5_c4
    core/coolant_channel_r4_c0
    core/coolant_channel_r4_c1
    core/coolant_channel_r4_c3
    core/coolant_channel_r4_c4
    core/coolant_channel_r3_c1
    core/coolant_channel_r3_c2
    core/coolant_channel_r3_c3
    core/coolant_channel_r2_c0
    core/coolant_channel_r2_c1
    core/coolant_channel_r2_c3
    core/coolant_channel_r2_c4
    core/coolant_channel_r1_c1
    core/coolant_channel_r1_c2
    core/coolant_channel_r1_c3
    core/coolant_channel_r0_c1
    core/coolant_channel_r0_c3'
    mat_prop = Hw
    value_type = min
  []

  ################################################################################################################
  ################  ********************    SECONDARY LOOP POSTPROCESSORS   *******************  #################
  ################################################################################################################

  ##########################
  # Heat exchanger
  ##########################

  ######### temperature

  [hx_sec_T_in]
    type = SideAverageValue
    boundary = hx_sec_inlet
    variable = T
  []
  [hx_sec_T_out]
    type = SideAverageValue
    boundary = hx_sec_outlet
    variable = T
  []
  [hx_sec_T_difference]
    type = DifferencePostprocessor
    value1 = hx_sec_T_out
    value2 = hx_sec_T_in
  []

  ######### pressure

  [hx_sec_p_in]
    type = SideAverageValue
    boundary = hx_sec_inlet
    variable = p
  []
  [hx_sec_p_out]
    type = SideAverageValue
    boundary = hx_sec_outlet
    variable = p
  []
  [hx_sec_p_difference]
    type = DifferencePostprocessor
    value1 = hx_sec_p_out
    value2 = hx_sec_p_in
  []

  ######### flow rate

  [hx_sec_m_dot_in]
    type = ADFlowBoundaryFlux1Phase
    boundary = hx_sec_inlet
    equation = mass
  []
  [hx_sec_m_dot_out]
    type = ADFlowBoundaryFlux1Phase
    boundary = hx_sec_outlet
    equation = mass
  []
  [hx_sec_m_dot_difference]
    type = DifferencePostprocessor
    value1 = hx_sec_m_dot_out
    value2 = hx_sec_m_dot_in
  []

  ######### power

  [hx_sec_power_in]
    type = ADFlowBoundaryFlux1Phase
    boundary = hx_sec_inlet
    equation = energy
  []
  [hx_sec_power_out]
    type = ADFlowBoundaryFlux1Phase
    boundary = hx_sec_outlet
    equation = energy
  []
  [hx_sec_power_difference]
    type = DifferencePostprocessor
    value1 = hx_sec_power_out
    value2 = hx_sec_power_in
  []
  [hx_sec_power_flux]
    type = ADHeatRateConvection1Phase
    P_hf = '${fparse pi * hx_dia_sec * hx_nb_channels }'
    block = hx/sec
  []
[]

##########################
##########################
####                  ####
####        8         ####
####                  ####
####    EXECUTION     ####
####                  ####
####                  ####
##########################
##########################

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'

  end_time = 1000000
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    growth_factor = 1.2
    cutback_factor = 0.8
  []
  dtmin = 1e-7
  dtmax = 100000
  steady_state_detection = true
  steady_state_start_time = 2000
  solve_type = NEWTON
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  nl_max_its = 15

  l_tol = 1e-4
  l_max_its = 10

  # automatic_scaling = true

  petsc_options_iname = '-pc_type'
  petsc_options_value = ' lu     '
[]

[Outputs]
  [e]
    type = Exodus
    file_base = 'htgr_thm_primary'
  []

  [csv]
    type = CSV
    file_base = 'htgr_thm_primary'
  []

  checkpoint = true

  [console]
    type = Console
    max_rows = 10
    show = 'core_T_in core_T_out  core_p_in core_p_out core_p_difference
    hx_pri_T_in hx_pri_T_out hx_pri_p_in hx_pri_p_out hx_pri_p_difference hx_pri_m_dot_in hx_pri_m_dot_difference hx_pri_power_in hx_pri_power_out hx_pri_power_difference


    Hw_core_max Hw_core_min'

    # pump_T_in pump_T_out pump_p_in pump_p_out pump_p_difference
  []
[]
