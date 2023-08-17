# parameters of the coolant channels

radius_coolant = 0.00635 # m
length_channel = 2 # m
channel_n_elems = 50

# numbers of channels and assemblies

nb_assembly = 55
nb_coolant_per_assembly = 18
nb_fuel_per_assembly = 42
nb_coolant_tot = '${fparse nb_assembly * nb_coolant_per_assembly}'

# other parameters of the assembly

lattice_pitch = 0.022 # m
pincell_apothem = '${fparse lattice_pitch / 2}'

# Apothem of PolygonConcentricCircleMeshGenerator:
#
#                 #
#             #       #
#         #               #
#         #<------        #
#         #               #
#             #       #
#                 #

assembly_radius = '${fparse 5 * lattice_pitch}'
assembly_apothem = '${fparse sqrt(3) / 2 * assembly_radius}'
radius_fuel = 0.00794 # m

# calculus of the equivalent parameters of a cylindrical heat structure around the coolant channel

section_assembly = '${fparse 3 * ( sqrt(3) / 2) * 5 * lattice_pitch * 5 * lattice_pitch}' # calculus of the assembly section (hexagonal) using the lattice pitch (1/5 of a side of the hexagon)
section_fuel_channel = '${fparse pi * radius_fuel * radius_fuel}'
section_moderator_and_coolant = '${fparse section_assembly - ( nb_fuel_per_assembly * section_fuel_channel )}'
# radius_equivalent_graphite = '${fparse sqrt(section_moderator_and_coolant / pi)}'

# radius_equivalent_channel = '${fparse sqrt(nb_coolant_per_assembly) * radius_coolant}'

L_plenum_inlet = 0.2
L_upcomer_out = 0.5

x_core_origin = ${fparse L_plenum_inlet + L_upcomer_out}

[Mesh]

  #################################      # This parameter allows us to execute the file but stop at this block so we can see intermediate output.
  # final_generator = extrude # User: Change this based on build step
  ################################

  ### Step 1. Create Pin Unit Cells
  # There are 3 unique pin in the fuel assembly

  [moderator_pincell]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6 # must be six to use hex pattern
    num_sectors_per_side = '2 2 2 2 2 2 '
    background_intervals = 2.
    background_block_ids = '1'
    polygon_size = ${pincell_apothem}
    polygon_size_style = 'apothem'
    ring_radii = '${fparse pincell_apothem / 2 }'
    ring_intervals = '2'
    ring_block_ids = '103 101' # 103 is tri mesh
    preserve_volumes = on
    quad_center_elements = false
  []
  [coolant_pincell]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6 # must be six to use hex pattern
    num_sectors_per_side = '2 2 2 2 2 2 '
    background_intervals = 2.
    background_block_ids = '2'
    polygon_size = ${pincell_apothem}
    polygon_size_style = 'apothem'
    ring_radii = '${radius_coolant}'
    ring_intervals = '2'
    ring_block_ids = '203 201' # 203 is tri mesh
    preserve_volumes = on
    quad_center_elements = false
  []
  [fuel_pincell_ring1]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6 # must be six to use hex pattern
    num_sectors_per_side = '2 2 2 2 2 2 '
    background_intervals = 2.
    background_block_ids = '3'
    polygon_size = ${pincell_apothem}
    polygon_size_style = 'apothem'
    ring_radii = '${radius_fuel}'
    ring_intervals = '5'
    ring_block_ids = '303 301' # 303 is tri mesh
    preserve_volumes = on
    quad_center_elements = false
  []

  [fuel_pincell_ring2]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6 # must be six to use hex pattern
    num_sectors_per_side = '2 2 2 2 2 2 '
    background_intervals = 2.
    background_block_ids = '4'
    polygon_size = ${pincell_apothem}
    polygon_size_style = 'apothem'
    ring_radii = '${radius_fuel}'
    ring_intervals = '5'
    ring_block_ids = '403 401' # 403 is tri mesh
    preserve_volumes = on
    quad_center_elements = false
  []

  [fuel_pincell_ring3]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6 # must be six to use hex pattern
    num_sectors_per_side = '2 2 2 2 2 2 '
    background_intervals = 2.
    background_block_ids = '5'
    polygon_size = ${pincell_apothem}
    polygon_size_style = 'apothem'
    ring_radii = '${radius_fuel}'
    ring_intervals = '5'
    ring_block_ids = '503 501' # 503 is tri mesh
    preserve_volumes = on
    quad_center_elements = false
  []

  [fuel_pincell_ring4]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6 # must be six to use hex pattern
    num_sectors_per_side = '2 2 2 2 2 2 '
    background_intervals = 2.
    background_block_ids = '6'
    polygon_size = ${pincell_apothem}
    polygon_size_style = 'apothem'
    ring_radii = '${radius_fuel}'
    ring_intervals = '5'
    ring_block_ids = '603 601' # 603 is tri mesh
    preserve_volumes = on
    quad_center_elements = false
  []

  ### Step 2. Create Fuel assembly

  [fuel_assembly]
    type = PatternedHexMeshGenerator
    inputs = 'coolant_pincell fuel_pincell_ring1 fuel_pincell_ring2 fuel_pincell_ring3 fuel_pincell_ring4 moderator_pincell'
    # Pattern ID     0                 1                  2                  3                  4                   5
    hexagon_size = ${assembly_apothem}
    background_block_id = 10
    background_intervals = 1
    pattern = '4 4 0 4 4;
              4 0 3 3 0 4;
             0 3 2 0 2 3 0;
            4 3 0 1 1 0 3 4;
           4 0 2 1 5 1 2 0 4;
            4 3 0 1 1 0 3 4;
             0 3 2 0 2 3 0;
              4 0 3 3 0 4;
               4 4 0 4 4;'
  []

  ### Step 3. Create dummy assembly and use it with fuel assembly to create the pattern of the core

  [dummy_assembly]
    type = HexagonConcentricCircleAdaptiveBoundaryMeshGenerator
    num_sectors_per_side = '4 4 4 4 4 4'
    hexagon_size = ${assembly_apothem}
    background_intervals = 2
    background_block_ids = '20 21'
    # external_boundary_id = 9998
  []

  [core]
    type = PatternedHexMeshGenerator
    inputs = 'fuel_assembly dummy_assembly'
    # Pattern ID     0        1
    pattern_boundary = none
    generate_core_metadata = true
    pattern = '1 0 0 0 1;
              0 0 0 0 0 0;
             0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0;
           1 0 0 0 0 0 0 0 1;
            0 0 0 0 0 0 0 0;
             0 0 0 0 0 0 0;
              0 0 0 0 0 0;
               1 0 0 0 1'
  []

  ### Step 4. Delete Dummy Assemblies

  [core_dummy_del]
    type = BlockDeletionGenerator
    block = '20 21'
    input = core
  []

  ### Step 5. Add Core Periphery

  [outer_shield]
    type = PeripheralRingMeshGenerator
    input = core_dummy_del
    peripheral_layer_num = 1
    peripheral_ring_radius = 0.9
    input_mesh_external_boundary = 10000
    peripheral_ring_block_id = 250
    peripheral_ring_block_name = outer_shield
  []

  [gap_with_vessel]
    type = PeripheralRingMeshGenerator
    input = outer_shield
    peripheral_layer_num = 1
    peripheral_ring_radius = 0.95
    input_mesh_external_boundary = 10000
    peripheral_ring_block_id = 350
    peripheral_ring_block_name = gap_with_vessel
  []

  [vessel]
    type = PeripheralRingMeshGenerator
    input = gap_with_vessel
    peripheral_layer_num = 1
    peripheral_ring_radius = 1.15
    input_mesh_external_boundary = 10000
    peripheral_ring_block_id = 450
    peripheral_ring_block_name = vessel
  []

  [vessel_bis]
    type = PeripheralRingMeshGenerator
    input = vessel
    peripheral_layer_num = 1
    peripheral_ring_radius = 1.20
    input_mesh_external_boundary = 10000
    peripheral_ring_block_id = 550
    peripheral_ring_block_name = vessel_bis
  []

  # ### Step 8. Slice to 1/6 Core

  # [coreslice_1]
  #   type = PlaneDeletionGenerator
  #   point = '0 0 0'
  #   normal = '10 17.32 0'
  #   input = outer_shield
  #   new_boundary = 147
  # []

  # [coreslice_2]
  #   type = PlaneDeletionGenerator
  #   point = '0 0 0'
  #   normal = '10 -17.32 0'
  #   input = coreslice_1
  #   new_boundary = 147
  # []

  ### Step 6. Extrude to 3D

  [extrude]
    type = AdvancedExtruderGenerator
    input = vessel_bis
    # heights = '0.2 1.6 0.2'
    # num_layers = '1 8 1'

    heights = '2'
    num_layers = '10'

    direction = '0 0 1'
    # subdomain_swaps = '10 1000 100 1000 101 1000 103 1003 200 1000 201 1000 203 1003 301 1000 303 1003;
    #                     10 10   100 100  101 101  103 103  200 200  201 201  203 203  301 301  303 303;
    #                     10 1000 100 1000 101 1000 103 1003 200 200  201 201  203 203  301 1000 303 1003'
    subdomain_swaps = '10 10   100 100  101 101  103 103  200 200  201 201  203 203  301 301  303 303  401 401  403 403  501 501  503 503  601 601  603 603'

    # 10 = graphite matrix
    # 100 = "graphite channel"
    # 200 = coolant channel
    # 300, 400, 500, 600 = fuel channel

    top_boundary = 2000
    bottom_boundary = 3000
  []

  ### Step 7. Create the boundary with the flow channel

  [add_coolant_boundary]
    type = SideSetsBetweenSubdomainsGenerator
    input = extrude
    primary_block = '2'
    paired_block = '201'
    new_boundary = coolant_boundary
  []

  ### Step 8. Delete coolant channels

  [delete_coolant_channels]
    type = BlockDeletionGenerator
    block = '201 203'
    input = add_coolant_boundary
  []

  ### Step 9. Add a boundary around the reflector

  [add_reflector_boundary]
    type = SideSetsBetweenSubdomainsGenerator
    input = delete_coolant_channels
    primary_block = '250'
    paired_block = '350'
    new_boundary = reflector_boundary
  []

  ### Step 10. Add the boundaries (in and out) of the vessel

  [add_vessel_boundary_in]
    type = SideSetsBetweenSubdomainsGenerator
    input = add_reflector_boundary
    primary_block = '450'
    paired_block = '350'
    new_boundary = vessel_boundary_in
  []

  [add_vessel_boundary_out]
    type = SideSetsBetweenSubdomainsGenerator
    input = add_vessel_boundary_in
    primary_block = '450'
    paired_block = '550'
    new_boundary = vessel_boundary_out
  []

  ### Step 11. Delete dummy blocks

  [delete_gap]
    type = BlockDeletionGenerator
    block = '350'
    input = add_vessel_boundary_out
  []

  [delete_vessel_bis]
    type = BlockDeletionGenerator
    block = '550'
    input = delete_gap
  []


  ### Step 12. Rotate the core to get the same orientation than the primary loop

  [core_rotate]
    type = TransformGenerator
    input = delete_vessel_bis
    transform = ROTATE
    vector_value = '90. 90. 90.'
  []

  [core_translate]
    type = TransformGenerator
    input = core_rotate
    transform = TRANSLATE
    vector_value = '${x_core_origin} 0. 0.'
  []

[]
