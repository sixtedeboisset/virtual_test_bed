# Description of the system


To complete the 1D thermal-hydraulic simulation of the HTGR, a 3D model of the core is built. The geometry is presented in [!citep](Duchnowski2022): the core contains 55 fuel assemblies. Each of them is composed of a graphite matrix, with a hexagonal section of 0.11 m radius and is 2 meters high. It contains 18 coolant channels, whose radius is 6.35 mm, and 42 fuel channels, whose radius is 7.94 mm. The lattice pitch between the channels is 0.022 m. The following pitches are used for the complete core and for each assembly and presented in [core_pattern] and [assembly_pattern]:

!media media/gcmr/3D_core_balance_of_plant/description_core_pattern.png
      style=display: block;margin-left:auto;margin-right:auto;width:40%;
      id=core_pattern
      caption= Core pattern

!media media/gcmr/3D_core_balance_of_plant/description_assembly_pattern.png
      style=display: block;margin-left:auto;margin-right:auto;width:40%;
      id=assembly_pattern
      caption= Assembly pattern

The coolant is Helium, the fuel channel contains TRISO particles in a graphite matrix. It is H-451 graphite as specified in [!citep](Duchnowski2022), and SA508 is used as the steel type for the vessel. It is one of the most considered alloys considered to make the reactor pressure vessel of HTGRs.

The assemblies are surrounded by an outer shield. The core is in a steel vessel. Between this vessel and the outer shield, a gap of 5 cm surrounds the core and brings the primary fluid to the core inlet plenum. In this gap and before it enters in the coolant channels of the core, the fluid gets some power. A fraction of the thermal power is also lost by radiation of the vessel, heated by the primary fluid flow in the gap and by radiative exchange from the graphite outer shield.

When helium escapes from the core, it goes in the core outlet plenum and then in the pipes of the primary loop.

These elements are summarized in [core_diagram]:

!media media/gcmr/3D_core_balance_of_plant/core_with_plena.png
      style=display: block;margin-left:auto;margin-right:auto;width:40%;
      id=core_diagram
      caption= Core diagram

The primary thermal hydraulic system is the same than presented for the 1D model of the balance of plant, except the fact that a more precise model of the flow in the core is considered (as presented above). The secondary side is reduced to the heat exchanger with the operating parameters measured in the 1D simulation.
