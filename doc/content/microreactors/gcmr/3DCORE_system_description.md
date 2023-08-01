# Description of the system


To complete the 1D thermal-hydraulic simulation of the HTGR, a 3D model of the core is built. The geometry is presented in [!citep](Duchnowski2022): the core contains 55 fuel assemblies. Each of them is composed of a graphite matrix, with a hexagonal section of 0.11 m radius and is 2 meters high. It contains 18 coolant channels, whose radius is 6.35 mm, and 42 fuel channels, whose radius is 7.94 mm. The lattice pitch between the channels is 0.022 m. The following pitches are used for the complete core and for each assembly and presented in [core_pattern] and [assembly_pattern]:

!media media/gcmr/balance_of_plant/description_core_pattern.png
      style=display: block;margin-left:auto;margin-right:auto;width:40%;
      id=core_pattern
      caption= Core pattern

!media media/gcmr/balance_of_plant/description_assembly_pattern.png
      style=display: block;margin-left:auto;margin-right:auto;width:40%;
      id=assembly_pattern
      caption= Assembly pattern

The coolant is Helium, the fuel channels contain TRISO particles in a graphite matrix.

The system obeys to the following heat conduction equation:

\begin{equation}
    \rho c \frac{\partial{T}}{\partial{t}} = div( k * \vec{grad}(T)) + q
\end{equation}

where $c$ is the specific heat, $\rho$ the volumic mass, $k$ the thermal conductivity, $T$ the temperature and $q$ the power sources.

The thermal hydraulic system is the same than presented in the Balance of plant section.
