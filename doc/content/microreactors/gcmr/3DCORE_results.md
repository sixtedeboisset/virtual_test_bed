Some differences appear between the results of the complete core coupled with the 1D model of the primary loop and those of the 1D model of the loops. The channels going through the core can experience different conditions (particularly the power provided to each one). Thus, the mass flow rate and consequently the fluid temperature and heat transfer coefficient can differ.

The mass flow rate decreases from 14 kg/s to its steady state value (which is 8.5 kg/s). As shown in [FIG–PLOT-m_dot_primary], it takes approximately 400 s to reach the steady state. It is comparable in some extent with the 1D simulation. The steady state is reach after a longer time, because in that case the complete secondary loop is also considered. The coupling changes also a bit the operating values, which are lower than with the 1D simulation using a simplified core.

[FIG–PLOT-m_dot_primary]

As imagined, the primary pressures shown in [FIG–PLOT-p_primary] stay close to 90 bar. This is imposed as an initial condition in the whole primary loop and during the operations by the pressurizer between the core and the heat exchanger. The secondary pressure is imposed at the heat exchanger outlet at 5.5 bar, and the pressure drop in this component is about 0.1 bar.

[FIG–PLOT-p_primary]

The fluid temperatures are a bit higher at the core and heat exchanger (secondary side) outlets. They are plotted in [FIG–PLOT–T_fluid]. It is due to the lower mass flow rate: the provided thermal power is the same, consequently the temperature difference must be higher. The core inlet temperature stays nevertheless close to 900 K and the higher temperatures are not too far from 1200 K, which is the 1D simulation core outlet temperature. The 3D representation, given in [FIG–3D–T_fluid], shows that the fluid takes power from the core when it goes through it.

[FIG–PLOT–T_fluid]
[FIG–3D–T_fluid]

The solid temperatures are consequently a bit higher for the maximum and average values than in the 1D simulation, but they stay close to what was got in that case. They are plotted in [FIG–PLOT–T_solid]. A 3D representation is also provided in [FIG–3D–T_solid] which makes appear clearly that the solid temperature increases progressively in the core.

[FIG–PLOT–T_solid]
[FIG–3D–T_solid]

The core heat transfer coefficient is also a bit lower than what was got in the 1D simulation. Its evolution is plotted in [FIG–PLOT–Hw_core] and a 3D representation of the steady state is provided in [FIG–3D–Hw_core] (in this case, the colors are put on the solid surrounding the boundaries where the heat transfer coefficient is defined). This effect is due to the lowest mass flow rate, because the Dittus-Boelter correlation, which is used here, gives:

[h = 0.023 lambda/Dh * Re^(4/5) * Pr^(0.4)], with [Re = (rho u Dh)/mu = 4/pi * m_dot/(Dh * mu)]
Consequently, the lower mass flow rate got in this simulation induces a lower heat transfer coefficient. It is possible to expect this type of results, because the heat transfer is logically easier when the mass flow rate is bigger.

[FIG–PLOT–Hw_core]
[FIG–3D–Hw_core]

The last step is to check the power transfer in the system, plotted in [FIG–PLOT–power]. The power transferred through the coolant channels boundaries is quickly the same than the power difference in the fluid between the upcomer inlet and core outlet. A longer transient occurs in the secondary side of the heat exchanger, but the transferred steady state power is the same.

[FIG–PLOT–power]

