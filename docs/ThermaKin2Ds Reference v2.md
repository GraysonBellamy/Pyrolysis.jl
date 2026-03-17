

# **Technical Reference for Numerical Pyrolysis Model ThermaKin2Ds**

**May 2019 Edition**

Stanislav I. Stoliarov

University of Maryland, Department of Fire Protection Engineering, College Park, MD 20742, USA

## **Summary**

This document provides a description of the governing equations, solution methodology and applications of a numerical pyrolysis solver, ThermaKin2Ds. ThermaKin2Ds is a result of further development of ThermaKin2D solver [1], which, in turn, is an extension of ThermaKin [2]. The main purpose of ThermaKin2Ds is to enable a detailed analysis and quantitative prediction of the processes that take place inside and at the surface of a burning material in a wide range of fire scenarios. The solver accounts for chemical reactions (chemical mechanisms of essentially arbitrary complexity can be modeled), transient conduction (governed by temperature dependent properties), in-depth thermal radiation absorption and emission, mass transport of gaseous species through the condensed phase and expansion or contraction of the pyrolyzing object.

ThermaKin2Ds can be used to model both anaerobic and oxidative pyrolysis decoupled from the gas-phase (flaming) combustion. ThermaKin2Ds can also be used to model interaction between pyrolyzing material and gas-phase flame, including various modes of flame spread. The gas-phase flame is represented in the ThermaKin2Ds through an empirical submodel that relates the rate of production of combustible gases by a material object to the spatially resolved heat feedback that this object receives from the flame fueled by these gases. In addition to the solver description, this document contains examples of ThermaKin2Ds input and output used to model thermal analysis and standard material flammability tests. The section of the document containing these examples is continuously updated to include new results.

## 1. Governing Equations

In ThermaKin2Ds, the material is represented by a mixture of components. By default, the maximum number of components is restricted to 50 (this restriction is used to optimize computer memory management and can be lifted by adjusting corresponding settings in the code). Every component is characterized by density, heat capacity, thermal conductivity, gas transfer coefficient, emissivity, and radiation absorption (or extinction) coefficient. The first 4 properties in this list are defined by a flexible function of temperature ( $T$ ),

$$\text{property} = p_0 + p_1 T + p_2 T^n \quad (1)$$

where  $p_0$ ,  $p_1$ ,  $p_2$ , and  $n$  are user-specified parameters. Emissivity and absorption coefficient are defined by single (constant) values. All components are divided into 3 categories: solids, liquids, and gases. This categorization is used in the calculation of material density as explained below.

The heat capacity or specific heat of material ( $c$ ) is calculated as

$$c = \sum_{i=1}^N mf_i c_i \quad (2)$$

where  $mf_i$  and  $c_i$  are mass fraction and heat capacity of the  $i$ -th component; and  $N$  is the number of components. The density of material ( $\rho$ ) is defined by

$$\rho = \frac{1}{\sum_{s=1}^{N_s} \frac{mf_s}{\rho_s} + \sum_{l=1}^{N_l} \frac{mf_l}{\rho_l} + \gamma \sum_{g=1}^{N_g} \frac{mf_g}{\rho_g}} \quad (3)$$

where  $\rho$  with a subscript designates component density. Subscripts  $s$ ,  $l$ , and  $g$  are used to refer to solid, liquid and gaseous components, respectively. Swelling factor  $\gamma$ , which may assume a value between 0 and 1, describes reaction of a volume of material to the presence of gases. When  $\gamma=0$ , the presence of gases has no effect on the volume. When  $\gamma=1$ , gases contribute to the volume of material in accordance with their respective densities.  $\gamma$  is calculated by volume-weighted averaging of the swelling factor specified for solids ( $\gamma_s$ ) and liquids ( $\gamma_l$ ):

$$\gamma = \frac{\gamma_s \sum_{s=1}^{N_s} \frac{mf_s}{\rho_s} + \gamma_l \sum_{l=1}^{N_l} \frac{mf_l}{\rho_l} + \tau \sum_{g=1}^{N_g} \frac{mf_g}{\rho_g}}{\sum_{s=1}^{N_s} \frac{mf_s}{\rho_s} + \sum_{l=1}^{N_l} \frac{mf_l}{\rho_l} + \tau \sum_{g=1}^{N_g} \frac{mf_g}{\rho_g}} \quad (4)$$

$\tau$  is a parameter, which is typically very small ( $\tau \ll 1$ ), that is used to ensure that, at the limit of very high gas content, the volume of material converges to that defined by gas densities.

### 1.1 Chemical Reactions

Components may undergo chemical reactions, which can also be used to represent physical transitions (such as melting). Up to 50 reactions can be specified (this restriction is used to optimize computer memory management and can be lifted by adjusting corresponding settings in the code). Each reaction may have 1 to 2 reactants and 0 to 2 products:

![](3c11688359c620619db37c4a1c01a46d_img.jpg)

$$\theta_1 \text{ COMP1} + \theta_2 \text{ COMP2} \rightarrow \theta_3 \text{ COMP3} + \theta_4 \text{ COMP4} + h \quad (5)$$

$\theta_i$  are stoichiometric coefficients and  $h$  is the heat of reaction. Different reactions may involve the same components (i.e., reactions can be coupled in parallel or consecutive fashion). The rate of reaction ( $r$ ) taking place in a unit volume of material is defined by

$$r = A \exp\left(-\frac{E}{RT}\right) \xi_{\text{COMP1}} \xi_{\text{COMP2}} \quad (6)$$

where  $A$  and  $E$  are the Arrhenius pre-exponential factor and activation energy; and  $R$  is the molar gas constant.  $\xi$  is the concentration of a given component expressed in the units of mass per unit volume ( $\xi_{\text{COMP1}} = mf_{\text{COMP1}} \rho$ ). In the absence of the second reactant,  $\xi_{\text{COMP2}}$  is set to 1.

The rate of consumption/formation of a reactant/product is calculated by multiplying  $r$  by the corresponding stoichiometric coefficient. The rate of production of heat is calculated by multiplying  $r$  by  $h$ . Following a convention frequently adopted in pyrolysis studies, a positive  $h$  is used to represent an exothermic reaction.  $h$  is defined by the same type of temperature dependence as that used for component properties (Equation 1). It should be noted that the reaction parameter definition is not constrained by either mass or energy conservation. It is a responsibility of the user to define reaction parameters in a way that ensures that masses of reactants and products are equal and that the heat capacities of reactants and products are related to the temperature dependence of  $h$  [3].

The reaction description also includes specification of a lower or upper temperature limit. If temperature decreases below the lower limit or increases above the upper limit, the rate of reaction is set to 0. Application of this limit increases computational efficiency (reaction rates are evaluated only within the temperature range where they are important) and facilitates usage of reactions for the description of thermodynamically-controlled phase transitions.

### 1.2 Heat and Mass Transfer

The conduction of heat is described by Fourier's law:

$$q_x = -k \frac{\partial T}{\partial x} \quad (7)$$

where  $q_x$  is the heat flux in the direction of Cartesian coordinate  $x$ ; and  $k$  is the thermal conductivity of material. The value of conductivity depends on relative amounts and spatial distribution of components [4]. If components are stacked in uniform layers that are normal to the direction of heat flow, the thermal conductivity is

$$k_n = \frac{1}{\sum_{i=1}^N \frac{v_i}{k_i}} \quad (8)$$

where  $k_i$  and  $v_i$  are the thermal conductivity and fractional volume of the  $i$ -th component (note that, when  $\gamma=0$ , gaseous components do not contribute to the thermal conductivity). On the other hand, if the layers are parallel to the direction of heat flow, the thermal conductivity is

$$k_p = \sum_{i=1}^N k_i v_i \quad (9)$$

For an arbitrary spatial distribution of components, the exact analytical expression of thermal conductivity is not available. However, under the assumption that components do not affect thermal conductivities of each other, Equations 8 and 9 provide lower and upper limits for the value of  $k$ . This means that the thermal conductivity of a multicomponent material can be represented as

$$k = \beta k_p + (1 - \beta) k_n \quad (10)$$

where  $\beta$  is a parameter that may assume a value between 0 and 1. In ThermaKin2Ds, this representation is used in conjunction with the assumption that a pyrolyzing material can be characterized by a single value of  $\beta$ .

Thermal radiation from an external source is absorbed inside the material according to a generalized version of the Beer-Lambert law:

$$\frac{\partial I_{ex}}{\partial x} = -I_{ex} \sum_{i=1}^N \alpha_i \xi_i \quad (11)$$

where  $I_{ex}$  is the flux of the radiation in  $x$  direction; and  $\alpha_i$  is the absorption (or extinction) coefficient of the  $i$ -th component. The radiation is not spectrally resolved. To comply with the second law of thermodynamics, the material is prescribed to re-radiate energy to the environment according to

$$\frac{\partial I_{rr}}{\partial x} = \frac{\sigma T^4 \sum_{i=1}^N \varepsilon_i v_i}{I_{ex}^0} \frac{\partial I_{ex}}{\partial x} \quad (12)$$

where  $I_{rr}$  is the heat flux radiated toward and through a material object boundary normal to  $x$ ;  $I_{ex}^0$  is the external radiation flux through the boundary (incident flux minus reflected); and  $\sigma$  is the Stefan-Boltzmann constant.  $\varepsilon_i$  is the emissivity of the  $i$ -th component, which in the current framework equals to the complement of reflectivity of the same component. In the case where no external radiation is applied,  $I_{ex}^0$  value used in Equation 12 is set to unity to produce a meaningful calculation of radiative loss.

Equations 11 and 12 describe radiative exchange between a material object and environment. The radiative transfer inside the object is modeled using the conduction equation (Equation 7) combined with the thermal conductivity expressed as the third power of temperature (using Equation 1). This approach is referred to as radiation diffusion approximation; it has been shown to be accurate for optically thick medium [5].

The transfer of mass is assumed to be driven by the gradient of volumetric fraction expressed through concentration:

$$J_g^x = -\rho_g \lambda \frac{\partial \left( \frac{\xi_g}{\rho_g} \right)}{\partial x} \quad (13)$$

Here,  $J_g^x$  is the mass flux of gas  $g$  in  $x$  direction. Only gaseous components are assumed to undergo this transfer.  $\lambda$  is the gas transfer coefficient of material. It is calculated from the corresponding component coefficients using the same approach as that utilized for the thermal conductivity (see Equations 8-10). Note that  $\lambda$  does not depend on the nature of gas that is being transferred (i.e., on volumetric basis, all gases subjected to the same volumetric fraction gradient are transferred at the same rate).

Application of Boyle's law, which states that the product of pressure and volume of a fixed amount of gas is constant, transforms Equation 13 into

$$J_g^x = -\frac{\rho_g \lambda}{P^{def}} \frac{\partial (\phi P_g)}{\partial x} \quad (14)$$

where  $\phi$  is the volume fraction of material occupied by gases;  $P_g$  is the partial pressure of gas  $g$ ; and  $P^{def}$  is the pressure at which the gas density ( $\rho_g$ ) is specified. This transformation helps clarify the transport model represented by Equation 13. If material expands proportionally to the volume of added gases (i.e.,  $\gamma > 0$  and pressure inside the material is constant), the flow of gas is a diffusion-like process. On the other hand, if material is rigid and does not expand with addition of gases (i.e.,  $\gamma = 0$  and  $\phi$  is constant), the transport model assumes the form of Darcy's law, which is frequently used to describe the flow of fluids through porous medium [6].

### 1.3 Conservation Equations

The chemical reaction and transport descriptions are connected through mass and energy conservation statements. These statements are formulated in two dimensions, either in Cartesian coordinates,  $x$  and  $y$ , or cylindrical coordinates,  $x$  and  $z$ . In the latter case,  $x$  is the axis of an axisymmetric object and  $z$  is the radial coordinate. The conservation of mass of component  $j$  is given by

$$\frac{\partial \xi_j}{\partial t} = \sum_{i=1}^{N_r} \theta_i^j r_i - \frac{\partial J_j^x}{\partial x} - \frac{\partial J_j^y}{\partial y} + \frac{\partial}{\partial x} \left( \xi_j \int_0^x \frac{1}{\rho} \frac{\partial \rho}{\partial t} dx \right) \quad (15)$$

and

$$\frac{\partial \xi_j}{\partial t} = \sum_{i=1}^{N_r} \theta_i^j r_i - \frac{\partial J_j^x}{\partial x} - \frac{1}{z} \frac{\partial (z J_j^z)}{\partial z} + \frac{\partial}{\partial x} \left( \xi_j \int_0^x \frac{1}{\rho} \frac{\partial \rho}{\partial t} dx \right) \quad (16)$$

in Cartesian and cylindrical coordinates, respectively. Here,  $t$  is time and  $N_r$  is the number of reactions.  $\theta_i^j$  is a stoichiometric coefficient, which is negative when component  $j$  is  $i$ -th reaction reactant and positive when it is this reaction's product. The second and third terms on the right-hand-side of Equations 15 and 16 account for gas transfer; they are present only if component  $j$  is a gas. The last terms on the right-hand-side of both equations arise due to an application of the Eulerian (stationary) coordinate framework to a medium that contracts or expands in response to density changes. Here it is assumed that the contraction or expansion occurs in one dimension ( $x$ ), with respect to a stationary plane defined by  $x = 0$ .

The conservation of energy is given by

$$\begin{aligned}
\sum_{j=1}^N \xi_j c_j \frac{\partial T}{\partial t} = \\
\sum_{i=1}^{N_r} h_i r_i - \frac{\partial q_x}{\partial x} - \frac{\partial q_y}{\partial y} - \frac{\partial I_{ex}}{\partial x} + \frac{\partial I_{rr}}{\partial x} - \sum_{g=1}^{N_g} c_g \left( J_g^x \frac{\partial T}{\partial x} + J_g^y \frac{\partial T}{\partial y} \right) + c \rho \frac{\partial T}{\partial x} \int_0^x \frac{1}{\rho} \frac{\partial \rho}{\partial t} dx
\end{aligned} \tag{17}$$

and

$$\begin{aligned}
\sum_{j=1}^N \xi_j c_j \frac{\partial T}{\partial t} = \\
\sum_{i=1}^{N_r} h_i r_i - \frac{\partial q_x}{\partial x} - \frac{1}{z} \frac{\partial(zq_z)}{\partial z} - \frac{\partial I_{ex}}{\partial x} + \frac{\partial I_{rr}}{\partial x} - \sum_{g=1}^{N_g} c_g \left( J_g^x \frac{\partial T}{\partial x} + J_g^z \frac{\partial T}{\partial z} \right) + c \rho \frac{\partial T}{\partial x} \int_0^x \frac{1}{\rho} \frac{\partial \rho}{\partial t} dx
\end{aligned} \tag{18}$$

in Cartesian and cylindrical coordinates, respectively. These equations rely on the assumptions that the radiative exchange between the material object and environment (defined by  $I_{ex}$  and  $I_{rr}$ ) is one-dimensional and that the rate of local heat transfer between gases and condensed phase components is infinitely fast. The latter assumption is implemented in the convective heat flow terms (sixth right-hand-side terms in Equations 17 and 18). An approximate finite rate correction to this assumption can be made by scaling gas heat capacities ( $c_g$ ) downward from their actual values. The last right-hand-side terms in Equations 17 and 18 account for the heat transfer associated with contraction or expansion of the material object.

## 2. Boundary Conditions

The boundary condition specification in ThermaKin2Ds is subdivided into 3 object geometries 1D, 2D and 2Dax. 1D object can be used to setup a one-dimensional pyrolysis or combustion simulation. This geometry can also be used to setup zero-dimensional, infinitely fast condensed-phase transport simulation mimicking conditions in thermal analysis instruments such as thermogravimetric analyzer (TGA). 2D object can be used to setup a two-dimensional pyrolysis or combustion simulation in Cartesian coordinates. This geometry was developed to simulate flame spread on a solid surface. 2Dax object can be used to setup a two-dimensional pyrolysis or combustion simulation in cylindrical coordinates. This geometry was developed to simulate conditions in the Controlled Atmosphere Pyrolysis Apparatus II (CAPA II) [7].

All these geometries share the same description of mass transport at the boundary. The mass flux of component  $i$  out of a material object is expressed as

$$J_i^B = \begin{cases} a_i \rho_i^B \left( \frac{\xi_i^B}{\rho_i^B} - b_i \right) \\ \text{or} \\ a_i \exp\left(-\frac{b_i}{RT^B}\right) \end{cases} \quad (19)$$

where B superscript indicates physical quantities at the boundary; and  $a_i$  and  $b_i$  are component-specific constants. The primary function of the linear expression is to remove/introduce gases from/to a pyrolyzing material. In this case,  $b_i$  represents volumetric fraction of gas  $i$  in the environment. The exponential expression is added to enable simulation of a surface reaction (such as oxidation) or material loss through dripping. In that case,  $b_i$  is the effective activation energy of this process.

### 2.1 1D Object

The 1D object is defined using the Cartesian coordinate  $x$  and 2 boundaries, top and bottom. The boundaries are defined separately using the same mathematical framework.  $x = 0$  corresponds to the position of the bottom boundary. This boundary serves as a reference point for the expansion or contraction of 1D object.

The convective heat flux out of the object boundary is expressed as

$$q^B = h_c (T^B - T^e) \quad (20)$$

where  $h_c$  is the convection coefficient and  $T^e$  is the environmental temperature that is defined through a time dependent heating rate:

$$T^e = T_0^e + \int_0^t \frac{dT^e}{dt'} dt' \quad (21)$$

$$\frac{dT^e}{dt'} = h_T^0 \left( 1 - \exp\left(-h_T^1 t'\right) \left( \cos\left(h_T^2 t'\right) + h_T^3 \sin\left(h_T^2 t'\right) \right) \right) \quad (22)$$

Here,  $T_0^e$  and  $h_T$  are user-specified parameters. Note that the numerical superscripts in Equation 22 do not indicate power; they are used to distinguish between different  $h_T$  parameters. This expression for  $T^e$  has been shown to capture sample temperature histories observed in thermal analysis experiments with a high degree of accuracy [8].

The external radiative heat flux incident onto a boundary,  $I_{ex}^B$ , is expressed through a piecewise-linear function of time. This function has 2 time segments, each defined by the initial heat flux

value at the beginning of the segment, rate of change and duration of the segment. When the time of the simulation exceeds the total time of both segments, the heat flux assumes the value attained at the end of the second segment. It is also possible to specify this piecewise-linear heat flux dependence to be repeated in time until the end of the simulation. The external radiative heat flux through the boundary ( $I_{ex}^0$ ) is computed by taking into account material reflectivity:

$$I_{ex}^0 = I_{ex}^B \frac{\int_0^X \sum_{i=1}^N \varepsilon_i v_i I_{ex} \sum_{j=1}^N \alpha_j \xi_j dx}{\int_0^X I_{ex} \sum_{j=1}^N \alpha_j \xi_j dx} \quad (23)$$

Here,  $X$  represents the thickness of the object.

A simple surface flame model is implemented as an additional 1D object boundary condition. This model relies on the following ignition/extinction criterion:

$$CI = \sum_{i=1}^N \frac{J_i^B}{J_i^{CI}} \quad (24)$$

where  $J_i^{CI}$  is the user-specified critical (or ignition) mass flux of component  $i$ . Note that only positive  $J_i^B$  are counted (i.e., component flows into the object do not contribute to the criterion value). If  $CI$  reaches or exceeds 1, a constant user-defined value is added to  $I_{ex}^B$  to simulate radiation from the flame. Simultaneously, the current values of  $h_c$  and  $T^e$  in Equation 20 are replaced with new user-defined constant values simulating the flame's convective heat feedback.

### 2.2 2D Object

The 2D object is defined using the Cartesian coordinates  $x$  and  $y$ , and 2 boundaries, front and back. These boundaries are parallel to each other and orthogonal to the coordinate  $x$ .  $x = 0$  corresponds to the position of the back boundary. This boundary serves as a reference point for the expansion or contraction of the 2D object. Two other boundaries, top and bottom, which are orthogonal to the coordinate  $y$ , are set to be impenetrable to mass or heat flow.  $y = 0$  corresponds to the position of the bottom boundary. It should be noted that, either due to initial geometry specification and/or due to material expansion or contraction during pyrolysis, the front boundary may assume a non-uniform profile (i.e., the object may have a variable thickness). In the 2D object boundary definition, this non-uniformity is ignored; it is assumed that the front boundary surface area is defined by its projection onto  $x = 0$  plane. This is not the case for the 2Dax object as explained in the next subsection.

For the 2D object, the external heat flux (radiative or convective) is defined using a piecewise-linear spatial distribution function

$$D = \begin{cases} D_0^1 + D_y^1 y & y < Y_1 \\ D_0^2 + D_y^2 y & Y_1 \leq y < Y_2 \\ D_0^3 + D_y^3 y & Y_2 \leq y < Y_3 \\ 0 & y \geq Y_3 \end{cases} \quad (25)$$

where all  $D$  and  $Y$  parameters on the right-hand-side are user-specified constants. Note that the numerical superscripts do not indicate power; they are used to distinguish between different  $D$  parameters. The function  $D$  can be specified to represent radiative heat flux incident onto a boundary ( $I_{ex}^B$ ) or environmental temperature ( $T^e$ ), which defines convective heat flux out of the object in accordance with Equation 20. The latter definition requires that a constant convection coefficient ( $h_c$ ) is also provided.

The spatial distribution of the heat flux defined by the function  $D$  can also be varied in time. The user is required to specify the times when the external heat flux exposure begins and ends. The value of the heat flux can also be either ramped up linearly in time (from 0 to the value defined by Equation 25 and, in the case of convection,  $h_c$ ), ramped down (from the value defined by Equation 25 to 0) or held steady (at the value defined by Equation 25). Three separate modules for this external heat flux specification are provided for each 2D object boundary (front and back) to gain an ability to simulate complex thermal environments.

The surface flame heat feedback formulation implemented for the 2D object is based on detailed measurements conducted on buoyancy-driven laminar flames spreading upward on vertical flat plates [9] and buoyancy-driven turbulent flames spreading upward in a corner [10]. It is believed that this formulation is sufficiently flexible to accommodate other flame spread scenarios. According to the measurement results, the values of parameters used in the flame heat feedback formulation depend strongly on the geometry of the problem and weakly on the composition of the burning solid.

At the core of this formulation is an expression for the height of the flame

$$Y_f = Y_f^0 + Y_f^{CI} \left( \int_0^Y CI dy \right)^{Y_f^{PW}} \quad (26)$$

where all  $Y_f$  parameters on the right-hand-side are user-specified constants;  $CI$  is the ignition/extinction criterion defined by Equation 24; and  $Y$  and is the distance from the bottom to

the top boundary. The flame height value is used not only as the main length scale for the heat feedback but also as a global ignition/extinction criterion. The value of  $Y_f^0$  is usually negative [9][10] and the flame heat feedback is added only when the computed flame height becomes positive. Additional flexibility is achieved by specifying the earliest time when the flame heat feedback can be turned on. This feature helps define boundary conditions for the scenarios where the surface flame is initiated by an external burner exposure and the burner and surface flames coexist for a certain period of time [10].

The flame heat flux (radiative or convective) is defined using a spatial distribution of the form

$$D_f = \begin{cases} D_f^1 & 0 \leq (y - y_b) < Y_f \text{ and } (y - y_b) < Y_f^s \\ D_f^2 & 0 \leq (y - y_b) < Y_f \text{ and } (y - y_b) \geq Y_f^s \\ D_f^1 \exp(-D_f^3(y - y_b)^2) & (y - y_b) < 0 \\ D_f^1 D_f^4 \exp\left(-\ln(D_f^4) \left(\frac{y - y_b + Y_f^t}{Y_f + Y_f^t}\right)^2\right) & (y - y_b) \geq Y_f \text{ and } (y - y_b) < Y_f^s \\ D_f^2 D_f^4 \exp\left(-\ln(D_f^4) \left(\frac{y - y_b + Y_f^t}{Y_f + Y_f^t}\right)^2\right) & (y - y_b) \geq Y_f \text{ and } (y - y_b) \geq Y_f^s \end{cases} \quad (27)$$

where  $D_f$  with numerical superscripts and  $Y_f$  with letter superscripts are used-specified constants. Here,  $y_b$  designates the position of the base of the flame that is defined as the closest point to the bottom boundary where  $CI \geq 1$ . As in the case of external heat flux definition, function  $D_f$  can be specified to represent incident radiative heat flux or environmental temperature ( $T^e$ ), which defines convective heat flux in conjunction with a constant value of  $h_c$  specified for the flame. The flame heat flux and all external heat fluxes are added together. The total incident radiation from the flame and external sources ( $I_{ex}^B$ ) is converted to the radiation through the boundary ( $I_{ex}^0$ ) in accordance with Equation 23.

The last parameter used in the definition of the 2D object is background temperature,  $T^b$ .  $T^b$  represents equilibrium temperature of the object in the absence of external heating or flame. The value of this parameter is added to all values of  $D$  and  $D_f$  when these functions are used to define convective heat transfer. The value of this parameter is also used to compute background radiative heat flux,  $\sigma(T^b)^4$ , which is added to  $I_{ex}^B$ .

### 2.3 2Dax Object

The 2Dax object is defined using the cylindrical coordinates  $x$  and  $z$ , and 2 boundaries, front and back.  $x = 0$  corresponds to the position of the back boundary, which is orthogonal to the coordinate  $x$ . This boundary serves as a reference point for the expansion or contraction of the 2Dax object.  $z = 0$  defines to the center axis of the object. It should be noted that, either due to initial geometry specification and/or due to material expansion or contraction during pyrolysis, the front boundary may assume a non-uniform profile (i.e., the object may have thickness that varies with the radial coordinate  $z$ ). In the 2Dax object boundary definition, the impacts of this non-uniformity on the surface area and heat exposure of the front boundary are taken into account.

Analogous to the 2D object, the external heat flux (radiative or convective) is defined using a piecewise-linear spatial distribution function

$$L = \begin{cases} L_0^1 + L_z^1 z & z < Z_1 \\ L_0^2 + L_z^2 z & Z_1 \leq z < Z_2 \\ L_0^3 + L_z^3 z & Z_2 \leq z < Z_3 \\ 0 & z \geq Z_3 \end{cases} \quad (28)$$

where all  $L$  and  $Z$  parameters on the right-hand-side are user-specified constants. Note that the numerical superscripts do not indicate power; they are used to distinguish between different  $L$  parameters. The function  $L$  can be specified to represent radiative heat flux incident onto a boundary ( $I_{ex}^B$ ) or environmental temperature ( $T^e$ ), which defines convective heat flux out of the object in accordance with Equation 20. The latter definition requires that a constant convection coefficient ( $h_c$ ) is also provided.

The spatial distribution of the heat flux defined by the function  $L$  can also be varied in time. The user is required to specify the times when the external heat flux exposure begins and ends. The value of the heat flux can also be either ramped up linearly in time (from 0 to the value defined by Equation 28 and, in the case of convection,  $h_c$ ), ramped down (from the value defined by Equation 28 to 0) or held steady (at the value defined by Equation 28). Three separate modules for this external heat flux specification are provided for each 2Dax object boundary to gain an ability to simulate complex thermal environments.

An additional external heat flux module is provided to enable modeling of pyrolysis experiments in the CAPA II [7]. The parameters for the convective heat flux expressed by Equation 20 are defined in this module as follows:

$$T^e = T_1^e \exp(T_2^e t) + T_3^e \exp(T_4^e t) \quad (29)$$

$$h_c = h_c^0 + h_c^1 z \quad (30)$$

Here, all  $T^e$  with numerical subscripts and  $h_c$  with numerical superscripts represent user-specified constants.

The incident radiative heat flux is defined in this module as a function of the object surface distance from a reference plane  $x = x_r$ , object surface radial position ( $z$ ) and absolute value  $\phi$  of the angle, expressed in radians, between the vector normal to the object surface and vector normal to the  $x = x_r$  plane:

$$\begin{aligned} I_{ex}^{B-2Dax} = I_0 & \left( 1 + I_1^x (x - x_r) + I_2^x (x - x_r)^2 \right) \times \\ & \left( 1 + \left( I_0^{xz} + I_1^{xz} (x - x_r) + I_2^{xz} (x - x_r)^2 \right) z + \left( I_0^{xzz} + I_1^{xzz} (x - x_r) + I_2^{xzz} (x - x_r)^2 \right) z^2 \right) \times \\ & \left( 1 + \left( I_0^{x\phi} + I_1^{x\phi} (x - x_r) + I_2^{x\phi} (x - x_r)^2 \right) \phi + \left( I_0^{x\phi\phi} + I_1^{x\phi\phi} (x - x_r) + I_2^{x\phi\phi} (x - x_r)^2 \right) \phi^2 \right) \end{aligned} \quad (31)$$

All  $I$  parameters on the right-hand-side represent user-specified constants.  $I_0$  is the heat flux set point. This is the heat flux measured at a physical location defined by  $x = x_r$ ,  $z = 0$  and  $\phi = 0$ . Note that the numerical superscripts used with  $I$  parameters do not indicate power; they are used to distinguish between different constant values.

Equation 31 is valid for  $x \geq x_r$ . For the cases where  $x < x_r$ , the incident radiative heat flux becomes

$$I_{ex}^{B-2Dax} = I_0 \left( 1 + I_0^{xz} z + I_0^{xzz} z^2 \right) \left( 1 + I_0^{x\phi} \phi + I_0^{x\phi\phi} \phi^2 \right) \quad (32)$$

The heat fluxes defined by the CAPA II module can be turned on and off at user-specified times. All external and CAPA II heat fluxes are added together. The total incident radiative heat flux is converted to the heat flux through the boundary ( $I_{ex}^0$ ) using Equation 23.

As in the case of 2D object, specification of the heat flow through each 2Dax boundary requires definition of background temperature,  $T^b$ .  $T^b$  represents equilibrium temperature of the object in the absence of external heating. The value of this parameter is added to  $T^e$  (defined by Equation 29) and  $L$  function, when this function represents the environmental temperature. The value of this parameter is also used to compute background radiative heat flux,  $\sigma(T^b)^4$ , which is added to the sum of the incident radiative fluxes.

## 3. Solution Methodology

To solve the conservation equations, a material object is divided into finite volumes (or elements). For the 1D object, the elements are characterized by size  $\Delta x$ . For two-dimensional objects, the elements are characterized by two sizes:  $\Delta x$  and  $\Delta y$ , in the case of 2D object, and  $\Delta x$  and  $\Delta z$ , in the case of 2Dax object. These sizes represent discretization with respect to the corresponding coordinates.

Each element is also characterized by masses of components and temperature. These parameters serve as primary object descriptors. Changes in these descriptors with time are computed using a small timestep,  $\Delta t$ . For  $x$  dimension of any object, the time integration is based on the Crank-Nicolson scheme [11]. A detailed description of this integration procedure can be found in an earlier publication [2]. For  $y$  and  $z$  dimensions, a simple explicit integration is used. The  $x$  dimension of the 2D and 2Dax objects is expected to be associated with higher temperature and concentration gradients; therefore, a more stable, semi-implicit integration technique is applied to this coordinate.  $x$  is also the dimension associated with the radiative heat transport and object deformation (see Equations 15-18), which add “stiffness” to the conservation equations.

Previously performed simulations [1] indicate that the integration in the  $x$  dimension is stable at a wide range of selected values of  $\Delta x$  and  $\Delta t$ . The accuracy of the results gradually decreases with increasing value of either of these parameters. Therefore, to ensure that a given simulation is converged, a sensitivity of the results to the values of the integration parameters should be always examined.  $\Delta x = 5 \times 10^{-5}$  m and  $\Delta t = 0.01$  s represent a good initial guess for a typical pyrolysis problem. The stability of the solution in the explicit integration dimension ( $y$  or  $z$ ) is a subject of a well-defined condition:

$$\frac{(\Delta y)^2}{\Delta t} > 2A \quad (33)$$

where  $A$  represents either the gas transfer coefficient ( $\lambda$ ) or thermal diffusivity,  $\frac{k}{\rho c}$ , of the object, whichever is greater.

It should be noted that, as a result of object deformation,  $\Delta x$  of individual elements changes with time. These changes, accumulated over time, may have substantial negative effects on the accuracy and stability of the solution process. To minimize these effects, element sizes are adjusted before every timestep. If an element is larger than a pre-set value of  $\Delta x$ , it is split in two. If it is smaller, a fraction of the following element is added to bring it to the pre-set value. The temperature and composition of the mixed element are recalculated to ensure the conservation of energy and species.

To simplify the solution process, the radiation-related terms (terms 4 and 5 on the right-hand-side of Equations 17 and 18) are incorporated into to the energy conservation using either the maximum absorption or random absorption algorithm. In the case of the maximum absorption algorithm, the element that, according to Equation 11, absorbs most of the external radiation is assumed to absorb all of it (corrected for reflection). In the case of the random absorption algorithm, a Monte Carlo

approach is used to distribute the radiative energy. The absorbing element is selected at random using the distribution of energy obtained from the solution of Equation 11 as a probability density guiding this selection. In both cases, the element selection is performed at every timestep. Re-radiation of energy to the environment is computed using temperature and emissivity of the same element that is selected to absorb (and reflect) external radiation.

ThermaKin2Ds is implemented using ANSI/ISO C++ and its standard library. A single instance of this program is designed to run on a single core of a single processor. Each subroutine of ThermaKin2Ds has been carefully verified; the results of these verification exercises are reported elsewhere [1][2][12].

## 4. Examples of Simulations: Input and Output

ThermaKin2Ds is compiled into a console application as a single executable file, ThermaKin2Ds.exe. When the program starts, it asks for 3 file names. The information on components and reactions is read from the components (.cmp) file. The information on the initial state of a material object, boundary conditions and integration parameters is read from the conditions (.cnd) file. The last name is that of a new file where the calculation results are stored. All files are written in ASCII text format. The input is case sensitive. An example of components file is provided in Figure 1. This figure also contains an explanation of the input parameters. Examples of conditions files specific to particular testing scenarios are provided in the following subsections.

### 4.1 *Thermogravimetric Analysis, Differential Scanning Calorimetry, and Microscale Combustion Calorimetry*

Thermal analysis experiments, including thermogravimetric analysis (TGA), differential scanning calorimetry (DSC), and microscale combustion calorimetry (MCC) [13], can be simulated using 1D object geometry. The simulation is setup by representing a material sample with a single element. The environmental temperature is defined using a time dependent heating rate measured in a given experiment and a high convection coefficient ( $h_c = 1 \times 10^5 \text{ W m}^{-2} \text{ K}^{-1}$ ) to ensure that the simulated sample temperature perfectly follows the experimental temperature profile. The mass transport at the boundary is prescribed to quickly remove all gases forming inside the sample. These settings are an idealization of the actual experiment. In this idealization, the rate of both mass and heat transport through the sample are defined to be effectively infinitely fast. An example of conditions file used to simulate a thermal analysis experiment in an inert gaseous environment is provided in Figure 2.

```

//Extruded clear poly(methylmethacrylate) or PMMA produced by Evonik
//Leventon I.T., Li J., Stollari S.I., Combustion and Flame, vol. 162, pp.
3884-3895 (2015)

COMPONENT:      PMMA
STATE:          S
DENSITY:        1155 0 0 0
HEAT CAPACITY: 601 3.63 0 0
CONDUCTIVITY:   0.45 -3.8e-4 0 0
TRANSPORT:      2e-5 0 0 0
EMISSIVITY & ABSORPTION: 0.95 1.94

COMPONENT:      PMMA_glass
STATE:          S
DENSITY:        1155 0 0 0
HEAT CAPACITY: 601 3.63 0 0
CONDUCTIVITY:   0.27 -2.4e-4 0 0
TRANSPORT:      2e-5 0 0 0
EMISSIVITY & ABSORPTION: 0.95 1.94

COMPONENT:      PMMA_res
STATE:          S
DENSITY:        1155 0 0 0
HEAT CAPACITY: 601 3.63 0 0
CONDUCTIVITY:   0.27 -2.4e-4 0 0
TRANSPORT:      2e-5 0 0 0
EMISSIVITY & ABSORPTION: 0.95 1.94

COMPONENT:      PMMA_gas
STATE:          G
DENSITY:        1155 0 0 0
HEAT CAPACITY: 1800 0 0 0
CONDUCTIVITY:   0.27 -2.4e-4 0 0
TRANSPORT:      2e-5 0 0 0
EMISSIVITY & ABSORPTION: 0.95 1.94

MIXTURES
S SWELLING:      0
L SWELLING:      0
G SWELLING LIMIT: 1e-30
PARALL CONDUCTIVITY: 0.5
PARALL TRANSPORT: 0.5

REACTION:      PMMA + NOCOMP -> PMMA_glass + NOCOMP
STOICHIOMETRY: 1 0 1 0
ARRHENIUS:    1 0
HEAT:          0 0 0 0
TEMP LIMIT:    L 378

REACTION:      PMMA_glass + NOCOMP -> PMMA_res + PMMA_gas
STOICHIOMETRY: 1 0 0.015 0.985
ARRHENIUS:    8.6e12 188100
HEAT:          -846000 0 0 0
TEMP LIMIT:    L 378

```

← Optional comments on the file contents can be provided at the beginning or the end

← Each component description starts with its name followed by its state or category: S – for solid, L – for liquid or G – for gas

← Each of these component properties is defined using a set of coefficients,  $p_0, p_1, p_2$ , and  $n$ , of Equation 1. These coefficients are defined so that Equation 1 yields corresponding properties in SI units:  $\text{kg m}^{-3}$  – for density,  $\text{J kg}^{-1} \text{K}^{-1}$  – for heat capacity,  $\text{W m}^{-1} \text{K}^{-1}$  – for thermal conductivity, and  $\text{m}^2 \text{s}^{-1}$  – for gas transfer coefficient (TRANSPORT)

← Component emissivity (equal to 1 – reflectivity) and absorption coefficient in  $\text{m}^2 \text{kg}^{-1}$ , which is normalized by component density

← Indicates lower (L) temperature limit in K below which the reaction is turned off. Used in conjunction with  $A = 1 \text{ s}^{-1}$  and  $E = 0$  to define a fast (thermodynamically-controlled) phase transition: at temperature above 378 K, all PMMA is quickly converted to PMMA\_glass. The reverse transition can also be defined using the upper (U) temperature limit

← These inputs define dimensionless parameters,  $\beta$  in Equation 10, used to compute thermal conductivity and gas transfer coefficient for a mixture of a given composition

← Dimensionless parameters,  $\gamma_s, \gamma_l$  and  $\tau$ , specifying reaction of material volume to the presence of gases. See Equations 3 and 4. Current (default) parameter values represent the case where gases do not contribute to the volume

← Reaction is defined by specifying names of reactants and products and absolute values of their stoichiometric coefficients. NOCOMP indicates that the reaction does not have second reactant or product

← Pre-exponential factor in  $\text{s}^{-1}$  or  $\text{m}^3 \text{kg}^{-1} \text{s}^{-1}$ , activation energy in  $\text{J mol}^{-1}$ , and heat of reaction in  $\text{J kg}^{-1}$  defined using a set of coefficients,  $p_0, p_1, p_2$ , and  $n$ , of Equation 1. Negative heat indicates endotherm

**Figure 1.** PMMA.cmp – an example of components file used to define material properties.

```

//Example of Thermal Analysis Experiment Simulation (nominal heating rate is 10
K/min)

OBJECT TYPE: 1D ← Define object geometry:
                1D, 2D or 2Dax

OBJECT STRUCTURE
*****

THICKNESS: 1e-5 ← For 1D object, the initial object structure is specified by defining
TEMPERATURE: 313      layers of specified thickness in m, temperature in K, and composition
MASS FRACTIONS:         in component mass fractions. The first layer specified is adjacent to
PMMA 1.0                the top boundary; the last – to the bottom boundary. Here, only one
                        layer consisting of PMMA component is specified

OBJECT BOUNDARIES
*****

TOP BOUNDARY
MASS TRANSPORT: YES ← YES indicates that gases are allowed to flow through this boundary.
PMMA_gas LIN 0.05 0 ← NO at both boundaries completely turns off gas transport (across the
                        boundaries and inside the object)
OUTSIDE INIT TEMP: 313
OUTSIDE HEAT RATE: 0.167 0.00302 0.00458 -0.992 ← LIN indicates that the linear expression of Equation 19 is used. EXP
CONVECTION COEFF: 1e5                           should be specified to switch to the exponential expression. The
                        following values are the expression coefficients,  $a_i$  and  $b_i$ , used for
                        the gaseous component PMMA_gas. These are typical (default)
                        coefficient values employed to specify rapid gas removal at the
                        boundary

EXTERNAL RADIATION: NO ← These are parameters of Equation 22 specifying heating rate
                        evolution in time. The first value is the nominal or set heating
                        rate in K s-1. The rest of the parameters account for the fact that
                        this heating rate is not achieved instantaneously due to thermal
                        inertia of the system. The values of these parameters are
                        instrument and heating rate specific. To discard effects of thermal
                        inertia, the last 3 parameters should be set to 1010, 0, and 0

FLAME: NO ← The external radiation and flame model are turned off at both
                        boundaries

BOTTOM BOUNDARY
MASS TRANSPORT: NO
OUTSIDE INIT TEMP: 313
OUTSIDE HEAT RATE: 0 0 0 0
CONVECTION COEFF: 0
EXTERNAL RADIATION: NO
FLAME: NO

INTEGRATION PARAMETERS
*****

ELEMENT SIZE: 1e-5 ← The element size (m) is the same as the object size. I.e., the object
TIME STEP: 0.01 ← is represented by one element. The integration time step is 0.01 s
DURATION: 2800 ← (default value). The duration of the simulated experiment, also in
                        s, is defined here to reach a prescribed maximum temperature (≈
                        780 K)

OUTPUT FREQUENCY:
ELEMENTS: 1 ← Data output will be provided for each element, once per second
TIME STEPS: 100

```

**Figure 2.** TA\_PMMA.cnd – an example of conditions file used to simulate a thermal analysis experiment.

To perform this type of simulation correctly, all radiative losses from the object surfaces should be turned off. This is accomplished by setting emissivity values of all components (specified in the components file) to 0. A portion of the output file generated by ThermaKin2Ds using the input listed in Figure 2 is shown in Figure 3. The format of the output is self-explanatory. This output

can be easily imported into a spreadsheet program, such as Microsoft Excel, by identifying the source as a tab-separated text file.

```

ThermaKin Program Version 7

Components file: TA_PMMA.cmp
Conditions file: TA_PMMA.cnd

Number of components: 4
Number of reactions: 2
Mixture rules assigned: yes

Object type: 1D

Top Boundary
External radiation: off
Mass transport: on
Ignition: off

Bottom Boundary
External radiation: off
Mass transport: off
Ignition: off

Started on: Wed Mar 13 14:09:26 2019

*****
Time [s] = 0.000000e+00

BOUNDARY      AREA [m^2]      HEAT FLOW IN [W]      MASS FLOW OUT [kg/s]:      PMMA      PMMA_glass      PMMA_res      PMMA_gas
TOP           1.000000e+00      0.000000e+00      0.000000e+00      0.000000e+00      0.000000e+00      0.000000e+00      0.000000e+00
BOTTOM        1.000000e+00      0.000000e+00      0.000000e+00      0.000000e+00      0.000000e+00      0.000000e+00      0.000000e+00

FROM TOP [m]  TEMPERATURE [K]      CONCENTRATION [kg/m^3]:      PMMA      PMMA_glass      PMMA_res      PMMA_gas
5.000000e-06      3.130000e+02      1.155000e+03      0.000000e+00      0.000000e+00      0.000000e+00

Total thickness [m] = 1.000000e-05
Total mass [kg/m^2] = 1.155000e-02

*****
Time [s] = 1.000000e+00

BOUNDARY      AREA [m^2]      HEAT FLOW IN [W]      MASS FLOW OUT [kg/s]:      PMMA      PMMA_glass      PMMA_res      PMMA_gas
TOP           1.000000e+00      2.531193e-02      0.000000e+00      0.000000e+00      0.000000e+00      0.000000e+00      0.000000e+00
BOTTOM        1.000000e+00      0.000000e+00      0.000000e+00      0.000000e+00      0.000000e+00      0.000000e+00      0.000000e+00

FROM TOP [m]  TEMPERATURE [K]      CONCENTRATION [kg/m^3]:      PMMA      PMMA_glass      PMMA_res      PMMA_gas
5.000000e-06      3.130006e+02      1.155000e+03      0.000000e+00      0.000000e+00      0.000000e+00      0.000000e+00

Total thickness [m] = 1.000000e-05
Total mass [kg/m^2] = 1.155000e-02
...

```

**Figure 3.** A small portion of ThermaKin2Ds output file, TA\_PMMA.txt, containing results of a thermal analysis experiment simulation.

A MathWorks MATLAB R2017a script, ThermaKin\_Visualizer\_v2.m, is also available to quickly display the results of the thermal analysis simulations, compare them with experimental data and/or save a summary of these results in an Excel file. A screenshot of this interactive script is shown in Figure 4. The summary file generated using this script provides a range of data relevant to the thermal analysis experiments including the rates of production of individual gaseous components normalized by the initial sample mass. Multiplying these rates by the corresponding heats of combustion and adding them together produces heat release rate (HRR) [8] – a simulated result of MCC experiment, which can be presented either as a function of sample temperature or time. ThermaKin\_Visualizer\_v2.m can be used to convert ThermaKin2Ds output to MCC HRR curve and save the curve in an Excel file or compare it to experimental data. To perform this

conversion, the user needs switch to the MCC pane, load the ThermaKin2Ds output and enter the heats of combustion of the gaseous components.

![Figure 4: A screenshot of the ThermaKin_Visualizer_v2.m MATLAB script interface. The interface is divided into two main sections: 'Inputs' on the left and 'Graphs' on the right. The 'Inputs' section includes a dropdown menu for 'TGA/DSC', buttons for 'Experimental Data Input' (Browse), 'Thermakin Output' (Browse), and 'TA_PMMA.txt was loaded', a 'Data Cursor' button, a 'Plot' button, and a 'Save Thermakin Output As' button. The 'Graphs' section displays four plots: 1. 'Mass Loss Rate' (y-axis: Mass Loss Rate [1/s] x 10^-3, x-axis: Temperature [K]) showing a sharp peak around 650 K. 2. 'Total Mass Normalized by Initial Mass' (y-axis: Mass / Initial Mass, x-axis: Temperature [K]) showing a step-like decrease starting around 600 K. 3. 'Heat Flow (endotherm is positive)' (y-axis: Heat Flow [W/g], x-axis: Temperature [K]) showing a peak around 650 K. 4. 'Integral Heat Flow (endotherm is positive)' (y-axis: Integral Heat Flow [J/g], x-axis: Temperature [K]) showing a step-like increase starting around 600 K. All plots show a single red line representing the simulation results.](3ae74a33759ae31781f484406db4feed_img.jpg)

Figure 4: A screenshot of the ThermaKin\_Visualizer\_v2.m MATLAB script interface. The interface is divided into two main sections: 'Inputs' on the left and 'Graphs' on the right. The 'Inputs' section includes a dropdown menu for 'TGA/DSC', buttons for 'Experimental Data Input' (Browse), 'Thermakin Output' (Browse), and 'TA\_PMMA.txt was loaded', a 'Data Cursor' button, a 'Plot' button, and a 'Save Thermakin Output As' button. The 'Graphs' section displays four plots: 1. 'Mass Loss Rate' (y-axis: Mass Loss Rate [1/s] x 10^-3, x-axis: Temperature [K]) showing a sharp peak around 650 K. 2. 'Total Mass Normalized by Initial Mass' (y-axis: Mass / Initial Mass, x-axis: Temperature [K]) showing a step-like decrease starting around 600 K. 3. 'Heat Flow (endotherm is positive)' (y-axis: Heat Flow [W/g], x-axis: Temperature [K]) showing a peak around 650 K. 4. 'Integral Heat Flow (endotherm is positive)' (y-axis: Integral Heat Flow [J/g], x-axis: Temperature [K]) showing a step-like increase starting around 600 K. All plots show a single red line representing the simulation results.

**Figure 4.** Graphs generated by ThermaKin\_Visualizer\_v2.m MATLAB script used to post-process the thermal analysis simulation results stored in TA\_PMMA.txt file.

### 4.2 Cone Calorimetry

Modeling cone calorimetry [14] requires specification of the heat flux from the radiant heater and the flame that forms on the top surface of the sample upon ignition. Recently performed measurements of the flame heat flux indicate that its spatial distribution across the sample surface is highly non-uniform [15]. It was shown that, to model the heat feedback from the flame accurately, the surface should be portioned into 2 zones. Here we present a simplified and less accurate approach where the surface-area-averaged value of the flame heat flux is used within the framework of a single 1D simulation. Details of the averaged (single-zone) flame model can be found in the Appendix A (Supplementary data) of [15].

Cone calorimetry samples are typically placed onto a thick thermal insulator, which provides nearly adiabatic conditions at the back boundary. In practice, the insulating materials are never perfect and even small heat losses to them may affect the outcome of the test. Therefore, in the presented simulation, the insulator is modeled explicitly by defining 2 additional components, KAWOOL and KAWOOL\_back. The first component represents the properties of Kaowool PM, a thermally insulating ceramic fiber board manufactured to close tolerances. The second component has the same properties with the exception of emissivity; the emissivity is set to 0 to prevent radiative losses. A thin layer of this component is added at the bottom of the insulation

layer to create a perfectly adiabatic boundary. These component definitions, which are shown in Figure 5, are added to the components file. The gas transfer coefficients for these components are defined to be effectively zero ( $1 \times 10^{-30} \text{ m}^2 \text{ s}^{-1}$ ) to prevent movement of gaseous components into the insulation.

```
...  
  
COMPONENT:      KAOWOOL  
STATE:          S  
DENSITY:        256 0 0 0  
HEAT CAPACITY:  1070 0 0 0  
CONDUCTIVITY:   0.0519 -4e-5 1e-7 2  
TRANSPORT:      1e-30 0 0 0  
EMISSIVITY & ABSORPTION: 0.95 10000  
  
COMPONENT:      KAOWOOL_back  
STATE:          S  
DENSITY:        256 0 0 0  
HEAT CAPACITY:  1070 0 0 0  
CONDUCTIVITY:   0.0519 -4e-5 1e-7 2  
TRANSPORT:      1e-30 0 0 0  
EMISSIVITY & ABSORPTION: 0 10000  
  
...
```

**Figure 5.** Definitions of components representing thermal insulation used in cone calorimetry. These definitions were added to CC\_PMMA.cmp file.

An example of conditions file used to simulate a cone calorimetry experiment performed at  $49.7 \text{ kW m}^{-2}$  radiant heater heat flux (this is the average heat flux across the sample top surface) on a 6 mm thick sample of PMMA is provided in Figure 6. The results of this simulation can be processed and visualized using the Cone pane of ThermaKin\_Visualizer\_v2.m. A screenshot demonstrating this process is shown in Figure 7.

As in the case of the MCC simulation visualization, the user is required to enter the heats of combustion of the gaseous components. In the example shown on the figure, the script is employed to compare the results of the simulations with the experimental HRR and mass loss rate (MLR) data uploaded from CC\_PMMA\_ExperimentalData.xlsx file. Note that the script also allows to extract up to 5 temporal temperature profiles (TC1 through TC5) at fixed distances from the bottom boundary. If the specified distance exceeds the thickness of the object, the temperature at the top boundary (corresponding, in this case, to the sample surface facing the heater) is displayed and recorded.

```

//Example of Cone Calorimetry Simulation (cone heater is set at 49.7 kW/m2)
//Heat of combustion of PMMA_gas is 2.55e7 J/kg
//PMMA_gas ignition mass flux is 2.1e4 W/m2 / 2.55 J/kg = 8.24e-4 kg/m2-s
//For further information, see McCoy C. G. et al., Fire Safety Journal, vol. 103,
pp. 38-48 (2019)

OBJECT TYPE: 1D

OBJECT STRUCTURE
*****
THICKNESS: 0.006
TEMPERATURE: 300
MASS FRACTIONS:
PMMA 1.0

THICKNESS: 0.0171
TEMPERATURE: 300
MASS FRACTIONS:
KAOWOOL 1.0

THICKNESS: 0.0019
TEMPERATURE: 300
MASS FRACTIONS:
KAOWOOL_back 1.0

OBJECT BOUNDARIES
*****

TOP BOUNDARY

MASS TRANSPORT: YES
PMMA_gas LIN 0.05 0

OUTSIDE INIT TEMP: 300
OUTSIDE HEAT RATE: 0 0 0 0
CONVECTION COEFF: 9.5

EXTERNAL RADIATION: YES
TIME PROG1: 49.7e3 0 1e4
TIME PROG2: 0 0 0
REPEAT: NO
ABSORPTION MODE: RAND

FLAME: YES
IGNITION MASS FLUXS:
PMMA_gas 8.24e-4
OUTSIDE TEMP: 2150
CONVECTION COEFF: 15.3
RADIATION: 2.8e3

BOTTOM BOUNDARY
...

```

0.006 m thick PMMA sample at the top; 0.0171 m thick ceramic fiber insulation at the bottom. A thin layer of KAOWOOL\_back at the very bottom has zero emissivity and is used to prevent radiative losses from the bottom boundary

The external radiation model is turned on and used to represent the radiation from the cone heater

RAND means that the random absorption algorithm is used to handle radiation. MAX should be used to switch to the maximum absorption algorithm

The model of the heat feedback from the surface flame is also turned on

Constant environmental temperature in K, which matches the initial temperature of the sample, and convection coefficient in  $\text{W m}^{-2} \text{K}^{-1}$  used to define convective losses from the top sample surface prior to ignition

The incident radiative heat flux is held constant at its initial value,  $4.97\text{e}4 \text{ W m}^{-2}$ , for  $1\text{e}4 \text{ s}$  or duration of the simulation

Rate of change of the heat flux in  $\text{W m}^{-2} \text{s}^{-1}$

Second time segment of the heat flux program that starts immediately after the first one

Critical mass flux in  $\text{kg m}^{-2} \text{s}^{-1}$ ; defines ignition according to Equation 24

Environmental temperature and convection coefficient that define convective heat flow at the top boundary after ignition

Radiative heat flux that is added to one specified above after ignition

The bottom boundary and integration parameter specifications are not shown on this figure. They are similar to those explained in Figure 2

**Figure 6.** CC\_PMMA.cnd – an example of conditions file used to simulate a cone calorimetry experiment.

![Figure 7: A screenshot of the ThermaKin_Visualizer_v2.m MATLAB interface. The interface includes an 'Inputs' panel on the left with options for 'Cone' simulation, 'Experimental Data Input' (loading CC_PMMA_ExperimentalData.xlsx), 'Thermakin Output' (loading CC_PMMA.txt), and 'Data Cursor'. The 'Graphs' panel displays four plots: 1) 'Heat Release Rate' (kW/m²) vs. 'Time [s]' (0-300), showing 'Simulation Total' (red solid line) and 'Experimental' (black dashed line) data. A note states 'Ignition is assumed to occur when virtual HRR exceeds 21 kW/m²'. 2) 'Mass Loss Rate' (kg/m²/s) vs. 'Time [s]' (0-300), comparing 'Simulation Total' (red solid) and 'Experimental' (black dashed). 3) 'In-depth Temperatures' (K) vs. 'Time [s]' (0-300), showing 'Model TC1' through 'Model TC5' (solid lines) and 'Exp TC1' through 'Exp TC5' (dashed lines). 4) 'Thickness' (m) vs. 'Time [s]' (0-300), showing 'Simulation' (red solid) and 'Experimental' (black dashed). A table on the left shows 'dHc (kJ/g)' for 'PMMA_gas' as 25.500. Another table lists 'Position from bottom [m]' for models TC1 to TC5: TC1=0.0250, TC2=0.0190, TC3=0, TC4=0, TC5=0.](60e9207be66a64332619bb4b667fe67b_img.jpg)

Figure 7: A screenshot of the ThermaKin\_Visualizer\_v2.m MATLAB interface. The interface includes an 'Inputs' panel on the left with options for 'Cone' simulation, 'Experimental Data Input' (loading CC\_PMMA\_ExperimentalData.xlsx), 'Thermakin Output' (loading CC\_PMMA.txt), and 'Data Cursor'. The 'Graphs' panel displays four plots: 1) 'Heat Release Rate' (kW/m²) vs. 'Time [s]' (0-300), showing 'Simulation Total' (red solid line) and 'Experimental' (black dashed line) data. A note states 'Ignition is assumed to occur when virtual HRR exceeds 21 kW/m²'. 2) 'Mass Loss Rate' (kg/m²/s) vs. 'Time [s]' (0-300), comparing 'Simulation Total' (red solid) and 'Experimental' (black dashed). 3) 'In-depth Temperatures' (K) vs. 'Time [s]' (0-300), showing 'Model TC1' through 'Model TC5' (solid lines) and 'Exp TC1' through 'Exp TC5' (dashed lines). 4) 'Thickness' (m) vs. 'Time [s]' (0-300), showing 'Simulation' (red solid) and 'Experimental' (black dashed). A table on the left shows 'dHc (kJ/g)' for 'PMMA\_gas' as 25.500. Another table lists 'Position from bottom [m]' for models TC1 to TC5: TC1=0.0250, TC2=0.0190, TC3=0, TC4=0, TC5=0.

**Figure 7.** Graphs generated by ThermaKin\_Visualizer\_v2.m MATLAB script used to post-process the cone calorimetry simulation results stored in CC\_PMMA.txt file and compare these results with the experimental data stored in CC\_PMMA\_ExperimentalData.xlsx file.

## References

1. Stoliarov S. I., Leventon I. T., Lyon R. E., *Two-dimensional Model of Burning for Pyrolyzable Solids*, Fire and Materials, vol. 38, pp. 391-408 (2014).
2. Stoliarov S. I., Lyon R. E., *Thermo-Kinetic Model of Burning*, Federal Aviation Administration Technical Note, DOT/FAA/AR-TN08/17 (2008).
3. Blinder S. M., *Advanced Physical Chemistry: A Survey of Modern Theoretical Principles*, The Macmillan Company, London, UK, 1969.
4. Holman J. P., *Heat Transfer, Ninth Edition*, McGraw-Hill, Boston, MA, 2002.
5. Siegel R., Howell J., *Thermal Radiation Heat Transfer, Fourth Edition*, Taylor & Francis, New York, NY, 2002.
6. Scheidegger A. E., *The Physics of Flow through Porous Media*, University of Toronto Press, Toronto, Canada, 1974.
7. Swann J. D., Ding Y., McKinnon M. B., Stoliarov S. I., *Controlled Atmosphere Pyrolysis Apparatus II (CAPA II): A New Tool for Analysis of Pyrolysis of Charring and Intumescent Polymers*, Fire Safety Journal, vol. 91, pp. 130-139 (2017).

8. McKinnon M. B., Stoliarov S. I., *Pyrolysis Model Development for a Multilayer Floor Covering, Materials*, vol. 8, pp. 6117-6153 (2015).
9. Leventon I. T., Korver K. T., Stoliarov S. I., *A Generalized Model of Flame to Surface Heat Feedback for Laminar Wall Flames*, *Combustion and Flame*, vol. 179, pp. 338-353 (2017).
10. Lannon C. M., Stoliarov S. I., Lord J. M., Leventon I. T., *A Methodology for Predicting and Comparing the Full-scale Fire Performance of Similar Materials based on Small-scale Testing*, *Fire and Materials*, vol. 42, pp. 710-724 (2018).
11. Press W. H., Teukolsky S. A., Vetterling W. T., Flannery B. P., *Numerical Recipes in C++: The Art of Scientific Computing*, Cambridge University Press, Cambridge, UK, 2002.
12. Swann J. D., Ding Y., Stoliarov S. I., *Characterization of Pyrolysis and Combustion of Rigid Poly(vinyl chloride) using Two-dimensional Modeling*, *International Journal of Heat and Mass Transfer*, vol. 132, pp. 347-361 (2019).
13. Lyon R. E., Walters R. N., Stoliarov S. I., Safronava N., *Principles and Practice of Microscale Combustion Calorimetry*, Federal Aviation Administration Report, DOT/FAA/TC-12/53 (2013).
14. ASTM E1354-16a, *Standard Test Method for Heat and Visible Smoke Release Rates for Materials and Products Using an Oxygen Consumption Calorimeter*, ASTM International, West Conshohocken, 2016.
15. McCoy C. G., Tilles J. L., Stoliarov S. I., *Empirical Model of Flame Heat Feedback for Simulation of Cone Calorimetry*, *Fire Safety Journal*, vol. 103, pp. 38-48 (2019).