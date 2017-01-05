Below are some examples of visualizations created with Co-Blade:


![alt tag](https://raw.githubusercontent.com/nnmrec/Co-Blade/master/Documentation/images/Co-Blade logo.PNG)

Figure 1. Replica of the Sandia SNL100-00 100-meter 13.2MW wind turbine blade, modeled in the Co-Blade software. Co-Blade is able to model a large variety of composite layups and blade topologies.


![alt tag](https://raw.githubusercontent.com/nnmrec/Co-Blade/master/Documentation/images/Co-Blade post process GUI.png)

Figure 2. This graphical user interface (GUI) allows the user to plot a 3D representation of the blade and to select via a spreadsheet which panels and/or laminas are shown in the plot. Also, the user can select which output parameter (i.e. Effective modulus, panel thickness, strain, stress, etc.) are visualized in the selected panels and/or laminas. Many types of visualizations can be created by using this GUI, and some examples are illustrated in Figure 3 below.


| (a)  | ![alt tag](https://raw.githubusercontent.com/nnmrec/Co-Blade/master/Documentation/images/Co-Blade post process Blade Geometry.png)            |
| (b)  | ![alt tag](https://raw.githubusercontent.com/nnmrec/Co-Blade/master/Documentation/images/Co-Blade post process Max Stress FC blade shell.png) |
| (c)  | ![alt tag](https://raw.githubusercontent.com/nnmrec/Co-Blade/master/Documentation/images/Co-Blade post process Max Stress FC blade root.png)  |
| (d)  | ![alt tag](https://raw.githubusercontent.com/nnmrec/Co-Blade/master/Documentation/images/Co-Blade post process Max Stress FC spar cap.png)    |
| (e)  | ![alt tag](https://raw.githubusercontent.com/nnmrec/Co-Blade/master/Documentation/images/Co-Blade post process Max Stress FC shear webs.png)  |

Figure 3. Some examples of figures that can be created to visualize Co-Blade results using the graphical user interface (GUI) described in Figure 2. This particular series of figures visualizes stresses in multiple layers of the composite blade. (a) showing only the external shape of the blade, (b) max stress failure criteria in the "blade-shell" material covering the entire top surface of the blade, (c) max stress failure criteria in the root build-up "blade-root" material, which lies directly under the "blade-shell" material, (d) max stress failure criteria in the spar cap "spar-uni" material, which lies directly under the "blade-root" material, (e) max stress failure criteria in the shear web "web-shell" material. Values greater than 1 for the max stress failure criteria indicate that the material has exceeded its maximum allowable stress.


| (a)  | ![alt tag](https://raw.githubusercontent.com/nnmrec/Co-Blade/master/Documentation/images/Co-Blade applied loads and body forces 3D plot.png) |
| (b)  | ![alt tag](https://raw.githubusercontent.com/nnmrec/Co-Blade/master/Documentation/images/Co-Blade applied loads and body forces 2D plot.png) |

Figure 4. (a) The external blade shape (outlined in black), the displaced external blade shape (outlined in red), and the magnitude and direction of the applied forces are plotted in the view of the global coordinate system. (b) Forces (px, py, pz) and moments (qz) are applied per unit length. Aerodynamic forces and moments are applied along the blade pitch axis, while body forces act along the blade inertial axis. Additional coupling between extensional, bending, and torsion is accounted for though offsets of the blade center-of-mass, tension center, and shear center from the blade pitch axis.


| (a)  | ![alt tag](https://raw.githubusercontent.com/nnmrec/Co-Blade/master/Documentation/images/Co-Blade effective Youngs modulus.png)       |
| (b)  | ![alt tag](https://raw.githubusercontent.com/nnmrec/Co-Blade/master/Documentation/images/Co-Blade effective shear modulus.png)        |
| (c)  | ![alt tag](https://raw.githubusercontent.com/nnmrec/Co-Blade/master/Documentation/images/Co-Blade effective bending stiffnesses.png)  |

Figure 5. Subplots (a) and (b) show the different laminates in each cross section, colored by thier effective Young's modulus and effective shear modulus. Subplot (c) shows the principal flapwise and edgewise bending stiffnesses (EI). Co-Blade computes structural properties with respect to several points of reference (e.g. the mass/tension/shear center, and global axes) because different aeroelastic simulation codes (such as FAST, ADAMS, BLADED, etc.) require structural properties to be input with respect to specific coordinate systems.


| (a)  | ![alt tag](https://raw.githubusercontent.com/nnmrec/Co-Blade/master/Documentation/images/Co-Blade effective normal stress.png) |
| (b)  | ![alt tag](https://raw.githubusercontent.com/nnmrec/Co-Blade/master/Documentation/images/Co-Blade effective shear stress.png)  |
| (c)  | ![alt tag](https://raw.githubusercontent.com/nnmrec/Co-Blade/master/Documentation/images/Co-Blade buckling criteria.png)       |

Figure 6. Subplots (a) and (b) show the effective normal stress and effective shear stress in the cross section laminates. These effective stresses are converted into equivalent extensional, bending, and shear loads on a laminated plate within classical lamination theory to recover the 2D lamina-level strains and stresses. Subplot (c) shows the panel buckling criteria within the blade--panels subjected to compression and shear are prone to buckle or wrinkle, and a buckling criteria of R greater than 1 indicates that the panel has exceeded its critical buckling stress (see the user's guide for additional explanation on buckling).


| (a)  | ![alt tag](https://raw.githubusercontent.com/nnmrec/Co-Blade/master/Documentation/images/Co-Blade lamina strain e-11.png) |
| (b)  | ![alt tag](https://raw.githubusercontent.com/nnmrec/Co-Blade/master/Documentation/images/Co-Blade lamina stress s-11.png) |

Figure 7. Subplots (a) and (b) show the span variant and chordwise variant normal stresses and normal strains in direction-1 (the principal direction) within each individual lamina of the top, bottom, and web panels. The coordinate z/L is the distance measured along the blade pitch axis starting at the blade root and normalized by the blade length. The coordinate x/c is the distance measured along the chord line starting at the leading edge and normalized by the chord. The coordinate s/Lw is the distance measured along the mid-wall of the shear web starting at the end of the web connected to the top surface and normalized by the total length of the web panel. A line is plotted for every single lamina within the entire blade, and the lines are color coded by their material name. This plot can become easily cluttered, but it provides a quick reference for the range of lamina stress and strain magnitudes within the entire blade. In subplot (b), note that the stress in the "core" materials (which are structural foam) is approximately zero due to their relatively low stiffness compared to the other materials (which are high stiffness E-glass).


![alt tag](https://raw.githubusercontent.com/nnmrec/Co-Blade/master/Documentation/images/Co-Blade flapwise modal displacements 1st 3 modes.png)

Figure 8. The BModes code has been integrated with Co-Blade to predict coupled mode shapes and frequencies. In this example, the flapwise modal displacements along the blade length for the first 3 modes are shown.