# OpenFoam-Generalised-Wall-Function
Implementation of the Generalised Wall Function by Popovac in OpenFoam v9

##Test case 2-D NACA Airfoil

##Usage:
Execute the make in the popovacGWallFunction directory

./Allwmake

Include the following in the "controlDict" to load it:

libs ("libmynutkGWF.so");

Set the following in 0/epsilon and 0/nut for wall condition type:
type            epsilonGWF;

type            nutkGWF;

Further validation, verification and work in progress...

References
--------------------------------------------
(1) M. Popovac. “Modelling and Simulation of Turbulence and Heat Transfer in Wall-
Bounded Flows”. PhD thesis. Delft, Netherlands: Technische Universiteit Delft, Oct.
2006 (cited on page i).

(2) K. H. M. Popovac. “A combined WF and ItW treatment of wall boundary conditions
for turbulent convective heat transfer”. In: 9th UK National Heat Transfer Conference,
Manchester, UK (2005) (cited on pages i, 2).

