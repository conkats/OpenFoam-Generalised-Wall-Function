# OpenFoam-Generalised-Wall-Function
Implementation of the Generalised Wall Function by Popovac

-Usage:
-Compile the GAWF using ./Allwmake to add the library in the path
-Include the following:
    In Controldict: libs ("libmynutkGWF.so");
    For the wall patches definition in the dictionaries of
    0/nut: type nutkGWF;
    0/epsilon: type            epsilonGWF;

References
--------------------------------------------
(1) M. Popovac. “Modelling and Simulation of Turbulence and Heat Transfer in Wall-
Bounded Flows”. PhD thesis. Delft, Netherlands: Technische Universiteit Delft, Oct.
2006 (cited on page i).

(2) K. H. M. Popovac. “A combined WF and ItW treatment of wall boundary conditions
for turbulent convective heat transfer”. In: 9th UK National Heat Transfer Conference,
Manchester, UK (2005) (cited on pages i, 2).

