# v1.0.0
Cite as: International Journal of Heat and Mass Transfer 137 (2019) 835-846
 (https://www.sciencedirect.com/science/article/pii/S0017931019304508) doi: 10.1016/j.ijheatmasstransfer.2019.03.152

# rotatingElectrodesFoam
It is described how to simulate current distribution in electrochemical reactors with rotating electrodes under mass-transfer controlled conditions, by the developed solver: concentrationPimpleFoam. It is shown how to pre-process, run and post-process a basic case involving a 2D axisimetric domain. 
The proposed strategy allows to compute the wall share stress and flux of a single specie towards electrodes.

# Disclaimer
This offering is not approved or endorsed by OpeFOAM Foundation, producer and distributor of the OpenFOAM software via www.openfoam.org.

# Usage
In applications (A) you will find the scripts to compile the solver, turbulence model and post-processing utility in order to obtain the wall share stress and mass-transfer profile.
In tutorials (B) you will find two examples for an axisymmetric rotating cylinder electrode (RCE) and an axisymmetric rotating disc electrode (RDE). 
In more information (C) you will find a brief description of the proposed tool

# #  A) Applications
**1.**  
_A)_ Paste applications/utilities/Solvers/concentrationPimpleFoam inside OpenFOAM user directory (Applications/Utilities/Solvers).  
_B)_ Open a terminal inside concentrationPimpleFoam.  
_C)_ Run wmake.  
**2.**  
_A)_ Paste applications/utilities/turbulenceModels/kOmegaSSTCC and SpalartAllmarasRCC inside OpenFOAM user directory (Applications/Utilities/turbulenceModels).  
_B)_ Open a terminal inside kOmegaSSTCC and SpalartAllmarasRCC.  
_C)_ Run wmake.  
**3.**  
_A)_ Paste applications/utilities/postProcessing/wallFlux and /wallSS inside OpenFOAM user directory (Applications/Utilities/postProcessing).  
_B)_ Open a terminal inside wallFlux and wallSS.  
_C)_ Run wmake.  

# #  B) Tutorials
**1-** Paste tutorials inside OpenFOAM user directory (Run/Tutorials).  
**2-** Enter to RDE or RCE and open a Terminal.  
**3-** Modify properties (diffusion coefficient, turbulent Schmith number and kinematic viscosity) inside constant/transportProperties.  
**4-** Modify turbulence properties (laminar or turbulent fluid flow with the corresponding turbulence model) inside constant/turbulenceProperties. 
**5-** Run ./Allrun.  

# # C) More information
**Geometry and Mesh**  
The geometry is defined and then meshed using the OpenFOAM blockMesh tool, after running the blockMesh utility. 

**Files setup by the user**  

_0 directory:_ Boundary ans initial conditions. 

_constant directory:_ As is any standard OpenFOAM case, the constant folder must contain a standard polymesh directory, generated by the standard blockMeshDict-file, which defines the full domain and its mesh and is located in the system directory. Also, must contain transport properties and the desirable turbulence model.

_system directory:_ In the system directory, there is information about the discretization procedure, and control variables of the numerical scheme. 

**Post-processing**  
By the command lines >> postProcess -func flux and >> concentrationPimpleFoam -postProcess -func shearStress it is possible to obtain by means of the post processing utility wallFlux and wallSS (supplied here) the total flux and wall shear stress per electrode.

**Running the case**  
After following each step as defined in the above tutorial. The case can be executed by using the command in the terminal ./Allrun. An Allrun can be explained as a script file which contains all the commands used to execute the case.

**Monitoring flux (mol/m^2/s) vs. iteration**  
By runnung the command "gnuplot scripts/plot_res" in a new terminal it is possible to monitor the flux over the working electrode vs. time (iteration).

