# CREPE

## CREPE - Cosmic Ray Electron Propogation codE

CREPE is a numerical model that models the cosmic ray electron propagation in radial coordinates.
CREPE solves the diffusion-loss equation in radial coordinates using a 4th order Runga Kutta method for time stepping.
4th order central differences scheme to solve the spatial derivatives and an upwind scheme to solve for the energy derivatives. Details of such algorithms can be found in the textbook; numerical recipes.

Full details of the code and physics can be found in Mulcahy et al. 2016 Astronomy and Astrophysics
http://www.aanda.org/articles/aa/abs/2016/08/aa28446-16/aa28446-16.html

At the moment this code is very focused on the needs of the author.
I will be planning to make it more user friendly in time so it can used by wider community.

To do:
- Delete unneeded sections  
- General tidy up needs to be performed
- Courant conditions work properly
- inputs are via a parset --> done but needs refinement
- adaptive length
- Energy dependence of diffusion coefficient can be switched on
