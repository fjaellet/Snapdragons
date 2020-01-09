c     Parameters for stellar splitting and mock observations
c     generation from cosmological simulations
c     Santi Roca-Fabrega 31-08-2016
      
c ***  Parameters for the IMF ***
c ** Flag for the IMF to use
      integer imftype

      parameter (imftype=1) ! 1=salpeter

c ** Salpeter xIMF=1.35, Msmin=0.1d0 **
c ** Theoretical work for salpeter http://astro1.physics.utoledo.edu/~megeath/ph6820/lecture12_eqn.pdf

       double precision xIMF, Msmin

       parameter (xIMF=1.35d0) !IMF Slope
       parameter (Msmin=0.1d0) ! M1 parameter

c *** conversion units factors **

       double precision VUKMS,LUKPC

       parameter (VUKMS=207.4d0)
       parameter (LUKPC=100.0d0)

