c *** Mock generation options ***

      integer Binary,error,EXTINCT,Schlegel,idum
      double precision ChosenLim,xsun,ysun,zsun,vlsr
      double precision usun,vsun,wsun,R0,usol,vsol,wsol

c *** Output to binary (1) or ASCII (0)? ***
      Binary=1
c *** Minimum magnitude necessary ***
      ChosenLim=12.0d0
c *** Add Gaia errors? (1=yes) ***
      error=1
c *** Add Milky Way Extinction (1=yes) ***
      EXTINCT=1
c *** Extintion map (1=Schlegel, 0=2d analytic) ***
      Schlegel=1
c *** Random seed ***
      idum=-1314
c *** Set Sun position / velocity (kpc from center)***
      xsun=0.0d0
      ysun=-8.0d0
      zsun=0.0d0
      vlsr=1.1d0*VUKMS
      R0=dsqrt(xsun**2.0+ysun**2.0+zsun**2.0)
      xsun=xsun/100.0d0
      ysun=ysun/100.0d0
      zsun=zsun/100.0d0
c *** Maximum radius around solar position to generate stars
      rad_limit_kpc=10.0d0

c *** 1.1d0*VUKMS=228.14 km/s ***
      Usol = 10.00d0 ! solar velocity(Dehnen i Binney 1998)
      Vsol = 5.25d0
      Wsol = 7.17d0
      usun=Usol
      vsun=Vsol
      wsun=Wsol
