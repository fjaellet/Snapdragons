c----------------------------------------------------------------------
c  This is 'const_ast.h'
c  Astronomical and physical constants (SI units)
c  ep0 is the reference epoch for the astrometric parameters,
c  reckoned from J2000.0
  
      double precision pc,Gte,time
c en cm
      PARAMETER (pc = 3.086d18,
c en  dines cm2/g2 (cm3/(s2 g)) (Bowers)
     &           Gte = 6.670d-8,
c canvi de km/s/kpc a 1/(100 Myr)  
     &           time=0.10225d0)

      double precision kt,a1,a2,Ag,Dg
      parameter (kt=4.7404705d0)
      parameter (a1=(62.87124882d0*deg))
      parameter (a2=(32.93680516d0*deg))
      parameter (Ag=(282.8594813d0*deg))
      parameter (Dg=(27.12825111d0*deg))


c  Constants related to galactic potential

      REAL*8 Mb,Md,Mh,K

c values of the potential
      PARAMETER (ZMsol = 1.989d33,! Msol in gm
     +           abu = 0.387d0,! en Kpc
     +           ad = 5.3178d0,
     +           ah = 12.0d0,
     +           bd = 0.25d0,
     +           Mb = 1.406d10,!in Msol
     +           Md = 8.561d10,!in Msol
     +           Mh = 1.071d11)!in Msol

c other variables:
       PARAMETER(K=Gte*ZMsol/pc/1.d13,
     +           abuabu=abu*abu,
     +           bdbd=bd*bd,
     +           rKMb=K*Mb,
     +           q9=-3.d0*rKMb,
     +           rKMd=K*Md,
     +           q8=-3.d0*rKMd,
     +           rKMh=K*Mh)
