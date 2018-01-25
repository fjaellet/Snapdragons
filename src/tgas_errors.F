c ---------------------------------------------------------------------
c
c Subroutine that computes and assign TGAS errors for TYCHO-2 only subset and HIPPARCOS subset
c to observables. Errors are computed following the paper Michalik et al. 2015
c 
c Note: No correlation is assumed between astrometric errors. No dependency with scaning law 
c is assumed. The error for mu_alpha, mu_delta is imposed to be the same as mu as an 
c approximation.
c
c The code allows to compute the astrometric errors.
c
c Universitat de Barcelona
c Contact: Roger Mor, Merce Romero-Gomez
c email: rmor@am.ub.es
c
c This is the First Version and work is still in progress:
c Work is in progress to include the correlations
c Work is in progress to add errors as a function of G magnitude
c Future Work: Include dependencies on the scanning law
c
c ---------------------------------------------------------------------
c
c Updated: September 2015
c
c Input values:
c -------------
c     a(i): Astrometric true values for the star (not affected by errors)
c         a(1): Equatorial heliocentric right ascension (units: radians)
c         a(2): Equatorial heliocentric declination: delta (units: radians) 
c         a(3): Parallax: pi (units: mas)
c         a(4): Equatorial proper motions in right ascension in true arcs on the sky: mua_{*}=mua*cos(delta) (units: mas/yr)
c         a(5): Equatorial proper motions: mud (units: mas/yr)
c
c     V_T: V tycho magnitude. (TGAS is not going to contain new magnitudes for the stars)
c*********************************************************************
c     subset: INTEGER Variable to select the TGAS subset             *
c--------------------------------------------------------------------*
c                                                                    *     
c            subset=1 : For Tycho Only subset of TGAS                *
c            subset=2 : For HIPPARCOS subset of TGAS                 * 
c*********************************************************************
c
c Output values: 
c -------------
c Parameters of the star expected in TGAS (affected by errors)
c TGAS errors assigned to each parameter:  
c         TGASER(1): Standard deviation in Equatorial heliocentric right ascension in true arcs on the sky:: alpha_{*}=alpha*cos(delta)  (units: mas)
c         TGASER(2): Standard deviation in Equatorial heliocentric declination: delta (units: mas) 
c         TGASER(3): Standard deviation in Parallax: pi (units: mas)
c         TGASER(4): Standard deviation in Equatorial proper motions in right ascension in true arcs on the sky: mua_{*}=mua*cos(delta) (units: mas/yr)
c         TGASER (5): Standard deviation in Equatorial proper motions: mud (units: mas/yr)
c   Astrometric data affected by errors  
c     ao(i): Astrometric values for the star affected by Gaia errors
c         TGASOB(1): Observed Equatorial heliocentric right ascension (units: radians)
c         TGASOB(2): Observed Equatorial heliocentric declination: delta (units: radians) 
c         TGASOB(3): Observed Parallax: pi (units: mas)
c         TGASOB(4): Observed Equatorial proper motions in right ascension in true arcs on the sky: mua_{*}=mua*cos(delta) (units: mas/yr)
c         TGASOB(5): Observed Equatorial proper motions: mud (units: mas/yr)

      Subroutine TGASErrors(V_T,a,TGASOB,TGASER, subset,Dseed)
  
      Double Precision V_T,frac
      Double Precision a(5),TGASER(5), TGASOB(5)
      Integer Dseed, subset
      external gasdev
      double precision gasdev

      INCLUDE 'const_math.h'
      INCLUDE 'const_ast.h'

c      Dseed=time() ! number of seconds from 1970 (If you prefere a fixed seed use the line below)
c      Dseed=214749.d0 
c# Previous line can be substituted by this one if desired.

c *** Computation of the observed astrometric quantities 
c *** The error in right ascension TGASER(1) denotes true arc on the sky
c *** So the right ascension shall be converted to that before the 
c *** random error is assigned 
c *** alpha_{*}=alpha*cos(delta) 
      a(1)=a(1)*dcos(a(2))
c *** Conversion of (alpha*,delta) from radians to mas
      a(1)=a(1)/mas
      a(2)=a(2)/mas
c Here we assign the values for the standard deviation depending on magnitude and on subset.
      if(subset.eq.1) then
         if(V_T.lt.7.0d0) then
            TGASER(1)=0.244d0
            TGASER(2)=0.244d0
            TGASER(3)= 0.399d0
            TGASER(4)= 0.198d0
            TGASER(5)= 0.198d0
         elseif(V_T.ge.7.0d0.and.V_T.lt.7.5d0) then
            frac=(V_T-7.0d0)/(7.5d0-7.0d0)
            TGASER(1)=0.244d0-frac*0.244d0+frac*0.198d0
            TGASER(2)=0.244d0-frac*0.244d0+frac*0.198d0
            TGASER(3)=0.399d0-frac*0.399d0+frac*0.348d0
            TGASER(4)=0.198d0-frac*0.198d0+frac*0.264d0
            TGASER(5)=0.198d0-frac*0.198d0+frac*0.264d0
         elseif(V_T.ge.7.5d0.and.V_T.lt.8.5d0) then 
            frac=(V_T-7.5d0)
            TGASER(1)=0.198d0-frac*0.198d0+frac*0.191d0
            TGASER(2)=0.198d0-frac*0.198d0+frac*0.191d0
            TGASER(3)=0.348d0-frac*0.348d0+frac*0.327d0
            TGASER(4)=0.264d0-frac*0.264d0+frac*0.403d0
            TGASER(5)=0.264d0-frac*0.264d0+frac*0.403d0
         elseif(V_T.ge.8.5d0.and.V_T.lt.9.5d0) then
            frac=(V_T-8.5d0)
            TGASER(1)=0.191d0-frac*0.191d0+frac*0.230d0
            TGASER(2)=0.191d0-frac*0.191d0+frac*0.230d0
            TGASER(3)=0.327d0-frac*0.327d0+frac*0.407d0
            TGASER(4)=0.403d0-frac*0.403d0+frac*0.680d0
            TGASER(5)=0.403d0-frac*0.403d0+frac*0.680d0
         elseif(V_T.ge.9.5d0.and.V_T.lt.10.5d0) then 
            frac=(V_T-9.5d0)
            TGASER(1)=0.230d0-frac*0.230d0+frac*0.329d0
            TGASER(2)=0.230d0-frac*0.230d0+frac*0.329d0
            TGASER(3)=0.407d0-frac*0.407d0+frac*0.601d0
            TGASER(4)=0.680d0-frac*0.680d0+frac*1.145d0
            TGASER(5)=0.680d0-frac*0.680d0+frac*1.145d0
         elseif(V_T.ge.10.5d0.and.V_T.lt.11.5d0) then 
            frac=(V_T-10.5d0)
            TGASER(1)=0.329d0-frac*0.329d0+frac*0.379d0
            TGASER(2)=0.329d0-frac*0.329d0+frac*0.379d0
            TGASER(3)=0.601d0-frac*0.601d0+frac*0.722d0
            TGASER(4)=1.145d0-frac*1.145d0+frac*1.522d0
            TGASER(5)=1.145d0-frac*1.145d0+frac*1.522d0
         elseif(V_T.ge.11.5d0.and.V_T.lt.12.0d0) then
            frac=(V_T-11.5d0)/(12.0d0-11.5d0)
            TGASER(1)=0.379d0-frac*0.379d0+frac*0.349d0
            TGASER(2)=0.379d0-frac*0.379d0+frac*0.349d0
            TGASER(3)=0.722d0-frac*0.722d0+frac*0.702d0
            TGASER(4)=1.522d0-frac*1.522d0+frac*1.615d0
            TGASER(5)=1.522d0-frac*1.522d0+frac*1.615d0
         elseif(V_T.ge.12.0d0) then
            TGASER(1)=0.349d0
            TGASER(2)=0.349d0
            TGASER(3)=0.702d0
            TGASER(4)=1.615d0
            TGASER(5)=1.615d0
         endif  
      elseif(subset.eq.2) then
         if(V_T.lt.7.0d0) then
            TGASER(1)=0.116d0
            TGASER(2)=0.116d0
            TGASER(3)=0.180d0
            TGASER(4)=0.017d0
            TGASER(5)=0.017d0
         elseif(V_T.ge.7.0d0.and.V_T.lt.7.5d0) then
            frac=(V_T-7.0d0)/(7.5d0-7.0d0)
            TGASER(1)=0.116d0-frac*0.116d0+frac*0.120d0
            TGASER(2)=0.116d0-frac*0.116d0+frac*0.120d0
            TGASER(3)=0.180d0-frac*0.180d0+frac*0.192d0
            TGASER(4)=0.017d0-frac*0.017d0+frac*0.021d0
            TGASER(5)=0.017d0-frac*0.017d0+frac*0.021d0
         elseif(V_T.ge.7.5d0.and.V_T.lt.8.5d0) then 
            frac=(V_T-7.5d0)
            TGASER(1)=0.120d0-frac*0.120d0+frac*0.125d0
            TGASER(2)=0.120d0-frac*0.120d0+frac*0.125d0
            TGASER(3)=0.192d0-frac*0.192d0+frac*0.198d0
            TGASER(4)=0.021d0-frac*0.021d0+frac*0.029d0
            TGASER(5)=0.021d0-frac*0.021d0+frac*0.029d0
         elseif(V_T.ge.8.5d0.and.V_T.lt.9.5d0) then
            frac=(V_T-8.5d0)
            TGASER(1)=0.125d0-frac*0.125d0+frac*0.133d0
            TGASER(2)=0.125d0-frac*0.125d0+fraC*0.133d0
            TGASER(3)=0.198d0-frac*0.198d0+frac*0.217d0
            TGASER(4)=0.021d0-frac*0.021d0+frac*0.039d0
            TGASER(5)=0.021d0-frac*0.021d0+frac*0.039d0
         elseif(V_T.ge.9.5d0.and.V_T.lt.10.5d0) then 
            frac=(V_T-9.5d0)
            TGASER(1)=0.133d0-frac*0.133d0+frac*0.154d0
            TGASER(2)=0.133d0-frac*0.133d0+frac*0.154d0
            TGASER(3)=0.217d0-frac*0.217d0+frac*0.253d0
            TGASER(4)=0.039d0-frac*0.039d0+frac*0.058d0
            TGASER(5)=0.039d0-frac*0.039d0+frac*0.058d0
         elseif(V_T.ge.10.5d0.and.V_T.lt.11.5d0) then 
            frac=(V_T-10.5d0)
            TGASER(1)=0.154d0-frac*0.154d0+frac*0.128d0
            TGASER(2)=0.154d0-frac*0.154d0+frac*0.128d0
            TGASER(3)=0.253d0-frac*0.253d0+frac*0.211d0
            TGASER(4)=0.058d0-frac*0.058d0+frac*0.087d0
            TGASER(5)=0.058d0-frac*0.058d0+frac*0.087d0
         elseif(V_T.ge.11.5d0.and.V_T.lt.12.0d0) then
            frac=(V_T-11.5d0)/(12.0d0-11.5d0)
            TGASER(1)=0.128d0-frac*0.128d0+frac*0.151d0
            TGASER(2)=0.128d0-frac*0.128d0+frac*0.151d0
            TGASER(3)=0.211d0-frac*0.211d0+frac*0.248d0
            TGASER(4)=0.087d0-frac*0.087d0+frac*0.135d0
            TGASER(5)=0.087d0-frac*0.087d0+frac*0.135d0
         elseif(V_T.ge.12.0d0) then
            TGASER(1)=0.151d0
            TGASER(2)=0.151d0
            TGASER(3)=0.248d0
            TGASER(4)=0.135d0
            TGASER(5)=0.135d0
         endif  
      endif
      do i=1,5
         TGASOB(i)=gasdev(a(i),TGASER(i),Dseed)
      enddo
      
c *** The true and observed (alpha*,delta) are converted from mas to radians 
      a(1)=a(1)*mas
      a(2)=a(2)*mas
      TGASOB(1)=TGASOB(1)*mas
      TGASOB(2)=TGASOB(2)*mas
c *** and the alpha_{*} is converted to alpha
      TGASOB(1)=TGASOB(1)/dcos(a(2))
      a(1)=a(1)/dcos(a(2))

      return
      end

c ****************************************************************
c Subrountine Tcoord
c ****************************************************************
      subroutine Tcoord(Tao,Ta,xo_p,yo_p,zo_p,lo_p,bo_p
     &     ,vlo_p,vbo_p,l_p,b_p,vl_p,vb_p)

      implicit none
      double precision SI_tmateqgal(0:2,0:2),SI_tmatgaleq(0:2,0:2)
      double precision d,xeq,yeq,zeq,rxy,xo_p,yo_p,zo_p,temp
     &     ,vxo_p,vyo_p,vzo_p,Uo,Vo,Wo,A(3),r(3),Tao(5),b(6),vl_p
     &     ,vb_p,vlo_p,vbo_p,l1,l2,l3,c1,c2,c3,vl2,vlo2,Tao4,Tao5
     &     ,sinphi,phi,lo_p,bo_p,l_p,b_p,Ta(5),Ta4,Ta5
     &     ,cosphi,angp,dngp
      INCLUDE 'const_math.h'

c galactic (Cartesian, Solar origin) -> heliocentlic (Cartesian)
      SI_tmateqgal(0,0)=-0.054875290d0
      SI_tmateqgal(1,0)=0.494109674d0
      SI_tmateqgal(2,0)=-0.867666026d0
      SI_tmateqgal(0,1)=-0.873437003d0
      SI_tmateqgal(1,1)=-0.444830323d0
      SI_tmateqgal(2,1)=-0.198075198d0
      SI_tmateqgal(0,2)=-0.483834719d0
      SI_tmateqgal(1,2)=0.746981441d0
      SI_tmateqgal(2,2)=0.455985407d0
c heliocentlic (Cartesian) -> galactic (Cartesian, Solar origin)
      SI_tmatgaleq(0,0)=-0.054877639d0
      SI_tmatgaleq(1,0)=-0.873437219d0
      SI_tmatgaleq(2,0)=-0.483835031d0
      SI_tmatgaleq(0,1)=0.494110185d0
      SI_tmatgaleq(1,1)=-0.444829280d0
      SI_tmatgaleq(2,1)=0.746982180d0
      SI_tmatgaleq(0,2)=-0.867665586d0
      SI_tmatgaleq(1,2)=-0.198076592d0
      SI_tmatgaleq(2,2)=0.455983866d0
c *** Calculate catesian error added positions ***
      d=1.0d0/(Tao(3)*1.0e5)
      rxy=d*dcos(Tao(2))
      xeq=rxy*dcos(Tao(1))
      yeq=rxy*dsin(Tao(1))
      zeq=d*dsin(Tao(2))
      xo_p=SI_tmatgaleq(0,0)*xeq+SI_tmatgaleq(1,0)*yeq
     &     +SI_tmatgaleq(2,0)*zeq
      yo_p=SI_tmatgaleq(0,1)*xeq+SI_tmatgaleq(1,1)*yeq
     &     +SI_tmatgaleq(2,1)*zeq
      zo_p=SI_tmatgaleq(0,2)*xeq+SI_tmatgaleq(1,2)*yeq
     &     +SI_tmatgaleq(2,2)*zeq

      c1=dcos(Tao(5)*1000.0d0*mas)*dcos(Tao(4)*1000.0d0*mas)
      c2=dcos(Tao(5)*1000.0d0*mas)*dsin(Tao(4)*1000.0d0*mas)
      c3=dsin(Tao(5)*1000.0d0*mas)
      l1=SI_tmatgaleq(0,0)*c1+SI_tmatgaleq(1,0)*c2+
     &     SI_tmatgaleq(2,0)*c3
      l2=SI_tmatgaleq(0,1)*c1+SI_tmatgaleq(1,1)*c2+
     &     SI_tmatgaleq(2,1)*c3
      l3=SI_tmatgaleq(0,2)*c1+SI_tmatgaleq(1,2)*c2+
     &     SI_tmatgaleq(2,2)*c3
c *** Jo's equations for mu(a,d) to mu(l,b) ***
      angp=192.85948d0*0.0174533d0
      dngp=27.12825d0*0.0174533d0
      call equtogal(Ta(1),Ta(2),l_p,b_p)
      if (l_p.lt.0.0d0) then
         l_p=l_p+2.0d0*pi
      endif
      cosphi=(dsin(dngp)-dsin(Ta(2))*dsin(b_p))/
     &     (dcos(Ta(2))*dcos(b_p))
      sinphi=(dsin(Ta(1)-angp)*dcos(dngp))/(dcos(b_p))
      vl_p=cosphi*Ta(4)+sinphi*Ta(5)
      vb_p=cosphi*Ta(5)-sinphi*Ta(4)
c *** observed mu(a,d) -> mu(l,b) calculation ***
      call equtogal(Tao(1),Tao(2),lo_p,bo_p)
      if (lo_p.lt.0.0d0) then
         lo_p=lo_p+2.0d0*pi
      endif
      cosphi=(dsin(dngp)-dsin(Tao(2))*dsin(bo_p))/
     &     (dcos(Tao(2))*dcos(bo_p))
      sinphi=(dsin(Tao(1)-angp)*dcos(dngp))/(dcos(bo_p))
      vlo_p=cosphi*Tao(4)+sinphi*Tao(5)
      vbo_p=cosphi*Tao(5)-sinphi*Tao(4)

c      Ta4=Ta(4)*1000.0d0*mas
c      Ta5=Ta(5)*1000.0d0*mas
c      call equtogal(Ta(1),Ta(2),l_p,b_p)
c      if (l_p.lt.0.0d0) then
c         l_p=l_p+2.0d0*pi
c      endif
c      sinphi=0.40893064d0*dsin(pi/2.0d0+Ta(1))/(pi/2.0d0-b_p)
c      phi=dasin(sinphi)
cc *** Quadrant loop for phi ***
c      if (dsin(sinphi).lt.0.0d0) then
c         if (dcos(sinphi).gt.0.0d0) then
c            phi=2.0d0*pi+dasin(sinphi)
c         endif
c      endif   
c      if (dcos(sinphi).lt.0.0d0) then
c         phi=pi-dasin(sinphi)
c      endif
c      
c      vl_p=(Ta4*dcos(phi)+Ta5*dsin(phi))/(1000.0d0*mas)
c      vb_p=(Ta5*dcos(phi)-Ta4*dsin(phi))/(1000.0d0*mas)
c      write(6,*) 'Mul=',vlo_p,'Mub=',vbo_p

c *** vbo_p quadrant loop ***
c      if (dsin(l3).lt.0.0d0) then
c         if (dcos(l3).gt.0.0d0) then   
c            vbo_p=vbo_p+2.0d0*pi
c         else if (dcos(l3).lt.0.0d0) then
c            vbo_p=pi-vbo_p
c         endif
c      else if (dsin(l3).gt.0.0d0) then
c         if (dcos(l3).lt.0.0d0) then
c            vbo_p=pi-vbo_p
c         endif
c      endif
cc *** vlo_p quadrant loop ***
c      temp=l1/dcos(vbo_p)
c      if (dsin(temp).lt.0.0d0) then
c         vlo_p=2.0d0*pi-dacos(l1/dcos(vbo_p))
c      else if (dsin(temp).gt.0.0d0) then
c         vlo_p=dacos(l1/dcos(vbo_p))
c      endif   
cc *** vlo2 quadrant loop (should match vlo_p)
c      temp=l2/dcos(vbo_p)
c      if (dsin(temp).lt.0.0d0) then
c         if (dcos(temp).lt.0.0d0) then
c            vlo2=pi-dasin(temp)
c         else if (dcos(temp).gt.0.0d0) then
c            vlo2=2.0d0*pi+dasin(temp)
c         endif
c      else if (dsin(temp).gt.0.0d0) then
c         if (dcos(temp).gt.0.0d0) then
c            vlo2=dasin(temp)
c         else if (dcos(temp).lt.0.0d0) then
c            vlo2=pi-dasin(temp)
c         endif
c      endif

c      write(6,*) vlo_p,vbo_p,vlo2
c      vlo_p=(vlo_p/mas)/1000.0d0
c      vbo_p=(vbo_p/mas)/1000.0d0

c      write(6,*) vlo_p,vlo2,vbo_p,l1,l2,l3,dcos(vbo_p),
c     &     l2/dcos(vbo_p),l1/dcos(vbo_p)
c It'll be something to do with the dacos2 or whatever . . .
c incorrect treatment of signs *sigh* again
     
      return
      end

c subroutine equatorial_to_galactic(alfa,delta,gl,gb)
c input: equatorial coordinates (alfa,delta)
c output: galactic coordinates (xl,xb)
c Units: all angles l,b,alpha,delta in radians
c -------------------------------------------------------------

      subroutine equtogal(alfa,delta,gl,gb)
      implicit real*8(a-h,o-z)
      real*8 lq,lpl
      include 'const_math.h'
      include 'const_ast.h'
c      INCLUDE 'const_pot.h'

      T6=60.d0
      T36=3600.d0
      RD=PI/180.d0
      LQ=122.92851d0*RD
      ALFAG=(12.d0+51.d0/T6+26.2755d0/T36)*15.d0*RD
      DECG=(27.d0+7.d0/T6+41.704d0/T36)*RD
      SDG=dSIN(DECG)
      CDG=dCOS(DECG)
 
      sd=dsin(delta)
      cd=dcos(delta)

      sb=sdg*sd+cdg*cd*dcos(alfa-alfag)
      cb=sqrt(1.d0-sb**2.d0)
      gb=datan2(sb,cb)
      slpl=(cd*dsin(alfa-alfag))/cb
      clpl=(cdg*sd-sdg*cd*dcos(alfa-alfag))/cb
      lpl=datan2(slpl,clpl)
      if(lpl.le.0.d0)lpl=lpl+2.d0*pi
      gl=lq-lpl

      return
      end
