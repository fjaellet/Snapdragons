      program PopSynth
      implicit none
      integer MNSc,MN2D,MN3D,MN,MNiso
      parameter (MN=5000000) ! *** Max number of particles in input file
      parameter (MNSc=4096) !Number of grid points in extinction map Schegel
      parameter (MN3D=1024) !Number of grid points in extinction map 3D
      parameter (MN2D=1024) !Number of grid points in extinction map 2D
      parameter (MNiso=34)  !Number of isochrones in isochromes file
c from  IsoNo.txt (name of isochrone as function of metallicity)      
      integer i,niso(0:MNiso-1),Ns,j,k,jindex(0:6200),npr
     &     ,IsoNo(0:MNiso-1),l,q,inputdata,GenSel,per3
     &     ,IsoStart(0:MNiso-1),ageindex,np,m,per2,lbegin,lend,ngen
      integer dummyint,nparts
      double precision MiMax,Afac,z,PAge,vgen,igen,avgen,aigen,z2
      double precision frac,Vigen,dummy,dmod,mag2,rd,circ,z1,mtemp,mp
      double precision Age(0:MN-1),logL_Lo(0:MN-1),logG(0:MN-1),zsh,
     &     Mi(0:MN-1),Mc(0:MN-1),logT(0:MN-1),
     &     C_O(0:MN-1),VMU(0:MN-1),
     &     VMB(0:MN-1),VMV(0:MN-1),VMR(0:MN-1),VMI(0:MN-1),
     &     VMJ(0:MN-1),VMH(0:MN-1),VMK(0:MN-1),
     &     mbol(0:MN-1)
      double precision vl_p,vb_p,vlo_p,vbo_p,vrots,mtn,Nsf,z3,z4
      double precision M_hec,period,pmode,logMdot,intIMF,vcir
      double precision mtotalneed,rad_limit_kpc
      double precision d_p,Milb,mag,mgen,sumgen,R,tss,month
      double precision x_p,y_p,z_p,vx_p,vy_p,vz_p,m_p,ex,ex31,ex148
      integer il,ib,ie1,ie2,ie3,ie4,ie5,ie6,ie7,ie8
      integer il3,ib3,id3,jflag
      double precision dl,db,wl1,wl2,wb1,wb2,wd1,wd2,exS,ex3
      double precision dl3,db3,dd3,CAfactor
      double precision leS(0:MNSc-1),beS(0:MNSc-1),le3(0:MN3D-1),
     &     be3(0:MN3D-1),de3(0:99),datS(0:16777215),dat3(0:104857599)
      double precision l_p,b_p,a(6),ao(6),ae(6),p(2),po(2),pe(2), 
     &     ap(4),apo(4),ape(4)
      double precision M_PI,t0,uh,vh,wh,makms,mdkms,makmso,mdkmso
      parameter (M_PI=4.0d0*datan(1.0d0)) !pi value
      double precision alpha,delta,g,gb,gl,GRVS
      double precision Uo,Vo,Wo,xo_p,yo_p,zo_p,vxo_p,vyo_p,vzo_p,
     &     vro_p,vroto_p,rop,od_p,rp,vr_p,vrot_p
      real*8 ran2
      character*40 FN
      character*25 ONAME
      character*1 check
      character*1 hash
      include 'parameters.h'
      include 'Mockinioptions.h'
c *** GEnrating MOck DAta From Sims GEMODAFS version 1.0 31/08/16 ***

c *** For ASCII input file ini/input.dat ***

c *** TEMP particle ***
      open(20,FILE='ini/UBV/IsoNo.txt', STATUS='OLD')
      do i=0,MNiso-1
         read(20,'(I10)') IsoNo(i) ! Read isochrones metallicity*1e6
      enddo
      close(20)

      per2=0
      per3=0
      ngen=0
      GenSel=0
c *** Read extinction arrays ***
      if (EXTINCT.eq.1) then
         if (Schlegel.eq.1) then !Schlegel extinction array
c table from http://irsa.ipac.caltech.edu/applications/DUST/
            open(51,file='ini/ScEx1.bin',status='old'
     &           ,form='unformatted')
            write(6,*) 'Reading Extinction arrays'
            do i=0,16777215
               read(51) datS(i)
            enddo

c Define londitud, latitud bins (gal coords) for extinction
            dl=(2.0d0*M_PI)/4095.0d0 
            db=M_PI/4095.0d0   
            write(6,*) 'Schlegel extinction array Read'
         else
            open(51,file='ini/2dEx1.bin',status='old'
     &           ,form='unformatted') !2D analytic extinction
            write(6,*) 'Reading Extinction arrays'
            do i=0,1048575
               read(51) datS(i)
            enddo
             write(6,*) '2D analytical extinction array Read'

c Define londitud, latitud bins (gal coords) for extinction
             dl=(2.0d0*M_PI)/1023.0d0 
             db=M_PI/1023.0d0   
         endif
         close(51)
         open(51,file='ini/3dEx1.bin',status='old',form='unformatted')
         do i=0,104857599
            read(51) dat3(i)
         enddo
         close(51)
         write(6,*) '3D factor array Read'
         write(6,*) 'Extinction running'

c Define longitud, latitud and distance bins (gal coords) for extinction
c     0<l<2*pi 
         dl3=(2.0d0*M_PI)/1023.0d0  
c     -pi/2<b<pi/2, with +pi/2 for +ve
         db3=M_PI/1023.0d0    
         dd3=3.47712123/99.0d0 
      else
         write(6,*) '*** No Extinction Added ***'
      endif
c *** Read in (Padova) Isochrone Files ***   
      j=0
      k=1
      jindex(0)=0
      do i=0,MNiso-1
         write(FN,'(a18,i5.5,a4)')'ini/UBV/output_0.0',IsoNo(i),'.dat'
         open(20,FILE=FN,STATUS='OLD')
c *** 10 headers read and discarded ***
         do l=1,9
            read(20,'(A150)')
         enddo
         hash='#'
c *** Read data ***
         niso(i)=0
         write(ONAME,'(a12,i5.5,a4)')'ini/UBV/V-VI',IsoNo(i),'.dat'
         do l=0,30000
            read(20,*,end=600) check
            if (index(check,hash).ne.1) then
               backspace(20)
               read(20,*) Age(j),Mi(j),Mc(j),logL_Lo(j),logT(j)
     &              ,logG(j),mbol(j),VMU(j),VMB(j),VMV(j),VMR(j)
     &              ,VMI(j),VMJ(j),VMH(j),VMK(j),C_O(j),M_hec
     &              ,period,pmode,logMdot,intIMF
               j=j+1
               niso(i)=niso(i)+1
            else
               read(20,*,end=600) check
               jindex(k)=j
               k=k+1
            endif
         enddo 
  600    close(20)
         IsoStart(i)=j
      enddo
      write(6,*) 'Isochrones Read'
 1940 format(17(1pE13.5),I10,3(1pE13.5))
c *** Open files to read & write ***
c *** open output file (binary or ASCII) ***
      if (Binary.eq.0) then
         open(10,file='output/GeneratedStars.dat',status='unknown')
c         open(15,file='ini/input.dat',status='unknown')
      else
      open(10,file='output/GeneratedStars.bin',status='unknown'
     &     ,form='unformatted')
      endif
c *** Open input file ***
c *** For example file ***
c      open(30,file='ini/lbselsl90d5.dat',status='unknown')
c 2030 format(16(1pE13.5),I10)
c *** For simpler file ***
      open(30,file='ini/input.dat',status='unknown')
c 2030 format(9(1pE13.5))
      np=50000000 !max number of particles
      sumgen=0.0d0
      mtotalneed=0.0d0
      nparts=0
c   16 format(9(1pE13.5))  
      do q=0,np-1
         mtemp=0.0d0
c *** For example file ***
c         read(30,2030,end=603) dummy,dummy,dummy,x_p,y_p,z_p,dummy,
c     &        dummy,vx_p,vy_p,vz_p,dummy,dummy,z,tss,m_p,dummyint
c *** x,y,z(kpc),vx,vy,vz(km/s),metallicity(z),age(Gyr),mass(solar)
         read(30,*,end=603) x_p,y_p,z_p,vx_p,vy_p,vz_p,m_p,z1,z2,tss,
     +z3,z4,vcir
         z=10**z1+10**z2*10**z1
c         write(15,16) x_p,y_p,z_p,vx_p,vy_p,vz_p,z,tss,m_p
         nparts=nparts+1
         x_p=x_p/100.0d0
         y_p=y_p/100.0d0
         z_p=z_p/100.0d0
         vx_p=(vx_p+usun)/VUKMS
         vy_p=(vy_p+vlsr+vsun)/VUKMS
         vz_p=(vz_p+wsun)/VUKMS
         z=z*1e6               ! Metallicity in isochrone filename is in
c z*10^6

c *** if z is lower than the one in the first isochrone file, assigns the
c lower metallicity in files, if higher, the highest ***

         if (z.lt.100.0d0) then
            z=100.0d0
         else if (z.gt.30000.0d0) then
            z=30000.0d0
         endif 

c *** definitions of usefull variables and changes of units ***

         rp=dsqrt((x_p+xsun)**2+(y_p+ysun)**2)
         vr_p=((vx_p*(x_p+xsun)+vy_p*(y_p+ysun)+vz_p*(z_p+zsun))/rp)*
     +VUKMS
         vrot_p=(-1.0d0*((x_p+xsun)*vy_p-(y_p+ysun)*vx_p)/rp)*VUKMS
         rp=dsqrt((x_p+xsun)**2+(y_p+ysun)**2+(z_p+zsun)**2)
         d_p=dsqrt(x_p**2+y_p**2+z_p**2)
         dmod=5.0d0*(dlog10(d_p*100000.0d0)-1.0d0)! distance term for 
c        V-Mv=5(log(r)-1)
         PAge=log10(tss*1000000000.0d0) !log(age) in yrs

c *** Cut in maximum radius

       if(rp.lt.rad_limit_kpc/100.) then

c *** Find approproate isochrone segment ***
c *** Random age & metallicity if unknown ***
c     z=(ran1(idum)*0.9d0+0.1d0)*19000.0d0
c     PAge=log10(ran1(idum)*12000000000.0d0)
         ageindex=nint((PAge-6.6d0)/0.02d0)

c limits in max and min allowed ages (depend on isochrones)
         if (ageindex.lt.0) then 
            ageindex=0
         else if (ageindex.gt.181) then
            ageindex=181
         endif

c *** extinction (CTIO) ***
         if (EXTINCT.eq.1) then
c *** Make grids l,b,d ***
            if (Schlegel.eq.1) then
               do i=0,MNSc-1
                  leS(i)=i*dl
                  beS(i)=i*db-M_PI/2.0d0
               enddo
            else
            do i=0,MN2D-1
                  leS(i)=i*dl
                  beS(i)=i*db-M_PI/2.0d0
            enddo
            endif
            do i=0,MN3D-1
               le3(i)=i*dl3
               be3(i)=i*db3-M_PI/2.0d0
            enddo
            do i=0,99
               de3(i)=i*dd3-2.0d0
            enddo

c *** l,b,d from cartesian input ***
            l_p=-datan2(y_p,-x_p)+M_PI
            b_p=dasin(z_p/d_p)
            vrots=(-1.0d0*(x_p*vy_p-y_p*vx_p)/d_p)*VUKMS
            vl_p=vrots/(4.74d0*(d_p*1.0e5))
            vb_p=vz_p*VUKMS/(4.74d0*(d_p*1.0e5))

c *** Find il,ib,id (position in the extinction map grid)
c     for each particle in the input ***
            il=int(l_p/dl)
            ib=int((b_p+(M_PI/2.0d0))/db)
            il3=int(l_p/dl3)
            ib3=int((b_p+(M_PI/2.0d0))/db3)
            id3=int((dlog10(d_p*100.0d0)+2.0d0)/dd3)
            if (float(id3).gt.99.0) then
               id3=99
            else if (float(id3).lt.0.0) then
               id3=0
            endif
c *** Find weighting (extinction interpolated between contiguous
c     extinction bins to the real position of particle in 
c     extinction Schegel or 2D bin) ***
            wl2=(l_p-leS(il))/dl
            wl1=1-wl2
            wb2=(b_p-beS(ib))/db
            wb1=1-wb2
c *** Find position in ScEx array ***
            if (Schlegel.eq.1) then
               ie1=(ib)*MNSC+il
               ie2=(ib+1)*MNSC+il
               ie3=(ib)*MNSC+il+1
               ie4=(ib+1)*MNSC+il+1
            else
               ie1=(ib)*MN2D+il
               ie2=(ib+1)*MN2D+il
               ie3=(ib)*MN2D+il+1
               ie4=(ib+1)*MN2D+il+1
            endif
c *** Schlegel (or 2D) extinction value ***
            exS=datS(ie1)*wl1*wb1+datS(ie2)*wl1*wb2
     &         +datS(ie3)*wl2*wb1+datS(ie4)*wl2*wb2
c *** Weighting for 3D (extinction interpolated between contiguous
c     extinction bins to the real position of particle in
c     3D extinction bin) ***
            wl2=(l_p-le3(il3))/dl3
            wl1=1-wl2
            wb2=(b_p-be3(ib3))/db3
            wb1=1-wb2
            wd2=(dlog10(d_p*100.0d0)-de3(id3))/dd3
            wd1=1-wd2
c *** Find position in 3DEx array ***
            ie1=(MN3D*MN3D*id3)+(MN3D*ib3)+il3
            ie2=(MN3D*MN3D*(id3+1))+(MN3D*ib3)+il3
            ie3=(MN3D*MN3D*(id3+1))+(MN3D*ib3+1)+il3
            ie4=(MN3D*MN3D*(id3+1))+(MN3D*(ib3+1))+il3+1
            ie5=(MN3D*MN3D*(id3+1))+(MN3d*ib3)+il3+1
            ie6=(MN3D*MN3D*id3)+(MN3D*(ib3+1))+il3
            ie7=(MN3D*MN3D*id3)+(MN3D*(ib3+1))+il3+1
            ie8=(MN3D*MN3D*id3)+(MN3D*ib3)+il3+1
c *** 3D extinction factor ***
            ex3=dat3(ie1)*wl1*wb1*wd1+dat3(ie2)*wl1*wb1*wd2+
     &          dat3(ie3)*wl1*wb2*wd2+dat3(ie4)*wl2*wb2*wd2+
     &          dat3(ie5)*wl2*wb1*wd2+dat3(ie6)*wl1*wb2*wd1+
     &          dat3(ie7)*wl2*wb2*wd1+dat3(ie8)*wl2*wb1*wd1
c *** Final extinciton value ***
            ex=ex3*exS  !Schlegel+98 or 2D + 3D correction

c *** Schlegel Modified as per Sharma et al (2014) ***
c            ex=ex*(0.6d0+0.2d0*(1.0d0-dtanh((ex-0.15d0)/0.1d0)))
            if (q.eq.0) then
c               write(6,*) 'Using correction from Sharma et al. (2014)'
               write(6,*) 'Not using correction from Sharma
     & et al. (2014)'
            endif
         else
            ex=0.0
         endif
         ex31=ex*3.24d0 !Extinction to absorption in V band
         ex148=ex*1.962d0 ! extinction to absorption in I band
c AI/E(B − V) = 1.962 and AV /E(B − V) = 3.24, Savage & Mathis 1979; Schlegel et al. 1998)
c *** End Extinction ***


c *** Assigning isochrones to particle depending on metals
c *** total age bins = 181
c All isochrones has been stored in same array, now select the entire
c isochrone for the proper metals (lbegin to lend)
         MiMax=0.0d0
         do i=1,MNiso-1
            if (z.le.dble(IsoNo(i))) then
               if (dble(IsoNo(i))-z.lt.z-dble(IsoNo(i-1))) then
                  lbegin=jindex(181*(i)+ageindex)
                  lend=jindex(181*(i)+ageindex+1)-1
               else
                  lbegin=jindex(181*(i-1)+ageindex)
                  lend=jindex(181*(i-1)+ageindex+1)-1
               endif
               do l=lbegin,lend
c *** Calc MiMax: maximum stellar mass to take in the isochrone 
c                 depending on the age of the population (sim particle) ***
c WARNING: We are taking only the initial mass of the star, not accounting
c for mass loss in stellar evolution (Mc(l)=mass at present).
                  if (Mi(l).gt.MiMax) then
                     MiMax=Mi(l)
                  endif
               enddo
c *** Calc constants for the IMF (imftype=1 is Salpeter)
               if(imftype.eq.1) then
               Afac=-m_p*0.35d0/(MiMax**(-0.35d0)-Msmin**(-0.35d0))
               endif
c *** Calc Milb that is the minimum star mass limit for stars that will contribute
c to the total magnitud (not interested into fainter stars that
c     will not be seen.
c We take advantage that in the isochrone library star contribution is 
c ordered by mass from low to high.
               Milb=10.0d0**12.0d0 !A very high number just in case
               do l=lbegin,lend
                  mag=VMV(l)+dmod+ex31 !computing aparent magnitude
                  mag2=VMV(l-1)+dmod+ex31
                  if (mag.lt.ChosenLim) then
                     frac=(ChosenLim-mag2)/(mag-mag2) !Interpolation between
                     Milb=(Mi(l)-Mi(l-1))*frac+Mi(l-1)!present and previous to
c the limit, mass.
                     goto 667
                  endif

c If star particle is too far to be observed (The more massive 
c stars when splitting the stellar population have Vmag>limit)
c selected min mass will be higher than selected max mass, then
c no stars are generated.

               enddo
            
c *** Calc number to generate *** (new part)
  667          continue
             if(imftype.eq.1) then
               Nsf=((Afac/(-1.35d0))*(MiMax**(-1.35d0)-Milb
     &              **(-1.35d0)))
               Ns=int(Nsf)
c select randomly if we choose one more or one less star (Nsf may not be
c an integer)
               if (ran2(idum).lt.(Nsf-Ns)) then
                  Ns=Ns+1
               endif

c Control of the total stellar mass generated if Nstars=Ns
               mtn=(Afac/(-0.35d0))
     &              *(MiMax**(-0.35d0)-Milb**(-0.35d0))
c Control for Ns=0 (if error on m_p or particle too far, too faint)
               if (Ns.lt.0.0d0) then
                  Ns=0
                  mtn=0.0d0
               endif
              endif
               mtotalneed=mtotalneed+mtn
               if (Ns.eq.0) then
                  goto 601
               endif
c *** Generate stars ***
               do m=0,Ns-1
                  R=ran2(idum)
c *** Calculate stellar mass ***
                  mgen=(R*(MiMax**(-1.35))+(1-R)*(Milb**(-1.35)))
     &                 **(1/(-1.35))
                  do l=lbegin,lend
c Again taking advantage of that mass is ordered from less to more 
c massive stars in the Padova files.
                     if (mgen.lt.Mi(l)) then
c *** Interpolate for V and V-I ***
                        frac=(mgen-Mi(l-1))/(Mi(l)-Mi(l-1))
                        avgen=frac*(VMV(l)-VMV(l-1))+VMV(l-1)
                        aigen=frac*(VMI(l)-VMI(l-1))+VMI(l-1)
c *** Add stellar extinction (and absolute -> apparent) ***
                        vgen=avgen+dmod+ex31
c *** If satisfies max Vmag criteria ***
                        if (vgen.lt.ChosenLim) then
                           igen=aigen+dmod+ex148
                           VIgen=vgen-igen
c *** Write stellar output ***
c All stars are at the same position as the mother particle
                           x_p=x_p*LUKPC
                           y_p=y_p*LUKPC
                           z_p=z_p*LUKPC
                       call carte_to_equatorial(x_p,y_p,z_p,alpha,delta)
                           a(1)=alpha
                           a(2)=delta 
                           a(3)=1.0d0/(d_p*100000.0d0)
                           uh=vx_p*VUKMS-usun
                           vh=vy_p*VUKMS-vlsr-vsun
                           wh=vz_p*VUKMS-wsun
                           call UVWH_to_equatorials(uh,vh,wh,a)
                           jflag=-1 
                           do j=1,4
                              ap(j)=0.0d0
                           enddo
                           do j=3,5
                              a(j)=a(j)*1000.0d0
                           enddo
c *** Add Gaia errors ***
                           if (error.eq.1) then
c *** mission length in months, 60=5 years ***
                              month=60.0d0
c *** Calibration astrometric factor:
c *** CAfactor=1.0d0 means no calibration errors added ***
c *** CAfactor=1.2d0 gives 20% calibration error (on the safe side) ***
                              CAfactor=1.2d0
                              call Gaia_errors(month,CAfactor,jflag,Vgen
     &                             ,VIgen,a,ao,ae,p,po,pe,ap,apo,ape
     &                             ,idum)  
                           endif  
                           do j=3,5
                              a(j)=a(j)/1000.0d0
                              ao(j)=ao(j)/1000.0d0
                           enddo
                           do j=1,5
                              ae(j)=ae(j)/1000.0d0
                           enddo

                          call equatorials_TO_UVWH(ao,Uo,Vo,Wo)
                 call equatorial_to_carte(1.d3/ao(3),ao(1),ao(2),
     +xo_p,yo_p,zo_p)
                           vxo_p=(Uo+usun)/VUKMS
                           vyo_p=(Vo+vsun+VLSR)/VUKMS
                           vzo_p=(Wo+wsun)/VUKMS
                           rop=dsqrt((xo_p+xsun)**2+(yo_p+ysun)**2)
                           od_p=1.0d0/(ao(3)*1.0e3)
                           vro_p=((vxo_p*(xo_p+xsun)+vyo_p*(yo_p+ysun)
     &                          +vzo_p*(zo_p+zsun))/rop)*VUKMS
                          vroto_p=(-1.0d0*((xo_p+xsun)*vyo_p-(yo_p+ysun)
     &                          *vxo_p)/rop)*VUKMS
                           vrots=(-1.0d0*(xo_p*vyo_p-yo_p*vxo_p)/
     &                          (od_p/100.0d0))*VUKMS
                           rop=dsqrt((xo_p+xsun)**2+(yo_p+ysun)**2
     &                          +(zo_p+zsun)**2)
                           vlo_p=vrots/(4.74d0*(od_p*1000.0d0))
                           vbo_p=vzo_p*VUKMS/(4.74d0*(od_p*1000.0d0))
                           x_p=x_p/LUKPC
                           y_p=y_p/LUKPC
                           z_p=z_p/LUKPC
                           d_p=d_p*LUKPC
                           GRVS=Vgen-0.0119d0-1.2092d0*VIgen+0.0188d0
     &                          *(VIgen**2)+0.0005d0*(VIgen**3)
                           makmso=ao(4)*4.74d0*od_p*1000.0d0
                           mdkmso=ao(5)*4.74d0*od_p*1000.0d0
                           makms=a(4)*4.74d0*d_p*1000.0d0
                           mdkms=a(5)*4.74d0*d_p*1000.0d0
c **** HER YOU CAN ADD CUSTOM REQUIREMENTS FOR OUTPUTS ****
c                           if ((od_p.gt.0.0d0).and.
c     &                          (od_p.lt.40.0d0)) then
c                           if (GRVS.lt.16.5d0) then
                              if (Binary.eq.0) then
c *** Ascii ***
                                write(10,102) (a(j),j=1,6),(ao(j),j=1,6)
     &                             ,(ae(j),j=1,6),(p(j),j=1,2)
     &                             ,(po(j),j=1,2),(pe(j),j=1,2)
     &                             ,vgen,d_p,x_p,y_p,z_p,vx_p,vy_p,vz_p
     &                             ,mgen,GRVS,VIgen,ex31,avgen
     &                             ,avgen-aigen,z,PAge,od_p,xo_p,yo_p
     &                             ,zo_p,vxo_p,vyo_p,vzo_p,vro_p,vroto_p
     &                             ,rp,vr_p,vrot_p,rop,makms,mdkms
     &                             ,makmso,mdkmso,vl_p,vb_p,vlo_p,vbo_p
     &                             ,l_p,b_p
c *** Binary ***
                             else
                                write(10) (a(j),j=1,6),(ao(j),j=1,6)
     &                             ,(ae(j),j=1,6),(p(j),j=1,2)
     &                             ,(po(j),j=1,2),(pe(j),j=1,2)
     &                             ,vgen,d_p,x_p,y_p,z_p,vx_p,vy_p,vz_p
     &                             ,mgen,GRVS,VIgen,ex31,avgen
     &                             ,avgen-aigen,z,PAge,od_p,xo_p,yo_p
     &                             ,zo_p,vxo_p,vyo_p,vzo_p,vro_p,vroto_p
     &                             ,rp,vr_p,vrot_p,rop,makms,mdkms
     &                             ,makmso,mdkmso,vl_p,vb_p,vlo_p,vbo_p
     &                             ,l_p,b_p
                              endif
                              GenSel=GenSel+1
c                           endif
c                           endif
c *** END CUSTOM SELECTIONS ***
                           d_p=d_p/LUKPC
                        endif
                        goto 602
                     endif
                  enddo
  602             sumgen=sumgen+mgen
                  mtemp=mtemp+mgen
                  ngen=ngen+1
               enddo
               goto 601
            endif
         enddo
c *** Percentage done ***
  601    if (int((float(q)/float(nparts))*100.0d0).gt.per2) then
            per2=int((float(q)/float(nparts))*100.0d0)
            if (inputdata.lt.2) then
               write(6,*) per2,'% complete'
            endif
         endif
       endif
      enddo 
  102 format(63(1pE13.5))
      close(11)
c *** Finish ***
  603 write(6,*) '100 % complete'
      write(6,*) 'You generated',sumgen,' solar masses of stars!'
      write(6,*) 'You needed',mtotalneed,' solar masses of stars!'
      if (sumgen.gt.mtotalneed) then
         write(6,*) 'You have',sumgen-mtotalneed,' excess'
      else
         write(6,*) 'You have',mtotalneed-sumgen,' deficit'
      endif
      write(6,*) 'You generated',ngen,' stars',GenSel,' selected'
      close(30)
      close(10)

      END
