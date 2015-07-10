! This program was written to generate Real Gas Property (RGP) files to be used
! with ANSYS CFX (Current version 14.5).  The user is required to set certain
! input variables which generates an RGP file that follows the ANSYS guidelines
! see http://www.kxcad.net/ansys/ANSYS_CFX/help/help/Modelling/bk05ch10s07s01.html

!Current Features:
! Works for CO2
! Output file to ANSYS CFX version 14.5 specifications

PROGRAM RGP

      !implicit double precision (a-h,o-z)
      !implicit integer (i-k,m,n)
      !implicit character*30 (l)
      implicit NONE
      integer ncmax
      parameter (ncmax=20)   !max number of components in mixture
      double precision, dimension (ncmax) :: x(ncmax),xliq(ncmax),xvap(ncmax)
      double precision, allocatable :: t(:),p(:),h(:,:),w(:,:),v(:,:),cv(:,:),cp(:,:),s(:,:),dPdv(:,:),eta(:,:),tcx(:,:), &
      &tsat(:),hsat(:),wsat(:),vsat(:),cvsat(:),cpsat(:),ssat(:),dPdvsat(:),etasat(:),tcxsat(:), &
      &tsat2(:),hsatl(:),wsatl(:),rhosatl(:),cvsatl(:),cpsatl(:),ssatl(:),dPdrhosatl(:),etasatl(:),tcxsatl(:), &
      &psat2(:),hsatv(:),wsatv(:),rhosatv(:),cvsatv(:),cpsatv(:),ssatv(:),dPdrhosatv(:),etasatv(:),tcxsatv(:), &
      &blank1(:),blank2(:),arr(:)
      double precision :: acf,conp,conv,conpl,convl,conpv,convv,cpcrit,cvcrit,density,dpdrho,q,rgas,rhocrit,rhol,rhov, &
      &dc,dip,dl,dm,dmax,dpdrhol,dpdrhov,dv,e,enthalpy,enthalpyl,enthalpyv,pstep,pstepsat,scrit,speedcrit, &
      &entropy,entropyl,entropyv,hcrit,htyp,pc,pm,pmax,pmaxsat,pminsat,psat,preslower,upperPres,speedOsound, &
      &speedOsoundl,speedOsoundv,supercool,tc,templower,uppertemp,TminSat,TmaxSat,thermalconductivity,thermalconductivityl, &
      &thermalconductivityv,thermcrit,tliq,tm,tmax,tmin,tnbp,tpp,ttp,tstep,tvap,visccrit,viscosity,viscosityl,viscosityv,wm, &
      &y,zc,dbl
      integer :: i,j,ierr,kq,np,nt,nsat,npsat,ilng
      character :: hrf*3, herr*255, lfluid*30, loutfile*30
      character :: lFileName*255, hnam*12, hn80*80, hstr*255
      character :: hcasn*12
      character(LEN=255) :: hf(ncmax)
      character*7 hfmix
      integer :: prevStatus
      integer :: currStatus
      integer :: tenPercent
      integer :: AllocStatus
      i=1
      !PROMPT FOR OUTPUT FILE
      WRITE(*,*) '                       Carleton University'
      WRITE(*,*) '        Department of Mechanical and Aerospace Engineering' 
      WRITE(*,*)
      WRITE(*,*) '               88888888ba    ,ad8888ba,  88888888ba   '
      WRITE(*,*) '               88      "8b  d8"      "8b 88      "8b  '
      WRITE(*,*) '               88      .8P d8            88       8P  '
      WRITE(*,*) '               88aaaaaa8P  88            88aaaaaa8P   '
      WRITE(*,*) '               88""""88    88      88888 88""""""     '
      WRITE(*,*) '               88    `8b   Y8,       `88 88           '
      WRITE(*,*) '               88     `8b   Y8a.    .a88 88           '
      WRITE(*,*) '               88      `8b   `"Y88888P"  88           '
      WRITE(*,*)
      WRITE(*,*) '                           FILE GENERATOR'
      WRITE(*,*)
      WRITE(*,'(A)') 'Output Data File Name (i.e. filename.txt):'
      READ (*,*) lOutfile
      OPEN(UNIT=10,FILE=lOutFile,STATUS='NEW')
      lfluid='SCO2'
      WRITE(*,'(A)') 'Fluid file name (i.e. filename.fld):'
      READ (*,*) lFileName
      hf(1)=lFileName
      hfmix='hmx.bnc'
      hrf='DEF'
      CALL SETUP (i,hf,hfmix,hrf,ierr,herr)
      IF (ierr.ne.0) write (*,*) herr
      CALL INFO (1,wm,ttp,tnbp,tc,pc,dc,zc,acf,dip,rgas)
      CALL NAME (1,hnam,hn80,hcasn)
      Call PASSCMN ('ptpn',0,1,0, hstr,ilng,tpp,arr,ierr,herr)
      dc=wm*dc!Convert density to kg/m^3
      Dmax=wm*Dmax
      !tpp = 517.95 !kPa *********** ADDED FROM CO2.FLD *******
      WRITE(*,*)
      WRITE(*,*) 'hnam',hnam 
      WRITE(*,*) 'hn80',hn80
      WRITE(*,*) 'Molecular Weight:                 ',wm,'g/mol'
!     WRITE(*,*) 'Accentric Factor [-]: ',acf SETUP error
!     WRITE(*,*) 'Dipole Moment [debye]: ',dip
      WRITE(*,*) 'Triple Point Temperature:         ',ttp,'K'
      WRITE(*,*) 'Triple Point Pressure:            ',tpp,'kPa'
      WRITE(*,*) 'Normal Boiling Point:             ',tnbp,'K'
      
      WRITE(*,*) 'Critical Temperature:             ',tc,'K'
      WRITE(*,*) 'Critical Pressure:                ',pc,'kPa'
      WRITE(*,*) 'Critical Density:                 ',dc,'kg/m^3'
      WRITE(*,*) 'Compressibility at Critical Point: ',zc
      WRITE(*,*) 'Gas Constant:                     ',rgas,'J/mol-K'
      CALL MAXT (x,tm,pm,Dm,ierr,herr)
      CALL LIMITS(htyp,x,tmin,tmax,Dmax,pmax)
      WRITE(*,*)
      WRITE(*,*)'Limits for this fluid:'
      WRITE(*,*)
      WRITE(*,*)'Min Temperature:                  ',tmin,'K'
      WRITE(*,*)'Max Temperature:                  ',tmax,'K'
      WRITE(*,*)'Max Saturation temperature:        ',tm,'K'
      WRITE(*,*)'Max Density:                      ',Dmax,'kg/m^3'
      WRITE(*,*)'Max Pressure:                     ',pmax,'kPa'
      
      ! PROMPT FOR BOUNDS AND STEP SIZE
      
      WRITE(*,*)'--------------------------------------------'
      WRITE(*,*)
      WRITE(*,*)'Enter - Lower bound for Temperature [K]:'
      READ(*,*)TempLower
      WRITE(*,*)'Enter - Upper bound for Temperature [K]:'
      READ(*,*)UpperTemp
      WRITE(*,*)'Enter - Number of Steps for Temperature:'
      READ(*,*)nt
      WRITE(*,*)'Enter - Lower bound for Pressure [Pa]:'
      READ(*,*)PresLower
      WRITE(*,*)'Enter - Upper bound for Pressure [Pa]:'
      READ(*,*)UpperPres
      WRITE(*,*)'Enter - Number of Steps for Pressure:'
      READ(*,*)np
      WRITE(*,*)'Enter - Minimum Temperature for Saturation Tables:'
      READ(*,*)TminSat
      WRITE(*,*)'Enter - Maximum Temperature for Saturation Tables:'
      READ(*,*)TmaxSat
      WRITE(*,*)'Enter - Number of Steps for Saturation Tables:'
      READ(*,*)npsat
      
      SuperCool=0
      
      WRITE(*,*)
      WRITE(*,*)'Generating RGP File Please Wait...'
      WRITE(*,*)
      
      ! SETUP OUTPUT FILE 
      WRITE(10,901)' $$$$HEADER'
      WRITE(10,901)' $$$',hnam!,lfluid !character*8 (key into the RGP file see NOTE 3)  if this program is used with a fluid other than CO2, this line will need to be changed, along with the corresponding line in the $$$DATA header
      WRITE(10,*) 1!enter an integer (this line is ignored in ANSYS CFX)
      WRITE(10,901)'$$PARAM'! See NOTE 4
      WRITE(10,*) 26!enter an integer(number of parameters)
      WRITE(10,901)'DESCRIPTION'
      WRITE(10,901)'Carbon Dioxide from NIST'!character*50(description of the material)
      WRITE(10,901)'NAME'
      WRITE(10,901)hnam!lfluid!character*8(material name, same as $$$<component>)
      WRITE(10,901)'INDEX'
      WRITE(10,901)'SCO2'!lfluid!character*50 (index into clients RGDB program)
      WRITE(10,901)'DATABASE'                                      
      WRITE(10,901)'NIST REAL GAS PROPERTY DATABASE'                                      
      WRITE(10,901)'MODEL'                                      
      WRITE(10,*)3!integer(level of property info available 1,2,or 3)
      WRITE(10,901)'UNITS'                                      
      WRITE(10,*)1!integer(unit system of 1 ,2,3,4 or 5)
      WRITE(10,901)'PMIN_SUPERHEAT'!	at Triple Point Pressure?
      WRITE(10,900)PresLower!real (Pa)
      WRITE(10,901)'PMAX_SUPERHEAT'!  at Critical Pressure?
      WRITE(10,900)UpperPres!real (Pa)
      WRITE(10,901)'TMIN_SUPERHEAT'! at  Triple Point Temperature?
      WRITE(10,900)TempLower!real (K)
      WRITE(10,901)'TMAX_SUPERHEAT'! at Critical Temperature?
      WRITE(10,900)UpperTemp!real (K)
      WRITE(10,901)'SUPERCOOLING'
      WRITE(10,900)SuperCool!real(supercooling level in superheat tables)(Optional)
      WRITE(10,901)'P_CRITICAL'
      WRITE(10,900)Pc*1000!real
      WRITE(10,901)'P_TRIPLE'!        See NOTE 7
      WRITE(10,900)tpp*1000!real
      WRITE(10,901)'T_CRITICAL'
      WRITE(10,900)Tc!real
      WRITE(10,901)'T_TRIPLE'
      WRITE(10,900)Ttp!real
      WRITE(10,901)'GAS_CONSTANT'
      WRITE(10,900)rgas*1000/wm!real
      
      WRITE(10,901)'TABLE_1'!See NOTE 8
      WRITE(10,*)nt,np!integer
      WRITE(10,901)'TABLE_2'
      WRITE(10,*)nt,np!integer
      WRITE(10,901)'TABLE_3'
      WRITE(10,*)nt,np!integer
      WRITE(10,901)'TABLE_4'
      WRITE(10,*)nt,np!integer
      WRITE(10,901)'TABLE_5'
      WRITE(10,*)nt,np!integer
      WRITE(10,901)'TABLE_6'
      WRITE(10,*)nt,np!integer
      WRITE(10,901)'TABLE_7'
      WRITE(10,*)nt,np!integer
      WRITE(10,901)'TABLE_8'!		See NOTE 9
      WRITE(10,*)nt,np!integer
      WRITE(10,901)'TABLE_9'
      WRITE(10,*)nt,np!integer
      WRITE(10,901)'SAT_TABLE'
      WRITE(10,*)npsat,4,9 !integer
      
      WRITE(10,901) '$$$$DATA'
      WRITE(10,901) '$$$SCO2'!,lfluid !character*8 (key into the RGP file see NOTE 3)
      WRITE(10,*) 1!enter an integer (this line is ignored in ANSYS CFX)
      WRITE(10,901)'$$PARAM'! See NOTE 4
      WRITE(10,*) 26!enter an integer(number of parameters)
      WRITE(10,901)'DESCRIPTION'
      WRITE(10,901)'Carbon Dioxide from NIST'!character*50(description of the material)
      WRITE(10,901)'NAME'
      WRITE(10,901)'SCO2'!lfluid!character*8(material name, same as $$$<component>)
      WRITE(10,901)'INDEX'
      WRITE(10,901)'SCO2'!lfluid!character*50 (index into clients RGDB program)
      WRITE(10,901)'DATABASE'                                      
      WRITE(10,901)'NIST REAL GAS PROPERTY DATABASE'                                      
      WRITE(10,901)'MODEL'                                      
      WRITE(10,*)3!integer(level of property info available 1,2,or 3)
      WRITE(10,901)'UNITS'                                      
      WRITE(10,*)1!integer(unit system of 1 ,2,3,4 or 5)
      WRITE(10,901)'PMIN_SUPERHEAT'!	at Triple Point Pressure?
      WRITE(10,900)PresLower!real (Pa)
      WRITE(10,901)'PMAX_SUPERHEAT'!  at Critical Pressure?
      WRITE(10,900)UpperPres!real (Pa)
      WRITE(10,901)'TMIN_SUPERHEAT'! at  Triple Point Temperature?
      WRITE(10,900)TempLower!real (K)
      WRITE(10,901)'TMAX_SUPERHEAT'! at Critical Temperature?
      WRITE(10,900)UpperTemp!real (K)
      WRITE(10,901)'SUPERCOOLING'
      WRITE(10,900)SuperCool!real(supercooling level in superheat tables)(Optional)
      WRITE(10,901)'P_CRITICAL'
      WRITE(10,900)Pc*1000!real
      WRITE(10,901)'P_TRIPLE'!        See NOTE 7
      WRITE(10,900)tpp*1000!real
      WRITE(10,901)'T_CRITICAL'
      WRITE(10,900)Tc!real
      WRITE(10,901)'T_TRIPLE'
      WRITE(10,900)Ttp!real
      WRITE(10,901)'GAS_CONSTANT'
      WRITE(10,900)rgas*1000/wm!real
      
      WRITE(10,901)'TABLE_1'!See NOTE 8
      WRITE(10,*)nt,np!integer
      WRITE(10,901)'TABLE_2'
      WRITE(10,*)nt,np!integer
      WRITE(10,901)'TABLE_3'
      WRITE(10,*)nt,np!integer
      WRITE(10,901)'TABLE_4'
      WRITE(10,*)nt,np!integer
      WRITE(10,901)'TABLE_5'
      WRITE(10,*)nt,np!integer
      WRITE(10,901)'TABLE_6'
      WRITE(10,*)nt,np!integer
      WRITE(10,901)'TABLE_7'
      WRITE(10,*)nt,np!integer
      WRITE(10,901)'TABLE_8'!		See NOTE 9
      WRITE(10,*)nt,np!integer
      WRITE(10,901)'TABLE_9'
      WRITE(10,*)nt,np!integer
      WRITE(10,901)'SAT_TABLE'
      WRITE(10,*)npsat,4,9 !integer
      
      WRITE(10,901)'$$SUPER_TABLE'
      WRITE(10,*)9 !integer (number of superheat tables, nn = 9)
      tstep = (UpperTemp-TempLower)/REAL(nt-1)
      pstep = (UpperPres-PresLower)/REAL(np-1)

      ! allocate 2D array
      WRITE (*,*) 'Allocating superheat array...'
      ALLOCATE (t(nt), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 001'
      ALLOCATE (p(np), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 002'
      ALLOCATE (h(nt,np), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 003'
      ALLOCATE (w(nt,np), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 004'
      ALLOCATE (v(nt,np), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 005'
      ALLOCATE (cv(nt,np), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 006'
      ALLOCATE (cp(nt,np), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 007'
      ALLOCATE (dPdv(nt,np), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 008'
      ALLOCATE (s(nt,np), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 009'
      ALLOCATE (eta(nt,np), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 010'
      ALLOCATE (tcx(nt,np), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 011'
      
      ! allocate 1D array
      nsat = np
      ALLOCATE (tsat(nsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 012'
      ALLOCATE (hsat(nsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 013'
      ALLOCATE (wsat(nsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 014'
      ALLOCATE (vsat(nsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 015'
      ALLOCATE (cvsat(nsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 016'
      ALLOCATE (cpsat(nsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 017'
      ALLOCATE (dPdvsat(nsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 018'
      ALLOCATE (ssat(nsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 019'
      ALLOCATE (etasat(nsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 020'
      ALLOCATE (tcxsat(nsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 021'
      
      CALL TPFLSH (tc,pc,x,rhoCrit,dl,dv,xliq,xvap,q,e,Hcrit,Scrit,cvCrit,cpCrit,SpeedCrit,ierr,herr)
      IF (ierr.ne.0) write (*,*) herr
      CALL DPDD (tc,rhoCrit,x,dPdrhov)
      IF (ierr.ne.0) write (*,*) herr
      CALL TRNPRP (tc,rhoCrit,x,viscCrit,ThermCrit,ierr,herr)
      IF (ierr.ne.0) write (*,*) herr
      
      !Fill superheat arrays with RGP data
      !For 3D table:
      !any point in the liquid area gets filled with saturation properties at that pressure
      !any point below the crit temp but above the crit pressure gets filled with phi(T, P) (ie normally)
      !For 2D tables:
      !for every pressure in the 2D table there is an entry in the 1D table
      !for pressures above Pcrit, the entry is of phi(Tmin, P) instead of phi(Tsat, P)
      write (*,*) 'Filling superheat table...'
      prevStatus = 0
      currStatus = 0
      tenPercent = INT(UpperPres-PresLower)/10
      j=1
      p(j) = PresLower
      DO WHILE (j <= np)
          i=1
          t(i) = TempLower
          
          ! show progress every 10%
          currStatus = INT(p(j) - PresLower)/tenPercent
          IF (prevStatus /= currStatus) THEN
              write (*,*) (currStatus*10), '%'
              prevStatus = currStatus
          END IF
          
          DO WHILE (i <= nt)
              IF (p(j) <= pc*1000) THEN
                  ! find the sat temp for this pressure, the if we are above it, work normal
                  ! if we are below the sat temp, then fill 2D table with phi(Tsat, p)
                  PSat= p(j)/1000
                  q=1
                  CALL PQFLSH (PSat,q,x,kq,tvap,rhov,Dl,Dv,xliq,xvap,e,enthalpyv,entropyv,convv,conpv,speedOsoundv,ierr,herr)
                  IF (ierr.ne.0) write (*,*) herr
                  
                  IF (t(i) <= tvap) THEN
                      ! use Tsat properties at the given p
                      h(i,j)  = enthalpyv*1000/wm
                      w(i,j)  = speedOsoundv
                      v(i,j)  = 1/(rhov*wm)
                      cv(i,j) = convv*1000/wm
                      cp(i,j) = conpv*1000/wm
                      s(i,j)  = entropyv*1000/wm
                      CALL DPDD (tvap,rhov,x,dPdrhov)
                      IF (ierr.ne.0) write (*,*) herr
                      dPdv(i,j) = (dPdrhov/wm*1000*(-(rhov*wm)**2))
                      CALL TRNPRP (tvap,rhov,x,viscosityv,thermalconductivityv,ierr,herr)
                      IF (ierr.ne.0) write (*,*) herr
                      eta(i,j) = viscosityv/1000000
                      tcx(i,j) = thermalconductivityv
                  ELSE
                      ! get properties normally
                      CALL TPFLSH (t(i),(p(j)/1000),x,density,dl,dv,xliq,xvap,q,e,enthalpy,entropy,conv,conp,speedOsound,ierr,herr)
                      IF (ierr.ne.0) write (*,*) herr
                      h(i,j) = enthalpy*1000/wm
                      w(i,j) = speedOsound
                      v(i,j) =  1/(density*wm)
                      cv(i,j) = conv*1000/wm
                      cp(i,j) = conp*1000/wm
                      s(i,j) = entropy*1000/wm
                      CALL DPDD (t(i),density,x,dPdrho)
                      IF (ierr.ne.0) write (*,*) herr
                      dPdv(i,j) = (dPdrho/wm*1000*(-(density*wm)**2))
                      CALL TRNPRP (t(i),density,x,viscosity,thermalconductivity,ierr,herr)
                      IF (ierr.ne.0) write (*,*) herr
                      eta(i,j) = viscosity/1000000
                      tcx(i,j) = thermalconductivity
                  END IF
              ELSE
                  ! P > Pcrit
                  
                  ! run normally
                  CALL TPFLSH (t(i),(p(j)/1000),x,density,dl,dv,xliq,xvap,q,e,enthalpy,entropy,conv,conp,speedOsound,ierr,herr)
                  IF (ierr.ne.0) write (*,*) herr
                  h(i,j) = enthalpy*1000/wm
                  w(i,j) = speedOsound
                  v(i,j) =  1/(density*wm)
                  cv(i,j) = conv*1000/wm
                  cp(i,j) = conp*1000/wm
                  s(i,j) = entropy*1000/wm
                  CALL DPDD (t(i),density,x,dPdrho)
                  IF (ierr.ne.0) write (*,*) herr
                  dPdv(i,j) = (dPdrho/wm*1000*(-(density*wm)**2))
                  CALL TRNPRP (t(i),density,x,viscosity,thermalconductivity,ierr,herr)
                  IF (ierr.ne.0) write (*,*) herr
                  eta(i,j) = viscosity/1000000
                  tcx(i,j) = thermalconductivity
              END IF
              
              i = i + 1
              IF (i <= nt) THEN
                  t(i) = t(i-1)+tstep
              END IF
          END DO
          
          ! Fill 1D Table (Sat line)
          IF (p(j)<(Pc*1000))THEN
              ! fill with phi(Tsat, P)
              PSat= p(j)/1000
              q=1
              CALL PQFLSH (PSat,q,x,kq,tvap,rhov,Dl,Dv,xliq,xvap,e,enthalpyv,entropyv,convv,conpv,speedOsoundv,ierr,herr)
              IF (ierr.ne.0) write (*,*) herr
              tsat(j) = tvap
              hsat(j) = enthalpyv*1000/wm
              wsat(j) = speedOsoundv
              vsat(j) = 1/(rhov*wm)
              cvsat(j) = convv*1000/wm
              cpsat(j) = conpv*1000/wm
              ssat(j) = entropyv*1000/wm
              CALL DPDD (tvap,rhov,x,dPdrhov)
              IF (ierr.ne.0) write (*,*) herr
              dPdvsat(j) = (dPdrhov/wm*1000*(-(rhov*wm)**2))
              CALL TRNPRP (tvap,rhov,x,viscosityv,thermalconductivityv,ierr,herr)
              IF (ierr.ne.0) write (*,*) herr
              etasat(j) = viscosityv/1000000
              tcxsat(j) = thermalconductivityv
          ELSE
              ! fill with phi(Tmin, P)
              PSat= p(j)/1000
              q=1
              CALL TPFLSH (TempLower,(p(j)/1000),x,density,dl,dv,xliq,xvap,q,e,enthalpy,entropy,conv,conp,speedOsound,ierr,herr)
              IF (ierr.ne.0) write (*,*) herr
              tsat(j) = TempLower
              hsat(j) = enthalpy*1000/wm
              wsat(j) = speedOsound
              vsat(j) = 1/(density*wm)
              cvsat(j) = conv*1000/wm
              cpsat(j) = conp*1000/wm
              ssat(j) = entropy*1000/wm
              CALL DPDD (TempLower,density,x,dPdrho)
              IF (ierr.ne.0) write (*,*) herr
              dPdvsat(j) = (dPdrho/wm*1000*(-(density*wm)**2))
              CALL TRNPRP (TempLower,density,x,viscosity,thermalconductivity,ierr,herr)
              IF (ierr.ne.0) write (*,*) herr
              etasat(j) = viscosity/1000000
              tcxsat(j) = thermalconductivity
          END IF
          
          j = j+1
          IF (j <= np) THEN
              p(j) = p(j-1)+pstep
          END IF
      END DO
      
      q=1
      CALL TQFLSH (TminSat,q,x,kq,PminSat,density,Dl,Dv,x,y,e,enthalpy,entropy,conv,conp,speedOsound,ierr,herr)
      IF (ierr.ne.0) write (*,*) herr
      CALL TQFLSH (TmaxSat,q,x,kq,PmaxSat,density,Dl,Dv,x,y,e,enthalpy,entropy,conv,conp,speedOsound,ierr,herr)
      IF (ierr.ne.0) write (*,*) herr
      pstepsat = (Pmaxsat-Pminsat)/REAL(npsat-1)
      
      WRITE (*,*) 'Allocating saturation table...'
      ALLOCATE (tsat2(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 022'
      ALLOCATE (psat2(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 023'
      ALLOCATE (blank1(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 024'
      ALLOCATE (blank2(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 025'
      
      ALLOCATE (hsatl(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 026'
      ALLOCATE (wsatl(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 027'
      ALLOCATE (rhosatl(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 028'
      ALLOCATE (cvsatl(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 029'
      ALLOCATE (cpsatl(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 030'
      ALLOCATE (dPdrhosatl(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 031'
      ALLOCATE (ssatl(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 032'
      ALLOCATE (etasatl(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 033'
      ALLOCATE (tcxsatl(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 034'
      
      ALLOCATE (hsatv(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 035'
      ALLOCATE (wsatv(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 036'
      ALLOCATE (rhosatv(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 037'
      ALLOCATE (cvsatv(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 038'
      ALLOCATE (cpsatv(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 039'
      ALLOCATE (dPdrhosatv(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 040'
      ALLOCATE (ssatv(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 041'
      ALLOCATE (etasatv(npsat), STAT = AllocStatus)
      IF (AllocStatus /= 0) write (*,*)'out of memory 042'
      ALLOCATE (tcxsatv(npsat), STAT = AllocStatus) 
      IF (AllocStatus /= 0) write (*,*)'out of memory 043'
      
      ! fill saturation table 
      WRITE (*,*) 'Filling saturation table...'
      j=1
      psat2(j) = PminSat
      DO WHILE (j <= npsat)
          q=1
          CALL PQFLSH (psat2(j),q,x,kq,tvap,rhov,Dl,Dv,xliq,xvap,e,enthalpyv,entropyv,convv,conpv,speedOsoundv,ierr,herr)
          IF (ierr.ne.0) write (*,*) herr
          tsat2(j) = tvap
          blank1(j)=0
          blank2(j)=0
          hsatv(j) = enthalpyv*1000/wm
          wsatv(j) = speedOsoundv
          rhosatv(j) = (rhov*wm)
          cvsatv(j) = convv*1000/wm
          cpsatv(j) = conpv*1000/wm
          ssatv(j) = entropyv*1000/wm
          CALL DPDD (tvap,rhov,x,dPdrhov)
          IF (ierr.ne.0) write (*,*) herr
          dPdrhosatv(j) = wm/(dPdrhov*1000)
          CALL TRNPRP (tsat2(j),rhov,x,viscosityv,thermalconductivityv,ierr,herr)
          IF (ierr.ne.0) write (*,*) herr
          etasatv(j) = viscosityv/1000000
          tcxsatv(j) = thermalconductivityv
          
          q=0
          CALL PQFLSH (psat2(j),q,x,kq,tliq,rhol,Dl,Dv,xliq,xvap,e,enthalpyl,entropyl,convl,conpl,speedOsoundl,ierr,herr)
          IF (ierr.ne.0) write (*,*) herr
          hsatl(j) = enthalpyl*1000/wm
          wsatl(j) = speedOsoundl
          rhosatl(j) = (rhol*wm)
          cvsatl(j) = convl*1000/wm
          cpsatl(j) = conpl*1000/wm
          ssatl(j) = entropyl*1000/wm
          CALL DPDD (tliq,rhol,x,dPdrhol)
          IF (ierr.ne.0) write (*,*) herr
          dPdrhosatl(j) = wm/(dPdrhol*1000)
          CALL TRNPRP (tsat2(j),rhol,x,viscosityl,thermalconductivityl,ierr,herr)
          IF (ierr.ne.0) write (*,*) herr
          etasatl(j) = viscosityl/1000000
          tcxsatl(j) = thermalconductivityl
          
          j = j+1
          IF (j <= npsat) THEN
              psat2(j) = psat2(j-1) + pstepsat
          END IF
      END DO
      
      WRITE (*,*) 'Writing data to file...'
      !Print data to the text file in the correct format
      WRITE (*,*) 'Table 1'
      WRITE(10,901)'$TABLE_1'
      WRITE(10,*)nt,np!integer(size of superheat arrays, See NOTE 10)
      WRITE(10,900) (t(i),i=1,nt) !|all following real
      WRITE(10,900) (p(j),j=1,np) !|
      WRITE(10,900) ((h(i,j),i=1,nt),j=1,np) !|
      WRITE(10,900) (tsat(i),i=1,nsat) !|
      WRITE(10,900) (hsat(i),i=1,nsat) !|
      
      WRITE (*,*) 'Table 2'
      WRITE(10,901)'$TABLE_2'
      WRITE(10,*)nt,np!integer(size of superheat arrays, See NOTE 10)
      WRITE(10,900) (t(i),i=1,nt) !|all following real
      WRITE(10,900) (p(j),j=1,np) !|
      WRITE(10,900) ((w(i,j),i=1,nt),j=1,np) !|
      WRITE(10,900) (tsat(i),i=1,nsat) !|
      WRITE(10,900) (wsat(i),i=1,nsat) !|
      
      WRITE (*,*) 'Table 3'
      WRITE(10,901)'$TABLE_3'
      WRITE(10,*)nt,np!integer(size of superheat arrays, See NOTE 10)
      WRITE(10,900) (t(i),i=1,nt) !|all following real
      WRITE(10,900) (p(j),j=1,np) !|
      WRITE(10,900) ((v(i,j),i=1,nt),j=1,np) !|
      WRITE(10,900) (tsat(i),i=1,nsat) !|
      WRITE(10,900) (vsat(i),i=1,nsat) !|
      
      WRITE (*,*) 'Table 4'
      WRITE(10,901)'$TABLE_4'
      WRITE(10,*)nt,np!integer(size of superheat arrays, See NOTE 10)
      WRITE(10,900) (t(i),i=1,nt) !|all following real
      WRITE(10,900) (p(j),j=1,np) !|
      WRITE(10,900) ((cv(i,j),i=1,nt),j=1,np) !|
      WRITE(10,900) (tsat(i),i=1,nsat) !|
      WRITE(10,900) (cvsat(i),i=1,nsat) !|
      
      WRITE (*,*) 'Table 5'
      WRITE(10,901)'$TABLE_5'
      WRITE(10,*)nt,np!integer(size of superheat arrays, See NOTE 10)
      WRITE(10,900) (t(i),i=1,nt) !|all following real
      WRITE(10,900) (p(j),j=1,np) !|
      WRITE(10,900) ((cp(i,j),i=1,nt),j=1,np) !|
      WRITE(10,900) (tsat(i),i=1,nsat) !|
      WRITE(10,900) (cpsat(i),i=1,nsat) !|
      
      WRITE (*,*) 'Table 6'
      WRITE(10,901)'$TABLE_6'
      WRITE(10,*)nt,np!integer(size of superheat arrays, See NOTE 10)
      WRITE(10,900) (t(i),i=1,nt) !|all following real
      WRITE(10,900) (p(j),j=1,np) !|
      WRITE(10,900) ((dPdv(i,j),i=1,nt),j=1,np) !|
      WRITE(10,900) (tsat(i),i=1,nsat) !|
      WRITE(10,900) (dPdvsat(i),i=1,nsat) !|
      
      WRITE (*,*) 'Table 7'
      WRITE(10,901)'$TABLE_7'
      WRITE(10,*)nt,np!integer(size of superheat arrays, See NOTE 10)
      WRITE(10,900) (t(i),i=1,nt) !|all following real
      WRITE(10,900) (p(j),j=1,np) !|
      WRITE(10,900) ((s(i,j),i=1,nt),j=1,np) !|
      WRITE(10,900) (tsat(i),i=1,nsat) !|
      WRITE(10,900) (ssat(i),i=1,nsat) !|
      
      WRITE (*,*) 'Table 8'
      WRITE(10,901)'$TABLE_8'
      WRITE(10,*)nt,np!integer(size of superheat arrays, See NOTE 10)
      WRITE(10,900) (t(i),i=1,nt) !|all following real
      WRITE(10,900) (p(i),i=1,np) !|
      WRITE(10,900) ((eta(i,j),i=1,nt),j=1,np) !|
      WRITE(10,900) (tsat(i),i=1,nsat) !|
      WRITE(10,900) (etasat(i),i=1,nsat) !|
      
      WRITE (*,*) 'Table 9'
      WRITE(10,901)'$TABLE_9'
      WRITE(10,*)nt,np!integer(size of superheat arrays, See NOTE 10)
      WRITE(10,900) (t(i),i=1,nt) !|all following real
      WRITE(10,900) (p(j),j=1,np) !|
      WRITE(10,900) ((tcx(i,j),i=1,nt),j=1,np) !|
      WRITE(10,900) (tsat(i),i=1,nsat) !|
      WRITE(10,900) (tcxsat(i),i=1,nsat) !|
      
      WRITE (*,*) 'Sat Table'
      WRITE(10,901)'$$SAT_TABLE'
      WRITE(10,*)npsat,4,9 !(size of saturation arrays, see NOTE 11)
      WRITE(10,900) (psat2(i)*1000,i=1,npsat)
      WRITE(10,900) (tsat2(i),i=1,npsat)
      WRITE(10,900) (blank1(i),i=1,npsat)
      WRITE(10,900) (blank2(i),i=1,npsat)
      WRITE(10,900) (hsatl(i),i=1,npsat)
      WRITE(10,900) (cpsatl(i),i=1,npsat)
      WRITE(10,900) (rhosatl(i),i=1,npsat)
      WRITE(10,900) (dPdrhosatl(i),i=1,npsat)!
      WRITE(10,900) (ssatl(i),i=1,npsat) !
      WRITE(10,900) (cvsatl(i),i=1,npsat)!
      WRITE(10,900) (wsatl(i),i=1,npsat) 
      WRITE(10,900) (etasatl(i),i=1,npsat)
      WRITE(10,900) (tcxsatl(i),i=1,npsat)
      WRITE(10,900) (hsatv(i),i=1,npsat)
      WRITE(10,900) (cpsatv(i),i=1,npsat)!
      WRITE(10,900) (rhosatv(i),i=1,npsat)
      WRITE(10,900) (dPdrhosatv(i),i=1,npsat)!
      WRITE(10,900) (ssatv(i),i=1,npsat)
      WRITE(10,900) (cvsatv(i),i=1,npsat) 
      WRITE(10,900) (wsatv(i),i=1,npsat)
      WRITE(10,900) (etasatv(i),i=1,npsat)
      WRITE(10,900) (tcxsatv(i),i=1,npsat)
      
      CLOSE(10)
      
      WRITE(*,*)
      WRITE(*,*)'RGP File Generated!'  
      READ(*,*)
  900   format (5ES17.7E3)
  901   format (A)
      
      
END PROGRAM RGP
