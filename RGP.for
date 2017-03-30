! This program was written to generate Real Gas Property (RGP) files to be used
! with ANSYS CFX (Current version 14.5).  The user is required to set certain
! input variables which generates an RGP file that follows the ANSYS guidelines
! 
! See Anysys manual section "12.7.1.1. Detailed .rgp File Format"
! this is easily found on Google.

! Current Features:
! Output file to ANSYS CFX version 14.5 specifications
! Input using arguments of interactively
! Any fluid :D

! Notes:
! - Must be compiled with -std=f2008 -ffree-form
! - Other flags can be whatever 
! - Must be linked to refprop objects

module RGP_DATA
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
    character(LEN=7) :: hfmix
    integer :: prevStatus
    integer :: currStatus
    integer :: tenPercent
    integer :: AllocStatus
end module

program RGP
    USE RGP_DATA
    implicit none
    
    write(*,*) '                       Carleton University'
    write(*,*) '        Department of Mechanical and Aerospace Engineering' 
    write(*,*) '                                 &'
    write(*,*) '                Queensland University of Technology'
    write(*,*) '    Department of Chemistry Physics and Mechanical Engineering' 
    write(*,*)
    write(*,*) '               88888888ba    ,ad8888ba,  88888888ba   '
    write(*,*) '               88      "8b  d8"      "8b 88      "8b  '
    write(*,*) '               88      .8P d8            88       8P  '
    write(*,*) '               88aaaaaa8P  88            88aaaaaa8P   '
    write(*,*) '               88""""88    88      88888 88""""""     '
    write(*,*) '               88    `8b   Y8,       `88 88           '
    write(*,*) '               88     `8b   Y8a.    .a88 88           '
    write(*,*) '               88      `8b   `"Y88888P"  88           '
    write(*,*)
    write(*,*) '                           FILE GENERATOR'
    write(*,*)

    call checkArg();
    call getFilesConfig();
    call setupREFPROP();
    call getRGPBounds();

    write(*,*)'Writing header...'
    call writeRGPHeader();

    write(*,*)'Allocating arrays...'
    call allocateArrays();

    write(*,*)'Calculating sat table...'
    call calculateSatTable();

    write(*,*)'Calculating super table...'
    call calculateSuperTable();

    write(*,*)'Writing super table...'
    call writeSuperTable();

    write(*,*)'Writing sat table...'
    call writeSatTable();
    
    write(*,*)'Done! Yay!'
    close(10);

contains
    subroutine checkArg()

        if (COMMAND_ARGUMENT_COUNT() .gt. 1) then
            write(*,*) 'Using arguments for inputs'
            call get_command_argument(1, lOutfile)
            call get_command_argument(2, lFileName)
        end if

        if (COMMAND_ARGUMENT_COUNT() .eq. 2+9) then
            call getDblArg(3,  TempLower)
            call getDblArg(4,  UpperTemp)
            call getIntArg(5,  nt)
            call getDblArg(6,  PresLower)
            call getDblArg(7,  UpperPres)
            call getIntArg(8,  np)
            call getDblArg(9,  TminSat)
            call getDblArg(10, TmaxSat)
            call getIntArg(11, npsat)
        end if

    end subroutine checkArg

    subroutine getFilesConfig()
        USE RGP_DATA
        implicit none

        ! Check if lOutfils is set, or get from cmd line
        i = SCAN(lOutfile,'.')
        if (i == 0) then
            write(*,'(A)') 'Output Data File Name (i.e. filename.txt):'
            read (*,*) lOutfile
        else
            write(*,*) 'Output Data File: ', lOutfile
        end if

        ! Check if lFileName is set, or get from cmd line
        i = SCAN(lFileName,'.')
        if (i == 0) then
            write(*,'(A)') 'Fluid file name (i.e. filename.fld):'
            read (*,*) lFileName
        else
            write(*,*) 'Fluid file name:  ', lFileName
        end if

        ! Trim excess spaces
        lOutfile = TRIM(lOutfile)
        lFileName = TRIM(lFileName)

        ! Set lfluid from the file name
        i = SCAN(lFileName,'.')
        lfluid = lFileName(1:i-1)

        ! Create write unit for output file
        OPEN(UNIT=10,FILE=lOutFile)
    end subroutine getFilesConfig

    subroutine setupREFPROP()
        USE RGP_DATA
        implicit none

        ! Set fluid to load
        i=1
        hf(1)=lFileName
        
        ! default fluid path is fluids/
        call SETPATH ('fluids/')

        ! Use default mixtures setup and
        ! fluids preferred (DEFault) EOS
        hfmix='hmx.bnc'
        hrf='DEF'

        ! Call refprop setup
        call SETUP (i,hf,hfmix,hrf,ierr,herr)
        call checkRPError(ierr, herr)

        write(*,*) 'REFPROP loaded'

        call INFO (1,wm,ttp,tnbp,tc,pc,dc,zc,acf,dip,rgas)
        call NAME (1,hnam,hn80,hcasn)
        call PASSCMN ('ptpn',0,1,0,hstr,ilng,tpp,arr,ierr,herr)
        call MAXT (x,tm,pm,Dm,ierr,herr)
        call LIMITS(htyp,x,tmin,tmax,Dmax,pmax)

        dc   = wm * dc
        Dmax = wm * Dmax
        hnam = TRIM(hnam);

        write(*,*)
        write(*,*) 'Short name:                        ',hnam 
        write(*,*) 'Molecular Weight:                  ',wm,'g/mol'
        write(*,*) 'Triple Point Temperature:          ',ttp,'K'
        write(*,*) 'Triple Point Pressure:             ',tpp,'kPa'
        write(*,*) 'Normal Boiling Point:              ',tnbp,'K'

        write(*,*) 'Critical Temperature:              ',tc,'K'
        write(*,*) 'Critical Pressure:                 ',pc,'kPa'
        write(*,*) 'Critical Density:                  ',dc,'kg/m^3'
        write(*,*) 'Compressibility at Critical Point: ',zc
        write(*,*) 'Gas Constant:                      ',rgas,'J/mol-K'

        write(*,*)
        write(*,*)'Limits for this fluid:'
        write(*,*)'Min Temperature:                    ',tmin,'K'
        write(*,*)'Max Temperature:                    ',tmax,'K'
        write(*,*)'Max Saturation temperature:         ',tm,'K'
        write(*,*)'Max Density:                        ',Dmax,'kg/m^3'
        write(*,*)'Max Pressure:                       ',pmax,'kPa'
    end subroutine setupREFPROP

    subroutine getRGPBounds()
        USE RGP_DATA
        implicit none

        ! PROMPT FOR BOUNDS AND STEP SIZE
        write(*,*)'--------------------------------------------'
        write(*,*)
        call askQuestion('Lower bound for Temperature [K]:', TempLower, tc, tmax)
        call askQuestion('Upper bound for Temperature [K]:', UpperTemp, TempLower, tmax)
        call askQuestion_int('Number of Steps for Temperature:', nt)

        call askQuestion('Lower bound for Pressure [Pa]:', PresLower)
        call askQuestion('Upper bound for Pressure [Pa]:', UpperPres)
        call askQuestion_int('Number of Steps for Pressure :', np)

        PresLower = max(PresLower, DBLE(0.0));
        UpperPres = min(UpperPres, pmax*1000);

        call askQuestion('Minimum Temperature for Saturation Tables [k]:', TminSat, tmin, tc)
        call askQuestion('Maximum Temperature for Saturation Tables [k]:', TmaxSat, TminSat, tc)
        call askQuestion_int('Number of Steps for Saturation Tables        :', npsat)
    
        ! Tmin and Tmax ranges are fairly strict, if violated do the following
        TminSat = max(TminSat, tmin+DBLE(0.01));
        TmaxSat = min(TmaxSat, tm-DBLE(0.01));
        
        ! set t and p step
        tstep = (UpperTemp-TempLower)/REAL(nt-1)
        pstep = (UpperPres-PresLower)/REAL(np-1)
    end subroutine getRGPBounds

    !---------------------------------------------------------------------!
    !                       DATASET CALCULATION                           !
    !                                                                     !
    ! getSatLineData   : Gets information on the saturation line          !
    ! getSuperDataTP   : Gets information at any vapour point             !
    ! getSuperDataPQ   : calls getSatLineData for super table             !
    !                                                                     !
    !---------------------------------------------------------------------!


    subroutine calculateSuperTable()
        USE RGP_DATA
        implicit none

        prevStatus = 0
        currStatus = 0
        tenPercent = int(UpperPres-PresLower)/10
        j=1
        p(j) = PresLower
        do while (j <= np)
            i=1
            t(i) = TempLower
            
            ! show progress every 10%
            currStatus = int(p(j) - PresLower)/tenPercent
            if (prevStatus /= currStatus) then
                write (*,*) (currStatus*10), '%'
                prevStatus = currStatus
            end if
            
            do while (i <= nt)

                ! if p(j) <= pc (below critical pressure), calc Phi saturated line
                if (p(j) <= pc*1000) then
                    PSat = p(j)/1000
                    q    = 1

                    ! find the sat temp for this pressure
                    call PQFLSH (PSat,q,x,kq,tvap,rhov,Dl,Dv,xliq,xvap,e,enthalpyv,entropyv,convv,conpv,speedOsoundv,ierr,herr)
                    call checkRPError(ierr, herr)
                    
                    ! if t(i) <= tvap (below sat temp), use phi(Tsat, p), q=1 line
                    if (t(i) <= tvap) then
                        call getSuperDataPQ(PSat, h(i,j), w(i,j), v(i,j), cv(i,j), cp(i,j),&
                        &                   s(i,j), dPdv(i,j),eta(i,j), tcx(i,j))
                    
                    ! if t(i) > tvap (above sat temp), get vapour properties
                    else
                        call getSuperDataTP(t(i), p(j), h(i,j), w(i,j), v(i,j), cv(i,j), cp(i,j),&
                        &                   s(i,j), dPdv(i,j), eta(i,j), tcx(i,j))
                    end if
                
                ! if p(j) > pc (above critical pressure), get vapour properties
                else
                    call getSuperDataTP(t(i), p(j), h(i,j), w(i,j), v(i,j), cv(i,j), cp(i,j),&
                    &                   s(i,j), dPdv(i,j),eta(i,j), tcx(i,j))
                end if
                
                i = i + 1
                if (i <= nt) then
                    t(i) = t(i-1)+tstep
                end if
            end do
            
            ! Fill 1D Table (Sat line)
            if (p(j)<(Pc*1000))then
                q    = 1
                PSat = p(j)/1000
                call getSuperDataPQ(PSat, hsat(j), wsat(j), vsat(j), cvsat(j), cpsat(j),&
                &                   ssat(j), dPdvsat(j),etasat(j), tcxsat(j), tsat(j))
            else
                ! fill with phi(Tmin, P)
                q=1
                PSat    = p(j)/1000
                tsat(j) = TempLower
                call getSuperDataTP(TempLower, p(j), hsat(j), wsat(j), vsat(j), cvsat(j), cpsat(j),&
                &                   ssat(j), dPdvsat(j),etasat(j), tcxsat(j))
            end if
            
            j = j+1
            if (j <= np) then
                p(j) = p(j-1)+pstep
            end if
        end do

    end subroutine calculateSuperTable


    subroutine calculateSatTable()
        write (*,*) 'Filling saturation table...'

        ! Get pressure step
        q=1
        call TQFLSH (TminSat,q,x,kq,PminSat,density,Dl,Dv,x,y,e,enthalpy,entropy,conv,conp,speedOsound,ierr,herr)
        call checkRPError(ierr, herr)

        call TQFLSH (TmaxSat,q,x,kq,PmaxSat,density,Dl,Dv,x,y,e,enthalpy,entropy,conv,conp,speedOsound,ierr,herr)
        call checkRPError(ierr, herr)

        pstepsat = (Pmaxsat-Pminsat)/REAL(npsat-1)

        ! Iterate
        j=1
        psat2(j) = PminSat
        do while (j <= npsat)

            q=0
            call getSatLineData(psat2(j), hsatl(j), wsatl(j), rhosatl(j), cvsatl(j), cpsatl(j),&
            &                   ssatl(j), dPdrhosatl(j), etasatl(j), tcxsatl(j))
            
            q=1
            call getSatLineData(psat2(j), hsatv(j), wsatv(j), rhosatv(j), cvsatv(j), cpsatv(j),&
            &                   ssatv(j), dPdrhosatv(j), etasatv(j), tcxsatv(j), tsat2(j))

            blank1(j)=0
            blank2(j)=0

            j = j+1
            if (j <= npsat) then
                psat2(j) = psat2(j-1) + pstepsat
            end if
        end do
    end subroutine calculateSatTable

    !---------------------------------------------------------------------!
    !                    FLID STATE CALCULATION                           !
    !                                                                     !
    ! getSatLineData   : Gets information on the saturation line          !
    ! getSuperDataTP   : Gets information at any vapour point             !
    ! getSuperDataPQ   : calls getSatLineData for super table             !
    !                                                                     !
    !---------------------------------------------------------------------!

    subroutine getSatLineData(psat2_j, hsat_q, wsat_q, rhosat_q, cvsat_q, cpsat_q, ssat_q, dPdrhosat_q, etasat_q, tcxsat_q, tsat_q)
        USE RGP_DATA
        implicit none

        double precision, intent(in) :: psat2_j
        double precision, intent(out) :: hsat_q, wsat_q, rhosat_q, cvsat_q, cpsat_q, ssat_q, dPdrhosat_q, etasat_q, tcxsat_q
        double precision, optional, intent(out) :: tsat_q
        double precision :: rho

        call PQFLSH (psat2_j,q,x,kq,tvap,rho,Dl,Dv,xliq,xvap,e,enthalpy,entropy,conv,conp,speedOsound,ierr,herr)
        call checkRPError(ierr, herr)

        call DPDD (tvap,rho,x,dPdrho)
        call checkRPError(ierr, herr)

        call TRNPRP (tvap,rho,x,viscosity,thermalconductivity,ierr,herr)
        call checkRPError(ierr, herr)

        rhosat_q    = rho*wm
        hsat_q      = enthalpy*1000/wm
        ssat_q      = entropy*1000/wm
        cvsat_q     = conv*1000/wm
        cpsat_q     = conp*1000/wm
        wsat_q      = speedOsound

        dPdrhosat_q = wm/(dPdrho*1000)
        etasat_q    = viscosity/1000000
        tcxsat_q    = thermalconductivity

        if (present(tsat_q)) tsat_q = tvap
    end subroutine getSatLineData


    subroutine getSuperDataPQ(p_s, h_s, w_s, v_s, cv_s, cp_s, s_s, dPdv_s, eta_s, tcx_s, tvap_s)
        USE RGP_DATA
        implicit none

        double precision, intent(in)  :: p_s
        double precision, intent(out) :: h_s, w_s, v_s, cv_s, cp_s, s_s, dPdv_s, eta_s, tcx_s
        double precision, optional, intent(out) :: tvap_s
        double precision :: rho, dPdrhosat_q

        call getSatLineData(p_s, h_s, w_s, rho, cv_s, cp_s, s_s, dPdrhosat_q, eta_s, tcx_s, tvap_s)

        v_s    =  1/rho
        dPdv_s = -((rho*wm)**2)/dPdrhosat_q  ! (dPdrho/wm*1000*(-(rho*wm)**2))
    end subroutine getSuperDataPQ


    subroutine getSuperDataTP(t_s, p_s, h_s, w_s, v_s, cv_s, cp_s, s_s, dPdv_s, eta_s, tcx_s)
        USE RGP_DATA
        implicit none

        double precision, intent(in)  :: t_s, p_s
        double precision, intent(out) :: h_s, w_s, v_s, cv_s, cp_s, s_s, dPdv_s, eta_s, tcx_s
        double precision :: rho

        call TPFLSH (t_s,p_s/1000,x,rho,dl,dv,xliq,xvap,q,e,enthalpy,entropy,conv,conp,speedOsound,ierr,herr)
        call checkRPError(ierr, herr)

        call DPDD (t_s,rho,x,dPdrho)
        call checkRPError(ierr, herr)

        call TRNPRP (t_s,rho,x,viscosity,thermalconductivity,ierr,herr)
        call checkRPError(ierr, herr)

        h_s    = enthalpy*1000/wm
        w_s    = speedOsound
        v_s    = 1/(rho*wm)
        cv_s   = conv*1000/wm
        cp_s   = conp*1000/wm
        s_s    = entropy*1000/wm
        dPdv_s = (dPdrho/wm*1000*(-(rho*wm)**2))
        eta_s  = viscosity/1000000
        tcx_s  = thermalconductivity
    end subroutine getSuperDataTP


    !---------------------------------------------------------------------!
    !                   TABLE WRITING FUNCTIONS                           !
    !                                                                     !
    ! writeRGPHeader    : Writes RGP main and sub header                  !
    !                                                                     !
    ! See below for other table writing functions                         !
    !---------------------------------------------------------------------!


    subroutine writeRGPHeader()
        USE RGP_DATA
        implicit none

        ! SETUP outPUT FILE 
        write(10,901)'$$$$HEADER'

        do i=1,2
            if (i .eq. 2) write(10,901)'$$$$DATA'

            write(10,901)'$$$'//TRIM(lfluid)       !character*8 (key into the RGP file see NOTE 3)  if this program is used with a fluid other than CO2, this line will need to be changed, along with the corresponding line in the $$$DATA header
            write(10,*) 1                          !enter an integer (this line is ignored in ANSYS CFX)
            write(10,901)'$$PARAM'                 ! See NOTE 4
            write(10,*) 28                         !enter an integer(number of parameters)
            write(10,901)'DESCRIPTION'             ! 1
            write(10,901)TRIM(lfluid)//' from NIST'!character*50(description of the material)
            write(10,901)'NAME'                    ! 2
            write(10,901)TRIM(lfluid)              !character*8(material name, same as $$$<component>)
            write(10,901)'INDEX'                   ! 3
            write(10,901)TRIM(lfluid)              !character*50 (index into clients RGDB program)     
            write(10,901)'DATABASE'                ! 4
            write(10,901)'NIST REFPROP'
            write(10,901)'MODEL'                   ! 5                 
            write(10,*)3                           !integer (level of property info available 1,2,or 3)
            write(10,901)'UNITS'                   ! 6                   
            write(10,*)1                           !integer(unit system of 1 ,2,3,4 or 5)
            write(10,901)'PMIN_SUPERHEAT'          ! 7
            write(10,900)PresLower                 !real (Pa) at Triple Point Pressure?
            write(10,901)'PMAX_SUPERHEAT'          ! 8
            write(10,900)UpperPres                 !real (Pa) at Critical Pressure?
            write(10,901)'TMIN_SUPERHEAT'          ! 9
            write(10,900)TempLower                 !real (K) at  Triple Point Temperature?
            write(10,901)'TMAX_SUPERHEAT'          ! 10
            write(10,900)UpperTemp                 !real (K) at Critical Temperature?
            write(10,901)'TMIN_SATURATION'         ! 11
            write(10,900)TminSat
            write(10,901)'TMAX_SATURATION'         ! 12
            write(10,900)TmaxSat
            write(10,901)'SUPERCOOLING'            ! 13
            write(10,900)SuperCool                 !real(supercooling level in superheat tables)(Optional)
            write(10,901)'P_CRITICAL'              ! 14
            write(10,900)Pc*1000                   !real
            write(10,901)'P_TRIPLE'                ! 15 See NOTE 7
            write(10,900)tpp*1000                  !real
            write(10,901)'T_CRITICAL'              ! 16
            write(10,900)Tc                        !real
            write(10,901)'T_TRIPLE'                ! 17
            write(10,900)Ttp                       !real
            write(10,901)'GAS_CONSTANT'            ! 18
            write(10,900)rgas*1000/wm              !real
            
            write(10,901)'TABLE_1'                 ! 19 See NOTE 8
            write(10,*)nt,np                       !integer
            write(10,901)'TABLE_2'                 ! 20
            write(10,*)nt,np                       !integer
            write(10,901)'TABLE_3'                 ! 21
            write(10,*)nt,np                       !integer
            write(10,901)'TABLE_4'                 ! 22
            write(10,*)nt,np                       !integer
            write(10,901)'TABLE_5'                 ! 23
            write(10,*)nt,np                       !integer
            write(10,901)'TABLE_6'                 ! 24
            write(10,*)nt,np                       !integer
            write(10,901)'TABLE_7'                 ! 25
            write(10,*)nt,np                       !integer
            write(10,901)'TABLE_8'                 ! 26 See NOTE 9
            write(10,*)nt,np                       !integer
            write(10,901)'TABLE_9'                 ! 27
            write(10,*)nt,np                       !integer
            write(10,901)'SAT_TABLE'               ! 28
            write(10,*)npsat,4,9                   !integer
            
            if (i .eq. 2) write(10,901)'$$SUPER_TABLE'
            if (i .eq. 2) write(10,*)9 !integer (number of superheat tables, nn = 9)

            ! Output format units
            900   format (5ES17.7E3)
            901   format (A)
        end do

    end subroutine writeRGPHeader


    !---------------------------------------------------------------------!
    !                           SUPERHEAT TABLE                           !
    !                                                                     !
    ! writeSuperTable   : Writes all super table items                    !
    ! writeSupTableHdr  : Adds info about the super table                 !
    ! writeTable(...)   : Writes the given data to file                   !
    !                                                                     !
    !---------------------------------------------------------------------!

    subroutine writeSuperTable()
        call writeTable1()
        call writeTable2()
        call writeTable3()
        call writeTable4()
        call writeTable5()
        call writeTable6()
        call writeTable7()
        call writeTable8()
        call writeTable9()
    end subroutine writeSuperTable

    subroutine writeSupTableHdr(tableno)
        USE RGP_DATA
        implicit none
        character(len=1), intent(in) :: tableno

        write(10,901)'$TABLE_'//tableno 
        write(10,*)nt,np!integer(size of superheat arrays, See NOTE 10)
        write(10,900) (t(i),i=1,nt) !|all following real
        write(10,900) (p(j),j=1,np) !|

        ! Output format units
        900   format (5ES17.7E3)
        901   format (A)
    end subroutine writeSupTableHdr

    subroutine writeTable1()
        USE RGP_DATA
        implicit none

        call writeSupTableHdr('1')
        write(10,900) ((h(i,j),i=1,nt),j=1,np) !|
        write(10,900) (tsat(i),i=1,nsat) !|
        write(10,900) (hsat(i),i=1,nsat) !|

        ! Output format units
        900   format (5ES17.7E3)
    end subroutine writeTable1

    subroutine writeTable2()
        USE RGP_DATA
        implicit none

        call writeSupTableHdr('2')
        write(10,900) ((w(i,j),i=1,nt),j=1,np) !|
        write(10,900) (tsat(i),i=1,nsat) !|
        write(10,900) (wsat(i),i=1,nsat) !|

        ! Output format units
        900   format (5ES17.7E3)
    end subroutine writeTable2

    subroutine writeTable3()
        USE RGP_DATA
        implicit none

        call writeSupTableHdr('3')
        write(10,900) ((v(i,j),i=1,nt),j=1,np) !|
        write(10,900) (tsat(i),i=1,nsat) !|
        write(10,900) (vsat(i),i=1,nsat) !|

        ! Output format units
        900   format (5ES17.7E3)
    end subroutine writeTable3

    subroutine writeTable4()
        USE RGP_DATA
        implicit none

        call writeSupTableHdr('4')
        write(10,900) ((cv(i,j),i=1,nt),j=1,np) !|
        write(10,900) (tsat(i),i=1,nsat) !|
        write(10,900) (cvsat(i),i=1,nsat) !|

        ! Output format units
        900   format (5ES17.7E3)
    end subroutine writeTable4

    subroutine writeTable5()
        USE RGP_DATA
        implicit none

        call writeSupTableHdr('5')
        write(10,900) ((cp(i,j),i=1,nt),j=1,np) !|
        write(10,900) (tsat(i),i=1,nsat) !|
        write(10,900) (cpsat(i),i=1,nsat) !|

        ! Output format units
        900   format (5ES17.7E3)
    end subroutine writeTable5

    subroutine writeTable6()
        USE RGP_DATA
        implicit none

        call writeSupTableHdr('6')
        write(10,900) ((dPdv(i,j),i=1,nt),j=1,np) !|
        write(10,900) (tsat(i),i=1,nsat) !|
        write(10,900) (dPdvsat(i),i=1,nsat) !|

        ! Output format units
        900   format (5ES17.7E3)
    end subroutine writeTable6

    subroutine writeTable7()
        USE RGP_DATA
        implicit none

        call writeSupTableHdr('7')
        write(10,900) ((s(i,j),i=1,nt),j=1,np) !|
        write(10,900) (tsat(i),i=1,nsat) !|
        write(10,900) (ssat(i),i=1,nsat) !|

        ! Output format units
        900   format (5ES17.7E3)
    end subroutine writeTable7

    subroutine writeTable8()
        USE RGP_DATA
        implicit none

        call writeSupTableHdr('8')
        write(10,900) ((eta(i,j),i=1,nt),j=1,np) !|
        write(10,900) (tsat(i),i=1,nsat) !|
        write(10,900) (etasat(i),i=1,nsat) !|

        ! Output format units
        900   format (5ES17.7E3)
    end subroutine writeTable8

    subroutine writeTable9()
        USE RGP_DATA
        implicit none

        call writeSupTableHdr('9')
        write(10,900) ((tcx(i,j),i=1,nt),j=1,np) !|
        write(10,900) (tsat(i),i=1,nsat) !|
        write(10,900) (tcxsat(i),i=1,nsat) !|

        ! Output format units
        900   format (5ES17.7E3)
    end subroutine writeTable9

    !---------------------------------------------------------------------!
    !                           SATURATION TABLE                          !
    ! Note: 12.7 last paragraph:                                          !
    ! The $$SAT_TABLE section does not have to exist under a $$$Database  !
    ! access key, that the SUPERCOOLING parameter can be zero and that    !
    ! only the $$SUPER_TABLE section is necessary.                        !
    !                                                                     !
    ! writeSatTable     : Adds the sat table to the RGP file              !
    ! writeSatData(...) : Writes the given data to file                   !
    !                                                                     !
    !---------------------------------------------------------------------!
    
    subroutine writeSatTable()
        USE RGP_DATA
        implicit none

        write (*,*) 'Sat Table'
        write(10,901)'$$SAT_TABLE'
        write(10,*)npsat,4,9 !(size of saturation arrays, see NOTE 11)
        write(10,900) (psat2(i)*1000,i=1,npsat)
        write(10,900) (tsat2(i),i=1,npsat)
        write(10,900) (blank1(i),i=1,npsat)
        write(10,900) (blank2(i),i=1,npsat)

        ! Liquid group
        call writeSatData(hsatl, cpsatl, rhosatl, dPdrhosatl, ssatl, cvsatl,&
        &                 wsatl, etasatl, tcxsatl)
        
        ! Vapour group
        call writeSatData(hsatv, cpsatv, rhosatv, dPdrhosatv, ssatv, cvsatv, &
        &                 wsatv, etasatv, tcxsatv)

        ! Output format units
        900   format (5ES17.7E3)
        901   format (A)
    end subroutine writeSatTable

    subroutine writeSatData(H, Cp, Rho, dPdRho_T, S, Cv, SpSnd, eta, tcx)
        implicit none
        double precision, allocatable :: H(:), Cp(:), Rho(:), dPdRho_T(:), S(:), Cv(:), SpSnd(:), eta(:), tcx(:)

        write(10,900) (H(i),i=1,npsat)        ! enthalpy
        write(10,900) (Cp(i),i=1,npsat)       ! Cp
        write(10,900) (Rho(i),i=1,npsat)      ! density
        write(10,900) (dPdRho_T(i),i=1,npsat) ! change density with pressure const T
        write(10,900) (S(i),i=1,npsat)        ! entropy 
        write(10,900) (Cv(i),i=1,npsat)       ! Cv
        write(10,900) (SpSnd(i),i=1,npsat)    ! Speed of Sound
        write(10,900) (eta(i),i=1,npsat)      ! molar viscosity
        write(10,900) (tcx(i),i=1,npsat)      ! thermal conductivity

        ! Output format units
        900   format (5ES17.7E3)
    end subroutine writeSatData


    !---------------------------------------------------------------------!
    !                           HELPER FUNCTIONS                          !
    !---------------------------------------------------------------------!
    subroutine askQuestion(question, response, minimum, maximum)
        implicit none
        character(len=*), intent(in) :: question
        double precision, intent(inout) :: response
        double precision, optional, intent(in) :: minimum, maximum

        if (response .eq. 0) then
            write(*,*) 'Enter - ', question
            read(*,*) response
        else
            write(*,*) question, response
        end if
        
        if (present(minimum) .and. response .lt. minimum) then
            write(*,*) response, 'below below minimum of', minimum, 'this will be used instead'
            response = minimum
        end if
        
        if (present(maximum) .and. response .gt. maximum) then
            write(*,*) response, 'above maximum of', maximum, 'this will be used instead'
            response = maximum
        end if
    end subroutine askQuestion

    subroutine askQuestion_int(question, response)
        implicit none
        character(len=*), intent(in) :: question
        integer, intent(inout) :: response

        if (response .eq. 0) then
            write(*,*) 'Enter - ', question
            read(*,*) response
        else
            write(*,*) question, response
        end if
    end subroutine askQuestion_int

    subroutine getDblArg(position, value)
        implicit none
        character(len=255)     :: arg
        integer, intent(in)  :: position
        double precision, intent(out) :: value
        
        call get_command_argument(position,  arg)
        read(arg, *)value
    end subroutine getDblArg

    subroutine getIntArg(position, value)
        implicit none
        character(len=255)     :: arg
        integer, intent(in)  :: position
        integer, intent(out) :: value
        
        call get_command_argument(position,  arg)
        read(arg, *)value
    end subroutine getIntArg

    subroutine checkRPError(ierr, herr)
        implicit none;
        integer, intent(in) :: ierr
        character(len=255), intent(in) :: herr

!       ierr value : flash resulted in : the results are
!       ierr > 0   : error             : calculations not possible,
!       ierr < 0   : warning           : results may be questionable
        
        if (ierr .gt. 0) then
            write (*,*) herr
            call exit(-1)
        end if
    end subroutine checkRPError

    subroutine allocateArrays()
        USE RGP_DATA
        implicit none

        ! allocate 2D array
        write (*,*) 'Allocating superheat array...'
        ALLOCATE (t(nt), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 001'
        ALLOCATE (p(np), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 002'
        ALLOCATE (h(nt,np), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 003'
        ALLOCATE (w(nt,np), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 004'
        ALLOCATE (v(nt,np), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 005'
        ALLOCATE (cv(nt,np), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 006'
        ALLOCATE (cp(nt,np), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 007'
        ALLOCATE (dPdv(nt,np), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 008'
        ALLOCATE (s(nt,np), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 009'
        ALLOCATE (eta(nt,np), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 010'
        ALLOCATE (tcx(nt,np), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 011'

        ! allocate 1D array
        nsat = np
        ALLOCATE (tsat(nsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 012'
        ALLOCATE (hsat(nsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 013'
        ALLOCATE (wsat(nsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 014'
        ALLOCATE (vsat(nsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 015'
        ALLOCATE (cvsat(nsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 016'
        ALLOCATE (cpsat(nsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 017'
        ALLOCATE (dPdvsat(nsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 018'
        ALLOCATE (ssat(nsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 019'
        ALLOCATE (etasat(nsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 020'
        ALLOCATE (tcxsat(nsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 021'

        ALLOCATE (tsat2(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 022'
        ALLOCATE (psat2(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 023'
        ALLOCATE (blank1(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 024'
        ALLOCATE (blank2(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 025'

        ALLOCATE (hsatl(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 026'
        ALLOCATE (wsatl(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 027'
        ALLOCATE (rhosatl(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 028'
        ALLOCATE (cvsatl(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 029'
        ALLOCATE (cpsatl(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 030'
        ALLOCATE (dPdrhosatl(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 031'
        ALLOCATE (ssatl(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 032'
        ALLOCATE (etasatl(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 033'
        ALLOCATE (tcxsatl(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 034'

        ALLOCATE (hsatv(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 035'
        ALLOCATE (wsatv(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 036'
        ALLOCATE (rhosatv(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 037'
        ALLOCATE (cvsatv(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 038'
        ALLOCATE (cpsatv(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 039'
        ALLOCATE (dPdrhosatv(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 040'
        ALLOCATE (ssatv(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 041'
        ALLOCATE (etasatv(npsat), STAT = AllocStatus)
        if (AllocStatus /= 0) write (*,*)'out of memory 042'
        ALLOCATE (tcxsatv(npsat), STAT = AllocStatus) 
        if (AllocStatus /= 0) write (*,*)'out of memory 043'
    end subroutine allocateArrays

end program RGP
