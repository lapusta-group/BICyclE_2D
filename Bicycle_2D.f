c--------------------------------------------------------------------------
c BICyclE (Boundary Integral for Cycles of Earthquakes)
c 2D version parallelized with MPI
c
c This code was originally written by Nadia Lapusta,
c and subsequently modified by Hiroyuki Noda and Valere Lambert
c
c Things to consider modifying:
c - Parameter variables 
c - Time related parameters 
c - Output related parameters 
c - Output points (controls where the "H###" files are recorded)
c - Values in Set_initials -> you can check that these are working properly with the "initialSetup" file 
c
c Important Notes:
c - N1/NPU must be an integer value! <- make sure to make consistent with number of CPUs in run.sh if changed here
c - *pr : previous values
c - *s  : predicted values
c - b*  : initial value
c
c-------------------------------------------------------------------------
      PROGRAM Code
      IMPLICIT NONE
      INCLUDE 'mpif.h'						        
      INCLUDE 'fftw3.f'    
c
c-----Definition of the variables-----------------------------------------
      CHARACTER(80) :: jobname,ach	         	        ! jobname used for files output for the whole sim, ach used to hold the numbering scheme for filenames related to specific events
      CHARACTER(200) :: wdir                                    ! output directory
      LOGICAL cyclecont						! tells the basic cycle whether to run again or not
      INTEGER i,effi,ierror,ihconv,ihconv2,j,jbeta,jshift, ! ierror is error flag, 0 = no error, 1 = Nlnew is too small, 2 = NRSearch error
     &     jdtmax,jshiftpr,jTw,N1,N1fr,Nlnew,ndt,ndtmax,	
     &     k,vi,niter,NRiter,iiter,mode,FS,	        
     &     N1FULL,NPU,N1X,NSUBX,Enumax,outStepSize              
      REAL*8 cs,cp,delay,dt,dtelas,dtincf,dtmax,dtmin,dtpr,
     &     divN1,eta,gamma,hcell,pi,q,nu,rau,rho,	 
     &     rhoexpr,t,tmax,xLam,xLam0,xfr,xfrvw,xmu,Vo,Vpl,betas,	         
     &     Hseis,Sfac,x0fr,x0frvw,tol
      INTEGER ierr,myrank,size
      EXTERNAL bej1sx,bessj0,xker2,integral,w
c
c-----Parameter variables-----CONTROL THIS AREA---------------------------
c  Parameters governing grid size and frictional domain

      PARAMETER(xLam0=50.d3,xfr=40.d3,xfrvw=10.d3,		! Total domain, frictional domain, velocity weakening domain
     &     FS=1,xLam=FS*xLam0,                                  ! FS : 1 = whole space, 2 = half space - for free surface ( Mode III )
     &     x0fr=0.5d0*(xLam0-xfr)*(FS-1),                       ! Absolute shift for frictional and VW region (free surface)
     &     x0frvw=0.5d0*(xLam0-xfrvw)*(FS-1),                   ! Set to zero when full space
     &     N1=2000,N1fr=N1*nint(xfr/xLam),		        ! N1 = Number of cells, N1fr = number of cells in frictional region
     &     N1FULL=FS*N1,                                        ! N1FULL = total number of cells, for wave numbers in kernels
     &     jbeta=3,						! jbeta controls number of timesteps needed for shear wave to cross a single cell
     &     jdtmax=1.d3,						! jdtmax is approximately the number of timesteps we want to guarantee in the interseismic period - affects dtmax
     &     delay=0.5d0,						! offset needed for mid-point rule integration
     &     gamma = 1.5,						! Factor to adjust jTw (elastodynamic time window)
     &     ndtmax = 1.d9,					! Max number of timesteps for the simulation
     &     Enumax = 5,  					! Max number of events for the simulation
     &     tmax = 1200.d0*3.1536e7,				! Max number of total seconds in simulation (divide number by ~3*10^7 to get years)
     &     jTw=N1*jbeta*gamma*nint(xfrvw/xLam),		  	! Number of Fourier coefficients [Time window used for fourier transforms, good limit is at least N1*jbeta*(VWDomain/TotalDomain)
     &     Vpl=1.d-9,						! Plate rate
     &     NRiter=1000,tol=1.d-6,                               ! Number of iterations and tolerance for Newton-Raphson solver
     &     niter=1,mode=3,Sfac=2.d0,                            ! niter = max number of iterations for determining better guess (of slip, stress, etc.)
								! mode = Mode of Simulation (2,3, or crustal plane), Sfac - adjust timesteps smaller or larger -> timestep*(1/Sfac)
c  mode = 2 : in-plane
c  mode = 3 : anti-plane
c  mode = 4 : crustal plane
     &     NPU=50,NSUBX=N1FULL/NPU,N1X=N1/NPU,                  ! Number of computer nodes used for sim, NSUBX = N1FULL/NPU,  MUST BE AN INTEGER!!!!!
     &     Nlnew=jTw*NSUBX,					! Holds the total number of convolution time windows and is the number of elements that memory has to provide
     &     Hseis = 5.d3,					! Seismogenic depth, only applicable to crustal plane model (mode 4) 
     &     outStepSize = 10)                                    ! spatial interval for "line" output files

      PARAMETER (pi=3.141592653589793d0)
c
c-------Physical parameters------------------------------------------------ 
c Physical properties of the elastic material (rocks)
      PARAMETER(betas=1.d0,                                      ! Scaling parameter for radiation damping term in Quasi-Dyn simulations (See Lapusta + Liu 2009)
     &     cs = 3464.d0*betas,                                  ! Shear wave speed (m/s)
     &     nu = 0.25,						! Poisson's ratio (unitless)
     &     cp = dsqrt(2.d0*(1.d0-nu)/(1.d0-2.d0*nu))*cs,	! Pressure wave speed (m/s)
     &     rau = 2670.d0,					! Density (kg/m^3)
     &     xmu = cs*cs*rau/(betas*betas), 			! Shear Modulus (Pa)
     &     eta = 0.50*xmu/cs,					! Radiation damping term in NR Search
     &     rho = (cp/cs)**2,					! used in kernel computations
     &     rhoexpr = 2.d0*(1.d0-1.d0/rho))			! used in convolutions
c
c------Physical Domain----------------------------------------------------
      PARAMETER(divN1 = 1.d0/dfloat(N1FULL),			! 1/N1
     &     hcell = xLam*divN1,					! cell size
     &     dtelas = hcell/cs/jbeta,				! we want the minimum timestep to be such that it takes "jbeta" timesteps for a shear wave to cross one cell in our simulation
     &     dtmin = dtelas)					! set the minimum timestep equal to the above statement
c
c-----Fault slips, slip rates, and stresses-------------------------------
c   *pr : previous values
c   *s  : predicted values
      REAL*8 del(N1X),dels(N1X),delpr(N1X),delV(N1X),
     &     tau(N1X),taupr(N1X),
     &     V(N1X),Vs(N1X),VV(N1X),Vpr(N1X),
     &     dLmax(N1X),func(N1X),x(N1X),tauo(N1X),
     &     dfdlV(N1X),dfssdlV(N1X)		! dfdlV = delta friction, delta log of slip rate; dfssdlV = delta friction, steady state, delta log of slip rate
c   Corresponding global arrays
      REAL*8, DIMENSION(:), ALLOCATABLE :: delall,
     &     delallpr,funcall,tauall,xall,tauoall,
     &     Vall,dfdlVall,dfssdlVall

      REAL*8 An(NSUBX),dtlocal(NPU)
c   Arrays storing kernels
      REAL*8 ww(Nlnew)
c   Fourier Coeff. of the slip rate history
      complex(kind(0d0)) Vk(Nlnew)
c   Fourier Coeff.s
      complex(kind(0d0)) FCdelV(NSUBX),FCfunc(NSUBX),		! FCdelV = Fourier coeff. of slip, FCfunc = Fourier coeff of functional
     &     FCVV(NSUBX),FCVpr(NSUBX),FCVk(NSUBX)			! FCVV = Fourier coeff. of slip rate, FCVpr = Fourier coeff of previous slip rate, FCVk = Fourier coeff of stored Vk kernel
c   Snh : Head of the convolution bw. wwn and FCVm.
c   Snt : Tail of the convolution bw. wwn and FCVm.
c   termdyn : Dynamic terms in the stress transfer
c   termst  : Static terms in the stress transfer
c   termtot : termdyn + termst
      complex(kind(0d0)) Snh(NSUBX),Snt(NSUBX),Snhss(NSUBX),
     &     termdyn,termst,termtot(NSUBX)
c
c   FFTW-related variables
      REAL*8, DIMENSION(:), ALLOCATABLE :: delVall,VVall      
      COMPLEX*16, DIMENSION(:), ALLOCATABLE :: FCdelVall,
     &     FCVVall,FCfuncall
      COMPLEX*16, DIMENSION(:), ALLOCATABLE :: inFFFTW,
     &     outFFFTW,inIFFTW,outIFFTW
      INTEGER*8 planF,planI
c
c-----Frictional constitutive parameters and variables--------------------
      REAL*8 L(N1X),cca(N1X),ccb(N1X),
     &     psi(N1X),psipr(N1X),psis(N1X),fo(N1X),
     &     Seff(N1X)
c     Corresponding global arrays
      REAL*8, DIMENSION(:), ALLOCATABLE :: Lall,ccaall,ccball,
     &     psiall,foall,Seffall

c     Minimum value of L
      REAL*8 minL,LLmin

c     Material number for choosing friction law
      INTEGER Mnum(N1X),nparam
      PARAMETER(nparam=5)
      REAL*8 param(nparam)
c
      INTEGER jTwn(NSUBX),jnew(3,NSUBX)
      INTEGER jTwnall(N1FULL)
      REAL*8 wwsum(NSUBX)
c
c-----Output-related variables------------------------------------
      INTEGER noutH						! Number of observation points
      PARAMETER(noutH=16)					
      INTEGER P1(noutH),P1r(noutH),Num_procP(noutH)
      REAL*8 tvsx,tvsxinc,Vmax,Vmaxpr,tevne,tevneinc,
     &     tevneb,Vmaxlocal,Temaxlocal,Temax
      REAL*8 Vmaxin(2), Vmaxout(2)
      INTEGER Vmaxloc, Temaxloc
      INTEGER maxVIndex, SaveTime 
      REAL*8 tvkininc,tvkin
      REAL*8 Vevne
      INTEGER idelevne
c
c-----Event catalogue output -------------------------------------
      INTEGER Endtin,Endten
      REAL*8 Etin,Eten
      INTEGER Enum
c     Event descriptors
      Enum = 0
      idelevne = 0
      Etin = 0.d0
      Eten = 0.d0
      Endtin = 0
      Endten = 0

      ALLOCATE(FCdelVall(N1FULL),FCVVall(N1FULL),FCfuncall(N1FULL), 
     &        stat=ierr) 
      ALLOCATE(inFFFTW(N1FULL),outFFFTW(N1FULL),inIFFTW(N1FULL),
     &        outIFFTW(N1FULL),stat=ierr)
      ALLOCATE(delVall(N1),VVall(N1))

c     
c-----------END OF DECLARATION STATEMENTS--------------------------------------
c BELOW THIS LIES EXECUTABLE CODE, DON'T MOVE ANY OF THIS UP INTO THE
c DECLARATION STATEMENT AREA, FORTRAN WON'T LIKE THIS
c
c-----Time related parameters-------CONTROL THIS AREA--------------------------
      dtincf = 2.0d0						! controls the rate at which the timestep increases (by a factor of dtincf)
      Vevne = 1.d-3                                             ! defines boundary between aseismic and seismic slip (m/s), put in order such that 
      dtmax = (Vevne/Vpl)*xLam/cs/jdtmax		        ! max timestep should depend on ratio of interseismic to seismic time
c								! depends on ~(Vseis/Vpl) = ~10^9 and time for shear waves to cross domain
								! jdtmax is how many timesteps we want to "guarantee" in the interseismic time

      IF (dtmax.GT.dtmin*dfloat(jTw)) then			! double check to make sure that jshift is created correctly given the size of dtmax
         jshift=jTw						! if dtmax is bigger than dtmin*jTw then we're going to shift by the entire elastodynamic time window
      ELSE							! otherwise, we shift the elastodynamic time window by a fraction of dtmax/dtmin
         jshift=dtmax/dtmin
         IF(jshift.eq.0) jshift=1				! check to make sure that jshift >= 1
         dtmax=jshift*dtmin					! set dtmax = dtmin (it's not allowed to be smaller)
      ENDIF

c------Output related parameters-------CONTROL THIS AREA-----------------------

      tvsxinc = 2*3.1536d7					! line_c file is output with this time period
      tvkininc = 8.64d4                                         ! Cycle VTD is output with this time period
      tvsx = tvsxinc						! current time for line_c output
      tvkin = tvkininc                                           
      Vmax = 0.d0						! maximum velocity on fault
      Vmaxpr = 0.d0
      tevne = 0.d0						! current time in event (used for line_e output)
      tevneinc = 0.05d0						! line_e files are output every "tevneinc" seconds
      tevneb = 0.d0						! records time at beginning of event (so we can have relative event times rather than since the start of the sim) 
      ierror = 0						! Let's start off with no errors      

c------------------------------------------------------------------------------
c     Initialize MPI                                                                                                                                                                                       
      CALL MPI_Init(ierr)
      CALL MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
      CALL MPI_Comm_size(MPI_COMM_WORLD,size,ierr)

      CALL readinput(wdir,myrank)
      CALL FLUSH(6)


c-------- FFTW Initiation-----------------------------------------------------
      CALL dfftw_plan_dft_1d(planF,N1FULL,inFFFTW,outFFFTW,
     &                                     FFTW_FORWARD, FFTW_ESTIMATE)
      CALL dfftw_plan_dft_1d(planI,N1FULL,inIFFTW,outIFFTW,
     &                                     FFTW_BACKWARD,FFTW_ESTIMATE)
c
c---------Allocate arrays and create output files for whole simulation--------
c Files here use ones (1_) in their names so they are sorted
c to the front of the files and don't get mixed in with files for each event

      IF(myrank.EQ.0)THEN
c     Allocate arrays for head node
         ALLOCATE(Lall(N1),ccaall(N1),ccball(N1),
     &        psiall(N1),foall(N1),
     &        Seffall(N1),stat=ierr)
         
         ALLOCATE(delall(N1),delallpr(N1),funcall(N1),tauall(N1),
     &        xall(N1),tauoall(N1),Vall(N1),dfdlVall(N1),dfssdlVall(N1),
     &        stat=ierr)

         IF (0/=ierr) STOP "could not allocate memory"
         
         OPEN(100, FILE = trim(wdir)//'/1_line_c'//jobname) ! line file for interseismic period
         OPEN(101, FILE = trim(wdir)//'/1_vmax'//jobname) ! maximum slip rate and temperature
         OPEN(102, FILE = trim(wdir)//'/1_initialSetup'//jobname) ! output of initial setup of fault to check and see if you're doing the simulation you want to 
         OPEN(103, FILE = trim(wdir)//'/1_ErReport'//jobname) ! Error report
         OPEN(104, FILE = trim(wdir)//'/1_HpointLocations'//jobname) ! Output locations of H points 
         OPEN(105, FILE = trim(wdir)//'/1_settings'//jobname) ! Output physical settings
         OPEN(106, FILE = trim(wdir)//'/1_Ecatalogue'//jobname) ! Event information
      ENDIF
c------Physical Domain--------------------------------------------------
      DO i=1,N1X
         x(i)=(dfloat(myrank*N1X+i-N1/2) - 0.5d0)*hcell
      ENDDO
      
      CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
      CALL MPI_Gather(x,   N1X,MPI_REAL8,
     &                xall,N1X,MPI_REAL8,
     &                0,MPI_COMM_WORLD,ierr)

c------Initial conditions------------------------------------------------
      CALL Set_initials(Seff,tauo,V,del,			! Set the initial conditions on the fault
     &     psi,Mnum,cca,ccb,L,fo,Vo,eta,
     &     nparam,x,N1X,xfr,xfrvw,xLam,x0fr,x0frvw)


      DO i = 1,N1X						! Set initial values on each nodes
         delpr(i) = del(i)
         dels(i) = del(i)
         psipr(i) = psi(i)
         psis(i) = psi(i)
c   Solve for tau and V consistent with tauo and psi
         IF(dabs(x(i)-x0fr).le.(0.5d0*xfr))then			! Solving only on the frictionally resolved regime
            CALL const_param(param,nparam,fo(i),Vo,
     &           cca(i),ccb(i),L(i))

            CALL NRsearch(V(i),tau(i),
     &           tauo(i),func(i),
     &           Seff(i),Vpl,eta,psi(i),
     &           param,nparam,Mnum(i),ierror,x(i),
     &           NRiter,tol)
         ENDIF
         Vpr(i) = V(i)
         Vs(i) = V(i)
      ENDDO
      DO i=1,NSUBX
         FCVpr(i) = 0.0d0
         FCVk(i) = 0.0d0
         Snt(i) = 0.d0
         Snh(i) = 0.d0
         Snhss(i) = 0.d0
      ENDDO
      DO i = 1, Nlnew
         Vk(i) = 0.d0
      ENDDO
      DO i = 1, NPU
         dtlocal(i) = 0
      ENDDO


c
c-----Output points--------------------------------------
      DO i=1,noutH
         P1(i) = nint(N1*(xLam0 - (i-1)*2.4d3)/xLam0)
      ENDDO

c     Location of the H### points
      IF (myrank.eq.0) then
         DO i = 1,noutH
            WRITE(104,517) P1(i), xall(P1(i)) 
         ENDDO
c     Output key settings
         WRITE(105,518) N1, xLam0, xfr, xfrvw, cs, xmu
      ENDIF
c
      DO i = 1, noutH
         CALL Mode_convertor(P1(i),N1X,Num_procP(i),P1r(i))	! Divide the output points amongst the nodes
         IF (myrank.EQ.Num_procP(i)) THEN			! Num_procP gives the node number for the point P1(i)
            IF (i.LE.9) THEN					! and P1r(i) is the point # out of "NSUBX" points on that node
               WRITE(ach,400) 0,0,0,i
            ELSEIF (i.LE.99) THEN
               WRITE(ach,401) 0,0,i
            ELSE
               WRITE(ach,402) 0,i
            ENDIF
            OPEN(10+i, FILE = trim(wdir)//'/H'//ach//jobname)
         ENDIF
      ENDDO
c
c------Kernel Computations--------------------------------
      CALL conv_storage(N1FULL,NSUBX,jTw,Nlnew,jobname,jnew,jTwn,
     &     ierror,myrank,wdir)
      IF (ierror.eq.1) then
         STOP
         CALL MPI_FINALIZE(ierr)
      ENDIF
c
      IF(mode.eq.2) then 					! in-plane (This is kernel computation for mode II. Hseis is not used. HIRO 20100903)
         Call kern2(N1FULL,NSUBX,Nlnew,jbeta,delay,jobname,rho,dtelas,
     &     bej1sx,w,bessj0,xker2,integral,jTwn,jnew,ww,wwsum,cs,
     &     xLam,myrank)
         DO i=1,NSUBX
            effi=myrank*NSUBX+i
            IF (effi.le.(N1FULL/2+1)) then
               k = effi - 1
            ELSE
               k = effi - 1 - N1FULL
            ENDIF
            An(i)=xmu*pi*dabs(k/xLam)
         ENDDO
      ELSEIF (mode.eq.3) then 					! anti-plane (This is kernel computation for mode III. Hseis is not used. HIRO 20100903)
         CALL kern3(N1FULL,NSUBX,Nlnew,jbeta,delay,jobname,jTwn,
     &     jnew,ww,wwsum,cs,xLam,dtelas,myrank)
         DO i=1,NSUBX
            effi=myrank*NSUBX+i
            IF (effi.le.(N1FULL/2+1)) then
               k = effi - 1
            ELSE
               k = effi - 1 - N1FULL
            ENDIF
            An(i)=xmu*pi*dabs(k/xLam)
         ENDDO
      ELSEIF (mode.eq.4) then 					! crustal-plane (This is kernel computation for the crustal plane model. Hseis matters here only. HIRO 20100903)
         CALL kern4(N1FULL,NSUBX,Nlnew,jbeta,delay,jobname,jTwn,
     &     jnew,ww,wwsum,cs,nu,xLam,dtelas,Hseis,myrank)
         DO i=1,NSUBX
            effi=myrank*NSUBX+i
            IF (effi.le.(N1FULL/2+1)) then
               k = effi - 1
            ELSE
               k = effi - 1 - N1FULL
            ENDIF
            An(i)=0.5d0*xmu*dsqrt((2.d0*pi*dabs(k/xLam)/(1.d0-nu))**2
     &                                            +(4.d0/pi/Hseis)**2)
         ENDDO
      ENDIF

c     Initialize timing
      dt = dtmin
      jshift = 1
      t = 0
      ndt = 0
      dtpr=dt        
      jshiftpr=jshift
c
c---------- Output initial information ---------------------
      CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
      CALL MPI_Gather(del,   N1X,MPI_REAL8,
     &     delall,N1X,MPI_REAL8,
     &     0,MPI_COMM_WORLD,ierr)
      CALL MPI_Gather(V,   N1X,MPI_REAL8,
     &     Vall,N1X,MPI_REAL8,
     &     0,MPI_COMM_WORLD,ierr)
      CALL MPI_Gather(tau,   N1X,MPI_REAL8,
     &     tauall,N1X,MPI_REAL8,
     &     0,MPI_COMM_WORLD,ierr)
      CALL MPI_Gather(psi,   N1X,MPI_REAL8,
     &     psiall,N1X,MPI_REAL8,
     &     0,MPI_COMM_WORLD,ierr)
      CALL MPI_Gather(L,   N1X,MPI_REAL8,
     &     Lall,N1X,MPI_REAL8,
     &     0,MPI_COMM_WORLD,ierr)
c
      IF(myrank.eq.0)then
        CALL out_vsx(100,N1,					! line_c
     &        delall,Vall,tauall,psiall,
     &        ndt,t,hcell,outStepSize)

      ENDIF


      DO i=1,noutH						! H### files
        IF (myrank.eq.Num_procP(i)) then
          WRITE(10+i,430)
     &           ndt,t,
     &           del(P1r(i)),V(P1r(i)),
     &           tau(P1r(i)),psi(P1r(i))
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     Basic Cycle
c-----------------------------------------------------------------------

      cyclecont = .true.
      DO WHILE(cyclecont)       ! continue basic cycle while cyclecont = 1

c     Find first prediction of slip
         DO i=1,N1X
            dels(i)=delpr(i)+Vpr(i)*dt
            delV(i)=dels(i)-Vpl*(t+dt)
c     State evolution
            IF(dabs(x(i)-x0fr).le.(0.5d0*xfr))then ! if the point is in the frictionally resolved region...
               CALL const_param(param,nparam,fo(i),Vo,
     &              cca(i),ccb(i),L(i))
               CALL s_evolution(psis(i),psipr(i),psipr(i),dt, ! psipr(i) is added as the referred state variable in the second half of the timestep. HIRO 20100903
     &              Vpr(i),Vpr(i),Seff(i),
     &              param,nparam,Mnum(i))
            ENDIF
         ENDDO                  ! end 1,N1X determination of state variables
c     
c     Fourier transform of slip difference
         CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
         CALL MPI_Gather(delV,   N1X,MPI_REAL8,
     &        delVall,N1X,MPI_REAL8,
     &        0,MPI_COMM_WORLD,ierr)
         IF(myrank.eq.0)then
            CALL OneFFTF(delVall,FCdelVall,planF,N1,N1FULL,FS)
         ENDIF
         
         CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
         CALL MPI_Scatter(FCdelVall,NSUBX,MPI_COMPLEX16,
     &        FCdelV,   NSUBX,MPI_COMPLEX16,
     &        0,MPI_COMM_WORLD,ierr)
c     
c     Fix the range of convolution.
         DO i=1,NSUBX
            IF(jshift.ge.jTwn(i)) then
               jnew(3,i)=jnew(1,i)
            ELSE
               ihconv=jnew(3,i)-jnew(1,i)
               IF(ihconv.ge.jshift) then
                  jnew(3,i)=jnew(3,i)-jshift
               ELSE
                  jnew(3,i)=jnew(2,i)-(jshift-ihconv)
               ENDIF
            ENDIF
         ENDDO
c     
c     "Tail" means convolution contribution from the previous timesteps.
         DO i=1,NSUBX
!     Zero tails of the convolution.
            Snt(i)=0.d0
            IF(jshift.LT.jTwn(i)) then 
!     Nonzero tail only in this case, depending on Fourier mode.
               IF(jshiftpr.ge.jTwn(i)) then
                  DO k=jnew(1,i)+1,jnew(2,i) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     Vk(k)=FCVk(i)
                  ENDDO
               ENDIF
               ihconv=jnew(3,i)-jnew(1,i)
               ihconv2=jnew(2,i)-jnew(3,i)
               IF(ihconv2.ge.(jshift+1)) then
                  DO k=jnew(3,i)+jshift+1,jnew(2,i)
                     Snt(i)=Snt(i)+ww(k-ihconv)*Vk(k)
                  ENDDO
                  DO k=jnew(1,i)+1,jnew(3,i)
                     Snt(i)=Snt(i)+ww(k+jTwn(i)-ihconv)*Vk(k)
                  ENDDO
               ELSE
                  DO k=jnew(1,i)+1+(jshift-ihconv2),jnew(3,i)
                     Snt(i)=Snt(i)+ww(k+jTwn(i)-ihconv)*Vk(k)
                  ENDDO
               ENDIF
            ENDIF	
         ENDDO
c     "Head" means the contribution from the current time step.
         DO i=1,NSUBX
            IF(jshift.ge.jTwn(i)) then
               Snh(i)=wwsum(i)
            ELSE
               Snh(i)=0.d0
               ihconv=jnew(3,i)-jnew(1,i)
               IF(jnew(3,i)+jshift.le.jnew(2,i)) then
                  DO k=jnew(3,i)+1,jnew(3,i)+jshift
                     Snh(i)=Snh(i)+ww(k-ihconv)
                  ENDDO
               ELSE
                  DO k=jnew(3,i)+1,jnew(2,i)
                     Snh(i)=Snh(i)+ww(k-ihconv)
                  ENDDO
                  DO k=jnew(1,i)+1,jnew(1,i)+
     &                 jshift-(jnew(2,i)-jnew(3,i))
                     Snh(i)=Snh(i)+ww(k+jTwn(i)-ihconv)
                  ENDDO
               ENDIF
            ENDIF
            Snh(i)=Snh(i)*FCVpr(i)
         ENDDO
c     
         DO i=1,NSUBX
            IF(mode.eq.2)then
               termst=-An(i)*rhoexpr*FCdelV(i)
            ELSEIF((mode.eq.3).or.(mode.eq.4))then
               termst=-An(i)*FCdelV(i)
            ENDIF
            termdyn=An(i)*(Snh(i)+Snt(i))*dtelas
            FCfunc(i)=termst+termdyn
         ENDDO
         
c     Inverse Fourier transform to obtain stress transfer functional
         CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
         CALL MPI_Gather(FCfunc,   NSUBX,MPI_COMPLEX16,
     &        FCfuncall,NSUBX,MPI_COMPLEX16,
     &        0,MPI_COMM_WORLD,ierr)
         IF(myrank.eq.0)then
            call OneFFTB(FCfuncall,funcall,planI,N1,N1FULL)
         ENDIF
         
         CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
         CALL MPI_Scatter(funcall,N1X,MPI_REAL8,
     &        func,   N1X,MPI_REAL8,
     &        0,MPI_COMM_WORLD,ierr)
c     
c     Find first prediction of velocity
         DO i = 1,N1X
            IF(dabs(x(i)-x0fr).le.(0.5d0*xfr))then
               CALL const_param(param,nparam,fo(i),Vo,
     &              cca(i),ccb(i),L(i))
               CALL NRsearch(Vs(i),tau(i),
     &              tauo(i),func(i),
     &              Seff(i),Vpl,eta,psis(i),
     &              param,nparam,Mnum(i),ierror,x(i),
     &              NRiter,tol)
c     Modifiable: Vs,V3s,tau,tau3,ierror
            ENDIF
         ENDDO
c     Compute FCV using the predicted slip rates.
         DO i =1,N1X
            VV(i) = Vs(i) - Vpl
         ENDDO
         
c     Fourier transform slip rate difference
         CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
         CALL MPI_Gather(VV,   N1X,MPI_REAL8,
     &        VVall,N1X,MPI_REAL8,
     &        0,MPI_COMM_WORLD,ierr)
         IF(myrank.eq.0)then
            CALL OneFFTF(VVall,FCVVall,planF,N1,N1FULL,FS)
         ENDIF
         
         CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
         CALL MPI_Scatter(FCVVall,NSUBX,MPI_COMPLEX16,
     &        FCVV,   NSUBX,MPI_COMPLEX16,
     &        0,MPI_COMM_WORLD,ierr)
c     
         DO i=1,NSUBX           ! normalized Fourier coef-ts
            FCVk(i)=0.5d0*(FCVpr(i)+FCVV(i))
         ENDDO
c     
c     If niter > 0, iterate for a better prediction.
         IF(niter.eq.0)then
            DO i=1,N1X
               V(i)=Vs(i)
               del(i)=dels(i)
               psi(i)=psis(i)
c     Store the history
               IF(jshift.LT.jTwn(i))then
                  IF(jnew(3,i)+jshift.le.jnew(2,i))then
                     DO k=jnew(3,i)+1,jnew(3,i)+jshift
                        Vk(k)=FCVk(i)
                     ENDDO
                  ELSE
                     DO k=jnew(3,i)+1,jnew(2,i)
                        Vk(k)=FCVk(i)
                     ENDDO
                     DO k=jnew(1,i)+1,jnew(1,i)+jshift-
     &                    (jnew(2,i)-jnew(3,i))
                        Vk(k)=FCVk(i)
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
         ELSE
            DO iiter=1,niter
ccccc 
c     Begining of iteration. Not indented for this block.
c     Find better prediction of Slip
               DO i=1,N1X
                  dels(i)=delpr(i)+dt*(Vpr(i)+Vs(i))/2.d0
                  delV(i)=dels(i)-Vpl*(t+dt)
c     State evolution
                  IF(dabs(x(i)-x0fr).le.(0.5d0*xfr))then
                     CALL const_param(param,nparam,fo(i),Vo,
     &                    cca(i),ccb(i),L(i))
                     CALL s_evolution(psis(i),psipr(i),psis(i),dt, ! psis(i) is added as the referred state variable in the second half of the timestep. HIRO 20100903
     &                    Vpr(i),Vs(i),Seff(i),
     &                    param,nparam,Mnum(i))
                  ENDIF
               ENDDO
c     
c     Fourier trasnform of slip difference
               CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
               CALL MPI_Gather(delV,   N1X,MPI_REAL8,
     &              delVall,N1X,MPI_REAL8,
     &              0,MPI_COMM_WORLD,ierr)
               IF(myrank.eq.0)then
                  CALL OneFFTF(delVall,FCdelVall,planF,N1,N1FULL,FS)
               ENDIF
               
               CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
               CALL MPI_Scatter(FCdelVall,NSUBX,MPI_COMPLEX16,
     &              FCdelV,   NSUBX,MPI_COMPLEX16,
     &              0,MPI_COMM_WORLD,ierr)
c     
c     Better guess for the functional term. 
c     Only the "head" portion is updated.
               DO i=1,NSUBX
                  IF(jshift.ge.jTwn(i)) then 
!     Using presummed values of kernel
                     Snh(i)=wwsum(i)*FCVk(i)
                  ELSE
                     Snhss(i)=0.d0
                     ihconv=jnew(3,i)-jnew(1,i)
                     IF(jnew(3,i)+jshift.le.jnew(2,i))then
c     !all needed values fit inbetween jnew(3,i) and jnew(2,i)
c     !only one cycle is needed
                        DO k=jnew(3,i)+1,jnew(3,i)+jshift
                           Snhss(i)=Snhss(i)+ww(k-ihconv)
                        ENDDO
                     ELSE
                        DO k=jnew(3,i)+1,jnew(2,i)
                           Snhss(i)=Snhss(i)+ww(k-ihconv)
                        ENDDO
                        DO k=jnew(1,i)+1,jnew(1,i)+jshift
     &                       -(jnew(2,i)-jnew(3,i))
                           Snhss(i)=Snhss(i)+ww(k+jTwn(i)-ihconv)
                        ENDDO
                     ENDIF
                     Snhss(i)=Snhss(i)*FCVk(i)
                  ENDIF
               ENDDO
c     
c     Calculation of the convolution term using the same "Tail" contribution.
               DO i=1,NSUBX
                  IF(mode.eq.2)then
                     termst=-An(i)*rhoexpr*FCdelV(i)
                  ELSEIF((mode.eq.3).or.(mode.eq.4))then
                     termst=-An(i)*FCdelV(i)
                  ENDIF
                  termdyn=An(i)*(Snhss(i)+Snt(i))*dtelas 
                  FCfunc(i)=termst+termdyn
               ENDDO
c     
c     Inverse Fourier transform to get stress transfer functional
               CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
               CALL MPI_Gather(FCfunc,   NSUBX,MPI_COMPLEX16,
     &              FCfuncall,NSUBX,MPI_COMPLEX16,
     &              0,MPI_COMM_WORLD,ierr)
               IF(myrank.eq.0)then
                  CALL OneFFTB(FCfuncall,funcall,planI,N1,N1FULL)
               ENDIF
               
               CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
               CALL MPI_Scatter(funcall,N1X,MPI_REAL8,
     &              func,   N1X,MPI_REAL8,
     &              0,MPI_COMM_WORLD,ierr)
c     
c     Find a better prediction of velocity
               DO i = 1,N1X
                  IF (dabs(x(i)-x0fr).le.(0.5d0*xfr)) then
                     CALL const_param(param,nparam,fo(i),Vo,
     &                    cca(i),ccb(i),L(i))
                     CALL NRsearch(Vs(i),tau(i),
     &                    tauo(i),func(i),
     &                    Seff(i),Vpl,eta,psis(i),
     &                    param,nparam,Mnum(i),ierror,x(i),
     &                    NRiter,tol)
c     Modifiable: Vs,V3s,tau,tau3,ierror
                  ENDIF
               ENDDO
c     Compute FCVV,FCVV3 using the predicted slip rates.
               DO i =1,N1X
                  VV(i) = Vs(i) - Vpl
               ENDDO
               
c     Fourier transform of slip rate difference
               CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
               CALL MPI_Gather(VV,   N1X,MPI_REAL8,
     &              VVall,N1X,MPI_REAL8,
     &              0,MPI_COMM_WORLD,ierr)
               IF(myrank.eq.0)then
                  CALL OneFFTF(VVall,FCVVall,planF,N1,N1FULL,FS)
               ENDIF
               
               CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
               CALL MPI_Scatter(FCVVall,NSUBX,MPI_COMPLEX16,
     &              FCVV,   NSUBX,MPI_COMPLEX16,
     &              0,MPI_COMM_WORLD,ierr)
c     
               DO i=1,NSUBX     ! normalized Fourier coef-ts
                  FCVk(i)=0.5d0*(FCVpr(i)+FCVV(i))
               ENDDO
ccccc 
c     End of iteration
ccccc 
               IF(iiter.eq.niter)then
                  DO i=1,N1X
                     V(i)=Vs(i)
                     del(i)=dels(i)
                     psi(i)=psis(i)
		  END DO
    		  DO i=1,NSUBX
c     Store the history.
                     IF(jshift.LT.jTwn(i))then
                        IF(jnew(3,i)+jshift.le.jnew(2,i))then
                           DO k=jnew(3,i)+1,jnew(3,i)+jshift
                              Vk(k)=FCVk(i)
                           ENDDO
                        ELSE
                           DO k=jnew(3,i)+1,jnew(2,i)
                              Vk(k)=FCVk(i)
                           ENDDO
                           DO k=jnew(1,i)+1,jnew(1,i)+jshift-
     &                          (jnew(2,i)-jnew(3,i))
                              Vk(k)=FCVk(i)
                           ENDDO
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF

c     Update Time 
         t=t+dt
         ndt=ndt+1

c-------------------------------------------------------------------------
c     Below are output statements until the end of the main loop.
c     These outputs occur while the event is ongoing.
c-------------------------------------------------------------------------
c     Output slip info at points
         DO i=1,noutH
            IF (myrank.eq.Num_procP(i)) then
               WRITE(10+i,430)
     &              ndt,t,
     &              del(P1r(i)),V(P1r(i)),
     &              tau(P1r(i)),psi(P1r(i))
             ENDIF
            
         ENDDO

c     Output initial setup 
         IF (ndt.eq.1) then
            CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
            CALL MPI_Gather(V,   N1X,MPI_REAL8,
     &           Vall,N1X,MPI_REAL8,
     &           0,MPI_COMM_WORLD,ierr)
            CALL MPI_Gather(Seff,   N1X,MPI_REAL8,
     &           Seffall,N1X,MPI_REAL8,
     &           0,MPI_COMM_WORLD,ierr)
            CALL MPI_Gather(fo,   N1X,MPI_REAL8,
     &           foall,N1X,MPI_REAL8,
     &           0,MPI_COMM_WORLD,ierr)
            CALL MPI_Gather(ccb,   N1X,MPI_REAL8,
     &           ccball,N1X,MPI_REAL8,
     &           0,MPI_COMM_WORLD,ierr)
            CALL MPI_Gather(cca,   N1X,MPI_REAL8,
     &           ccaall,N1X,MPI_REAL8,
     &           0,MPI_COMM_WORLD,ierr)
            CALL MPI_Gather(psi,   N1X,MPI_REAL8,
     &           psiall,N1X,MPI_REAL8,
     &           0,MPI_COMM_WORLD,ierr)
            CALL MPI_Gather(L,   N1X,MPI_REAL8,
     &           Lall,N1X,MPI_REAL8,
     &           0,MPI_COMM_WORLD,ierr)
            CALL MPI_Gather(tauo,   N1X,MPI_REAL8,
     &           tauoall,N1X,MPI_REAL8,
     &           0,MPI_COMM_WORLD,ierr)
            IF(myrank.eq.0)then
               DO i = 1,N1
                  WRITE(102,418) xall(i),foall(i),ccball(i),ccaall(i),
     &                 tauoall(i),Seffall(i),Vall(i),Lall(i),psiall(i)
                  
               ENDDO       
            ENDIF
         ENDIF
         
         
c     Find the maximum slip rate and temperature on each node and solve for potency rate in VW region
         Vmaxpr = Vmax
         Vmaxlocal = 0.d0
         DO i = 1,N1X
            IF (dabs(V(i)).GT.Vmaxlocal)THEN ! find Vmax on the fault
               Vmaxlocal = V(i)
               maxVIndex = myrank*N1X+i
            ENDIF
         ENDDO

c     Fill node array for Vmax and Temax
         Vmaxin(1) = Vmaxlocal
         Vmaxin(2) = maxVIndex
         
c     Find the maximum slip rate among all of the nodes and share this
         CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
         CALL MPI_AllREDUCE(Vmaxin,Vmaxout,1,MPI_2DOUBLE_PRECISION,
     &        MPI_MAXLOC,MPI_COMM_WORLD,ierr)
         Vmax = Vmaxout(1)
         VmaxLoc = Vmaxout(2)

c     We will use the main node for outputting
         IF(myrank.eq.0)then    ! 1_vmax
            WRITE(101,408) t,dt,ndt,Vmax,VmaxLoc
         ENDIF
         

         IF (idelevne.eq.0) then ! If we're not currently in an event...check to see if Vmax > threshold for an event

c     Starting event !!!!!
            IF (Vmax.ge.Vevne) then 

c     Set event timing info for different thresholds
               IF (Vmax.ge.Vevne) THEN
                  idelevne = 1  ! set event flag to "true" 
                  Endtin = ndt  ! set event initial timestep to the current timestep 
                  Etin = 
     &                 t - dt*dlog(Vmax/Vevne)/dlog(Vmax/Vmaxpr)
                  Eten = t      ! event finish time (just set it to the start time for now)
               ENDIF
               
               Enum = Enum + 1  ! increment the event number
               tevneb = t       ! set the event beginning time to the current time  
               tevne = tevneinc ! set next output time initially to the increment size ( line_e files ) 
              
c     Call state information to main node for outputting
               CALL MPI_Gather(V,   N1X,MPI_REAL8,
     &              Vall,N1X,MPI_REAL8,
     &              0,MPI_COMM_WORLD,ierr)
               CALL MPI_Gather(del,   N1X,MPI_REAL8,
     &              delall,N1X,MPI_REAL8,
     &              0,MPI_COMM_WORLD,ierr)
               CALL MPI_Gather(tau,   N1X,MPI_REAL8,
     &              tauall,N1X,MPI_REAL8,
     &              0,MPI_COMM_WORLD,ierr)
               CALL MPI_Gather(psi,   N1X,MPI_REAL8,
     &              psiall,N1X,MPI_REAL8,
     &              0,MPI_COMM_WORLD,ierr)
               
               IF(myrank.EQ.0)THEN
c     Create output for slip distribution 
                  IF (Enum.LE.9) THEN ! Store numbering scheme into "ach" for naming of output files
                     WRITE(ach,400) 0,0,0,Enum
                  ELSEIF (Enum.LE.99) THEN
                     WRITE(ach,401) 0,0,Enum
                  ELSEIF (Enum.LE.999) THEN
                     WRITE(ach,402) 0,Enum
                  ELSE
                     WRITE(ach,403) Enum
                  ENDIF
                  
c     Output for individual events (numbered in the 200s) 
c     Use '//ach' with filenames in order to get proper numbering for events 
c     Should be closed at the end of each event 
                  OPEN(200,FILE = trim(wdir)//'/line_e_'//ach) ! point, slip, velocity, stress, state var, temperature, effective normal stress, for each timestep in event 
                  OPEN(201,FILE = trim(wdir)//'/slipbefore_'//ach) ! accumulated slip along the fault right before the event starts 
                  OPEN(202,FILE = trim(wdir)//'/slipafter_'//ach) ! accumulated slip along the fault right after the event ends 

c     Output line_e ( first time step of events)
                  CALL out_vsx(200,N1, 
     &                 delall,Vall,tauall,psiall,
     &                 ndt,t,hcell,outStepSize)

c     Slip BEFORE an individual event
                  CALL out_vsx(201,N1, 
     &                 delall,Vall,tauall,psiall,
     &                 ndt,t,hcell,outStepSize)
               ENDIF            ! (myrank.eq.0)               
c     Otherwise we're not in an event (Vmax < Vevne)
            ELSE   
c     Check to see if we save quasi-static time step information
               IF(t.GT.tvsx)then                  
c     Call everything to main node for outputting
                  CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
                  CALL MPI_Gather(V,   N1X,MPI_REAL8,
     &                 Vall,N1X,MPI_REAL8,
     &                 0,MPI_COMM_WORLD,ierr)
                  CALL MPI_Gather(del,   N1X,MPI_REAL8,
     &              delall,N1X,MPI_REAL8,
     &                 0,MPI_COMM_WORLD,ierr)
                  CALL MPI_Gather(tau,   N1X,MPI_REAL8,
     &                 tauall,N1X,MPI_REAL8,
     &                 0,MPI_COMM_WORLD,ierr)
                  CALL MPI_Gather(psi,   N1X,MPI_REAL8,
     &                 psiall,N1X,MPI_REAL8,
     &                 0,MPI_COMM_WORLD,ierr)
                  IF(myrank.eq.0)THEN
                     CALL out_vsx(100,N1,
     &                    delall,Vall,tauall,psiall,
     &                    ndt,t,hcell,outStepSize)

                  ENDIF
                  tvsx = tvsx + tvsxinc
               ENDIF               
            ENDIF               ! Vmax > Vevne
c-------------------------------------------------------------------------
c---------End initialization of event parameters--------------------------
c-------------------------------------------------------------------------
c     Event is already ongoing, we keep track of information until an event end time is determined and
c     the waves leave the system in case the velocity is oscillatory or some other rupture is triggered
         ELSEIF((idelevne.EQ.1).or.
     &           (t.le.(Eten+dble(jTw)*dtelas)))THEN
            
c     If the max slip rate falls beneath threshold then event it tentatively ending
            IF((idelevne.EQ.1).and.(Vmax.LT.Vevne))THEN
               idelevne = 2
               Eten = t
               Endten = ndt
            ELSEIF((idelevne.EQ.2).and.(Vmax.ge.Vevne))THEN
               idelevne = 1
            ENDIF
            
                                    
c     Check to see if we need to gather things for the line_e file
            IF ((t-tevneb).ge.tevne)then
               CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
               CALL MPI_Gather(V,   N1X,MPI_REAL8,
     &              Vall,N1X,MPI_REAL8,
     &              0,MPI_COMM_WORLD,ierr)
               CALL MPI_Gather(del,   N1X,MPI_REAL8,
     &              delall,N1X,MPI_REAL8,
     &              0,MPI_COMM_WORLD,ierr)
               CALL MPI_Gather(tau,   N1X,MPI_REAL8,
     &              tauall,N1X,MPI_REAL8,
     &              0,MPI_COMM_WORLD,ierr)
               CALL MPI_Gather(psi,   N1X,MPI_REAL8,
     &              psiall,N1X,MPI_REAL8,
     &              0,MPI_COMM_WORLD,ierr)

               IF(myrank.eq.0)then
                  CALL out_vsx(200,N1,
     &                 delall,Vall,tauall,psiall,
     &                 ndt,t,hcell,outStepSize) 
                  
               ENDIF
               tevne = tevne + tevneinc ! update the event output time
            ENDIF               ! IF ((t-tevneb).ge.tevne) then 
c            
c-------------------------------------------------------------------------
c------------End of an Event----------------------------------------------
c-------------------------------------------------------------------------
         ELSE      
c     Set event flag to zero for end of event
            idelevne = 0 

            CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
c     Gather arrays for line_e file (event end)
            CALL MPI_Gather(V,   N1X,MPI_REAL8,
     &           Vall,N1X,MPI_REAL8,
     &           0,MPI_COMM_WORLD,ierr)
            CALL MPI_Gather(del,   N1X,MPI_REAL8,
     &           delall,N1X,MPI_REAL8,
     &           0,MPI_COMM_WORLD,ierr)
            CALL MPI_Gather(tau,   N1X,MPI_REAL8,
     &           tauall,N1X,MPI_REAL8,
     &           0,MPI_COMM_WORLD,ierr)
            CALL MPI_Gather(psi,   N1X,MPI_REAL8,
     &           psiall,N1X,MPI_REAL8,
     &           0,MPI_COMM_WORLD,ierr)
            
            IF(myrank.EQ.0)THEN                     
c     Catalogue of event information
               WRITE(106,404) Enum,Endtin,Endten,Etin,Eten 
C     Slip AFTER an individual event
               CALL out_vsx(202,N1,
     &              delall,Vall,tauall,psiall,
     &              ndt,t,hcell,outStepSize)

c     Close event files
               CLOSE(200)       ! line_e
               CLOSE(201)       ! slipbefore_
               CLOSE(202)       ! slipafter_
            ENDIF               ! (myrank.eq.0)
c-------------------------------------------------------------------
         ENDIF                  !IF (idelevne.eq.0) ELSE  (event over now)

c     Store the information in this time step and determine next time step
         dtpr=dt
         jshiftpr=jshift
         DO i=1,N1X
            Vpr(i)=V(i)
            taupr(i)=tau(i)
            delpr(i)=del(i)
            psipr(i)=psi(i)
         ENDDO
     
         DO i=1,NSUBX
            FCVpr(i)=FCVV(i)
         ENDDO
         
c     Update frictional properties and estimate dLmax
         CALL Set_friction(dfdlV,dfssdlV,
     &        V,Seff,nparam,Mnum,psi,cca,ccb,L,
     &        fo,Vo,N1X,x,xfr,x0fr)
         CALL Set_dLmax(dLmax,dfdlV,dfssdlV,hcell,L,Seff,
     &        N1X,xmu,x,xfr,nu,mode,Sfac,x0fr)         

c     Determine the next time step
         CALL dtevol(N1X,x,jTw,dLmax,dtmax,dtmin,dtincf,dtlocal,
     &        dt,V,xfr,NPU,x0fr)

         CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
         CALL MPI_ALLREDUCE(dtlocal,dt,1,MPI_REAL8,MPI_MIN,
     &        MPI_COMM_WORLD,ierr)
         
c     dt must be INTEGER*dtmax if dt < jTw*dtmin .
         IF(dt.GT.(dble(jTw)*dtmin))THEN
            jshift = jTw
         ELSE
            jshift = dt/dtmin
            IF(jshift.eq.0)THEN
               jshift = 1
               dt = dtmin
            ELSE
               dt = dtmin*jshift
            ENDIF
         ENDIF

 2000    cyclecont= (t.le.tmax).and.(ndt.le.ndtmax)
     &        .and.(Enum.le.Enumax).and.(ierror.eq.0)
      ENDDO                     !End of Main(Time) Cycle
      

c     close remaining output files
      CLOSE(100)                ! line_c
      CLOSE(101)                ! vmax
      CLOSE(102)                ! initialSetup 
      CLOSE(103)                ! ErReport
      CLOSE(104)                ! HpointLocations  
      CLOSE(105)                ! settings
      CLOSE(106)                ! Ecatalogue

c     
c     Error report, if any
      IF(ierror.eq.2) THEN
         WRITE(103,*) 'N-R search!'
         WRITE(103,*) myrank
      ENDIF
      
c---------FORMATTING FOR OUTPUTS------------------------------------------
c  Key:
c 
c   w - the number of positions to be used (should at least be 6 greater than d, need to have space for 0.E+##)
c   d - number of digits to the right of the decimal point
c
c   Aw - text string
c   Dw - double precision numbers, exponent notation
c   Ew.d - real numbers, exponent notation
c   Fw - real numbers, fixed point format
c   Iw - INTEGER
c   nX - horizontal skip (space)
c   / - vertical skip (newline)
    
c  Please using numbering beginning in the 400s for easier searching (won't accidentally match a file number) 
 400  format (I1,I1,I1,I1)                                                 ! used for formatting "ach"
 401  format (I1,I1,I2)                                                    ! used for formatting "ach"
 402  format (I1, I3)                                                      ! used for formatting "ach"
 403  format (I4)                                                          ! used for formatting "ach"
 408  format (1X,E22.14,1X,E15.7,1X,I8,1X,E15.7,1X,I10)                    ! vmax
 404  format (1X,I3,2(2X,I10),2(2X,E22.14))                                ! EventCatalogue
 418  format (9(1X,E15.7))						   ! initialSetup
 430  format (I9,1X,5(1X,E22.14))					   ! H### files
 517  format (I7,1X,E15.7)                                                 ! 1_HpointLocations
 518  format (I8,5(1X,E15.7))                                              ! 1_settings

      STOP
      CALL MPI_FINALIZE(ierr)
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccmccccccccccc
c                                                                      c
c eeeee  n   n  dddd                                                   c
c e      nn  n  d   d                                                  c
c eeeee  n n n  d   d                                                  c
c e      n  nn  d   d                                                  c
c eeeee  n   n  dddd                                                   c
c                                                                      c
c                ooo   fffff                                           c
c               o   o  f                                               c
c               o   o  ffff                                            c
c               o   o  f                                               c
c                ooo   f                                               c
c                                                                      c
c                      m   m    a     iii   n   n   !!                 c
c                      mm mm   a a     i    nn  n   !!                 c
c                      m m m  aaaaa    i    n n n   !                  c
c                      m   m  a   a    i    n  nn                      c
c                      m   m  a   a   iii   n   n   !                  c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c-----------------------------------------------------------------------                            
c------------- FFT Functions to incorporate Mirroring ------------------                              
c--------------- Added by Valere Lambert 03-11-18 ----------------------    
c-----------------------------------------------------------------------        
c-----------------------------------------------------------------------
c     subroutine OneFFTF computes Fourier transform accounting for mirror       
c                                                                           
c     INPUT:                                                                              
c     FS: Symmetry for free surface 
c     1 = no free surface, 2 = free surface mirror         
c-----------------------------------------------------------------------
      Subroutine OneFFTF(Input,Output,planF,N1,N1FULL,FS)
      IMPLICIT NONE
      INTEGER N1,N1FULL,i,FS
      INTEGER*8 planF
      REAL*8 Input(N1)
      complex*16 Output(N1FULL)

      complex*16 inFFTW(N1FULL),outFFTW(N1FULL)
      complex*16 temp(N1)

c  copy input to local buffer                                                        
      IF (FS.eq.1) THEN
         DO i=1,N1
            inFFTW(i)=dcmplx(Input(i))
         ENDDO
      ELSEIF (FS.eq.2) THEN     ! Mirror image                           
         DO i=1,N1
            temp(i)=dcmplx(Input(i))
            inFFTW(i)=temp(i)
            inFFTW(N1FULL+1-i)=temp(i)
         ENDDO
      ENDIF

c  forward Fourier transform                                           
      call dfftw_execute_dft(planF,inFFTW,outFFTW)

      DO i=1,N1FULL
         Output(i)=outFFTW(i)/N1FULL
      ENDDO

      RETURN
      END

c--------------------------------------------------------------------
c     subroutine OneFFTB performs the inverse Fourier transform  
c--------------------------------------------------------------------
      Subroutine OneFFTB(Input,Output,planI,N1,N1FULL)
      IMPLICIT NONE
      INTEGER N1,N1FULL,i
      INTEGER*8 planI
      REAL*8 Output(N1)
      complex*16 Input(N1FULL),outIFFTW(N1FULL)
      call dfftw_execute_dft(planI,Input,outIFFTW)
      DO i=1,N1
         Output(i) = dble(outIFFTW(i))
      ENDDO

      RETURN
      END

c-------------------------------------------------------------------------
c     Parse data from input line
c-------------------------------------------------------------------------
      SUBROUTINE getdata(unit,line)
      IMPLICIT NONE
c     First implemented in Potsdam, Feb, 1999
c     Last modified: Valere Lambert, 2018
c     This subroutine extracts data from an input line
      INTEGER unit
      CHARACTER(200) :: line
      CHARACTER(1) :: char
      INTEGER i
      LOGICAL linecont
c     this subroutine reads over all comment lines starting with "#".
      char='#'
      linecont = .TRUE.
      DO WHILE(linecont)
         IF(char.EQ.'#')THEN
            read(unit,'(a)')line
            i=1
            char=line(1:1)
         ELSEIF(char.EQ.' ')THEN
            i=i+1
            char=line(i:i)
         ELSE
            linecont = .FALSE.
         ENDIF
      ENDDO
      RETURN
      END
c-------------------------------------------------------------------------
c     Read output directory from submission file
c-----------------------------------------------------------------------      
      SUBROUTINE readinput(wdir,myrank)
c     Last modified: Valere Lambert, 2018
c     Read directories from input
c     UNIT = 5 for input from keyboard
      INCLUDE 'mpif.h'
      INTEGER :: myrank,position,ierr
      CHARACTER(200) :: dataline
      INTEGER, PARAMETER :: psize=512
      CHARACTER, DIMENSION(psize) :: packed
      CHARACTER(200) :: wdir

      position = 0
      IF(0.EQ.myrank)THEN
c     Get output directory
         CALL getdata(5,dataline)
         READ (dataline,'(a)') wdir
      ENDIF

      CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(wdir,200,MPI_CHARACTER,
     &     0,MPI_COMM_WORLD,ierr)
      RETURN
      END 

c-----------------------------------------------------------------------
c-----------------------ORIGINAL SUBROUTINES----------------------------
c-----------------------------------------------------------------------
      Subroutine out_vsx(nfile,N1,
     &              delall,Vall,tauall,psiall,
     &              ndt,t,hcell,stepsize)
      ! Output physical variables along fault

      IMPLICIT NONE
      INTEGER i,nfile,N1,ndt,stepsize
      REAL*8 Vall(N1),tauall(N1),delall(N1),psiall(N1),
     &     Lall(N1),t,xcord,hcell
c
      DO i = 1, N1, stepsize                                 			! go through all points in steps of "stepsize" in order to reduce file size of "line_*" files
         xcord = (dfloat(i)-N1/2)*hcell						! find the xcoordinate of the points in space, given their indices from 1,N1
            WRITE(nfile,337) xcord,delall(i),Vall(i),				! output xcoord and it's physical variables
     &           tauall(i),psiall(i)
      ENDDO

      WRITE(nfile,365) '# **',' above is for ndt, t = ',ndt,t			! record the time and timestep below all of the points
      WRITE(nfile,'()')
      endfile nfile
      backspace nfile
c
337   format (1X,E13.5,4(1X,E15.7))						! formatting for output (see above formatting table for help)
365   format (T1,A,A,I8,1X,E25.14)
      RETURN
      END
c
c--------------------------------------------------------------------
      Subroutine dtevol(N1X,x,jTw,dLmax,dtmax,dtmin,dtincf,
     &     dtlocal,dt,V,xfr,NPU,x0fr)
      ! Calculate local timestep, used to determine global
      ! optimal/minimum timestep
c     Modifiable: dtlocal

      IMPLICIT NONE
      INTEGER i,N1X,jTw,temp,NPU
      REAL*8 dtmax,dtmin,dtincf,dtcell,dtnx,Vmax,dtin,dt,xfr,x0fr
      REAL*8 V(N1X),dLmax(N1X),dtlocal(NPU),x(N1X)
c
      dtcell=dtmax								! Set the timestep for this cell to the max
      DO i=1,N1X								! Loop through all points on this node
         IF(dabs(x(i)-x0fr).le.(0.5d0*xfr))THEN					! If the point is in the frictional region, we evaluate it
            dtin=dLmax(i)/dabs(V(i))					!Calculates timestep based on dLmax and current velocity at this point
            IF(dtin.LT.dtcell) THEN						! If calculated timestep is smaller than dtcell, make it dtcell
               dtcell=dtin
            ENDIF
         ENDIF
      ENDDO
c
      IF(dtcell.GT.dtmax) THEN							! Checks to make sure dtcell<=dtmax and dtcell>=dtmin
         dtnx=dtmax
      ELSEIF(dtcell.LT.dtmin) THEN
         dtnx=dtmin
      ELSE
         dtnx=dtcell
      ENDIF
c
      IF(dtnx.GT.(dt*dtincf))   dtnx=dt*dtincf					! We limit the timestep to increasing only as fast as dt*dtincf
      dtlocal(1)=dtnx								! Set the timestep for points on this node 
c										! (will be evaluated against other nodes to determine overall timestep)
c
      RETURN
      END		
c
c-----------------------------------------------------------------------    
      SUBROUTINE NRsearch(V,tau,
     &           tauo,func,
     &           Seff,Vpl,eta,psi,
     &           param,nparam,Mnum,ierror,locx,
     &           NRiter,tol)
      ! This subroutine solves an equation, (strength) = (stress)
      ! using Newton-Raphson method

      IMPLICIT NONE
c
      ! IO parameters
      ! INTEGERs
      ! nparam : Number of the frictional parameters
      ! Mnum   : Material number
      ! ierror : Error flag!
      !
      ! REAL*8
      ! V      : Slip rate [to be updated]
      ! tau    : Shear stress [to be updated]
      ! tauo   : Loading in
      ! func   : Stress transfer
      ! Seff   : Effective normal stress
      ! Vpl    : Plate rate (Elasticity sees difference from it)
      ! eta    : Radiation damping term, mu/2cs
      ! psi    : State variable, f-(direct term)
      ! param  : An array containing the frictional parameters
      INTEGER nparam,Mnum,ierror,NRiter
      REAL*8 V,tau,tauo,func,Seff,
     &     Vpl,eta,psi,param(nparam),locx,eta0
c
      ! Other parameters
      ! j      : Number of updates in N-R search
      ! V      : Absolute value of the slip rate
      ! tau    : Absolute value of the shear stress
      ! tpotential: tau at V = Vpl
      ! ratio  : Tangent of the frictional resistance (upward plunge)
      ! sc_str : Stress scale
      ! sc_vel : Velocity scale
      ! tol    : Tolerance
      ! dntau  : Normalized update in tau
      ! dnV    : Normalized update in V
      ! F1     : tau-(tauo+f-eta*V)
      ! F2     : tau-Seff*f
      ! F11,F12,F21,F22
      !        : Partial drivative of F1 and F2 w/r.to tau (1) and V (2)
      ! det    : Determinant of a matrix, (Fij)
      INTEGER j
      REAL*8 Vab,tauab,tpotential,ratio,sc_str,sc_vel,tol,dntau,dnV,
     &     F1,F2,F11,F12,F21,F22,det,dV
c
      ! Externally defined functions
      REAL*8 fric_f,fric_dfdV
c
      Vab=dabs(V)
      tauab=dabs(tau)
      tpotential=dabs(tauo+func+eta*Vpl)
c
      IF(Seff.le.0.d0)THEN							! Check to see if Seff is non-negative
         V = (tauo+func+eta*Vpl)/eta						! If Seff < 0 (non physical), THEN prescribe V...
         tau=0.d0								! and set tau = 0
         RETURN
      ENDIF
c
ccc This is a 2D Newton-Raphson scheme 
      ! Tolerance and initial scales
      IF(tau.eq.0.d0)then
         sc_str=1.d0
      ELSE
         sc_str = tauab
      ENDIF
      sc_vel = Vab
      eta0 = eta
c
      ! Initialize the update 
      dntau = tol*10.d0
      dnV   = tol*10.d0 
c
      DO j=1,NRiter
         ! Newton-Raphson search
         ! These functions are reduced to zero.
         F1 = (tauab-tpotential+eta0*Vab)/sc_str
         F2 = (tauab-Seff*fric_f(Vab,psi,param,nparam,Mnum))
     &        /sc_str
c
         IF ((dabs(dntau).le.tol).and.(dabs(dnV).le.tol).and.			! if the updates to tau and V are below tolerance...
     &       (dabs(F1).le.tol).and.(dabs(F2).le.tol)) THEN			! and the functions F1 and F2 have been minimized to below tolerance...
            IF(Vab.ge.1000.d0)THEN						! if the solution we got has a velocity > 1000 m/s...
              ierror = 2							! we consider it unphysical and throw an error
               RETURN
            ENDIF
            V=Vab*(tauo+func+eta0*Vpl)/tpotential				! ...otherwise set V and tau
            tau=tauab*(tauo+func+eta0*Vpl)/tpotential
            RETURN
         ENDIF
c
         ! Jacobian, indep. variables are 1:tau/sc_str 2:V/sc_vel
         F11 = 1.d0
         F21 = 1.d0
         F12 = eta0*sc_vel/sc_str
         F22 = -Seff*fric_dfdV(Vab,psi,param,nparam,Mnum)
     &        *sc_vel/sc_str
         det = F11*F22-F12*F21
         ! Shift in normalized tau and V
         dntau = -( F22*F1 - F12*F2)/det
         dnV   = -(-F21*F1 + F11*F2)/det 
         IF(dnV.le.-1.d0)THEN							! if updates in V are too large...
            dntau = dntau * (-0.9d0/dnV)
            dnV = -0.9d0
         ENDIF
         IF(dntau.le.-1.d0)THEN							! if updates in tau are too large...
            dnV = dnV * (-0.9d0/dntau)
            dntau = -0.9d0
         ENDIF
         ! Update
         tauab = tauab + dntau*sc_str						! do the updates
         Vab = Vab + dnV*sc_vel

         ! Choose a proper scale. 
         ! V may be changing by orders in this step!
         sc_str = tauab
         sc_vel = Vab
      ENDDO									! do j=1,1000
c
      ierror = 2								! if the NR-search doesn't converge after 1000 iterations, output an error
      RETURN
      END
c
c-----------------------------------------------------------------------
      Subroutine kern2(N1,NSUBX,Nlnew,jbeta,delay,jobname,rho,
     &     dtelas,bej1sx,w,bessj0,xker2,integral,jTwn,jnew,ww2,ww2sum,
     &     cs,xLam,myrank)
      ! kernel calculation for mode 2
c     Modifiable: ww2,ww2sum
c     The values are computed for the mid-point rule integration 
c     (requires delay = 0.5)

      IMPLICIT NONE
      CHARACTER*80 jobname
      INTEGER N1,NSUBX,Nlnew,jbeta,i,k1,k,ihelp,myrank,effi
      REAL*8 delay,s1,s2,ww2base,pi,rho,b5,dtelas,q,xLam
      REAL*8 bej1sx,w,bessj0,xker2,integral,cs
      PARAMETER (pi=3.141592653589793d0)
      INTEGER jnew(3,NSUBX),jTwn(NSUBX)
      REAL*8 ww2(Nlnew),ww2sum(NSUBX)
c
c     ! to save kernel values for some modes:
      DO i=1,NSUBX
         effi=myrank*NSUBX+i
         IF (effi.le.(N1/2+1)) THEN
            k1 = effi - 1
         ELSE
            k1 = effi - 1 - N1
         ENDIF
         q=2.d0*pi*dabs(k1/xLam)   					! modal frequency 
         ihelp = jTwn(i)
cc!  Assign kernel value for current time, delayed by the delay:
         ww2base = 2.d0*(1.d0-1.d0/rho)
         s1 = 0.d0
         s2 = q*cs*dtelas*(0.d0+delay)
c     b5 = wwbase - 0.5d0*(xker2(s1,rho)+xker2(s2,rho))*An*delay
         b5 = ww2base - integral(xker2,10,s1,s2,rho)
         k = 1
         ww2(jnew(1,i)+k)=b5    
c
cc!  Assign all the other kernel values:
         DO k=2,ihelp
c     ! delay all modes for mid-point scheme (with delay = 0.5);
            s1 = q*cs*dtelas*(dfloat(k-2)+delay)
            s2 = q*cs*dtelas*(dfloat(k-1)+delay)
c     b5 = b5 - 0.5d0*(xker2(s1,rho)+xker2(s2,rho))*An
            b5 = b5 - integral(xker2,10,s1,s2,rho)
            ww2(jnew(1,i)+k)=b5    
         ENDDO
      ENDDO
c     ! Presummed kernel for long time steps:
      DO i=1,NSUBX
         ww2sum(i)=0.d0
         DO k=jnew(1,i)+1,jnew(2,i)
            ww2sum(i)=ww2sum(i)+ww2(k)
         ENDDO
      ENDDO
      RETURN
      END

c---------------------------------------------------------------------    
      Subroutine kern3(N1,NSUBX,Nlnew,jbeta,delay,jobname,
     &     jTwn,jnew,ww3,ww3sum,cs,xLam,dtelas,myrank)
      ! kernel calculation for mode 3
c     Modifiable: ww3,ww3sum

      IMPLICIT NONE
      CHARACTER*80 jobname
      INTEGER N1,NSUBX,jTw,Nlnew,jbeta,i,k1,k,ihelp,myrank,effi
      REAL*8 delay,slimit,pi,q,cs,xLam,dtelas
      PARAMETER (pi=3.141592653589793d0)
      INTEGER jnew(3,NSUBX),jTwn(NSUBX)
      REAL*8 ww3(Nlnew),ww3sum(NSUBX)
c
      REAL*8 AA0,AA1,AA2,AA3,AA4,AA5,AA6,AA7,AA8,
     &     AA9,AA10,AA11,AA12,AA13,AA14,
     &     BB0,BB1,BB2,BB3,BB4,BB5,BB6,BB7,
     &     CC1,CC2,CC3,CC4,CC5,CC6,CC7,
     &     b0,b1,b2,b3,b4,b5
c
      DATA AA0,AA1,AA2,AA3,AA4,AA5,AA6,AA7,AA8,AA9,AA10,
     &     AA11,AA12,AA13,AA14
     &     /4.D0,-1.066666666666666D1,1.706666666666666D1,
     &     -1.625396825396825D1,1.011358024691358D1,
     &     -4.413198653198653D0,1.422569529236196D0,
     &     -3.522553120013437D-1,6.906966901987132D-2,
     &     -1.098652045228363D-2,1.445845115798365D-3,
     &     -1.600144792187913D-4,1.509880214167057D-5,
     &     -1.229043723375708D-6,8.718339712615857D-8/
      DATA BB0,BB1,BB2,BB3,BB4,BB5,BB6,BB7
     &     /1.25D-1,-5.26410175D-3,1.273027121D-3,-6.73431837D-4,
     &     5.005321894D-4,-3.452926821D-4,1.591757834D-4,
     &     -3.364922862D-5/
      DATA CC1,CC2,CC3,CC4,CC5,CC6,CC7
     &     /1.757799D-2,-2.274954641D-3,8.472097978D-4,-5.141838574D-4,
     &     3.188374766D-4,-1.387793491D-4,2.835472838D-5/
c
c     ! to save kernel values for some modes:
      DO i=1,NSUBX     ! compute kernel for modes n and (-n) simultaneously
         effi=myrank*NSUBX+i
         IF (effi.le.(N1/2+1)) THEN
            k1 = effi - 1
         ELSE
            k1 = effi - 1 - N1
         ENDIF
         q=2.d0*pi*dabs(k1/xLam)
         ihelp=jTwn(i)
         DO k=1,ihelp
c        ! delay all modes for mid-point scheme (with delay = 0.5);
            slimit=q*cs*dtelas*(dfloat(k-1)+delay)
            IF (slimit.LT.8.d0) THEN
               b0=slimit/8.d0
               b1=b0*b0 
               b2=b1*(AA10+b1*(AA11+b1*(AA12+b1*(AA13+b1*AA14))))
               b2=b1*(AA5+b1*(AA6+b1*(AA7+b1*(AA8+b1*(AA9+b2)))))
               b2=b0*(AA0+b1*(AA1+b1*(AA2+b1*(AA3+b1*(AA4+b2)))))
               b5=1.d0-b2
            ELSE
               b0=8.d0/slimit
               b1=b0*b0
               b2=slimit-0.785398163d0
               b3=b0*(BB0+b1*(BB1+b1*(BB2+b1*(BB3+b1*(BB4
     &              +b1*(BB5+b1*(BB6+b1*BB7)))))))
               b4=b1*(CC1+b1*(CC2+b1*(CC3+b1*(CC4
     &              +b1*(CC5+b1*(CC6+b1*CC7))))))
               b5=DSQRT(0.636619772d0/slimit)*(DCOS(b2)*b3+DSIN(b2)
     &              *b4)
            ENDIF
c
            ww3(jnew(1,i)+k)=b5
         ENDDO
      ENDDO
c     ! Presummed kernel for long time steps:
      DO i=1,NSUBX
         ww3sum(i)=0.d0
         DO k=jnew(1,i)+1,jnew(2,i)
            ww3sum(i)=ww3sum(i)+ww3(k)
         ENDDO
      ENDDO
      RETURN
      END
c
c---------------------------------------------------------------------
      Subroutine kern4(N1,NSUBX,Nlnew,jbeta,delay,jobname,
     &     jTwn,jnew,ww4,ww4sum,cs,nu,xLam,dtelas,Hseis,myrank)
      ! kernel calculation for crustal plane model
c     Modifiable: ww4,ww4sum

      IMPLICIT NONE
      CHARACTER*80 jobname
      INTEGER N1,NSUBX,jTw,Nlnew,jbeta,i,k1,k,ihelp,myrank,effi
      REAL*8 delay,slimit,pi,q,cs,xLam,dtelas,Hseis,nu
      PARAMETER (pi=3.141592653589793d0)
      INTEGER jnew(3,NSUBX),jTwn(NSUBX)
      REAL*8 ww4(Nlnew),ww4sum(NSUBX)
c
      REAL*8 AA0,AA1,AA2,AA3,AA4,AA5,AA6,AA7,AA8,
     &     AA9,AA10,AA11,AA12,AA13,AA14,
     &     BB0,BB1,BB2,BB3,BB4,BB5,BB6,BB7,
     &     CC1,CC2,CC3,CC4,CC5,CC6,CC7,
     &     b0,b1,b2,b3,b4,b5
c
      DATA AA0,AA1,AA2,AA3,AA4,AA5,AA6,AA7,AA8,AA9,AA10,
     &     AA11,AA12,AA13,AA14
     &     /4.D0,-1.066666666666666D1,1.706666666666666D1,
     &     -1.625396825396825D1,1.011358024691358D1,
     &     -4.413198653198653D0,1.422569529236196D0,
     &     -3.522553120013437D-1,6.906966901987132D-2,
     &     -1.098652045228363D-2,1.445845115798365D-3,
     &     -1.600144792187913D-4,1.509880214167057D-5,
     &     -1.229043723375708D-6,8.718339712615857D-8/
      DATA BB0,BB1,BB2,BB3,BB4,BB5,BB6,BB7
     &     /1.25D-1,-5.26410175D-3,1.273027121D-3,-6.73431837D-4,
     &     5.005321894D-4,-3.452926821D-4,1.591757834D-4,
     &     -3.364922862D-5/
      DATA CC1,CC2,CC3,CC4,CC5,CC6,CC7
     &     /1.757799D-2,-2.274954641D-3,8.472097978D-4,-5.141838574D-4,
     &     3.188374766D-4,-1.387793491D-4,2.835472838D-5/
c
c     ! to save kernel values for some modes:
      DO i=1,NSUBX     ! compute kernel for modes n and (-n) simultaneously
         effi=myrank*NSUBX+i
         IF (effi.le.(N1/2+1)) THEN
            k1 = effi - 1
         ELSE
            k1 = effi - 1 - N1
         ENDIF
         q=dsqrt((2.d0*pi*dabs(k1/xLam)/(1.d0-nu))**2
     &                            +(4.d0/pi/Hseis)**2)
         ihelp=jTwn(i)
         DO k=1,ihelp
c        ! delay all modes for mid-point scheme (with delay = 0.5);
            slimit=q*cs*dtelas*(dfloat(k-1)+delay)
            IF (slimit.LT.8.d0) THEN
               b0=slimit/8.d0
               b1=b0*b0 
               b2=b1*(AA10+b1*(AA11+b1*(AA12+b1*(AA13+b1*AA14))))
               b2=b1*(AA5+b1*(AA6+b1*(AA7+b1*(AA8+b1*(AA9+b2)))))
               b2=b0*(AA0+b1*(AA1+b1*(AA2+b1*(AA3+b1*(AA4+b2)))))
               b5=1.d0-b2
            ELSE
               b0=8.d0/slimit
               b1=b0*b0
               b2=slimit-0.785398163d0
               b3=b0*(BB0+b1*(BB1+b1*(BB2+b1*(BB3+b1*(BB4
     &              +b1*(BB5+b1*(BB6+b1*BB7)))))))
               b4=b1*(CC1+b1*(CC2+b1*(CC3+b1*(CC4
     &              +b1*(CC5+b1*(CC6+b1*CC7))))))
               b5=DSQRT(0.636619772d0/slimit)*(DCOS(b2)*b3+DSIN(b2)
     &              *b4)
            ENDIF
c
            ww4(jnew(1,i)+k)=b5
         ENDDO
      ENDDO
c     ! Presummed kernel for long time steps:
      DO i=1,NSUBX
         ww4sum(i)=0.d0
         DO k=jnew(1,i)+1,jnew(2,i)
            ww4sum(i)=ww4sum(i)+ww4(k)
         ENDDO
      ENDDO
      RETURN
      END
c
c---------------------------------------------------------------
      SUBROUTINE conv_storage(N1,NSUBX,jTw,Nlnew,jobname,
     &     jnew,jTwn,ierror,myrank,wdir)
      ! Check if array is suitable for specified convolution window
c     Modifiable: jnew,jTwn,ierror
c
      IMPLICIT NONE
      CHARACTER(200) :: wdir
      INTEGER N1,NSUBX,jTw,Nlnew,i,k,ierror,myrank,effi
      INTEGER jnew(3,NSUBX),jTwn(NSUBX)
      REAL*8 pi,ratio
      CHARACTER*80 jobname
      PARAMETER (pi=3.141592653589793d0)
c
      DO i=1,NSUBX
         jTwn(i) = jTw
      ENDDO
c
      DO i=1,NSUBX
         IF(i.eq.1) THEN
            jnew(1,1)=0
            jnew(2,1)=jTwn(1)
         ELSE
            jnew(1,i) = jnew(2,i-1)
            jnew(2,i) = jnew(2,i-1)+jTwn(i)
         ENDIF
      ENDDO
c
      DO i=1,NSUBX
            jnew(3,i) = jnew(1,i)
c     in the beginning the most current 
c     value will be the first one in storage
      ENDDO
c
      IF(myrank.eq.0) THEN
         IF (jnew(2,NSUBX).GT.Nlnew) THEN ! if the Nlnew is too small
            OPEN (1, FILE = trim(wdir)//'conv_ERROR'//jobname)
            WRITE(1,*)'   ERROR: Nlnew should be'
            WRITE(*,*)'   ERROR: Nlnew should be'
            WRITE(1,*)jnew(2,NSUBX)
            WRITE(*,*)jnew(2,NSUBX)
            WRITE(1,*)'   Right now Nlnew is'
            WRITE(*,*)'   Right now Nlnew is'
            WRITE(1,*)Nlnew
            WRITE(*,*)Nlnew
            WRITE(1,*)'Aborting the run...'
            WRITE(*,*)'Aborting the run...'
            CLOSE(1)
            ierror = 1
         ENDIF
         IF (jnew(2,NSUBX).LT.Nlnew) THEN ! if the Nlnew is too large
            OPEN (1, FILE = trim(wdir)//'conv_WARNING'//jobname)
            WRITE(1,*)'   Warning: could save some memory using Nlnew'
            WRITE(*,*)'   Warning: could save some memory using Nlnew'
            WRITE(1,*)jnew(2,NSUBX)
            WRITE(*,*)jnew(2,NSUBX)
            WRITE(1,*)'   Right now Nlnew is'
            WRITE(*,*)'   Right now Nlnew is'
            WRITE(1,*)Nlnew
            WRITE(*,*)Nlnew
            CLOSE(1)
         ENDIF
      ENDIF
      RETURN
      END
c
c-----------------------------------------------------------------------
      Subroutine Set_friction(dfdlV,dfssdlV,
     &     V,Seff,nparam,Mnum,psi,cca,ccb,L,
     &     fo,Vo,N1X,x,xfr,x0fr)
      ! This subroutine is called near the end of the time loop, and 
      ! updates the frictional constitutive parameters refering 
      ! updated value of V and T.
      IMPLICIT NONE
c
      ! IO parameters
      !    REAL*8
      ! V     : Slip rate
      ! cca   : a-value
      ! ccb   : b-value
      ! L     : State evolution distance
      ! fo    : Reference friction coefficient [scalar]
      ! Vo    : Reference slip rate
      ! x     : point along domain 
      ! xfr   : size of frictional domain (meters)
      !    INTEGER
      ! Mnum  : Material number
      ! N1X : Number of grids in the x-direction
      INTEGER N1X,Mnum(N1X),nparam
      REAL*8 V(N1X),psi(N1X),
     &     x(N1X),xfr,x0fr,
     &     dfdlV(N1X),dfssdlV(N1X),L(N1X),cca(N1X),ccb(N1X),
     &     fo(N1X),Vo,Seff(N1X)
c
      ! Externally defined functions
      REAL*8 fric_dfdV,fric_dfssdV,fric_f,fric_fss

      ! Other parameters
      ! i     : Index in the x-direction
      ! V     : Absolute value of the slip rate
      ! cca_b : "a-b"-value
      INTEGER i
      REAL*8 param(nparam)
c
      DO i = 1,N1X
         IF((dabs(x(i)-x0fr).le.(0.5d0*xfr)))THEN
            CALL const_param(param,nparam,fo(i),Vo,
     &           cca(i),ccb(i),L(i))
            dfdlV(i) = fric_dfdV(V(i),psi(i),
     &           param,nparam,Mnum(i))*dabs(V(i))
            dfssdlV(i) = fric_dfssdV(V(i),Seff(i),
     &           param,nparam,Mnum(i))*dabs(V(i))

         ENDIF
      ENDDO
c
      RETURN
      END
c
c-----------------------------------------------------------------------
      Subroutine Set_initials(Seff,tauo,V,del,
     &     psi,Mnum,cca,ccb,L,fo,Vo,eta,
     &     nparam,x,N1X,xfr,xfrvw,xLam,x0fr,x0frvw)
      ! This subroutine is called prior to the main time loop, and 
      ! sets the initial condition (effective normal stress, loading, 
      ! slip rate, state variable, and frictional 
      ! constitutive parameters).
      IMPLICIT NONE
c
      ! IO parameters
      !    REAL*8
      ! Seff  : Effective normal stress
      ! tauo  : Loading
      ! V     : Slip rate
      ! del   : Slip
      ! psi   : State variable
      ! cca   : a-value, (a constant)
      ! ccb   : b-value, (a constant)
      ! L     : State evolution distance (characteristic slip distance)
      ! fo    : Reference friction coefficient [scalar]
      ! Vo    : Reference slip rate
      ! x     : x-coordinate
      ! xfr   : Fault dimension in x-direction
      ! xfrvw : VW area dimension in x-direction
      !    INTEGER
      ! Mnum : Material number
      ! N1X   : Number of grids in the x-direction on each node
      INTEGER N1X
      INTEGER Mnum(N1X),nparam,inFlag
      REAL*8 Seff(N1X),tauo(N1X),V(N1X),
     &     del(N1X),x(N1X),psi(N1X),
     &     cca(N1X),ccb(N1X),L(N1X),fo(N1X),Vo,eta

      REAL*8 xfr,xfrvw,xLam,x0fr,x0frvw
c
      ! Other parameters
      ! i     : Index in the x-direction
      ! Vini  : Initial slip rate
      ! Seff0 : Effective normal stress used in the initialization
      INTEGER i
      REAL*8 Vini,Vpl,Seff0,tau0,param(nparam),pi
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Note on the fault constitutive laws
c
c   Material number (Mnum) = m*100+n
c  
c  m : State-evolution
c  1 : Regularized aging law
c  2 : Non-regularized aging law
c  3 : Regularized slip law
c  4 : Non-regularized slip law
c
c  n : Steady-state
c  1 : Usual logarithmic
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ! Common properties. Some of them are not used in certain Mnums
      Seff0 = 50.d6
      Vo = 1.d-6
      pi=3.141592653589793d0
c
      Vpl = 1.d-9
      Vini = Vpl

c     Field properties
      DO i = 1,N1X
            del(i) = 0.d0					! Slip (m)
            Seff(i)=Seff0				        ! Effective Normal Stress (Pa)
            V(i) = Vini					        ! Slip rate (m/s)            

c     Rate-and-State properties
            fo(i) = 0.6d0           				! Reference friction coefficient
            L(i) = 4.d-3                                        ! Characteristic slip weakening distance (m)
            cca(i) = 0.025d0				        ! a-value for R+S friction
            ccb(i) = 0.015d0				        ! b-value for R+S friction
            Mnum(i)=101   					! Material number (see above table)

c     Set VW zone in the middle
            IF(abs(x(i)-x0frvw).LT.(0.5d0*xfrvw))THEN
               cca(i) = 0.010d0
            ENDIF

c     Initial shear stress
            tauo(i)=cca(i)*Seff(i)*dasinh( (V(i)/(2*Vo))*dexp((fo(i)+
     &               ccb(i)*dlog(Vo/V(i)))/cca(i)))

            IF(((x(i)-x0frvw).GT.(-0.5d0*xfrvw)).AND.
     &           (x(i)-x0frvw).LT.(-0.5d0*xfrvw+1.5d3))THEN
               tauo(i)=cca(i)*Seff(i)*dasinh((1.d0/(2*Vo))*dexp((fo(i)+
     &               ccb(i)*dlog(Vo/Vpl))/cca(i)))
            ENDIF

c     Initial state variable
            psi(i) = (cca(i)/ccb(i))*dlog(2*Vo/V(i) * 
     &            dsinh((tauo(i) - eta*V(i))/(cca(i)*Seff(i))))
     &            - fo(i)/ccb(i)

      ENDDO
      RETURN
      END
c
c-------------------------------------------------------------
      Subroutine Mode_convertor(Input,Divider,Integ,Rema)		
      ! Used for dividing up the output points ("H###") amongst the nodes
      ! and determining the points local index
      IMPLICIT NONE
      INTEGER Input,Divider,Integ,Rema
      IF (Mod(Input,Divider).eq.0) THEN
         Integ = INT(Input/Divider) - 1
         Rema = Divider
      ELSE
         Integ = INT(Input/Divider)
         Rema = Mod(Input,Divider)
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
      Subroutine Set_dLmax(dLmax,dfdlV,dfssdlV,hcell,L,Seff,
     &        N1X,xmu,x,xfr,nu,mode,Sfac,x0fr)
c     Determine the maximum slip distance, used to determine maximum timestep
c     Modifiable: dLmax
      IMPLICIT NONE
      INTEGER i,N1X,myrank,mode
      REAL*8 dLmax(N1X),dfdlv(N1X),dfssdlV(N1X),
     &     L(N1X),Seff(N1X),
     &     x(N1X),xfr,nu,Sfac,
     &     hcell,gamma,ro,pi,expr1,expr2,xmu,Xith,
     &     hcell_,xmu_,x0fr
      PARAMETER (pi=3.141592653589793d0)
c
      IF(mode.eq.2)THEN
         xmu_= xmu/(1.d0-nu)
         hcell_= hcell
      ELSEIF(mode.eq.3)THEN
         xmu_= xmu
         hcell_= hcell
      ELSEIF(mode.eq.4)THEN
         xmu_= xmu
         hcell_= hcell*(1.d0-nu)
      ENDIF
      gamma = pi/4.d0
      DO i=1,N1X
         IF(dabs(x(i)-x0fr).le.(0.5d0*xfr))THEN
            expr1=-dfssdlV(i)/dfdlv(i)
            expr2=gamma*xmu_/hcell_*L(i)/(dfdlv(i)*Seff(i))
            ro=expr2-expr1
            IF ((0.25d0*ro*ro-expr2).ge.0.d0) THEN
               Xith=1.d0/ro
            ELSE
               Xith=1.d0-expr1/expr2
            ENDIF
c
            Xith = Xith / Sfac
c
            IF (Xith.GT.(0.33d0)) THEN
               dLmax(i)=0.33d0*L(i)
            ELSE
               dLmax(i)=Xith*L(i)  
            ENDIF
         ENDIF
      ENDDO
      RETURN
      END
c
c-----------------------------------------------------------------
      REAL*8 function bej1sx(x)						
      ! bessel function
      implicit REAL*8(a-h,o-z)
c     
c     computes J1(x)/x for any x (from Numerical Recipes)
c     
      DATA R1,R2,R3,R4,R5,R6/72362614232.D0,-7895059235.D0,
     *     242396853.1D0,-2972611.439D0,15704.48260D0,
     *     -30.16036606D0/
      DATA S1,S2,S3,S4,S5,S6/144725228442.D0,2300535178.D0, 
     *     18583304.74D0,99447.43394D0,376.9991397D0,1.D0/ 
      DATA P1,P2,P3,P4,P5/1.D0,.183105D-2,-.3516396496D-4,
     *     .2457520174D-5,-.240337019D-6/
      DATA Q1,Q2,Q3,Q4,Q5/.04687499995D0,-.2002690873D-3,
     *     .8449199096D-5,-.88228987D-6,.105787412D-6/
      IF(dabs(X).LT.8.d0)THEN 
         Y=X**2 
         BEJ1SX=(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))) 
     *        /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))) 
      ELSE 
         AX=dabs(X) 
         Z=8.d0/AX 
         Y=Z**2 
         XX=AX-2.356194491d0
         BEJ1SX=dsqrt(.636619772d0/AX)*(dcos(XX)*(P1+Y*
     *        (P2+Y*(P3+Y*(P4+Y*P5))))-Z*dsin(XX)*
     *        (Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))/X 
      ENDIF 
      RETURN 
      END 
c
c------------------------------------------------------------------
      REAL*8 FUNCTION W(X)
      ! WOFX
      implicit REAL*8 (a-h,o-z)
c     computes directly 1-int(from 0 to x of J1(x)/x dx)
c     by John Morrissey
C
      DATA A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14
     & /4.D0,-1.066666666666666D1,1.706666666666666D1,
     & -1.625396825396825D1,1.011358024691358D1,-4.413198653198653D0,
     & 1.422569529236196D0,-3.522553120013437D-1,6.906966901987132D-2,
     & -1.098652045228363D-2,1.445845115798365D-3,-1.600144792187913D-4,
     & 1.509880214167057D-5,-1.229043723375708D-6,8.718339712615857D-8/
	DATA B0,B1,B2,B3,B4,B5,B6,B7
     &  /1.25D-1,-5.26410175D-3,1.273027121D-3,-6.73431837D-4,
     &  5.005321894D-4,-3.452926821D-4,1.591757834D-4,-3.364922862D-5/
	DATA C1,C2,C3,C4,C5,C6,C7
     &  /1.757799D-2,-2.274954641D-3,8.472097978D-4,-5.141838574D-4,
     &  3.188374766D-4,-1.387793491D-4,2.835472838D-5/
C
	AX=dabs(X)
	IF (AX.LT.8.d0)THEN
         Y=AX/8.d0
         Z=Y*Y
         W=Y*(A0+Z*(A1+Z*(A2+Z*(A3+Z*(A4+Z*(A5+Z*(A6+Z*(A7+Z*(A8+Z*(A9+
     &   Z*(A10+Z*(A11+Z*(A12+Z*(A13+Z*A14))))))))))))))
	 W=1.d0-W
	ELSE
	 Y=8.d0/AX
	 Z=Y*Y
	 XX=X-.785398163d0
	 W=dsqrt(.636619772/AX)*(COS(XX)*Y*(B0+Z*(B1+Z*(B2+Z*(B3+Z*(B4+
     &   Z*(B5+Z*(B6+Z*B7)))))))+SIN(XX)*Z*(C1+Z*(C2+Z*(C3+Z*(C4+Z*(C5+
     &   Z*(C6+Z*C7)))))))
	ENDIF
	W=W*SIGN(1.,X)
	RETURN
	END
c
c--------------------------------------------------------------------
	REAL*8 function bessj1(x)
        ! Bessel function
	implicit REAL*8 (a-h,o-z)
c	computes J1(x) for any x (from Numerical Recipes)
c
      	SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     #       s1,s2,s3,s4,s5,s6
	DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,
     *	   242396853.1d0,-2972611.439d0,15704.48260d0,
     *     -30.16036606d0/
	data s1,s2,s3,s4,s5,s6/144725228442.d0,2300535178.d0,
     *     18583304.74d0,99447.43394d0,376.9991397d0,1.d0/
      	DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,
     *	   .2457520174d-5,-.240337019d-6/
	DATA q1,q2,q3,q4,q5/.04687499995d0,-.2002690873d-3,
     *     .8449199096d-5,-.88228987d-6,.105787412d-6/
      	IF(dabs(x).LT.8.d0)THEN
          y=x**2
          bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/
     *        (s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      	ELSE
          ax=dabs(x)
          z=8.d0/ax
          y=z**2
          xx=ax-2.356194491d0
          bessj1=dsqrt(.636619772/ax)*(dcos(xx)*
     *        (p1+y*(p2+y*(p3+y*(p4+y*p5))))-z*sin(xx)*
     *        (q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1.,x)
      	ENDIF
      	RETURN
      	END
c
c------------------------------------------------------------------
	REAL*8 function bessj0(X)
        ! Bessel function
	implicit REAL*8 (a-h,o-z)
c	computes J0(x) for any x (from Numerical Recipes)
c
	DATA P1,P2,P3,P4,P5/1.D0,-.1098628627D-2,.2734510407D-4,
     *	  -.2073370639D-5,.2093887211D-6/
	DATA Q1,Q2,Q3,Q4,Q5/-.1562499995D-1,.1430488765D-3,
     *    -.6911147651D-5,.7621095161D-6,-.934945152D-7/
	DATA R1,R2,R3,R4,R5,R6/57568490574.D0,-13362590354.D0,
     *     651619640.7D0,
     *	  -11214424.18D0,77392.33017D0,-184.9052456D0/
	DATA S1,S2,S3,S4,S5,S6/57568490411.D0,1029532985.D0,
     *	  9494680.718D0,59272.64853D0,267.8532712D0,1.D0/
	IF(ABS(X).LT.8.)THEN
	  Y=X**2
	  BESSJ0=(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))
     *		/(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
	ELSE
	  AX=dabs(X)
	  Z=8.d0/AX
	  Y=Z**2
	  XX=AX-.785398164d0
	  BESSJ0=dsqrt(.636619772/AX)*(COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y
     *		*P5))))-Z*SIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
	ENDIF
	RETURN
	END
c
c----------------------------------------------------------------------
	REAL*8 function xker2(x,rho)
        ! computes the kernel of the convlution for the mode II crack problem
        ! Note rho=(cp/cs)**2

	implicit REAL*8 (a-h,o-z)

	xa=dsqrt(rho)*x
	xker2=bej1sx(x)+4.d0*x*(w(xa)-w(x))
     #       -4.d0/dsqrt(rho)*bessj0(xa)+3.d0*bessj0(x)
	RETURN
	end
c
c-----------------------------------------------------------------------
        REAL*8 function integral(xker2,n,s1,s2,rho)
        implicit REAL*8 (a-h,o-z)

        ds = (s2-s1)/n
        integral = 0.5*xker2(s1,rho)
        DO i = 1,n-1
          integral = integral + xker2(s1+i*ds,rho)
        ENDDO
        integral = ds*(integral+0.5*xker2(s2,rho))
        RETURN
        end

c-----------------------------------------------------------------------
c   Fault constitutive properties
c-----------------------------------------------------------------------
c
c   Material number (Mnum) = m*100+n
c  
c  m : State-evolution
c  1 : Regularized aging law
c  2 : Non-regularized aging law
c  3 : Regularized slip law
c  4 : Non-regularized slip law
c
c  n : Steady-state
c  1 : Usual logarithmic
c
c-----------------------------------------------------------------------
      Subroutine const_param(param,nparam,fo,Vo,
     &                 cca,ccb,L)
      ! This subroutine construct an array, param to simplify the 
      ! structure of other subroutines and functions.
      IMPLICIT NONE
c
      ! IO parameters
      INTEGER nparam
      REAL*8 param(nparam),fo,Vo,cca,ccb,L
c
      param(1)=fo
      param(2)=Vo
      param(3)=cca
      param(4)=ccb
      param(5)=L
c
      RETURN
      END
c
c-----------------------------------------------------------------------
      Subroutine s_evolution(psin,psipr,psis,dt,Vpr,Vs,Seff,
     $                 param,nparam,Mnum)
      ! This subroutine evolve the state variable, psi.
      IMPLICIT NONE
c
      ! IO parameters
      ! INTEGER
      ! nparam : Number of the frictional parameters
      ! Mnum   : Material number
      !
      ! REAL*8
      ! psin   : New value of the state variable, psi
      ! psipr  : Previous value of the state variable, psi
      ! dt     : Length of timestep
      ! Vpr    : Previous slip rate
      ! Vs     : Predicted slip rate
      ! param  : Array containing a set of frictional parameters
      !  param(1) : Ref. friction
      !  param(2) : Ref. slip rate
      !  param(3) : a, not necessarily df/dln(V)|_{state}
      !  param(4) : b, not necessarily df/dln(V)|_{state} - dfss/dln(V)
      !                This is used to define a state variable, theta.
      !  param(5) : L

      INTEGER nparam,Mnum
      REAL*8 psin,psipr,psis,dt,Vpr,Vs,param(nparam) !            Definition of psis is added. HIRO 20100903
c
      ! Other parameters
      ! INTEGER
      ! Snum    : A part of Mnum for the formulation of state-evolution
      !
      ! REAL*8
      ! psiss   : Steady state psi
      ! thetapr : Previous state in theta(=exp(psi/b)), nondimensional
      ! thetass : Steady state theta
      ! thetan  : New value of the state variable, theta
      INTEGER Snum
      REAL*8 psiss,thetapr,thetass,thetan,thetas,zhp,zhpo,Vpr_,Vs_, ! Definition of thetas is added. HIRO 20100903
     &     Seff
      ! External functions
      REAL*8 fric_fss,fric_f
c
      Snum = Mnum-mod(Mnum,100)
      Vpr_ = dabs(Vpr)
      Vs_ = dabs(Vs)
c
      IF(param(4).eq.0.d0)THEN
      ! If b=0, the const. law is purely direct. 
         psin = psipr
      ELSEIF(Snum.eq.100)THEN
      ! Regularized aging law. 
      ! f = a asinh (v/2vo exp((fo+bpsi)/a))
      ! Not suitable for a variable a-value.
         ! Using an analytic expression, second order
         thetapr = dexp(psipr)
         thetass = param(2)/Vpr_
c
         IF((1.d0-dexp(-0.5d0*Vpr_*dt/param(5))).le.1.d-6)THEN
            thetan = thetapr + 
     &           0.5d0*Vpr_*dt/param(5)*(thetass-thetapr)
         ELSE
            thetan=thetass+(thetapr-thetass)
     &           *dexp(-0.5d0*Vpr_*dt/param(5))
         ENDIF
c
         thetas = dexp(psis) !                                           Added. HIRO 20100903
         thetass = param(2)/Vs_
c     
         IF((1.d0-dexp(-0.5d0*Vs_*dt/param(5))).le.1.d-6)THEN
            thetan = thetan + 
     &           0.5d0*Vs_*dt/param(5)*(thetass-thetas) !                Modified. thetan => thetas . This is needed to make the integration second order accurate interseismically. See Noda and Lapusta, 2010. HIRO 20100903 
         ELSE
            thetan=thetass+(thetan-thetass)
     &           *dexp(-0.5d0*Vs_*dt/param(5))
         ENDIF
         psin = dlog(thetan)
c     
      ELSEIF(Snum.eq.200)THEN
      ! Non-regularized aging law.
      ! f = Direct + b*psi
      ! Constant b-value (log-t healing rate) is needed.
         ! Using an analytic expression, 2nd order
         thetapr = dexp(psipr)
         psiss = (fric_fss(Vpr_,Seff,param,nparam,Mnum)
     &          - fric_f(Vpr_,0.d0,param,nparam,Mnum))/param(4)
         thetass = dexp(psiss)
         IF((1.d0-dexp(-0.5d0*Vpr_*dt/param(5))).le.1.d-6)THEN
            thetan = thetapr + 
     &           0.5d0*Vpr_*dt/param(5)*(thetass-thetapr)
         ELSE
            thetan=thetass+(thetapr-thetass)
     &           *dexp(-0.5d0*Vpr_*dt/param(5))
         ENDIF
         psiss = (fric_fss(Vs_,Seff,param,nparam,Mnum)
     &          - fric_f(Vs_,0.d0,param,nparam,Mnum))/param(4)
         thetass = dexp(psiss)
         IF((1.d0-dexp(-0.5d0*Vs_*dt/param(5))).le.1.d-6)THEN
            thetan = thetan + 
     &           0.5d0*Vs_*dt/param(5)*(thetass-thetan)
         ELSE
            thetan=thetass+(thetan-thetass)
     &           *dexp(-0.5d0*Vs_*dt/param(5))
         ENDIF
         psin = dlog(thetan)
c
      ELSEIF(Snum.eq.300)THEN
      ! Regularized slip law. 
      ! f = a asinh (V/2Vo exp((fo+bpsi)/a))
      ! Not suitable for a variable so-called a-value.
         ! Using analytic expression, 2nd order
         psiss = (param(3)*dlog(2.d0*param(2)/Vpr_*dsinh(
     &        fric_fss(Vpr_,Seff,param,nparam,Mnum)/param(3)))-param(1))
     &        /param(4)
         IF((1.d0-dexp(-Vpr_*0.5d0*dt/param(5))).le.1.d-6)THEN
            psin = psipr + Vpr_*0.5d0*dt/param(5)*(psipr-psiss)
         ELSE
            psin = psiss + (psipr-psiss)*dexp(-Vpr_*0.5d0*dt/param(5))
         ENDIF
         psiss = (param(3)*dlog(2.d0*param(2)/Vpr_*dsinh(
     &        fric_fss(Vs_,Seff,param,nparam,Mnum)/param(3)))-param(1))
     &        /param(4)
         IF((1.d0-dexp(-Vs_*0.5d0*dt/param(5))).le.1.d-6)THEN
            psin = psin + Vs_*0.5d0*dt/param(5)*(psin-psiss)
         ELSE
            psin = psiss + (psin-psiss)*dexp(-Vs_*0.5d0*dt/param(5))
         ENDIF
c
      ELSEIF(Snum.eq.400)THEN
      ! Non-regularized slip law.
      ! f = Direct + psi
      ! Only L (param(5)) is needed.
         ! Using analytic expression, 2nd order
         psiss = fric_fss(Vpr_,Seff,param,nparam,Mnum)
     &          - fric_f(Vpr_,0.d0,param,nparam,Mnum)
         IF((1.d0-dexp(-Vpr_*0.5d0*dt/param(5))).le.1.d-6)THEN
            psin = psipr + Vpr_*0.5d0*dt/param(5)*(psipr-psiss)
         ELSE
            psin = psiss + (psipr-psiss)*dexp(-Vpr_*0.5d0*dt/param(5))
         ENDIF
         psiss = fric_fss(Vs_,Seff,param,nparam,Mnum)
     &          - fric_f(Vs_,0.d0,param,nparam,Mnum)
         IF((1.d0-dexp(-Vs_*0.5d0*dt/param(5))).le.1.d-6)THEN
            psin = psin + Vs_*0.5d0*dt/param(5)*(psin-psiss)
         ELSE
            psin = psiss + (psin-psiss)*dexp(-Vs_*0.5d0*dt/param(5))
         ENDIF
c

      ENDIF
      RETURN
      END
c
c-----------------------------------------------------------------------
      REAL*8 function fric_dfdV(V,psi,param,nparam,Mnum)
      ! This function gives df/dV with fixed state
      IMPLICIT NONE
c
      ! IO parameters
      ! INTEGER
      ! nparam : Number of the frictional parameters
      ! Mnum   : Material number
      !
      ! REAL*8
      ! V      : Slip rate
      ! psi    : State variable, psi
      ! param  : An array containing the frictional parameters
      INTEGER nparam,Mnum
      REAL*8 V,psi,param(nparam)
c
      ! Other paramters
      ! INTEGER
      ! Snum : Number for the formulation of the state evoltuion
      INTEGER Snum
      REAL*8 ps
c
      Snum = Mnum-mod(Mnum,100)

      IF((Snum.eq.100).or.(Snum.eq.300))THEN 
c     Regularized rate- and state- law
         ps = param(1)+param(4)*psi
         fric_dfdV = param(3)
     &        / dsqrt(1.d0+(V/(2.d0*param(2))*dexp(ps/param(3)))**2)
     &        / (2.d0*param(2))*dexp(ps/param(3))
      ELSEIF((Snum.eq.200).or.(Snum.eq.400))THEN
      ! Non-regularized rate- and state- law
         fric_dfdV = param(3)/V
      ENDIF
c
      RETURN
      END
c
c-----------------------------------------------------------------------
      REAL*8 function fric_f(V,psi,param,nparam,Mnum)
      ! This function gives the value of friction coefficient.
      IMPLICIT NONE
c
      ! IO parameters
      ! INTEGER
      ! nparam : Number of the frictional parameters
      ! Mnum   : Material number
      ! REAL*8
      ! V     : Slip rate
      ! psi   : State variable, psi
      ! param : An array containing the frictional parameters
      INTEGER nparam,Mnum
      REAL*8 V,psi,param(nparam)
c
      ! Other parameters
      ! Snum : Number for the formulation of the state evoltuion

      INTEGER Snum
      REAL*8 lz,lzy,R,Q,ps,help
c
      Snum = Mnum-mod(Mnum,100)

      IF((Snum.eq.100).or.(Snum.eq.300))THEN 
c     Regularized rate-and-state
         ps = param(1)+param(4)*psi
         fric_f = param(3)*dasinh(V/(2.d0*param(2))
     &        *dexp(ps/param(3)))
      ELSEIF((Snum.eq.200).or.(Snum.eq.400))THEN
c     Non-regularized rate-and-state
         fric_f = param(1) + param(3)*dlog(V/param(2))+param(4)*psi  
      ENDIF
c
      RETURN
      END
c
c-----------------------------------------------------------------------
      REAL*8 function fric_dfssdV(V,Seff,param,nparam,Mnum)
      ! This function gives dfss/dV
      IMPLICIT NONE
c
      ! IO parameters
      ! INTEGER
      ! nparam : Number of the frictional constitutive parameters
      ! Mum    : Material number
      !
      ! REAL*8
      ! V     : Slip rate
      ! param : An array containing the frictional constitutve parameters
      INTEGER nparam,Mnum
      REAL*8 V,param(nparam),Seff
c
      ! Other paramters
      ! help : helper function #1
      REAL*8 help
c
      ! External functions
      REAL*8 fric_dfdV
c
      help = 0.5d0 * dexp(param(1)/param(3)) 
     &     * (V/param(2))**(1-param(4)/param(3)) 
      fric_dfssdV = param(3)/sqrt(1.d0+help*help)*help 
     &     *(1.d0-param(4)/param(3))/V 
      

      RETURN
      END
c
c-----------------------------------------------------------------------
      REAL*8 function fric_fss(V,Seff,param,nparam,Mnum) 
      ! This function gives the steady state friction coefficient.
      IMPLICIT NONE
c
      ! IO parameters
      ! INTEGER
      ! nparam   : Number of the frictional parameters
      ! Mnum     : Material number
      ! REAL*8
      ! V        : Slip rate
      ! Seff     : Effective confining stress
      ! param    : An array containing the frictional parameters
      INTEGER nparam,Mnum
      REAL*8 V,param(nparam),Seff
c
      ! External functions
      REAL*8 fric_f
c
      fric_fss = fric_f(V,dlog(param(2)/V),param,nparam,Mnum)
c
      RETURN
      END
c
c-----------------------------------------------------------------------
c     End of the fault constitutive properties
c-----------------------------------------------------------------------
