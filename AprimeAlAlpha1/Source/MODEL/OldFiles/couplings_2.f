c----------------------------------------------------------------------
C  couplings.f 
c----------------------------------------------------------------------
c  This files takes the inputs of the standard model from a Les Houches 
c  file (param_card.dat) and calculates all the couplings that end up
c  in the Feynman rules, i.e., in the HELAS routines.
c   
c  With the convention of the New Web Generation MadEvent in this
c  part no non-trivial computation is made. The SM model secondary
c  parameters have been all calculated by the SM calculator, SMCalc
c  which produced param_card.dat.
c
c  The only exception to the above rule is for alpha_S. In the case
c  where a pdf is used in the initial state, then the values as(MZ)
c  set in the les houches card is superseeded.
c  Scales and pdf's are set in the run_card.dat.
c
c This file contains the following routines:
c 
c- madgraph original call to set the parameters
c- lh_readin is called from here.
c  It talks to the lh_readin through the common block values.
c      subroutine setpara
c
c-routine to read the LH parameters in
c      subroutine lh_readin
c
c-to set
c      subroutine set_it(npara,ivalue,value,name,id,
c      subroutine case_trap(string,length)
c      subroutine no_spaces(buff,nchars)
c---------------------------------------------------------------------- 


      subroutine setpara(param_name,readlha)
c***********************************************************************
c This subroutine sets up the HELAS couplings of the STANDARD MODEL.
c***********************************************************************
      implicit none
c
c local
c
      character*(*) param_name
      logical readlha
      logical coupl_debug
      integer i
      real*8 dum
c
c     common file with the couplings
c
      include 'coupl.inc'
      include 'input.inc'
c
c     local
c
      double precision  v
      double precision  ee, ee2, ez, ey, sw, cw, sc2, sin2w, wm
      double precision  gwne, gwud, lambda, lam4, xt, rew, rqcd
      double precision  alphas, alfa, alfaw, mfrun
c     following line added by Rouven
      double precision aval, dval, apval, tval, Anuc, fullcoupling 
c
      external          alphas, alfa, alfaw, mfrun
c
c     Common to lh_readin and printout
c
      double precision  alpha, gfermi, alfas
      double precision  mtMS,mbMS,mcMS,mtaMS !MSbar masses
      double precision  Vud,Vus              !CKM matrix elements
      common/values/    alpha,gfermi,alfas,   
     &                  mtMS,mbMS,mcMS,mtaMS,
     &                  Vud
c
c constants
c
      double complex  ci
      parameter( ci = ( 0.0d0, 1.0d0 ) )
      double precision  Zero, One, Two, Three, Four, Half, Rt2
      parameter( Zero = 0.0d0, One = 1.0d0, Two = 2.0d0 )
      parameter( Three = 3.0d0, Four = 4.0d0, Half = 0.5d0 )
      parameter( Rt2   = 1.414213562d0 )
      double precision  Pi, Fourpi
      parameter( Pi = 3.14159265358979323846d0 )
      parameter( Fourpi = Four * Pi )

      integer Nin, Nout
      parameter (Nin=2,Nout=4)  ! for A' this time
!      parameter (Nin=2,Nout=6)  ! for A'

c
c     alfas and run
c************
c Uncomment the following lines in order to use alphas from the PDF
c      include '../alfas.inc'
c      include '../run.inc'
c***********
c-------------------------------------------
c Rouven added the following block to implement form 
c factors on 8/11/09:
c Momenta of particles in event
      include '../genps.inc'
      double precision pp(0:3,max_particles)
      common/momenta_pp/pp 

      coupl_debug=.false.

c------------------------------------------
c Start calculating the couplings for HELAS
c------------------------------------------
c
      if(readlha) then 
         call lh_readin(param_name)
         G = DSQRT(4d0*PI*ALFAS) ! use setting of the param_card.dat @ NLO
      endif
c     

      GG(1) = -G
      GG(2) = -G     
c
c auxiliary local values
c
      wm = sqrt(zmass**2/Two+
     $     sqrt(zmass**4/Four-Pi/Rt2*alpha/gfermi*zmass**2))
      sin2w  = One-(wm/zmass)**2
      cw  = sqrt( One - sin2w )
      ee2 = alpha * Fourpi
      sw  = sqrt( sin2w )
      ee  = sqrt( ee2 )
      ez  = ee/(sw*cw)
      ey  = ee*(sw/cw)
      sc2 = sin2w*( One - sin2w )
      v   = Two*wm*sw/ee   ! the wmass is used to calculate v
      lambda = hmass**2 / (Two * v**2)
c
c vector boson couplings
c
      gw   = ee/sw
      gwwa = ee
      gwwz = ee*cw/sw
c
c gauge & higgs boson coupling constants
c
      gwwh  = dcmplx( ee2/sin2w*Half*v, Zero )
      gzzh  = dcmplx( ee2/sc2*Half*v, Zero )
      ghhh  = dcmplx( -hmass**2/v*Three, Zero )
      gwwhh = dcmplx( ee2/sin2w*Half, Zero )
      gzzhh = dcmplx( ee2/sc2*Half, Zero)
      ghhhh = ghhh/v
c
c fermion-fermion-vector couplings
c
      gal(1) = dcmplx(  ee          , Zero )
      gal(2) = dcmplx(  ee          , Zero )
      gau(1) = dcmplx( -ee*Two/Three, Zero )
      gau(2) = dcmplx( -ee*Two/Three, Zero )
      gad(1) = dcmplx(  ee/Three    , Zero )
      gad(2) = dcmplx(  ee/Three    , Zero )

      gwf(1) = dcmplx( -ee/sqrt(Two*sin2w), Zero )
      gwf(2) = dcmplx(  Zero              , Zero )

      gzn(1) = dcmplx( -ez*Half                     , Zero )
      gzn(2) = dcmplx(  Zero                        , Zero )
      gzl(1) = dcmplx( -ez*(-Half + sin2w)          , Zero )
      gzl(2) = dcmplx( -ey                          , Zero )
      gzu(1) = dcmplx( -ez*( Half - sin2w*Two/Three), Zero )
      gzu(2) = dcmplx(  ey*Two/Three                , Zero )
      gzd(1) = dcmplx( -ez*(-Half + sin2w/Three)    , Zero )
      gzd(2) = dcmplx( -ey/Three                    , Zero )

c fermion-fermion-Higgs couplings (complex) hff(2)
c
c NOTE: the running mass is evaluated @ the same order 
c nloop of alpha_s set by the PDF choice
c 

      if(mtMS.gt.1d0) then
         ghtop(1) = dcmplx( -mtMS/v, Zero )
      else
         ghtop(1) = dcmplx( Zero,Zero)
      endif
      ghtop(2) = ghtop(1)

      if(mbMS.gt.1d0) then
         ghbot(1) = dcmplx( -mbMS/v, Zero )
      else
         ghbot(1) = dcmplx( Zero, Zero )
      endif
      ghbot(2) = ghbot(1)
      
      if(mcMS.gt.1d0) then
         ghcha(1) = dcmplx( -mcMS/v, Zero )
      else
         ghcha(1) = dcmplx( Zero, Zero )
      endif
      ghcha(2) = ghcha(1)

      ghtau(1) = dcmplx( -mtaMS/v, Zero )
      ghtau(2) = ghtau(1)

c
c     CKM matrix: 
c     symmetric 3x3 matrix, Vud=Vcs, Vus=Vcd Vcb=Vub=0
c
c     >>>>>>>>>>>>>>>***** NOTE****<<<<<<<<<<<<<<<<<<<<<<<<<
c     these couplings matter only when interaction_CKM.dat
c     is used to generate all the diagrams with off-diagonal
c     couplings. The default of MadEvent is a diagonal
c     CKM matrix.

	  Vus=DSQRT(1d0-Vud**2)
      do i=1,2
         gwfc(i) = gwf(i)*Vud
         gwfs(i) = gwf(i)*Vus
         gwfm(i) =-gwf(i)*Vus
      enddo

c---------------------------------------------------------
c Set Photon Width to Zero, used by symmetry optimization
c---------------------------------------------------------

      awidth = 0d0
c************************************            
c UserMode couplings
c************************************

      Anuc = 26.98
      aval = 111.0/( elemass*Znuc**(1.0/3.0) )
c      aval = 106.0/( elemass*Znuc**(1.0/3.0) ) ! Beryllium
      dval = 0.164/Anuc**(2.0/3.0)
      apval = 773.0/( elemass*Znuc**(2.0/3.0) ) 
c      apval = 571.4/( elemass*Znuc**(2.0/3.0) ) ! Beryllium

      tval = -(pp(0,Nin)-pp(0,Nout))**2
     $       +(pp(1,Nin)-pp(1,Nout))**2
     $       +(pp(2,Nin)-pp(2,Nout))**2
     $       +(pp(3,Nin)-pp(3,Nout))**2


      fullcoupling = ee * (Znuc**2 *aval**4 *
     $     tval**2/((1+aval**2*tval)*(1+tval/dval))**2 +
     $     Znuc * apval**4 * tval**2 * (1+1.9276*tval)**2 /
     $     ((1+apval**2*tval)*(1+1.40845*tval)**4)**2 )**0.5

      if(coupl_debug) then
         do i=1,7
            write(*,*) "P" , i,"= ", 
     $           " {" , pp(0,i), ", " , pp(1,i),
     $           ", " , pp(2,i), ", " , pp(3,i), "}"
         enddo
         
         write(*,*) "   (t=(P",Nin,"-P",Nout,")^2=", tval, ")",
     $        "...   coupl/eZ=",fullcoupling/(ee*Znuc)
      endif

      GAN(1)=dcmplx(fullcoupling, Zero) 
      GAN(2)=dcmplx(fullcoupling, Zero) 
      
c      write(*,*) " fullcoupling= ",fullcoupling

      GAF(1)=dcmplx( ee      ,Zero)
      GAF(2)=dcmplx( ee      ,Zero)
c      GAN(1)=dcmplx( ee*Znuc ,Zero) ! add formfactors
c      GAN(2)=dcmplx( ee*Znuc ,Zero) ! add formfactors
      GEAP(1)=dcmplx( ee*epsilon ,Zero)
      GEAP(2)=dcmplx( ee*epsilon ,Zero)
      GEAPX(1)=dcmplx( 11.7047*ee*epsilon*Zero ,Zero)
      GEAPX(2)=dcmplx( 11*7047*ee*epsilon*Zero ,Zero)  

c----------------------------
c end subroutine coupsm
c----------------------------


      return
      end

      
