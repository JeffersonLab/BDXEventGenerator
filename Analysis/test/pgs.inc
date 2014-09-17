c ------------------------------------------------------------------------------------------
c PGS.INC - include file for PGS package 
c
c           Version 3.3
c ------------------------------------------------------------------------------------------
c
c   This include file is to be used with the PGS package.  See
c   http://www.physics.rutgers.edu/~jconway/soft/pgs/pgs.html
c   for documentation and details.
c
c ------------------------------------------------------------------------------------------


c STDHEP common blocks (see L. Garren, STDHEP 4.10, PM0091)  -------------------------------

c generated particle list

      integer    nmxhep
      parameter (nmxhep=4000)

      integer nevhep,nhep,isthep,idhep,jmohep,jdahep

      double precision phep,vhep

      common /hepevt/ nevhep,           ! event number
     .                nhep,             ! number of entries in record
     .                isthep(nmxhep),   ! status code
     .                idhep(nmxhep),    ! particle ID (PDG standard)
     .                jmohep(2,nmxhep), ! index to first and second particle mothers
     .                jdahep(2,nmxhep), ! index to first and last daughter particles
     .                phep(5,nmxhep),   ! 4-vector and mass
     .                vhep(4,nmxhep)    ! (x,y,z) of production, and production time (mm/c)

      integer numhep                    ! for symmetry
      equivalence(numhep,nhep)          ! with other lists

      integer numgen                    ! for symmetry
      equivalence(numgen,nhep)          ! with other lists

c table of particle properties

      integer    nmxln

      parameter (nmxln=2000)

      integer idt
      real*4 stmerr,stwerr

      real*8 stmass,stwidth

      character*21 stname

      common /stdtbl/ idt(nmxln),       ! particle ID (PDG standard)
     .                stmerr(2,nmxln),  ! positive (1) and negative (2) error on mass
     .                stwerr(2,nmxln),  ! positive (1) and negative (2) error on width
     .                stmass(nmxln),    ! particle mass (GeV/c^2 for this, and rest)
     .                stwidth(nmxln),   ! particle width
     .                stname(nmxln)     ! particle name


c STDHEP logical unit numbers and I/O information

      integer lnhwrt,lnhrd,lnhout,lnhdcy,lnhrdm

      common /heplun/ lnhwrt,           ! logical unit number for writing
     .                lnhrd,            ! logical unit number for reading
     .                lnhout,           ! logical unit number for output
     .                lnhdcy,           ! logical unit number for decay table (NYI)
     .                lnhrdm(15)        ! logical unit number for reading multiple streams

      real stdecom,stdxsec
      double precision stdseed1,stdseed2
      integer nevtreq,nevtgen,nevtwrt

      common /stdcm1/ stdecom,  ! center of mass energy
     .                stdxsec,  ! cross section
     .                stdseed1, ! random number seed 1
     .                stdseed2, ! random number seed 2
     .                nevtreq,  ! events to generate
     .                nevtgen,  ! actually generated
     .                nevtwrt   ! written to output file
      save /stdcm1/


c ---------------------------------------------------------------------------------------------

c PGS functions for ET, pt, etc. of generated, reconstructed, clusters and tracks

      real*8 et_gen,pt_gen,p_gen,eta_gen,phi_gen
      real*8 et_obj,pt_obj,p_obj,eta_obj,phi_obj
      real*8 et_clu,pt_clu,p_clu,eta_clu,phi_clu
      real*8 et_trk,pt_trk,p_trk,eta_trk,phi_trk

      real*8 cos_theta,cos_del_phi,del_phi,pgs_cluster_width

      real*8 v4mass,v4mass2,v4et,v4p,v4pt,v4eta,v4phi

      logical pgs_cal_fid_cut

c PGS detector parameters

      integer netacal,nphical

      real*8 deta,dphi,crack_frac,rmax,bfield,sagitta_res
      real*8 ecal_res,hcal_res,met_res
      real*8 seed_thresh,shoulder_thresh,conesize
      real*8 eff_track,eta_max_track,min_track_pt


      common /pgspar/ netacal,		! number of eta cells in calorimeter
     .                nphical,		! number of phi cells in calorimeter
     .                deta,		! eta width of calorimeter cells
     .                dphi,		! phi width of calorimeter cells  (2*pi/nphical)
     .                ecal_res,         ! electromagnetic calorimeter resolution * sqrt(E)
     .                hcal_res,         ! hadronic calorimeter resolution * sqrt(E)
     .                met_res,          ! MET resolution
     .                crack_frac,       ! calorimeter cell edge "crack" fraction
     .                seed_thresh,      ! calorimeter cluster seed tower threshold
     .                shoulder_thresh,  ! calorimeter shoulder tower threshold
     .                conesize,         ! cluster finder cone size
     .                rmax,             ! outer radius of tracking (m)
     .                bfield,           ! magnetic field (T)
     .                sagitta_res,      ! sagitta resolution (m)
     .                eff_track,	! track finding efficiency
     .                min_track_pt,     ! minimum track pt (GeV/c)
     .                eta_max_track     ! maximum tracking coverage in eta

c PGS event header and control information

      integer nevpgs,nprpgs

      character*6 optpgs
      character*80 pgs_input_file,pgs_output_file,pgs_param_file
      character*80 pgs_isajet_decay,pgs_isajet_cards,pgs_pythia_cards

      common /pgsevt/ nevpgs,            ! number of events to generate/read
     .                nprpgs,            ! number of events to print out 
     .                optpgs,            ! type of run: 'PYTHIA', 'ISAJET', 'FILE', ...
     .                pgs_input_file,    ! input file
     .                pgs_output_file,   ! output file
     .                pgs_param_file,    ! detector parameter file
     .                pgs_isajet_decay,  ! ISAJET decay table file name
     .                pgs_isajet_cards,  ! ISAJET card file name
     .                pgs_pythia_cards   ! PYTHIA card file name

c PGS track list

      integer ntrkmx
      parameter (ntrkmx=500)

      integer numtrk,dumtrk,indtrk

      real*8 ptrk,qtrk

      common /pgstrk/ numtrk,dumtrk,          ! number of tracks
     .                indtrk(ntrkmx),         ! index to HEPEVT particle
     .                ptrk(3,ntrkmx),         ! track 3-vector 
     .                qtrk(ntrkmx)            ! track charge


c PGS calorimeter tower arrays

      real*8 pi
      parameter(pi=3.141592654)      

      integer nphimax,netamax   
      parameter (nphimax=600)   
      parameter (netamax=600)   

      real*8 ecal,hcal,met_cal,phi_met_cal

      common /pgscal/ ecal(netamax,nphimax),  ! electromagnetic energy in each tower
     .                hcal(netamax,nphimax),  ! hadronic energy in each tower
     .                met_cal,                ! calorimeter missing ET
     .                phi_met_cal             ! calorimeter missing ET phi


c PGS calorimeter cluster list

      integer nclumx

      parameter (nclumx=50)

      integer cclu,numclu,dumclu,etaclu,phiclu,mulclu

      real*8 pclu,etclu,emclu,ehclu,efclu,widclu

      common /pgsclu/ cclu(netamax,nphimax),  ! map of cluster indices
     .                numclu,dumclu,          ! number of clusters in list
     .                pclu(5,nclumx),         ! cluster 4 vector and mass
     .                etaclu(nclumx),         ! cluster seed tower eta
     .                phiclu(nclumx),         ! cluster seed tower phi
     .                emclu(nclumx),          ! cluster electromagnetic energy
     .                ehclu(nclumx),          ! cluster hadronic energy
     .                efclu(nclumx),          ! cluster electromagnetic fraction
     .                widclu(nclumx),         ! cluster width sqrt(deta**2+dphi**2)
     .                mulclu(nclumx)          ! cluster tower multiplicity


c PGS trigger object list

      integer    ntrgmx
      parameter (ntrgmx=500)

      integer numtrg,dumtrg,indtrg,typtrg

      real*8 vectrg

      common /pgstrg/ numtrg,dumtrg,          ! number of trigger objects
     .                indtrg(ntrgmx),         ! index to HEPEVT particle (where relevant)
     .                typtrg(ntrgmx),         ! reconstructed type:  0=photon
                                              !                      1=electron
                                              !                      2=muon
                                              !                      3=tau (hadronic)
                                              !                      4=jet
                                              !                      5=detached vertex
                                              !                      6=MET
     .                vectrg(10,ntrgmx)       ! trigger object vector: 1 = eta
                                              !                        2 = phi
                                              !                        3 = ET of cluster
                                              !                        4 = cluster #
                                              !                        5 = pt of track (if any)
                                              !                        6 = track #

c PGS reconstructed object list

      integer    nobjmx
      parameter (nobjmx=500)

      integer numobj,dumobj,indobj,typobj
      real*8 pobj,qobj,vecobj

      common /pgsrec/ numobj,dumobj,          ! number of reconstructed objects
     .                indobj(nobjmx),         ! index to HEPEVT particle (where relevant)
     .                typobj(nobjmx),         ! reconstructed type:  0 = photon
                                              !                      1 = electron
                                              !                      2 = muon
                                              !                      3 = tau (hadronic)
                                              !                      4 = jet
                                              !                      5 = heavy charged
     .                pobj(4,nobjmx),         ! four vector of reconstructed object
     .                qobj(nobjmx),           ! charge of reconstructed object
     .                vecobj(10,nobjmx)       ! interesting object quantities

c            --------------------------
c            table of vecobj quantities
c            --------------------------
c
c    -------------------------------------------------------------------------------------
c     type            1           2          3        4        5        6        7
c    -------------------------------------------------------------------------------------
c    0  photon     EM energy  HAD energy  track E   N(trk)   width      -        -
c    1  electron    "   "      "     "       "        "        -        -        -
c    2  muon        "   "      "     "       "        "     trk iso E   -        -
c    3  tau         "   "      "     "       "        "      width    mtau     ptmax
c    4  jet         "   "      "     "       "        "        "      flavor   c,b tags ->
c    -------------------------------------------------------------------------------------
c
c  b, c tagging: vecobj(7,iobj) non-zero if charm tag (jet prob. alg.)
c                vecobj(8,iobj) non-zero if b     tag ( "    "    "  )
c                vecobj(9,iobj) non-zero if b     tag (impact method)
c
c    --> all algorithms include rates for gluon, uds, c and b jets
c
