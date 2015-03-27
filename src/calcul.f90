!###############################################################################

!                            subroutine calcul

!###############################################################################


!*******************************************************************************
subroutine calcul
  !*******************************************************************************

  !===============================================================================
  !modules
  !===============================================================================
  !-------------------------------------------------------------------------------
  !Modules de definition des structures et variables => src/mod/module_...f90
  !-------------------------------------------------------------------------------
  use module_canopee
  use module_commun
  use module_connexions
  use module_connexions_lim
  use module_coordonnees
  use module_elevation_surf_libre
  use module_equations
  use module_grille
  use module_impressions
  use module_intercom
  use module_lecture_fichier
  use module_metriques
  use module_mpi
  use module_mpi2
  use module_mpi3
  use module_mpif
  use module_mumps
  use module_objet
  use module_reprise
  use module_statistiques
  use module_struct_psm
  use module_tests_arret
  use module_thermophysique
  use module_unites
  use module_utilitaires
  use module_intercom
  !-------------------------------------------------------------------------------
  !Modules de subroutine utilitaires generiques => src/mod/subroutines_generiques.f90
  !-------------------------------------------------------------------------------
  use module_generique_desalloc
  use module_generique_initialise
  use module_generique_lecval
  use module_generique_mpimax
  use module_generique_mpimaxtab
  use module_generique_mpimin
  use module_generique_mpimintab
  use module_generique_mpisum
  !-------------------------------------------------------------------------------
  !Modules de subroutines => src/preparation.f90 et src/utilitaires_calcul.f90
  !-------------------------------------------------------------------------------
  use module_subroutines_preparation
  use module_subroutines_utilitaires_calcul
  !-------------------------------------------------------------------------------
  !Modules de subroutines de la librairie
  !-------------------------------------------------------------------------------
  use bib_module_caracteristiques_thermophys
  use bib_module_conditions_dom_scalaire
  use bib_module_connexions
  use bib_module_impression
  use bib_module_initialisation_vague
  use bib_module_interpolation_viscosite
  use bib_module_maillage
  use bib_module_metriques
  use bib_module_moyenne
  use bib_module_resolution_advection
  use bib_module_resolution_energie
  use bib_module_resolution_helmholtz
  use bib_module_resolution_navier
  use bib_module_resolution_transport
  use bib_module_resolution_turbulence
  use bib_module_reprise
  use bib_module_tests_d_arret
  use bib_module_viscosite_les
  use bib_module_tab_local2global
  !-------------------------------------------------------------------------------

  !===============================================================================
  !declarations des variables
  !===============================================================================
  implicit none

  !-------------------------------------------------------------------------------
  !maillage
  !-------------------------------------------------------------------------------
  integer, dimension(:  ),         allocatable :: typlt,typlk,typle,typlv2
  integer, dimension(:,:),         allocatable :: typlv,typlo,typlh
  !variables multibloc
  integer, dimension(:),           allocatable :: njr
  integer, dimension(:,:),         allocatable :: nbloc,lmg
  integer, dimension(:,:),         allocatable :: jbloc
  character(len=33), dimension(:), allocatable :: linom
  integer, dimension(:),           allocatable :: cnode,typnvb,typntb
  integer, dimension(:),           allocatable :: milieut,milieuv,milieuc
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !energie
  !-------------------------------------------------------------------------------
  integer, dimension(:),           allocatable :: jpent
  double precision,  dimension(:), allocatable :: tep,tp0,tp1,smc,smd,cal,cal0,com,&
       emi,tin,smcr,pro,ent,ent0,ent1,rcp,compliq,com0
  double precision,  dimension(:,:),   allocatable :: frs,frs0,bcoeft
  double precision,  dimension(:,:,:), allocatable :: cot,cot0
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !navier-stokes, is_darcy
  !-------------------------------------------------------------------------------
  integer, dimension(:),           allocatable :: jpenv,jpenp,loiper,loiperv
  double precision,  dimension(:), allocatable :: pre,div,rov,rov0,vts,vtset,qmv,slv,&
       pep,pep0,rho,rho0,rod,xit,beg,smv,smc0,vt0,vt1,vt2, &
       bip,pin,vim,vil,vik,vic,vir, &
       qm0,qm1,vin,gra,xiv,for,         &
       smcor,pr0,pr1,ssm,phi,phi0,phi1,por,pov,cpep
  double precision,  dimension(:,:),   allocatable :: bcoefv
  double precision,  dimension(:,:,:), allocatable :: per,per0
  !-------------------------------------------------------------------------------
  !

  !-------------------------------------------------------------------------------
  !turbulence
  !-------------------------------------------------------------------------------
  !k-e, rng, v2-f
  integer, dimension(:),           allocatable :: jpenk,jpene,jpenv2
  double precision,  dimension(:), allocatable :: cin,dis,duv,duf,cin0,cin1,dis0,dis1, &
       duv0,duv1,din,ein,duvin,sml,smm,smn, &
       smo,smi,smj,vit
  double precision,  dimension(:,:), allocatable :: bcoefk,bcoefe,bcoefv2
  logical, dimension(:),             allocatable :: lobstur
  !statisitiques
  double precision,  dimension(:),   allocatable :: cinmy,vitmy,tepmy,premy,prep,tepp, &
       tepec,preec,prevar,tepvar,rhomy,&
       vtsp,vtsmy,vtsec,vtsvar
  double precision,  dimension(:,:), allocatable :: vormy,vorec,vorvar,cisvar,vtpvar,&
       vorp,coumy,qdmmy
  double precision,  dimension(:,:), allocatable :: tempsm_nc,tempsm0_nc
  double precision,  dimension(:,:), allocatable :: vtsmyf,vtspf,vtsvard0,vtsvard1
  double precision,  dimension(:,:), allocatable :: vtsvarf,vtsecf,vtsecd1,vtsecd0
  double precision,  dimension(:,:), allocatable :: vtsmyd1,vtsmyd0,vtspd1,vtspd0
  double precision,  dimension(:,:), allocatable :: cisvard1,cisvard0,cisvarf
  !inra
  double precision,  dimension(:),   allocatable :: z0,hd,tsurf,hc,cd,cda,zdc,resusp,  &
       mais,orge,ble,chaume,tetaf,af                            
  !lim
  double precision,  dimension(:,:,:), allocatable :: bfidi,rx,ry,rz
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !advection
  !-------------------------------------------------------------------------------
  integer,           dimension(:),       allocatable :: jpeno,ftc,interpar
  double precision,  dimension(:),       allocatable :: ret
  double precision,  dimension(:,:),     allocatable :: cou,boo,ton,dist,dism,bik,kin
  double precision,  dimension(:),       allocatable :: mou
  double precision,  dimension(:,:,:,:), allocatable :: ftp 
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !transport
  !-------------------------------------------------------------------------------
  integer,           dimension(:),     allocatable :: jpenh
  double precision,  dimension(:),     allocatable :: mas,mas0
  double precision,  dimension(:,:),   allocatable :: esp,esp0,esp1,smp,smq,difc,dfvc,hin
  double precision,  dimension(:,:),   allocatable :: esph,btrh,ttrh
  double precision,  dimension(:,:,:), allocatable :: bcoefh
  !inra
  double precision,  dimension(:,:),   allocatable :: ssimpmais,sssedmais,ssimp,sssed,dts
  double precision,  dimension(:,:),   allocatable :: vch,vdp
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !helmholtz
  !-------------------------------------------------------------------------------
  integer,           dimension(:), allocatable :: jpenz
  double precision,  dimension(:), allocatable :: bil,lin,psi,smb,smt,dic
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! vof - sm
  !-------------------------------------------------------------------------------
  integer,           dimension(:),   allocatable :: ktp,kte,ktt,cori2,cork2,corj2
  double precision,  dimension(:),   allocatable :: wee,foe
  double precision,  dimension(:,:), allocatable :: mtp,mte,mtt,wep,fop,wet,fot
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !utilitaires
  !-------------------------------------------------------------------------------
  double precision,  dimension(:),   allocatable :: volf
  double precision,  dimension(:,:), allocatable :: tau,couini
  !-------------------------------------------------------------------------------

  !===============================================================================
  !variables temporaires de travail
  !===============================================================================
  double precision :: tpotv,tetotv0,tecinv0,tepotv0
  integer :: ltp,ltn,lvp,ic,l,i,j,k,n,lvb,lvf,lmp,lvn,e
  double precision :: rho2,kprime,taup,tdeb,tfin,tfin2
  character(len=1) :: a
  logical :: is_util_vague
  double precision,  dimension(:,:), allocatable :: vtspp
  double precision :: norv
  !-------------------------------------------------------------------------------

  ! Deewakar Variables
  ! How to give R by reading from .don file
	
  double precision :: Gas_constant_air
  Gas_constant_air=287.D0

  !===============================================================================
  !fin des declarations
  !===============================================================================

  !-------------------------------------------------------------------------------
  if (is_debug) write(*,*) 'entree calcul'
  !-------------------------------------------------------------------------------

  !===============================================================================
  !initialisation des structures et des variables scalaires definies dans les modules
  !de src/mod;
  !couplage avec l'interface graphique (lecture du fichier .thetis) ou bien
  !couplage avec l'interface texte (lecture du fichier .don)
  !la subroutine 'initialisation_structures' se trouve dans le fichier thetis.f90
  !===============================================================================
  call initialisation_structures
  !===============================================================================

  if (utilitaires%is_test_perf) tdeb=mpi_wtime()

  !===============================================================================
  !creation du maillage (calcul des coordonnées des noeuds)
  !===============================================================================
  call maillage (nbloc,lmg,linom,cnode,njr,jbloc,milieut,mpips,mpivs)
  !===============================================================================

  !===============================================================================
  ! pour la lecture d'un mot clé utilisateur, décommenter la ligne suivante :
  ! call lecval('mon_motcle',var)
  ! la dupliquer pour la lecture de plusieurs mots-clés.
  ! mon_motcle est un mot-clé choisi par l'utilisateur, présent dans le
  ! fichier de données, precedé de user et suivi d'un réel, d'un entier ou de oui/non :
  ! user mon_motcle 1.d0
  ! var doit être déclaré dans le bloc de déclaration de la subroutine calcul, être
  ! de type double precision, integer ou logical selon le type de la valeur lue. apres lecture
  ! du fichier de données, var contiendra la valeur du réel, de l'entier ou aura
  ! pour valeur .true. ou .false.
  ! exemple
  !===============================================================================



  !===============================================================================
  !calcul des matrices de connexions pour les différentes grilles
  !===============================================================================
  call connexions (jcoeft,jcoefv,jcoefm,jcoefc,jelem,jelec,jpenv,jpenp,jpent,jpenk,&
       jpene,jpenv2,jpenh,jpeno,jpenz,jpert,jperv,lmg,indvl,ilimt,&
       jlimt,ilimt2,jlimt2,typnt,klt,cnode,ktj,kvj,kmj,kcj,kro, &
       krm,krj,nbloc,jnode)
  !===============================================================================

  if (utilitaires%is_test_perf) then
     tfin=mpi_wtime()
     tfin=tfin - tdeb
     call mpimax(tfin)
     if (rang==0) print*,'tps ini',tfin
  endif

  !===============================================================================
  !calcul des métriques
  !===============================================================================
  call metriques(det,dev,dem,nbloc,cnode)
  !===============================================================================


  !===============================================================================
  !maillage grille decalee
  !===============================================================================
  call initialise(milieuv,kkv)
  call initialise(mtv,kkv,3)
  call coordonnees_grille_v(mtv,mtg)
  !===============================================================================


  !===============================================================================
  !preparation des équations : allocation des tableaux, conditions initiales,
  !conditions aux limites, penalisation, etc.
  !===============================================================================

  !-------------------------------------------------------------------------------
  !allocation de tableaux si navier-stokes n'est pas résolu
  !-------------------------------------------------------------------------------
  if (.not.is_navier .and. .not.is_darcy) then
     call initialise(vts,kkv)
     call initialise(vt0,kkv) 
     call initialise(vim,kkt,fluide%viscosite)
     call initialise(pep,kkt)
     call initialise(div,kkt)            
     call initialise(rho,kkt,fluide%masse_vol)
     call initialise(rov,kkv,fluide%masse_vol)
     if (is_obstacle.or.is_poreux) then
        call initialise(rho0,kkt,fluide%masse_vol)
        call initialise(rov0,kkv,fluide%masse_vol)
        call initialise(pep0,kkt)
        call initialise(per,kkv,1,1)
        call initialise(per0,kkv,1,1)
        call initialise(loiper,kkt,val=0)
        call initialise(loiperv,kkv,val=0)
     endif
     call initialise(pre,kkt)
     call initialise(xit,kkt,fluide%coef_compressibilite)
     call initialise(xiv,kkv,fluide%coef_compressibilite)
     if (is_compressible) call initialise(rod,kkt)
  endif
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  !préparation de navier-stokes, is_darcy
  !-------------------------------------------------------------------------------
  if (is_navier.or.is_darcy) then
     call preparation_navier (pre,vts,vtset,qmv,div,rho,rho0,rov,rov0,rod,smv,smc0,&
          slv,vt0,vt1,vt2,qm0,&
          qm1,pep,pep0,por,pov,cpep,loiper,loiperv,vim,xit,xiv,bcoefv,vin,&
          bip,pin,per,per0,vil,vik,vic,vir,                       &
          tep,beg,pro,for,gra,smcor,pr0,pr1,ssm,                  &
          tau,phi,phi0,phi1,milieut,milieuv,jpenv,                &
          jpenp,linom,lmg,nbloc,typlv)
  endif
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  !préparation de l'équation de l'énergie
  !-------------------------------------------------------------------------------
  if (is_energie) then
     call preparation_energie (tep,smc,smd,tp0,tp1,cal,cal0,com,com0,beg,ent,ent0,&
          ent1,vts,rho,rho0,rov,rov0,rcp,tin,frs,frs0,compliq,kte,mte,foe,wee,  &
          tsurf,cou,cot,cot0,pro,bcoeft,milieut,milieuv,jpent,typlt,&
          linom,lmg,nbloc)
  endif
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  !préparation de l'advection
  !-------------------------------------------------------------------------------
  if (is_advection) then
     call preparation_advection (cou,dist,dism,mou,boo,ton,bik,kin,ret,vts,pep,&
          milieut,milieuc,milieuv,interpar,jpeno,linom,lmg,nbloc,ftp,ftc,ktp,mtp,&
          wep,fop,typlo,mpiad,mpire)
  endif
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  !préparation des modèles de turbulence rans et les
  !-------------------------------------------------------------------------------
  if (is_turbulence) then
     call preparation_turbulence (cin,dis,duv,duf,cin0,cin1,dis0,dis1,duv0,duv1,&
          din,ein,duvin,sml,smm,smn,smo,smi,smj,vit,z0,hd,hc,cd,cda,zdc,resusp,&
          mais,orge,ble,chaume,tetaf,af,bcoefk,bcoefe,bcoefv2,jpenk,jpene,jpenv2,&
          milieut,lobstur,linom,lmg,nbloc,typlk,typle,typlv2,typlv,vin)
  endif
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  !préparation de l'équation de transport
  !-------------------------------------------------------------------------------
  if (is_transport) then
     call preparation_transport (esp,esph,smp,smq,esp0,esp1,hin,difc,dfvc,vts,vch,&
          vdp,ktt,mtt,wet,fot,ssimpmais,sssedmais,ssimp,sssed,dts,milieut,bcoefh,&
          jpenh,btrh,ttrh,linom,lmg,nbloc,typlh)
  endif
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  !préparation de is_helmholtz
  !-------------------------------------------------------------------------------
  if (is_helmholtz) then
     call preparation_helmholtz (psi,smb,smt,dic,bil,lin,jpenz)
  endif
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  !préparation des domaines fictifs
  !-------------------------------------------------------------------------------
  if (psm%pent.or.psm%penv) then
     call predomfic (jcoeft,jcoefv,ktj,kvj,kkt,kkv,ndim,is_debug)
  endif
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  !tableaux  de moyennes, écarts-types et variances
  !-------------------------------------------------------------------------------
  if (statistiques%is_statistiques) &
       call initialisation_moyennes (cinmy,vitmy,tepmy,premy,coumy,qdmmy,rhomy,vtsmy,&
       vormy,tepec,preec,vtsec,vorec,&
       prevar,tepvar,vtsvar,vorvar,cisvar,vtpvar,vtsp,vorp,tepp,prep, &
       vtsmyd1,vtsmyd0,vtsmyf,&
       vtspd1,vtspd0,vtspf,   &
       vtsvard1,vtsvard0,vtsvarf,vtsecf,vtsecd1,vtsecd0,&
       cisvard1,cisvard0,cisvarf,tempsm_nc,tempsm0_nc)
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  !autres initialisations
  !-------------------------------------------------------------------------------
  call initialisation_utilitaires(vts,vt0,vt1,rho,cou,dism,esp,cin,cin0,cin1,&
       dis,dis0,dis1,tep,tp0,tp1,z0,hd,tsurf,hc,cd,cda,zdc,volf,  &
       mais,orge,ble,chaume,tetaf,af,vin,ein,din,tin,&
       rod,xiv,xit,pre,pr0,rov,vim,div,pep,smb,& 
       jpenv,jpene,jpenk,jpent,rho2,kprime,mas,mas0,&
       bcoefv,bcoefe,bcoefk,bcoeft,couini,bil,lin,jpenz)
  !-------------------------------------------------------------------------------


  !===============================================================================
  !fin de la preparation des équations
  !===============================================================================


  !===============================================================================
  !lecture du fichier de reprise
  !===============================================================================
  if (reprise%is_lecture) then
     call reprise_lecture (vts,vt0,vt1,pre,div,pr0,pr1,tep,tp0,tp1,cin,cin0,cin1,dis,dis0,&
          dis1,duv,duv0,duv1,duf,vit,dist,dism,ftp,interpar,esp,esp0, & 
          esp1,frs,cou,cinmy,vitmy,tepmy,premy,vtsmy,vormy,coumy,&
          qdmmy,rhomy,vtsvar,cisvar,vtpvar,tepvar,prevar,vorvar, &
          tempsm_nc,tempsm0_nc,vtsmyf,vtsvard0,vtsvard1,vtsvarf, &
          vtsecf,vtsecd1,vtsecd0,vtsmyd1,vtsmyd0,cisvard1,cisvard0,&
          cisvarf,smc,smd,smv,slv,smp,smq,sml,smn,smm,smo,smi,smj, &
          ktp,kte,ktt,wee,foe,mtp,mte,mtt,wep,fop,wet,fot,ent,ent0,ent1,&
          phi,phi0,phi1,vtset)
  endif
  !===============================================================================


  !-------------------------------------------------------------------------------
  !calcul des caracteristiques thermophysiques de l'ecoulement
  !-------------------------------------------------------------------------------
  if (is_turbulence.and.turbulence_modele==turbulence_modele_les) &
       call viscosite_les(vit,vts,vim,cou,pro,turbulence,navier,mpips,typnt,&
       typlv,jpent,jpenv)

  ! Deewakar : All the thermophysical properties are calculated here

  call caracteristiques_thermophys(pre,tep,rho,rho0,rov,rov0,vim,vil,vik,vic,vir,&
       esp,com,com0,cot,cot0,cal,cal0,xit,xiv,beg,rcp,cou,dism,bik,kin,boo,ton,pep, &
       pep0,per,per0,por,pov,cpep,loiper,loiperv,vit,ret,frs,&
       compliq,difc,dfvc,milieut,milieuv,milieuc,jpeno)
  !initialisation de la pression
  !=> ici plutot que dans prens (prepa.f90) car peut dependre de la masse volumique
  if (iteration_temps==0.and.is_navier) then
     call conditions_domaine_scalaire('initialise','pression',1,kkt,mtg,tab=pre,&
          cond_dom=navier%cond_ini_pression,masse_vol=rho)
     pr0=pre
  endif
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  !autres initialisations 2
  !-------------------------------------------------------------------------------
  call initialisation_utilitaires_2(vts,vt0,vt1,rho,cou,dism,esp,cin,cin0,cin1,&
       dis,dis0,dis1,tep,tp0,tp1,z0,hd,tsurf,hc,cd,cda,zdc,volf,  &
       mais,orge,ble,chaume,tetaf,af,vin,ein,din,tin,&
       rod,xiv,xit,pre,pr0,rov,vim,div,pep,smb,& 
       jpenv,jpene,jpenk,jpent,rho2,kprime,mas,mas0,&
       bcoefv,bcoefe,bcoefk,bcoeft,couini,bil,lin,jpenz)
  !-------------------------------------------------------------------------------



  !===============================================================================
  !impression des conditions initiales
  !===============================================================================
  if (impressions%is_impressions) then
!     call mpi_barrier(comut,code)
!     tdeb=mpi_wtime()

     call impression(vts,pre,tep,rho,div,cin,dis,duv,duf,vit,vim,pep,por,&
          rcp,cou,tin,din,ein,duvin,com,ent,dist,dism,&
          cal,smc,smv,esp,frs,esph,psi,&
          gra,pro,cda,ssimp,sssed,cinmy,vitmy,tepmy,premy,&
          vtsmy,vormy,coumy,tepec,preec,vtsec,vorec,&
          prevar,tepvar,vtsvar,vtsp,vorvar,cisvar,vtpvar,&
          vtsmyf,vtsmyd1,vtsmyd0,vtsvard0,vtsvard1,vtsvarf,vtsecf,&
          vtsecd1,vtsecd0,cisvard1,cisvard0,cisvarf)

!     tfin=mpi_wtime()
!     tfin2=tfin - tdeb
!     call mpimax(tfin2)
!     tfin2=dble(kktg*0.000044)/tfin2
!     if (rang==0) print*,'tps impressions min_SAR',tfin2,dble(kktg*0.000044) !=> min_SAR
!     tfin=tfin - tdeb
!     tfin=dble(kkp*4*11/1000000)/tfin
!     call mpisum(tfin) !=>max_SAR
!     if (rang==0) print*,'tps impressions max_SAR',tfin,dble(kkp*4*11/1000000) !=> min_SAR
     
     ! impression possible de tableaux supplémentaires de dimension kkt (jusqu'à 8),
     ! en rajoutant après le dernier argument usr1=tab1,usr2=tab2,usr3=tab3,...,
     ! ou tab1, tab2... sont les tableaux de dimension kkt a imprimer. il faut rajouter
     ! imprime usr1 nom1, imprime usr2 nom2, etc. dans le fichier de donnees
  endif
  !===============================================================================

!!$  call initialise(vtspp,kkt,3)
!!$  call intervitesse (vts,vtspp)
!!$  norv=0.d0
!!$  do ltp=1,kkt
!!$   norv=max(norv,sqrt(vtspp(ltp,1)**2  + vtspp(ltp,2)**2 + vtspp(ltp,3)**2))
!!$  enddo
!!$  call mpimax(norv)
!!$  if (rang==0) print*,'norv',norv 
!!$  !-------------------------------------------------------------------------------



  !*******************************************************************************
  !boucle en temps   
  ! Deewakar : Begining of iteration in time
  !*******************************************************************************
  do while (iteration_temps<iterations_temps_max)
     !*******************************************************************************

  !===============================================================================
  if (utilitaires%is_test_perf.and.iteration_temps==1) then !en fait a partir de la seconde iter en temps
     temps_hypresca=0.d0
     temps_hyprevec=0.d0
     temps_thetis=0.d0
     call mpi_barrier(comut,code)
     temps_thetis_debut=mpi_wtime()
  endif
  !===============================================================================

     !-------------------------------------------------------------------------------
     !preparation temps
     !-------------------------------------------------------------------------------
     call preparation_temps(vts,vt0,rho,cou,cin,z0,hd,tsurf,hc,cda,zdc,vin,ein,din,&
          tin,jpenv,jpene,jpenk,jpent,jpenz,kprime,bcoefv,bcoefe,&
          bcoefk,bcoeft,cal,com,taup,rov,smv,vit,slv,ssm,&
          smb,bil,lin,typlv,dism,pre,vim,vtsmy,premy)
     !-------------------------------------------------------------------------------


     !===============================================================================
     !résolution des équations
     !===============================================================================

     !-------------------------------------------------------------------------------
     !résolution de navier-stokes
     !-------------------------------------------------------------------------------
     if ((is_navier.or.is_darcy) .and. &
          iteration_temps>=navier%activation%debut .and. &
          iteration_temps<=navier%activation%fin) then
        !-------------------------------------------------------------------------------
        call resolution_navier (pre,vts,vtset,div,cin,rho,cou,dism,slv,smv,smc0,vin,bip,       &
             pin,vt0,vt1,vt2,pr0,pr1,ssm,vim,vil,vik,vic,vir,xit,&
             beg,per,rov,for,phi,phi0,phi1,pov,bcoefv,       &
             navier,energie,turbulence,advection,ft,ftp,gra,smcor,frs,&
             jpenv,jpenp,typlv,typntb,typnvb,nbloc,&
             mpins,mpips,interpar,lg,vl,milieuv,mumps_par_ns,mumps_par_ps,tep,cal,cot)
        !-------------------------------------------------------------------------------
        if (is_navier) call imprim_param(imp,'navier')
        if (is_darcy)  call imprim_param(imp,'darcy')
        !-------------------------------------------------------------------------------
     endif
     !------------------------------------------------------------------------------- 

!!$  call intervitesse (vts,vtspp)
!!$  norv=0.d0
!!$  do ltp=1,kkt
!!$   norv=max(norv,sqrt(vtspp(ltp,1)**2  + vtspp(ltp,2)**2 + vtspp(ltp,3)**2))
!!$  enddo
!!$  call mpimax(norv)
!!$  if (rang==0) print*,'norv',norv


     !-------------------------------------------------------------------------------
     !résolution de l'équation de l'énergie
     !-------------------------------------------------------------------------------
     if (is_energie .and.&
          iteration_temps>=energie%activation%debut .and. &
          iteration_temps<=energie%activation%fin) then
        !-------------------------------------------------------------------------------
        call resolution_energie(tep,rho,vts,smc,smd,smcr,tp0,tp1,tin,cal,ent,ent0,ent1, &
             vim,pep,rcp,rov,com,cot,frs,cou,esp,compliq,bcoeft,jpent,typlt,&
             jpert,vt0,mte,kte,wee,foe,energie,navier,advection,&
             mpien,mumps_par_en,beg,xit)
        call imprim_param(imp,'energie')
        !-------------------------------------------------------------------------------
     endif
     !-------------------------------------------------------------------------------


     !-------------------------------------------------------------------------------
     !résolution des modèles de turbulence rans
     !-------------------------------------------------------------------------------
     if (is_turbulence .and. &
          .not.turbulence_modele==turbulence_modele_les .and. &
          iteration_temps>=turbulence%activation%debut .and. &
          iteration_temps<=turbulence%activation%fin) then
        !-------------------------------------------------------------------------------
        call resolution_turbulence(cin,dis,duv,duf,cin0,cin1,dis0,dis1,duv0,duv1,din,ein,duvin,&
             sml,smm,smn,smo,smi,smj,hd,lobstur,vts,tep,rho,rov,cal,beg, &
             pep,vim,vit,gra,cda,bcoefk,bcoefe,bcoefv2,jpert,jpenk,jpene,&
             jpenv2,mpitu,mpips,canop,turbulence,navier,energie,mumps_par_ke,mumps_par_v2,   &
             mumps_par_kl)
        call imprim_param(imp,'turbulence')
        !-------------------------------------------------------------------------------
     endif
     !-------------------------------------------------------------------------------

     !-------------------------------------------------------------------------------
     !résolution de l'equation d'advection
     !-------------------------------------------------------------------------------
     if (is_advection .and. &
          iteration_temps>=advection%activation%debut .and. &
          iteration_temps<=advection%activation%fin) then
        !-------------------------------------------------------------------------------
        call resolution_advection(cou,dist,vts,vt0,boo,ton,ret,ftp,pov,jpeno,ktp,wep,&
             fop,mtp,advection,lg,vl,ft,mpiad,interpar)
        call imprim_param(imp,'advection')
        !-------------------------------------------------------------------------------
     endif
     !-------------------------------------------------------------------------------


     !-------------------------------------------------------------------------------
     !résolution de l'équation de transport
     !-------------------------------------------------------------------------------
     if (is_transport .and. &
          iteration_temps>=advection%activation%debut .and. &
          iteration_temps<=advection%activation%fin) then 
        !-------------------------------------------------------------------------------
        call resolution_transport (esp,esph,vts,smp,smq,esp0,esp1,difc,dfvc,bcoefh,hin,btrh ,  &
             ttrh,vit,tep,vch,vdp,lobstur,rho,ktt,mtt,wet,fot,vt0,       &
             transport,turbulence,ch,energie,jpert,nbloc,vim,jpenh,typlh,&
             mpitr,mpins,hc,ssimp,sssed,ssimpmais,sssedmais,resusp,mais,&
             ble,chaume,dts,canop,cin,dis,compliq,por,mumps_par_tr) 
        call imprim_param(imp,'transport')
     endif
     !-------------------------------------------------------------------------------


     !-------------------------------------------------------------------------------
     !résolution de l'équation d'helmholtz
     !-------------------------------------------------------------------------------
     if (is_helmholtz) then
        call resolution_helmholtz(psi,dic,smb,smt,bil,lin,det,dev,helmholtz,mpihe,&
             mumps_par_he,cort,jcoeft,jpenz,jpert,ktj,kkt,bzp,kkv,ndim,is_debug) 
        call imprim_param(imp,'helmholtz')
     endif
     !-------------------------------------------------------------------------------


     !===============================================================================
     !fin de la résolution des équations
     !===============================================================================


     !===============================================================================
     !calcul des caractéristiques (loi d'état, rhéologie, interpolations...)
     !===============================================================================
     if (is_turbulence.and.turbulence_modele==turbulence_modele_les) &
          call viscosite_les (vit,vts,vim,cou,pro,turbulence,navier,mpips,typnt,&
          typlv,jpent,jpenv)
     !-------------------------------------------------------------------------------
     ! Deewakar Comment : Thermophysical properties are calculated for next iteration here
     call caracteristiques_thermophys(pre,tep,rho,rho0,rov,rov0,vim,vil,vik,vic,vir,&
       esp,com,com0,cot,cot0,cal,cal0,xit,xiv,beg,rcp,cou,dism,bik,kin,boo,ton,pep, &
       pep0,per,per0,por,pov,cpep,loiper,loiperv,vit,ret,frs,&
       compliq,difc,dfvc,milieut,milieuv,milieuc,jpeno)
     !===============================================================================
     
     ! Deewakar Equations Modified
      !rho=(pre+101325.)/(Gas_constant_air*tep)
      !beg=1.0/tep
     ! xit=1/(pre+101325.)

     !===============================================================================
     !utilitaires
     !===============================================================================

     !-------------------------------------------------------------------------------
     !calcul des moyennes, écarts-types et variances
     !-------------------------------------------------------------------------------
     if (statistiques%is_statistiques) then
        call moyenne (vit,tep,vts,pre,cou,rho, &
             tempsm_nc,tempsm0_nc,&
             cinmy,vitmy,tepmy,premy,vtsmy,vormy,coumy,qdmmy,rhomy, &
             tepec,preec,vtsec,vorec, &
             prevar,tepvar,vtsvar,vorvar,cisvar,vtpvar, &
             vtsmyd1,vtsmyd0,vtsmyf,&
             vtspd1,vtspd0,vtspf,&
             vtsvard1,vtsvard0,vtsvarf,vtsecf,vtsecd1,vtsecd0,&
             cisvard1,cisvard0,cisvarf,&
             vtsp,vorp,tepp,prep,navier%pas_de_temps,energie%pas_de_temps)
     endif
     !-------------------------------------------------------------------------------


     !-------------------------------------------------------------------------------
     !utilitaires et validation
     !-------------------------------------------------------------------------------
     call calcul_utilitaires (cou,dism,frs,esp,vts,rho,cal,com,vin,tep,psi,cin,dis,duv,vt0,&
          vt1,smv,pep,pre,vim,mas0,milieut,jpent,jpenv,rho2,&
          difc,hin,volf,ssimp,sssed,vdp,vch,smq,dts,hc,vit,vtsmy,premy,&
          couini,tin)
     !-------------------------------------------------------------------------------


     !===============================================================================
     !test d'arrêt
     !===============================================================================
     if ( mod(iteration_temps,tests_arret%frequence)==0 .and. tests_arret%is_test) then
        call tests_d_arret(vt0,vts,tep,tp0,cin,cin0,dis,dis0,duv,duv0,esp,esp0,div)
     endif
     !===============================================================================



     !===============================================================================
     !impressions
     !===============================================================================
     if (impressions%is_impressions .or. tests_arret%is_arret) then
        !-------------------------------------------------------------------------------
        call impression(vts,pre,tep,rho,div,cin,dis,duv,duf,vit,vim,pep,por,      &
             rcp,cou,tin,din,ein,duvin,com,ent,dist,dism,                   &
             cal,smc,smv,esp,frs,esph,psi, &
             gra,pro,cda,ssimp,sssed,cinmy,vitmy,tepmy,premy,          &
             vtsmy,vormy,coumy,tepec,preec,vtsec,vorec,                &
             prevar,tepvar,vtsvar,vtsp,vorvar,cisvar,vtpvar,        &
             vtsmyf,vtsmyd1,vtsmyd0,vtsvard0,vtsvard1,vtsvarf,vtsecf,  &
             vtsecd1,vtsecd0,cisvard1,cisvard0,cisvarf)
        ! impression possible de tableaux supplémentaires de dimension kkt (jusqu'à 8),
        ! en rajoutant après le dernier argument usr1=tab1,usr2=tab2,usr3=tab3,...,
        ! ou tab1, tab2... sont les tableaux de dimension kkt a imprimer. il faut rajouter
        ! imprime usr1 nom1, imprime usr2 nom2, etc. dans le fichier de donnees
     endif
     !===============================================================================



     !===============================================================================
     !pour la validation des cas test => impression dans 'verification.out'
     if (is_validation) call validation
     !===============================================================================


     !-------------------------------------------------------------------------------
     !passage tn=>tn+1
     !-------------------------------------------------------------------------------
     call passage_temps_n_plus_1(vts,vt0,vt1,vt2,pre,pr0,pr1,tep,tp0,tp1,esp,esp0,esp1,&
       cin,cin0,cin1,dis,dis0,dis1,duv,duv0,duv1,qmv,qm0,qm1,ent,ent0,ent1,frs,frs0,&
       phi,phi0,phi1)
     !-------------------------------------------------------------------------------


     !===============================================================================
     !ecriture du fichier de reprise
     !===============================================================================
     if (reprise%ecriture_freq/=0 .and.&
          ((reprise%is_ecriture.and.mod(iteration_temps,reprise%ecriture_freq)==0).or.&
          (reprise%is_ecriture.and.tests_arret%is_arret))) then 
        call reprise_ecriture (vts,vt0,vt1,pre,pr1,tep,tp0,tp1,cin,cin0,cin1,div,dis,&
             dis0,dis1,duv,duv0,duv1,duf,vit,dist,dism,ftp,interpar,esp,& 
             esp0,esp1,frs,cou,cinmy,vitmy,tepmy,premy,vtsmy,vormy,&
             coumy,qdmmy,rhomy,vtsvar,cisvar,vtpvar,tepvar,prevar, &
             vorvar,tempsm_nc,tempsm0_nc,vtsmyf,vtsvard0,vtsvard1, &
             vtsvarf,vtsecf,vtsecd1,vtsecd0,vtsmyd1,vtsmyd0,cisvard1, &
             cisvard0,cisvarf,smc,smd,smv,slv,smp,smq,sml,smn,smm,smo,&
             smi,smj,ktp,kte,ktt,wee,foe,mtp,mte,mtt,wep,fop,wet,fot,&
             ent,ent0,ent1,phi,phi0,phi1,vtset)
     endif
     !===============================================================================


     !===============================================================================
     !arret
     !===============================================================================
     if (tests_arret%is_arret) then
        call imprim_param(imp,'arret')
        exit
     endif
     !===============================================================================
     !call flush(imp)
     !call flush(ifih2)
     !===============================================================================



     !*******************************************************************************
     ! fin de la boucle en temps
     !*******************************************************************************
  enddo
  !*******************************************************************************
  ! ENd of Iteration

  !===============================================================================
  if (utilitaires%is_test_perf) then
!     temps_thetis_fin=mpi_wtime()
!     temps_thetis_fin=temps_thetis_fin - temps_thetis_debut
!     call mpimax(temps_thetis_fin)
!     temps_thetis_fin=temps_thetis_fin/real(iterations_temps_max-1)
!     temps_hyprevec_final=temps_hyprevec_final/real(iterations_temps_max-1)
!     temps_hypresca_final=temps_hypresca_final/real(iterations_temps_max-1)
!     if (rang==0) then
!        open(unit=1000,file='temps_thetis.txt',position='append')
        !write(1000,'(i8,5(1x,e13.5))') nb_proc, temps_thetis_fin, temps_hyprevec_final, temps_hypresca_final,&
        !     temps_thetis_fin - temps_hyprevec_final - temps_hypresca_final, tfin
!        write(1000,'(i8,3(1x,e13.5))') nb_proc,tfin,tfin2
!        close(1000)
!     endif
  endif
  !===============================================================================


  !===============================================================================
  !desallocation de la memoire
  !===============================================================================

  !-------------------------------------------------------------------------------
  !maillage
  !-------------------------------------------------------------------------------
  if (is_debug) print*,'deallocation maillage',rang
  call desalloc(mtg)    ; call desalloc(mtv)
  call desalloc(det)    ; call desalloc(dev)    ; call desalloc(dem)
  call desalloc(klt)    ; call desalloc(indvl)  ; call desalloc(ilimt) ; call desalloc(ilimt2)
  call desalloc(jlimt)  ; call desalloc(jlimt2) ; call desalloc(typnt)
  call desalloc(pdom)   ; call desalloc(pdomo)  ; call desalloc(cort)  ; call desalloc(cori)
  call desalloc(cork)   ; call desalloc(corj)
  call desalloc(jcoeft) ; call desalloc(ktj)
  call desalloc(corv)   ; call desalloc(kro)    ; call desalloc(jcoefv); call desalloc(kvj)
  call desalloc(corh)   ; call desalloc(jcoefc) ; call desalloc(kcj)   ; call desalloc(jelec)
  call desalloc(corm)   ; call desalloc(krm)    ; call desalloc(kmj);    call desalloc(jcoefm)
  call desalloc(jelem)  ; call desalloc(krj)
  call desalloc(core)   ; call desalloc(coreh)
  if (allocated(grille%sous_domaines_x1)) deallocate(grille%sous_domaines_x1)
  if (allocated(grille%sous_domaines_x2)) deallocate(grille%sous_domaines_x2)
  if (allocated(grille%sous_domaines_x3)) deallocate(grille%sous_domaines_x3)
  if (allocated(impressions%segments))    deallocate(impressions%segments)
  if (allocated(impressions%points))      deallocate(impressions%points)
  if (allocated(impressions%coord))       deallocate(impressions%coord)
  if (allocated(impressions%vscoef))      deallocate(impressions%vscoef)
  if (allocated(impressions%lscoef))      deallocate(impressions%lscoef)
  if (allocated(gages))       deallocate(gages)
  if (allocated(ecoord))      deallocate(ecoord)
  if (allocated(evscoef))     deallocate(evscoef)
  if (allocated(elscoef))     deallocate(elscoef)
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  if (nb_proc>1) then
     !-------------------------------------------------------------------------------
     if (is_debug) print*,'deallocation //',rang
     call desalloc(dims)  ; call desalloc(coords) ; call desalloc(vois)  ; call desalloc(lmn)
     call desalloc(mmn)   ; call desalloc(nmn)    ; call desalloc(ideb)  ; call desalloc(ifin)
     call desalloc(kdeb)  ; call desalloc(kfin)   ; call desalloc(jdeb)  ; call desalloc(jfin)
     call desalloc(kktm)  ; call desalloc(kkcm)   ; call desalloc(lml)   ; call desalloc(nml)
     call desalloc(mml)   ; call desalloc(idebh)  ; call desalloc(ifinh) ; call desalloc(kdebh) 
     call desalloc(kfinh) ; call desalloc(jdebh)  ; call desalloc(jfinh) ; call desalloc(lmc)
     call desalloc(nmc)   ; call desalloc(mmc)
     call desalloc(recvh) ; call desalloc(envvh)  ; call desalloc(nrvh)
     deallocate(mpins%rec,mpins%env,mpins%pin,mpins%dom,mpins%kic,mpins%irec,&
          mpins%ienv,mpins%prenv,mpins%prrec)
     deallocate(mpips%rec,mpips%env,mpips%pin,mpips%dom,mpips%kic,mpips%irec,&
          mpips%ienv,mpips%prenv,mpips%prrec)
     if (is_energie)   deallocate(mpien%rec,mpien%env,mpien%pin,mpien%dom,mpien%kic,mpien%irec,&
          mpien%ienv,mpien%prenv,mpien%prrec)
     if (is_transport) deallocate(mpitr%rec,mpitr%env,mpitr%pin,mpitr%dom,mpitr%kic,mpitr%irec,&
          mpitr%ienv,mpitr%prenv,mpitr%prrec)
     if (is_helmholtz) deallocate(mpihe%rec,mpihe%env,mpihe%pin,mpihe%dom,mpihe%kic,mpihe%irec,&
          mpihe%ienv,mpihe%prenv,mpihe%prrec)
     if (is_turbulence.and.(&
          turbulence_modele==turbulence_modele_ke.or.&
          turbulence_modele==turbulence_modele_rng.or.&
          turbulence_modele==turbulence_modele_v2f.or.&
          turbulence_modele==turbulence_modele_kl)) &
          deallocate(mpitu%rec,mpitu%env,mpitu%pin,mpitu%dom,mpitu%kic,mpitu%irec,&
          mpitu%ienv,mpitu%prenv,mpitu%prrec)
     if (is_advection) deallocate(mpiad%rec,mpiad%env,mpiad%irec,&
          mpiad%ienv,mpiad%prenv,mpiad%prrec)
     if (advection%is_regularisation) deallocate(mpire%rec,mpire%env,mpire%pin,mpire%dom,mpire%kic,mpire%irec,&
          mpire%ienv,mpire%prenv,mpire%prrec)
  endif
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !matrices de connexions
  !-------------------------------------------------------------------------------
  if (is_debug) print*,'deallocation matrice connexion',rang
  call desalloc(njr)     ; call desalloc(jpert)   ; call desalloc(jperv) ; call desalloc(nbloc)
  call desalloc(lmg)     ; call desalloc(typlt)   ; call desalloc(typlv) ; call desalloc(typlk)
  call desalloc(typle)   ; call desalloc(typlv2)  ; call desalloc(typlo) ; call desalloc(typlh)
  call desalloc(jnode)   ; call desalloc(jbloc)   ; call desalloc(linom) ; call desalloc(cnode)
  call desalloc(typnvb)  ; call desalloc(typntb)
  call desalloc(milieut) ; call desalloc(milieuv) ; call desalloc(milieuc)
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !caracteristiques des fluides et solides
  !-------------------------------------------------------------------------------
  if (is_debug) print*,'deallocation carac',rang
  if (allocated(fluides)) deallocate(fluides)
  if (allocated(solides)) deallocate(solides)
  if (allocated(especes)) deallocate(especes)
  if (allocated(phases))  deallocate(phases)
  if (allocated(tcesp))   deallocate(tcesp)
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !energie
  !-------------------------------------------------------------------------------
  if (is_debug) print*,'deallocation energie',rang
  call desalloc(jpent) ; call desalloc(tep)     ; call desalloc(tp0)    ; call desalloc(tp1)
  call desalloc(smc)   ; call desalloc(smd)     ; call desalloc(cal)    ; call desalloc(cal0)
  call desalloc(com)   ; call desalloc(emi)     ; call desalloc(tin)    ; call desalloc(smcr)
  call desalloc(pro)   ; call desalloc(ent)     ; call desalloc(ent0)   ; call desalloc(ent1)
  call desalloc(rcp)   ; call desalloc(compliq) ; call desalloc(com0)
  call desalloc(frs)   ; call desalloc(frs0)    ; call desalloc(bcoeft) 
  call desalloc(cot)   ; call desalloc(cot0)  
  if (is_energie) then
     if (is_debug) print*,'deallocation energie',rang
  endif
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !navier-stokes
  !-------------------------------------------------------------------------------
  if (is_debug) print*,'deallocation navier',rang
  call desalloc(jpenp)
  call desalloc(jpenv); call desalloc(pre);    call desalloc(div);   call desalloc(rov)
  call desalloc(vts);   call desalloc(vtset);  call desalloc(qmv);   call desalloc(slv)
  call desalloc(pep);   call desalloc(rho);    call desalloc(rod);   call desalloc(xit)
  call desalloc(beg);   call desalloc(smv);    call desalloc(vt0);   call desalloc(vt1)
  call desalloc(vt2);   call desalloc(smc0)
  call desalloc(bip);   call desalloc(pin);    call desalloc(vim)
  call desalloc(vil);   call desalloc(vik);    call desalloc(vic);   call desalloc(vir)
  call desalloc(qm0);   call desalloc(qm1);    call desalloc(vin);   call desalloc(gra)
  call desalloc(xiv);   call desalloc(for);    call desalloc(smcor); call desalloc(pr0)
  call desalloc(pr1)
  call desalloc(ssm);   call desalloc(bcoefv); call desalloc(per)
  call desalloc(pep0);  call desalloc(per0)  ; call desalloc(cpep)
  call desalloc(phi) ;  call desalloc(phi0);   call desalloc(phi1)
  call desalloc(pov) ;  call desalloc(rho0);   call desalloc(rov0) ; call desalloc(por) 
  call desalloc(loiper) ; call desalloc(loiperv)
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  !turbulent (rans)
  !-------------------------------------------------------------------------------
  if (is_debug) print*,'deallocation turbulent',rang
  call desalloc(jpenk);     call desalloc(jpene);      call desalloc(jpenv2);   call desalloc(cin)
  call desalloc(dis);       call desalloc(duv);        call desalloc(duf);      call desalloc(cin0)
  call desalloc(cin1);      call desalloc(dis0);       call desalloc(dis1);     call desalloc(duv0)
  call desalloc(duv1);      call desalloc(din);        call desalloc(ein);      call desalloc(duvin)
  call desalloc(sml);       call desalloc(smm);        call desalloc(smn);      call desalloc(smo)
  call desalloc(smi);       call desalloc(smj);        call desalloc(vit)
  call desalloc(bcoefk);    call desalloc(bcoefe);     call desalloc(bcoefv2)
  if (allocated(lobstur))  deallocate(lobstur)
  call desalloc(cinmy);      call desalloc(vitmy);    call desalloc(tepmy)
  call desalloc(premy);     call desalloc(prep);       call desalloc(tepp);     call desalloc(tepec)
  call desalloc(preec);     call desalloc(prevar);     call desalloc(tepvar);   call desalloc(rhomy)
  call desalloc(vtsp);      call desalloc(vtsmy);      call desalloc(vtsec);    call desalloc(vtsvar)
  call desalloc(vormy);     call desalloc(vorec);      call desalloc(vorvar);   call desalloc(cisvar)
  call desalloc(vtpvar);    call desalloc(vorp);       call desalloc(coumy);    call desalloc(qdmmy)
  call desalloc(tempsm_nc); call desalloc(tempsm0_nc); call desalloc(vtsmyf);   call desalloc(vtspf)
  call desalloc(vtsvard0);  call desalloc(vtsvard1);   call desalloc(vtsvarf);  call desalloc(vtsecf)
  call desalloc(vtsecd1);   call desalloc(vtsecd0);    call desalloc(vtsmyd1);  call desalloc(vtsmyd0)
  call desalloc(vtspd1);    call desalloc(vtspd0);     call desalloc(cisvard1); call desalloc(cisvard0)
  call desalloc(cisvarf);   call desalloc(z0);         call desalloc(hd);       call desalloc(tsurf)
  call desalloc(hc);        call desalloc(cd);         call desalloc(cda);      call desalloc(zdc)
  call desalloc(resusp);    call desalloc(mais);       call desalloc(orge);     call desalloc(ble)
  call desalloc(chaume);    call desalloc(tetaf);      call desalloc(af);       call desalloc(bfidi)
  call desalloc(rx);        call desalloc(ry);         call desalloc(rz)
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  !advection
  !-------------------------------------------------------------------------------
  if (is_debug) print*,'deallocation advection',rang
  call desalloc(jpeno) ; call desalloc(ftc) ; call desalloc(interpar) ; call desalloc(ret)
  call desalloc(cou)   ; call desalloc(boo) ; call desalloc(ton)      ; call desalloc(dist)
  call desalloc(dism)
  call desalloc(bik)   ; call desalloc(kin) ; call desalloc(mou)
  if (allocated(ftp)) deallocate(ftp)
   !-------------------------------------------------------------------------------



  !-------------------------------------------------------------------------------
  !transport
  !-------------------------------------------------------------------------------
  if (is_debug) print*,'deallocation transport',rang
  call desalloc(jpenh)     ; call desalloc(mas)  ; call desalloc(mas0) ; call desalloc(esp)
  call desalloc(esp0)      ; call desalloc(esp1) ; call desalloc(smp)  ; call desalloc(smq)
  call desalloc(difc)      ; call desalloc(dfvc) 
  call desalloc(hin)       ; call desalloc(esph)
  call desalloc(btrh)      ; call desalloc(ttrh)      ; call desalloc(bcoefh)
  call desalloc(ssimpmais) ; call desalloc(sssedmais) ; call desalloc(ssimp)  ; call desalloc(sssed)
  call desalloc(dts)       ; call desalloc(vch)       ; call desalloc(vdp)
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  !is_helmholtz
  !-------------------------------------------------------------------------------
  if (is_debug) print*,'deallocation helmholtz',rang
  call desalloc(bil); call desalloc(lin); call desalloc(psi); call desalloc(smb)
  call desalloc(smt); call desalloc(dic); call desalloc(jpenz)
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  ! vof - sm
  !-------------------------------------------------------------------------------
  if (is_debug) print*,'deallocation vof-sm',rang
  call desalloc(ktp); call desalloc(kte); call desalloc(ktt); call desalloc(wee)
  call desalloc(foe); call desalloc(mtp); call desalloc(mte); call desalloc(mtt)
  call desalloc(wep); call desalloc(fop); call desalloc(wet); call desalloc(fot)
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !vof-lag
  !-------------------------------------------------------------------------------
  if (is_debug) print*,'deallocation structures vof-lag',rang
  if (is_advection) then
     if (advection%equations(1)%methode==advection_methode_vof_lag) call endvl(vl)
  endif
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  !utilitaires
  !-------------------------------------------------------------------------------
  if (is_debug) print*,'deallocation utilitaires volume',rang
  if (allocated(utilitaires%volume_phases)) deallocate(utilitaires%volume_phases)
  if (allocated(utilitaires%volume_phases_ini)) deallocate(utilitaires%volume_phases_ini)
  if (allocated(utilitaires%erreur_volume_phases)) deallocate(utilitaires%erreur_volume_phases)
  call desalloc(tau)
  if (is_debug) print*,'deallocation objets',rang
  if (allocated(obj)) then
     do l=1,size(obj)
        if (allocated(obj(l)%lag))    deallocate(obj(l)%lag)
        if (allocated(obj(l)%lagtr))  deallocate(obj(l)%lagtr)
        if (allocated(obj(l)%lagtru)) deallocate(obj(l)%lagtru)
        if (allocated(obj(l)%lagtrw)) deallocate(obj(l)%lagtrw)
        if (allocated(obj(l)%lagtrv)) deallocate(obj(l)%lagtrv)
        if (allocated(obj(l)%hea))    deallocate(obj(l)%hea)
        if (allocated(obj(l)%heav))   deallocate(obj(l)%heav)
        if (allocated(obj(l)%dist))   deallocate(obj(l)%dist)
        if (allocated(obj(l)%distv))  deallocate(obj(l)%distv)
        if (allocated(obj(l)%mtgtr))  deallocate(obj(l)%mtgtr)
        if (allocated(obj(l)%mtvtr))  deallocate(obj(l)%mtvtr)
     enddo
     deallocate(obj) 
  endif
  if (allocated(par))     deallocate(par)    
  if (allocated(forc))    deallocate(forc)
  if (allocated(bdlim))   deallocate(bdlim)
  if (allocated(bdini))   deallocate(bdini)
  if (allocated(impressions%champs_user)) deallocate(impressions%champs_user)
  if (allocated(impressions%is_champs_user)) deallocate(impressions%is_champs_user)
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  !ns
  !-------------------------------------------------------------------------------
  if (is_debug) print*,'deallocation structures navier',rang
  if (allocated(navier%lagrangien_div)) deallocate(navier%lagrangien_div)
  if (allocated(navier%solveur_quantite_de_mvt_out)) &
       deallocate(navier%solveur_quantite_de_mvt_out)

  call desalloc_limite(vitesse=navier%cond_lim_vitesse)
  call desalloc_struct_partielles_domaine(navier%cond_ini_x1%partielles)
  call desalloc_struct_partielles_domaine(navier%cond_ini_x2%partielles)
  call desalloc_struct_partielles_domaine(navier%cond_ini_x3%partielles)
  call desalloc_struct_partielles_domaine(navier%cond_ini_pression%partielles)

  call desalloc_limite(dirichlet=navier%cond_lim_pression)
  call desalloc_struct_partielles_domaine(navier%penalisation_pression%partielles)

  !-------------------------------------------------------------------------------
  !en
  !-------------------------------------------------------------------------------
  if (is_debug) print*,'deallocation structures energie',rang
  call desalloc_limite(scalaire=energie%cond_lim)
  call desalloc_struct_partielles_domaine(energie%terme_source%partielles)
  call desalloc_struct_partielles_domaine(energie%terme_lineaire%partielles)
  call desalloc_struct_partielles_domaine(energie%cond_ini%partielles)
  call desalloc_struct_partielles_domaine(energie%penalisation%partielles)

  !-------------------------------------------------------------------------------
  !tr
  !-------------------------------------------------------------------------------
  if (is_debug) print*,'deallocation structures transport',rang
  if (allocated(transport%equations)) then
     do e=1,nesp
        call desalloc_limite(scalaire=transport%equations(e)%cond_lim)
        call desalloc_struct_partielles_domaine(transport%equations(e)%terme_source%partielles)
        call desalloc_struct_partielles_domaine(transport%equations(e)%terme_lineaire%partielles)
        call desalloc_struct_partielles_domaine(transport%equations(e)%cond_ini%partielles)
        call desalloc_struct_partielles_domaine(transport%equations(e)%penalisation%partielles)
     enddo
     deallocate(transport%equations)
  endif

  !-------------------------------------------------------------------------------
  !ad
  !-------------------------------------------------------------------------------
  if (is_debug) print*,'deallocation structures advection',rang
  if (allocated(advection%equations)) then
     do ic=1,nc
        call desalloc_limite(dirichlet=advection%equations(ic)%cond_lim)        
        call desalloc_struct_partielles_domaine(advection%equations(ic)%cond_ini%partielles)
        call desalloc_struct_partielles_domaine(advection%equations(ic)%penalisation%partielles)
     enddo
     deallocate(advection%equations)
  endif

  !-------------------------------------------------------------------------------
  !tu
  !-------------------------------------------------------------------------------
  if (is_debug) print*,'deallocation structures turbulence',rang
  call desalloc_limite(scalaire=turbulence%k%cond_lim)
  call desalloc_struct_partielles_domaine(turbulence%k%terme_source%partielles)
  call desalloc_struct_partielles_domaine(turbulence%k%terme_lineaire%partielles)
  call desalloc_struct_partielles_domaine(turbulence%k%cond_ini%partielles)
  call desalloc_struct_partielles_domaine(turbulence%k%penalisation%partielles)

  call desalloc_limite(scalaire=turbulence%epsilon%cond_lim)
  call desalloc_struct_partielles_domaine(turbulence%epsilon%terme_source%partielles)
  call desalloc_struct_partielles_domaine(turbulence%epsilon%terme_lineaire%partielles)
  call desalloc_struct_partielles_domaine(turbulence%epsilon%cond_ini%partielles)
  call desalloc_struct_partielles_domaine(turbulence%epsilon%penalisation%partielles)

  call desalloc_limite(scalaire=turbulence%v2%cond_lim)
  call desalloc_struct_partielles_domaine(turbulence%v2%terme_source%partielles)
  call desalloc_struct_partielles_domaine(turbulence%v2%terme_lineaire%partielles)
  call desalloc_struct_partielles_domaine(turbulence%v2%cond_ini%partielles)
  call desalloc_struct_partielles_domaine(turbulence%v2%penalisation%partielles)

  !-------------------------------------------------------------------------------
  !he
  !-------------------------------------------------------------------------------
  if (is_debug) print*,'deallocation structures helmholtz',rang
  call desalloc_limite(scalaire=helmholtz%cond_lim)
  call desalloc_struct_partielles_domaine(helmholtz%terme_source%partielles)
  call desalloc_struct_partielles_domaine(helmholtz%terme_lineaire%partielles)
  call desalloc_struct_partielles_domaine(helmholtz%cond_ini%partielles)
  call desalloc_struct_partielles_domaine(helmholtz%penalisation%partielles)

  !-------------------------------------------------------------------------------
  !mumps
  !-------------------------------------------------------------------------------
  call desalloc_mumps
  !-------------------------------------------------------------------------------

  !===============================================================================
  !fin de la désallocation de la memoire
  !===============================================================================
  if (is_debug) print*,'sortie calcul'
end subroutine calcul
!*******************************************************************************


