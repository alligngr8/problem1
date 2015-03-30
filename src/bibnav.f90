!###############################################################################

!                                navier-stokes

!###############################################################################

!===============================================================================
module bib_module_resolution_navier
contains
  !===============================================================================
  !*******************************************************************************
  subroutine resolution_navier (pre,vts,vtset,div,cin,rho,cou,dism,slv,smv,smc0,vin,bip,       &
       pin,vt0,vt1,vt2,pr0,pr1,ssm,vim,vil,vik,vic,vir,xit,&
       beg,per,rov,for,phi,phi0,phi1,pov,bcoefv,       &
       navier,energie,turbulence,advection,ft,ftp,gra,smcor,frs,&
       jpenv,jpenp,typlv,typntb,typnvb,nbloc,&
       mpins,mpips,interpar,lg,vl,milieuv,mumps_par_ns,mumps_par_ps,tep,cal,cot) 
    !*******************************************************************************
    use bib_module_cond_lim_libre_vitesse
    use bib_module_conexnav
    use bib_module_pen_nav
    use bib_module_prepavol
    use bib_module_projection_scalaire
    use bib_module_projection_vectorielle
    use bib_module_resol_mumps_dist
    use bib_module_resol_hypre_v
    use bib_module_resol_thetis
    use bib_module_tensions_superficielles
    use bib_module_pen_nav_sca
    use module_generique_mpimax
    !-------------------------------------------------------------------------------
    use module_mpif
    use module_commun
    use module_connexions
    use module_connexions_lim
    use module_coordonnees
    use module_grille
    use module_intercom
    use module_metriques
    use module_mpi
    use module_objet 
    use module_utilitaires, only : utilitaires
    use module_thermophysique
    !-------------------------------------------------------------------------------
    use module_struct_advection
    use module_struct_energie
    use module_struct_front_tracking
    use module_struct_mumps
    use module_struct_navier
    use module_struct_psm
    use module_struct_turbulence
    use module_struct_vof_lag
    !-------------------------------------------------------------------------------
    use module_generique_initialise
    use module_generique_mpisum
    use module_generique_mpimaxtab
    use module_generique_mpimintab
    !-------------------------------------------------------------------------------
    implicit none
    !-------------------------------------------------------------------------------
    !variables globales
    !-------------------------------------------------------------------------------
    integer, dimension(:),   allocatable   :: interpar
    integer, dimension(:,:), allocatable   :: nbloc
    !-------------------------------------------------------------------------------
    double precision,  dimension(:),   allocatable   :: pre,rho,div,pov
    double precision,  dimension(:),   allocatable   :: vim,vik,vic,vir
    double precision,  dimension(:),   allocatable   :: vts,vtset,rov,vt0
    double precision,  dimension(:),   allocatable   :: vin
    double precision,  dimension(:,:), allocatable   :: bcoefv,frs
    !-------------------------------------------------------------------------------
    integer, dimension(:), allocatable             :: jpenv,jpenp,typntb,typnvb,&
         milieuv
    integer, dimension(:,:), allocatable           :: typlv
    double precision,  dimension(:), allocatable   :: cin,xit,for,gra,vt1,vt2,slv,smv,&
         smc0,beg,vil,bip,pin,phi,phi0,phi1,smcor,pr0,pr1,ssm
    double precision,  dimension(:,:),     allocatable :: cou,dism
    double precision,  dimension(:,:,:,:), allocatable :: ftp
    double precision,  dimension(:,:,:),   allocatable :: per
    !-------------------------------------------------------------------------------
    type(struct_navier)                 :: navier
    type(struct_energie)                :: energie
    type(struct_turbulence_systeme)     :: turbulence
    type(struct_mpieq)                  :: mpins,mpips
    type(struct_advection_systeme)      :: advection
    type(struct_front_tracking)                        :: ft
    type(struct_vof_lag)                        :: lg
    type(struct_vof_lag_2)                       :: vl
    type(dmumps_struc)                  :: mumps_par_ns,mumps_par_ps
    !-------------------------------------------------------------------------------
    !declaration des variables locales
    !-------------------------------------------------------------------------------
    double precision, dimension(:,:,:), allocatable  :: mob
    integer                            :: n,n3,lvp
    !-------------------------------------------------------------------------------
    integer, dimension(:), allocatable     :: kic
    integer, dimension(:), allocatable     :: jcof,icof
    double precision,  dimension(:), allocatable     :: coef
    double precision,  dimension(:), allocatable     :: rcn,grc,grd1,smc,sol,tet,rep
    double precision,  dimension(:), allocatable     :: skn,rovb,rhob,frlv,gphi
    !-------------------------------------------------------------------------------
    integer                                      :: k,kdv,la,kt,npt,nps,ndd,nic,no
    double precision, dimension(kkv)             :: grd0
    double precision, dimension(ndim)            :: pm,pp
    double precision                             :: res,rho2,pvg,pvd
    !-------------------------------------------------------------------------------
    integer                            :: i,impp,lts,l,ltp
    integer                            :: kvi,kmi,kmv,kmm
    integer                            :: pout,proc,lax
	
	!-------------------------------------------------------------------------------
	! Deewakar variables
	double precision, dimension(:),    allocatable :: tep,cal,term1,numerator,denominator
	double precision, dimension(:),    allocatable :: add_smt,rep_new
	double precision, dimension(:),    allocatable :: diverge_V,diverge_flux,flux,grad_T
	double precision			       :: scalar
	double precision, dimension(:,:,:),allocatable  :: cot
	
    !-------------------------------------------------------------------------------
    character(len=100)                 :: fmt1
    !-------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    if (is_debug) write (*,*) 'entree resolution_navier'
    !-------------------------------------------------------------------------------

    !===============================================================================
    !boussinesq
    !===============================================================================
    if (navier%is_boussinesq) then
       allocate(rovb(kkv)) ; allocate(rhob(kkt))
       rovb=rov
       rhob=rho
       rov=fluide%masse_vol
       rho=fluide%masse_vol
    endif
    !===============================================================================


    !===============================================================================
    !calcul de la fraction liquide pour le changement de phase de type alliage
    !===============================================================================
    if (energie%is_changement_etat.and.energie%change_modele==change_modele_alliage) then
       call initialise(frlv,kkv)
       call interpolations(1.d0-frs(:,1),frlv,.true.,.false.)
    endif
    !===============================================================================


    !===============================================================================
    !calcul de la mobilite mu/k
    !===============================================================================
    if (navier%is_brinkman) then
       !-------------------------------------------------------------------------------
       allocate (mob(kkv,npd,npe))
       call permob (per,mob,vik,navier%is_iso,navier%is_ort,navier%is_ani)
       if (is_poreux)        where(milieuv==2) mob(:,1,1)=mob(:,1,1)*(pov + 1.d-20)
       if (energie%is_changement_etat.and.energie%change_modele==change_modele_alliage)&
            mob(:,1,1)=mob(:,1,1)*frlv
       !-------------------------------------------------------------------------------
    endif
    !===============================================================================




    !===============================================================================
    allocate(sol(kkv),smc(kkv))
    !-------------------------------------------------------------------------------
    if (schema_temps==schema_temps_gear.and.navier%methode_type==navier_methode_lagrangien)&
         allocate (grd1(kkv))
    !-------------------------------------------------------------------------------
    !schema theta
    allocate (tet(kkv))
    tet=0.5d0
    !===============================================================================



    !===============================================================================
    !calcul de grad(rho*k) si rans
    !===============================================================================
    if (is_turbulence .and. .not.turbulence_modele==turbulence_modele_les) then
       !-------------------------------------------------------------------------------
       allocate (rcn(kkt),grc(kkv))
       rcn=rho*cin
       call gradient(rcn,grc)
       deallocate(rcn)
    end if
    !===============================================================================


    !===============================================================================
    ! calcul de la tension superficielle
    !===============================================================================
    if (navier%is_tension_sup) then
       call tensions_superficielles(cou,skn,dism,ftp,rho,rov,navier,advection,ft,mpins,mpips)
    endif
    !===============================================================================


    !===============================================================================
    ! gestion des collisions pour les particules
    !===============================================================================
    if (allocated(advection%equations)) then
       if (advection%equations(1)%methode==advection_methode_vof_lag) then
          call sourcepartns (smv,cou,par,det,dev,mtg,mtv,kro,interpar,lg,kkt,kkv,navier%pas_de_temps,nc, &
               ndim,is_periodique_x1,is_periodique_x2,is_periodique_x3,is_debug)
       endif
    endif
    !===============================================================================



    !-------------------------------------------------------------------------------
    ! calcul des connectivites
    !-------------------------------------------------------------------------------
    call conexnav (coef,jcof,icof,jcoefv,kvj,kro,navier%schema,kkv,kdv,npt,nps,ndim,ndd)
    !-------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    ! calcul de kic
    !-------------------------------------------------------------------------------
    ! npt est l ordre du systeme lineaire
    ! nic est le nombre d'inconnues a resoudre
    ! kic est un tableau qui contient les inconnues a resoudre effectivement
    !-------------------------------------------------------------------------------
    npt=kkv
    allocate (kic(npt))
    !-------------------------------------------------------------------------------
    ! l utilisateur est cense remplir ici kic avec les valeurs des inconnues
    ! ou laisser kic a zero si l inconnue n est pas a resoudre
    !-------------------------------------------------------------------------------
    nic = npt
    kic = 0
    do lvp=1,npt
       kic(lvp)= lvp
    enddo
    !-------------------------------------------------------------------------------


    !-------------------------------------------------------------------------------
    ! lagrangien automatique
    !-------------------------------------------------------------------------------
    navier%is_viscosite_compression = .false.
    allocate(rep(kkt))
    rep = navier%lagrangien_param
    !-------------------------------------------------------------------------------


    !-------------------------------------------------------------------------------
    if (navier%schema==navier_schema_weno) then
       !-------------------------------------------------------------------------------
       ! schema weno inertie
       !-------------------------------------------------------------------------------
       ! call nsinerexp (vtw,dev,vin,bcoefv,kro,jpenv,bvp,nbfv,jcoefv,typlv,indvl, &
       !                       kvj,navier%pas_de_temps,kkv,ndim,bvpl,is_debug)
       call nsinerexp3 (vts,dev,vin,bcoefv,kro,jpenv,bvp,nbfv,jcoefv,typlv,indvl,kvj,navier%pas_de_temps, &
            kkv,ndim,bvpl,mpins,is_debug)
       !centre
       ! call nsinerexp4 (vts,dev,vin,bcoefv,kro,jpenv,bvp,nbfv,jcoefv, &
       !                        kvj,dv,kkv,ndim,is_debug)
       if (schema_temps==schema_temps_gear) then
          vt1=vt0
       endif
       vt0=vts 
       navier%is_udiv=.false.
       !-------------------------------------------------------------------------------
    endif
    !-------------------------------------------------------------------------------



    !===============================================================================
    !===============================================================================
    !equation de navier-stokes - lagrangien augmente ou etape de prediction
    !===============================================================================
    !===============================================================================

    !-------------------------------------------------------------------------------
    do la=1,navier%lagrangien_iterations
       !-------------------------------------------------------------------------------


       !-------------------------------------------------------------------------------
       ! calcul du second membre
       !-------------------------------------------------------------------------------
       smc=0.d0
       !-------------------------------------------------------------------------------
       !gravite
       if (navier%is_boussinesq) then
          if (navier%is_gravite.and.navier%pression_type==pression_type_motrice) smc=(rovb - fluide%masse_vol)*gra
          if (navier%is_gravite.and.navier%pression_type==pression_type_totale)  smc=rovb*gra
       else
          if (navier%is_gravite.and.navier%pression_type==pression_type_motrice) smc=(rov - fluide%masse_vol)*gra
          if (navier%is_gravite.and.navier%pression_type==pression_type_totale)  smc=rov*gra
       endif
       !-------------------------------------------------------------------------------
       !porosite
       if (is_poreux) where(milieuv==2) smc=smc*pov
       if (energie%is_changement_etat.and.energie%change_modele==change_modele_alliage) &
            smc=smc*frlv
       !-------------------------------------------------------------------------------
       !autre termes sources
       if (navier%is_terme_source) smc=smc + smv
       if (navier%is_tension_sup) smc=smc + skn
       if (navier%is_coriolis) smc=smc - rov*smcor
       !-------------------------------------------------------------------------------
       !turbulent RANS
       if (is_turbulence .and. .not.turbulence_modele==turbulence_modele_les) then
          smc = smc - 2.d0/3.d0*grc
       endif
       !-------------------------------------------------------------------------------

       
       !===============================================================================
       !Correction de vitesse ; projection scalaire
       !===============================================================================
       !-------------------------------------------------------------------------------
       if (navier%is_correction_vitesse) then
          !-------------------------------------------------------------------------------
          call projection_scalaire(phi,phi0,phi1,pre,pr0,vts,vik,div,rov,rho,bip,pin,ssm,xit,mob,&
            per,frlv,pov,slv,jpenp,jpenv,bcoefv,vin,vtset,mpips,mpins,&
            typlv,energie,typntb,typnvb,milieuv,mumps_par_ps,navier,vt0,vt1,vt2,smc,smc0)
       endif
       !===============================================================================

       !-------------------------------------------------------------------------------
       !gradient de pression
       if (nb_proc>1) call mpiech(pre,kkt,mpips)
       call gradient(pre,grd0)
       if (is_poreux) where(milieuv==2) grd0=grd0*pov
       if (energie%is_changement_etat.and.energie%change_modele==change_modele_alliage) &
            grd0=grd0*frlv
       if (schema_temps==schema_temps_gear.and.navier%methode_type==navier_methode_lagrangien) &
          call gradient(pr1,grd1)
       if (navier%methode_type==navier_methode_lagrangien) then
          smc = smc + coef_temps_2/coef_temps_1*grd0 
          if (schema_temps==schema_temps_gear) smc = smc + coef_temps_3/coef_temps_1*grd1
       else
          smc = smc - grd0 
       endif
       !-------------------------------------------------------------------------------
       !temps
       smc=smc - rov*coef_temps_2*vt0/navier%pas_de_temps
       if (schema_temps==schema_temps_gear) smc=smc - rov*coef_temps_3*vt1/navier%pas_de_temps
       !-------------------------------------------------------------------------------
       !limites
       do k = 1,nbfv
          do l=1,bvp
             lvp = jpenv(l)
             smc(lvp) = smc(lvp) + bcoefv(l,k)*vin(l)
          enddo
       enddo
       !-------------------------------------------------------------------------------
       do lvp = 1,kkv
          smc(lvp) = smc(lvp) + 1.d-28*sin(dble(lvp)/20.d0)
       enddo
       !-------------------------------------------------------------------------------


       !-------------------------------------------------------------------------------
       kvi = 0
       kmi = 0
       kmv = 0
       kmm = 0
       !-------------------------------------------------------------------------------

       !-------------------------------------------------------------------------------
       ! schema hybride
       !-------------------------------------------------------------------------------
       if (navier%schema==navier_schema_hybride) then
          tet=1.d0
          do lvp=1,kkv
             do n=1,ndim
                pp(n)=1.d0
                pm(n)=1.d0

                pvg=(vts(jcoefv(lvp,kvj(2,n)))-vts(jcoefv(jcoefv(lvp,kvj(2,n)),kvj(2,n))))*&
                     (vts(lvp)-vts(jcoefv(lvp,kvj(2,n))))
                if (abs(pvg)>1.d-15) pm(n)=sign(1.d0,pvg)

                pvd=(vts(jcoefv(lvp,kvj(3,n)))-vts(lvp))*&
                     (vts(jcoefv(jcoefv(lvp,kvj(3,n)),kvj(3,n)))-vts(jcoefv(lvp,kvj(3,n))))
                if (abs(pvd)>1.d-15) pp(n)=sign(1.d0,pvd)

                if((pm(n)<0.d0).or.(pp(n)<0.d0)) tet(lvp)=0.d0
             enddo
          enddo
       endif
       !-------------------------------------------------------------------------------
		! Deewakar addition
		
		allocate (term1(kkt))
		allocate(numerator(kkt))
		allocate(denominator(kkt))
		allocate(diverge_V(kkt))
		allocate(rep_new(kkt))
		
		term1=1./xit
		numerator=((beg)**2)*tep
		denominator=(rho)*cal*((xit)**2)
		!call diverge(vts,diverge_V,scalar)   
		
		rep_new=navier%pas_de_temps*(term1+(numerator/denominator)) !*diverge_V ! remove div V
		
		! Modyfying for source term
		
		allocate(grad_T(kkv))
		allocate(flux(kkv))
		allocate(diverge_flux(kkt))
		allocate(add_smt(kkt))
		
		call gradient(tep,grad_T)
		flux=-cot(:,1,1)*grad_T
		call diverge(flux,diverge_flux,scalar)
		add_smt=(beg*diverge_flux)/(rho*cal*xit)
		
		!print *,'max val of xit an dbeg', maxval(xit),maxval(beg)
		!rep_new=rep_new-add_smt
       !-------------------------------------------------------------------------------
       ! remplissage de la matrice
       !-------------------------------------------------------------------------------
       if (ecoulement==ecoulement_compressible_non_conservatif) &
            rep=rep_new
       lax=la
       if (.not.schema_temps==schema_temps_gear.or.la>1) then
          sol = vt0
       else
          sol = 2.d0*vt0 - vt1
       endif
       call prepavol (coef,jcof,icof,sol,vil,vik,vic,vir,rep,slv,rho,rov,for, &
            xit,gra,beg,frlv,pov,mob,ssm,navier,energie,tet,bcoefv,kdv,&
            nps,jpenv,lax,utilitaires%is_onde,mpips)

!rho2=0.d0
!do l=1,kkv
!rho2=rho2+smc(l)
!enddo
!call mpisum(rho2)
!print*,'smv',rho2

!rho2=0.d0
!do lvp=1,kkv
!do i=icof(lvp),icof(lvp+1)-1
!if (coef(i)<1.D20) rho2=rho2+coef(i)
!enddo
!enddo
!call mpisum(rho2)
!print*,'coef',rho2
!!corv ; mpins
!l=0
!do lvp=1,kkv
!l=l+corv(lvp)+mpins%kic(lvp)
!enddo
!print*,'l',rang,l
       !-------------------------------------------------------------------------------


       !-------------------------------------------------------------------------------
       !periodique
       !-------------------------------------------------------------------------------
       if (.not.(navier%solveur_quantite_de_mvt%type==solveur_type_iteratif.and.&
            navier%solveur_quantite_de_mvt%iteratif==solveur_iteratif_hypre)) then
          do l=1,size(jperv,1)
             if (jperv(l)/=0) then
                lvp=jpenv(l)
                do i=icof(lvp),icof(lvp+1)-1
                   coef(i)=0.d0
                enddo
                coef(icof(lvp))  =1.d0
                coef(icof(lvp)+1)=-1.d0
                jcof(icof(lvp)+1)=jperv(l)
                smc(lvp)=0.d0
             endif
          enddo
       endif
       !-------------------------------------------------------------------------------


       !-------------------------------------------------------------------------------
       !condition limites orlanski, contraintes nulle, etc.
       if (navier%is_ouvert) call cond_lim_libre_vitesse(coef,smc,rov,vts,pre,pr1,vim,vir,rho,&
            vin,rep,icof,jcof,jpenv,typlv,nps,kkt,npt,kkv,ndim,bvp,is_debug)
       !-------------------------------------------------------------------------------


       !-------------------------------------------------------------------------------
       !impression de la matrice et de la connectique
       !-------------------------------------------------------------------------------
       impp=0
       if (impp==1) then

          if (nb_proc>1) then
             if (rang==0) then
                do ltp=0,nb_proc-1
                   if (rang==ltp) then
                      write(*,*) 'rangmat',rang
                      fmt1='(a,10i12)'
                      write(fmt1(4:5),'(i2.2)') kdv
                      write(*,fmt1) 'coefv        ',(i,i=1,kdv)
                      do lvp=1,npt
                         !!coef
                         fmt1='(1x,i4,e20.10,09e22.13)'
                         !fmt1='(1x,i4,e20.10,09e18.9)'
                         if (mpins%kic(lvp)/=0) then
                            write(fmt1(15:16),'(i2.2)') icof(lvp+1)-icof(lvp)
                            !!write(*,fmt1) lvp,smc(lvp),(coef(lts),lts=icof(lvp),icof(lvp+1)-1)
                            write(*,fmt1) corv(lvp),smc(lvp),(coef(lts),lts=icof(lvp),icof(lvp+1)-1)
                         endif
                      enddo
                   endif
                enddo
             endif

          else

             fmt1='(a,10i12)'
             write(fmt1(4:5),'(i2.2)') kdv
             write(*,fmt1) 'coefv        ',(i,i=1,kdv)
             do lvp=1,npt
                !!coef
                fmt1='(1x,i5,e12.3,14e12.3)'
                !fmt1='(1x,i4,e20.10,14e20.10)'
                !!if (mpins%kicmu(lvp)/=0) then
                !write(fmt1(15:16),'(i2.2)') icof(lvp+1)-icof(lvp)
                write(fmt1(14:15),'(i2.2)') icof(lvp+1)-icof(lvp)
                write(*,fmt1) lvp,smc(lvp),(coef(lts),lts=icof(lvp),icof(lvp+1)-1)
                !!endif
             enddo
          endif
       endif
       !-------------------------------------------------------------------------------


    
       !-------------------------------------------------------------------------------
       !modifcation matrice pour domaines fictifs
       !-------------------------------------------------------------------------------
       if (psm%penv) then
          if(navier%methode_type==navier_methode_correc_pression.or.&
           navier%methode_type==navier_methode_correc_pression_rot)then
          call pen_nav_sca(npt,nps,kkv,kic,nic,icof,jcof,smc,vts,coef,jcoefv,ndim)
        else
          call pen_nav(npt,nps,kkv,kic,nic,icof,jcof,smc,vts,coef,jcoefv,ndim)
        endif
       endif
       !-------------------------------------------------------------------------------



       !-------------------------------------------------------------------------------
       if (navier%solveur_quantite_de_mvt%type==solveur_type_direct) then
          !-------------------------------------------------------------------------------

          if (navier%solveur_quantite_de_mvt%direct==solveur_direct_mumps) then

             !mumps
             if (nb_proc==1) then
                call resol_mumps(mumps_par_ns,coef,jcof,icof,sol,smc,res,npt,nps,kic,nic,&
                     navier%solveur_quantite_de_mvt%mumps,is_debug)
             else
                call resol_mumps_dist(mumps_par_ns,coef,jcof,icof,sol,smc,res,npt,nps,corv,kkvg,mpins,&
                     navier%solveur_quantite_de_mvt%mumps,ndim,is_debug)
             endif
             kt=0

          else if (navier%solveur_quantite_de_mvt%direct==solveur_direct_pardiso) then

             !pardiso
             if (nb_proc==1) then
                pout = 0
                proc = navier%solveur_quantite_de_mvt%pardiso%threads
                call resol_pardiso(coef,jcof,icof,sol,smc,res,npt,nps,kic,nic,proc,pout)
                kt=0
             else
                print*,'pardiso fonctionne uniquement dans la version sequentielle'
                call sortie()
             endif

          else
             write(6,*) 'cette bibliotheque n existe pas pour navier'
             write(4,*) 'cette bibliotheque n existe pas pour navier'
             call sortie
          endif

          !-------------------------------------------------------------------------------
       elseif (navier%solveur_quantite_de_mvt%type==solveur_type_iteratif) then
          !-------------------------------------------------------------------------------

          if (navier%solveur_quantite_de_mvt%iteratif==solveur_iteratif_thetis) then

             if (utilitaires%is_test_perf) temps_hyprevec_debut=mpi_wtime()

             call resol_thetis(coef,jcof,icof,sol,smc,navier%pas_de_temps,&
                  navier%solveur_quantite_de_mvt%thetis,&
                  res,kt,npt,nps,ndd,kic,nic,comut,mpins,nb_proc,rang,ndim,is_debug)

             if (utilitaires%is_test_perf.and.iteration_temps>1) then
                temps_hyprevec_fin=mpi_wtime()
                temps_hyprevec_fin=temps_hyprevec_fin - temps_hyprevec_debut
                call mpimax(temps_hyprevec_fin)
                temps_hyprevec_final=temps_hyprevec_final + temps_hyprevec_fin
             endif


          else if (navier%solveur_quantite_de_mvt%iteratif==solveur_iteratif_hypre) then

             if (nb_proc>1) then
                call resol_hypre_v (coef,jcof,icof,sol,smc,res,npt,npu,npw,npuh,npvh,npwh,nps,&
                     kt,navier%solveur_quantite_de_mvt,ilower,iupper,iloweru,iupperu,ilowerv,iupperv,&
                     ilowerw,iupperw,mpins,comut,rang,nb_proc,ndim,kdv,is_debug)
             else
                print*,'hypre fonctionne uniquement dans la version parallele'
                call sortie
             endif

          else
             write(6,*) 'cette bibliotheque n existe pas pour navier'
             write(4,*) 'cette bibliotheque n existe pas pour navier'
             call sortie
          endif

          !-------------------------------------------------------------------------------
       endif
       !-------------------------------------------------------------------------------


       !-------------------------------------------------------------------------------
       vts = sol
       !-------------------------------------------------------------------------------

       !===============================================================================
       !fin de la resolution
       !===============================================================================

       !-------------------------------------------------------------------------------
       call diverge(vts,div,navier%quantite_de_mvt_div)
       !-------------------------------------------------------------------------------
		
       
	   !-------------------------------------------------------------------------------
       ! calcul de la pression
	   
	   !------------------------------------------------------------------------------
	   ! Deewakar change pressure
	   
	   rep=rep_new ! includes the compressibility as well as thermal effect as per the desired setting of rep_new .
					!This will be passed below as rep.
       !-------------------------------------------------------------------------------
       if (navier%methode_type==navier_methode_lagrangien) then
          !-------------------------------------------------------------------------------
          if (navier%is_viscosite_compression) then
             pre = pre - vil * div
             if (utilitaires%is_onde) pre = pre + vil * ssm
          else
             pre = pre - rep * div
             if (utilitaires%is_onde) pre = pre + rep * ssm
          endif
          if (navier%is_penalisation_pression) then
             do l=1,bpp
                ltp=jpenp(l)
                pre(ltp)=(pre(ltp) + bip(l)*pin(l))/(1.d0 + bip(l))
             enddo
          endif
          if (navier%is_extrapolp) call extrapol_p_coins(pre,jpenv,typlv)
          !-------------------------------------------------------------------------------
       endif
       !-------------------------------------------------------------------------------

       !-------------------------------------------------------------------------------
       navier%lagrangien_iterations_out=la
       navier%solveur_quantite_de_mvt_out(la)%residu=res
       navier%solveur_quantite_de_mvt_out(la)%iterations=kt
       navier%lagrangien_div(la)=navier%quantite_de_mvt_div
       if (navier%quantite_de_mvt_div < navier%lagrangien_tolerance) exit

       !-------------------------------------------------------------------------------
    enddo
    !-------------------------------------------------------------------------------
    !=============================================================================== 
    ! fin de la boucle du lagrangien augmente ou de l etape de prediction
    !===============================================================================


    !-------------------------------------------------------------------------------
    if (nb_proc>1) call mpiech(pre,kkt,mpips)
    !-------------------------------------------------------------------------------



    !===============================================================================
    !projection scalaire
    !===============================================================================
    !-------------------------------------------------------------------------------
    if (navier%is_correction_pression.and.abs(navier%quantite_de_mvt_div)>1.d-14) then
       !-------------------------------------------------------------------------------
       call projection_scalaire(phi,phi0,phi1,pre,pr0,vts,vik,div,rov,rho,bip,pin,ssm,xit,mob,&
            per,frlv,pov,slv,jpenp,jpenv,bcoefv,vin,vtset,mpips,mpins,&
            typlv,energie,typntb,typnvb,milieuv,mumps_par_ps,navier,vt0,vt1,vt2,smc,smc0)
       !-------------------------------------------------------------------------------
    elseif (navier%is_correction_pression) then
       !-------------------------------------------------------------------------------
       navier%solveur_projection_sca_out%residu=0.d0
       navier%solveur_projection_sca_out%iterations=0
    endif
    !-------------------------------------------------------------------------------


    !===============================================================================
    !projection vectorielle
    !===============================================================================
    !------------------------------------------------------------------------------- 
    if (navier%is_projection_vectorielle.and.abs(navier%quantite_de_mvt_div)>1.d-14) then
       !-------------------------------------------------------------------------------
       !operateur grad ( div u' ) = - grad ( div v* )
       !-------------------------------------------------------------------------------
       call projection_vectorielle (vts,rov,div,bcoefv,jpenv,jperv,jlimt,jlimt2,&
            ilimt,ilimt2,typnt,typlv,indvl,vin,navier,coef_temps_1,coef_temps_2,&
            coef_temps_3,mpins,nbloc) 
       !------------------------------------------------------------------------------ 
    elseif (navier%is_projection_vectorielle) then
       navier%solveur_projection_vec_out%residu=0.d0
       navier%solveur_projection_vec_out%iterations=0
    endif
    !-----------------------------------------------------------------------------
    !===============================================================================



    !-------------------------------------------------------------------------------
    !modifcation solution pour domaines fictifs
    !-------------------------------------------------------------------------------
    if (psm%penv) then
       if (.not.allocated(psm%vts)) allocate(psm%vts(kkv))
       psm%vts=vts
       do no=1,nbo    
          if (btest(obj(no)%pid,2).and.(size(obj(no)%heav,1).gt.1)) then
             if (obj(no)%id.eq.0) then
                do lvp=1,kkv
                   if (obj(no)%heav(lvp).ge.0.5d0) psm%vts(lvp)=0.d0
                enddo
                do ltp=1,kkt
                   if (obj(no)%hea(ltp).ge.0.01d0) div(ltp)=0.d0
                enddo
             else
                do lvp=1,kkv
                   if (obj(no)%heav(lvp).le.0.5d0) psm%vts(lvp)=0.d0
                enddo
                do ltp=1,kkt
                   if (obj(no)%hea(ltp).le.0.99d0) div(ltp)=0.d0
                enddo
             endif
          endif
       enddo
    endif
    !-------------------------------------------------------------------------------


    !===============================================================================
    !desallocation
    !===============================================================================
    !-------------------------------------------------------------------------------
    deallocate(kic)
    deallocate(sol)
    deallocate(smc)
    deallocate(tet)
    deallocate(jcof)
    deallocate(icof)
    deallocate(coef)
    deallocate(rep)

    ! Deewakar deallocation
    
    
	deallocate (term1)
	deallocate(numerator)
	deallocate(denominator)
	deallocate(diverge_V)
	deallocate(rep_new)

	deallocate(grad_T)
	deallocate(flux)
	deallocate(diverge_flux)
	deallocate(add_smt)

    if (is_turbulence .and. .not.turbulence_modele==turbulence_modele_les) &
         deallocate(grc)
    if (navier%is_tension_sup) &
         deallocate(skn)
    if (schema_temps==schema_temps_gear.and.navier%methode_type==navier_methode_lagrangien) &
         deallocate(grd1)
    if (navier%is_brinkman) deallocate(mob)
    if (energie%is_changement_etat.and.energie%change_modele==change_modele_alliage) &
         deallocate(frlv)
    !-------------------------------------------------------------------------------


    !===============================================================================
    !boussinesq
    !===============================================================================
    if (navier%is_boussinesq) then
       rov=rovb
       rho=rhob
       deallocate(rovb,rhob)
    endif
    !===============================================================================


    !-------------------------------------------------------------------------------
    if (is_debug) write (*,*) 'sortie resolution_navier'
    !-------------------------------------------------------------------------------
  end subroutine resolution_navier
  !*******************************************************************************
  !===============================================================================
end module bib_module_resolution_navier
!===============================================================================


!*******************************************************************************
subroutine extrapol_p_coins(pre,jpenv,typlv)
  !*******************************************************************************
  use module_grille, only : kkt,ndim,kkv,bvp,bvpl
  use module_coordonnees
  use module_metriques_v
  use module_connexions_t
  use module_connexions_lim
  implicit none
  !-------------------------------------------------------------------------------
  double precision,  dimension(kkt)    :: pre
  integer,           dimension(*)      :: jpenv
  integer,           dimension(bvpl,*) :: typlv
  !-------------------------------------------------------------------------------
  integer                 :: d,ln,i,k,m,n,lvp,ltp,l1,l2,ideb,ifin,l
  integer, dimension(kkv) :: jgrvp
  double precision                  :: p,a,b
  !-------------------------------------------------------------------------------

  jgrvp=-1
  do ln=1,bvp
     jgrvp(jpenv(ln))=ln
  enddo

  if (ndim==2) then
     ideb=ilimt(6,1)
     ifin=ilimt(10,1)
  else
     ideb=ilimt(20,1)
     ifin=ilimt(28,1)
  endif

  do i=ideb,ifin-1

     ltp=jlimt(i)
     k=0
     do n=1,ndim
        do m=6,7
           lvp=jcoeft(ltp,ktj(m,n))
           l=jgrvp(lvp)
           if (typlv(l,2)==2.or.typlv(l,2)==3) k=k+1
        enddo
     enddo
     if (k>=2*(ndim-1)+1) then
        p=0.d0
        do n=1,ndim
           d=0
           if (jcoeft(jcoeft(ltp,ktj(3,n)),ktj(2,n))==ltp) then
              d=3
           elseif (jcoeft(jcoeft(ltp,ktj(2,n)),ktj(3,n))==ltp) then
              d=2
           endif
           if (d==2.or.d==3) then
              l1=jcoeft(ltp,ktj(d,n))
              l2=jcoeft(l1,ktj(d,n))
              a=dev(jcoeft(ltp,ktj(4+d,n)),n)
              b=dev(jcoeft(l1,ktj(4+d,n)),n)
              p = p + (a + b)/b*pre(l1) - a*pre(l2) 
           endif
        enddo
        pre(ltp)=p/dble(ndim)
     endif

  enddo

  !-------------------------------------------------------------------------------
end subroutine extrapol_p_coins
!*******************************************************************************


!*******************************************************************************
!fin
!*******************************************************************************
