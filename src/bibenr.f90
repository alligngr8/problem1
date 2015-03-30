

!###############################################################################

!                                   energie

!###############################################################################


!===============================================================================
module bib_module_resolution_energie
contains
  !===============================================================================
  !*******************************************************************************
  subroutine resolution_energie(tep,rho,vts,smc,smd,smcr,tp0,tp1,tin,cal,ent,ent0,ent1, &
       vim,pep,rcp,rov,com,cot,frs,cou,esp,compliq,bcoeft,jpent,typlt,&
       jpert,vt0,mte,kte,wee,foe,energie,navier,advection,&
       mpien,mumps_par_en,beg,xit)
    !*******************************************************************************
    use bib_module_cond_lim_libre_scalaire
    use bib_module_conextra
    use bib_module_pen_enr
    use bib_module_pen_mat_trans
    use bib_module_prepatoc
    use bib_module_resistance_thermique_contact
    use bib_module_resol_mumps_dist
    use bib_module_resol_hypre_t
    use bib_module_resol_thetis
    use bib_module_tracer_mpi2
    !-------------------------------------------------------------------------------
    use module_commun
    use module_connexions
    use module_coordonnees
    use module_generique_mpimax
    use module_grille
    use module_impressions
    use module_metriques
    use module_mpi
    use module_mpif
    use module_objet, only : obj,nbo
    use module_thermophysique
    !-------------------------------------------------------------------------------
    use module_struct_advection
    use module_struct_energie
    use module_struct_mumps
    use module_struct_navier
    use module_struct_psm
    !-------------------------------------------------------------------------------
    implicit none
    !-------------------------------------------------------------------------------
    !variables globales
    !-------------------------------------------------------------------------------
    integer, dimension(:), allocatable :: kte
    integer, dimension(:), allocatable :: jpent,typlt,jpert
    !-------------------------------------------------------------------------------
    double precision, dimension(:),    allocatable  :: tep,tp0,tin,cal,rho,pep,rcp,vim,ent, &
         ent0,ent1,vts,rov,com,vt0,compliq,xit,beg
    double precision, dimension(:,:),  allocatable  :: mte
    double precision, dimension(:,:),  allocatable  :: bcoeft
    double precision, dimension(:),    allocatable  :: wee,foe
    double precision, dimension(:,:,:),allocatable  :: cot
    !-------------------------------------------------------------------------------
    double precision, dimension(:),     allocatable  :: smcr,tp1,smc,smd
    double precision, dimension(:,:),   allocatable  :: cou,frs,esp
    !-------------------------------------------------------------------------------
    type(struct_energie)           :: energie
    type(struct_navier)            :: navier
    type(struct_advection_systeme) :: advection
    type(struct_mpieq)             :: mpien
    type(dmumps_struc)             :: mumps_par_en
    !-------------------------------------------------------------------------------
    !variables locales
    !-------------------------------------------------------------------------------
    integer, dimension(:),  allocatable  :: jcof,icof
    integer, dimension(:),  allocatable  :: kic
    double precision,  dimension(:),  allocatable  :: coef,coef_diag
    !-------------------------------------------------------------------------------
    double precision, dimension(:), allocatable :: vts2,pro,gts,bio,sou,ret,frl,frl0,&
         frl1,frl2,bnul,tnul
    double precision, dimension(:), allocatable :: tpp0,smt,nul,rhob,rovb,sol
    double precision, dimension(:,:), allocatable :: vtsp
    !-------------------------------------------------------------------------------
    integer :: kdc,lts,npt,ltp,ndd,nps,lcp,n,bit,solveur_iterations_moy
    integer :: i,j,l,impp,nic,pout,proc
    integer :: kdt = 0
    !-------------------------------------------------------------------------------
    integer :: nbb,ntt,nbcfl
    integer :: nes,npes,npars
    integer, dimension(:), allocatable :: ktes
    double precision,  dimension(:,:),  allocatable  :: mtes
    double precision,  dimension(:),    allocatable  :: wees,foes
    double precision :: vcar,solveur_residu_moy,dtcfl,dtrest,tstart,tend,tolerance_save
    logical :: coupler,is_reduced_tolerance
    !-------------------------------------------------------------------------------

	! Deewakar variables additions--------------------------------------------------
	double precision,  dimension(:), allocatable   :: diverge_V,new_term
	double precision 			       :: scalar
    !-------------------------------------------------------------------------------
    if (is_debug) write (*,*) 'entree energie'
    !-------------------------------------------------------------------------------
    coupler = .false.
    allocate(smt(kkt))
    allocate(sol(kkt))
    if (energie%is_changement_etat) allocate(coef_diag(kkt))

    !-------------------------------------------------------------------------------
    if (energie%is_enthalpie) then
       sol=ent-energie%enthalpie_ref
    else
       sol=tep-energie%temperature_ref
    endif
    !-------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    !modifcation vitesses pour domaines fictifs
    !-------------------------------------------------------------------------------
    if (psm%penv) then
       allocate(vts2(kkv)); vts=0.d0
       vts2=vts
       vts=psm%vts
    endif
    !-------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    !initialisation changement de phase
    !-------------------------------------------------------------------------------
    if (energie%is_changement_etat) then
       solveur_iterations_moy=0
       solveur_residu_moy=0.d0
       energie%change_iterations_out=0
       energie%change_residu=1.d20
       allocate(tpp0(kkt))
       allocate(nul(kkt)) ; nul=0.d0
       allocate(frl(kkt))
       allocate(frl0(kkt))
       allocate(frl1(kkt))
       allocate(frl2(kkt))

       if (advection%is_grille_duale) THEN
          do ltp=1,kkt
             lcp=jcoeft(ltp,ktc)              
             if (frs(ltp,1)>cou(lcp,1)) frs(ltp,1)=cou(lcp,1)
          enddo
       else
          do ltp=1,kkt                         
             if (frs(ltp,1)>cou(ltp,1)) frs(ltp,1)=cou(ltp,1)
          enddo
       endif

       frl(:)=1.d0 - frs(:,1)
       frl0=frl
       frl2=frl

       allocate(sou(kkt))
       sou=0.d0
       tpp0 = tep
       is_reduced_tolerance=.false.
       if (energie%change_modele==change_modele_franc_diffus) then
          is_reduced_tolerance=.true.
          if (energie%solveur%iteratif==solveur_iteratif_hypre) tolerance_save=energie%solveur%hypre%tolerance
          if (energie%solveur%iteratif==solveur_iteratif_thetis) tolerance_save=energie%solveur%thetis%tolerance
       endif

    endif
    if (.not.energie%is_changement_etat) then
       energie%change_residu = 1.d10
       energie%change_tolerance = -1.d0
       energie%change_iterations=0
       energie%change_iterations_out=0
    endif
    !-------------------------------------------------------------------------------


    !-------------------------------------------------------------------------------
    ! calcul des connectivites et de la matrice
    !-------------------------------------------------------------------------------
    call conextra (coef,jcof,icof,jcoeft,ktj,energie%schema,energie_schema_centre,&
         energie_schema_hybride,energie_schema_upwind,energie_schema_quick,energie_schema_double,&
         energie_schema_tvd,energie_schema_weno,energie_schema_vof_sm,energie%is_ani,kkt,kdt,kdc,&
         npt,nps,ndim,ndd)
    !-------------------------------------------------------------------------------
    ! calcul de kic
    !-------------------------------------------------------------------------------
    ! npt est l ordre du systeme lineaire
    ! nic est le nombre d'inconnues a resoudre
    ! kic est un tableau qui contient les inconnues a resoudre effectivement
    !-------------------------------------------------------------------------------
    allocate (kic(npt))
    !-------------------------------------------------------------------------------


    !-------------------------------------------------------------------------------
    !boussinesq
    !-------------------------------------------------------------------------------
    if (navier%is_boussinesq) then ! attention pas fait en diphasique
       allocate(rovb(kkv)) ; allocate(rhob(kkt))
       rovb=rov
       rhob=rho
       rov=fluide%masse_vol
       rho=fluide%masse_vol
    endif
    !-------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    if (energie%is_dissipation.and.is_navier) allocate(pro(kkt))
    if (energie%schema==energie_schema_tvd.or.&
         energie%schema==energie_schema_weno.or.&
         energie%schema==energie_schema_vof_sm) then
       allocate(bio(btp)) 
       if (energie%schema==energie_schema_vof_sm) then
         allocate(bnul(btp)) ; bnul=0.d0
         allocate(tnul(btp)) ; tnul=0.d0
       endif
       allocate(ret(kkt))
       allocate(gts(kkv))
    endif
    !-------------------------------------------------------------------------------


    !-------------------------------------------------------------------------------
    if (energie%schema==energie_schema_tvd.or.&
         energie%schema==energie_schema_weno.or.&
         energie%schema==energie_schema_vof_sm) then
       !------------------------------------------------------------------------------
       if (.not.is_poreux) then
          bio(:) = bcoeft(:,1)
          ret=1.d0
          call vitgrad (vts,gts)
          if (energie%schema==energie_schema_tvd) then
             call voftvd (tep,gts,vts,bio,tin,ret,jcoeft,&
                  ktj,energie%pas_de_temps,kkt,kkv,ndim,det,&
                  energie%iterations_advection_out,is_debug,.false.,jpent,btp,mpien,.false.)
          elseif(energie%schema==energie_schema_weno) then
             call vofweno (tep,gts,vts,bio,tin,ret,jcoeft,&
                  ktj,energie%pas_de_temps,kkt,kkv,ndim,det,   &
                  energie%iterations_advection_out,is_debug,.false.,jpent,btp,mpien,.false.)
          elseif(energie%schema==energie_schema_vof_sm) then
             call tracer_mpi2(tep,vts,vt0,mte,kte,wee,foe,dev,det,bio,tin,energie%pas_de_temps, &
                      jcoeft,jcoefv,jpent,kro,ktj,kvj,energie%iterations_advection_out,2,&
                      kkt,kkv,iteration_temps,ndim,energie%npe,energie%ne,energie%npar,btp,&
                      is_periodique_x1,is_periodique_x2,is_periodique_x3,mpien,&
                      .false.,nbb,ntt,nbcfl,dtcfl,dtrest,energie%enrichissement)

          endif
          if (schema_temps==schema_temps_gear) tp1 = tp0
          tp0 = tep
       else
          print*,'schema tvd/weno/vsm pas fait pour l energie si poreux'
          call sortie
       endif
       !-------------------------------------------------------------------------------
    endif
    !-------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    if (energie%is_changement_etat.and.energie%change_modele==change_modele_franc_diffus.and.&
         energie%change_delta_t_tolerance_auto) then
       if (nb_proc>1) then
          tstart=mpi_wtime()
       else 
          call cpu_time(tstart)
       endif
    endif
    !-------------------------------------------------------------------------------


    !-------------------------------------------------------------------------------
    !boucle principale du point fixe, active si changement d'etat
    !-------------------------------------------------------------------------------
    do while (energie%change_residu >= energie%change_tolerance.and.&
         energie%change_iterations_out <= energie%change_iterations)
       !-------------------------------------------------------------------------------

       !-------------------------------------------------------------------------------
       if (energie%is_changement_etat) then
          !-------------------------------------------------------------------------------
          energie%change_iterations_out=energie%change_iterations_out + 1
          energie%change_residu = 0.d0
          frl1=frl
          tpp0=tep !temperature at previous iteration
          !-------------------------------------------------------------------------------
       endif
       !-------------------------------------------------------------------------------


       !-------------------------------------------------------------------------------
       !second membre
       !-------------------------------------------------------------------------------
       !temps
       !-------------------------------------------------------------------------------
       if(energie%is_enthalpie) then
          !-------------------------------------------------------------------------------
          if (.not.is_poreux) then
             smt=- coef_temps_2*(ent0 - energie%enthalpie_ref)/energie%pas_de_temps
             if (schema_temps==schema_temps_gear)&
                  smt=smt - coef_temps_3*(ent1 - energie%enthalpie_ref)/energie%pas_de_temps
          else
             print*,'equation de l energie en enthalpie et en poreux pas faite'
             call sortie()
          endif
          !-------------------------------------------------------------------------------
       else
          !-------------------------------------------------------------------------------
          if (.not.is_poreux) then
             smt=- coef_temps_2*(tp0 - energie%temperature_ref)*rho*cal/energie%pas_de_temps
             if (schema_temps==schema_temps_gear)&
                  smt=smt - coef_temps_3*(tp1 - energie%temperature_ref)*rho*cal/energie%pas_de_temps
          else
             smt=- coef_temps_2*(tp0-energie%temperature_ref)*rcp/energie%pas_de_temps
             if (schema_temps==schema_temps_gear)&
                  smt=smt - coef_temps_3*(tp1-energie%temperature_ref)*rcp/energie%pas_de_temps
          endif
          !-------------------------------------------------------------------------------
       endif
       !limite
       if (.not.energie%is_enthalpie) then
          do l=1,btp
             ltp=jpent(l)
             smt(ltp)=smt(ltp) + bcoeft(l,1)*(tin(l) - energie%temperature_ref)
          enddo
       elseif (energie%is_enthalpie.and.energie%is_limite_enthalpie) then
          do l=1,btp
             ltp=jpent(l)
             smt(ltp)=smt(ltp) + bcoeft(l,1)*(tin(l) - energie%enthalpie_ref)
          enddo
       elseif (energie%is_enthalpie.and.energie%is_limite_temperature) then
          do l=1,btp
             ltp=jpent(l)
             smt(ltp)=smt(ltp) + bcoeft(l,1)*(tin(l)*rho(ltp)*cal(ltp) - energie%enthalpie_ref)
          enddo
       endif
       !sources
       if (energie%is_terme_source) smt=smt + smc
       if (energie%is_rayonnement) smt=smt + smcr
       if (energie%is_changement_etat) smt=smt - sou
       !-------------------------------------------------------------------------------

       !-------------------------------------------------------------------------------
       !terme de dissipation visqueuse 
       !-------------------------------------------------------------------------------
       if (energie%is_dissipation) then
          !-------------------------------------------------------------------------------
          if (is_navier) then
             call termpro (pro,vts)
             smt = smt + vim*pro
          elseif (is_darcy) then
             allocate(vtsp(kkt,ndim))
             call intervitesse(vts,vtsp)
             do ltp=1,kkt
                vcar=0.d0
                do n=1,ndim
                   vcar=vcar + vtsp(ltp,n)**2
                enddo
                smt = smt + vim/pep*vcar
             enddo
             deallocate(vtsp)
          endif
       endif
       !-------------------------------------------------------------------------------
	   !-------------------------------------------------------------------------------
       if (energie%is_resistance_th) call resistance_thermique_contact(cot,cou)
       !-------------------------------------------------------------------------------
! mod
       !-------------------------------------------------------------------------------
       

	 ecoulement=ecoulement_incompressible
	 ! Deewakar Energy modifcation
		allocate(diverge_V(kkt))
		allocate(new_term(kkt))
		  
		call diverge(vts,diverge_V,scalar)
		new_term=(beg*tep)*diverge_V*(1.d0/xit)
		smt=smt - new_term
		  
	  ! end if 
		
	   !remplissage de la matrice
       !-------------------------------------------------------------------------------
       call prepatoc (coef,jcof,icof,tep,vts,rho,cal,com,cot,rov,smd,det,dev,dem,     &
            bcoeft,jpent,jcoeft,jcoefv,ktj,kvj,kro,              &
            energie%schema,energie_schema_centre,energie_schema_hybride,&
            energie_schema_upwind,energie_schema_quick,energie_schema_double,&
            energie%is_diffusion,energie%is_advection,                 &
            coupler,is_compressible,energie%is_iso,energie%is_ort,energie%is_ani,&
            energie%pas_de_temps,coef_temps_1,kkt,kkv,kkm,kdt,kdc,nps,ndim,nbft,ncd,&
            nce,energie%is_terme_lineaire,btp,energie%is_enthalpie,psm%pent)
			
		
		ecoulement=ecoulement_compressible_non_conservatif
       !-------------------------------------------------------------------------------


       !-------------------------------------------------------------------------------
       !poreux
       !-------------------------------------------------------------------------------
       if (is_poreux) then
          do ltp=1,kkt
             l=icof(ltp)
             coef(l)=coef(l) - coef_temps_1*rho(ltp)*cal(ltp)/energie%pas_de_temps + &
                  coef_temps_1*rcp(ltp)/energie%pas_de_temps
          enddo
       endif
       !-------------------------------------------------------------------------------

       !-------------------------------------------------------------------------------
       !periodique
       !-------------------------------------------------------------------------------
       if (.not.(energie%solveur%type==solveur_type_iteratif.and.&
            energie%solveur%iteratif==solveur_iteratif_hypre)) then
          do l=1,size(jpert,1)
             if (jpert(l)/=0) then
                ltp=jpent(l)
                do i=icof(ltp),icof(ltp+1)-1
                   coef(i)=0.d0
                enddo
                coef(icof(ltp))  = 1.d0
                coef(icof(ltp)+1)=-1.d0
                jcof(icof(ltp)+1)=jpert(l)
                !   print*,'period',rang,ltp,jpert(l),cort(jpert(l))
                smt(ltp)=0.d0
             endif
          enddo
       endif
       !-------------------------------------------------------------------------------

       !-------------------------------------------------------------------------------
       !cl orlanski
       !-------------------------------------------------------------------------------
       if (energie%is_ouvert) &
            call cond_lim_libre_scalaire(coef,smt,tep - energie%temperature_ref,vts,&
            tin,icof,jcof,typlt,jpent,energie%pas_de_temps,nps,kkt,npt,kkv,btp,btpl,ndim,is_debug)
       !-------------------------------------------------------------------------------


       !-------------------------------------------------------------------------------
       if (energie%is_changement_etat) then
          do ltp=1,kkt
             coef_diag(ltp)=coef(icof(ltp))
          enddo
          !-------------------------------------------------------------------------------
          if (energie%change_modele==change_modele_franc_diffus) then
             !-------------------------------------------------------------------------------
             if (is_reduced_tolerance) then
                if (energie%solveur%iteratif==solveur_iteratif_hypre) &
                     energie%solveur%hypre%tolerance=energie%change_delta_t_tolerance
                if (energie%solveur%iteratif==solveur_iteratif_thetis) &
                     energie%solveur%thetis%tolerance=energie%change_delta_t_tolerance
             endif
             sol=0.d0  !initial guess=0

             do j=1,btp
                ltp=jpent(j)
                l=icof(ltp)                                                                            
                lcp=ltp                 
                if (advection%is_grille_duale) lcp=jcoeft(ltp,ktc)    
                if (coef(l) > 1.d35) then                    
                   tep(ltp)=smt(ltp)/coef(l)+energie%temperature_ref
                   smt(ltp)=0.d0                 
                endif
             enddo

             if (nb_proc>1) call mpiech(tep,kkt,mpien)

             do ltp=1,kkt  
                l=icof(ltp)                                                                            
                lcp=ltp                 
                if (advection%is_grille_duale) lcp=jcoeft(ltp,ktc)    
                if (coef_diag(ltp) < 1.d35) then                    
                   smt(ltp)=- smt(ltp) + phases(1)%masse_vol*energie%change_chaleur_latente*(frl(ltp) - frl0(ltp))/energie%pas_de_temps
                   do i=0,kdt-1 
                      smt(ltp)=smt(ltp) + coef(l+i)*(tep(jcof(l+i)) - energie%temperature_ref)
                   enddo
                   coef(l)=coef(l) + phases(1)%masse_vol*energie%change_chaleur_latente*5.d-1/energie%change_delta_t*&
                        (1.d0-dtanh((energie%change_temp_fusion - tep(ltp))/energie%change_delta_t)**2)/energie%pas_de_temps*cou(lcp,1)                                
                endif
             enddo

          endif
       endif
       !-------------------------------------------------------------------------------


       !-------------------------------------------------------------------------------
       !impression matrice
       !-------------------------------------------------------------------------------
       impp=1
       if(impp.eq.0) then
          !do l=0,nb_proc-1
             if (rang==0) then
                write(*,*) 'rang',rang,'coeft',kdt,kdc,kkt
                write(*,601) (ltp,ltp=1,kdt+kdc)
                do ltp=1,kkt

                   ! if (mpienergie%kic(ltp)>0) write(*,600) ltp,smt(ltp), &
                   write(*,600) ltp,smt(ltp), &
                        (coef(lts),lts=1+(ltp-1)*(kdt+kdc), &
                                !                                      1+(ltp-1)*(kdt+kdc)+4)
                        (ltp-1)*(kdt+kdc)+(kdt+kdc))
                enddo
             endif
          !enddo
600       format(1x,i5,e24.15,' ',9e24.15,i5)
601       format(1x,16x,9i11)
          !-------------------------------------------------------------------------------
       endif
       !-------------------------------------------------------------------------------

       !-------------------------------------------------------------------------------
       ! l utilisateur est cense remplir ici kic avec les valeurs des inconnues
       ! ou laisser kic a zero si l inconnue n est pas a resoudre
       !-------------------------------------------------------------------------------
       nic = npt
       kic = 0
       do ltp=1,npt
          kic(ltp)= ltp
       enddo
       !-------------------------------------------------------------------------------

       !-------------------------------------------------------------------------------
       !modifcation matrice pour domaines fictifs
       !-------------------------------------------------------------------------------
       if (psm%pent) then
          bit=0
          do n=1,nbo
             if (btest(obj(n)%pid,0)) bit=1
             if (btest(obj(n)%pid,4)) bit=4    
          enddo
          if (bit.eq.1) then
             call pen_enr(npt,nps,kkt,kic,nic,icof,jcof,smt,sol,coef,jcoeft,ndim)
          else if (bit.eq.4) then
             call pen_mat_trans(npt,nps,kkt,kic,nic,icof,jcof,smt,sol,com,coef,det,psm%altert,&
                  psm%idpt,psm%mtg,jcoeft,psm%btp,psm%ncoft,psm%hasht,ndim)
          endif
       endif
       !-------------------------------------------------------------------------------

       !-------------------------------------------------------------------------------
       !resolution
       !-------------------------------------------------------------------------------

       !-------------------------------------------------------------------------------
       if (energie%solveur%type==solveur_type_iteratif) then
          !-------------------------------------------------------------------------------

          if (energie%solveur%iteratif==solveur_iteratif_thetis) then
             !thetis
             call resol_thetis(coef,jcof,icof,sol,smt,energie%pas_de_temps,&
                  energie%solveur%thetis,&
                  energie%solveur_out%residu,energie%solveur_out%iterations,&
                  npt,nps,ndd,kic,nic,comut,mpien,nb_proc,rang,ndim,is_debug)
          else if (energie%solveur%iteratif==solveur_iteratif_hypre) then
             !hypre
             if (energie%solveur%hypre%type==solveur_hypre_cg.and.is_navier) &
                  energie%solveur%hypre%type=solveur_hypre_bicg
             if (nb_proc>1) then
                call resol_hypre_t (coef,jcof,icof,sol,smt,kkt,npt,nps,kdt,energie%solveur,&
                     energie%solveur_out,nb_proc,ilower,iupper,comut,rang,ndim,mpien,is_debug)
             else
                print*,'hypre fonctionne uniquement dans la version parallele'
                call sortie()
             endif

          else
             write (6,*) 'cette bibliotheque n existe pas pour l energie'
             write (4,*) 'cette bibliotheque n existe pas pour l energie'
             call sortie()
          endif

          !-------------------------------------------------------------------------------
       elseif (energie%solveur%type==solveur_type_direct) then
          !-------------------------------------------------------------------------------

          if (energie%solveur%direct==solveur_direct_mumps) then
             !mumps
             if (nb_proc==1) then
                call resol_mumps(mumps_par_en,coef,jcof,icof,sol,smt,energie%solveur_out%residu,&
                     npt,nps,kic,nic,energie%solveur%mumps,is_debug)
             else
                call resol_mumps_dist(mumps_par_en,coef,jcof,icof,sol,smt,&
                     energie%solveur_out%residu,npt,nps,cort,kktg,mpien,&
                     energie%solveur%mumps,ndim,is_debug)
             endif
             energie%solveur_out%iterations=0

          else if (energie%solveur%direct==solveur_direct_pardiso) then
             !pardiso
             if (nb_proc==1) then
                pout = 0
                proc = energie%solveur%pardiso%threads
                call resol_pardiso(coef,jcof,icof,sol,smt,energie%solveur_out%residu,npt,nps,kic,nic,proc,pout)
                energie%solveur_out%iterations=0
             else
                print*,'pardiso fonctionne uniquement dans la version sequentielle'
                call sortie()
             endif

          else
             write (6,*) 'cette bibliotheque n existe pas pour l energie'
             write (4,*) 'cette bibliotheque n existe pas pour l energie'
             call sortie()
          endif

          !-------------------------------------------------------------------------------
       endif
       !-------------------------------------------------------------------------------

       !-------------------------------------------------------------------------------
       if (energie%is_enthalpie) then
          ent=sol(1:kkt) + energie%enthalpie_ref 
          tep=ent/(rho*cal) !approx. 
       else
          if (energie%change_modele/=change_modele_franc_diffus) then
             tep=sol + energie%temperature_ref
          else
             tep=tep - sol
          endif
       endif
       !-------------------------------------------------------------------------------
       !print*,'tep',tep

       !-------------------------------------------------------------------------------
       if (energie%is_changement_etat) then
          !-------------------------------------------------------------------------------
          if (energie%change_modele/=change_modele_franc_diffus) then
             !-------------------------------------------------------------------------------
             if (energie%change_modele==change_modele_alliage) then
                do ltp=1,kkt
                   compliq(ltp)=esp(ltp,1)/(energie%change_allia_coef_partition + &
                        (1.d0 - energie%change_allia_coef_partition)*frl(ltp))
                   if(compliq(ltp)>energie%change_allia_c_eutectique) &
                        compliq(ltp)=energie%change_allia_c_eutectique
                enddo
             endif
             do ltp = 1,kkt
                lcp=ltp
                if (advection%is_grille_duale) lcp=jcoeft(ltp,ktc)
                !!if (cou(lcp,1) >= 0.5d0) then
                if (cou(lcp,1) >= 1.d-13) then
                   if (energie%change_modele==change_modele_alliage) then
                      energie%change_temp_fusion=energie%change_allia_tm_zero + &
                           energie%change_allia_pente_liquidus*compliq(ltp)
                      if(compliq(ltp).eq.energie%change_allia_c_eutectique) &
                           energie%change_temp_fusion=energie%change_allia_t_eutectique
                   endif

                   if (coef_diag(ltp) >1.d35) then
                      sou(ltp)=0.d0
                      frl(ltp)=1.d0-cou(lcp,1) 
                      if (tep(ltp)>energie%change_temp_fusion) frl(ltp)=1.d0 

                   else    

                      if ((tep(ltp)-energie%change_temp_fusion)*(tpp0(ltp)-energie%change_temp_fusion)<0.d0 .and.&
                           abs(tep(ltp)-tpp0(ltp))>1.d-15 .and. energie%change_iterations_out>1 ) then 
                         frl(ltp)=(frl2(ltp)*(tep(ltp)-energie%change_temp_fusion)-frl1(ltp)*(tpp0(ltp)-energie%change_temp_fusion))/(tep(ltp)-tpp0(ltp))   !=SECANT METHOD (locally)
                      else     
                         frl(ltp) = frl1(ltp) +  (tep(ltp)-energie%change_temp_fusion)*coef_diag(ltp)/&
                              (phases(1)%masse_vol*energie%change_chaleur_latente)*energie%pas_de_temps     !="NEWTON" METHOD  (locally and only diagonal of Jacobian Matrix)

                         if (frl(ltp) > 1.d0) then
                            frl(ltp)=1.d0
                         elseif (frl(ltp)<1.d0-cou(lcp,1)) then
                            frl(ltp)=1.d0-cou(lcp,1)                                                   
                         endif
                      endif

                      sou(ltp) = (frl(ltp) - frl0(ltp))*phases(1)%masse_vol*energie%change_chaleur_latente/energie%pas_de_temps
                      if (abs(frl(ltp)-frl1(ltp)) > energie%change_residu) then
                         energie%change_residu = abs(frl(ltp)-frl1(ltp))                                                    
                      endif

                   endif
                endif
             enddo
             frl2=frl1
             if (nb_proc>1) then
                call mpimax(energie%change_residu)
                call mpiech(tep,kkt,mpien)
             endif
             !-------------------------------------------------------------------------------
          elseif (energie%change_modele==change_modele_franc_diffus) then
             !-------------------------------------------------------------------------------
             do ltp = 1,kkt
                lcp=ltp                
                if (advection%is_grille_duale) lcp=jcoeft(ltp,ktc)
                !!if (cou(lcp,1) >= 0.5d0) then
                if (cou(lcp,1) >= 1.d-13) then
                   frl(ltp) = 5.d-1*(1.d0-dtanh( (energie%change_temp_fusion-tep(ltp))/energie%change_delta_t))*cou(lcp,1)+(1.d0-cou(lcp,1))                    
                   if (coef_diag(ltp)<1.d35) then    

                      if (abs(frl(ltp)-frl1(ltp)) > energie%change_residu) then
                         energie%change_residu = abs(frl(ltp)-frl1(ltp))
                      endif

                      if (energie%pas_de_temps<1.d35 .and. (abs(frl(ltp)-frl1(ltp))>5.d-1 .or. &
                           (tep(ltp)-energie%change_temp_fusion)*(tpp0(ltp)-energie%change_temp_fusion)<0.d0)) then
                         !(frl(ltp)-(1.d0-cou(lcp,1))-5.d-1*cou(lcp,1))*(frl1(ltp)-(1.d0-cou(lcp,1))-5.d-1*cou(lcp,1))<-1.d-13)) then    !stabilization of the newton method / dont limit for steady case
                         tep(ltp)=energie%change_temp_fusion
                         frl(ltp)=5.d-1*cou(lcp,1) +(1.d0-cou(lcp,1))                                                                    
                      endif

                   endif
                endif
             enddo
             if (nb_proc>1) then
                call mpimax(energie%change_residu)
                call mpiech(tep,kkt,mpien)
             endif
             if (is_reduced_tolerance) then
                if (mod(energie%change_iterations_out,5)==0.and.energie%change_residu>0.49d0) then
                   energie%change_delta_t_direction=1
                   energie%change_delta_t_tolerance=max(tolerance_save,energie%change_delta_t_tolerance/energie%change_delta_t_multiplicateur)                                                                        
                endif
                if(energie%solveur%iteratif==solveur_iteratif_hypre) then 
                   if (energie%change_residu<=energie%change_tolerance.and.energie%change_iterations_out==1) then                         
                      is_reduced_tolerance=.false.
                      solveur_residu_moy=0.d0
                      energie%change_iterations_out=0 
                      energie%change_residu=1.d0
                   endif
                   energie%solveur%hypre%tolerance=tolerance_save
                else
                   if (energie%change_residu<=energie%change_tolerance.and.energie%change_iterations_out==1) then   
                      is_reduced_tolerance=.false.
                      solveur_residu_moy=0.d0
                      energie%change_iterations_out=0 
                      energie%change_residu=1.d0
                   endif
                   energie%solveur%thetis%tolerance=tolerance_save
                endif
             endif
             !-------------------------------------------------------------------------------
          endif
          !-------------------------------------------------------------------------------

          !-------------------------------------------------------------------------------

          if (energie%pas_de_temps>1.d35) energie%change_residu=0.d0 !limit the steady case to a single iteration 

          solveur_iterations_moy = solveur_iterations_moy + energie%solveur_out%iterations
          solveur_residu_moy = solveur_residu_moy + energie%solveur_out%residu
          if (rang==0.and.impressions%is_impression_solveur) &
               write (6,603) '   Changement de phase - iter. et res. du point fixe :',&
               energie%change_iterations_out,energie%change_residu

          !-------------------------------------------------------------------------------
       else !not changement etat
          !-------------------------------------------------------------------------------
          energie%change_tolerance = 1.d+40
          !-------------------------------------------------------------------------------
       endif
       !-------------------------------------------------------------------------------
603    format (a,i8,f24.15)

       !-------------------------------------------------------------------------------
    enddo
    !-------------------------------------------------------------------------------


    !-------------------------------------------------------------------------------
    if (energie%is_changement_etat.and.energie%change_modele==change_modele_franc_diffus.and.&
         energie%change_delta_t_tolerance_auto) then
       !-------------------------------------------------------------------------------        
       if (nb_proc>1) then
          tend=mpi_wtime()
       else 
          call cpu_time(tend)
       endif
       tend=tend-tstart
       if (nb_proc>1) call mpimax(tend)

       if (is_reduced_tolerance) then
          if (rang==0) print*,'Newton_eps',energie%change_delta_t_tolerance,'Time iteration',tend
          if (tend>energie%change_delta_t_temps) energie%change_delta_t_direction=-energie%change_delta_t_direction
          energie%change_delta_t_tolerance=max(min(1.d-1,energie%change_delta_t_tolerance*&
               energie%change_delta_t_multiplicateur**(energie%change_delta_t_direction)),tolerance_save)        
          energie%change_delta_t_temps=tend        
       endif
    endif
    !-------------------------------------------------------------------------------



    !-------------------------------------------------------------------------------
    !advection de frs
    !-------------------------------------------------------------------------------
    if (energie%is_changement_etat) then
       frs(:,1) = 1.d0 - frl(:)
       if (energie%schema==energie_schema_tvd) then
          call voftvd (frs(:,1),gts,vts,nul,nul,ret,jcoeft,ktj,energie%pas_de_temps,kkt,kkv,ndim,det,&
               energie%iterations_advection_out,is_debug,.false.,jpent,btp,mpien,.false.)
       elseif(energie%schema==energie_schema_weno) then
          call vofweno (frs(:,1),gts,vts,nul,nul,ret,jcoeft,ktj,energie%pas_de_temps,kkt,kkv,ndim,det,   &
               energie%iterations_advection_out,is_debug,.false.,jpent,btp,mpien,.false.)
       elseif(energie%schema==energie_schema_vof_sm) then
           allocate(mtes(size(mte,1),size(mte,2)),ktes(size(kte,1)),wees(size(wee,1)),foes(size(foe,1)))
           mtes=mte; ktes=kte; wees=wee; foes=foe
           npars=energie%npar
           nes=energie%ne
           npes=energie%npe
           call tracer_mpi2(frs(:,1),vts,vt0,mtes,ktes,wees,foes,dev,det,bnul,tnul,energie%pas_de_temps,&
                      jcoeft,jcoefv,jpent,kro,ktj,kvj,energie%iterations_advection_out,2,&
                      kkt,kkv,iteration_temps,ndim,npes,nes,npars,btp,&
                      is_periodique_x1,is_periodique_x2,is_periodique_x3,mpien,&
                      .false.,nbb,ntt,nbcfl,dtcfl,dtrest,energie%enrichissement)
           deallocate(mtes,ktes,wees,foes)
       endif
    endif
    !-------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    ! vof-sm
    !-------------------------------------------------------------------------------
    if(energie%schema==energie_schema_vof_sm) then
       !call projtrace (tep-tp0,kte,foe,wee,energie%npe,kkt)
       call projtrace (tep-tp0,tep,kte,foe,wee,energie%npe,kkt)
    endif
    !-------------------------------------------------------------------------------


    !-------------------------------------------------------------------------------
    if (energie%is_changement_etat) then
       if (energie%change_iterations_out/=0) then
          energie%solveur_out%iterations = solveur_iterations_moy/energie%change_iterations_out
          energie%solveur_out%residu = solveur_residu_moy/dble(energie%change_iterations_out) 
       endif
    endif
    !-------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    !desallocation

    !################################################################################
    deallocate(diverge_V)
    deallocate(new_term)

    !##################################################################################
    !-------------------------------------------------------------------------------
    if (psm%penv) then
       vts=vts2
       deallocate(vts2)
    endif
    deallocate(coef,jcof,icof,sol)
    if (energie%change_modele==change_modele_franc_diffus) deallocate(coef_diag)
    deallocate(smt,kic)
    if (energie%is_changement_etat) deallocate(tpp0,nul,sou,frl,frl0,frl1,frl2)
    if (navier%is_boussinesq) then
       rov=rovb
       rho=rhob
       deallocate(rhob,rovb)
    endif
    !-------------------------------------------------------------------------------
    if (energie%is_dissipation.and.is_navier) deallocate(pro)
    if (energie%schema==energie_schema_tvd.or.&
         energie%schema==energie_schema_weno.or.&
         energie%schema==energie_schema_vof_sm) then
        deallocate(ret,gts,bio)
        if (energie%schema==energie_schema_vof_sm) deallocate(bnul,tnul)
    endif
    !-------------------------------------------------------------------------------


    !-------------------------------------------------------------------------------
    if (is_debug) write (*,*) 'sortie energie'
    !-------------------------------------------------------------------------------
  end subroutine resolution_energie
  !*******************************************************************************


  !===============================================================================
end module bib_module_resolution_energie
!===============================================================================

!*******************************************************************************
!fin
!*******************************************************************************
