

!############################################################################################################
!misc_subs.F90:
      subroutine zcoor(itag,inode,kbpl,ztmp)
!-------------------------------------------------------------------------------
!     Calculate z-coord. at a _wet_ node

!############################################################################################################
!misc_subs.F90:
      subroutine levels0(iths,it)
!-------------------------------------------------------------------------------
! Routine to update level indices and wetting and drying.
! Use levels1() for better inundation if resolution is fine enough.
!-------------------------------------------------------------------------------

!############################################################################################################

!     Do upwind and TVD transport
      subroutine do_transport_tvd_imp(it,ntr,difnum_max_l) !,nvrt1,npa1,dfh1)
!schism_step 4773        call do_transport_tvd_imp(it,ntracers,difnum_max_l) !,nvrt,npa,dfh)

! it       - timestep-counter  #; info only
! ntracers - number of tracers
!      real(rkind), intent(out) :: difnum_max_l !max. horizontal diffusion number reached by this process (check stability)

! local_to_global_ schism_init: write(10,'(1000(1x,i10))')ns_global,ne_global,np_global,nvrt,nproc,ntracers,ntrs(:) !global info


!#ifdef USE_MPIMODULE
!      use mpi
!#endif
#ifdef USE_QSIM
! flux_adv_vface  !unmodified vertical fluxes (positive upward)
! =0
! dfhm if itur==5 !Tsinghua group:0825 !diff in transport Eq. at nodes & whole levels 1007 | itur=3 is default
! bdy_frc   !body force at prism center Q_{i,k}
! flx_sf    !surface b.c. \kappa*dC/dz = flx_sf (at element center)
! flx_bt    !bottom b.c.
!dfh(nvrt,npa)   - vertikale? Diffusion Turbulenzmodell
!hdif(nvrt,npa)  - horizontale Diffusion nur von hdif.gr3 ??
! delj(ns), 
! iwsett(ntracers), 
! kbs(nsa)

      use schism_glbl, only : ! fixed in
                              !integer,parameter :: rkind = 8      ! Default real datatype
                              rkind,                                 &    = 8
                              ! settling Velocities computation
                              iwsett,wsett,                          &    = 0,0 no settling
                              ! clipping values for salt and temp
                              saltmax,saltmin,tempmax,tempmin,       &    = 0:35 ; 0:35 set by QSim
                              ! Initial tracer conc. at nodes
                              tr_nd0,                                &    = 0 set by QSim
                              !coefficients
                              max_iadjust_mass_consv,                &
                              dtb_min_transport,                     &
                              ! misc_subs.f90: weno1_coef,weno2_coef
                              mnweno1,mnweno2,nweno1,nweno2,         &
                              !weno stuff
                              itvd_e,isten1,isten2,fwts2,            &
                              wts1,wts2,wmat1,wmat2,                 &
                              iside_table,                           &   !a record of all interface sides within the current rank
                              isten_qual2,                           &   !stencil quality, check if at least 1 stencil is on one side of an element side
                              !boundaries
                              itrtype,irange_tr,                     &
                              isbnd,isbe,trth,trobc,                 &
                              ! integer,parameter :: natrm=12 !# of _available_ tracer models at the moment (including T,S)
                              natrm,                                 &
                              ! call aquire_hgrid(.true.)    -------------------
                              ! mesh
                              iegl2(:,:),                            &   ! added to global_to_local.prop
                              ! local_to_global 
                              kbe,ic3(:,:),                          &   ! elements  added
                              elside(:,:),ssign(:,:),area,           &   ! elements  added
                              isdel(:,:),isbs,kbs,                   &   ! sides added
                              distj,delj,snx,sny,                    &   ! sides added
                              nvrt,ntrs(natrm),                      &
                              ne,ielg,ns,i34,isidenode,elnode,       &
                              nsa,nea,                               &   ! added to local_to_global
                              xnd(m),ynd(m),real(dp00(m)),kbp00(m)   &   ! not needed for transport
                              !param.nml &CORE
                              dt,                                    &
                              !param.nml &OPT
                              ihdif,ihconsv,isconsv,itur,itr_met,    &
                              h_tvd,eps1_tvd_imp,eps2_tvd_imp,       &
                              ip_weno,courant_weno,ntd_weno,nquad,   &
                              epsilon1,epsilon2,ielm_transport       &
                              ,                                      &
                              ! timestep in
                              zs,ze,                                       &   levels1 or levels0 calling zcoor
                              ! ELM
                              sdbt,                                        &   avoid by ielm_transport=0
                              !  writeout_nc
                              su2,sv2,eta2,                                &   ++
                              dfh,hdif,                                    &   ++
                              flux_adv_vface,                              &   ++
                              dp,                                          &   ! depth, writeout_nc
                              ! wet/dry flag
                              idry_s,idry_e,idry_e_2t,                     &   ++
                              !=0
                              dfhm,bdy_frc,flx_sf,flx_bt                   &
                              ,                                            &
                              ! timestep out
                              total_mass_error,                   tr_el,         &  ++ temp+salt
                              ! tr_nd !schism_step ab zeile 7603 aus tr_el berechnet
                                                                  tr_nd          &
                              ,                                                  &
                              ! output
                              errmsg
                              

      use schism_msgp, only : nproc,myrank,comm,ierr,rtype                            &
                              parallel_abort,exchange_e2di_2t,exchange_e3d_2t_tr,     &
                              exchange_s3dw,exchange_s3d_tr3,exchange_s3d_tr2
#else
      use schism_glbl
      use schism_msgp
#endif
      use misc_modules


Wyrwa@voss-cln-preprocess:~/test/schism_modelle/sc21_elbe_n17f.vls/outputs> ncdump -h schout_000005_1.nc
        float time(time) ;
                time:i23d = 0 ;
        float wetdry_node(time, nSCHISM_hgrid_node) ;
                wetdry_node:i23d = 1 ;
                wetdry_node:ivs = 1 ;
        float wetdry_elem(time, nSCHISM_hgrid_face) ;
                wetdry_elem:i23d = 4 ;
                wetdry_elem:ivs = 1 ;
        float wetdry_side(time, nSCHISM_hgrid_edge) ;
                wetdry_side:i23d = 7 ;
                wetdry_side:ivs = 1 ;
        float elev(time, nSCHISM_hgrid_node) ;
                elev:i23d = 1 ;
                elev:ivs = 1 ;
        float diffusivity(time, nSCHISM_hgrid_node, nSCHISM_vgrid_layers) ;
                diffusivity:i23d = 2 ;
                diffusivity:ivs = 1 ;
        float hvel_side(time, nSCHISM_hgrid_edge, nSCHISM_vgrid_layers, two) ;
                hvel_side:i23d = 8 ;
                hvel_side:ivs = 2 ;
        float temp_elem(time, nSCHISM_hgrid_face, nSCHISM_vgrid_layers) ;
                temp_elem:i23d = 6 ;
                temp_elem:ivs = 1 ;
        float salt_elem(time, nSCHISM_hgrid_face, nSCHISM_vgrid_layers) ;
                salt_elem:i23d = 6 ;
                salt_elem:ivs = 1 ;
        float flux_adv_vface(time, nSCHISM_hgrid_face, nSCHISM_vgrid_layers) ;
                flux_adv_vface:i23d = 6 ;
                flux_adv_vface:ivs = 1 ;
}
