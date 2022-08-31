

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
                              rkind,                                 &
                              ! settling Velocities computation
                              iwsett,wsett,                          &
                              ! clipping values for salt and temp
                              saltmax,saltmin,tempmax,tempmin,       &
                              ! Initial tracer conc. at nodes
                              tr_nd0,                                &
                              !coefficients
                              max_iadjust_mass_consv,                &
                              dtb_min_transport,                     &
                              ! misc_subs.f90: weno1_coef,weno2_coef
                              mnweno1,mnweno2,nweno1,nweno2,         &
                              !weno stuff
                              itvd_e,isten1,isten2,fwts2,            &
                              wts1,wts2,wmat1,wmat2,                 &
                              !boundaries
                              itrtype,irange_tr,                     &
                              isbnd,isbe,trth,trobc,                 &
                              ! integer,parameter :: natrm=12 !# of _available_ tracer models at the moment (including T,S)
                              natrm,                                 &
                              ! schism_init
                              delj,                                  &
                              ! call aquire_hgrid(.true.) 
                              distj,area,isbs,snx,sny,               &
                              ! mesh
                              iegl2,kbs,kbe,elside,ic3,nsa,nea,dp,   &
                              ssign,isten_qual2,isdel,iside_table,   &
                              ! local_to_global 
                              nvrt,ntrs(natrm),                      &
                              ne,ielg,ns,i34,isidenode,elnode,       &
                              !param.nml &CORE
                              dt,                                    &
                              !param.nml &OPT
                              ihdif,ihconsv,isconsv,itur,itr_met,    &
                              h_tvd,eps1_tvd_imp,eps2_tvd_imp,       &
                              ip_weno,courant_weno,ntd_weno,nquad,   &
                              epsilon1,epsilon2,ielm_transport       &
                              ,                                      &
                              ! timestep in
                              ! wet/dry flag                               &
                              idry_s,idry_e,idry_e_2t,                     &
                              ! ELM
                              sdbt,                                        &
                              !  writeout_nc
                              su2,sv2,eta2,zs,ze,dfh,hdif,                 &
                              flux_adv_vface,                              &
                              !=0
                              dfhm,bdy_frc,flx_sf,flx_bt                   &
                              ,                                            &
                              ! timestep out
                              total_mass_error,                   tr_el,         &
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
