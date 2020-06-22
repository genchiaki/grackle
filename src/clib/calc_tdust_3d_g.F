#include "phys_const.def"

!=======================================================================
!\////////////////////  SUBROUTINE CALC_TDUST_3D  \\\\\\\\\\\\\\\\\\\\\\

      subroutine calc_tdust_3d_g(
     &     d, de, HI, HII, 
     &     HeI, HeII, HeIII,
     &     HM, H2I, H2II, 
     &     in, jn, kn, 
     &     nratec, iexpand,
     &     ispecies, idim,
     &     is, js, ks, 
     &     ie, je, ke, 
     &     aye, temstart, temend,
     &     fgr, gasgra,
     &     gamma_isrfa, isrf,
     &     utem, uxyz, uaye,
     &     urho, utim,
     &     gas_temp, dust_temp,
     &     iisrffield, isrf_habing
#ifdef SMBH_RAD
     &   , dust2_temp
     &   , ibhrffield, bhrf, S_Silicate, S_Graphite
     &   , kp_gr_Size, kp_gr_N, kp_gr_dTd, kp_gr_Td
     &   , kp_Silicate, kp_Graphite
#endif
     & )

!  COMPUTE THE DUST TEMPERATURE
!
!  written by: Britton Smith
!  date: July, 2011
!  modified1: 
!
!  PURPOSE:
!    Calculate dust heat balance to get the dust temperature.
!
!  INPUTS:
!
!  PARAMETERS:
!
!-----------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"
#ifdef _OPENMP
#include "omp_lib.h"
#endif

!  Arguments

      integer in, jn, kn, is, js, ks, ie, je, ke, nratec,
     &        iexpand, ispecies, idim, iisrffield
#ifdef SMBH_RAD
     &      , ibhrffield
#endif
      real*8  aye, temstart, temend,
     &        utem, uxyz, uaye, urho, utim, fgr, isrf
      R_PREC  d(in,jn,kn),
     &     de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &     HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &     HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &     gas_temp(in,jn,kn), dust_temp(in,jn,kn),
     &     isrf_habing(in,jn,kn)
#ifdef SMBH_RAD
     &   , dust2_temp(in,jn,kn)
     &   , bhrf(in,jn,kn)
      real*8  S_Silicate, S_Graphite
      integer kp_gr_Size, kp_gr_N(1)
      real*8  kp_gr_dTd, kp_gr_Td(kp_gr_N(1))
     &      , kp_Silicate(kp_gr_Size), kp_Graphite(kp_gr_Size)
#endif

!  Chemistry tables
      real*8  gasgra(nratec), gamma_isrfa

!  Parameters

      real*8 mh
      parameter (mh = mass_h)

!  Locals

      integer i, j, k
      integer t, dj, dk
      real*8  trad, zr, logtem0, logtem9, dlogtem,
     &        coolunit, dbase1, tbase1, xbase1

!  Slice locals
 
      integer indixe(in)
      R_PREC  t1(in), t2(in), logtem(in), tdef(in), 
     &     tgas(in), tdust(in), nh(in), gasgr(in),
     &     myisrf(in)
#ifdef SMBH_RAD
     &   , mybhrf(in)
     &   , tdust2(in)
#endif
          
!  Iteration mask for multi_cool

      logical itmask(in)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Set log values of start and end of lookup tables

      logtem0  = log(temstart)
      logtem9  = log(temend)
      dlogtem  = (log(temend) - log(temstart))/real(nratec-1)

!     Set units

      tbase1   = utim
      xbase1   = uxyz/(aye*uaye)    ! uxyz is [x]*a      = [x]*[a]*a'        '
      dbase1   = urho*(aye*uaye)**3 ! urho is [dens]/a^3 = [dens]/([a]*a')^3 '
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
      zr       = 1._RKIND/(aye*uaye) - 1._RKIND

!     Set CMB temperature

      trad = 2.73_RKIND * (1._RKIND + zr)

!  Loop over zones, and do an entire i-column in one go
      dk = ke - ks + 1
      dj = je - js + 1

! parallelize the k and j loops with OpenMP
! flat j and k loops for better parallelism
#ifdef _OPENMP
!$omp parallel do schedule(runtime) private(
!$omp&   i, j, k,
!$omp&   indixe,
!$omp&   t1, t2, logtem, tdef,
!$omp&   tgas, tdust, nh, gasgr, myisrf,
!$omp&   itmask )
#endif
      do t = 0, dk*dj-1
        k = t/dj      + ks+1
        j = mod(t,dj) + js+1

        do i = is+1, ie+1

!     Set itmask to all true

           itmask(i) = .true.

!     Compute interstellar radiation field

           if (iisrffield .gt. 0) then
              myisrf(i) = isrf_habing(i,j,k)
           else
              myisrf(i) = isrf
           endif
#ifdef SMBH_RAD
           if (ibhrffield .gt. 0) then
              mybhrf(i) = bhrf(i,j,k)
           else
              mybhrf(i) = 0._DKIND
           endif
#endif

!     Compute hydrogen number density

           nh(i) = HI(i,j,k) + HII(i,j,k)
           if (ispecies .gt. 1) then
              nh(i) = nh(i) + H2I(i,j,k) + H2II(i,j,k)
           endif

!     We have not converted to proper, so use urho and not dom

           nh(i) = nh(i) * urho / mh

!     Compute log temperature and truncate if above/below table max/min

           tgas(i)   = gas_temp(i,j,k)
           logtem(i) = log(tgas(i))
           logtem(i) = max(logtem(i), logtem0)
           logtem(i) = min(logtem(i), logtem9)

!     Compute index into the table and precompute parts of linear interp

           indixe(i) = min(nratec-1,
     &          max(1,int((logtem(i)-logtem0)/dlogtem,DIKIND)+1))
           t1(i) = (logtem0 + (indixe(i) - 1)*dlogtem)
           t2(i) = (logtem0 + (indixe(i)    )*dlogtem)
           tdef(i) = (logtem(i) - t1(i)) / (t2(i) - t1(i))

!     Lookup values and do a linear temperature in log(T)
!     Convert back to cgs

           gasgr(i) = gasgra(indixe(i)) + tdef(i)
     &          *(gasgra(indixe(i)+1) -gasgra(indixe(i)))
#ifdef SMBH_RAD
           gasgr(i) = gasgr(i) * coolunit / mh
#else
           gasgr(i) = gasgr(i) * fgr * coolunit / mh
#endif

        enddo

#ifdef SMBH_RAD
!     --- Compute silicate temperature in a slice ---

        call calc_tdust_1d_g(tdust, tgas, nh, gasgr,
     &       gamma_isrfa, myisrf,
     &       itmask, trad, in, is, ie, j, k
     &     , mybhrf, S_Silicate
     &     , kp_gr_Size, kp_gr_N, kp_gr_dTd, kp_gr_Td, kp_Silicate)

!     --- Compute graphite temperature in a slice ---

        call calc_tdust_1d_g(tdust2, tgas, nh, gasgr,
     &       gamma_isrfa, myisrf,
     &       itmask, trad, in, is, ie, j, k
     &     , mybhrf, S_Graphite
     &     , kp_gr_Size, kp_gr_N, kp_gr_dTd, kp_gr_Td, kp_Graphite)
#else
!     --- Compute dust temperature in a slice ---

        call calc_tdust_1d_g(tdust, tgas, nh, gasgr,
     &       gamma_isrfa, myisrf,
     &       itmask, trad, in, is, ie, j, k)
#endif

!     Copy slice values back to grid

        do i = is+1, ie+1
           dust_temp(i,j,k) = tdust(i)
#ifdef SMBH_RAD
           dust2_temp(i,j,k) = tdust2(i)
!!         write(*,*) 'hhh', tdust(i), tdust2(i)
#endif
        enddo

      enddo

      return
      end
