!-------------------------------------------------------------------------------

!                      Code_Saturne version 7.1
!                      ------------------------
! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2021 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

!===============================================================================
! Function:
! ---------

! Basic example of cs_user_boundary_conditions subroutine.f90
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[out]    icodcl        boundary condition code:
!>                               - 1 Dirichlet
!>                               - 2 Radiative outlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!> \param[in]     itrifb        indirection for boundary faces ordering
!> \param[in,out] itypfb        boundary face types
!> \param[in,out] izfppp        boundary face zone number
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradt \, \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
!_______________________________________________________________________________

subroutine cs_f_user_boundary_conditions &
 ( nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   dt     ,                                                       &
   rcodcl )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use atincl
use ctincl
use cs_fuel_incl
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor), itypfb(nfabor)
integer          izfppp(nfabor)

double precision dt(ncelet)
double precision rcodcl(nfabor,nvar,3)

! Local variables

!< [loc_var_dec]
integer          ifac, iel, ii
integer          ilelt, nlelt
integer          f_id_rough, f_id_t_rough
double precision uref2
double precision rhomoy, xdh
double precision xitur

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer :: bfpro_rom
double precision, dimension(:), pointer :: bpro_roughness
double precision, dimension(:), pointer :: bpro_roughness_t
!< [loc_var_dec]

! inlet profiles 
integer            ny,ny2,ny3,ibc,ix
double precision   dmin,dmin2,dmin3
integer            ituser(5000)
double precision   rtuser(5000)
double precision   uent,vent,went,xkent,xeent,d2s3
! inlet profiles 



!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

!< [init]
allocate(lstelt(nfabor))  ! temporary array for boundary faces selection

call field_get_val_s(ibrom, bfpro_rom)
!< [init]
!
!!===============================================================================
!! Assign boundary conditions to boundary faces here
!
!! For each subset:
!! - use selection criteria to filter boundary faces of a given subset
!! - loop on faces from a subset
!!   - set the boundary condition for each face
!!===============================================================================
!
!! Assign an inlet to boundary faces of group '2' and x < 0.01,
!
!!< [example_1]
!call getfbr('2 and x < 0.01', nlelt, lstelt)
!
!
!do ilelt = 1, nlelt
!
!  ifac = lstelt(ilelt)
!  iel = ifabor(ifac)
!
!  itypfb(ifac) = ientre
!
!  rcodcl(ifac,iu,1) = 1.1d0
!  rcodcl(ifac,iv,1) = 1.1d0
!  rcodcl(ifac,iw,1) = 1.1d0
!
!  uref2 = rcodcl(ifac,iu,1)**2  &
!        + rcodcl(ifac,iv,1)**2  &
!        + rcodcl(ifac,iw,1)**2
!  uref2 = max(uref2,1.d-12)
!
!  !   Turbulence example computed using equations valid for a pipe.
!
!  !   We will be careful to specify a hydraulic diameter adapted
!  !     to the current inlet.
!
!  !   We will also be careful if necessary to use a more precise
!  !     formula for the dynamic viscosity use in the calculation of
!  !     the Reynolds number (especially if it is variable, it may be
!  !     useful to take the law from 'usphyv'. Here, we use by default
!  !     the 'viscl0" value.
!  !   Regarding the density, we have access to its value at boundary
!  !     faces (romb) so this value is the one used here (specifically,
!  !     it is consistent with the processing in 'usphyv', in case of
!  !     variable density)
!
!  ! Hydraulic diameter: to be defined in the notebook GUI
!  xdh = notebook_parameter_value_by_name('hydraulic_diam')
!
!  !   Calculation of turbulent inlet conditions using
!  !     standard laws for a circular pipe
!  !     (their initialization is not needed here but is good practice).
!  rhomoy  = bfpro_rom(ifac)
!
!  call turbulence_bc_inlet_hyd_diam(ifac, uref2, xdh, rhomoy, viscl0,  &
!                                    rcodcl)
!
!  ! Handle scalars
!  if (nscal.gt.0) then
!    do ii = 1, nscal
!      rcodcl(ifac,isca(ii),1) = 1.d0
!    enddo
!  endif
!
!enddo
!< [example_1]

! Assign an inlet to boundary faces of group '3'
call getfbr('VFLAME', nlelt, lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)
  iel  = ifabor(ifac)

  itypfb(ifac) = ientre

  rcodcl(ifac,iu,1) = 0.0d0
  rcodcl(ifac,iv,1) = 0.0d0
  rcodcl(ifac,iw,1) = 0.0d0

  uref2 = rcodcl(ifac,iu,1)**2   &
        + rcodcl(ifac,iv,1)**2   &
        + rcodcl(ifac,iw,1)**2
  uref2 = max(uref2,1.d-12)

  ! Turbulence example computed using turbulence intensity data.

  ! We will be careful to specify a hydraulic diameter adapted
  !   to the current inlet.

  ! Hydraulic diameter: to be defined in the notebook GUI
!  xdh = notebook_parameter_value_by_name('hydraulic_diam')
  xdh = 10.12
  ! Turbulence intensity: to be defined in the notebook GUI
!  xitur = notebook_parameter_value_by_name('turb_intensity')
  xitur = 0.0

  ! Calculation of turbulent inlet conditions using
  !   the turbulence intensity and standard laws for a circular pipe
  !   (their initialization is not needed here but is good practice)

  call turbulence_bc_inlet_turb_intensity(ifac, uref2, xitur, xdh,  &
                                          rcodcl)

  ! --- Handle scalars
  if (nscal.gt.0) then

    ii=1
    rcodcl(ifac,isca(ii),1) = 1.d0
    ii=2
    rcodcl(ifac,isca(ii),1) = 1.d0
    ii=3
    rcodcl(ifac,isca(ii),1) = 0.d0
    ii=4
    rcodcl(ifac,isca(ii),1) = 0.d0

  endif

enddo

!< [example_2]
call getfbr('INLET', nlelt, lstelt)

ibc = 33
ny= 101
ny2= 101
ny3= 101
d2s3=2.d0/3.d0 


if(ntcabs.eq.ntpabs+1)then
  open(file='U_mean.txt',unit=ibc,FORM='formatted')
        do iel=1,ny
          read(ibc,*) rtuser(iel),rtuser(iel+ny)
        enddo
  close(ibc)

  open(file='k.txt',unit=ibc,FORM='formatted')
        do iel=1,ny2
          read(ibc,*) rtuser(iel+(3*ny2)),rtuser(iel+(4*ny2))
        enddo
  close(ibc)

  open(file='epsilon.txt',unit=ibc,FORM='formatted')
        do iel=1,ny3
          read(ibc,*) rtuser(iel+(6*ny3)),rtuser(iel+(7*ny3))
        enddo
  close(ibc)

endif

!PRINT *, "The value of nfabor is:", nfabor

!dmin = 1.d13
do ifac = 1, nlelt
  ituser(ifac) = 0
  ituser(ifac+(5*ny2))=0
  ituser(ifac+(8*ny3))=0
enddo

!PRINT *, "The value of nlelt is:", nlelt


do ilelt = 1, nlelt
  
  dmin  = 1.d13
  dmin2 = 1.d13
  dmin3 = 1.d13
  ifac = lstelt(ilelt)
  iel = ifabor(ifac)
!  PRINT *, "The value of ifac is:",ilelt, ifac

!------for reading mean velocity-------------
  do ix = 1, ny
    if(abs(cdgfbo(2,ifac)-rtuser(ny+ix)).le.dmin)then
      dmin = abs(cdgfbo(2,ifac)-rtuser(ny+ix))
      ituser(ifac) = ix
    endif
  enddo
!------for reading urms---------------------------
  do ix = 1, ny2
    if(abs(cdgfbo(2,ifac)-rtuser((4*ny2)+ix)).le.dmin2)then
      dmin2 = abs(cdgfbo(2,ifac)-rtuser((4*ny2)+ix))
      ituser(ifac+(5*ny2)) = ix
    endif
  enddo
!------for reading epsilon---------------------------
  do ix = 1, ny3
    if(abs(cdgfbo(2,ifac)-rtuser((7*ny3)+ix)).le.dmin3)then
      dmin3 = abs(cdgfbo(2,ifac)-rtuser((7*ny3)+ix))
      ituser(ifac+(8*ny3)) = ix
    endif
  enddo

enddo


do ilelt = 1, nlelt

  ifac = lstelt(ilelt)
  iel = ifabor(ifac)


    itypfb(ifac) = ientre

     uent = rtuser(ituser(ifac))
     vent = 0.d0
     went = 0.d0
    
    rcodcl(ifac,iu,1) = uent
    rcodcl(ifac,iv,1) = vent
    rcodcl(ifac,iw,1) = went

!    uref2 = rcodcl(ifac,iu,1)**2  &
!           +rcodcl(ifac,iv,1)**2  &
!           +rcodcl(ifac,iw,1)**2
!    uref2 = max(uref2,1.d-12)

!    xitur= rtuser(ituser(ifac+(5*ny2))+(3*ny2))!0.02d0

    xkent= rtuser(ituser(ifac+(5*ny2))+(3*ny2))!0.02d0

!    xkent  = 1.5d0*uref2*xitur**2!epzero
!    xeent  = ((sqrt((2.d0/3.d0)*xkent))**(3.d0))/(2.d0)!epzero
!    xeent  = ((xkent)**(3.d0/2.d0))/(0.0254d0)!epzero
     xeent= rtuser(ituser(ifac+(8*ny3))+(6*ny3))!0.02d0



     ! itytur is a flag equal to iturb/10
    if    (itytur.eq.2) then

      rcodcl(ifac,ik,1)  = xkent
      rcodcl(ifac,iep,1) = xeent

    elseif(itytur.eq.3) then

      rcodcl(ifac,ir11,1) = d2s3*xkent
      rcodcl(ifac,ir22,1) = d2s3*xkent
      rcodcl(ifac,ir33,1) = d2s3*xkent
      rcodcl(ifac,ir12,1) = 0.d0
      rcodcl(ifac,ir13,1) = 0.d0
      rcodcl(ifac,ir23,1) = 0.d0
      rcodcl(ifac,iep,1)  = xeent

!    elseif(iturb.eq.50) then

!      rcodcl(ifac,ik,1)   = xkent
!      rcodcl(ifac,iep,1)  = xeent
!      rcodcl(ifac,iphi,1) = d2s3
!      rcodcl(ifac,ifb,1)  = 0.d0

    elseif(iturb.eq.60) then

      rcodcl(ifac,ik,1)   = xkent
      rcodcl(ifac,iomg,1) = xeent/cmu/xkent

    endif

  ! --- Handle scalars
  if (nscal.gt.0) then
      ii = 1
      rcodcl(ifac,isca(ii),1) = 0.d0
      ii = 2
      rcodcl(ifac,isca(ii),1) = 0.d0
      ii = 3
      rcodcl(ifac,isca(ii),1) = 0.d0
      ii = 4
      rcodcl(ifac,isca(ii),1) = 0.d0

  endif

enddo
!< [example_2]

! Assign an outlet to boundary faces of group 'outlet'

!< [example_3]
call getfbr('OUTLET', nlelt, lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Outlet: zero flux for velocity and temperature, prescribed pressure
  !         Note that the pressure will be set to P0 at the first
  !         free outlet face (isolib)

  itypfb(ifac) = isolib

enddo
!< [example_3]

! Assign a wall to boundary faces of group '5'

!< [example_4]
call getfbr('TOP_WALL', nlelt, lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Wall: zero flow (zero flux for pressure)
  !       friction for velocities (+ turbulent variables)
  !       zero flux for scalars

  itypfb(ifac) = iparoi

  ! If sliding wall with velocity u = 1:
  ! rcodcl(ifac, iu, 1) = 1.d0

  ! If sliding wall with velocity u = 0: nothing to do

  if (nscal.gt.0) then
!!!
!!!    ! If temperature prescribed to 20 with wall law (scalar ii=1):
!!!    ii = 1
!!!    icodcl(ifac, isca(ii))    = 5
!!!    rcodcl(ifac, isca(ii), 1) = 20.d0
!!!
!!!    ! If temperature prescribed to 50 with no wall law (simple Dirichlet)
!!!    !   with exchange coefficient 8 (scalar ii=2):
!!!    ii = 2
!!!    icodcl(ifac, isca(ii))    = 1
!!!    rcodcl(ifac, isca(ii),1)  = 50.d0
!!!    rcodcl(ifac, isca(ii), 2) = 8.d0
!!!
! If flux prescribed to 4.d0 (scalar ii=3):
    ii=1
    icodcl(ifac, isca(ii))    = 3
    rcodcl(ifac, isca(ii), 3) = 0.d0

! If temperature prescribed to 20 with wall law (scalar ii=1):
    ii = 2
    icodcl(ifac, isca(ii))    = 5
    rcodcl(ifac, isca(ii), 1) = 0.d0

! If RPV_variance prescribed to 0 (scalar ii=3):
    ii = 3
    icodcl(ifac, isca(ii))    = 3
    rcodcl(ifac, isca(ii), 3) = 0.d0

! If Theta_variance prescribed to 0 (scalar ii=4):
    ii = 4
    icodcl(ifac, isca(ii))    = 5
    rcodcl(ifac, isca(ii), 1) = 0.d0

   endif
enddo
!< [example_4]

call getfbr('BOTTOM_WALL', nlelt, lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Wall: zero flow (zero flux for pressure)
  !       friction for velocities (+ turbulent variables)
  !       zero flux for scalars

  itypfb(ifac) = iparoi

  ! If sliding wall with velocity u = 1:
  ! rcodcl(ifac, iu, 1) = 1.d0

  ! If sliding wall with velocity u = 0: nothing to do

  if (nscal.gt.0) then
    
    ii = 1
    icodcl(ifac, isca(ii))    = 3
    rcodcl(ifac, isca(ii), 3) = 0.d0

    ii = 2
    icodcl(ifac, isca(ii))    = 5
    rcodcl(ifac, isca(ii), 1) = 0.d0

! If RPV_variance prescribed to 0 (scalar ii=3):
    ii = 3
    icodcl(ifac, isca(ii))    = 3
    rcodcl(ifac, isca(ii), 3) = 0.d0

! If Theta_variance prescribed to 0 (scalar ii=4):
    ii = 4
    icodcl(ifac, isca(ii))    = 5
    rcodcl(ifac, isca(ii), 1) = 0.d0


!!!
!!!    ! If temperature prescribed to 50 with no wall law (simple Dirichlet)
!!!    !   with exchange coefficient 8 (scalar ii=2):
!!!    ii = 2
!!!    icodcl(ifac, isca(ii))    = 1
!!!    rcodcl(ifac, isca(ii),1)  = 50.d0
!!!    rcodcl(ifac, isca(ii), 2) = 8.d0
!!!
!!!    ! If flux prescribed to 4.d0 (scalar ii=3):
!!!    ii = 3
!!!    icodcl(ifac, isca(ii))    = 3
!!!    rcodcl(ifac, isca(ii), 3) = 4.d0
!!!
  endif
enddo
!< [example_4]


!< [finalize]
deallocate(lstelt)  ! temporary array for boundary faces selection
!< [finalize]

return
end subroutine cs_f_user_boundary_conditions
