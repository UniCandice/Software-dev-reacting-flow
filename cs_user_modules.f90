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
! Purpose:
! -------

!> \file cs_user_modules-user-arrays.f90
!>
!> \brief User-defined module example
!>
!> This file is compiled before all other user Fortran files.
!> To ensure this, it must not be renamed.
!>
!> The user may define an arbitrary number of modules here, even though
!> only one is defined in the example.
!
!> \cond DOXYGEN_SHOULD_SKIP_THIS

!-------------------------------------------------------------------------------

module user_module

  !=============================================================================

  implicit none

  !=============================================================================
!< [variables for combustion]
  !nll93 options for combustion models 
  ! 0: no model
  ! 1: EBU model
  ! 2: BML model  
  ! 3: SDR standard 
  ! 31: SDR_FWI 
  integer, parameter :: ico_model = 2
!EBU
  double precision   :: c_ebu     = 2.5d0
!BML
  double precision   :: g_bml     = 1.5d0
  double precision   :: sigmay_bml= 0.5d0
  double precision   :: cl_bml    = 1.0d0
  double precision   :: n_bml     = 1.0d0   

!SDR
  double precision   :: K_c_tcp    = 0.76d0*2.3d0  ! 4.5 the heat release parameter 
  double precision   :: beta_prime = 6.7d0
  double precision   :: c_m        = 0.844d0

  integer, parameter :: ico_variance =0
  integer, parameter :: ico_countergrad = 0

!< [variables for combustion]


contains

! subroutine for calculating the C_w-T_w at the wall and linked it to the values of the cell 
subroutine  cs_f_user_C_T_discrepancy_wall &
 ( nvar   , nscal  ,                                              &
   icodcl ,  rcodcl )

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

double precision rcodcl(nfabor,nvar,3)

! Local variables
integer       ii,ifac,iel
integer       ilelt,nlelt


integer, allocatable, dimension(:) :: lstelt

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

!< [init]
allocate(lstelt(nfabor))  ! temporary array for boundary faces selection

!< [init]

call getfbr('TOP_WALL', nlelt, lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)


  if (nscal.gt.0) then
    ii=1
    print *,'C_w=', rcodcl(ifac, isca(ii), 1)

    ii = 2
    print *,'T_w=', rcodcl(ifac, isca(ii), 1)



   endif
enddo



        
end subroutine cs_f_user_C_T_discrepancy_wall



end module user_module

!-------------------------------------------------------------------------------

!> (DOXYGEN_SHOULD_SKIP_THIS) \endcond
