!===============================================================================
! User source terms definition.
!
! 1) Momentum equation (coupled solver)
! 2) Species transport
! 3) Turbulence (Rij-epsilon)
!===============================================================================

!-------------------------------------------------------------------------------

!                      code_saturne version 7.1
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
!> \file cs_user_source_terms.f90
!>
!> \brief User subroutines for additional right-hand side source terms
!>
!> See \ref cs_user_source_terms and
!> \ref cs_user_source_terms-scalar_in_a_channel for examples.
!>
!===============================================================================
!> \brief Additional right-hand side source terms for velocity components equation
!> (Navier-Stokes)
!>
!> \deprecated Use \ref cs_user_source_terms instead.
!>
!> \section ustsnv_use  Usage
!>
!> The additional source term is decomposed into an explicit part (\c crvexp) and
!> an implicit part (\c crvimp) that must be provided here.
!> The resulting equation solved by the code for a velocity is:
!> \f[
!>  \rho \norm{\vol{\celli}} \DP{\vect{u}} + ....
!>   = \tens{crvimp} \cdot \vect{u} + \vect{crvexp}
!> \f]
!>
!> Note that \c crvexp and \c crvimp are defined after the Finite Volume integration
!> over the cells, so they include the "volume" term. More precisely:
!>   - crvexp is expressed in kg.m/s2
!>   - crvimp is expressed in kg/s
!>
!> The \c crvexp and \c crvimp arrays are already initialized to 0
!> before entering the
!> the routine. It is not needed to do it in the routine (waste of CPU time).
!>
!> \remark The additional force on \f$ x_i \f$ direction is given by
!>  \c crvexp(i, iel) + vel(j, iel)* crvimp(j, i, iel).
!>
!> For stability reasons, Code_Saturne will not add -crvimp directly to the
!> diagonal of the matrix, but Max(-crvimp,0). This way, the crvimp term is
!> treated implicitely only if it strengthens the diagonal of the matrix.
!> However, when using the second-order in time scheme, this limitation cannot
!> be done anymore and -crvimp is added directly. The user should therefore test
!> the negativity of crvimp by himself.
!>
!> When using the second-order in time scheme, one should supply:
!>   - crvexp at time n
!>   - crvimp at time n+1/2
!>
!> The selection of cells where to apply the source terms is based on a
!> \ref getcel command. For more info on the syntax of the \ref getcel command,
!> refer to the user manual or to the comments on the similar command
!> \ref getfbr in the routine \ref cs_user_boundary_conditions.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss terms
!> \param[in]     ncesmp        number of cells with mass source terms
!> \param[in]     ivar          index number of the current variable
!> \param[in]     icepdc        index number of cells with head loss terms
!> \param[in]     icetsm        index number of cells with mass source terms
!> \param[in]     itypsm        type of mass source term for each variable
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        head loss coefficient
!> \param[in]     smacel        value associated to each variable in the mass
!>                              source terms or mass rate
!> \param[out]    crvexp        explicit part of the source term
!> \param[out]    crvimp        implicit part of the source term
!_______________________________________________________________________________

subroutine ustsnv &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel ,                                              &
   crvexp , crvimp )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          ivar

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)
double precision crvexp(3,ncelet), crvimp(3,3,ncelet)

! Local variables
integer          iel, ipcrom, ipp
double precision ckp, qdm
integer, allocatable, dimension(:) :: lstelt

!===============================================================================


!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate a temporary array for cells selection
allocate(lstelt(ncel))

ckp  = 0.d0
qdm  = 0.d0

!do iel = 1, ncel
!   crvimp(1,1,iel) = - volume(iel)*propce(iel,ipcrom)*ckp
!enddo

do iel = 1, ncel
   crvexp(1,iel) =   volume(iel)*qdm
enddo
!--------
! Formats
!--------

 1000 format(' User source terms for variable ', a8,/)

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine ustsnv


!===============================================================================

!===============================================================================
!>    User subroutine.
!> \brief    Additional right-hand side source terms for scalar equations (user
!>     scalars and specific physics scalars).
!>
!> \deprecated Use \ref cs_user_source_terms instead.
!>
!> Usage
!> -----
!> The routine is called for each scalar, user or specific physisc. It is
!> therefore necessary to test the value of the scalar number iscal to separate
!> the treatments of the different scalars (if (iscal.eq.p) then ....).
!>
!> The additional source term is decomposed into an explicit part (crvexp) and
!> an implicit part (crvimp) that must be provided here.
!> The resulting equation solved by the code for a scalar f is:
!>
!>   \f[ \rho*volume*\frac{df}{dt} + .... = crvimp*f + crvexp \f]
!>
!>
!> Note that crvexp and crvimp are defined after the Finite Volume integration
!> over the cells, so they include the "volume" term. More precisely:
!>   - crvexp is expressed in kg.[scal]/s, where [scal] is the unit of the scalar
!>   - crvimp is expressed in kg/s
!>
!>
!> The crvexp and crvimp arrays are already initialized to 0 before entering the
!> the routine. It is not needed to do it in the routine (waste of CPU time).
!>
!> For stability reasons, Code_Saturne will not add -crvimp directly to the
!> diagonal of the matrix, but Max(-crvimp,0). This way, the crvimp term is
!> treated implicitely only if it strengthens the diagonal of the matrix.
!> However, when using the second-order in time scheme, this limitation cannot
!> be done anymore and -crvimp is added directly. The user should therefore test
!> the negativity of crvimp by himself.
!>
!> When using the second-order in time scheme, one should supply:
!>   - crvexp at time n
!>   - crvimp at time n+1/2
!>
!>
!> The selection of cells where to apply the source terms is based on a getcel
!> command. For more info on the syntax of the getcel command, refer to the
!> user manual or to the comments on the similar command \ref getfbr in the routine
!> \ref cs_user_boundary_conditions.

!> WARNING: If scalar is the temperature, the resulting equation
!>          solved by the code is:
!>
!>  rho*Cp*volume*dT/dt + .... = crvimp*T + crvexp
!>
!>
!> Note that crvexp and crvimp are defined after the Finite Volume integration
!> over the cells, so they include the "volume" term. More precisely:
!>   - crvexp is expressed in W
!>   - crvimp is expressed in W/K
!>

!>
!> STEP SOURCE TERMS
!>===================
!> In case of a complex, non-linear source term, say F(f), for scalar f, the
!> easiest method is to implement the source term explicitely.
!>
!>   df/dt = .... + F(f(n))
!>   where f(n) is the value of f at time tn, the beginning of the time step.
!>
!> This yields :
!>   crvexp = volume*F(f(n))
!>   crvimp = 0
!>
!> However, if the source term is potentially steep, this fully explicit
!> method will probably generate instabilities. It is therefore wiser to
!> partially implicit the term by writing:
!>
!>   df/dt = .... + dF/df*f(n+1) - dF/df*f(n) + F(f(n))
!>
!> This yields:
!>   crvexp = volume*( F(f(n)) - dF/df*f(n) )
!>   crvimp = volume*dF/df
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss terms
!> \param[in]     ncesmp        number of cells with mass source terms
!> \param[in]     iscal         index number of the current scalar
!> \param[in]     icepdc        index number of cells with head loss terms
!> \param[in]     icetsm        index number of cells with mass source terms
!> \param[in]     itypsm        type of mass source term for each variable
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        head loss coefficient
!> \param[in]     smacel        value associated to each variable in the mass
!> \param[in]                   source terms or mass rate
!> \param[out]    crvexp        explicit part of the source term
!> \param[out]    crvimp        implicit part of the source term
!______________________________________________________________________________!


subroutine ustssc &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel ,                                              &
   crvexp , crvimp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use mesh
use field
use field_operator 
use cs_c_bindings
use user_module 
use cs_f_interfaces


!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)
double precision crvexp(ncelet), crvimp(ncelet)


!double precision ebu_c
double precision temp_u, temp_b, temp_q 
double precision Pe_Ql, delta_z, delta_th
double precision dist_q
! Local variables
!character(len=80) :: chaine
integer    ivar,iel
integer, allocatable, dimension(:) :: lstelt


double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, dimension(:), pointer ::  crom
double precision, dimension(:), pointer :: cvara_k, cvara_ep, cvara_omg
double precision, dimension(:,:), pointer :: cvara_rij
double precision, dimension(:), pointer :: cvara_scal
double precision, dimension(:), pointer :: omega_c
double precision, dimension(:), pointer :: theta
double precision, dimension(:), pointer :: w_dist
double precision, dimension(:), pointer :: temp             ! dimensional temperature 
double precision, dimension(:), pointer :: cpro_visct 




double precision, allocatable, dimension(:) :: cvara_scal_R ! Reynolds averaged progress variable
double precision, dimension(:), pointer :: ly_bml
double precision, dimension(:), pointer :: I_0          ! strech factor 
double precision, dimension(:), pointer :: Ka_L,Da_L,SDR_c,C3,C4,SDR_t
double precision, dimension(:), pointer :: diff_c_T
double precision, dimension(:), pointer :: RPV_variance, theta_variance

double precision, dimension(:), pointer :: countergrad_sum 




double precision tau_hr, lfs_sl,Re_tau


double precision Le
double precision B_SDR,PI_SDR 
double precision, allocatable, dimension(:) :: A1,A2,A3,A_eps
double precision, allocatable, dimension(:) :: psi
double precision, allocatable, dimension(:,:) :: gradc
double precision, allocatable, dimension(:,:) :: gradt 

double precision, allocatable, dimension(:)   :: sum_T_2c
double precision, allocatable, dimension(:)   :: sum_T_2t
!type(var_cal_opt) :: vcopt

integer   f_id, inc, iccocg, iprev 

double precision, allocatable, dimension(:,:,:) :: grad_cc
double precision, allocatable, dimension(:,:,:) :: grad_tt 



!===============================================================================



!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate a temporary array for cells selection
!allocate(lstelt(ncel))

! Allocate temporary arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))


! --- Numero du scalaire a traiter : ISCAL

! --- Numero de la variable associee au scalaire a traiter ISCAL
ivar = isca(iscal)

! --- Nom de la variable associee au scalaire a traiter ISCAL
!call field_get_label(ivarfl(ivar), chaine)

! --- Numero des grandeurs physiques (voir cs_user_boundary_conditions)
call field_get_val_s(icrom, crom)

!if ( ivar.eq.isca(1) ) then
!  call field_get_val_prev_s(ivarfl(isca(iscal)), cvara_scal)
!endif

call field_get_val_prev_s_by_name("RPV",cvara_scal)



if (itytur.eq.2.or.iturb.eq.50) then
  call field_get_val_prev_s(ivarfl(ik), cvara_k)
  call field_get_val_prev_s(ivarfl(iep), cvara_ep)
elseif (itytur.eq.3) then
  call field_get_val_prev_v(ivarfl(irij), cvara_rij)
  call field_get_val_prev_s(ivarfl(iep), cvara_ep)
elseif (iturb.eq.60) then
  call field_get_val_prev_s(ivarfl(ik), cvara_k)
  call field_get_val_prev_s(ivarfl(iomg), cvara_omg)
endif


!===============================================================================
! 2. PRISE EN COMPTE DES TERMES SOURCES
!===============================================================================

!if ( ivar.eq.isca(1) ) then

! ---> Terme source pour la fraction massique moyenne de gaz frais

!  if (vcopt%iwarni.ge.1) then
!    write(nfecra,1000) chaine(1:8)
!  endif

! ---> Calcul de K et Epsilon en fonction du modele de turbulence

  if (itytur.eq.2) then

    do iel = 1, ncel
      w1(iel) = cvara_k(iel)
      w2(iel) = cvara_ep(iel)
    enddo

  elseif (itytur.eq.3) then

    do iel = 1, ncel
      w1(iel) = 0.5d0 *( cvara_rij(1,iel)                    &
                        +cvara_rij(2,iel)                    &
                        +cvara_rij(3,iel) )
      w2(iel) = cvara_ep(iel)
    enddo

  elseif (iturb.eq.50) then

    do iel = 1, ncel
      w1(iel) = cvara_k(iel)
      w2(iel) = cvara_ep(iel)
    enddo

  elseif (iturb.eq.60) then

    do iel = 1, ncel
      w1(iel) = cvara_k(iel)
      w2(iel) = cmu*cvara_k(iel)*cvara_omg(iel)
    enddo

  endif



!===============================================================================
! 3. COMBUSTION MODEL CALCULATION
!===============================================================================
!NLL93 Hardcode change the varialbes as global SOON
Re_tau=1.0d0/viscl0
tau_hr=2.3d0
lfs_sl=0.7d0
temp_u=730.0d0
temp_b=2409.d0 
temp_q=780.0d0
Pe_Ql=2.2d0
delta_z=1.d0/Re_tau/0.7d0/lfs_sl  !TODO 110 Re_tau 0.7 Pr 0.7 lfs 
delta_th=0.04d0

Le=1.0d0
dist_q=0.01d0        ! 0.01d0

PI_SDR=Pe_Ql*erf((8.d0*Le-6.d0)/2.d0)
B_SDR = -6.d0*(Le-1.0d0)

!call field_get_val_s(ro0,cpro_ro0)
call field_get_val_prev_s_by_name("Theta",theta)

call field_get_val_s_by_name("source_term",omega_c)
call field_get_val_s_by_name("Temperature",temp)
call field_get_val_s_by_name("wall_distance",w_dist)
call field_get_val_s_by_name("countergrad_sum",countergrad_sum)


do iel =1,ncel
   temp(iel)=(1.d0-theta(iel))*temp_u+theta(iel)*temp_b
enddo 




! no model
if (ico_model.eq.0) then
  if ( isca(iscal).eq.isca(1) ) then
    do iel = 1, ncel
      crvexp(iel) =0.0d0
      crvimp(iel) =0.0d0
      omega_c(iel)=crvexp(iel)+crvimp(iel)*cvara_scal(iel) 
    enddo
  endif


! EBU
elseif (ico_model.eq.1) then
  if ( isca(iscal).eq.isca(1) ) then       
    do iel = 1, ncel
     
     if (temp(iel).le.temp_q) then 
      crvexp(iel)=0.0d0 
      crvimp(iel)=0.0d0
      omega_c(iel)=0.0d0 
     else
      w3(iel)=c_ebu*w2(iel)/w1(iel)*crom(iel)
      crvexp(iel)=w3(iel)*cvara_scal(iel)**2*volume(iel)
      crvimp(iel)=w3(iel)*(-2.d0*cvara_scal(iel)+1.d0)*volume(iel)
      omega_c(iel)=w3(iel)*cvara_scal(iel)*(1-cvara_scal(iel)) 
     endif 
    enddo
  endif


! BML
elseif(ico_model.eq.2) then
  allocate(cvara_scal_R(ncelet))
 
  call field_get_val_s_by_name("Ly",ly_bml)

  call field_get_val_s_by_name("I_0",I_0)

  do iel = 1, ncel
     I_0(iel)=1.0d0*(erf(w_dist(iel)/delta_z-1.4d0*Pe_Ql)+1.d0)
  enddo



  if ( isca(iscal).eq.isca(1) ) then
    do iel = 1, ncel
    
    if (temp(iel).le.temp_q .or.w_dist(iel).le.dist_q ) then
      crvexp(iel)=0.0d0 
      crvimp(iel)=0.0d0
      omega_c(iel)=0.0d0
    else 
      ly_bml(iel)=max(cl_bml*(lfs_sl/sqrt(2.0d0*w1(iel)/3.0d0))**n_bml &
                  *sqrt(w1(iel))**3/w2(iel), 1.d-8)
    
         
      cvara_scal_R(iel)=(1.d0+tau_hr)*cvara_scal(iel)/(1.d0+tau_hr*cvara_scal(iel))
!      crvexp(iel)=bml_g*(1.d0-cvara_scal_R(iel))*cvara_scal_R(iel)/bml_sigmay/bml_ly(iel)*volume(iel) ! Density
      w3(iel)=I_0(iel)*ro0*lfs_sl*g_bml/sigmay_bml/ly_bml(iel)
      crvexp(iel)=w3(iel)*cvara_scal_R(iel)**2.d0*volume(iel)+countergrad_sum(iel)*volume(iel)
      crvimp(iel)=w3(iel)*(-2.d0*cvara_scal_R(iel)+1.d0)*volume(iel) 
      omega_c(iel)=(crvexp(iel)+crvimp(iel)*cvara_scal(iel))/(volume(iel)+1.d-8)

   endif
   enddo
  endif
  ! Free memory 
  deallocate(cvara_scal_R)


! SDR 

elseif(ico_model.eq.3) then
  call field_get_val_s_by_name("Da_L",Da_L)
  call field_get_val_s_by_name("Ka_L",Ka_L)
  call field_get_val_s_by_name("C3",C3)
  call field_get_val_s_by_name("C4",C4)
  call field_get_val_s_by_name("SDR_c",SDR_c)

  if ( isca(iscal).eq.isca(1) ) then
    do iel = 1, ncel
    
    Ka_L(iel) = max(sqrt(delta_th*w2(iel)/(lfs_sl**3)),1.d-8) 
    Da_L(iel) = max(lfs_sl * w1(iel)/w2(iel)/delta_th,1.d-8)
    C3(iel) = 1.5d0*sqrt(Ka_L(iel))/(1.d0+sqrt(Ka_L(iel)))
    C4(iel) = 1.1d0/((1.d0+Ka_L(iel))**0.4d0)

    if (temp(iel).le.temp_q  .or.w_dist(iel).le.dist_q ) then ! 0.003
      crvexp(iel)=0.0d0 
      crvimp(iel)=0.0d0
      SDR_c(iel) =0.0d0
      omega_c(iel)=0.0d0
    else  

      SDR_c(iel)=max((2.d0*K_c_tcp*lfs_sl/delta_th+(w2(iel)/w1(iel))*(C3(iel)-tau_hr*C4(iel)*Da_L(iel)))/beta_prime , 1.d-8)         

      w3(iel)=2.d0*crom(iel)*SDR_c(iel)/(2*c_m-1.d0)
      crvexp(iel)=w3(iel)*cvara_scal(iel)**2.d0*volume(iel)
      crvimp(iel)=w3(iel)*(-2.d0*cvara_scal(iel)+1.d0)*volume(iel) 
     
      SDR_c(iel)=SDR_c(iel)*cvara_scal(iel)*(1.d0-cvara_scal(iel))
      omega_c(iel)=2.d0*crom(iel)*SDR_c(iel)/(2*c_m-1.d0)
      endif
    enddo
  endif

  ! SDR flame wall interaction

elseif(ico_model.eq.31) then
  call field_get_val_s_by_name("Da_L",Da_L)
  call field_get_val_s_by_name("Ka_L",Ka_L)
  call field_get_val_s_by_name("C3",C3)
  call field_get_val_s_by_name("C4",C4)
  call field_get_val_s_by_name("SDR_c",SDR_c)
  call field_get_val_s_by_name("Difference_c_T",diff_c_T)


  allocate(A1(ncelet), A2(ncelet), A3(ncelet), A_eps(ncelet))
  allocate(psi(ncelet))


  if ( isca(iscal).eq.isca(1) ) then
    do iel = 1, ncel
    Ka_L(iel) = max(sqrt(delta_th*w2(iel)/(lfs_sl**3)),1.d-8) 
    Da_L(iel) = max(lfs_sl * w1(iel)/w2(iel)/delta_th,1.d-8)
    C3(iel) = 1.5d0*sqrt(Ka_L(iel))/(1.d0+sqrt(Ka_L(iel)))
    C4(iel) = 1.1d0/((1.d0+Ka_L(iel))**0.4d0)



!    diff_c_T(iel)=0.0d0                        !TODO remove this after fix the parallel for the diff_c_T

!    psi(iel)=(max(5.d0*max(diff_c_T(iel),0.d0),1.0d0))**0.3d0
    psi(iel)=(max(5.d0*max(diff_c_T(iel),0.d0),1.0d0))**0.3d0

!    psi(iel)  = (max(0.d0,1.0d0))**0.3d0
    A1(iel)   = 0.5d0*(erf(3.0d0*(w_dist(iel)/delta_z-PI_SDR))+1.d0)
    A2(iel)   = 0.5d0*(erf(w_dist(iel)/delta_z-psi(iel)*PI_SDR)+1.0d0)
    A3(iel)   = 2.31d0*erf(2.6d0*(max(cvara_scal(iel)-theta(iel),1.d-8)))   
    A_eps(iel)= 0.5d0*(erf(w_dist(iel)/delta_z-PI_SDR)+1.0)

    if (temp(iel).le.temp_q  .or.w_dist(iel).le.dist_q ) then ! 0.003
      crvexp(iel)=0.0d0 
      crvimp(iel)=0.0d0
      SDR_c(iel)=0.0d0 
      omega_c(iel)=0.0d0
    else  

      w3(iel)=A_eps(iel)*exp(-1.2d0*Le*(diff_c_T(iel))**3.d0) &
              * max((2.d0*K_c_tcp*lfs_sl/delta_th+(w2(iel)/w1(iel))*(C3(iel)-tau_hr*C4(iel)*Da_L(iel)))/beta_prime , 1.d-8)         

      SDR_c(iel)=w3(iel)*cvara_scal(iel)**2.d0+w3(iel)*(-2.d0*cvara_scal(iel)+1.d0)*cvara_scal(iel)

      ! TODO fix the Re_tau and Pr 
      omega_c(iel)=(2.d0*crom(iel)*SDR_c(iel)/(2*c_m-1.d0)*A1(iel)*exp(Le*max(cvara_scal(iel)-theta(iel),1.d-8))+ &
                     ro0*lfs_sl*A2(iel)*A3(iel) &
                     *sqrt(crom(iel)*Re_tau*0.7d0*SDR_c(iel))*exp(-0.5*(w_dist(iel)/delta_z-PI_SDR)**2)/(Le**B_SDR)) 

      crvexp(iel)=omega_c(iel)*volume(iel)+countergrad_sum(iel)
      crvimp(iel)=0.d0 
     endif
    enddo
  endif
  deallocate(A1, A2, A3,A_eps)
  deallocate(psi)
endif


! add source term to theta, the second scalar 
if ( isca(iscal).eq.isca(2) ) then       
    do iel = 1, ncel    
      crvexp(iel)=omega_c(iel)*volume(iel)
      crvimp(iel)=0.0
    enddo
endif


! add  countergradient term as a source term



! add source term to the third scalar, RPV variance and the fourth scalar, theta variance 

if (ico_variance.eq.1) then 

iprev  = 1
inc    = 1
iccocg = 1

allocate(A_eps(ncelet))
allocate(gradc(3,ncelet))
allocate(gradt(3,ncelet))
allocate(sum_T_2c(ncelet))
allocate(sum_T_2t(ncelet))

call field_get_val_s(iprpfl(ivisct),cpro_visct)


call field_get_id('RPV',f_id)
call field_gradient_scalar(f_id,iprev,0,inc, iccocg ,gradc)
call field_get_id('Theta',f_id)
call field_gradient_scalar(f_id,iprev,0,inc, iccocg ,gradt)
call field_get_val_s_by_name("SDR_t",SDR_t)
call field_get_val_prev_s_by_name("RPV_variance",RPV_variance)
call field_get_val_prev_s_by_name("Theta_variance", theta_variance)


!const int ksigmas=cs_field_key_id("turbulent_schmidt")

!const cs_real_t turb_schmidt=cs_field_get_key_double(f,ksigmas); 




  do iel = 1, ncel 
    sum_T_2c(iel)=-crom(iel)*cpro_visct(iel)*(gradc(1,iel)/0.70d0)*gradc(1,iel) &
                  -crom(iel)*cpro_visct(iel)*(gradc(2,iel)/0.70d0)*gradc(2,iel) &
                  -crom(iel)*cpro_visct(iel)*(gradc(3,iel)/0.70d0)*gradc(3,iel) 

    sum_T_2t(iel)=-crom(iel)*cpro_visct(iel)*(gradt(1,iel)/0.70d0)*gradt(1,iel) &
                  -crom(iel)*cpro_visct(iel)*(gradt(2,iel)/0.70d0)*gradt(2,iel) &
                  -crom(iel)*cpro_visct(iel)*(gradt(3,iel)/0.70d0)*gradt(3,iel) 

  enddo 




  if ( isca(iscal).eq.isca(3) ) then       
    do iel = 1, ncel    
           
    if (temp(iel).le.temp_q  .or.w_dist(iel).le.dist_q ) then ! 0.003
      crvexp(iel)=0.0d0 
      crvimp(iel)=0.0d0
    else
!            print *,'omega_c', omega_c(iel)
!            print *, 'cvara_scal', cvara_scal(iel)
      crvexp(iel)=(2.0d0*omega_c(iel)*(c_m-cvara_scal(iel))&
                  -2.0d0*crom(iel)*SDR_c(iel)-2.0d0*sum_T_2c(iel))*volume(iel)
      crvimp(iel)=0.0d0
    endif 
    enddo
  endif

  if ( isca(iscal).eq.isca(4) ) then       

  do iel = 1, ncel 
    A_eps(iel)= 0.5d0*(erf(w_dist(iel)/delta_z-PI_SDR)+1.0) ! need to be consistent with the RPV
    if (temp(iel).le.temp_q  .or.w_dist(iel).le.dist_q ) then ! 0.003
      SDR_t(iel)=0.0d0 
    else

      w3(iel)=(1.0d0-erf(max(diff_c_T(iel),0.d0)**0.3))*A_eps(iel)*exp(-1.2d0*Le*(diff_c_T(iel))**3.d0) &
      * max((2.d0*K_c_tcp*lfs_sl/delta_th+(w2(iel)/w1(iel))*C3(iel)- tau_hr*C4(iel)*lfs_sl/delta_th)/beta_prime , 1.d-8) &     
      + erf(max(diff_c_T(iel),0.d0)**0.3)*theta_variance(iel)*w2(iel)/w1(iel)
      SDR_t(iel)=w3(iel)*theta(iel)**2.d0+w3(iel)*(-2.d0*theta(iel)+1.d0)*theta(iel)
    endif 
  enddo 


    do iel = 1, ncel    
    
     if (temp(iel).le.temp_q  .or.w_dist(iel).le.dist_q ) then ! 0.003
      crvexp(iel)=0.0d0 
      crvimp(iel)=0.0d0
    else
 
  
       crvexp(iel)=(2.0d0*omega_c(iel)*(c_m-cvara_scal(iel))*sqrt(max(RPV_variance(iel),1.d-2)/max(theta_variance(iel),1.d-2)) &
                 -2.0d0*crom(iel)*SDR_t(iel)-2.d0*sum_T_2t(iel))*volume(iel)
 
       crvimp(iel)=0.0
    endif 
    enddo
  endif

deallocate(gradc)
deallocate(gradt)
deallocate(sum_T_2c)
deallocate(sum_T_2t)
deallocate(A_eps)

endif 


! Free memory 
deallocate(w1, w2, w3)



!--------
! Formats
!--------

 1000 format(' User source terms for variable ',A8,/)

!----
! End
!----

! Deallocate the temporary array
!deallocate(lstelt)

return
end subroutine ustssc

!===============================================================================

!===============================================================================
!> \brief    Additional right-hand side source terms for vectorial equations
!>           (user vectors and specific physics vectors).
!>
!> \deprecated Use \ref cs_user_source_terms instead.
!>
!> Usage
!> -----
!> The routine is called for each vector, user or specific physisc. It is
!> therefore necessary to test the value of the vector number iscal to separate
!> the treatments of the different vectors (if (iscal.eq.p) then ....).
!>
!> The additional source term is decomposed into an explicit part (crvexp) and
!> an implicit part (crvimp) that must be provided here.
!> The resulting equation solved by the code for a vector f is:
!>
!>   \f[ \rho*volume*\frac{d\vect{f}}{dt} + .... = \tens{crvimp}*\vect{f} +
!>                                                 \vect{crvexp} \f]
!>
!>
!> Note that crvexp and crvimp are defined after the Finite Volume integration
!> over the cells, so they include the "volume" term. More precisely:
!>   - crvexp is expressed in kg.[scal]/s, where [scal] is vector unit
!>   - crvimp is expressed in kg/s
!>
!>
!> The crvexp and crvimp arrays are already initialized to 0 before entering the
!> the routine. It is not needed to do it in the routine (waste of CPU time).
!>
!> For stability reasons, Code_Saturne will not add -crvimp directly to the
!> diagonal of the matrix, but Max(-crvimp,0). This way, the crvimp term is
!> treated implicitely only if it strengthens the diagonal of the matrix.
!> However, when using the second-order in time scheme, this limitation cannot
!> be done anymore and -crvimp is added directly. The user should therefore test
!> the negativity of crvimp by himself.
!>
!> When using the second-order in time scheme, one should supply:
!>   - crvexp at time n
!>   - crvimp at time n+1/2
!>
!>
!> The selection of cells where to apply the source terms is based on a getcel
!> command. For more info on the syntax of the getcel command, refer to the
!> user manual or to the comments on the similar command \ref getfbr in the routine
!> \ref cs_user_boundary_conditions.

!> WARNING: If scalar is the temperature, the resulting equation
!>          solved by the code is:
!>
!>  rho*Cp*volume*dT/dt + .... = crvimp*T + crvexp
!>
!>
!> Note that crvexp and crvimp are defined after the Finite Volume integration
!> over the cells, so they include the "volume" term. More precisely:
!>   - crvexp is expressed in W
!>   - crvimp is expressed in W/K
!>

!>
!> STEEP SOURCE TERMS
!>===================
!> In case of a complex, non-linear source term, say F(f), for scalar f, the
!> easiest method is to implement the source term explicitely.
!>
!>   df/dt = .... + F(f(n))
!>   where f(n) is the value of f at time tn, the beginning of the time step.
!>
!> This yields :
!>   crvexp = volume*F(f(n))
!>   crvimp = 0
!>
!> However, if the source term is potentially steep, this fully explicit
!> method will probably generate instabilities. It is therefore wiser to
!> partially implicit the term by writing:
!>
!>   df/dt = .... + dF/df*f(n+1) - dF/df*f(n) + F(f(n))
!>
!> This yields:
!>   crvexp = volume*( F(f(n)) - dF/df*f(n) )
!>   crvimp = volume*dF/df
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss terms
!> \param[in]     ncesmp        number of cells with mass source terms
!> \param[in]     iscal         index number of the current scalar
!> \param[in]     icepdc        index number of cells with head loss terms
!> \param[in]     icetsm        index number of cells with mass source terms
!> \param[in]     itypsm        type of mass source term for each variable
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        head loss coefficient
!> \param[in]     smacel        value associated to each variable in the mass
!> \param[in]                   source terms or mass rate
!> \param[out]    crvexp        explicit part of the source term
!> \param[out]    crvimp        implicit part of the source term
!______________________________________________________________________________!


subroutine ustsvv &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel ,                                              &
   crvexp , crvimp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)
double precision crvexp(3,ncelet), crvimp(3,3,ncelet)

! Local variables

integer, allocatable, dimension(:) :: lstelt

!===============================================================================



!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate a temporary array for cells selection
allocate(lstelt(ncel))

!--------
! Formats
!--------

 1000 format(' User source terms for variable ',A8,/)

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine ustsvv

!===============================================================================

!===============================================================================
!> \brief Additional right-hand side source terms for turbulence models and
!>irijco =1
!>
!> \deprecated Use \ref cs_user_source_terms instead.
!>
!> \section cs_user_rij_source_terms_use  Usage
!>
!> The additional source term is decomposed into an explicit part (crvexp) and
!> an implicit part (crvimp) that must be provided here.
!> The resulting equations solved by the code are:
!> \f[
!>  \rho \norm{\vol{\celli}} \DP{\varia} + ....
!>   = \tens{crvimp} \varia + \vect{crvexp}
!> \f]
!> where \f$ \varia \f$ is the turbulence field of index \c f_id
!>
!> Note that crvexp, crvimp are defined after the Finite Volume
!> integration over the cells, so they include the "volume" term. More precisely:
!>   - crvexp is expressed in kg.m2/s2
!>   - crvimp is expressed in kg/s
!>
!> The crvexp, crvimp arrays are already initialized to 0 before
!> entering the routine. It is not needed to do it in the routine (waste of CPU time).
!>
!> For stability reasons, Code_Saturne will not add -crvimp directly to the
!> diagonal of the matrix, but Max(-crvimp,0). This way, the crvimp term is
!> treated implicitely only if it strengthens the diagonal of the matrix.
!> However, when using the second-order in time scheme, this limitation cannot
!> be done anymore and -crvimp is added directly. The user should therefore test
!> the negativity of crvimp by himself.
!>
!> When using the second-order in time scheme, one should supply:
!>   - crvexp at time n
!>   - crvimp at time n+1/2
!>
!> The selection of cells where to apply the source terms is based on a getcel
!> command. For more info on the syntax of the \ref getcel command, refer to the
!> user manual or to the comments on the similar command \ref getfbr in the routine
!> \ref cs_user_boundary_conditions.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss terms
!> \param[in]     ncesmp        number of cells with mass source terms
!> \param[in]     f_id          field index of the current turbulent variable
!> \param[in]     icepdc        index number of cells with head loss terms
!> \param[in]     icetsm        index number of cells with mass source terms
!> \param[in]     itypsm        type of mass source term for each variable
!> \param[in]     ckupdc        head loss coefficient
!> \param[in]     smacel        value associated to each variable in the mass
!>                              source terms or mass rate
!> \param[out]    crvexp        explicit part of the source term
!> \param[out]    crvimp        implicit part of the source term
!_______________________________________________________________________________

subroutine cs_user_turbulence_source_terms2 &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   f_id   ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ckupdc , smacel ,                                              &
   crvexp , crvimp )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use mesh
use field
use cs_f_interfaces
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          f_id

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)
double precision crvexp(6,ncelet), crvimp(6,6,ncelet)

! Local variables

integer, allocatable, dimension(:) :: lstelt

!===============================================================================



!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate a temporary array for cells selection
allocate(lstelt(ncel))

!--------
! Formats
!--------

 1000 format(' User source terms for turbulence model',/)

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine cs_user_turbulence_source_terms2
