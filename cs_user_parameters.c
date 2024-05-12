/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* code_saturne version 7.1 */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_parameters.c
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \ref parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select physical model options, including user fields.
 *
 * This function is called at the earliest stages of the data setup,
 * so field ids are not available yet.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_model(void)
{
/*--------------------------------------------------------------------------*/

  /* Example: Chose a turbulence model
   *   CS_TURB_NONE: no turbulence model (laminar flow)
   *   CS_TURB_MIXING_LENGTH: mixing length model
   *   CS_TURB_K_EPSILON: standard k-epsilon model
   *   CS_TURB_K_EPSILON_LIN_PROD: k-epsilon model with
   *     Linear Production (LP) correction
   *   CS_TURB_K_EPSILON_LS: Launder-Sharma low Re k-epsilon model
   *   CS_TURB_K_EPSILON_QUAD: Baglietto et al. low Re
   *   CS_TURB_RIJ_EPSILON_LRR: Rij-epsilon (LRR)
   *   CS_TURB_RIJ_EPSILON_SSG: Rij-epsilon (SSG)
   *   CS_TURB_RIJ_EPSILON_EBRSM: Rij-epsilon (EBRSM)
   *   CS_TURB_LES_SMAGO_CONST: LES (constant Smagorinsky model)
   *   CS_TURB_LES_SMAGO_DYN: LES ("classical" dynamic Smagorisky model)
   *   CS_TURB_LES_WALE: LES (WALE)
   *   CS_TURB_V2F_PHI: v2f phi-model
   *   CS_TURB_V2F_BL_V2K: v2f BL-v2-k
   *   CS_TURB_K_OMEGA: k-omega SST
   *   CS_TURB_SPALART_ALLMARAS: Spalart-Allmaras model */

  cs_turb_model_t *turb_model = cs_get_glob_turb_model();
  turb_model->iturb = CS_TURB_K_EPSILON_LS;

  /* Advanced choice of Wall function:
   *  Wall function for the dynamic, iwallf can be:
   *   CS_WALL_F_DISABLED: no wall function
   *   CS_WALL_F_1SCALE_POWER: deprecated 1 velocity scale power law
   *   CS_WALL_F_1SCALE_LOG: 1 velocity scale log law
   *   CS_WALL_F_2SCALES_LOG: 2 velocity scales log law (default for many models)
   *   CS_WALL_F_SCALABLE_2SCALES_LOG: Scalable log wll function with two
   *                                   velocity scales
   *   CS_WALL_F_2SCALES_VDRIEST: 2 velocity scales with Van Driest damping
   *   CS_WALL_F_2SCALES_SMOOTH_ROUGH: 2 velocity scales valid for smooth and
   *                                   rough regimes
   *   CS_WALL_F_2SCALES_CONTINUOUS: 2 velcoity scales (all y+), default for
   *                                 low Reynolds turbulence models.
   */
/*  {
    cs_wall_functions_t *wf = cs_get_glob_wall_functions();
     wf->iwallf = CS_WALL_F_SCALABLE_2SCALES_LOG;
  } */

   cs_parameters_add_variable("RPV",1);


   cs_parameters_add_variable("Theta",1); 

   cs_parameters_add_variable("RPV_variance",1);


   cs_parameters_add_variable("Theta_variance",1); 


   cs_parameters_add_property("Temperature",
                             1,
                             CS_MESH_LOCATION_CELLS);
   cs_parameters_add_property("RPV_bar",
                             1,
                             CS_MESH_LOCATION_CELLS);


   cs_parameters_add_property("source_term",
                             1,
                             CS_MESH_LOCATION_CELLS);


   cs_parameters_add_property("sflux_RPV",
                             3,
                             CS_MESH_LOCATION_CELLS);

   cs_parameters_add_property("sflux_Theta",
                             3,
                             CS_MESH_LOCATION_CELLS);


   int ico_model_c=2;
 


  if (ico_model_c == 2){
   cs_parameters_add_property("Ly",
                             1,
                             CS_MESH_LOCATION_CELLS);
 
   cs_parameters_add_property("I_0",
                             1,
                             CS_MESH_LOCATION_CELLS);

  }

  if (ico_model_c == 3 || ico_model_c ==31 ){
   cs_parameters_add_property("Da_L",
                             1,
                             CS_MESH_LOCATION_CELLS);
 
   cs_parameters_add_property("Ka_L",
                             1,
                             CS_MESH_LOCATION_CELLS);

   cs_parameters_add_property("C3",
                             1,
                             CS_MESH_LOCATION_CELLS);

   cs_parameters_add_property("C4",
                             1,
                             CS_MESH_LOCATION_CELLS);

   cs_parameters_add_property("SDR_c",
                             1,
                             CS_MESH_LOCATION_CELLS);
  }


  if (ico_model_c==31 ){

   cs_parameters_add_property("Difference_c_T",
		              1,
			      CS_MESH_LOCATION_CELLS);


  }

  if (ico_model_c==3 || ico_model_c==31 ){


   cs_parameters_add_property("SDR_t",
                             1,
                             CS_MESH_LOCATION_CELLS);

  }

 


   cs_parameters_add_property("countergrad_sum",
                             1,
                             CS_MESH_LOCATION_CELLS);

  
   cs_parameters_add_property("sflux_countergrad_scalar",
                             1,
                             CS_MESH_LOCATION_CELLS);



  /*--------------------------------------------------------------------------*/
}


 


/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify general numerical and physical user parameters.
 *
 * At the calling point of this function, most model-related most variables
 * and other fields have been defined, so specific settings related to those
 * fields may be set here.
 *
 * At this stage, the mesh is not built or read yet, so associated data
 * such as field values are not accessible yet, though pending mesh
 * operations and some fields may have been defined.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parameters(cs_domain_t   *domain)
{
 /* Time step type */

  cs_time_step_options_t *time_opt = cs_get_glob_time_step_options();
  time_opt->idtvar = CS_TIME_STEP_CONSTANT;

  /* Reference time step dt_ref
     The example given below is probably not adapted to your case. */

  cs_real_t dt_ref = 0.002;
  domain->time_step->dt_ref = dt_ref;

  /* Duration
     nt_max is absolute number of the last time step required
     if we have already run 10 time steps and want to
     run 10 more, nt_max must be set to 10 + 10 = 20 */

  domain->time_step->nt_max = 1000/*(int) (40./dt_ref)*/;

  /* Example: change Reference fluid properties options */
  /*----------------------------------------------------*/

  /* Members of the structure cs_fluid_properties_t
   */

  {
    cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();

    /*
      ro0        : density in kg/m3
      viscl0     : dynamic viscosity in kg/(m s)
      cp0        : specific heat in J/(Kelvin kg)
      t0         : reference temperature in Kelvin
      p0         : total reference pressure in Pascal
                   the calculation is based on a
                   reduced pressure P*=Ptot-ro0*g.(x-xref)
                   (except in compressible case)
      xyzp0(3)   : coordinates of the reference point for
                   the total pressure (where it is equal to p0)

      In general, it is not necessary to furnish a reference point xyz0.
      If there are outlets, the code will take the center of the
      reference outlet face.
      On the other hand, if we plan to explicitly fix Dirichlet conditions
      for pressure, it is better to indicate to which reference the
      values relate (for a better resolution of reduced pressure).

      Other properties are given by default in all cases.

      Nonetheless, we may note that:

      In the standard case (no combustion, electric arcs, compressibility):
      ---------------------
      ro0, viscl0 and cp0
          are useful and represent either the fluid properties if they
          are constant, either simple mean values for the initialization
          if properties are variable and defined in cs_user_physical_properties.
      t0  is not useful
      p0  is useful but is not used in an equation of state. p0
          is a reference value for the incompressible solver
          which will serve to set the (possible) domain outlet pressure.
          We may also take it as 0 or as a physical value in Pascals.

      With the electric module:
      ------------------------
      ro0, viscl0 and cp0
          are useful but simply represent mean initial values;
          the density, molecular dynamic viscosity, and specific
          heat are necessarily defined as fields (whether they are
          physically variable or not): see cs_user_physical_properties
          for the Joule effect
          module and the electric arcs dp_ELE data file.
      t0  is useful an must be in Kelvin (> 0) but represents a simple
          initialization value.
      p0  is useful bu is not used in the equation of state. p0
          is a reference value for the incompressible solver which
          will be used to calibrate the (possible) outlet pressure
          of the domain. We may take it as zero or as a physical
          value in Pascals.

      With gas combustion:
      --------------------
      ro0 is not useful (it is automatically recalculated by the
          law of ideal gases from t0 and p0).
      viscl0 is indispensable: it is the molecular dynamic viscosity,
          assumed constant for the fluid.
      cp0 is indispensable: it is the heat capacity, assumed constant,
          (modelization of source terms involving a local Nusselt in
          the Lagrangian module, reference value allowing the
          calculation of a radiative
          (temperature, exchange coefficient) couple).
      t0  is indispensible and must be in Kelvin (> 0).
      p0  is indispensable and must be in Pascal (> 0).

      With pulverized coal:
      ---------------------
      ro0 is not useful (it is automatically recalculated by the
          law of ideal gases from t0 and p0).
      viscl0 is indispensable: it is the molecular dynamic viscosity,
          assumed constant for the fluid (its effect is expected to
          be small compared to turbulent effects).
      cp0 is indispensable: it is the heat capacity, assumed constant,
          (modelization of source terms involving a local Nusselt in
          the coal or Lagrangian module, reference value allowing the
          calculation of a radiative
          (temperature, exchange coefficient) couple).
      t0  is indispensable and must be in Kelvin (> 0).
      p0  is indispensable and must be in Pascal (> 0).

      With compressibility:
      ---------------------
      ro0 is not useful, stricto sensu; nonetheless, as experience
          shows that users often use this variable, it is required
          to assign to it a strictly positive value (for example,
          an initial value).
      viscl0 is useful and represents the molecular dynamic viscosity,
          when it is constant, or a value which will be used during
          initializations (or in inlet turbulence conditions,
          depending on the user choice.
      cp0 is indispensable: it is the heat capacity, assumed constant
          in the thermodynamics available by default
      t0  is indispensable and must be in Kelvin (> 0).
      p0  is indispensable and must be in Pascal (> 0).
          With the thermodynamic law available by default,
      t0 and p0 are used for the initialization of the density.
      xyzp0 is not useful because the pressure variable directly
          represents the total pressure.
    */

    fp->ro0    = 1.0;
    fp->viscl0 = 1.0/110;
    fp->cp0    = 1.0;

    fp->t0 = 20. + 273.15;
    fp->p0 = 1.0;

    /* We only specify XYZ0 if we explicitely fix Dirichlet conditions
       for the pressure. */

    fp->xyzp0[0] = 0.;
    fp->xyzp0[1] = 0.;
    fp->xyzp0[2] = 0.;

    /* irovar, ivivar, icp: constant or variable density,
                              viscosity/diffusivity, and specific heat

         When a specific physics module is active
           (coal, combustion, electric arcs, compressible: see usppmo)
           we MUST NOT set variables 'irovar', 'ivivar', and 'icp' here, as
           they are defined automatically.
         Nonetheless, for the compressible case, ivivar may be modified
           in the uscfx2 user subroutine.

         When no specific physics module is active, we may specify if the
           density, specific heat, and the molecular viscosity
           are constant (irovar=0, ivivar=0, icp=-1), which is the default
           or variable (irovar=1, ivivar=1, icp=0)

         For those properties we choose as variable, the corresponding law
           must be defined in cs_user_physical_properties
           (incs_user_physical_properties.f90);
           if they are constant, they take values ro0, viscl0, and cp0.
    */

    fp->irovar = 1;
    fp->ivivar = 0;
    fp->icp = -1;

  }
  /* Example: Change physical constants
   *------------------------------------*/

  {
    cs_physical_constants_t *pc = cs_get_glob_physical_constants();

    pc->gravity[0] = 0.;
    pc->gravity[1] = 0.;
    pc->gravity[2] = 0.; /* gravity  (m/s2) in the z direction */
  }

  /* Example: Change options relative to the inner iterations
   * over prediction-correction.
   * - nterup: number of sub-iterations (default 1)
   * - epsup: relative precision (default 10e-5)
   *------------------------------------*/

  {
    cs_velocity_pressure_param_t *vp_param = cs_get_glob_velocity_pressure_param();
    vp_param->nterup = 4;
  }

  /* Example: activate the porous model
   * - 0 No porosity taken into account (Default)
   * - 1 Porosity taken into account
   *------------------------------------*/

  {
    cs_glob_porous_model = 0;
  }
  /* Example: log verbosity */
  /*------------------------*/

  {
    int n_fields = cs_field_n_fields();

    for (int f_id = 0; f_id < n_fields; f_id++) {

      cs_field_t  *f = cs_field_by_id(f_id);

      if (f->type & CS_FIELD_VARIABLE) {
        cs_equation_param_t *eqp = cs_field_get_equation_param(f);
        eqp->verbosity = 1;
      }
    }
  }

  /*Convective scheme

    blencv = 0 for upwind (order 1 in space, "stable but diffusive")
           = 1 for centered/second order (order 2 in space)
      we may use intermediate real values.
      Here we choose:
        for the velocity and user scalars:
          an upwind-centered scheme with 100% centering (blencv=1)
        for other variables
          the default code value (upwind standard, centered in LES)

    Specifically, for user scalars
      if we suspect an excessive level of numerical diffusion on
        a variable ivar representing a user scalar
        iscal (with ivar=isca(iscal)), it may be useful to set
        blencv = 1.0to use a second-order scheme in space for
        convection. For temperature or enthalpy in particular, we
        may thus choose in this case:

        cs_field_t *f = cs_thermal_model_field();
        cs_equation_param_t *eqp = cs_field_get_equation_param(f);
        eqp->blencv = 1.;

      For non-user scalars relative to specific physics
        implicitly defined by the model,
        the corresponding information is set automatically elsewhere:
        we do not modify blencv here. */

  {
    int n_fields = cs_field_n_fields();

    for (int f_id = 0; f_id < n_fields; f_id++) {

      cs_field_t *f = cs_field_by_id(f_id);

      if (f->type & CS_FIELD_VARIABLE) {
        cs_equation_param_t *eqp = cs_field_get_equation_param(f);
        eqp->blencv = 1;
      }
    }
  }

  /* Linear solver parameters (for each unknown)
     epsilo: relative precision for the solution of the linear system. */
  {
    int n_fields = cs_field_n_fields();

    for (int f_id = 0; f_id < n_fields; f_id++) {

      cs_field_t  *f = cs_field_by_id(f_id);

      if (f->type & CS_FIELD_VARIABLE) {
        cs_equation_param_t *eqp = cs_field_get_equation_param(f);
        eqp->epsilo = 1.e-6;
      }
    }
  }

  /* Dynamic reconstruction sweeps to handle non-orthogonlaities
     This parameter computes automatically a dynamic relax factor,
     and can be activated for any variable.
      - iswdyn = 0: no relaxation
      - iswdyn = 1: means that the last increment is relaxed
      - iswdyn = 2: (default) means that the last two increments are used
                    to relax.
  */
  {
    cs_equation_param_t *eqp = cs_field_get_equation_param(CS_F_(p));
    eqp->iswdyn = 2;
  }


  /* Example: Minimum and maximum admissible valuesfor each USER scalar
   * results are clipped at the end of each time step.
   * 
   * If min > max, we do not clip
   *
   * For a scalar jj representing the variance of another, we may
   * abstain from defining these values
   * (a default clipping is set in place).
   * This is the purpose of the test on iscavr(jj) in the example below.
   * 
   * For non-user scalars relative to specific physics (coal, combustion,
   * electric arcs: see usppmo) implicitly defined according to the
   * model, the information is automatically set elsewhere: we
   * do not set min or max values here. */

	 
   {
   cs_field_t*sca1 = cs_field_by_name("RPV");





   /* We define the min and max bounds */
   cs_field_set_key_double(sca1,
                           cs_field_key_id("min_scalar_clipping"),
                           0.0);
   cs_field_set_key_double(sca1,
                           cs_field_key_id("max_scalar_clipping"),
                           1.0);
   }


   {
   cs_field_t*sca1 = cs_field_by_name("RPV_variance");





   /* We define the min and max bounds */
   cs_field_set_key_double(sca1,
                           cs_field_key_id("min_scalar_clipping"),
                           0.0);
   cs_field_set_key_double(sca1,
                           cs_field_key_id("max_scalar_clipping"),
                           0.25);
   }




   {
   cs_field_t*sca1 = cs_field_by_name("Theta");





   /* We define the min and max bounds */
   cs_field_set_key_double(sca1,
                           cs_field_key_id("min_scalar_clipping"),
                           0.0);
   cs_field_set_key_double(sca1,
                           cs_field_key_id("max_scalar_clipping"),
                           1.0);
   } 


   {
   cs_field_t*sca1 = cs_field_by_name("Theta_variance");





   /* We define the min and max bounds */
   cs_field_set_key_double(sca1,
                           cs_field_key_id("min_scalar_clipping"),
                           0.0);
   cs_field_set_key_double(sca1,
                           cs_field_key_id("max_scalar_clipping"),
                           0.25);
   }


  /* Example: add boundary values for all scalars */
  /*----------------------------------------------*/

  {
    int n_fields = cs_field_n_fields();

    for (int f_id = 0; f_id < n_fields; f_id++) {

      cs_field_t  *f = cs_field_by_id(f_id);

      if (f->type & CS_FIELD_VARIABLE)
        cs_parameters_add_boundary_values(f);

    }
  }


}
/*----------------------------------------------------------------------------*/
/*!
 * \brief Define linear solver options.
 *
 * This function is called at the setup stage, once user and most model-based
 * fields are defined.
 *
 * Available native iterative linear solvers include conjugate gradient,
 * Jacobi, BiCGStab, BiCGStab2, and GMRES. For symmetric linear systems,
 * an algebraic multigrid solver is available (and recommended).
 *
 * External solvers may also be setup using this function, the cs_sles_t
 * mechanism allowing such through user-define functions.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_linear_solvers(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define time moments.
 *
 * This function is called at the setup stage, once user and most model-based
 * fields are defined, and before fine control of field output options
 * is defined.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_time_moments(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define internal coupling options.
 *
 * Options are usually defined using cs_internal_coupling_add_entity.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_internal_coupling(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify output user parameters.
 *
 * For CDO schemes, this function concludes the setup of properties,
 * equations, source terms...
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_finalize_setup(cs_domain_t   *domain)
{
  CS_UNUSED(domain);

  /* Frequency of log output */
  cs_glob_log_frequency = 1;

  cs_field_set_key_int(CS_F_(rho),cs_field_key_id("post_vis"),CS_POST_ON_LOCATION|CS_POST_MONITOR);
  cs_field_set_key_int(CS_F_(rho),cs_field_key_id("log"),CS_POST_ON_LOCATION|CS_POST_MONITOR);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
