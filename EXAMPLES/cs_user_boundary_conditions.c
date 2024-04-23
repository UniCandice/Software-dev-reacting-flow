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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_boundary_conditions.c
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \ref parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity inlet profile as an analytic function.
 *         elt_ids is optional. If not NULL, it enables to access to the coords
 *         array with an indirection. The same indirection can be applied to
 *         the val array if dense_output is set to false.
 *
 * \param[in]      time          when ?
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids (in coords and retval)
 * \param[in]      coords        where ? Coordinates array
 * \param[in]      dense_output  perform an indirection in res or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] val           resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_vel_profile(cs_real_t           time,
             cs_lnum_t           n_elts,
             const cs_lnum_t    *elt_ids,
             const cs_real_t    *coords,
             bool                dense_output,
             void               *input,
             cs_real_t          *val)
{
  CS_UNUSED(time);
  CS_UNUSED(input);

  const cs_real_3_t *elt_coords = (const cs_real_3_t *)coords;

  cs_real_3_t  *v = (cs_real_3_t *)val;

  for (cs_lnum_t i = 0; i < n_elts; i++) {

    const cs_lnum_t  elt_id = (elt_ids == NULL) ? i : elt_ids[i];
    const cs_lnum_t  j = dense_output ? i : elt_id;

    cs_real_t  y = elt_coords[elt_id][1] - 0.5;
    cs_real_t  z = elt_coords[elt_id][2] - 0.5;
    cs_real_t r = sqrt(y*y + z*z);

    v[j][0] = 1.5 * (1 - r);
    v[j][1] = 0;
    v[j][2] = 0;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the scalar inlet profile as an analytic function.
 *         elt_ids is optional. If not NULL, it enables to access to the coords
 *         array with an indirection. The same indirection can be applied to
 *         the val array if dense_output is set to false.
 *
 * \param[in]      time          when ?
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids (in coords and retval)
 * \param[in]      coords        where ? Coordinates array
 * \param[in]      dense_output  perform an indirection in res or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] val           resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_scalar_profile(cs_real_t           time,
                cs_lnum_t           n_elts,
                const cs_lnum_t    *elt_ids,
                const cs_real_t    *coords,
                bool                dense_output,
                void               *input,
                cs_real_t          *val)
{
  CS_UNUSED(time);

  const cs_real_3_t *elt_coords = (const cs_real_3_t *)coords;

  const cs_field_t *f = input;  /* field pointer passed as input
                                   upon assignment */

  if (strcmp(f->name, "scalar1") == 0) {

    for (cs_lnum_t i = 0; i < n_elts; i++) {

      const cs_lnum_t  elt_id = (elt_ids == NULL) ? i : elt_ids[i];
      const cs_lnum_t  j = dense_output ? i : elt_id;

      cs_real_t  y = elt_coords[elt_id][1] - 0.5;
      cs_real_t  z = elt_coords[elt_id][2] - 0.5;
      cs_real_t r = sqrt(y*y + z*z);

      val[j] = 4.0 * (1 - r);

    }
  }
}

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set boundary conditions to be applied.
 *
 * This function is called just before \ref cs_user_finalize_setup, and
 * boundary conditions can be defined in either of those functions,
 * depending on whichever is considered more readable or practical for a
 * given use.
 *
 * \param[in, out]  domain  pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions_setup(cs_domain_t  *domain)
{
  CS_UNUSED(domain);

  {
    cs_equation_param_t  *eqp = cs_equation_param_by_name("velocity");

    cs_real_t inlet_velocity[] = {1, 0, 0};

    cs_equation_add_bc_by_value(eqp,
                                CS_PARAM_BC_DIRICHLET,
                                "inlet",           // zone name
                                inlet_velocity);

  }

  {
    cs_equation_param_t  *eqp = cs_equation_param_by_name("velocity");

    cs_equation_add_bc_by_analytic(eqp,
                                   CS_PARAM_BC_DIRICHLET,
                                   "inlet",           // zone name
                                   _vel_profile,      // pointer to the function
                                   NULL);             // input structure
  }

  {
    cs_field_t  *f = cs_field_by_name("scalar1");
    cs_equation_param_t  *eqp = cs_equation_param_by_name("scalar1");

    cs_equation_add_bc_by_analytic(eqp,
                                   CS_PARAM_BC_DIRICHLET,
                                   "inlet",           // zone name
                                   _scalar_profile,   // pointer to the function
                                   f);                // input structure
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
