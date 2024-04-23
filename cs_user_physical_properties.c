/*============================================================================
 * User definition of physical properties.
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
 * \file cs_user_physical_properties.c
 *
 * \brief User definition of physical properties.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function called at each time step to define physical properties.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_physical_properties(cs_domain_t   *domain)
{
  CS_NO_WARN_IF_UNUSED(domain);
  cs_real_t tau_hrp = 2.3; /* heat release parameter*/
  cs_fluid_properties_t *fp =cs_get_glob_fluid_properties();

  /* Check fields exists */
/*  if (CS_F_(lambda) == NULL)
    bft_error(__FILE__, __LINE__, 0,_("error lambda not variable\n"));*/
  if (CS_F_(rho) == NULL)
    bft_error(__FILE__, __LINE__, 0,_("error rho not variable\n"));
/*  if (CS_F_(cp) == NULL)
    bft_error(__FILE__, __LINE__, 0,_("error cp not variable\n"));*/

  cs_real_t *cpro_rho = CS_F_(rho)->val;

  cs_field_t *f1 = cs_field_by_name_try("RPV");


  /* Impose BML to calculate the variable density */
  {
      const cs_mesh_t *z = domain->mesh;	  

    for (cs_lnum_t cell_id = 0; cell_id < z->n_cells; cell_id++) {
      cpro_rho[cell_id] = fp->ro0/(1.0+tau_hrp*f1->val[cell_id]); /* fp-ro0 unburned density*/
    }
  }
}
/*----------------------------------------------------------------------------*/

END_C_DECLS
