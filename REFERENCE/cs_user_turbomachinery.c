/*============================================================================
 * Definition of turbomachinery related options.
 *
 * Turbomachinery-related user functions (called in this order):
 *   1) Define rotor cells and associated axis
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
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_turbomachinery.c
 *
 * \brief Definition of turbomachinery related options.
 *
 * See \ref turbomachinery for examples.
*/
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define rotor/stator model.
*/
/*----------------------------------------------------------------------------*/

void
cs_user_turbomachinery(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define rotor axes, associated cells, and rotor/stator faces.
*/
/*----------------------------------------------------------------------------*/

void
cs_user_turbomachinery_rotor(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define rotation velocity of rotor.
*/
/*----------------------------------------------------------------------------*/

void
cs_user_turbomachinery_set_rotation_velocity(void)
{

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
