/*============================================================================
 * This function is called at the end of each time step, and has a very
 *  general purpose
 *  (i.e. anything that does not have another dedicated user function)
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

static void 

_geom_closest_point(cs_lnum_t         n_points,
                      const cs_real_t   point_coords[][3],
                      const cs_real_t   query_coords[3],
                      cs_lnum_t        *point_id,
                      int              *rank_id)
{
  cs_lnum_t id_min = -1;
  cs_real_t d2_min = HUGE_VAL;

  for (cs_lnum_t i = 0; i < n_points; i++) {
    cs_real_t d2 = cs_math_3_square_distance(point_coords[i], query_coords);
    if (d2 <= d2_min) {
      d2_min = d2;
      id_min = i;
    }
  }



    *point_id = id_min;
}







/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_extra_operations.c
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is called at the end of each time step.
 *
 * It has a very general purpose, although it is recommended to handle
 * mainly postprocessing or data-extraction type operations.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations(cs_domain_t     *domain)
{
  /* Calculation of the difference between the c and T */


  const cs_lnum_t n_cells     = domain->mesh->n_cells;
  const cs_lnum_t n_b_faces    = domain->mesh->n_b_faces;

  const cs_real_3_t *cell_cen
    = (const cs_real_3_t *)domain->mesh_quantities->cell_cen;
  const cs_real_3_t *cdgfbo
    = (const cs_real_3_t *)domain->mesh_quantities->b_face_cog;


  const cs_real_t    *f_b_RPV
    = (const cs_real_t *)cs_field_by_name("boundary_RPV")->val;
 
  const cs_real_t    *f_b_theta
    = (const cs_real_t *)cs_field_by_name("boundary_Theta")->val;

  cs_real_t    *diff_c_T =(cs_real_t *) cs_field_by_name("Difference_c_T")->val;


  cs_lnum_t nlelt;
  cs_lnum_t *lstelt;
  BFT_MALLOC(lstelt, n_b_faces, cs_lnum_t);

  cs_selector_get_b_face_list
    ("TOP_WALL or BOTTOM_WALL",
      &nlelt, lstelt);

    cs_real_3_t *xyz_wall = NULL;
/*    cs_real_3_t *xyz_cell = NULL;*/


    cs_real_t   *wall_RPV = NULL;
    cs_real_t   *wall_theta= NULL;

    BFT_MALLOC(xyz_wall, nlelt, cs_real_3_t);

    BFT_MALLOC(wall_RPV, nlelt, cs_real_t);
    BFT_MALLOC(wall_theta, nlelt, cs_real_t);
    
    for (cs_lnum_t ilelt = 0; ilelt < nlelt; ilelt ++) {
      cs_lnum_t face_id = lstelt[ilelt];
      xyz_wall[ilelt][0] = cdgfbo[face_id][0];
      xyz_wall[ilelt][1] = cdgfbo[face_id][1];
      xyz_wall[ilelt][2] = cdgfbo[face_id][2];

      wall_RPV[ilelt]=f_b_RPV[face_id];
      wall_theta[ilelt]=f_b_theta[face_id];
    }

/*
    BFT_MALLOC(xyz_cell, n_cells, cs_real_3_t);


    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id ++) {


      cs_real_t xyz_cell[cell_id][0] = cell_cen[cell_id][0],
      cs_real_t xyz_cell[cell_id][1] = cell_cen[cell_id][1],
      cs_real_t xyz_cell[cell_id][2] = cell_cen[cell_id][2];

    }
*/
/*
    cs_parall_allgather_ordered_r(nlelt, neltg, 1, xabs, xabs, xabsg);

    BFT_FREE(treloc);
    BFT_FREE(lstelt);
*/
   
    cs_lnum_t f_id_prev = -1;
    int irang1 = -1;
    int irangv;
  
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id ++) {

      cs_lnum_t f_id;


      cs_real_t xyz[3] = {cell_cen[cell_id][0],
                          cell_cen[cell_id][1],
                          cell_cen[cell_id][2]};

     _geom_closest_point(nlelt,
                            xyz_wall,
                            xyz,
                            &f_id,
			    &irangv);


            diff_c_T[cell_id] = wall_RPV[f_id] - wall_theta[f_id];

/*
      if (f_id != f_id_prev) {
          f_id_prev = f_id;
            diff_c_T[cell_id] = wall_RPV[f_id] - wall_theta[f_id];
        } */
       
      }

    BFT_FREE(lstelt);
    BFT_FREE(xyz_wall);
    BFT_FREE(wall_RPV);
    BFT_FREE(wall_theta);

}
/*----------------------------------------------------------------------------*/

END_C_DECLS
