/*============================================================================
 * Define postprocessing output.
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

#include "stdlib.h"
#include "string.h"

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local (user defined) function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Example function for advanced selection of interior faces.
 *
 * Selects interior faces separating cells of group "2" from those
 * of group "3", (assuming no cell has both colors).
 *
 * parameters:
 *   input    <-> pointer to input (unused here)
 *   n_faces  --> number of selected faces
 *   face_ids --> array of selected face ids (0 to n-1 numbering)
 *----------------------------------------------------------------------------*/

static void
_i_faces_select_example(void         *input,
                        cs_lnum_t    *n_faces,
                        cs_lnum_t   **face_ids)
{
  CS_UNUSED(input);

  cs_lnum_t i, face_id;
  int n_families = 0;
  int *family_list = NULL;
  int *family_mask = NULL;

  cs_lnum_t n_i_faces = 0;
  cs_lnum_t *i_face_ids = NULL;

  const cs_mesh_t *m = cs_glob_mesh;

  /* Allocate selection list */

  BFT_MALLOC(i_face_ids, m->n_i_faces, cs_lnum_t);

  /* Build mask on families matching groups "2" (1), "3" (2) */

  BFT_MALLOC(family_list, m->n_families, int);
  BFT_MALLOC(family_mask, m->n_families, int);

  for (i = 0; i < m->n_families; i++)
    family_mask[i] = 0;

  cs_selector_get_family_list("2",  &n_families, family_list);

  for (i = 0; i < n_families; i++)
    family_mask[family_list[i] - 1] += 1;

  cs_selector_get_family_list("3",  &n_families, family_list);

  for (i = 0; i < n_families; i++)
    family_mask[family_list[i] - 1] += 2;

  BFT_FREE(family_list);

  /* Now that mask is built, test for adjacency */

  for (face_id = 0; face_id < m->n_i_faces; face_id++) {

    /* Adjacent cells  and flags */

    cs_lnum_t c1 = m->i_face_cells[face_id][0];
    cs_lnum_t c2 = m->i_face_cells[face_id][1];

    int iflag1 = family_mask[m->cell_family[c1]];
    int iflag2 = family_mask[m->cell_family[c2]];

    /* Should the face belong to the extracted mesh ? */

    if ((iflag1 == 1 && iflag2 == 2) || (iflag1 == 2 && iflag2 == 1)) {
      i_face_ids[n_i_faces] = face_id;
      n_i_faces += 1;
    }

  }

  /* Free memory */

  BFT_FREE(family_mask);
  BFT_REALLOC(i_face_ids, n_i_faces, cs_lnum_t);

  /* Set return values */

  *n_faces = n_i_faces;
  *face_ids = i_face_ids;
}

/*----------------------------------------------------------------------------
 * Example function for selection of boundary faces.
 *
 * selects boundary faces of group "4".
 *
 * parameters:
 *   input    <-> pointer to input (unused here)
 *   n_faces  --> number of selected faces
 *   face_ids --> array of selected face ids (0 to n-1 numbering)
 *----------------------------------------------------------------------------*/

static void
_b_faces_select_example(void         *input,
                        cs_lnum_t    *n_faces,
                        cs_lnum_t   **face_ids)
{
  CS_UNUSED(input);

  cs_lnum_t n_b_faces = 0;
  cs_lnum_t *b_face_ids = NULL;

  const cs_mesh_t *m = cs_glob_mesh;

  /* Allocate selection list */

  BFT_MALLOC(b_face_ids, m->n_b_faces, cs_lnum_t);

  /* Use simple selection function */

  cs_selector_get_b_face_list("4", &n_b_faces, b_face_ids);

  /* Adjust array to final size (cleaner, but not required) */

  BFT_REALLOC(b_face_ids, n_b_faces, cs_lnum_t);

  /* Set return values */

  *n_faces = n_b_faces;
  *face_ids = b_face_ids;
}

/*----------------------------------------------------------------------------
 * Example function for selection of cells with scalar field values above
 * a certain threshold.
 *
 * In this example, the selection is base on the value of a scalar field
 * named "he_fraction" being above above 0.05.
 *
 * parameters:
 *   input    <-> pointer to input (unused here)
 *   n_cells  --> number of selected cells
 *   cell_ids --> array of selected cell ids (0 to n-1 numbering)
 *----------------------------------------------------------------------------*/

static void
_he_fraction_05_select(void        *input,
                       cs_lnum_t   *n_cells,
                       cs_lnum_t  **cell_ids)
{
  CS_UNUSED(input);

  cs_lnum_t _n_cells = 0;
  cs_lnum_t *_cell_ids = NULL;

  const cs_mesh_t *m = cs_glob_mesh;

  cs_field_t *f = cs_field_by_name_try("He_fraction"); /* Get access to field */

  if (f == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "No field with name \"He_fraction\" defined");

  /* Before time loop, field is defined, but has no values yet,
     so ignore that case (postprocessing mesh will be initially empty) */

  if (f->val != NULL) {

    BFT_MALLOC(_cell_ids, m->n_cells, cs_lnum_t); /* Allocate selection list */

    for (cs_lnum_t i = 0; i < m->n_cells; i++) {
      if (f->val[i] > 5.e-2) {
        _cell_ids[_n_cells] = i;
        _n_cells += 1;
      }
    }

    BFT_REALLOC(_cell_ids, _n_cells, cs_lnum_t); /* Adjust size (good practice,
                                                    but not required) */

  }

  /* Set return values */

  *n_cells = _n_cells;
  *cell_ids = _cell_ids;
}

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define post-processing writers.
 *
 * The default output format and frequency may be configured, and additional
 * post-processing writers allowing outputs in different formats or with
 * different format options and output frequency than the main writer may
 * be defined.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_writers(void)
{
  /* Set time plot file writer flush behavior defaults. */

  cs_time_plot_set_flush_default(1800, /* flush_wtime */
                                 -1);  /* n_buffer_steps */

  /* Default output format and options */

  /* Redefine default writer */
  /* ----------------------- */

  cs_post_define_writer(CS_POST_WRITER_DEFAULT,       /* writer_id */
                        "results",                    /* writer name */
                        "postprocessing",             /* directory name */
                        "EnSight Gold",               /* format_name */
                        "",                           /* format_options */
                        FVM_WRITER_FIXED_MESH,
                        true,                         /* output_at_start */
                        true,                         /* output_at_end */
                        10,                         /* frequency_n */
                        -1.0);                        /* frequency_t */

  /* Define additional writers */
  /* ------------------------- */
  cs_post_define_writer(1,       /* writer_id */
                        "BC_values",                    /* writer name */
                        "postprocessing",             /* directory name */
                        "EnSight Gold",               /* format_name */
                        "",                           /* format_options */
                        FVM_WRITER_FIXED_MESH,
                        true,                         /* output_at_start */
                        true,                         /* output_at_end */
                        10,                         /* frequency_n */
                        -1.0);                        /* frequency_t */

  /* Common parameters for all writers */

  double frequency_n = -1.0;
  double frequency_t = -1.0;

}


/*----------------------------------------------------------------------------*/
/*!
 * \brief Define post-processing meshes.
 *
 * The main post-processing meshes may be configured, and additional
 * post-processing meshes may be defined as a subset of the main mesh's
 * cells or faces (both interior and boundary).
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_meshes(void)
{
  /* Reconfigure predefined meshes (mesh_id -1 for volume, -2 for boundary */

  /* De-activate boundary mesh output by redefining it with no writer
     association (default is:
     int n_writers = 1;
     const int writer_ids[] = {CS_POST_WRITER_DEFAULT});
  */

  {
    int n_writers = 1;
    const int *writer_ids[] = {1};

    cs_post_define_surface_mesh(CS_POST_MESH_BOUNDARY,  /* mesh_id */
                                "boundary",  /* mesh name */
                                NULL,        /* interior face selection criteria */
                                "all[]",     /* boundary face selection criteria */
                                true,        /* add_groups */
                                true,        /* automatic variables output */
                                n_writers,
                                writer_ids);
  }
}







/*----------------------------------------------------------------------------*/
/*!
 * Override default frequency or calculation end based output.
 *
 * This allows fine-grained control of activation or deactivation,
 *
 * \param[in]  nt_max_abs  maximum time step number
 * \param[in]  nt_cur_abs  current time step number
 * \param[in]  t_cur_abs   absolute time at the current time step
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_activate(int     nt_max_abs,
                             int     nt_cur_abs,
                             double  t_cur_abs)
{
  CS_NO_WARN_IF_UNUSED(nt_cur_abs);
  CS_NO_WARN_IF_UNUSED(t_cur_abs);

  /* Use the cs_post_activate_writer() function to force the
   * "active" or "inactive" flag for a specific writer or for all
   * writers for the current time step.

   * the parameters for cs_post_activate_writer() are:
   *   writer_id <-- writer id, or 0 for all writers
   *   activate  <-- false to deactivate, true to activate */

  /* Example: deactivate all output before time step 1000 */

  if (nt_max_abs < 1000) {
    int writer_id = 0; /* 0: all writers */
    cs_post_activate_writer(writer_id, false);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
