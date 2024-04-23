/*============================================================================
 * User initialization prior to solving time steps.
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

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_initialization-base.c
 *
 * \brief Initialization prior to solving time steps.
 *        Basic examples
 *
 * See \ref cs_user_initialization for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize variables.
 *
 * This function is called at beginning of the computation
 * (restart or not) before the time step loop.
 *
 * This is intended to initialize or modify (when restarted)
 * variable and time step values.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_initialization(cs_domain_t     *domain)
{
  const cs_mesh_t *m = domain->mesh;
  cs_mesh_quantities_t *mq =domain->mesh_quantities;
  const cs_real_3_t *restrict cell_cen
	  =(const cs_real_3_t *restrict)mq->cell_cen;


  /* If this is restarted computation, do not reinitialize values */
  if (domain->time_step->nt_prev > 0)
    return;



  /*initialization k by channel flow */
 /* int i;
  int ny  = 101;
  int ny2 = 101;
  int ny3 = 101;
  int ix;
  double dmin;
  double dmin2;
  double dmin3;
  double xkent;
  double xeent;
  double d2s3 =2.0/3.0;
  int ituser[50000];
  double rtuser[50000];
  int ntcabs = domain->time_step->nt_cur;
  int ntpabs = domain->time_step->nt_prev;	  

  
  FILE *myFile;
  myFile = fopen ("U_mean.txt","r");

  if (myFile == NULL)
	  exit(1);

  for (i=0; i<ny;i++)
  { 
        if( fscanf(myFile, "%lf %lf", &rtuser[i], &rtuser[i+ny])!=2)
		exit(1);
  }
  fclose (myFile); */
/*
  FILE *myFile1;
  myFile1 = fopen ("k.txt","r");
  for (i=0; i<ny2;i++)
  { 
	  fscanf(myFile1, "%lf %lf", &rtuser[i+3*ny2], &rtuser[i+4*ny2]);

  }
  fclose (myFile1);


  FILE *myFile2;
  myFile2 = fopen ("epsilon.txt","r");
  for (i=0; i<ny3;i++)
  { 
	  fscanf(myFile2, "%lf %lf", &rtuser[i+6*ny3], &rtuser[i+7*ny3]);

  }
  fclose (myFile2);

  */
/* initialization for the index array */ 
/*  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++)
  {
  ituser[cell_id]=0;
  ituser[cell_id+5*ny2]=0;
  ituser[cell_id+8*ny3]=0;
  }
*/


/*---------- for reading mean velocity ----------------------*/
/*
  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++)
   {    
	  dmin  = 1000.0;
	  dmin2 = 1000.0;
	  dmin3 = 1000.0; 
          for(ix = 0; ix<ny; ix++)
	  {
	  if (abs(cell_cen[cell_id][1]-rtuser[ny+ix])<dmin)
	  {
	  dmin = abs(cell_cen[cell_id][1]-rtuser[ny+ix]);
	  ituser[cell_id] = ix;
          printf ("the value of cell_cen at cell id %d is %lf\n",cell_id,cell_cen[cell_id][1]); 

	  }
	  printf ("the value of rtuser at ix %d is %lf\n ", ix, rtuser[ny+ix]);

	  printf ("the value of ituser at ix %d is %d\n",ix, ituser[cell_id]); 
  
	  }
   }

*/


  /*initialize "velocity_x" field by channel flow */
  cs_field_t *f1 = cs_field_by_name_try("velocity");
  if (f1 != NULL) {
    for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++)
      f1->val[cell_id*3] = 0.0;
  }

  /*initialize "k" field by channel flow */
  cs_field_t *f2 = cs_field_by_name_try("k");

  if (f2 != NULL) {
    for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++)
      f2->val[cell_id] = 0.2;
  }

  
  /*initialize "epsilon" field by channel flow */
  cs_field_t *f3 = cs_field_by_name_try("epsilon");

  if (f3 != NULL) {
    for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++)
      f3->val[cell_id] = 0.2;
  }

  /* Initialize "RPV" field to 0.0 only if it exists  */
  cs_field_t *f4 = cs_field_by_name_try("RPV");

  if (f4 != NULL) {
    for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++)
      f4->val[cell_id] = 0.0;
  }

  /* Initialize "Theta" field to 0.0 only if it exists  */
  cs_field_t *f5 = cs_field_by_name_try("Theta");

  if (f5 != NULL) {
    for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++)
      f5->val[cell_id] = 0.0;
  }
  /* Initialize "RPV_variance" field to 0.0 only if it exists  */
  cs_field_t *f6 = cs_field_by_name_try("RPV_variance");

  if (f6 != NULL) {
    for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++)
      f6->val[cell_id] = f4->val[cell_id]*(1.0-f4->val[cell_id]);
  }

  /* Initialize "Theta_variance" field to 0.0 only if it exists  */
  cs_field_t *f7 = cs_field_by_name_try("Theta_variance");

  if (f7 != NULL) {
    for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++)
      f7->val[cell_id] = f5->val[cell_id]*(1.0-f5->val[cell_id]);
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
