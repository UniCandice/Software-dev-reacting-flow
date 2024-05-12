/*============================================================================
 * Examples for additional turbulence source terms for variable equations.
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


/*
#define VECTOR_DIM 3

// Function to calculate the norm of a vector
double vector_norm(double v[VECTOR_DIM]) {
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

// Function to normalize a vector
void normalize_vector(double v[VECTOR_DIM]) {
    double norm = vector_norm(v);
    
    if (norm == 0) {
        printf("Zero vector cannot be normalized.\n");
        return;
    }

    for (int i = 0; i < VECTOR_DIM; i++) {
        v[i] /= norm;
    }
}

// Function to normalize all vectors in a matrix
void normalize_all_vectors(double **matrix, int cell_number) {
    for (int i = 0; i < cell_number; i++) {
        normalize_vector(matrix[i]);
    }
}
*/


/*!
 * \file cs_user_source_terms-turbulence.c
 *
 * \brief Examples for additional turbulence source terms for
 *   variable equations.
 *
 * See the reference \ref cs_user_source_terms.c for documentation.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function called at each time step to define source terms.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in]       f_id     field id of the variable
 * \param[out]      st_exp   explicit source term
 * \param[out]      st_imp   implicit part of the source term
 */
/*----------------------------------------------------------------------------*/

void
cs_user_source_terms(cs_domain_t  *domain,
                     int           f_id,
                     cs_real_t    *st_exp,
                     cs_real_t    *st_imp)
{
  CS_NO_WARN_IF_UNUSED(domain);

  /* mesh quantities */
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_lnum_t  n_cells = cs_glob_mesh->n_cells;
  const cs_real_t  *cell_f_vol = cs_glob_mesh_quantities->cell_vol;

  /* Define pointer f to the current turbulent variable field */
  const cs_field_t  *f = cs_field_by_id(f_id);

  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);
  double C_eps1=1.44;
  double C_eps2=1.44;
 
  /* Define a k_st pointer to the source term to k equation 
   * the first is the pressure work and the second is the pressure dilitation term*/
  cs_real_2_t  *k_st; 
  /* Define a epsilon_st pointer to the source term to epsilon equation 
   * the first is the pressure work and the second is the pressure dilitaion term*/
  cs_real_2_t  *eps_st;



  /* Define a cpro_rom pointer to the density */
  const cs_real_t   *cpro_rom       = CS_F_(rho)->val;
   /* Define a cpro_visct pointer to the turbulent dynamic viscosity */
  const  cs_real_t  *cpro_visct     = CS_F_(mu_t)->val;
  /* Define a cvar_k pointer to the turbulent kenetic energy */
  const  cs_real_t  *cvar_k         = CS_F_(k)->val_pre;
  /* Define a cvar_eps pointer to the turbulent dissipation */
  const  cs_real_t  *cvar_eps       = CS_F_(eps)->val_pre;
  /* Define a cvar_rpv pointer to the reaction progress variable */
  const  cs_field_t *cvar_RPV       = cs_field_by_name_try("RPV");
 /* Define a cvar_rpv pointer to the reaction progress variable */
  const  cs_field_t *cvar_RPV_bar   = cs_field_by_name_try("RPV_bar");
  /* Define a cvar_rpv pointer to the dimensionless temperature */
  const  cs_field_t *cvar_theta       = cs_field_by_name_try("Theta");
  /* Define a cvar_omegac pointer to the source term of the reaction progress varaible */
  const  cs_field_t *cvar_omegac    = cs_field_by_name_try("source_term");
  /* Define a cvar_sflux_RPV pointer to the scalar flux of RPV */
  cs_field_t *cvar_sflux_RPV = cs_field_by_name_try("sflux_RPV");
  /* Define a cvar_sflux_RPV pointer to the scalar flux of dimensionless temperature */
  cs_field_t *cvar_sflux_theta = cs_field_by_name_try("sflux_Theta");
  /* Define a cvar_pr pointer to the pressure*/
  const  cs_field_t *cvar_pr        = CS_F_(p);

  /* gradc: reaction progress variable (RPV) gradient  
   * gradp: pressure gradient 
   * the data is stored ad gradc[n_cells][dim] as cs_real_3_t 
   * cvar_sflux_RPV: scalar flux of RPV stored in properties added by cs_user_parameters.c and the format is 
   * x direction: val[3*cell_id]
   * y direction: val[3*cell_id+1]
   * z direction: val[3*cell_id+2]*/  
  cs_real_3_t  *gradc,*gradp,*gradt;

  bool use_previous_t = true;

  BFT_MALLOC(gradc,m->n_cells_with_ghosts,cs_real_3_t);
     cs_field_gradient_scalar(cvar_RPV,
                              use_previous_t,
                              1,       /* inc */
                              true,    /* iccocg */ 
                              gradc);

   for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      cvar_sflux_RPV->val[3*cell_id]   = -cpro_visct[cell_id]* gradc[cell_id][0]/1.0; /*TODO turbulent Prandtl number Prt=1.0 */      
      cvar_sflux_RPV->val[3*cell_id+1] = -cpro_visct[cell_id]* gradc[cell_id][1]/1.0;        
      cvar_sflux_RPV->val[3*cell_id+2] = -cpro_visct[cell_id]* gradc[cell_id][2]/1.0;    
   }

   BFT_MALLOC(gradt,m->n_cells_with_ghosts,cs_real_3_t);
     cs_field_gradient_scalar(cvar_theta,
                              use_previous_t,
                              1,       /* inc */
                              true,    /* iccocg */ 
                              gradt);

   for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      cvar_sflux_theta->val[3*cell_id]   = -cpro_visct[cell_id]* gradt[cell_id][0]/1.0; /*TODO turbulent Prandtl number Prt=1.0 */      
      cvar_sflux_theta->val[3*cell_id+1] = -cpro_visct[cell_id]* gradt[cell_id][1]/1.0;        
      cvar_sflux_theta->val[3*cell_id+2] = -cpro_visct[cell_id]* gradt[cell_id][2]/1.0;    
   }

  BFT_MALLOC(gradp,m->n_cells_with_ghosts,cs_real_3_t);

     cs_field_gradient_scalar(cvar_pr,
                              use_previous_t,
                              1,       /* inc */
                              true,    /* iccocg */ 
                              gradp);

  BFT_MALLOC(k_st,m->n_cells_with_ghosts,cs_real_2_t);
  BFT_MALLOC(eps_st,m->n_cells_with_ghosts,cs_real_2_t);
  
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) { 
  k_st[cell_id][0]=-cvar_sflux_RPV->val[3*cell_id]*gradp[cell_id][0]-cvar_sflux_RPV->val[3*cell_id+1]*gradp[cell_id][1]-cvar_sflux_RPV->val[3*cell_id+2]*gradp[cell_id][2];
  k_st[cell_id][1]= 0.5*pow(2.3*0.7,2)*cvar_RPV->val_pre[cell_id]*cvar_omegac->val[cell_id]; /* TODO 2.3 heat release parameter 0.7 laminar flame speed */
  eps_st[cell_id][0]=  C_eps1*cvar_eps[cell_id]*k_st[cell_id][0]/cvar_k[cell_id]; 
  eps_st[cell_id][1]=  C_eps2*cvar_eps[cell_id]*k_st[cell_id][1]/cvar_k[cell_id]; 
  }

  BFT_FREE(gradc);
  BFT_FREE(gradp);
  BFT_FREE(gradt);


/*
   cs_field_t *cvar_Debug = cs_field_by_name_try("Debug");
   for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      cvar_Debug->val[cell_id] = k_st[cell_id][0]+k_st[cell_id][1];    
   }
*/

  /*  source term for turbulence models 
   * (Source term on the TKE 'k' here)
   *
   *  Source term for cvar_var:
   *  rho cell_f_vol d(cvar_var)/dt       = ...
   *                 ... - rho*cell_f_vol*ff - rho*cell_f_vol*cvar_var/tau
   *
   *  (Source term on the epsilon equation 'eps' here )
   */

  if (f == CS_F_(k)) {

    if (eqp->verbosity >= 1)
      bft_printf(" User source terms for turbulence variable %s\n",
                 cs_field_get_label(f));

    for (cs_lnum_t i = 0; i < n_cells; i++) {

      /* explicit source terms */
      st_exp[i] = (k_st[i][0]+k_st[i][1])*cell_f_vol[i];

      /* implicit source terms */
      st_imp[i] = 0.0;
    }
  }

  if (f == CS_F_(eps)) {

    if (eqp->verbosity >= 1)
      bft_printf(" User source terms for turbulence variable %s\n",
                 cs_field_get_label(f));


    for (cs_lnum_t i = 0; i < n_cells; i++) {

      /* explicit source terms */
      st_exp[i] = (eps_st[i][0]+eps_st[i][1])*cell_f_vol[i];

      /* implicit source terms */
      st_imp[i] = 0.0;
    }
  }

  BFT_FREE(k_st);
  BFT_FREE(eps_st);


  /*scalar flux for countergrad  */ 

  int ico_countergrad_c=1;

  if (ico_countergrad_c==1){

   /* Define a cvar_rpv pointer to the reaction progress variable */
  const  cs_field_t *countergrad_sum       = cs_field_by_name_try("countergrad_sum");
  /* Define a cvar_rpv pointer to the dimensionless temperature */
  const  cs_field_t *sflux_countergrad_scalar       = cs_field_by_name_try("sflux_countergrad_scalar");


  cs_real_3_t  *grad1,*grad2,*sflux_countergrad;

  bool use_previous_t = true;

  BFT_MALLOC(grad1,m->n_cells_with_ghosts,cs_real_3_t);
  BFT_MALLOC(grad2,m->n_cells_with_ghosts,cs_real_3_t);
  BFT_MALLOC(sflux_countergrad,m->n_cells_with_ghosts,cs_real_3_t);

     cs_field_gradient_scalar(cvar_RPV,
                              use_previous_t,
                              1,       
                              true,    
                              grad1);

//  normalize_all_vectors(gradc,m->n_cells_with_ghosts);

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
	    double norm = 0;
	    double tau_hrp= 2.3;
	    double ctilde = cvar_RPV->val_pre[cell_id];
	    double cbar = (1.0+tau_hrp)*ctilde/(1.0+tau_hrp*ctilde);
	    norm=pow(grad1[cell_id][0],2)+ pow (grad1[cell_id][1],2)+ pow(grad1[cell_id][2],2)+1e-8;
            norm=sqrt(norm);

            sflux_countergrad[cell_id][0]  = - grad1[cell_id][0]/norm*1.0*0.7*(ctilde-cbar);      
            sflux_countergrad[cell_id][1]  = - grad1[cell_id][1]/norm*1.0*0.7*(ctilde-cbar);        
            sflux_countergrad[cell_id][2]  = - grad1[cell_id][2]/norm*1.0*0.7*(ctilde-cbar);    
   }
  
   for (cs_lnum_t i = 0; i< 3; i++ ) {
	   for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
	   {
	   
	   sflux_countergrad_scalar->val[cell_id]= sflux_countergrad[cell_id][i];
	   
	   }
           cs_field_gradient_scalar(sflux_countergrad_scalar,
                              0,
                              1,       
                              true,    
                              grad2);

	   for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
	   {
	   
	   countergrad_sum->val[cell_id] -= grad2[cell_id][i] ;// sum the values in the diagonal 
	   
	   }

   }
// filtering out the gradient close to 0  
   for (cs_lnum_t cell_id=0; cell_id<n_cells;cell_id++)
   {
	   if (cvar_RPV->val[cell_id]<=0.1||cvar_RPV->val[cell_id]>=0.9)
	   {
		   countergrad_sum->val[cell_id] = 0.0;
	   
	   }
   
   }

// try other way 

//  normalize_all_vectors(gradc,m->n_cells_with_ghosts);
     cs_field_gradient_scalar(cvar_RPV,
                              use_previous_t,
                              1,       
                              true,    
                              grad1);

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
	    double norm = 0;
	    double tau_hrp= 2.3;
	    double ctilde = cvar_RPV->val_pre[cell_id];
            cvar_RPV_bar->val[cell_id] = (1.0+tau_hrp)*ctilde/(1.0+tau_hrp*ctilde);

	    norm=pow(grad1[cell_id][0],2)+ pow (grad1[cell_id][1],2)+ pow(grad1[cell_id][2],2)+1e-8;
            norm=sqrt(norm);

            sflux_countergrad[cell_id][0]  = - grad1[cell_id][0]/norm*1.0*0.7;      
            sflux_countergrad[cell_id][1]  = - grad1[cell_id][1]/norm*1.0*0.7;        
            sflux_countergrad[cell_id][2]  = - grad1[cell_id][2]/norm*1.0*0.7;    
   }

// filtering out the gradient close to 0  
   for (cs_lnum_t cell_id=0; cell_id<n_cells;cell_id++)
   {
	   if (cvar_RPV->val[cell_id]<=0.1||cvar_RPV->val[cell_id]>=0.9)
	   {
		   sflux_countergrad[cell_id][0] = 0.0;
	           sflux_countergrad[cell_id][1] = 0.0;
   		   sflux_countergrad[cell_id][2] = 0.0;

	   }
   
   }

 

      cs_field_gradient_scalar(cvar_RPV,
                              use_previous_t,
                              1,       
                              true,    
                              grad1);

   
      cs_field_gradient_scalar(cvar_RPV_bar,
                              0,
                              1,       
                              true,    
                              grad2);

    

	   for (cs_lnum_t i =0; i<1; i++){

           for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
	   {
	   
	   countergrad_sum->val[cell_id]=-sflux_countergrad[cell_id][0]*(grad1[cell_id][0]- grad2[cell_id][0])  
		                         -sflux_countergrad[cell_id][1]*(grad1[cell_id][1]- grad2[cell_id][1])
					 -sflux_countergrad[cell_id][2]*(grad1[cell_id][2]- grad2[cell_id][2]); // sum the values in the diagonal 
	   
	   } 
	   }

     

 
  BFT_FREE(sflux_countergrad);
  BFT_FREE(grad1);
  BFT_FREE(grad2);

  }




}

/*----------------------------------------------------------------------------*/

END_C_DECLS
