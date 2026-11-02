/* Version 21-10-2004 JvH: added memalloc for parallel-related params

 * This file contains the following functions:
 * void init_tracers()          - main init routine for tracers
 * void read_trac_input()       - read tracer input
 * void trac_memalloc()         - arrange mem.alloc. for tracer arrays
 * void trac_distr()            - distribute tracers over domain
 * void trac_arrange_funcinit() - main func. for trac.val. assignment
 * float trac_funcinit()        - calc. of tracer function value
 * void tracexit()              - clean prog end if error occurs
 */
#include <signal.h>
#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include <string.h>
#include <stdlib.h>
#include "element_definitions.h"
#include "global_defs.h"

void init_tracers_preliminary(E)
   struct All_variables *E;
{  int i,noz,lev;
   int input_int();
   int input_float();
   int input_int_vector();
   float trac_funcinit();

    /* Read if tracers are used. Default=not=0 */
  input_int("tracers_used",&(E->tracers.tracers_used),"0");

  if (E->tracers.tracers_used > 0)
  {  /* Tracers used: read more info */
   /* In case of active tracers, 1 or more function values are
    * assigned to each tracer (e.g. composition, strainrate, 
    * melt depletion etc). First read the number of functions 
    */
   E->tracers.ntracfunc=0;
   input_int("ntracfunc",&(E->tracers.ntracfunc),"0");
   fprintf(E->fp,"ntracfunc = # tracer functions = %d\n", 
                                      E->tracers.ntracfunc);
   if (E->tracers.ntracfunc > NTRACFUNCMAX)
   { fprintf(stderr,"Too many tracer function: truncating to %d\n",
             NTRACFUNCMAX);
     E->tracers.ntracfunc=NTRACFUNCMAX;
   }
   if (E->tracers.ntracfunc > 0)
     {
      /* Allocate memory for array in which tracer function types are
       * stored.
       */
      E->tracers.tracfunc =(int *) calloc(E->tracers.ntracfunc+1,
                                                          sizeof(int));
      input_int_vector("tracfunc",E->tracers.ntracfunc,
                                            (E->tracers.tracfunc));
      E->tracers.comporheo=1.;
      input_float("comporheo",&(E->tracers.comporheo),"1.");
      fprintf(E->fp,"Compositional rheology number Delta eta_C=%g\n", 
                   E->tracers.comporheo);
     }
  }

  if ( E->tracers.tracers_used > 0 &&
        E->tracers.ntracfunc > 0 &&
        E->tracers.comporheo!=1.)
   {  E->tracers.compoviscoused = 1;
      fprintf(E->fp, " Composition dependent rheology used\n");
      lev=E->mesh.levmax;
      noz=E->lmesh.NOZ[lev];

      E->tracers.init_at_inflow0=(float *) calloc(noz+1,sizeof(float));
      E->tracers.init_at_inflow1=(float *) calloc(noz+1,sizeof(float));

      if (E->viscosity.RHEOL==111)
       for (i=1;i<=noz;i++)
       {
	 E->tracers.init_at_inflow0[i]=
	   trac_funcinit(E,E->tracers.tracfunc[0],0.,0.,E->XX[lev][2][i]);
         E->tracers.init_at_inflow1[i]=
	   trac_funcinit(E,E->tracers.tracfunc[1],0.,0.,E->XX[lev][2][i]);
       }
    }
   else
    { E->tracers.compoviscoused = 0;
      fprintf(E->fp, " No composition dependent rheology used\n");
    }
  return;
}


void init_tracers(E)
   struct All_variables *E;
{
 
   void read_trac_input();
   void trac_memalloc();
   void parallel_tracer_setup();
   void trac_distr();
   float trac_funcinit();
   void tracexit();
   void katz_tables();
   void find_trac_in_elem();
   void update_trac_tables();
   void trac_arrange_funcinit();
   void tracinterp_t2e();
   void tracinterp_e2t();
   int input_int();


   if (E->tracers.tracers_used > 0)
   {  /* Tracers used: read more info */
      read_trac_input(E);
      /* Set defaults */
      E->tracers.type = 1; /* normal tracers: move in each dir */
  
      /* Tracer arrays memory allocation */
      trac_memalloc(E);

      /* Tracer distribution */
      trac_distr(E);
     
      /* After tracers have been positioned, they should be found in 
         the elems to allow for interpolation */
      if (E->tracers.type<2)
      {  find_trac_in_elem(E, E->tracers.xtrac, E->tracers.elemnr, 
                  &(E->tracers.lntrac),1);
         update_trac_tables(E, E->tracers.elemnr, 
                  E->tracers.ntrac_in_elem, E->tracers.tracs_in_elem,
                  E->tracers.lntrac);
      }

      if (E->tracers.ntracfunc>0)
      {  /* Tracer function assignment if tracers are 'active' */
         trac_arrange_funcinit(E); 

         /* Interpolate initial FEM fields to tracers if
          * necessary */
         tracinterp_e2t(E);

         /* tracer function values are interpolated to the FEM to add 
            Ra_c*C to the Stokes force term */
         tracinterp_t2e(E);

	 if (E->data.heat_of_melt !=0.)
	   /* katz tables are loaded to allow prediction of consumed Hmelt */
	   katz_tables(E);

      }

      /* Setup parallel array memallocs */
      parallel_tracer_setup(E);

      E->tracers.initialized=1;
     
   } 
   else 
      fprintf(stderr, "no tracers used\n");

   return;
}

void read_trac_input(E)
   struct All_variables *E;
{
   int input_float();
   int input_int();
   int input_int_vector();
   int input_float_vector();
   void tracexit();
   int ifunc,extraction;

   /* Error/Warning treatment:
      0 = Just warn, don't stop at problems
      1 = Stop at 'empty' elements in case of 'active' tracers
          (interp_t2e=1 for one of the functions)
    */
   E->tracers.Warn=0;
   input_int("tracer_Warn",&(E->tracers.Warn),"0");

   E->tracers.restart=0;
   input_int("tracer_restart",&(E->tracers.restart),"0");

   E->tracers.storage=100;
   input_int("tracer_store_every",&(E->tracers.storage),"100");

   /* Give the number of tracers per unit LENGTH */
   E->tracers.tracer_density=1.;
   input_float("tracer_density",&(E->tracers.tracer_density),"1.");
   input_float("thermochemical_layerd",&(E->tracers.thermochemical_layerd),".0");
   E->tracers.trapmelt_total=0.;
   fprintf(E->fp,"tracer density:%g\n", E->tracers.tracer_density);

   /* Buffer fraction added on top of lntrac when allocating tracer arrays.
    * Increase this if the model crashes with "ntrac>arraylen" during parallel
    * tracer exchange (e.g. set to 1.0 for a 100% buffer). Default: 0.5 */
   E->tracers.array_buf_factor = 0.5;
   input_float("tracer_array_buf_factor",&(E->tracers.array_buf_factor),"0.5");
   fprintf(E->fp,"tracer array buf factor:%g\n", E->tracers.array_buf_factor);

   /* Give the kind of tracer distribution:
    * 0     = random
    * 1|2|3 = fixed testcase positions 
    * 10    = subduction tracers
    */
   E->tracers.distr=0;
   input_int("tracer_distr",&(E->tracers.distr),"0");
   fprintf(E->fp,"tracer distr:%d\n", E->tracers.distr);

   // ntracfunc already read in in init_tracers_preliminary
   if (E->tracers.ntracfunc>0)
   {
      /* Read-in whether tracer function is to be interpolated
       * to/from FEM nodal points; 1st calloc & set defaults */
      E->tracers.interp_t2e =(int *) calloc(E->tracers.ntracfunc+1,
                                                          sizeof(int));
      E->tracers.interp_e2t =(int *) calloc(E->tracers.ntracfunc+2,
                                                          sizeof(int));
      for (ifunc=0; ifunc<E->tracers.ntracfunc; ifunc++)
      {    E->tracers.interp_t2e[ifunc] = 1;
           E->tracers.interp_e2t[ifunc] = 1;
      }
      input_int_vector("interp_t2e",E->tracers.ntracfunc,
                                            (E->tracers.interp_t2e));
      input_int_vector("interp_e2t",E->tracers.ntracfunc,
                                            (E->tracers.interp_e2t));
      E->tracers.interp_e2t[E->tracers.ntracfunc+1]=0;

      /* BMO_density_mode == 1: FeO (tracfunc 80) and H2O (tracfunc 81) must be
       * interpolated tracer->element so that felem[] arrays are populated and
       * available in compofuncbuoy() for the Du et al. (2019) density formula.
       * Override any user-supplied interp_t2e=0 for these two tracers here,
       * BEFORE the felem allocation loop below.                               */
      if (E->control.BMO_density_mode == 1) {
          for (ifunc = 0; ifunc < E->tracers.ntracfunc; ifunc++) {
              if (E->tracers.tracfunc[ifunc] == 80 ||
                  E->tracers.tracfunc[ifunc] == 81) {
                  if (E->tracers.interp_t2e[ifunc] != 1) {
                      E->tracers.interp_t2e[ifunc] = 1;
                      fprintf(E->fp,
                          "[BMO_density] forcing interp_t2e=1 for tracfunc %d"
                          " (ifunc=%d) required by BMO_density_mode=1\n",
                          E->tracers.tracfunc[ifunc], ifunc);
                  }
              }
          }
      }

      for (ifunc=0; ifunc<E->tracers.ntracfunc; ifunc++)
      {  fprintf(E->fp, "tracer function type %d = %d",
                  ifunc, E->tracers.tracfunc[ifunc]);
         fprintf(E->fp, ",   interpolation trac->FEM (0=no,1=yes):%d\n",
                  E->tracers.interp_t2e[ifunc]);
         fprintf(E->fp, ",   interpolation FEM->trac (0=no,1=yes):%d\n",
                  E->tracers.interp_e2t[ifunc]);
      }

      E->tracers.ntracbuoyfunc=0;
      input_int("ntracbuoyfunc",&(E->tracers.ntracbuoyfunc),"0");
      fprintf(E->fp,"ntracbuoyfunc = # tracer buoyancy fns = %d\n", 
              E->tracers.ntracbuoyfunc); 
      if (E->tracers.ntracbuoyfunc > E->tracers.ntracfunc)
      {  fprintf(E->fp," ntracbuoyfunc > ntracfunc is not possible\n");
         fprintf(E->fp," reduced ntracbuoyfunc to ntracfunc = %d\n",
                      E->tracers.ntracfunc);
         E->tracers.ntracbuoyfunc = E->tracers.ntracfunc;
      }

      E->tracers.comporayl =(float *)
               calloc(E->tracers.ntracbuoyfunc+1, sizeof(float));
      input_float_vector("comporayl",E->tracers.ntracbuoyfunc,
                                            (E->tracers.comporayl));
      input_float("lambda",&(E->tracers.lambda),".0");
      E->tracers.lambda*=1.e6*E->data.layer_km*E->data.layer_km/(3.1557e7*E->data.therm_diff);
      fprintf(E->fp,"Q0compo decay constant lambda: %g\n", E->tracers.lambda);
      
      /* ---- initial is wrong ----   */
      E->tracers.Q0compo_tracnr = -1;           
      E->tracers.Q0compo_U_Th_tracnr = -1;
      E->tracers.Q0compo_value = 0.0f;

      input_int("Q0compo_tracnr",&(E->tracers.Q0compo_tracnr),"-1");
      input_int("Q0compo_U_Th_tracnr",&(E->tracers.Q0compo_U_Th_tracnr),"-1");

      /* Cache the ifunc index for tracfunc 83 (HPE solid profile).
       * Used by the BMO_hpe_mantle_on=1 solid-mantle heating term in
       * element_residual (Advection_diffusion.c).                        */
      E->tracers.hpe_mantle_tracnr = -1;
      {
          int jf_hpe;
          for (jf_hpe = 0; jf_hpe < E->tracers.ntracfunc; jf_hpe++)
              if (E->tracers.tracfunc[jf_hpe] == 83)
              {   E->tracers.hpe_mantle_tracnr = jf_hpe;  break; }
      }
      if (E->control.BMO_hpe_mantle_on) {
          if (E->tracers.hpe_mantle_tracnr < 0)
              fprintf(E->fp,
                  "[BMO_hpe_mantle] WARNING: BMO_hpe_mantle_on=1 but"
                  " tracfunc 83 not found; mantle HPE heating disabled.\n");
          else
              fprintf(E->fp,
                  "[BMO_hpe_mantle] tracfunc 83 found at ifunc=%d;"
                  " H_mantle=%.4g, mult=%.4g, tau_hpe=%.4g\n",
                  E->tracers.hpe_mantle_tracnr,
                  (double)E->control.BMO_H_mantle,
                  (double)E->control.BMO_H_mantle_mult,
                  (double)E->control.BMO_tau_hpe);
      }
      //input_float("Q0compo_threshold",&(E->tracers.Q0compo_threshold),"0");
      input_float("Q0compo_excess_value",&(E->tracers.Q0compo_value),"0");
      //E->tracers.Q0compo_value-=E->control.Q0;
      if (E->tracers.Q0compo_tracnr>=E->tracers.ntracfunc)
        tracexit (E,"Q0compo_trac >= ntracfunc: Exit");
      if (E->tracers.Q0compo_U_Th_tracnr>=E->tracers.ntracfunc)
        tracexit (E,"Q0compo_U-Th_trac >= ntracfunc: Exit");
      else if (E->tracers.Q0compo_U_Th_tracnr>=0)
        fprintf(E->fp," Radioactive heating U-Th switched on trac#=%d\n", E->tracers.Q0compo_U_Th_tracnr);
      input_int("inflownode",&(E->tracers.bot_inflow),"0");
      if (E->tracers.bot_inflow<=1) 
	E->tracers.bot_inflow=E->lmesh.NOZ[E->mesh.levmax];

      input_int("melttracers_used",&(E->tracers.melttracers_used),"0");
      input_float("melttracer_density",&(E->tracers.melttracer_density),"32");
      input_float("deltarhodepl",&(E->data.deltarhodepl),"0.0"); /*[kg/m^3]*/
      input_float("deltarhomelt",&(E->data.deltarhomelt),"0.0"); /*[kg/m^3]*/
      /* Calculate nondimensional Comporayl out of deltarho(depl/melt) */
      if (E->tracers.ntracbuoyfunc>0)   
	{if (E->tracers.interp_e2t[1] == 6)
      /* if Katz' melt model and melt retention  */
 	    {
               if (E->data.deltarhomelt != 0.0) 
 		E->tracers.comporayl[1] = E->data.deltarhomelt*E->control.Atemp/(E->data.density*E->data.therm_exp*E->data.ref_temperature);
               fprintf(E->fp,"deltarhomelt[kg/m^3]: %e\n",E->data.deltarhomelt); 
               /*fprintf(E->fp,"Calculated out of deltarho(depl/melt): "); */
	    }
      if (E->tracers.interp_e2t[0] == 14 || E->tracers.interp_e2t[0] == 22)
        {
                E->tracers.comporayl[0]=E->control.Atemp/(E->data.density*E->data.therm_exp*E->data.ref_temperature);
        }
      else if (E->tracers.interp_e2t[0] == 3 || E->tracers.interp_e2t[0] == 4 ||  E->tracers.interp_e2t[0] == 20)
	{
               if (E->data.deltarhodepl != 0.0) 
 		E->tracers.comporayl[0]=E->data.deltarhodepl*E->control.Atemp/(E->data.density*E->data.therm_exp*E->data.ref_temperature);
               fprintf(E->fp,"deltarhodepl[kg/m^3]: %e\n",E->data.deltarhodepl); 
               fprintf(E->fp,"Calculated out of deltarho(depl/melt): "); 
 	}}
      fprintf(E->fp,"Compositional Rayleigh numbers:\n"); 
      for (ifunc=0; ifunc<E->tracers.ntracbuoyfunc; ifunc++)
          fprintf(E->fp,"  Ra_C[%d]=%g\n", ifunc,
                                        E->tracers.comporayl[ifunc]);

      E->tracers.mindiffts=-1.;
      input_float("mindiffts",&(E->tracers.mindiffts),"0.");
      fprintf(E->fp,"Minimum diffusion timestep=%g\n", 
                   E->tracers.mindiffts);


      input_int("extraction",&(E->tracers.extraction),"11");
      extraction=E->tracers.extraction;
      if (extraction==1 || extraction==11)
	{ // extraction < 11: water fractionation
          E->tracers.phiR = 0.005;
	  E->tracers.phiC = 0.02;
	}
      else if (extraction==2 || extraction==12)
	{ // extraction < 11: water fractionation
          E->tracers.phiR = 0.002;
	  E->tracers.phiC = 0.02;
	}
      else if (extraction==3 || extraction==13)
	{ // extraction < 11: water fractionation
          E->tracers.phiR = 0.0;
	  E->tracers.phiC = 0.02;
	}
      else if (extraction==4 || extraction==14)
	{ // extraction < 11: water fractionation
          E->tracers.phiR = 0.01;
	  E->tracers.phiC = 0.02;
	}
      else if (extraction==20)
        E->tracers.phiR = E->tracers.phiC = 100.0;
      else if (extraction==21)
	{
          E->tracers.phiR = 0.;
	  E->tracers.phiC = 100.0;
	}

      if (E->tracers.phiR>E->tracers.phiC)
	tracexit(E, "phiR larger than phiC: STOP");
   }

   /* Test output 
    * 1) print tracer location */
   E->tracers.tractest = 0;
   input_int("tractest",&(E->tracers.tractest),"0");
   fprintf(E->fp, "Tracer test nr. %d\n", E->tracers.tractest);

   return;
}

void trac_memalloc(E)
   struct All_variables *E;
{
   int lev, nox, noy, noz, a, b, i, j, k, nodel, len, nelem, ifunc, n, n2, ielem;
   float dens, size, xmin, xmax, ymin, ymax, zmin, zmax, small;
   float lbuffer, first_upper_bound, upper_bound;
   float c[5];
   char filename[100];
   char input_s[200];
   FILE *fp;

   small=1e-4;
   lev=E->mesh.levmax;

   void restart_melttracers();

   dens = E->tracers.tracer_density;
   if (E->mesh.nsd == 2)
   {  size = E->mesh.layer[1] * E->XX[lev][2][E->tracers.bot_inflow];
      E->tracers.ntrac = size * dens * dens;
   }
   else if (E->mesh.nsd == 3)
   {
   size=E->mesh.layer[1]*E->XX[lev][2][E->tracers.bot_inflow]*E->mesh.layer[3];
   E->tracers.ntrac = size * dens * dens * dens;
   }
   else
   {  fprintf(stderr, "Error in init_tracers: unknown nsd=%d\n",
                      E->mesh.nsd);
      return;
   }
   fprintf(E->fp,"# tracers:%d\n", E->tracers.ntrac);

   /* local x-, z-, and y-ranges of processor me */
   nox=E->lmesh.NOX[lev];
   noy=E->lmesh.NOY[lev];
   noz=E->lmesh.NOZ[lev];
   if (E->mesh.nsd==3)
   {  i=1; j=1; k=1;
      nodel = j + (i-1)*noz + (k-1)*noz*nox;
      xmin = E->XX[lev][1][nodel];
      zmin = E->XX[lev][2][nodel];
      ymin = E->XX[lev][3][nodel];
      i=nox; j=1; k=1;
      nodel = j + (i-1)*noz + (k-1)*noz*nox;
      xmax = E->XX[lev][1][nodel];
      i=1; j=noz; k=1;
      nodel = j + (i-1)*noz + (k-1)*noz*nox;
      zmax = E->XX[lev][2][E->tracers.bot_inflow];
      i=1; j=1; k=noy;
      nodel = j + (i-1)*noz + (k-1)*noz*nox;
      ymax = E->XX[lev][3][nodel];
      E->tracers.lxmin=xmin;
      E->tracers.lymin=ymin;
      E->tracers.lzmin=zmin;
      E->tracers.lxmax=xmax;
      E->tracers.lymax=ymax;
      E->tracers.lzmax=zmax;
   }
   else 
   {  i=1; j=1; 
      nodel = j + (i-1)*noz;
      xmin = E->XX[lev][1][nodel];
      zmin = E->XX[lev][2][nodel];
      i=nox; j=1;
      nodel = j + (i-1)*noz;
      xmax = E->XX[lev][1][nodel];
      i=1; j=noz; 
      nodel = j + (i-1)*noz;
      zmax = E->XX[lev][2][E->tracers.bot_inflow];
      E->tracers.lxmin=xmin;
      E->tracers.lzmin=zmin;
      E->tracers.lxmax=xmax;
      E->tracers.lzmax=zmax;
      ymin=ymax=0;
      fprintf(E->fp,"xmin=%g,xmax=%g,zmin=%g,zmax=%g\n",
               xmin,xmax,zmin,zmax);
   }

   /* determine initial # tracers for processor me 
    * - add a small number small to avoid problems with conversion
    *   from float to int 
    */
   if (E->mesh.nsd == 2)
   {  size = (xmax-xmin) * (zmax-zmin);
      if (size < small)
      {  /* In case of postprocessing: no mesh is built: 
            temporary shortcut: take global dimensions */
         size = E->mesh.layer[1] * E->tracers.lzmax;
      }
      E->tracers.lntrac = (int) (size * dens * dens + small);
   }
   else if (E->mesh.nsd == 3)
   {  size = (xmax-xmin) * (zmax-zmin) * (ymax-ymin);
      if (size < small)
      {  /* In case of postprocessing: no mesh is built: 
            temporary shortcut: take global dimensions */
         size = E->mesh.layer[1] * E->tracers.lzmax * E->mesh.layer[3];
      }
      E->tracers.lntrac = (int) (size * dens * dens * dens + small);
   }
   fprintf(E->fp,"local # tracers:%d\n", E->tracers.lntrac);
   E->tracers.lntrac_usual=E->tracers.lntrac;
   ielem=E->tracers.bot_inflow-1;
   E->tracers.lntrac_bottom_usual=E->tracers.lntrac
     *E->ECO[lev][ielem].size[2]/E->tracers.lzmax;

   /* Allocate memory for tracer arrays 
    * - if multiple processors are defined, then 
    *   make array somewhat longer to allow small variations in
    *   length (due to tracer exchange between processors) 
    */
   if (E->parallel.nproc>1)
      lbuffer = (int)(E->tracers.lntrac * E->tracers.array_buf_factor) + 10;
   else
      lbuffer = 1;
   E->tracers.larraylen = E->tracers.lntrac + lbuffer;
   fprintf(E->fp,"tracer larraylen=%d (lntrac=%d + lbuffer=%d)\n",
           E->tracers.larraylen, E->tracers.lntrac, lbuffer);
   len = E->tracers.larraylen;
   nelem=E->lmesh.NEL[lev];
   
   /* coord arrays (2nd order Runge Kutta requires 2 sets) */
   /* PC-scheme requires another backup set (xtrac2) */
   /* coordinate arrays */
   E->tracers.xtrac =malloc(sizeof(float*)*4);
   E->tracers.xtrac2=malloc(sizeof(float*)*4);
   E->tracers.xtrac[1]=(float *) calloc(len,sizeof(float));
   E->tracers.xtrac[2]=(float *) calloc(len,sizeof(float));
   E->tracers.xtrac1[1]=(float *) calloc(len,sizeof(float));
   E->tracers.xtrac1[2]=(float *) calloc(len,sizeof(float));
   E->tracers.xtrac2[1]=(float *) calloc(len,sizeof(float));
   E->tracers.xtrac2[2]=(float *) calloc(len,sizeof(float));
   if (E->mesh.nsd == 3)
     {
       E->tracers.xtrac[3]=(float *) calloc(len,sizeof(float));
       E->tracers.xtrac1[3]=(float *) calloc(len,sizeof(float));
       E->tracers.xtrac2[3]=(float *) calloc(len,sizeof(float));
     }

   /* velocity arrays (also 2 sets) */
   E->tracers.vtrac[1]=(float *) calloc(len,sizeof(float));
   E->tracers.vtrac[2]=(float *) calloc(len,sizeof(float));
   if (E->mesh.nsd == 3)
       E->tracers.vtrac[3]=(float *) calloc(len,sizeof(float));
   E->tracers.vtrac1[1]=(float *) calloc(len,sizeof(float));
   E->tracers.vtrac1[2]=(float *) calloc(len,sizeof(float));
   if (E->mesh.nsd == 3)
       E->tracers.vtrac1[3]=(float *) calloc(len,sizeof(float));

   /* elemnr arrays */
   E->tracers.elemnr=(int *) calloc(len,sizeof(int));
   E->tracers.elemnr1=(int *) calloc(len,sizeof(int));
   E->tracers.ntrac_in_elem=(int *) calloc(nelem+1,sizeof(int));
   E->tracers.ntrac_in_elem1=(int *) calloc(nelem+1,sizeof(int));
   E->tracers.tracs_in_elem=(int *) calloc(len+1,sizeof(int));
   E->tracers.tracs_in_elem1=(int *) calloc(len+1,sizeof(int));
   E->tracers.elemtraccounter=(int *) calloc(nelem+1,sizeof(int));
   E->tracers.trapmelt=(float *) calloc(nox*noy+1,sizeof(float));
   for (ifunc=0; ifunc<(nox*noy+1); ifunc++) E->tracers.trapmelt[ifunc]=0.;

   /* 2-D function array ftrac of size 
    *      E->tracers.ntracfunc X E->tracers.larraylen
    */
   if (E->tracers.ntracfunc>0)
   {  /* tracer function arrays */

     E->tracers.ftrac=malloc(sizeof(float*)*NTRACFUNCMAX);
     E->tracers.ftrac2=malloc(sizeof(float*)*NTRACFUNCMAX);
     for (ifunc=0; ifunc<E->tracers.ntracfunc; ifunc++)
      { 
	E->tracers.ftrac[ifunc]=(float *) calloc(len+1,sizeof(float));
        E->tracers.ftrac2[ifunc]=(float *) calloc(len+1,sizeof(float));
        if (E->tracers.interp_t2e[ifunc]==1)
           E->tracers.felem[ifunc]=(float *) 
               calloc(nelem+1,sizeof(float));
      }

      /* compositional buoyancy field, calc. per element */
      if (E->tracers.ntracbuoyfunc > 0.)
        E->tracers.compobuoy = (float *) calloc(nelem+1,sizeof(float));

      E->tracers.melting_at_trac=0.;
      if (E->tracers.restart != 0 && E->control.restart != 0)       
	{
	  sprintf(filename,"%s.ave.%d.%d",
		  E->control.data_file,E->parallel.me,E->control.first_tstep);
	  fp = fopen(filename,"r");
	  if (fp != NULL) 
	   {
	    fgets(input_s,200,fp);
	    sscanf(input_s,"%d %d %g %g %g %g %g %g ",
		  &a,&b,&c[0],&c[1],&c[2],&c[3],&c[4],
		  &E->tracers.melting_at_trac);
            fclose(fp);
	   }
	}

      /* melttracer memory allocation and setup                   *
       *                                                          * 
       * melttracers are used for trapping extracted melt in the  *
       * x/y-plane at z=0, they are advected simply with plate    *
       * velocity and collect all the melt extracted in the       *
       * column below.                                            */
 
      if (E->tracers.melttracers_used != 0)
       {
	if (E->tracers.restart != 0 && E->control.restart != 0)
          restart_melttracers(E);
        else 
	 {
	  dens=E->tracers.melttracer_density;
	  E->tracers.melt_columns = (int)(dens*E->mesh.layer[1]);
	  E->tracers.melt_rows=1;
	  E->tracers.melt_xdistance=E->mesh.layer[1]/(E->tracers.melt_columns);
	  if (E->mesh.nsd==3)
	    {
	     E->tracers.melt_rows = (int)(dens*E->mesh.layer[3]);
	    }
	  E->tracers.melt_ydistance=E->mesh.layer[3]/E->tracers.melt_rows;
	  E->tracers.melt_firstcolumn=0;

	  for (first_upper_bound=E->tracers.melt_xdistance/2; first_upper_bound<xmin; first_upper_bound+=E->tracers.melt_xdistance)
	    {;}
	  E->tracers.melt_lcolumns=1;
       // fprintf (E->fp,"first_upper_bound %f; \n\n\n\n",first_upper_bound);
	  for (upper_bound=first_upper_bound; upper_bound<xmax; upper_bound+=E->tracers.melt_xdistance)
	    E->tracers.melt_lcolumns++;

          E->tracers.melt_xbound=(double *) calloc(E->tracers.melt_lcolumns+2,sizeof(double));
	  E->tracers.melt_xbound[0]=first_upper_bound-E->tracers.melt_xdistance;
	  for (n=1; n<=E->tracers.melt_lcolumns; n++)
	    E->tracers.melt_xbound[n]=first_upper_bound+E->tracers.melt_xdistance*((float)(n-1));
	  E->tracers.melttracer_flag=0;
	  E->tracers.melt_lrows=1;
	  if (E->mesh.nsd==3)
	    {
	      for (first_upper_bound=0.;first_upper_bound<=ymin;first_upper_bound+=E->tracers.melt_ydistance)
		if (first_upper_bound==ymin) E->tracers.melttracer_flag=1;
	      for (upper_bound=first_upper_bound; upper_bound<ymax; upper_bound+=E->tracers.melt_ydistance)
		E->tracers.melt_lrows++;

	      E->tracers.melt_ybound=(double *)calloc(E->tracers.melt_lrows+2,sizeof(double));
	      E->tracers.melt_ybound[0]=first_upper_bound-E->tracers.melt_ydistance;
fprintf(E->fp,"ybound:%f\n",E->tracers.melt_ybound[0]);
	      for (n=1; n<=E->tracers.melt_lrows; n++)
		E->tracers.melt_ybound[n]=first_upper_bound+E->tracers.melt_ydistance*((float)(n-1));
	    }
	  else
	      E->tracers.melt_ybound=(double *)calloc(3,sizeof(double));

          E->tracers.melt=malloc(sizeof(float*)*(E->tracers.melt_lcolumns+2));
	  for (n=0; n<=E->tracers.melt_lcolumns; n++)
	    {
		 E->tracers.melt[n]=malloc(sizeof(float*)*(E->tracers.melt_lrows+1));
	         for (n2=0; n2<E->tracers.melt_lrows; n2++)
		  E->tracers.melt[n][n2]=0.;
	    }
	}
       }
   }

   return;
}

void restart_melttracers(E)
   struct All_variables *E;
{
    int dens;
    float first_bound;
    void tracexit();
    int columns_in,rows,firstcolumn_in,flag,lc,n;
    float x,y=0.;
    char filename[100];
    FILE *fp;

    /*    +------------------------------------+    *
     *    | restarting melttracers: reading in |    *
     *    | information from binary files.     |    *
     *    | does not yet work for ascii-files. |    *
     *    +------------------------------------+    *    MDB VII/2006 */

    if (E->tracers.melttracers_used==1) 
      {
       sprintf(filename,"%s.melttrac.%d.%d",E->control.data_file,
                   E->parallel.me,E->control.first_tstep);

       fp = fopen(filename,"r");
       if (fp == NULL)
         tracexit(E, "Error in restart_melttracers: Unable to open %s\n: STOP",filename);
       else 
	 {
	  fread(&columns_in,1,sizeof(int),fp);
	  fread(&rows,1,sizeof(int),fp);
	  fread(&firstcolumn_in,1,sizeof(int),fp);
	  fread(&lc,1,sizeof(int),fp);
	  fread(&x,1,sizeof(float),fp);

          E->tracers.melt_xbound=(double *) calloc(columns_in+2,sizeof(double));
	  E->tracers.melt_ybound=(double *)calloc(rows+2,sizeof(double));

	  E->tracers.melt_lcolumns=columns_in;
	  E->tracers.melt_lrows=rows;
	  E->tracers.melt_xbound[firstcolumn_in]=x;
	  E->tracers.melt_firstcolumn=firstcolumn_in;
	  E->tracers.melt_lastcolumn=lc;

	  if (E->mesh.nsd==3)
	    {
	      if (rows==1)
                tracexit(E, "Error in restart_melttracers: lrows=1 && E->mesh.nsd=3\n: STOP");
	      fread(&y,1,sizeof(float),fp);
	      fread(&flag,1,sizeof(int),fp);
	      E->tracers.melt_ybound[0]=y;
	      E->tracers.melttracer_flag=flag;
	    }
	  else if (rows>1)
            tracexit(E, "Error in restart_melttracers: lrows>1 && E->mesh.nsd!=3\n: STOP");

	  dens=E->tracers.melttracer_density;
	  E->tracers.melt_columns = (int)(dens*E->mesh.layer[1]);
	  E->tracers.melt_rows=1;
	  E->tracers.melt_xdistance=E->mesh.layer[1]/(E->tracers.melt_columns);
	  if (E->mesh.nsd==3)
	    {
	      E->tracers.melt_rows = (int)(dens*E->mesh.layer[3]);
	    }
	  E->tracers.melt_ydistance=E->mesh.layer[3]/E->tracers.melt_rows;
	  
          first_bound=E->tracers.melt_xbound[firstcolumn_in];
	  for (n=firstcolumn_in+1; n<=columns_in; n++)
	   {
	     first_bound+=E->tracers.melt_xdistance;
	     E->tracers.melt_xbound[n]=first_bound;
	   }
	  for (n=0; n<firstcolumn_in; n++)
	   {
	     first_bound+=E->tracers.melt_xdistance;
	     E->tracers.melt_xbound[n]=first_bound;
	   }

	  if (E->mesh.nsd==3)
	    {
	   // fprintf(E->fp,"ybound:%f\n",E->tracers.melt_ybound[0]);
              first_bound=E->tracers.melt_ybound[0];
	      for (n=1; n<=rows; n++)
		E->tracers.melt_ybound[n]=first_bound+E->tracers.melt_ydistance*((float)(n));
	    }
          E->tracers.melt=malloc(sizeof(float*)*(columns_in+2));
	  for (n=0; n<=columns_in; n++)
	    {
		 E->tracers.melt[n]=malloc(sizeof(float*)*(rows+1));
	         fread(E->tracers.melt[n],rows,sizeof(float),fp);
	    }
          fclose(fp); 
	 }
    }
    else if (E->tracers.melttracers_used==2)
      tracexit (E, "Error in restart_melttracers: restart of ascii-files not yet implemented\n: STOP");

    return;
}

void parallel_tracer_setup(E)
   struct All_variables *E;
{
   int ineighb, len, lev, nox, noy, noz, me, nprocx, nprocy, nprocz; 
   int ndim, iprocx, iprocy, iprocz;
   int nneighb_max;

   void alloc_sr();

   len    = E->tracers.larraylen;
   lev    = E->mesh.levmax;
   nox    = E->lmesh.NOX[lev];
   noy    = E->lmesh.NOY[lev];
   noz    = E->lmesh.NOZ[lev];
   me     = E->parallel.me;
   nprocx = E->parallel.nprocx;
   nprocy = E->parallel.nprocy;
   nprocz = E->parallel.nprocz;
   iprocx = E->parallel.me_loc[1];
   iprocz = E->parallel.me_loc[2];
   iprocy = E->parallel.me_loc[3];

   ndim = E->mesh.nsd;
   if (ndim==2) 
      nneighb_max=8;
   else 
      nneighb_max=26;

   if (E->parallel.nproc > 1)
   { E->tracers.crossborder=(int *) calloc(len,sizeof(int));
     /* determine if proc has neighbor, and if so, its proc number */
      
     /* determine if proc has neighbor, and if so:
      * - set used: 0=unused, 1=used
      * - determine size: large enough size of arrays in struct to 
      *   fit all tracs that travel to neighbor, assuming that 
      *   less than the average # tracs in particular boundary or
      *   corner elements will cross bnd (according to CFL).
      * - determine procnr: using iproc=iz + ix*nz + iy*nz*nx 
      */

     for (ineighb=1; ineighb<27; ineighb++) /*default=not used */
        E->tracers.sr[1].used=0;

     if (E->parallel.me_loc[1]!=0) /* xmin-neighbors */
     {  E->tracers.sr[1].used=1; 
        E->tracers.sr[1].size=len/nox;
        E->tracers.sr[1].procnr = me - nprocz;
        if (E->parallel.me_loc[2]!=0) /* zmin nghbr*/
        {  E->tracers.sr[5].used=1;
           E->tracers.sr[5].size=len/nox/noz;
           E->tracers.sr[5].procnr = me - nprocz - 1;
        }
        if (E->parallel.me_loc[2]!=E->parallel.nprocz-1) /* zmax nghbr*/
        {  E->tracers.sr[7].used=1;
           E->tracers.sr[7].size=len/nox/noz;
           E->tracers.sr[7].procnr = me - nprocz + 1;
        }
     }
     if (E->parallel.me_loc[1]!=E->parallel.nprocx-1) /* xmax-neighb */
     {  E->tracers.sr[2].used=1;
        E->tracers.sr[2].size=len/nox;
        E->tracers.sr[2].procnr = me + nprocz;
        if (E->parallel.me_loc[2]!=0) /* zmin nghbr*/
        {  E->tracers.sr[6].used=1;
           E->tracers.sr[6].size=len/nox/noz;
           E->tracers.sr[6].procnr = me + nprocz - 1;
        }
        if (E->parallel.me_loc[2]!=E->parallel.nprocz-1) /* zmax nghbr*/
        {  E->tracers.sr[8].used=1;
           E->tracers.sr[8].size=len/nox/noz;
           E->tracers.sr[8].procnr = me + nprocz + 1;
        }
     }
     if (E->parallel.me_loc[2]!=0) /* zmin-neighbors */
     {  E->tracers.sr[3].used=1;
        E->tracers.sr[3].size=len/noz;
        E->tracers.sr[3].procnr = me - 1;
     }
     if (E->parallel.me_loc[2]!=E->parallel.nprocz-1) /* zmax-neighb */
     {  E->tracers.sr[4].used=1;
        E->tracers.sr[4].size=len/noz;
        E->tracers.sr[4].procnr = me + 1;
     }

     if (ndim == 3)
     {  if (E->parallel.me_loc[3]!=0) /* ymin-neighbors */
        {  E->tracers.sr[9].used=1; 
           E->tracers.sr[9].size=len/noy;
           E->tracers.sr[9].procnr = me - nprocz*nprocx;
           if (E->parallel.me_loc[1]!=0) /* xmin-neighbors */
           {  E->tracers.sr[11].used=1; 
              E->tracers.sr[11].size=len/noy/nox;
              E->tracers.sr[11].procnr=me-nprocz-nprocz*nprocx;
              if (E->parallel.me_loc[2]!=0) /* zmin nghbr*/
              {  E->tracers.sr[19].used=1;
                 E->tracers.sr[19].size=len/noy/nox/noz;
                 E->tracers.sr[19].procnr=me-1-nprocz-nprocz*nprocx;
              }
              if (E->parallel.me_loc[2]!=E->parallel.nprocz-1) 
                                                        /* zmax nghbr*/
              {  E->tracers.sr[21].used=1;
                 E->tracers.sr[21].size=len/noy/nox/noz;
                 E->tracers.sr[21].procnr=me+1-nprocz-nprocz*nprocx;
              }
           }
           if (E->parallel.me_loc[1]!=E->parallel.nprocx-1) 
                                                        /* xmax-nghb */
           {  E->tracers.sr[12].used=1;
              E->tracers.sr[12].size=len/noy/nox;
              E->tracers.sr[12].procnr = me+nprocz-nprocz*nprocx;
              if (E->parallel.me_loc[2]!=0) /* zmin nghbr*/
              {  E->tracers.sr[20].used=1;
                 E->tracers.sr[20].size=len/noy/nox/noz;
                 E->tracers.sr[20].procnr=me-1+nprocz-nprocz*nprocx;
              }
              if (E->parallel.me_loc[2]!=E->parallel.nprocz-1) 
                                                        /* zmax nghbr*/
              {  E->tracers.sr[22].used=1;
                 E->tracers.sr[22].size=len/noy/nox/noz;
                 E->tracers.sr[22].procnr=me+1+nprocz-nprocz*nprocx;
              }
           }
           if (E->parallel.me_loc[2]!=0) /* zmin-neighbors */
           {  E->tracers.sr[15].used=1;
              E->tracers.sr[15].size=len/noy/noz;
              E->tracers.sr[15].procnr=me-1-nprocz*nprocx;
           }
           if (E->parallel.me_loc[2]!=E->parallel.nprocz-1) 
                                                       /* zmax-nghbr */
           {  E->tracers.sr[16].used=1;
              E->tracers.sr[16].size=len/noy/noz;
              E->tracers.sr[16].procnr=me+1-nprocz*nprocx;
           }
        }
        if (E->parallel.me_loc[3]!=E->parallel.nprocy-1) /* ymax-nghb */
        {  E->tracers.sr[10].used=1; 
           E->tracers.sr[10].size=len/noy;
           E->tracers.sr[10].procnr = me + nprocz*nprocx;
           if (E->parallel.me_loc[1]!=0) /* xmin-neighbors */
           {  E->tracers.sr[13].used=1; 
              E->tracers.sr[13].size=len/noy/nox;
              E->tracers.sr[13].procnr=me-nprocz+nprocz*nprocx;
              if (E->parallel.me_loc[2]!=0) /* zmin nghbr*/
              {  E->tracers.sr[23].used=1;
                 E->tracers.sr[23].size=len/noy/nox/noz;
                 E->tracers.sr[23].procnr=me-1-nprocz+nprocz*nprocx;
              }
              if (E->parallel.me_loc[2]!=E->parallel.nprocz-1) 
                                                        /* zmax nghbr*/
              {  E->tracers.sr[25].used=1;
                 E->tracers.sr[25].size=len/noy/nox/noz;
                 E->tracers.sr[25].procnr=me+1-nprocz+nprocz*nprocx;
              }
           }
           if (E->parallel.me_loc[1]!=E->parallel.nprocx-1) 
                                                        /* xmax-nghb */
           {  E->tracers.sr[14].used=1;
              E->tracers.sr[14].size=len/noy/nox;
              E->tracers.sr[14].procnr=me+nprocz+nprocz*nprocx;
              if (E->parallel.me_loc[2]!=0) /* zmin nghbr*/
              {  E->tracers.sr[24].used=1;
                 E->tracers.sr[24].size=len/noy/nox/noz;
                 E->tracers.sr[24].procnr=me-1+nprocz+nprocz*nprocx;
              }
              if (E->parallel.me_loc[2]!=E->parallel.nprocz-1) 
                                                        /* zmax nghbr*/
              {  E->tracers.sr[26].used=1;
                 E->tracers.sr[26].size=len/noy/nox/noz;
                 E->tracers.sr[26].procnr=me+1+nprocz+nprocz*nprocx;
              }
           }
           if (E->parallel.me_loc[2]!=0) /* zmin-neighbors */
           {  E->tracers.sr[17].used=1;
              E->tracers.sr[17].size=len/noy/noz;
              E->tracers.sr[17].procnr=me-1+nprocz*nprocx;
           }
           if (E->parallel.me_loc[2]!=E->parallel.nprocz-1) 
                                                       /* zmax-nghbr */
           {  E->tracers.sr[18].used=1;
              E->tracers.sr[18].size=len/noy/noz;
              E->tracers.sr[18].procnr=me+1+nprocz*nprocx;
           }
        }
     }
   
     /* allocate memory for used structs: */
     for (ineighb=1; ineighb<=nneighb_max; ineighb++)
        if (E->tracers.sr[ineighb].used==1) 
           alloc_sr(E,ineighb,E->tracers.sr[ineighb].size);

     /* the order in which neighbors are called to exchange tracer 
        info is either parallel_order1 (if iprocx+iprocz+iprocy=even) 
        or parallel_order2 (if iprocx+iprocz+iprocy=odd)
        (see tracers.h). */
     if ((iprocx + iprocz + iprocy) % 2 == 0)
        for (ineighb=0; ineighb<=nneighb_max; ineighb++)
           E->tracers.order[ineighb] = parallel_order1[ineighb];
     else
        for (ineighb=0; ineighb<=nneighb_max; ineighb++)
           E->tracers.order[ineighb] = parallel_order2[ineighb];
     fprintf(E->fp,"order: %d %d %d %d etc.\n",
            E->tracers.order[1], E->tracers.order[2],
            E->tracers.order[3], E->tracers.order[4]);
   }
   return;
}

void alloc_sr(E,ineighb, lbuffer)
   struct All_variables *E;
   int ineighb, lbuffer;
{  
   int ifunc;

   fprintf(E->fp,"Allocating %d memory for sendrec to neighbor %d\n", 
           lbuffer, ineighb);
   E->tracers.sr[ineighb].tracnr= 
         (int *) calloc(lbuffer,sizeof(int));
   E->tracers.sr[ineighb].elemnr= 
         (int *) calloc(lbuffer,sizeof(int));
   E->tracers.sr[ineighb].vec1[1]= 
         (float *) calloc(lbuffer,sizeof(float));
   E->tracers.sr[ineighb].vec1[2]= 
         (float *) calloc(lbuffer,sizeof(float));
   if (E->mesh.nsd == 3)
      E->tracers.sr[ineighb].vec1[3]=
           (float *) calloc(lbuffer,sizeof(float));
   E->tracers.sr[ineighb].vec2[1]=
        (float *) calloc(lbuffer,sizeof(float));
   E->tracers.sr[ineighb].vec2[2]=
        (float *) calloc(lbuffer,sizeof(float));
   if (E->mesh.nsd == 3)
      E->tracers.sr[ineighb].vec2[3]=
           (float *) calloc(lbuffer,sizeof(float));
   for (ifunc=0; ifunc<E->tracers.ntracfunc; ifunc++)
   {  E->tracers.sr[ineighb].func1[ifunc]= 
         (float *) calloc(lbuffer,sizeof(float));
      E->tracers.sr[ineighb].func2[ifunc]= 
         (float *) calloc(lbuffer,sizeof(float));
   }
   E->tracers.sr[ineighb].size = lbuffer;
   return;
}

/* Grow MPI send/recv buffers for one neighbor (must match alloc_sr layout). */
static void grow_sr(E, ineighb, new_size)
   struct All_variables *E;
   int ineighb, new_size;
{
   int ifunc, old_size;
   void tracexit();

#define GROW_INT(ptr, n, old_n) \
   do { \
      (ptr) = (int *) realloc((ptr), (n) * sizeof(int)); \
      if (!(ptr)) tracexit(E, "grow_sr: out of memory (int)"); \
      if ((n) > (old_n)) \
         memset((ptr) + (old_n), 0, ((n) - (old_n)) * sizeof(int)); \
   } while (0)

#define GROW_FLOAT(ptr, n, old_n) \
   do { \
      (ptr) = (float *) realloc((ptr), (n) * sizeof(float)); \
      if (!(ptr)) tracexit(E, "grow_sr: out of memory (float)"); \
      if ((n) > (old_n)) \
         memset((ptr) + (old_n), 0, ((n) - (old_n)) * sizeof(float)); \
   } while (0)

   if (new_size <= 0) return;
   old_size = E->tracers.sr[ineighb].size;
   if (old_size <= 0) {
      alloc_sr(E, ineighb, new_size);
      return;
   }
   if (new_size <= old_size) {
      E->tracers.sr[ineighb].size = new_size;
      return;
   }

   fprintf(E->fp,
           "grow_sr: neighbor %d buffer %d -> %d (proc %d)\n",
           ineighb, old_size, new_size, E->parallel.me);
   fprintf(stderr,
           "grow_sr: neighbor %d buffer %d -> %d (proc %d)\n",
           ineighb, old_size, new_size, E->parallel.me);

   GROW_INT(E->tracers.sr[ineighb].tracnr, new_size, old_size);
   GROW_INT(E->tracers.sr[ineighb].elemnr, new_size, old_size);
   GROW_FLOAT(E->tracers.sr[ineighb].vec1[1], new_size, old_size);
   GROW_FLOAT(E->tracers.sr[ineighb].vec1[2], new_size, old_size);
   GROW_FLOAT(E->tracers.sr[ineighb].vec2[1], new_size, old_size);
   GROW_FLOAT(E->tracers.sr[ineighb].vec2[2], new_size, old_size);
   if (E->mesh.nsd == 3) {
      GROW_FLOAT(E->tracers.sr[ineighb].vec1[3], new_size, old_size);
      GROW_FLOAT(E->tracers.sr[ineighb].vec2[3], new_size, old_size);
   }
   for (ifunc = 0; ifunc < E->tracers.ntracfunc; ifunc++) {
      GROW_FLOAT(E->tracers.sr[ineighb].func1[ifunc], new_size, old_size);
      GROW_FLOAT(E->tracers.sr[ineighb].func2[ifunc], new_size, old_size);
   }
   E->tracers.sr[ineighb].size = new_size;

#undef GROW_INT
#undef GROW_FLOAT
}

void trac_distr(E)
   /* This routine fills coordinate arrays of tracers by either:
      - random/non-random placement  (tracer_restart=0)
      - reading from ascii file (tracer_restart=1) 
      - reading from binary file (tracer_restart=2) */
   struct All_variables *E;
{
  int itrac, nfunc, ndim, n;
   float x, y, z;
   float xmin, xmax, ymin, ymax, zmin, zmax;
   float eps,dens;
   float dz,dx;
   float xx, yy, zz, dy, dzz, angle, tanangle, depth;
   int iy, iz, ny, nz, nx, ix, count;
   const int dims = E->mesh.nsd;   

   FILE *fp;
   char input_s[200], input_file[255];

   double drand48();
   void srand48();
   void read_trac_input();
   void vel_to_trac();
   void printtrac();
   void tracexit();
 
   xmin = E->tracers.lxmin;
   ymin = E->tracers.lymin;
   zmin = E->tracers.lzmin;
   xmax = E->tracers.lxmax;
   ymax = E->tracers.lymax;
   //zmax = E->tracers.lzmax;
   zmax = E->XX[E->mesh.levmax][2][E->tracers.bot_inflow];

   eps = 0.0001;

   if (E->tracers.restart == 1) 
   { /* Restart option: read tracer locations from file */
      sprintf(input_file,"%s.tracx.%d.%d",E->convection.old_T_file,
                        E->parallel.me, E->control.first_tstep);
      fprintf(stderr, "Reading tracer coords from restart file %s\n",
              input_file);
      fp = fopen(input_file,"r");
      fgets(input_s,200,fp);
      sscanf(input_s,"%d",&E->tracers.lntrac);
      if (E->tracers.lntrac >= E->tracers.larraylen)
      {  fprintf(stderr, "Error: proc %d has %d tracers, but ",
               E->parallel.me, E->tracers.lntrac);
         fprintf(stderr, "allocated memory is %d\n",
               E->tracers.larraylen);
         tracexit(E, "Error in tracdistr: lntrac>larraylen\n: STOP");
      }
      for (itrac=0; itrac<E->tracers.lntrac; itrac++)
      {  fgets(input_s,200,fp);
         sscanf(input_s,"%g %g",&x, &z);
         E->tracers.xtrac[1][itrac] = x;
         E->tracers.xtrac[2][itrac] = z;
      }
      fclose(fp);
   }
   else if (E->tracers.restart == 2) 
   { /* Restart option: read tracer locations and values from binary file */
      sprintf(input_file,"%s.tracval.%d.%d",E->convection.old_T_file,
                        E->parallel.me, E->control.first_tstep);
      fprintf(stderr, "Reading tracer coords from restart file %s\n",
              input_file);
      fp = fopen(input_file,"r");

      fread(&E->tracers.lntrac,1,sizeof(int),fp);
      fread(&nfunc,1,sizeof(int),fp);
      fread(&ndim,1,sizeof(int),fp);
 
      fprintf(stderr, "read in: ntrac: %d\n",
              E->tracers.lntrac);
                 
      if (E->tracers.lntrac >= E->tracers.larraylen)
      {  fprintf(stderr, "Error: proc %d has %d tracers, but ",
               E->parallel.me, E->tracers.lntrac);
         fprintf(stderr, "allocated memory is %d\n",
               E->tracers.larraylen);
         tracexit(E, "Error in tracdistr: lntrac>larraylen\n: STOP");
      }
      if (nfunc != E->tracers.ntracfunc)
	{
	  fprintf (stderr, "Error in tracdistr: nfunc (%i) does not correspond to restarted calculation (%i)\n",E->tracers.ntracfunc,nfunc);
	  tracexit (E, "Error in tracdistr: nfunc does not correspond to restarted calculation\n:STOP");
	}
      if (ndim != E->mesh.nsd)
	{
	  fprintf (stderr, "Error in tracdistr: dimensions (%i) do not correspond to restarted calculation (%i)\n",E->mesh.nsd,ndim);
         tracexit(E, "Error in tracdistr: dimensions do not correspond to restarted calculation\n: STOP");	  
	}
      for (n=1; n<=ndim; n++)
	fread(E->tracers.xtrac[n],E->tracers.lntrac,sizeof(float),fp);
      for (n=0; n<nfunc; n++)
	fread(E->tracers.ftrac[n],E->tracers.lntrac,sizeof(float),fp);
      fclose(fp);
   }

   // if (E->tracers.restart == 0):
   else
   {  /* Define new tracer locations */
      if (E->tracers.distr == 0)
      { /* Random tracer distribution */
         srand48(E->parallel.me);
         for (itrac=0; itrac<E->tracers.lntrac; itrac++)
         {  x = xmin + (xmax-xmin)*drand48();
            E->tracers.xtrac[1][itrac] = x;
            z = zmin + (zmax-zmin)*drand48();
            E->tracers.xtrac[2][itrac] = z;
            if (E->mesh.nsd == 3) 
            {  y = ymin + (ymax-ymin)*drand48();
               E->tracers.xtrac[3][itrac] = y;
            }
            /* Testing:
            if (E->mesh.nsd == 2) 
               fprintf(E->fp, "tracer x=%g, z=%g\n", 
                  E->tracers.xtrac[1][itrac], E->tracers.xtrac[2][itrac]);
            else 
               fprintf(E->fp, "tracer x=%g, y=%g, z=%g\n", 
                  E->tracers.xtrac[1][itrac], E->tracers.xtrac[3][itrac],  
                  E->tracers.xtrac[2][itrac]);
            */
         }
      }
      else if (E->tracers.distr == 1)
      { /* Non-random initial tracer location for testing purposes */
         for (itrac=0; itrac<E->tracers.lntrac; itrac++)
         {  E->tracers.xtrac[1][itrac] = 0.5;
            E->tracers.xtrac[2][itrac] = 0.25;
            if (E->mesh.nsd == 3) 
               E->tracers.xtrac[3][itrac] = 0.25;
         }
      }
      else if (E->tracers.distr == 2)
      { /* Non-random initial tracer location for testing purposes */
         for (itrac=0; itrac<E->tracers.lntrac; itrac++)
         {  E->tracers.xtrac[1][itrac] = 0.01;
            E->tracers.xtrac[2][itrac] = 0.01;
            if (E->mesh.nsd == 3) 
               E->tracers.xtrac[3][itrac] = 0.01;
         }
      }
      else if (E->tracers.distr == 3)
      { /* Non-random initial tracer location for testing purposes */
         for (itrac=0; itrac<E->tracers.lntrac; itrac++)
         {  E->tracers.xtrac[1][itrac] = 
                   max(min(0.25, xmax-eps),xmin+eps);
            E->tracers.xtrac[2][itrac] = 
                   max(min(0.25, zmax-eps),zmin+eps);
            if (E->mesh.nsd == 3) 
               E->tracers.xtrac[3][itrac] = 
                   max(min(0.25, ymax-eps),ymin+eps);
            if (E->mesh.nsd == 2) 
               fprintf(E->fp, "trac %d: (%g, %g)\n", itrac, 
                  E->tracers.xtrac[1][itrac],
                  E->tracers.xtrac[2][itrac]);
            else 
               fprintf(E->fp, "trac %d: (%g, %g %g)\n", itrac, 
                  E->tracers.xtrac[1][itrac],
                  E->tracers.xtrac[2][itrac],
                  E->tracers.xtrac[3][itrac]);
         }
      }
      else if (E->tracers.distr == 4)
        {/* Random tracer distribution, no tracers beyond the fence */
          srand48(E->parallel.me);
          for (itrac=0; itrac<E->tracers.lntrac; itrac++)
          {  x = xmin + (xmax-xmin)*drand48();
           E->tracers.xtrac[1][itrac] = x;
           z = zmin+eps + (zmax-zmin-2*eps)*drand48();
           E->tracers.xtrac[2][itrac] = z;
           if (E->mesh.nsd == 3) 
           {  y = ymin+eps + (ymax-ymin-2*eps)*drand48();
              E->tracers.xtrac[3][itrac] = y;
	   }
	  }
	}
      else if (E->tracers.distr == 10)
      { /* Non-random initial tracer location  */
         /* NB only 2-D yet */
         /* Set tracer type to 2: only move in hor. plane */
         E->tracers.type = 2;
         fprintf(E->fp, "distr 10: Subd fault tracers:\n");
         depth = 0.2;
         xx = 5.0;
         angle=45.; /* NB angle cannot be 0 !! */
         tanangle = 1./tan(angle/180.*3.1415926);
         dzz = depth/(E->tracers.lntrac-1);
         E->tracers.dzz=dzz;
         zz=0.;
         for (itrac=0; itrac<E->tracers.lntrac; itrac++)
         {  E->tracers.xtrac[1][itrac] = xx;
            E->tracers.xtrac[2][itrac] = zz;
            xx += tanangle*dzz;
            zz += dzz;
            fprintf(E->fp, "trac %d: (x,y)=(%g,%g)\n", itrac, 
                 E->tracers.xtrac[1][itrac], E->tracers.xtrac[2][itrac]);
         }
      }
      else if (E->tracers.distr == 11)
      { /* Non-random initial tracer location  */
         /* NB only 2-D yet */
         /* Set tracer type to 2: only move in hor. plane */
         E->tracers.type = 2;
         fprintf(E->fp, "distr 11: 3-D Subd fault surface tracers:\n");
         ymin = 0.0001;
         ymax = 1.9999;
         dy = (ymax-ymin)/(E->tracers.lntrac-1);
         yy=ymin;
         for (itrac=0; itrac<E->tracers.lntrac; itrac++)
         {  xx = 5.0-0.5*yy;
            E->tracers.xtrac[1][itrac] = xx;
            E->tracers.xtrac[2][itrac] = 0.;
            E->tracers.xtrac[3][itrac] = yy;
            yy += dy;
            fprintf(E->fp, "trac %d: (x,y,z)=(%g,%g,%g)\n", itrac, 
                 E->tracers.xtrac[1][itrac], 
                 E->tracers.xtrac[3][itrac],
                 E->tracers.xtrac[2][itrac]);
         }
      }
      
      else if (E->tracers.distr == 20)
        {  
         if (E->parallel.me !=0) E->tracers.lntrac=0;
   
         nz=10;
         dz=0.05;
         ny=E->tracers.lntrac/nz;
         dy=ymax/(ny+1);
         E->tracers.lntrac = nz*ny;
         iz=1; iy=1;
         for (itrac=0; itrac<E->tracers.lntrac; itrac++)
         {  E->tracers.xtrac[1][itrac] = 0.01;
            E->tracers.xtrac[2][itrac] = iz*dz;
            E->tracers.xtrac[3][itrac] = iy*dy;
            iz++;
            if (iz>nz)
            { iz-=nz;
              iy++;
            }
            fprintf(E->fp,"itr=%d, x=%g,z=%g,y=%g\n", itrac,
              E->tracers.xtrac[1][itrac], E->tracers.xtrac[2][itrac], 
              E->tracers.xtrac[3][itrac]);
         }
	}
      
      /*tad 4-17-06*/
      else if (E->tracers.distr == 1000)
      {
        count=0;
        dens = E->tracers.tracer_density;
        nx=(int) (dens*(xmax-xmin));
        //nx=(int) (dens*(xmax-(E->parallel.me_loc[1]*E->mesh.layer[1]/E->parallel.nprocx)));
        dx=(xmax-xmin)/nx;  //1/dens;
        nz=(int) (dens*(zmax-zmin));
        dz=(zmax-zmin)/nz;  //1/dens;
        ny=1;
        dy=0;
        if(dims==3){
          ny=(int) (dens*(ymax-ymin));
          dy=(ymax-ymin)/ny; //1/dens;
        }
        ix=0;
        iz=0;
        iy=0;
        for(iy=0; iy<ny; iy++)
          {
            for(iz=0; iz<nz; iz++)
              {
                for (ix=0; ix<nx; ix++)
                  {
                      /*ix + nx*iz + nx*nz*iy*/
                    E->tracers.xtrac[1][count] = dx*(ix+.5) +eps + xmin + dx*(drand48()-.5);
                     /*(E->parallel.me_loc[1]*E->mesh.layer[1]/E->parallel.nprocx) ;*/
                    E->tracers.xtrac[2][count] = dz*(iz+.5) +eps +zmin + dz*(drand48()-.5);
                    if (dims==3)
                      E->tracers.xtrac[3][count] = dy*(iy+.5) +eps +ymin + dy*(drand48()-.5);

                    count++;
                    /*if(E->parallel.me==1)
                      fprintf(stderr,"on proc %d, count=%d",E->parallel.me,count);*/
                  }
              }
          }
        /*fprintf(stderr,"on proc %d, count=%d and lntrac=%d, nx=%d dx=%.5e nz
=%d dz=%.5e X=%.10e Z=%.10e \n",
                E->parallel.me,count,E->tracers.lntrac,nx,dx,nz,dz,E->tracers.
xtrac[1][count],
                E->tracers.xtrac[2][count]);*/
        E->tracers.lntrac=count;
      }
   }

   /*
   fprintf(E->fp, "Allocated arrays with length %d\n", len);
   fprintf(E->fp, "xmin=%g, xmax=%g\n",xmin,xmax);
   fprintf(E->fp, "zmin=%g, zmax=%g\n",zmin,zmax);
   if (E->mesh.nsd == 3) 
      fprintf(E->fp, "ymin=%g, ymax=%g\n",ymin,ymax);
   fprintf(E->fp, "proc=%d, ntrac=%d, lntrac=%d\n",
      E->parallel.me, E->tracers.ntrac, E->tracers.lntrac);
   */

   /* Test tracer interpolation with V_node[idim] = 10*X_node[idim] 
   for (i=1;i<=E->lmesh.nno;i++)  
   {  E->V[1][i]=10*E->X[1][i];
      E->V[2][i]=10*E->X[2][i];
      if (E->mesh.nsd == 3) 
         E->V[3][i]=10*E->X[3][i];
   }
   vel_to_trac(E, E->V, E->tracers.xtrac, E->tracers.vtrac, 
                  E->tracers.elemnr, E->tracers.lntrac);
    */
   
   /* tracexit(E,"inittracers: TEST"); */
   return;
}

void trac_arrange_funcinit(E)
   struct All_variables *E;
{
   int ifunc, itrac, itracfunc;
   float x, y, z, val;
   FILE *fp;
   char filename[100], input_s[200];
   float trac_funcinit();

   /* For each function, fill array with tracer values: */

   if (E->tracers.restart == 1)
   { /* Restart option: read tracer values from file */
      if (E->parallel.nproc == 1) 
      {  /* No parallel job, so lntrac and tracer order are constant
            If ftrac doesn't change storing ftrac once is enough */
         for (ifunc=0; ifunc<E->tracers.ntracfunc; ifunc++)
         {  sprintf(filename,"%s.tracf.%d.%d",E->convection.old_T_file,
                           E->parallel.me, ifunc);
            fprintf(stderr, "Read tracer values from restart file %s\n",
                    filename);
            fp = fopen(filename,"r");
            for (itrac=0; itrac<E->tracers.lntrac; itrac++)
            {  fgets(input_s,200,fp);
               sscanf(input_s,"%g",&val);
               E->tracers.ftrac[ifunc][itrac] = val;
            }
            fclose(fp);
         }
      }
      else 
      {  /* Parallel job, so lntrac and tracer order are not constant
            ftrac is stored each output and latest should be read in */
         for (ifunc=0; ifunc<E->tracers.ntracfunc; ifunc++)
         {  sprintf(filename,"%s.tracf.%d.%d.%d",
                           E->convection.old_T_file,
                           E->parallel.me, ifunc, 
                           E->control.first_tstep);
            fprintf(stderr, "Read tracer values from restart file %s\n",
                    filename);
            fp = fopen(filename,"r");
            for (itrac=0; itrac<E->tracers.lntrac; itrac++)
            {  fgets(input_s,200,fp);
               sscanf(input_s,"%g",&val);
               E->tracers.ftrac[ifunc][itrac] = val;
            }
            fclose(fp);
         }
      }
   }
   else if (E->tracers.restart == 0) 
   {  /* Define trac-values functionally */
      for (ifunc=0; ifunc<E->tracers.ntracfunc; ifunc++)
      {  for (itrac=0; itrac<E->tracers.lntrac; itrac++)
         {  x = E->tracers.xtrac[1][itrac];
            z = E->tracers.xtrac[2][itrac];
            if (E->mesh.nsd == 3) 
               y = E->tracers.xtrac[3][itrac];
            else
               y = 0.;
            itracfunc = E->tracers.tracfunc[ifunc];
            val = trac_funcinit(E, itracfunc, x, y, z);
            E->tracers.ftrac[ifunc][itrac] = val;
         }
         if (E->tracers.tractest == 6 ) /* Print tracer func */
         {  sprintf(filename,"%s.tracf.%d.%d",E->control.data_file,
                         E->parallel.me,ifunc);
            fp=fopen(filename,"w");
            for (itrac=0; itrac < E->tracers.lntrac; itrac++)
            {  fprintf(fp,"%g\n", E->tracers.ftrac[ifunc][itrac]);
            }
            fclose(fp);
         }
      }
   }
   return;
}

void katz_tables(E)
   struct All_variables *E;
{
 float Temp, TsolP, TlherzLiq=1475.;
 float Tprime,deltaT,FcpxOut=1.;
 float pressure, F1, F2; 
 int p,T,i;
 float H2O;

 float calc_H2O();

 pressure= -0.025;
 for (p=0;p<=101;p++)
  {
  H2O=calc_H2O(E,0.);
  deltaT=43*pow(E->data.H2O,0.75);
  pressure+=0.05;
  F2=0.;
  FcpxOut = 0.15/(0.5+0.08*pressure);
  TsolP = 1085.7 + (132.9*pressure) - (5.1 * pressure*pressure);
  TlherzLiq = 1475.0+(80.0*pressure)-(3.2 * pressure*pressure);
  T=0;
  E->tracers.Fcorr[p][0]=1.+TsolP-deltaT;
  for (Temp=.9+TsolP-deltaT; F2<=FcpxOut; Temp+=2.)
   {
   T++;
    for (i=0;i<16;i++)
       { 
	deltaT = 43*pow(H2O,0.75);
        Tprime = (Temp - TsolP + deltaT) / (TlherzLiq - TsolP);
        F1 = Tprime*sqrt(Tprime);
        H2O+=calc_H2O(E,F1);
	H2O*=0.5;
       }
      for (i=0;i<16;i++)
       {
	deltaT = 43*pow(H2O,0.75);
        Tprime = (Temp+0.2 - TsolP + deltaT) / (TlherzLiq - TsolP);
        F2 = Tprime*sqrt(Tprime);
        H2O+=calc_H2O(E,F2);
	H2O*=0.5;
       }
     if (F2<=FcpxOut)
      {
       /*   1/(1 + dH/dT )  */
       E->tracers.Fcorr[p][T]=1/(1+((F2-F1)*E->data.heat_of_melt*5./E->data.Cp));
      }
   }          // T
  
  }           // p
}

float trac_funcinit(E, itracfunc, x, y, z)
   struct All_variables *E;
   int itracfunc;
   float x, y, z;
{
   static int imax,visits;
   FILE *fp;
   int i;
   char input_s[40],filename[140];
   static float rmd[52];
   static float depth[10003];
   static float density[10003];
   float val, val2, temp;
   void tracexit();
   const float factor = E->data.layer_km*E->data.density*E->data.grav_acc*1e-6;
   const float time_scale = E->data.layer_km*E->data.layer_km*
                            1.0e6/E->data.therm_diff;
   float katz_melting();
   void tracexit(); 

//   if  (((E->monitor.solution_cycles+1)%200)==50)
//     fprintf(stderr,"%d %g %g %g new shit\n",itracfunc, x, y, z);
   val = -12345.;
   if (itracfunc == 0) /* don't apply anything*/
      val = -12345.;
   else if (itracfunc == 1) /* test: ival=x */
      val = x;
   else if (itracfunc == 2) /* test: ival=z */
      val = z; 
   else if (itracfunc == 9) /* test: heavy stuff left */
   {  if (x > 0.5) 
         val=1.;
      else 
         val=-1.;
   }
   else if (itracfunc == 10) /* test: heavy stuff left */
   {  if (x < 0.5) 
         val=1.;
      else 
         val=-1.;
   }
   else if (itracfunc == 11) /* test: heavy stuff on top */
   {  if (z < (0.35+0.3*x)) 
         val=1.;
      else 
         val=0.;
   }
   else if (itracfunc == 12) /* test: heavy stuff on top */
   {  if (z < (0.65-0.3*x)) 
         val=1.;
      else 
         val=0.;
   }
   else if (itracfunc == 14) /* Magma Ocean Brown et al (2014) */
   {
     if (visits==0)
           { sprintf(filename,"%s.trac_init_14.input",E->control.data_file);
             visits=1;
             fp = fopen(filename,"r");
             if (fp==NULL)
               tracexit(E,"tracfunc=14: Unable to load trac_init_14.input: Exit");
             fgets(input_s,200,fp);
             sscanf(input_s,"%d",&imax);
             if (imax<5 || imax>1000) tracexit(E,"tracfunc=14: first entry in trac_init_14.input inappropriate: Exit");
             for (i=1; i<=imax; i++)
               { fgets(input_s,200,fp);
                 sscanf(input_s,"%g %g",&val, &val2);
                 depth[i] = val;
                 density[i] = val2;
               }
             fclose(fp);
             val=-12345;
           }
      if (z > depth[imax/2])
       for (i=imax/2+1; i<=imax; i++)
        {
           if (depth[i-1]<=z && depth[i]>=z)
            {
               val =(density[i-1]*(depth[i]-z) + density[i]*(z-depth[i-1]))/ (depth[i]-depth[i-1]);
               break;
            }
        }
      else
       for (i=2; i<=imax/2; i++)
        {
           if (depth[i-1]<=z && depth[i]>=z)
            {
               val =(density[i-1]*(depth[i]-z) + density[i]*(z-depth[i-1]))/ (depth[i]-depth[i-1]);
               break;
            }
        }
   }
  else if (itracfunc == 15) /* Magma Ocean regular-step Rayleigh fractionation */
   {
     if (E->control.init_age_distr == 160 || E->control.BMO_implicit_mode == 1)
       {
         /* Implicit BMO coupling (init_age_distr=160): no pre-computed
          * fractionation profile is needed.  Initialise tracers with a
          * uniform composition of 1 (no prior Rayleigh fractionation). */
         val = 1.0;
       }
     else
       {
         if (visits==0)
               { sprintf(filename,"%s.trac_init_15.input",E->control.data_file);
                 visits=1;
                 fp = fopen(filename,"r");
                 if (fp==NULL)
                   tracexit(E,"tracfunc=15: Unable to load trac_init_15.input: Exit");
                 fgets(input_s,200,fp);
                 sscanf(input_s,"%d",&imax);
                 if (imax<5 || imax>10001) tracexit(E,"tracfunc=15: first entry in trac_init_15.input inappropriate: Exit");
                 for (i=imax; i>=1; i--)
                   { fgets(input_s,200,fp);
                     sscanf(input_s,"%g",&val2);
                     density[i] = val2;
		     depth[i]=((float)(i-1))/((float)(imax-1));
		     fprintf(E->fp, "test output itracfunc 15 %d %g\n", i,depth[i]);
                   }
                 fclose(fp);
                 val=-12345;
                 imax--;
               }
            i = 2 + (int)(z*((float)(imax)));
            //     fprintf(E->fp, "Ttesto-- %d %g\n", i,z);
                   val =(density[i-1]*(depth[i]-z) + density[i]*(z-depth[i-1]))/ (depth[i]-depth[i-1]);
       }
   }
  else if (itracfunc == 16) /* Magma Ocean regular-step Rayleigh fractionation */
   {
     if (visits==0)
           { sprintf(filename,"%s.trac_init_16.input",E->control.data_file);
             visits=1;
             fp = fopen(filename,"r");
             if (fp==NULL)
               tracexit(E,"tracfunc=16: Unable to load trac_init_16.input: Exit");
             fgets(input_s,200,fp);
             sscanf(input_s,"%d",&imax);
             if (imax<5 || imax>10001) tracexit(E,"tracfunc=16: first entry in trac_init_16.input inappropriate: Exit");
             for (i=imax; i>=1; i--)
               { fgets(input_s,200,fp);
                 sscanf(input_s,"%g",&val2);
                 density[i] = val2;
                 depth[i]=((float)(i-1))/((float)(imax-1));
                 fprintf(E->fp, "test output itracfunc 16: %d %g %g\n", i,depth[i],density[i]);
               }
             fclose(fp);
             val=-12345;
             imax--;
           }
	if (z<E->control.protoM_depth) {
               z/=E->control.protoM_depth;
               i = 2 + (int)(z*((float)(imax)));
               val =(density[i-1]*(depth[i]-z) + density[i]*(z-depth[i-1]))/ (depth[i]-depth[i-1]); 
//             fprintf(E->fp, "Ttesto-- %d %g %g\n", i,z,val);
//             fprintf(E->fp, "Ttesto-- %d %g %g %g %g %g %g\n", i,z,depth[i],depth[i-1],density[i],density[i-1],val); 
           }
        else val=E->control.protoM_density;
   }
  else if (itracfunc == 17) /* BMO regular-step Rayleigh fractionation */
   {
     if (visits==0)
           { sprintf(filename,"%s.trac_init_17.input",E->control.data_file);
             visits=1;
             fp = fopen(filename,"r");
             if (fp==NULL)
               tracexit(E,"tracfunc=17: Unable to load trac_init_17.input: Exit");
             fgets(input_s,200,fp);
             sscanf(input_s,"%d",&imax);
             if (imax<5 || imax>10001) tracexit(E,"tracfunc=16: first entry in trac_init_16.input inappropriate: Exit");
             for (i=imax; i>=1; i--)
               { fgets(input_s,200,fp);
                 sscanf(input_s,"%g",&val2);
                 density[i] = val2;
                 depth[i]=((float)(i-1))/((float)(imax-1));
                 fprintf(E->fp, "test output itracfunc 16: %d %g %g\n", i,depth[i],density[i]);
               }
             fclose(fp);
             val=-12345;
             imax--;
           }
        if (z>=E->control.protoM_depth) {
             z=1.0-z;
             z/=(1-E->control.protoM_depth);
             i = 2 + (int)(z*((float)(imax)));
             val =(density[i-1]*(depth[i]-z) + density[i]*(z-depth[i-1]))/ (depth[i]-depth[i-1]);
//           fprintf(E->fp, "Ttesto-- %d %g %g\n", i,z,val);
//           fprintf(E->fp, "Ttesto-- %d %g %g %g %g %g %g\n", i,z,depth[i],depth[i-1],density[i],density[i-1],val); 
           }
        else val=E->control.protoM_density;
   }
  else if (itracfunc == 18) /* BMO regular-step trace element A */
   {
     if (visits==0)
           { sprintf(filename,"%s.trac_init_18.input",E->control.data_file);
             visits=1;
             fp = fopen(filename,"r");
             if (fp==NULL)
               tracexit(E,"tracfunc=17: Unable to load trac_init_17.input: Exit");
             fgets(input_s,200,fp);
             sscanf(input_s,"%d",&imax);
             if (imax<5 || imax>10001) tracexit(E,"tracfunc=16: first entry in trac_init_16.input inappropriate: Exit");
             for (i=imax; i>=1; i--)
               { fgets(input_s,200,fp);
                 sscanf(input_s,"%g",&val2);
                 density[i] = val2;
                 depth[i]=((float)(i-1))/((float)(imax-1));
                 fprintf(E->fp, "test output itracfunc 16: %d %g %g\n", i,depth[i],density[i]);
               }
             fclose(fp);
             val=-12345;
             imax--;
           }
        if (z>=E->control.protoM_depth) {
             z=1.0-z;
             z/=(1-E->control.protoM_depth);
             i = 2 + (int)(z*((float)(imax)));
             val =(density[i-1]*(depth[i]-z) + density[i]*(z-depth[i-1]))/ (depth[i]-depth[i-1]);
//           fprintf(E->fp, "Ttesto-- %d %g %g\n", i,z,val);
//           fprintf(E->fp, "Ttesto-- %d %g %g %g %g %g %g\n", i,z,depth[i],depth[i-1],density[i],density[i-1],val); 
           }
        else val=1.0;
   }
  else if (itracfunc == 19) /* BMO regular-step trace element B */
   {
     if (visits==0)
           { sprintf(filename,"%s.trac_init_17.input",E->control.data_file);
             visits=1;
             fp = fopen(filename,"r");
             if (fp==NULL)
               tracexit(E,"tracfunc=17: Unable to load trac_init_17.input: Exit");
             fgets(input_s,200,fp);
             sscanf(input_s,"%d",&imax);
             if (imax<5 || imax>10001) tracexit(E,"tracfunc=16: first entry in trac_init_16.input inappropriate: Exit");
             for (i=imax; i>=1; i--)
               { fgets(input_s,200,fp);
                 sscanf(input_s,"%g",&val2);
                 density[i] = val2;
                 depth[i]=((float)(i-1))/((float)(imax-1));
                 fprintf(E->fp, "test output itracfunc 16: %d %g %g\n", i,depth[i],density[i]);
               } 
             fclose(fp);
             val=-12345;
             imax--;
           }
        if (z>=E->control.protoM_depth) {
	     z=1.0-z;
             z/=(1-E->control.protoM_depth);
             i = 2 + (int)(z*((float)(imax)));
             val =(density[i-1]*(depth[i]-z) + density[i]*(z-depth[i-1]))/ (depth[i]-depth[i-1]);
//           fprintf(E->fp, "Ttesto-- %d %g %g\n", i,z,val);
//           fprintf(E->fp, "Ttesto-- %d %g %g %g %g %g %g\n", i,z,depth[i],depth[i-1],density[i],density[i-1],val); 
           }
        else val=1.0;
   }

   /* ================================================================
    * tracfunc == 80 : molar Fe# in solid or melt (tracfunc 80)
    *   - Solid mantle (z < MO_depth at t=0): Fe# = BMO_xi_s0
    *     (= D_eff(K_Fe, xi_m)*xi_m from phi0 bulk partition).
    *   - BMO melt (z >= MO_depth): melt Fe# = BMO_xi_l (from bulk partition
    *     at phi0 in BMO_init_compositions_from_bulk; Rayleigh evolution later).
    *   Use interp_e2t=0.  The 1e-5 z tolerance matches updatetrac().
    * ================================================================ */
   else if (itracfunc == 80)
   {
       float zBOMB = E->control.MO_depth;   /* dimensionless initial front */
       if (z < zBOMB - 1.0e-5f)
           val = E->control.BMO_xi_s0;   /* solid mantle: phi0 eq. solid */
       else
           val = E->control.BMO_xi_l;    /* BMO melt: liquid composition */
   }

   /* ================================================================
    * tracfunc == 81 : H2O mass fraction c_H2O in solid/liquid (Task 3)
    *   Solid mantle: Cw_s0 = Dw*Cw_melt0 (phi0 equilibrium solid).
    *   BMO melt: Cw_l from bulk partition (BMO_init_compositions_from_bulk).
    * ================================================================ */
   else if (itracfunc == 81)
   {
       float zBOMB = E->control.MO_depth;
       if (z < zBOMB - 1.0e-5f)
           val = E->control.BMO_Cw_s0;  /* solid mantle: phi0 eq. solid */
       else
           val = E->control.BMO_Cw_l;   /* BMO melt: liquid composition */
   }

   /* ================================================================
    * tracfunc == 82 : per-tracer "solidified" flag for the BMO
    * composition tracers (80/81).  0 = not yet crystallised, 1 = solid.
    * Initialised to 0 everywhere; updatetrac() sets it to 1 the first
    * step a tracer's depth lies below the crystallisation front, after
    * which the tracer's Fe#/H2O is no longer reset and only advects.
    * Use interp_e2t=0 and interp_t2e=0 so updatetrac() owns it entirely.
    * ================================================================ */
   else if (itracfunc == 82)
   {
       val = 0.0f;
   }

   /* ================================================================
    * tracfunc == 83 : HPE mass fraction c_HPE in solid/liquid
    *   Solid mantle: Chpe_s0 = Dhpe*Chpe_melt0 (phi0 equilibrium solid).
    *   BMO melt: Chpe_l from bulk partition.
    * ================================================================ */
   else if (itracfunc == 83)
   {
       float zBOMB   = E->control.MO_depth;
       if (z < zBOMB - 1.0e-5f)
           val = E->control.BMO_Chpe_s0;  /* solid mantle: phi0 eq. solid */
       else
           val = E->control.BMO_Chpe_l;   /* BMO melt: liquid composition */
   }

   else if (itracfunc == 20) /* test: PvK97-RT instability */
   {  if (z < (0.8-0.02*cos(3.436439131*x)))
         val=0.;
      else 
         val=1.;
   }
   else if (itracfunc == 21) /* test: home-made 3Dversion of PvK97-RT */
   {  if (z < (0.8-0.02*cos(3.436439131*x)*cos(3.436439131*y)))
         val=0.;
      else 
         val=1.;
   }
   else if (itracfunc == 22) /*Experiment 22, Bellahsen*/
   {  if (x > 3.2879 && x <= 6.9091)
      {  if (x > 3.397 && z < 0.1091)
          val = 1.;
         else if (x <= 3.397 && z < 0.3030)
          val = 1.;
         else
          val = 0.;
       }
       else
         val = 0.;
   }
   else if (itracfunc == 23) /*Experiment 22, Bellahsen,3D */
   {  if (x > 3.2879 && x <= 6.9091 && y>0.2 && y < 1.8)
      {  if (x > 3.397 && z < 0.1091)
          val = 1.;
         else if (x <= 3.397 && z < 0.3030)
          val = 1.;
         else
          val = 0.;
       }
       else
         val = 0.;
   }
/* else if (itracfunc==25)
   {
     if (been_here==0)
        {
         been_here=1;
         srand48(E->parallel.me);
        }
     if (z<E->tracers.pri_depth_ini) val = 0.;
     else
        {
         if (E->tracers.bridgmanite_LM>=(drand48())) val = 1.;
         else val = 0.;
        }
   }*/
   else if (itracfunc == 30) /* subduction benchm. mantle(0)/litho(1) */
   {  if (x > 1.333 && x <= 4.0 && z >= 0.0666)
      {  if ((x> 1.4666 && z < 0.2)) 
            val = 1.;
         else if (x<= 1.4666 && z < 0.333)
            val = 1.;
         else 
            val = 0.;
      }
      else 
         val=0.;
   }
   else if (itracfunc == 31) /* subduction benchm. air(0)/solid(1) */
   {  if (z < 0.0666)
         val = 1.;
      else 
         val = 0.;
   }
   else if (itracfunc == 32) 
       /* Newsubd subduction setup 
           0  1.33 1.4666     4
           v     v v          v
           +------------------+< 0 
           |     |1 ___1______|< 0.133
           |     |_|          |< 0.266
           |  0               |
           |      0       0   |
           +------------------+< 1
        */
   {  if (x > 1.333 && x <= 4.0)
      {  if ((x> 1.4666 && z < 0.133)) 
            val = 1.;
         else if (x<= 1.4666 && z < 0.266)
            val = 1.;
         else 
            val = 0.;
      }
      else 
         val=0.;
   }
   else if (itracfunc == 33) 
        /* As 32, but in 3-D
        */
     {  if (x > 1.333 && x <= 4.0 && y > 0.5 && y < 1.5)
      {  if ((x> 1.4666 && z < 0.133)) 
            val = 1.;
         else if (x<= 1.4666 && z < 0.266)
            val = 1.;
         else 
            val = 0.;
      }
     }
   else if (itracfunc == 34)
     {
	if (z>0.5) val=1.;
	else val=0.;
     }
   else if (itracfunc == 40)
     {
      if (z>=(1.-E->tracers.thermochemical_layerd/E->data.layer_km))
       val=1.;
      else
       val=0.;
     }

     else if (itracfunc==50)
       {
	 val=(x+z)/2;
       }
     else if (itracfunc==51)
       {
	 val=0.;
	 if (z<0.55 && z>0.45) val=.019;
	 
       }
     else if (itracfunc==52)
       {
	 val=0.;
       }
   /*
     else if (itracfunc==53)
       {
	 if (z<0.22) val=0.185*(0.22-z)/0.22;
	 else val=0;
       }
     else if (itracfunc==54)
       {
	 if (z<0.22) val=0.185*(0.22-z)/0.22;
	 else val=0;
	 if (val>0.02) val=0.0;
       }
   */
     /* ridge melt depletion */
     /* 20% depl up to 30km depth then linearly decreasing to 0% at 100km depth */
     /* works only for a box with a depth of 300km */
     else if (itracfunc==56)
       {
	 
	 z*=E->data.layer_km;
	 if (z<=30.) 
	   val=0.2;
	 else if (z<=100.)
	   val=(.2/70.)*(100.-z);
	 else val=0;
       }  
 
     /* ridge melt depletion 2 */
     /* 20% depl up to 45km depth then linearly decreasing to 0% at 120km depth */
     else if (itracfunc==57)
       {
	 z*=E->data.layer_km;
	 if (z<=45.) 
	   val=0.2;
	 else if (z<=100.)
	   val=(.2/55.)*(100.-z);
	 else val=0;
       }

     /* linear depletion */
     else if (itracfunc==55)
       { z*=E->data.layer_km;
	 if (z<=150.)
	   val=(150.-z)/1500.;
	 else val=0;
       }
   /* linear depletion, depletion at solidus for halfspace cooling in melting zone */
     else if (itracfunc==58)
       {
	 if (z<=150./E->data.layer_km) 
	   {
             temp = 0.5*z/sqrt(E->control.plate_age+x/E->control.plate_vel);
             temp = erf(temp);
	     val=katz_melting(E,temp,z*factor);
	     val=max((150.-z*E->data.layer_km)/1500.,val);
	   }
	 else val=0;
       }
   /* linear depletion,depletion at solidus for constant age in melting zone */
     else if (itracfunc==59)
       { temp=150./E->data.layer_km;
	 if (z<=temp) 
	   {
             temp = 0.5*z/sqrt(E->control.plate_age);
             temp = erf(temp);
	     val=katz_melting(E,temp,z*factor,E->data.H2O);
	     val=max((150.-z*E->data.layer_km)/1500.,val);
	   }
	 else val=0;
       }
   /* ridge melt depletion,depletion at solidus for halfspace cooling in melting zone */
     else if (itracfunc==60)
       {
	 y=30./E->data.layer_km;
	 if (z<=y)
	   {
             temp = 0.5*z/sqrt(E->control.plate_age+x/E->control.plate_vel);
             temp = erf(temp);
	     val=katz_melting(E,temp,z*factor,E->data.H2O);
	     val=max(0.2,val);
           }
	 else if (z<=100./E->data.layer_km) 
	   {
             temp = 0.5*z/sqrt(E->control.plate_age);
             temp = erf(temp);
	     val=katz_melting(E,temp,z*factor,E->data.H2O);
	     val=max(0.002857*(100.-z*E->data.layer_km),val);
           }
	 else val=0;
       }
   /* ridge melt depletion,depletion at solidus for constant age in melting zone */
     else if (itracfunc==61)
       {
	 if (z<=30./E->data.layer_km)
	   {
             temp = 0.5*z/sqrt(E->control.plate_age);
             temp = erf(temp);
	     val=katz_melting(E,temp,z*factor,E->data.H2O);	     
	     val=max(0.2,val);
           }
	 else if (z<=100./E->data.layer_km) 
	   {
             temp = 0.5*z/sqrt(E->control.plate_age);
             temp = erf(temp);
	     val=katz_melting(E,temp,z*factor,E->data.H2O);
	     val=max(0.002857*(100.-z*E->data.layer_km),val);
           }
	 else val=0;
       }
     else if (itracfunc==62)
       {
	 /* Katz ridge melt depletion for 1350K and wtH2O=0.02% */
	 z*=factor;
	 if (z<=0.6) 
	   val=0.21-z*0.025/0.6;
	 else if (z<=3.3)
	   {val=max((.185/2.1)*(2.7-z),(.025/.8)*(3.3-z));
	   val=max(val,0.0575*(3.-z));}
	 else val=0.;
       }  
     else if (itracfunc==63)
       {
	 // irregular melt inflow at about 1%
       }  
     else if (itracfunc==64)
       {
         // ridge melt depletion for age=8 Ma
         if ((z*factor)<5.)
           {
             temp = 0.5*z/sqrt(8.0*1.0e6*3.155734e7/time_scale);
             temp = erf(temp);
             val=katz_melting(E,temp,z*factor,E->data.H2O);
             if (z<=55./E->data.layer_km)
               val=max(val,0.15);
           }
         else val=0.;
       }
	else if (itracfunc==65)
       {
	 // ridge melt depletion for age=0 Ma
	 if ((z*factor)<5.) 
	   {
             temp=1.;
	     val=katz_melting(E,temp,z*factor,E->data.H2O);
	   }
	 else val=0.;
       }  
     else if (itracfunc==70)
       {
	 // ridge melt depletion loaded from file
	 if (visits==0) 
	   { sprintf(filename,"%s.ridge_melt_depl.input",E->control.data_file);
             visits=1;
	     fp = fopen(filename,"r");
             if (fp==NULL)
	       tracexit(E,"tracfunc=70: Unable to load A.ridge_melt_depl.input: Exit");
	     for (i=0; i<51; i++)
	       { fgets(input_s,200,fp);
		 sscanf(input_s,"%g",&val);
		 rmd[i] = val;
	       }
	     rmd[51]=val;
             fclose(fp);
	   }
	 z*=50.;
	 i=(int)(z);
	 z-=(float)(i);
	 val=(1-z)*rmd[i] + z*rmd[i+1];
       }  
     else if (itracfunc>100)
       {
	 // ridge melt depletion for age=x Ma
	 if ((z*factor)<5.) 
	   {
             temp = 0.5*z/sqrt((float)(itracfunc-100)*1.0e6*3.155734e7/time_scale);
             temp = erf(temp);
	     val=katz_melting(E,temp,z*factor,E->data.H2O);
	   }
	 else val=0.;
       }  
      else 
       tracexit(E,"Error in trac_funcinit: itracfunc not defined\n: STOP");
   return(val);
}
void tracexit(E,string)
   struct All_variables *E;
   char *string;
{  
   /* This routine ends the program if anything goes wrong with
    * the tracers. It attempts to do that in a clean way by first 
    * stopping all parallel stuff.
    */

   void parallel_process_termination();

   fflush(E->fp);
   fprintf(E->fp, "Program is exited with following messsage:\n");
   fprintf(E->fp, "%s\n", string);

   if (E->parallel.me == 0) 
   {  fprintf(stderr, "Program is exited with following messsage:\n");
      fprintf(stderr, "%s\n", string);
   }

   parallel_process_termination();

   exit(1);
}

/* ======================================================================
 * realloc_trac_arrays: expand all tracer arrays when larraylen is
 * about to be exceeded.  new_len must be > E->tracers.larraylen.
 * All arrays that were allocated with larraylen as their dimension
 * are reallocated here.  Newly added slots are zero-initialised.
 * ====================================================================== */
void realloc_trac_arrays(E, new_len)
   struct All_variables *E;
   int new_len;
{
   int ifunc, n, old_len;
   int ndim = E->mesh.nsd;

   old_len = E->tracers.larraylen;
   if (new_len <= old_len) return;

   fprintf(E->fp,
           "realloc_trac_arrays: expanding larraylen %d -> %d (proc %d)\n",
           old_len, new_len, E->parallel.me);
   fprintf(stderr,
           "realloc_trac_arrays: expanding larraylen %d -> %d (proc %d)\n",
           old_len, new_len, E->parallel.me);

#define REALLOC_FLOAT(ptr, n) \
   do { \
      (ptr) = (float *) realloc((ptr), (n) * sizeof(float)); \
      if (!(ptr)) { tracexit(E, "realloc_trac_arrays: out of memory"); } \
      memset((ptr) + old_len, 0, ((n) - old_len) * sizeof(float)); \
   } while(0)

#define REALLOC_INT(ptr, n) \
   do { \
      (ptr) = (int *) realloc((ptr), (n) * sizeof(int)); \
      if (!(ptr)) { tracexit(E, "realloc_trac_arrays: out of memory"); } \
      memset((ptr) + old_len, 0, ((n) - old_len) * sizeof(int)); \
   } while(0)

   /* coordinate arrays */
   for (n = 1; n <= ndim; n++)
   {  REALLOC_FLOAT(E->tracers.xtrac[n],  new_len);
      REALLOC_FLOAT(E->tracers.xtrac1[n], new_len);
      REALLOC_FLOAT(E->tracers.xtrac2[n], new_len);
   }

   /* velocity arrays */
   for (n = 1; n <= ndim; n++)
   {  REALLOC_FLOAT(E->tracers.vtrac[n],  new_len);
      REALLOC_FLOAT(E->tracers.vtrac1[n], new_len);
   }

   /* element index arrays */
   REALLOC_INT(E->tracers.elemnr,  new_len);
   REALLOC_INT(E->tracers.elemnr1, new_len);

   /* tracer-per-element lookup tables (allocated as len+1) */
   REALLOC_INT(E->tracers.tracs_in_elem,  new_len + 1);
   REALLOC_INT(E->tracers.tracs_in_elem1, new_len + 1);

   /* tracer function-value arrays (allocated as len+1) */
   for (ifunc = 0; ifunc < E->tracers.ntracfunc; ifunc++)
   {  REALLOC_FLOAT(E->tracers.ftrac[ifunc],  new_len + 1);
      REALLOC_FLOAT(E->tracers.ftrac2[ifunc], new_len + 1);
   }

   /* crossborder[itrac] is indexed by tracer id (0 .. ntrac-1) */
   if (E->parallel.nproc > 1 && E->tracers.crossborder != NULL)
   {  int *cb = (int *) realloc(E->tracers.crossborder, new_len * sizeof(int));
      if (!cb) { tracexit(E, "realloc_trac_arrays: crossborder OOM"); }
      memset(cb + old_len, 0, (new_len - old_len) * sizeof(int));
      E->tracers.crossborder = cb;
   }

#undef REALLOC_FLOAT
#undef REALLOC_INT

   /* Scale MPI send/recv buffers with larraylen (same ratios as parallel_tracer_setup). */
   if (E->parallel.nproc > 1 && old_len > 0)
   {  int ineighb, nneighb_max, new_sr;
      ndim = E->mesh.nsd;
      nneighb_max = (ndim == 2) ? 8 : 26;
      for (ineighb = 1; ineighb <= nneighb_max; ineighb++)
      {  if (E->tracers.sr[ineighb].used != 1) continue;
         new_sr = (int)((long)new_len * (long)E->tracers.sr[ineighb].size
                        / (long)old_len);
         if (new_sr < E->tracers.sr[ineighb].size)
            new_sr = E->tracers.sr[ineighb].size;
         grow_sr(E, ineighb, new_sr);
      }
   }

   E->tracers.larraylen = new_len;
}
