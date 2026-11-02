/*****************************************
 *   CC  III  TTTTT   CC   OO   MM MM    *
 *  C     I     T    C    O  O  M M M    *
 *  C     I     T    C    O  O  M   M    *
 *   CC  III    T     CC   OO   M   M    *
 *                                       *  
 * Developed at CIT for COnvection in    *
 * the Mantle by Louis Moresi 1992-today *
 *                                       *
 * You are free to use this code but it  * 
 * is distrubuted as BeWare i.e. it does *
 * not carry any guarantees or warranties *
 * of reliability.                       *
 *                                       *
 * Please respect all the time and work  *
 * that went into the development of the *
 * code.                                 *  
 *                                       *
 * LM                                    *
 *****************************************/
/* Heat-transport solver using Petrov-Galerkin streamline-upwind methods (after Brooks' thesis). */

#include <stdio.h>
#include <malloc.h>
#include <sys/types.h>
#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"
#include <mpi.h>

extern int Emergency_stop;

struct el { double gpt[9]; };

/* ============================================
   Generic adv-diffusion for temperature field.
   ============================================ */


void advection_diffusion_parameters(E)
     struct All_variables *E;

{

    void std_timestep();
    /* Set intial values, defaults & read parameters*/
    
    E->advection.temp_iterations = 2; /* petrov-galerkin iterations: minimum value. */
    E->advection.total_timesteps = 1; 
    E->advection.sub_iterations = 1;
    E->advection.last_sub_iterations = 1;
    E->advection.gamma = 0.5;
    E->advection.ADVECTION = 1;
    E->advection.dtmult = 1.; 
 
    input_boolean("ADV",&(E->advection.ADVECTION),"on");
   
    input_int("minstep",&(E->advection.min_timesteps),"1");
    input_int("maxstep",&(E->advection.max_timesteps),"1000");
    input_int("maxtotstep",&(E->advection.max_total_timesteps),"1000000");
    input_float("finetunedt",&(E->advection.fine_tune_dt),"0.9");
    input_float("fixed_timestep",&(E->advection.fixed_timestep),"0.0");
    input_int("adv_sub_iterations",&(E->advection.temp_iterations),"2,2,nomax");
    input_float("maxadvtime",&(E->advection.max_dimensionless_time),"10.0");
   
    input_float("sub_tolerance",&(E->advection.vel_substep_aggression),"0.005");  
    input_int("maxsub",&(E->advection.max_substeps),"25");
  
    input_float("liddefvel",&(E->advection.lid_defining_velocity),"0.01");
    input_float("sublayerfrac",&(E->advection.sub_layer_sample_level),"0.5");            
    input_float("timestepmax",&(E->control.timestepmax),"1.0e8");
    
    /* Multiplication of total timestep with a constant */
    input_float("dtmult",&(E->advection.dtmult),"1.0");            
	      
    /* inittreduct: dt scale factor ramping linearly from inittreduct to 1.0 over the
     * first n_early_steps steps (1.0 or n_early_steps=0 disables). */
    input_float("inittreduct",&(E->control.inittreduct),"0.1");
    input_int("n_early_steps",&(E->control.n_early_steps),"50");

    /* CFL_number: dt_adv <= CFL * min_element_size / |v| per element (0 disables). */
    input_float("CFL_number",&(E->control.CFL_number),"0.5");

    /* max_dT_frac: reserved (1.0 = disabled). */
    input_float("max_dT_frac",&(E->control.max_dT_frac),"1.0");
  /* allocate memory */

  return;
}

void advection_diffusion_allocate_memory(E)
     struct All_variables *E;

{ int i;

  E->Tdot= (float *)malloc((E->lmesh.nno+1)*sizeof(float));
    
  for(i=1;i<=E->lmesh.nno;i++)
    E->Tdot[i]=0.0;
return;
}


void PG_timestep(E)
     struct All_variables *E;
{    
    void timestep();
    //void time_display();
    void predictor();
    void corrector();
    void pg_solver();
    void remove_horiz_ave();
    void std_timestep();
    void temperatures_conform_bcs();
    void general_stokes_solver();
    void tracadv();
    void tracinterp_e2t();
    void updatetrac();    
    void tracinterp_t2e();
    void BMO_zone_temperature_override();
    float global_tdot();
    void vcopy();
    void strain_rate_2_inv();
    void apply_node_limiter();

    int get_v_estimate();
    int i,itrac,j,psc_pass,count,steps,ielem;
    float **ftrac_pointer;
    float **xtrac_pointer;
    float *DTdot;
    float temp;
    double eta_ave;
    int el,after8h_in;
    static int record_every,tracer_storage,after8h;
    const int vpts=vpoints[E->mesh.nsd];

    static int loops_since_new_eta = 0;
    static int been_here = 0;
    //static double *p_init_time;
    double CPU_time0();
   
    DTdot= (float *)malloc((E->lmesh.nno+1)*sizeof(float));


    if (been_here==0)    
     {  after8h=0;
       //p_init_time=&E->control.initial_time;
 	E->advection.timesteps=0;
	E->data.DioverRa=E->data.therm_exp*E->data.grav_acc*
	1.e3*E->data.layer_km/(E->control.Atemp*E->data.Cp);

      fprintf(E->fp,"Viscous dissipation is");
      if (E->control.shear_heating ==0) 
	fprintf(E->fp," not");
      fprintf(E->fp," activated\n");
     }
    been_here++;
    E->advection.timesteps++;
 
    std_timestep(E);

  
    // automatic output after 40 hrs. and 80 hrs.: Reset all
    if ((after8h==1 || after8h==3))
      { E->control.record_every=record_every;
	E->tracers.storage=tracer_storage;
        after8h++;
      }
    // Prepare to do automatic output after 40 hrs. and 80 hrs.
    if (after8h == 0 && E->parallel.me==0)
      {
	if ((CPU_time0()-E->control.initial_time)>144000.)
	  after8h++;  
      }
    else if (after8h == 2 && E->parallel.me==0)
      {
	if ((CPU_time0()-E->control.initial_time)>288000.)
	  after8h++;
      }

    MPI_Allreduce(&after8h, &after8h_in,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    after8h=after8h_in;

    if (after8h==1 || after8h==3)
      {	    record_every=E->control.record_every;
	    tracer_storage=E->tracers.storage;
	    E->control.record_every=E->monitor.solution_cycles;
	    E->tracers.storage=E->monitor.solution_cycles;     }


    //if (been_here==1) E->advection.timestep=.125e-6;

    /* In case of flux-controller: re-calculate plume-radius */
    if (E->control.fluxcontrol==1 &&
        E->advection.timesteps%E->control.updatefluxevery==0 &&
        E->advection.timesteps > 1)
    {  /* Update radius */
       plume_buoyancy_flux(E);
       if (E->data.buoy_flux < -1)
       E->segment.plume_radius *=
              (float)pow((double)(E->control.flux/-E->data.buoy_flux), 
                        (double)(E->control.relaxfluxupdate));
       if (E->parallel.me == 0)
          fprintf(stderr, "Flux=%g, Changing plume radius to %g\n", 
                  E->data.buoy_flux, E->segment.plume_radius);
       E->control.timdepTbc=1;
    }

    /* In case of timedep. T bc.: recall Tbc-routine */
    if (E->control.timdepTbc == 1) 
       temperature_boundary_conditions(E);

    count = 0;

    /* adjust latent heat of melt from updatetrac to length of tstep */
    if (E->data.heat_of_melt != 0.)
     for (ielem=1;ielem<=E->lmesh.nel;ielem++)
       E->Qmelt[ielem]/=E->advection.timestep;
 
      if (E->mesh.nsd==2)
       for (i=1;i<=E->lmesh.nnov;i++)
       {
	E->Vstar[1][i]=E->V[1][i];
	E->Vstar[2][i]=E->V[2][i];
       }
      else if (E->mesh.nsd==3)
       for (i=1;i<=E->lmesh.nnov;i++)
       {
	E->Vstar[1][i]=E->V[1][i];
	E->Vstar[2][i]=E->V[2][i];
	E->Vstar[3][i]=E->V[3][i];
       }


      /* viscous dissipation calculation */
    if (E->control.shear_heating > 0)
     {strain_rate_2_inv(E,E->eedot,0);
        /*if (E->viscosity.edot31in2Davail == 3)
          E->viscosity.calcstress=0; */
      for (el=1;el<=E->lmesh.nel;el++)
       {
        eta_ave=E->EVI[E->mesh.levmax][(el-1)*vpts+1];
        for(i=2;i<=vpts;i++)  
	  eta_ave*=E->EVI[E->mesh.levmax][(el-1)*vpts+i];
        eta_ave=(pow(eta_ave,1.0/(float)(vpts)));
        //fprintf(E->fp,"el %i, eps %g, eta %g\n",el,E->eedot[el],eta_ave);
        //eta_ave/=(float)(vpts);
       
        if (eta_ave<1000.)//(E->viscosity.max_value*.0001))
         E->visc_diss[el]=E->data.DioverRa*E->eedot[el]*eta_ave;
        else 
	 E->visc_diss[el]=0.;
        //fprintf(E->fp,"%i %g\n",el,E->visc_diss[el]);
       }
      }
    else if (been_here==1)
      for (el=1;el<=E->lmesh.nel;el++)
        E->visc_diss[el]=0.;

    if (E->advection.ADVECTION)
     {  
      predictor(E,E->T,E->Tdot);
      for(psc_pass=0;psc_pass<E->advection.temp_iterations; psc_pass++)
       { 
	 pg_solver(E,E->T,E->Tdot,DTdot,E->V,
              E->convection.heat_sources,1.0,1,E->TB,E->node);
	 corrector(E,E->T,E->Tdot,DTdot);

   //apply_node_limiter(E, E->T);   /* #### clamp nodal T to local neighbor range  */
       }
      temperatures_conform_bcs(E,E->T);
      /* Keep the BMO melt zone isothermal (T = MO_temperature for z >= MO_depth_initial);
       * only when BMO_T_override_mode == 1. */
      BMO_zone_temperature_override(E);
      thermal_buoyancy(E);
     }
    if (E->tracers.tracers_used!=0)
     {
      /* safe tracers for 2D */
      E->tracers.lntrac_stor=E->tracers.lntrac;
      if (E->mesh.nsd==2)
       for (itrac=0; itrac<E->tracers.lntrac; itrac++)
	{
	  E->tracers.xtrac2[1][itrac]=E->tracers.xtrac[1][itrac];
	  E->tracers.xtrac2[2][itrac]=E->tracers.xtrac[2][itrac];
	  for (i=0;i<E->tracers.ntracfunc;i++)
            E->tracers.ftrac2[i][itrac]=E->tracers.ftrac[i][itrac];
	}
      /* safe tracers for 3D */
      else if (E->mesh.nsd==3)
       for (itrac=0; itrac<E->tracers.lntrac; itrac++)
	{
	  E->tracers.xtrac2[1][itrac]=E->tracers.xtrac[1][itrac];
	  E->tracers.xtrac2[2][itrac]=E->tracers.xtrac[2][itrac];
	  E->tracers.xtrac2[3][itrac]=E->tracers.xtrac[3][itrac];
	  for (i=0;i<E->tracers.ntracfunc;i++)
            E->tracers.ftrac2[i][itrac]=E->tracers.ftrac[i][itrac];
	}
      /* swap pointers */
      ftrac_pointer=E->tracers.ftrac;
      E->tracers.ftrac=E->tracers.ftrac2;
      xtrac_pointer=E->tracers.xtrac;
      E->tracers.xtrac=E->tracers.xtrac2;


      /* update tracers */
      //time_display(E,p_init_time, 23);
      tracadv(E,0);
      //time_display(E,p_init_time, 24);
      tracinterp_e2t(E);
      //time_display(E,p_init_time, 25);
      updatetrac(E,0);
      //time_display(E,p_init_time, 26);    
      tracinterp_t2e(E);
      //if (E->data.heat_of_melt != 0.)
      // for (ielem=1;ielem<=E->lmesh.nel;ielem++)
      //  E->Qmelt[ielem]/=E->advection.timestep;
     }

    //time_display(E,p_init_time, 27);
    general_stokes_solver(E);
    //time_display(E,p_init_time, 28);
      /* reload tracers for 2D */
    if (E->tracers.tracers_used!=0)
     {
      E->tracers.lntrac=E->tracers.lntrac_stor;
      E->tracers.lntrac1=E->tracers.lntrac_stor;
      E->tracers.ftrac=ftrac_pointer;
      E->tracers.xtrac=xtrac_pointer;
     }

    if (E->advection.ADVECTION)
     {  
      for (i=1;i<=E->lmesh.nno;i++)
	{
	  E->T[i]=E->Tstar[i];
      	  E->Tdot[i]=0.;
	}
      
      for(psc_pass=0;psc_pass<E->advection.temp_iterations; psc_pass++)
       { 
	 pg_solver(E,E->T,E->Tdot,DTdot,E->V,
              E->convection.heat_sources,1.0,1,E->TB,E->node);
	 corrector(E,E->T,E->Tdot,DTdot);
       }     
     }

     if (E->mesh.nsd==2)
       for (i=1;i<=E->lmesh.nno;i++)
	{
	  E->V[1][i]+=E->Vstar[1][i];
	  E->V[2][i]+=E->Vstar[2][i];
  	  E->V[1][i]/=2.;
	  E->V[2][i]/=2.;
	  }
     else if (E->mesh.nsd==3)
       for (i=1;i<=E->lmesh.nno;i++)
	{
	  E->V[1][i]+=E->Vstar[1][i];
	  E->V[2][i]+=E->Vstar[2][i];
  	  E->V[1][i]/=2.;
	  E->V[2][i]/=2.;
	  E->V[3][i]+=E->Vstar[3][i];
	  E->V[3][i]/=2.;
	}

    E->monitor.elapsed_time += E->advection.timestep; 
    E->advection.total_timesteps++;
    if (been_here==1 && E->control.restart==0)
      E->monitor.elapsed_time /= 2.;  

    count++; 

    //apply_node_limiter(E, E->T);   /* #### clamp nodal T to local neighbor range  */

    temperatures_conform_bcs(E,E->T);
    BMO_zone_temperature_override(E);  /* isothermal BMO zone if mode==1 */
    thermal_buoyancy(E);
 
    E->advection.last_sub_iterations = count;
      
    if(E->monitor.solution_cycles < E->advection.max_timesteps)
        E->control.keep_going = 1;
    else
 	E->control.keep_going = 0;
    
    free((void *) DTdot );    /* free memory for vel solver */
    
return;  
}


/* ==============================
   predictor and corrector steps.
   ============================== */

void predictor(E,field,fielddot)
     struct All_variables *E;
     float *field,*fielddot;

{ 
    int node;
    float multiplier;

   multiplier = (1.0-E->advection.gamma) * E->advection.timestep;
      
    for(node=1;node<=E->lmesh.nno;node++)  {
	if(!(E->node[node] & (TBX | TBZ | TBY))) 
	    field[node] += multiplier * fielddot[node] ;
	E->Tstar[node]=field[node];  //save this Temperature for 2nd Iteration
	fielddot[node] = 0.0;
    }

   return; }

void corrector(E,field,fielddot,Dfielddot)
     struct All_variables *E;
     float *field,*fielddot,*Dfielddot;
    
{  int node;
   float multiplier;

   multiplier = E->advection.gamma * E->advection.timestep;
   for(node=1;node<=E->lmesh.nno;node++) {
       if(!(E->node[node] & (TBX | TBZ | TBY)))
   	   field[node] += multiplier * Dfielddot[node];
       fielddot[node] +=  Dfielddot[node];
   }
   
   return;  
 }

/* ===================================================
   pg_solver: advective-diffusive residual -> delta Tdot
   =================================================== */


void pg_solver(E,T,Tdot,DTdot,V,Q0,diff,bc,TBC,FLAGS)
     struct All_variables *E;
     float *T,*Tdot,*DTdot;
     float **V;
     struct SOURCES Q0;
     float diff;
     int bc;
     float **TBC;
     unsigned int *FLAGS;
{
    void get_global_shape_fn();
    void pg_shape_fn();
    void element_residual();
    void e_exchange_node_fc();
    void exchange_node_f20();
    int el,e,a,i,a1;
    double Eres[9];  /* correction to the (scalar) Tdot field */

    struct Shape_function PG;
 
    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int ends=enodes[dims];
  
    for(i=1;i<=E->lmesh.nno;i++)
 	  DTdot[i] = 0.0;

     
    for(el=1;el<=E->lmesh.nel;el++)    {

	  pg_shape_fn(E,el,&PG,V,diff);
	  element_residual(E,el,PG,V,T,Tdot,Q0,Eres,diff,TBC,FLAGS);

      for(a=1;a<=ends;a++) {
	    a1 = E->ien[el].node[a];
	    DTdot[a1] += Eres[a]; 
        }

      } /* next element */
   
      
   exchange_node_f20(E,DTdot,E->mesh.levmax);
   
    for(i=1;i<=E->lmesh.nno;i++) {
	  DTdot[i] *= E->Mass[i];         /* lumped mass matrix */
      
      }
    
    return;    
}



/* ===================================================
   Petrov-Galerkin shape functions for a given element
   =================================================== */

void pg_shape_fn(E,el,PG,V,diffusion)
     struct All_variables *E;
     int el;
     struct Shape_function *PG;
     float **V;
     float diffusion;

{ 
    int i,j;
    int *ienmatrix;

    double xsi1,xsi2,xsi3;
    double uc1,uc2,uc3;
    double aa1,aa2,aa3;
    double ah1,ah2,ah3;
    double size1,size2,size3;
    double neg1,neg2,neg3;
    double u1,u2,u3;
    double adiff,dx1,dx2,dx3,uxse,ueta,ufai,xse,eta,fai;

    double twodiff,prod,prod1;
  
    double K,unorm;
    
    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int ends=enodes[E->mesh.nsd];
    const int vpts=vpoints[E->mesh.nsd];
  
    ienmatrix=E->ien[el].node;

    twodiff = 2.0*diffusion;
 
   
    size1 = (double) E->eco[el].size[1]; 
    size2 = (double) E->eco[el].size[2]; 
    size3 = (double) E->eco[el].size[3]; 

    uc1 =  uc2 = uc3 = 0.0;
         
    if(3==dims) {
	for(i=1;i<=ENODES3D;i++) {
	    uc1 +=  E->N.ppt[GNPINDEX(i,1)]*V[1][ienmatrix[i]];
	    uc2 +=  E->N.ppt[GNPINDEX(i,1)]*V[2][ienmatrix[i]];
	    uc3 +=  E->N.ppt[GNPINDEX(i,1)]*V[3][ienmatrix[i]];
	}

    dx1 = 0.25*(E->X[1][ienmatrix[3]]+E->X[1][ienmatrix[4]]
               +E->X[1][ienmatrix[7]]+E->X[1][ienmatrix[8]]
               -E->X[1][ienmatrix[1]]-E->X[1][ienmatrix[2]]
               -E->X[1][ienmatrix[5]]-E->X[1][ienmatrix[6]]);
    dx2 = 0.25*(E->X[2][ienmatrix[3]]+E->X[2][ienmatrix[4]]
               +E->X[2][ienmatrix[7]]+E->X[2][ienmatrix[8]]
               -E->X[2][ienmatrix[1]]-E->X[2][ienmatrix[2]]
               -E->X[2][ienmatrix[5]]-E->X[2][ienmatrix[6]]);
    dx3 = 0.25*(E->X[3][ienmatrix[3]]+E->X[3][ienmatrix[4]]
               +E->X[3][ienmatrix[7]]+E->X[3][ienmatrix[8]]
               -E->X[3][ienmatrix[1]]-E->X[3][ienmatrix[2]]
               -E->X[3][ienmatrix[5]]-E->X[3][ienmatrix[6]]);
    uxse = fabs(uc1*dx1+uc2*dx2+uc3*dx3);

    dx1 = 0.25*(E->X[1][ienmatrix[2]]+E->X[1][ienmatrix[3]]
               +E->X[1][ienmatrix[6]]+E->X[1][ienmatrix[7]]
               -E->X[1][ienmatrix[1]]-E->X[1][ienmatrix[4]]
               -E->X[1][ienmatrix[5]]-E->X[1][ienmatrix[8]]);
    dx2 = 0.25*(E->X[2][ienmatrix[2]]+E->X[2][ienmatrix[3]]
               +E->X[2][ienmatrix[6]]+E->X[2][ienmatrix[7]]
               -E->X[2][ienmatrix[1]]-E->X[2][ienmatrix[4]]
               -E->X[2][ienmatrix[5]]-E->X[2][ienmatrix[8]]);
    dx3 = 0.25*(E->X[3][ienmatrix[2]]+E->X[3][ienmatrix[3]]
               +E->X[3][ienmatrix[6]]+E->X[3][ienmatrix[7]]
               -E->X[3][ienmatrix[1]]-E->X[3][ienmatrix[4]]
               -E->X[3][ienmatrix[5]]-E->X[3][ienmatrix[8]]);
    ueta = fabs(uc1*dx1+uc2*dx2+uc3*dx3);

    dx1 = 0.25*(E->X[1][ienmatrix[5]]+E->X[1][ienmatrix[6]]
               +E->X[1][ienmatrix[7]]+E->X[1][ienmatrix[8]]
               -E->X[1][ienmatrix[1]]-E->X[1][ienmatrix[2]]
               -E->X[1][ienmatrix[3]]-E->X[1][ienmatrix[4]]);
    dx2 = 0.25*(E->X[2][ienmatrix[5]]+E->X[2][ienmatrix[6]]
               +E->X[2][ienmatrix[7]]+E->X[2][ienmatrix[8]]
               -E->X[2][ienmatrix[1]]-E->X[2][ienmatrix[2]]
               -E->X[2][ienmatrix[3]]-E->X[2][ienmatrix[4]]);
    dx3 = 0.25*(E->X[3][ienmatrix[5]]+E->X[3][ienmatrix[6]]
               +E->X[3][ienmatrix[7]]+E->X[3][ienmatrix[8]]
               -E->X[3][ienmatrix[1]]-E->X[3][ienmatrix[2]]
               -E->X[3][ienmatrix[3]]-E->X[3][ienmatrix[4]]);
    ufai = fabs(uc1*dx1+uc2*dx2+uc3*dx3);

    xse = (uxse>twodiff)? (1.0-twodiff/uxse):0.0;
    eta = (ueta>twodiff)? (1.0-twodiff/ueta):0.0;
    fai = (ufai>twodiff)? (1.0-twodiff/ufai):0.0;

    unorm = uc1*uc1 + uc2*uc2 + uc3*uc3;

    adiff = (unorm>0.000001)?( (uxse*xse+ueta*eta+ufai*fai)/(2.0*unorm) ):0.0;

        for(i=1;i<=VPOINTS3D;i++) {
          u1 = u2 = u3 = 0.0;
          for(j=1;j<=ENODES3D;j++)  /* this line heavily used */ {
                u1 += V[1][ienmatrix[j]] * E->N.vpt[GNVINDEX(j,i)];
                u2 += V[2][ienmatrix[j]] * E->N.vpt[GNVINDEX(j,i)];
                u3 += V[3][ienmatrix[j]] * E->N.vpt[GNVINDEX(j,i)];
            }

          for(j=1;j<=ENODES3D;j++) {
        prod1 = (u1 * E->gNX[el].vpt[GNVXINDEX(0,j,i)] +
                 u2 * E->gNX[el].vpt[GNVXINDEX(1,j,i)] +
                 u3 * E->gNX[el].vpt[GNVXINDEX(2,j,i)] ) ;

                PG->vpt[GNVINDEX(j,i)] = E->N.vpt[GNVINDEX(j,i)] + adiff * prod1;
            }
          }
    }

  else if (2==dims)  {
        for(i=1;i<=ENODES2D;i++) {
	    uc1 +=  E->N.ppt[GNPINDEX(i,1)]*V[1][ienmatrix[i]];
	    uc2 +=  E->N.ppt[GNPINDEX(i,1)]*V[2][ienmatrix[i]];
          }
    dx1 = 0.5*(E->X[1][ienmatrix[3]]+E->X[1][ienmatrix[4]]
              -E->X[1][ienmatrix[1]]-E->X[1][ienmatrix[2]]);
    dx2 = 0.5*(E->X[2][ienmatrix[3]]+E->X[2][ienmatrix[4]]
              -E->X[2][ienmatrix[1]]-E->X[2][ienmatrix[2]]);
    uxse = fabs(uc1*dx1+uc2*dx2);

    dx1 = 0.5*(E->X[1][ienmatrix[2]]+E->X[1][ienmatrix[3]]
              -E->X[1][ienmatrix[1]]-E->X[1][ienmatrix[4]]);
    dx2 = 0.5*(E->X[2][ienmatrix[2]]+E->X[2][ienmatrix[3]]
              -E->X[2][ienmatrix[1]]-E->X[2][ienmatrix[4]]);
    ueta = fabs(uc1*dx1+uc2*dx2);

    xse = (uxse>twodiff)? (1.0-twodiff/uxse):0.0;
    eta = (ueta>twodiff)? (1.0-twodiff/ueta):0.0;

    unorm = uc1*uc1 + uc2*uc2;

    adiff = (unorm>0.000001)?( (uxse*xse+ueta*eta)/(2.0*unorm) ):0.0;

        for(i=1;i<=VPOINTS2D;i++) {
          u1 = u2 = 0.0;
          for(j=1;j<=ENODES2D;j++)  /* this line heavily used */ {
                u1 += V[1][ienmatrix[j]] * E->N.vpt[GNVINDEX(j,i)];
                u2 += V[2][ienmatrix[j]] * E->N.vpt[GNVINDEX(j,i)];
            }

          for(j=1;j<=ENODES2D;j++) {

        prod1 = (u1 * E->gNX[el].vpt[GNVXINDEX(0,j,i)] +
                 u2 * E->gNX[el].vpt[GNVXINDEX(1,j,i)]);
               
                PG->vpt[GNVINDEX(j,i)] = E->N.vpt[GNVINDEX(j,i)] + adiff * prod1;
            }
          }
    }

   return;
 }



/* ==========================================
   Residual force vector from heat-transport.
   Used to correct the Tdot term.
   =========================================  */

void element_residual(E,el,PG,vel,field,fielddot,Q0,Eres,diff,BC,FLAGS)
     struct All_variables *E;
     int el;
     struct Shape_function PG;
     float **vel;
     float *field,*fielddot;
     struct SOURCES Q0;
     double Eres[9];
     float diff;
     float **BC;
     unsigned int *FLAGS;

{
    int i,j,a,k,node,nodes[4],d,aid,back_front,onedfns;
    double Q;
    double adiab=E->data.adiabatic_gr*E->data.layer_km/E->data.ref_temperature;
    double dT[9];
    double tx1[9],tx2[9],tx3[9];
    double v1[9],v2[9],v3[9];
    double adv_dT,t2[4];
    double T,DT;
 
    register double prod,sfn;
    struct Shape_function1 GM;
    struct Shape_function1_dA dGamma;
    double temp;
    float visc_diss;
    float U_Th__withdecay,Q0comp_excess;
    float decay = E->tracers.Q0compo_value * exp(-E->tracers.lambda*E->monitor.elapsed_time); 
    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int ends=enodes[dims];
    const int vpts=vpoints[dims];
    const int diffusion = (diff != 0.0);
    const float nondim = 1.e6*E->data.layer_km*E->data.layer_km/(3.1557e7*E->data.therm_diff*E->data.Cp*E->data.ref_temperature);
    const float dimens_time = 1.e6*E->data.layer_km*E->data.layer_km/(3.1557e7*E->data.therm_diff); // years
    void get_global_1d_shape_fn();
    
    for(i=1;i<=vpts;i++)	{ 
	dT[i]=0.0;
	v1[i] = tx1[i]=  0.0;
	v2[i] = tx2[i]=  0.0;
	v3[i] = tx3[i]=  0.0;
    }
  
    U_Th__withdecay = nondim*(1.097738e-11*exp(-9.845841e-10*E->monitor.elapsed_time*dimens_time)   // U235
                           + 5.900101-12 * exp(-1.550665e-10*E->monitor.elapsed_time*dimens_time)   // U238              
			   + 4.10072e-12 * exp(-4.951051e-11*E->monitor.elapsed_time*dimens_time)); // Th

    for(j=1;j<=ends;j++)       {
	node = E->ien[el].node[j];
	T = field[node];  
	if(E->node[node] & (TBX | TBY | TBZ))
	    DT=0.0;
	else
	    DT = fielddot[node];
	
	if(3==dims) {
	     for(i=1;i<=vpts;i++)  {
		dT[i] += DT * E->N.vpt[GNVINDEX(j,i)];
		tx1[i] += E->gNX[el].vpt[GNVXINDEX(0,j,i)] * T; 
		tx2[i] += E->gNX[el].vpt[GNVXINDEX(1,j,i)] * T; 
		tx3[i] += E->gNX[el].vpt[GNVXINDEX(2,j,i)] * T; 
		sfn = E->N.vpt[GNVINDEX(j,i)];
		v1[i] += vel[1][node] * sfn;
		v2[i] += vel[2][node] * sfn;
		v3[i] += vel[3][node] * sfn;     
	     }
	}
	else {
	    for(i=1;i<=vpts;i++)  {
		dT[i] += DT * E->N.vpt[GNVINDEX(j,i)];
		tx1[i] +=  E->gNX[el].vpt[GNVXINDEX(0,j,i)] * T; 
		tx2[i] +=  E->gNX[el].vpt[GNVXINDEX(1,j,i)] * T;   
		sfn = E->N.vpt[GNVINDEX(j,i)];
		v1[i] += vel[1][node] * sfn;
		v2[i] += vel[2][node] * sfn;   } 
	}
    }
    if (E->advection.timesteps==0)
      Q0comp_excess = 0.;
    else
      {
       if (E->tracers.Q0compo_tracnr >= 0)
        {
          Q0comp_excess = E->tracers.felem[E->tracers.Q0compo_tracnr][el] * decay;
        }
       else Q0comp_excess = 0.; 
       if (E->tracers.Q0compo_U_Th_tracnr>=0)
	{
	  Q0comp_excess+= E->tracers.felem[E->tracers.Q0compo_U_Th_tracnr][el] * U_Th__withdecay;
	}
       /* HPE radiogenic heating of the solid mantle (BMO_hpe_mantle_on=1):
        * Q_HPE = H_mantle*mult*cHPE(el)*exp(-t*ln2/tau), cHPE = element average of
        * tracfunc 83.  No mask needed: the BMO melt has a Dirichlet T boundary condition. */
       if (E->control.BMO_hpe_mantle_on
           && E->tracers.hpe_mantle_tracnr >= 0
           && E->tracers.felem[E->tracers.hpe_mantle_tracnr] != NULL
           && E->control.BMO_tau_hpe > 1.0e-10f)
         {
           float hpe_decay = expf(-(float)M_LN2 * (float)E->monitor.elapsed_time
                                  / E->control.BMO_tau_hpe);
           Q0comp_excess += E->control.BMO_H_mantle * E->control.BMO_H_mantle_mult
                          * E->tracers.felem[E->tracers.hpe_mantle_tracnr][el]
                          * hpe_decay;
         }
      }

    // do not consider inflow elements (12 columns):
    // weird velocities there!
    if (E->control.shear_heating>0
	&& E->parallel.me_loc[1]==0){
      if (((el-1)%E->lmesh.ely)>=E->lmesh.elz*12) 
    	visc_diss=E->visc_diss[el];
      else visc_diss=0.;}
    else visc_diss=E->visc_diss[el];
    /* COMMENTCOMMENTCOMMENTCOMMENTCOMMENTCOMMENT
    Q=0.0;
    for(i=0;i<Q0.number;i++)
	Q += Q0.Q[i] * exp(-Q0.lambda[i] 
             * (E->monitor.elapsed_time+Q0.t_offset));
    COMMENTCOMMENTCOMMENTCOMMENTCOMMENTCOMMENT */

    /* construct residual from this information */
    if(diffusion)
    { if(3==dims) 
      { for(j=1;j<=ends;j++) 
        { Eres[j]=0.0;
          for(i=1;i<=vpts;i++) 
	    {
	      Q = E->control.Q0+(v2[i]*adiab)-(double) (E->Qmelt[el]-visc_diss-Q0comp_excess);
           // Q = E->control.Q0+(v2[i]*adiab)-E->Qmelt[el]+visc_diss;
              Eres[j] -= PG.vpt[GNVINDEX(j,i)] 
                       * E->gDA[el].vpt[i] 
                       * (dT[i] - Q + v1[i]*tx1[i] + v2[i]*tx2[i] + 
                                    v3[i]*tx3[i]) 
                     + diff*E->gDA[el].vpt[i] 
                       * (E->gNX[el].vpt[GNVXINDEX(0,j,i)]*tx1[i] +
                          E->gNX[el].vpt[GNVXINDEX(1,j,i)]*tx2[i] +
                          E->gNX[el].vpt[GNVXINDEX(2,j,i)]*tx3[i] ); 
	   }}
      }
      else  /* 2-D */
      { 
	for(j=1;j<=ends;j++) 
        { Eres[j]=0.0;
          for(i=1;i<=vpts;i++) 
	    {
              Q = E->control.Q0+(v2[i]*adiab)-(double) (E->Qmelt[el]-visc_diss-Q0comp_excess);

            Eres[j] -= PG.vpt[GNVINDEX(j,i)] 
                       * E->gDA[el].vpt[i] 
                       * (dT[i] - Q + v1[i] * tx1[i] + v2[i] * tx2[i]) 
                     + diff*E->gDA[el].vpt[i] 
                       * (E->gNX[el].vpt[GNVXINDEX(0,j,i)] * tx1[i] +
                          E->gNX[el].vpt[GNVXINDEX(1,j,i)] * tx2[i] ); 
	   }}
      }
    }
    else /* no diffusion term */
    { if(3==dims) 
      { for(j=1;j<=ends;j++) 
        { Eres[j]=0.0;
          for(i=1;i<=vpts;i++) 
	    {//Q = E->control.Q0+v2[i]*adiab-E->Qmelt[el]+visc_diss;
            Q = E->control.Q0+(v2[i]*adiab)-(double) (E->Qmelt[el]-visc_diss-Q0comp_excess);
            Eres[j] -= PG.vpt[GNVINDEX(j,i)] 
                       * E->gDA[el].vpt[i] 
                       * (dT[i] - Q + v1[i] * tx1[i] + v2[i] * tx2[i] 
                          + v3[i] * tx3[i]);
	    }}
      }
      else  /* 2-D */
      { for(j=1;j<=ends;j++) 
        { Eres[j]=0.0;
          for(i=1;i<=vpts;i++) 
	    {//Q = E->control.Q0+v2[i]*adiab-E->Qmelt[el]+visc_diss;
            Q = E->control.Q0+(v2[i]*adiab)-(double) (E->Qmelt[el]-visc_diss-Q0comp_excess);
            Eres[j] -= PG.vpt[GNVINDEX(j,i)] 
                       * E->gDA[el].vpt[i] 
                       * (dT[i] - Q + v1[i] * tx1[i] + v2[i] * tx2[i]);
	    }}
      }
    }        

    /* See brooks etc: the diffusive term is excused upwinding for 
       rectangular elements  */
  
    /* include BC's for fluxes at (nominally horizontal) edges 
       (X-Y plane) */


    /* COMMENTCOMMENTCOMMENTCOMMENTCOMMENTCOMMENTCOMMENTCOMMENT
       This part takes care of NON-ZERO BOTTOM heatflux, and doesn't
       work correctly. For ZERO bottom heatflux this part may be 
       skipped, because BC[2] (see below) is zero, and so it won't
       add to the residue (JvH)
    if(FLAGS!=NULL) 
    {  onedfns=0;
       for(a=1;a<=ends;a++)
          if (FLAGS[E->ien[el].node[a]] & FBZ) 
          {  if (!onedfns++) get_global_1d_shape_fn(E,el,&GM,&dGamma);

             nodes[1] = loc[loc[a].node_nebrs[0][0]].node_nebrs[2][0];
             nodes[2] = loc[loc[a].node_nebrs[0][1]].node_nebrs[2][0];
             nodes[4] = loc[loc[a].node_nebrs[0][0]].node_nebrs[2][2];
             nodes[3] = loc[loc[a].node_nebrs[0][1]].node_nebrs[2][2];
         
             for(aid=0,j=1;j<=onedvpoints[E->mesh.nsd];j++)
                if (a==nodes[j])
                    aid = j;
             if(aid==0)  
                printf("%d: mixed up in pg-flux int: looking for %d\n",
                    el,a);

             if (loc[a].plus[1] != 0)
                back_front = 0;
             else back_front = dims;

             for(j=1;j<=onedvpoints[dims];j++)
                for(k=1;k<=onedvpoints[dims];k++)
                   Eres[a] += dGamma.vpt[GMVGAMMA(1+back_front,j)] 
                              * E->M.vpt[GMVINDEX(aid,j)] 
                              * g_1d[j].weight[dims-1] 
                              * BC[2][E->ien[el].node[a]] 
                              * E->M.vpt[GMVINDEX(k,j)];
          }
    } 
    COMMENTCOMMENTCOMMENTCOMMENTCOMMENTCOMMENTCOMMENTCOMMENT */
    
    return; 
}




/* =====================================================
   Obtain largest possible timestep (no melt considered)
   =====================================================  */


void std_timestep(E)
     struct All_variables *E;

{ 
    static int been_here = 0;
    static int counter = 0;
    static float diff_timestep,root3,root2;
    int lev,i,d,n,nel;

    float global_fmin();
    int input_float();

    float adv_timestep;
    float ts,uc2,uc3,uc1,uc,size,step;
    float     MO_cool;
    
    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int ends=enodes[dims];
    
    nel=E->lmesh.nel;
    lev = E->mesh.levmax;
    MO_cool=0.7111*E->control.MO_coolingtime/((float)(E->lmesh.NOZ[lev])); // factor 0.711 is for geometry of MSfsvs*ut

    /* For decaying-speed boundary modes (157/158), tighten MO_cool based on
       the instantaneous front velocity, which is fastest at early times.     */
    if (E->control.init_age_distr == 157)
    {   float frac = 1.0f - (float)(E->monitor.elapsed_time / E->control.MO_coolingtime);
        if (frac < 0.01f) frac = 0.01f;
        float v_now = 2.0f * (1.0f - E->control.MO_depth_initial) * frac
                    / (float)(E->control.MO_coolingtime);
        if (v_now > 1.0e-10f)
            MO_cool = 0.7111f / ((float)(E->lmesh.NOZ[lev]) * v_now);
    }
    else if (E->control.init_age_distr == 158)
    {   float tau_frac = E->control.MO_decay_tau;
        if (tau_frac <= 0.0f) tau_frac = 0.33f;
        float tau   = tau_frac * (float)(E->control.MO_coolingtime);
        float t_now = (float)(E->monitor.elapsed_time);
        float T_tot = (float)(E->control.MO_coolingtime);
        float norm  = 1.0f - expf(-T_tot / tau);
        if (norm < 1.0e-12f) norm = 1.0e-12f;
        float v_now = ((1.0f - E->control.MO_depth_initial) / tau)
                    * expf(-t_now / tau) / norm;
        if (v_now > 1.0e-10f)
            MO_cool = 0.7111f / ((float)(E->lmesh.NOZ[lev]) * v_now);
    }
    else if (E->control.init_age_distr == 159)
    {   /* Stefan sqrt(t): v(t) = (1-d0)/(2*sqrt(t*T)); regularised with t_min = T/NOZ^2. */
        float t_now = (float)(E->monitor.elapsed_time);
        float T_tot = (float)(E->control.MO_coolingtime);
        float noz_f = (float)(E->lmesh.NOZ[lev]);
        float t_min = T_tot / (noz_f * noz_f);
        float t_eff = (t_now > t_min) ? t_now : t_min;
        float v_now = (1.0f - E->control.MO_depth_initial)
                    / (2.0f * sqrtf(t_eff * T_tot));
        if (v_now > 1.0e-10f)
            MO_cool = 0.7111f / (noz_f * v_now);
    }
    else if (E->control.init_age_distr == 160 || E->control.BMO_implicit_mode==1)
    {   /* Implicit Labrosse coupling: cap dt so the BMO front moves at most one local
           grid row per step, using the latest zdot (fallback: MO_coolingtime estimate). */
        float v_now = fabsf(E->control.BMO_zdot);
        if (v_now < 1.0e-10f)
            v_now = (1.0f - E->control.MO_depth_initial)
                  / ((E->control.MO_coolingtime > 1.0e-12f)
                     ? E->control.MO_coolingtime : 1.0e-12f);
        if (v_now > 1.0e-10f) {
            float dz_local;
            int   row_front = E->control.MO_depth_row;
            if (row_front >= 1 && row_front < E->lmesh.noz)
                dz_local = E->X[2][row_front+1] - E->X[2][row_front];
            else
                dz_local = 1.0f / (float)(E->lmesh.NOZ[lev]);
            if (dz_local < 1.0e-12f)
                dz_local = 1.0f / (float)(E->lmesh.NOZ[lev]);
            MO_cool = 0.7111f * dz_local / v_now;
        }
    }

    if(E->advection.fixed_timestep != 0.0) {
      E->advection.timestep = E->advection.fixed_timestep;
      return;
    }

    if (been_here == 0)
    {  diff_timestep = 1.0e8;
       /* Guard: if inittreduct was read as 0 (e.g. unset), treat as 1 (no reduction) */
       if (E->control.inittreduct <= 0.0f)
           E->control.inittreduct = 1.0f;
       for(i=1;i<=nel;i++)
          for(d=1;d<=dims;d++)
          {  ts = E->eco[i].size[d] * E->eco[i].size[d];
             if (diff_timestep > ts) diff_timestep = ts;
          }
       diff_timestep = 0.5*global_fmin(E,diff_timestep);
       if (E->tracers.mindiffts > 0.) 
           diff_timestep = max(diff_timestep, E->tracers.mindiffts);
       been_here++;
       fprintf(stderr,"Diffusion timestep = %g\n", diff_timestep);
    }
    
    adv_timestep = 1.0e8;
    if (E->control.timestepmax >0. && E->control.timestepmax<1.0e8) adv_timestep = E->control.timestepmax;
 
    /* --- Begin additional CFL / safety constraints --- */
    {
      const float CFL = (E->control.CFL_number > 0.0f) ? E->control.CFL_number : 0.5f;
      float adv_timestep_local = adv_timestep; /* original adv estimate */

      /* Per-element Courant limiter: dt <= CFL * min(size) / (|v| + eps). */
      {
        float min_dt_by_courant = 1.0e32f;
        const float eps = 1e-12f;
        int iel, nd;
        for(iel=1; iel<=nel; iel++) {
            /* compute element-local characteristic length (minimum of sizes) */
            float lmin = 1e30f;
            for(nd=1; nd<=dims; nd++) {
                if (E->eco[iel].size[nd] > 0.0f && E->eco[iel].size[nd] < lmin)
                   lmin = E->eco[iel].size[nd];
            }
            if (lmin <= 0.0f) continue;
            /* compute element-average velocity magnitude (consistent with other code) */
            float ux=0.0f, uy=0.0f, uz=0.0f;
            if (3==dims) {
              int n;
              for(n=1;n<=ENODES3D;n++) {
                int node = E->ien[iel].node[n];
                ux += E->N.ppt[GNPINDEX(n,1)] * E->V[1][node];
                uy += E->N.ppt[GNPINDEX(n,1)] * E->V[2][node];
                uz += E->N.ppt[GNPINDEX(n,1)] * E->V[3][node];
              }
            } else {
              int n;
              for(n=1;n<=ENODES2D;n++) {
                int node = E->ien[iel].node[n];
                ux += E->N.ppt[GNPINDEX(n,1)] * E->V[1][node];
                uy += E->N.ppt[GNPINDEX(n,1)] * E->V[2][node];
              }
            }
            float vmag = sqrt(ux*ux + uy*uy + uz*uz);
            float dt_local = CFL * lmin / (vmag + eps);
            if (dt_local < min_dt_by_courant) min_dt_by_courant = dt_local;
        } /* end element loop */

        if (min_dt_by_courant < adv_timestep_local) {
            adv_timestep_local = min_dt_by_courant;
            if (E->parallel.me==0) 
               fprintf(stderr,"[CFL guard] reducing adv_timestep -> %g (CFL=%g)\n",
                       (double)adv_timestep_local,(double)CFL);
        }
      }


    /* finally set adv_timestep to the min between original and local-limited */
    adv_timestep = adv_timestep_local;
    }
    /* --- End additional CFL / safety constraints --- */


    adv_timestep = 1.0e-32 +
                   min((E->advection.fine_tune_dt*adv_timestep),diff_timestep);
    
    if ((E->monitor.elapsed_time<E->control.MO_coolingtime) && (E->mesh.topvbc >= 2) && adv_timestep>MO_cool)
      adv_timestep=MO_cool;
    /* Implicit Labrosse mode (moving bottom boundary): cap dt to MO_cool until the front reaches the CMB. */
    if ((E->control.init_age_distr==160 || E->control.BMO_implicit_mode==1)
        && (E->mesh.botvbc >= 2)
        && (E->control.MO_depth < 1.0f - 1.0e-6f)
        && adv_timestep > MO_cool)
      adv_timestep = MO_cool;
 
    E->advection.timestep = E->advection.dtmult * 
                            global_fmin(E,adv_timestep);


    /* Early-step dt reduction: factor ramps linearly from inittreduct to 1.0 over n_early_steps. */
    if (E->control.inittreduct < 1.0f && E->control.n_early_steps > 0)
    {
        int cyc    = E->monitor.solution_cycles;
        int nsteps = E->control.n_early_steps;
        if (cyc < nsteps)
        {
            float factor = E->control.inittreduct
                         + (1.0f - E->control.inittreduct) * (float)cyc / (float)nsteps;
            E->advection.timestep *= factor;
            if (E->parallel.me == 0)
                fprintf(stderr,"[early-ramp] step %d/%d: dt *= %.4g -> %g\n",
                        cyc, nsteps, (double)factor,
                        (double)E->advection.timestep);
        }
    }
     
    if (E->parallel.me == 0)
       fprintf(stderr,"Counter=%d: timestep=%g, time=%g\n",
          counter, E->advection.timestep,E->monitor.elapsed_time);
    
    counter++;

    return; 
}


/* Extrapolate the velocity field from old solutions to take timesteps without
   re-solving Stokes (may be unsafe outside steady-state runs). */


int get_v_estimate(E)
    struct All_variables *E;
{

    int i;
    float alpha,beta,gamma,delta,t3;
    float dt1,dt2,dt3,delt;

    static int here;
    static float previous_delta_v;

    dt1=E->monitor.elapsed_time_vsoln-E->monitor.elapsed_time_vsoln1;
    delt = E->monitor.elapsed_time-E->monitor.elapsed_time_vsoln;

    if(delt < 1.0e-30) {
	alpha = 0.0;
	beta = 0.0;
	gamma = 0.0;
	if(3==E->mesh.dof) {
	    for(i=1;i<=E->lmesh.nno;i++)  {
		alpha += E->V[1][i] * E->V[1][i] + E->V[2][i] * E->V[2][i] + E->V[3][i] * E->V[3][i]; 
		beta += (E->V1[1][i] - E->V[1][i]) * (E->V1[1][i] - E->V[1][i]) +
		    (E->V1[2][i] - E->V[2][i]) * (E->V1[2][i] - E->V[2][i]) +
		    (E->V1[3][i] - E->V[3][i]) * (E->V1[3][i] - E->V[3][i]);  
		gamma += (E->Vest[1][i] - E->V[1][i]) * (E->Vest[1][i] - E->V[1][i]) +
		    (E->Vest[2][i] - E->V[2][i]) * (E->Vest[2][i] - E->V[2][i])+
		    (E->Vest[3][i] - E->V[3][i]) * (E->Vest[3][i] - E->V[3][i]);
	    }
	}
	else {
	    for(i=1;i<=E->lmesh.nno;i++) {
		alpha += E->V[1][i] * E->V[1][i] + E->V[2][i] * E->V[2][i];
		beta += (E->V1[1][i] - E->V[1][i]) * (E->V1[1][i] - E->V[1][i]) +
		  (E->V1[2][i] - E->V[2][i]) * (E->V1[2][i] - E->V[2][i]);
		gamma += (E->Vest[1][i] - E->V[1][i]) * (E->Vest[1][i] - E->V[1][i]) +
		  (E->Vest[2][i] - E->V[2][i]) * (E->Vest[2][i] - E->V[2][i]);
	    }
	}
	if(alpha != 0.0)
	    previous_delta_v = sqrt(max(beta,gamma)/alpha);
	else
	    previous_delta_v = 1.0e16;	

	/* If delt is so small, then a better estimate of the timescale for which 
	   this solution will be used can be found from the current timestep */
	
	delt = E->advection.timestep;
	
	
    }


    // first 30 steps: no extrapolation
    if(E->monitor.solution_cycles < 2001){
        for(i=1;i<=E->lmesh.nno;i++) { 
	    E->Vest[1][i] = E->V[1][i];
	    E->Vest[2][i] = E->V[2][i];
	    if(E->mesh.dof==3)
		E->Vest[3][i] = E->V[3][i];
        }
    }
    else {
	/*t3 = dt1*dt2*(dt2-dt1);
	  alpha = 1.0 + ( (dt1*dt1-dt2*dt2) + (dt2-dt1)*delt  ) * delt / t3;
	  beta =  (dt2 - delt) * dt2 * delt / t3;
	  gamma = (delt - dt1) * dt1 * delt / t3;
	  delta = 0.0;*/
	/*alpha = 1.0 + delt*(-delt*delt + delt*dt1 + delt*dt2 - dt1*dt2 + delt*dt3 - dt1*dt3 - dt2*dt3) / (dt1*dt2*dt3);
	  beta =  delt * (-delt+dt2) *(-delt+dt3) / (dt1*(-dt1+dt2)*(-dt1+dt3));
	  gamma = delt * (-delt+dt1) *(-delt+dt3) / (dt2*( dt1-dt2)*(-dt2+dt3));
	  delta = delt * (-delt+dt1) *(-delt+dt2) / (dt3*( dt1-dt3)*( dt2-dt3)); */

	beta = -delt/dt1;
	alpha = 1.0 - beta;
	for(i=1;i<=E->lmesh.nno;i++) {
	    E->Vest[1][i] = E->V[1][i] * alpha + beta * E->V1[1][i];
	    E->Vest[2][i] = E->V[2][i] * alpha + beta * E->V1[2][i];
	    if(E->mesh.dof==3)
		E->Vest[3][i] = E->V[3][i] * alpha + beta * E->V1[3][i];
	}

  /* --- clamp Vest to avoid extreme external extrapolation --- */
  float max_allowed_scale = 2.0f;  // est magnitude no more than 3x current V 
  int inode;
  for(inode=1; inode<=E->lmesh.nno; inode++) {
      float vmag = sqrt(E->V[1][inode]*E->V[1][inode] + E->V[2][inode]*E->V[2][inode]
                      + (E->mesh.dof==3 ? E->V[3][inode]*E->V[3][inode] : 0.0f));
      float vestmag = sqrt(E->Vest[1][inode]*E->Vest[1][inode] + E->Vest[2][inode]*E->Vest[2][inode]
                      + (E->mesh.dof==3 ? E->Vest[3][inode]*E->Vest[3][inode] : 0.0f));
      if (vmag < 1e-12f) vmag = 1e-12f;
      if (vestmag > max_allowed_scale * vmag) {
          float scale = (max_allowed_scale * vmag) / (vestmag + 1e-12f);
          E->Vest[1][inode] *= scale;
          E->Vest[2][inode] *= scale;
          if (E->mesh.dof==3) E->Vest[3][inode] *= scale;
            }
          }  
      }

    here ++;

 
    if((E->monitor.solution_cycles < 3) ) {
	return(0);
    }
    
    if(here >= E->advection.max_substeps) {
	here=0;
	return(0);
    }
    
    if((previous_delta_v * delt / dt1)  > E->control.accuracy*E->control.sub_stepping_factor)  {
	here=0;
	return(0);
    }
    
    return(1);
    
}

/* node_limiter.c -- node-based limiter to enforce local monotonicity of T
   Insert call: apply_node_limiter(E, E->T); */
void apply_node_limiter(struct All_variables *E, float *T)
{
    int i, el, j, node;
    const int nno = E->lmesh.nno;
    const int nel = E->lmesh.nel;
    /* allocate arrays for min/max */
    float *Tmin = (float*) malloc((nno+1) * sizeof(float));
    float *Tmax = (float*) malloc((nno+1) * sizeof(float));
    int *seen = (int*) malloc((nno+1) * sizeof(int));
    if (!Tmin || !Tmax || !seen) {
        if (Tmin) free(Tmin);
        if (Tmax) free(Tmax);
        if (seen) free(seen);
        return;
    }

    for (i=1;i<=nno;i++) { Tmin[i]=1e30f; Tmax[i]=-1e30f; seen[i]=0; }

    /* loop over elements and update neighbor extremes for nodes of each element */
    for (el=1; el<=nel; el++) {
        int *ienmatrix = E->ien[el].node;
        int ends = enodes[E->mesh.nsd]; /* number of nodes per element */
        int k;
        for (j=1; j<=ends; j++) {
            int ni = ienmatrix[j];
            /* record the node's temperature itself as neighbor (so it doesn't collapse to NaN) */
            if (T[ni] < Tmin[ni]) Tmin[ni] = T[ni];
            if (T[ni] > Tmax[ni]) Tmax[ni] = T[ni];
            seen[ni] = 1;
            /* also compare with other nodes in the same element */
            for (k=1; k<=ends; k++) {
                int nj = ienmatrix[k];
                if (T[nj] < Tmin[ni]) Tmin[ni] = T[nj];
                if (T[nj] > Tmax[ni]) Tmax[ni] = T[nj];
            }
        }
    }

    /* Now apply limiter: clamp T[i] into [Tmin[i], Tmax[i]] */
    for (i=1; i<=nno; i++) {
        if (!seen[i]) continue; /* isolated node? skip */
        float t = T[i];
        if (t < Tmin[i]) T[i] = Tmin[i];
        else if (t > Tmax[i]) T[i] = Tmax[i];
    }

    /* Finally, ensure global bounds [0,1] */
    for (i = 1; i <= nno; i++) {
    if (E->T[i] < 0.0)
        E->T[i] = 0.0;
    else if (E->T[i] > 1.0)
        E->T[i] = 1.0;
    }

    free(Tmin); free(Tmax); free(seen);
    return;
}