/* Assumes parameter list is opened and reads the things it needs. 
   Variables are initialized etc, default values are set */


#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "prem_density.h"
#include <stdlib.h> /* for "system" command */
#include <strings.h>

void set_convection_defaults(E)
     struct All_variables *E;
{
    void PG_timestep_with_melting();
    void PG_timestep();
    void PG_process();
    void read_convection_settings();
    void convection_derived_values();
    void convection_allocate_memory();
    void convection_boundary_conditions();
    void node_locations();
    void convection_initial_fields();
    void twiddle_thumbs();
 
    E->next_buoyancy_field = PG_timestep;
    E->special_process_new_velocity = PG_process; 
    E->special_process_new_buoyancy = twiddle_thumbs; 
    E->problem_settings = read_convection_settings;
    E->problem_derived_values = convection_derived_values;
    E->problem_allocate_vars = convection_allocate_memory;
    E->problem_boundary_conds = convection_boundary_conditions;
    E->problem_initial_fields = convection_initial_fields;
    E->problem_node_positions = node_locations;
    E->problem_update_node_positions = twiddle_thumbs;
    E->problem_update_bcs = twiddle_thumbs;

    sprintf(E->control.which_data_files,"Temp,Strf,Pres");
    sprintf(E->control.which_horiz_averages,"Temp,Visc,Vrms");
    sprintf(E->control.which_running_data,"Step,Time,");
    sprintf(E->control.which_observable_data,"Shfl");
 
return;
}

void read_convection_settings(E)
     struct All_variables *E;
    
{ 
    void advection_diffusion_parameters();
    
/* parameters */

    input_float("rayleigh",&(E->control.Atemp),"essential");
    
    input_boolean("halfspace",&(E->convection.half_space_cooling),"off");
    input_float("halfspage",&(E->convection.half_space_age),"nodefault");
    
    input_int("temperature_blobs",&(E->convection.temp_blobs),"0");
    input_float_vector("temperature_blobx",E->convection.temp_blobs,E->convection.temp_blob_x);
    input_float_vector("temperature_bloby",E->convection.temp_blobs,E->convection.temp_blob_y);
    input_float_vector("temperature_blobz",E->convection.temp_blobs,E->convection.temp_blob_z);
    input_float_vector("temperature_blobsize",E->convection.temp_blobs,E->convection.temp_blob_radius);
    input_float_vector("temperature_blobDT",E->convection.temp_blobs,E->convection.temp_blob_T);
    input_float_vector("temperature_blobbg",E->convection.temp_blobs,E->convection.temp_blob_bg);
    input_int_vector("temperature_blobsticky",E->convection.temp_blobs,E->convection.temp_blob_sticky);
    
    input_int("temperature_zones",&(E->convection.temp_zones),"0");
    input_float_vector("temperature_zonex1",E->convection.temp_zones,E->convection.temp_zonex1);
    input_float_vector("temperature_zonex2",E->convection.temp_zones,E->convection.temp_zonex2);
    input_float_vector("temperature_zonez1",E->convection.temp_zones,E->convection.temp_zonez1);
    input_float_vector("temperature_zonez2",E->convection.temp_zones,E->convection.temp_zonez2);
    input_float_vector("temperature_zoney1",E->convection.temp_zones,E->convection.temp_zoney1);
    input_float_vector("temperature_zoney2",E->convection.temp_zones,E->convection.temp_zoney2);
    input_float_vector("temperature_zoney2",E->convection.temp_zones,E->convection.temp_zoney2);
    input_float_vector("temperature_zoney2",E->convection.temp_zones,E->convection.temp_zoney2);
    input_float_vector("temperature_zonehw",E->convection.temp_zones,E->convection.temp_zonehw);
    input_float_vector("temperature_zonemag",E->convection.temp_zones,E->convection.temp_zonemag);
    input_int_vector("temperature_zonesticky",E->convection.temp_zones,E->convection.temp_zone_sticky);
    
    input_int("num_perturbations",&(E->convection.number_of_perturbations),"0,0,32");
    input_float_vector("perturbmag",E->convection.number_of_perturbations,E->convection.perturb_mag);
    input_float_vector("perturbk",E->convection.number_of_perturbations,E->convection.perturb_k);
    
    input_string("prevT",E->convection.old_T_file,"initialize");

  if (E->control.restart)  {
       input_int("restart_timesteps",&(E->control.first_tstep),"0");
       E->monitor.solution_cycles = E->control.first_tstep;
       input_string("oldfile",E->convection.old_T_file,"initialize");
       }
  else 
       E->control.first_tstep = 0;


    
    advection_diffusion_parameters(E);
    
    return;
}

/* =================================================================
   Any setup which relates only to the convection stuff goes in here
   ================================================================= */

void convection_derived_values(E)  
     struct All_variables *E;
 
{ 

return;
}

void convection_allocate_memory(E)
     struct All_variables *E;

{ void advection_diffusion_allocate_memory();

  advection_diffusion_allocate_memory(E);

return;
}

/* ============================================ */
    
void convection_initial_fields(E)
     struct All_variables *E;

{ 
    void convection_initial_temperature();

    report(E,"convection, initial temperature");
    convection_initial_temperature(E);

  return; }

/* =========================================== */

void convection_boundary_conditions(E)
     struct All_variables *E;

{
    void velocity_boundary_conditions();
    void temperature_boundary_conditions();
    void temperatures_conform_bcs();
 
    velocity_boundary_conditions(E);      /* universal */
    temperature_boundary_conditions(E);
    temperatures_conform_bcs(E);

    return;
}

/* ===============================
   Initialization of fields .....
   =============================== */

void convection_initial_temperature(E)
     struct All_variables *E;
{
   int i,j,k,l,p,node,ii;
   double temp,base,radius,radius2;
   const float factor = E->data.layer_km*E->data.density*E->data.grav_acc*1e-6;
   double adiab=E->data.adiabatic_gr*E->data.layer_km/E->data.ref_temperature;
   float TsolP,pressure;
   double drand48();
   int slabs();
   int crust();
   FILE *fp;
   void remove_horiz_ave();
   void vcopy();
   void temperatures_conform_bcs();
   void generic_data_storage();
   void thermal_buoyancy();
   void average_node_f20();
   void tracexit();
  
   int nfz,in1,in2,in3,instance,nox,noy,noz;
   char input_s[201],output_file[255],filename[255];
   float T_file[1003],depth[1003];
   float val2,x1,y1,z1,weight,para1,plate_velocity,delta_temp,age,onsetssc;
   float Xp, Yp, Rp2, r2, dT, pi;
  
   float randomnr;
   int rand();
 
   int read_previous_field();

   const int dims=E->mesh.nsd;

   pi = 3.141592635;

   noy=E->lmesh.noy;  
   noz=E->lmesh.noz;  
   nox=E->lmesh.nox;  
 
   if (E->parallel.me==0)
      fprintf(stderr,"[CIND-entry] convection_initial_temperature: init_age_distr=%d  restart=%d  MO_depth_initial=%g  MO_temperature=%g  BMO_TL_dPdim=%g\n",
              E->control.init_age_distr, E->control.restart,
              E->control.MO_depth_initial, E->control.MO_temperature,
              E->control.BMO_TL_dPdim);

   if (E->control.restart==0)
   { for(i=1;i<=noy;i++)  
        for(j=1;j<=nox;j++) 
           for(k=1;k<=noz;k++)
           {  node=k+(j-1)*noz+(i-1)*noz*nox;
              x1=E->X[1][node];
              z1=E->X[2][node];
              y1= (E->mesh.nsd==3) ? E->X[3][node] : 0.0; 
 
              if (E->control.init_age_distr<10) 
              {  temp = 0.5*E->X[2][k]/sqrt(E->control.plate_age);
                  E->segment.Tz_in[k] = erf(temp)+E->X[2][k]*adiab;
              }
              else if (E->control.init_age_distr==10) 
              {   E->segment.Tz_in[k] = z1;
                  E->segment.Tz_out[k] = z1;
              }
              else if (E->control.init_age_distr==11) 
                  if (k!=1)
                     E->segment.Tz_in[k] = 1.;
                  else 
                     E->segment.Tz_in[k] = 0.;
              else if (E->control.init_age_distr==12) 
              {  temp = 0.5*E->X[2][k]/sqrt(E->control.plate_age);
                  E->segment.Tz_in[k] = erf(temp)+z1*adiab;
              }
                          
              if (E->control.init_age_distr==0) 
              {  /* constant age at t=0 */
                 E->T[node] =  E->segment.Tz_in[k];
              }
              else if (E->control.init_age_distr==1) 
              { /* lin. age distr. at t=0 */
                age = E->control.plate_age + x1/E->control.plate_vel;
                temp = 0.5*z1/sqrt(age);
                E->T[node] = erf(temp)+z1*adiab;
              }
              else if (E->control.init_age_distr==3) 
              { /* as age_distr=1, but for x>xtrac0 then T=1*/
                /* Check if tracers are used and the right * type: */
                age = E->control.plate_age + x1/E->control.plate_vel;
                temp = 0.5*z1/sqrt(age);
                if (x1<5.0+z1)
                   E->T[node] = erf(temp)+z1*adiab;
                else
                   E->T[node] = 1.0;      
              }
	      else if (E->control.init_age_distr==4)
		{ /*  for a fracture zone at i=noy/2, middlemost node  */
		  /*  elements i=>noy/2 are concerned plate_age=plate_angle  */
		if (i < noy/2) 
                  age = E->control.plate_age + x1/E->control.plate_vel;
		else
		  age = E->control.plate_angle+x1/E->control.plate_vel;
                temp = 0.5*z1/sqrt(age);
                E->T[node] = erf(temp)+z1*adiab; 
		}
	      else if (E->control.init_age_distr==5)
		{ /*  for a fracture zone at i=noy/2 and i=noy/2+4  */
		  /*  elements noy/2=<i<=noy/2+4 are concerned plate_age=plate_angle  */
		  if (i <(noy/2) || i>(noy/2+4)) 
                  age = E->control.plate_age + x1/E->control.plate_vel;
		else
		  age = E->control.plate_angle+x1/E->control.plate_vel;
                temp = 0.5*z1/sqrt(age);
                E->T[node] = erf(temp)+z1*adiab; 
		}
	      else if (E->control.init_age_distr==6)
		{
		  E->T[node]=(x1+z1)/2;
		}
	      else if (E->control.init_age_distr==7)
		{
		  E->T[node]=1.;
		}
              else if (E->control.init_age_distr==8)
		{ /*    ssc-initial-profile, valid for ridge melt   *
		   *    depleted buoyant lid and EA=120 kJ/mol      */
		  age = E->control.plate_age + x1/E->control.plate_vel;
                  //intf(" %g ",age);
	          onsetssc = E->control.onset_ssc;
                  //intf(" %g ",onsetssc);
                  if (age > onsetssc)
                    age=onsetssc;
		  temp = 0.5*z1/sqrt(age);
                  E->T[node] = erf(temp)+z1*adiab;  
		}
              else if (E->control.init_age_distr==10) 
              { /* linear T-z relation: T=z, plus sinusoidal perturb */
                E->T[node] = z1;
                for(l=1;l<=E->convection.number_of_perturbations;l++)
                  E->T[node] += E->convection.perturb_mag[l-1]*
                    (float)
                    ( cos((double)(E->convection.perturb_k[l-1]*pi*x1))*
                      sin((double)(pi*z1)));
              }
              else if (E->control.init_age_distr==11) 
                  if (k!=1)
                     E->T[node] = 1.;
                  else 
                     E->T[node] = 0.;
              else if (E->control.init_age_distr==12) 
              { /* lin. age distr. at t=0 */
                age = E->control.plate_age;
                temp = 0.5*z1/sqrt(age);
                E->T[node] = erf(temp)+z1*adiab;
                
		pressure=z1*factor;
                TsolP = 1100.+102.04*pressure;
		/* if (z1<0.5) TsolP+=100.0*(0.5-z1); * real solidus due to initial linear depletion in Raddick-model *
		TsolP/=E->data.ref_temperature;
		if (E->T[node]>(TsolP))
		E->T[node]=TsolP;*/
                
		for(l=1;l<=E->convection.number_of_perturbations;l++)
                  E->T[node] += E->convection.perturb_mag[l-1]*
                    (float)
                    ( cos((double)(E->convection.perturb_k[l-1]*pi*x1))*
                      sin((double)(pi*z1)));
              }
              else if (E->control.init_age_distr==20) 
              {
                E->T[node] = 0.;
              }
              else if (E->control.init_age_distr>=100 && E->control.init_age_distr<130)
              {
                E->T[node] = 1.;
              }
              else if (E->control.init_age_distr>=130 && E->control.init_age_distr<150)
	      { if (z1>=E->control.MO_depth_initial)
		  E->T[node] = ((1.-z1)*E->control.MO_temperature+z1-E->control.MO_depth_initial)/(1.-E->control.MO_depth_initial);
		else E->T[node] = E->control.MO_temperature;
		if (j==1)
		  fprintf(stderr,"CIND 136: MO_T-final=%g, z=%g, T=%g\n",
                   E->control.MO_temp_final,z1, E->T[node]);
	      }
	      else if (E->control.init_age_distr == 160)
	      {
	        /* Implicit BMO coupling (Labrosse-2025): no external file needed.
	         * Solid mantle (z < z0): linear profile with scalar BMO_TL_dPdim
	         *   (shallow-mantle liquidus slope; differs from layers[0]).
	         * BMO melt (z >= z0): uniform T = MO_temperature.
	         * Runtime T_BOMB / Q_eq at the front still use layered mP(z0)
	         * via BMO_get_mP_at_z in update_BMO_implicit.                 */
	        float z0_init = E->control.BMO_z0_init;
	        if (z0_init < 1.0e-6f)
	            z0_init = E->control.MO_depth_initial;
	        if (z1 >= z0_init)
	            E->T[node] = E->control.MO_temperature;
	        else
	            E->T[node] = E->control.MO_temperature
	                       + E->control.BMO_TL_dPdim
	                         * (z1 - z0_init);
	        if (i==1 && j==1 && k==1)
	            fprintf(stderr,
	                "CIND 160: liquidus init, MO_T=%g, z0=%g,"
	                "  mP_shallow=%g, mP(z0)=%g, nlayers=%d\n",
	                E->control.MO_temperature, (double)z0_init,
	                (double)E->control.BMO_TL_dPdim,
	                (double)BMO_get_mP_at_z(E, z0_init),
	                E->control.BMO_TL_dPdim_nlayers);
	      }
	      else if (E->control.init_age_distr>160 && E->control.init_age_distr<170)
	      {
	        if (node==1)
	           { sprintf(filename,"%s.Temp_solidus.input",E->control.data_file);
	             fp = fopen(filename,"r");
	             if (fp==NULL)
                       fprintf(E->fp,"CIND 16X: file reading error\n");
	             for (ii=1; ii<=1001; ii++)
	               { fgets(input_s,200,fp);
	                 sscanf(input_s,"%g",&val2);
	                 T_file[ii] = val2/E->data.ref_temperature;
	                 depth[ii]=((float)(ii-1))/1000.;
	                 //fprintf(E->fp, "testoo %d %g %g\n",ii,depth[ii],T_file[ii]);
	               }
		     T_file[0]=T_file[1]/2.;  // mid mantle T T_file[0] could be higher than 0.5 if internal heating
	             if (E->control.MO_temperature != T_file[1]) tracexit(E, "Critical ERROR in CIND 16X: Exit (check MO_temperature)");;
		     fclose(fp);
	           }
		 if (z1>E->control.MO_depth_initial)
		  {
	           ii = 2 + ((int)((z1-E->control.MO_depth_initial)*1000));
	           E->T[node] = (T_file[ii-1]*(depth[ii]-z1) + T_file[ii]*(z1-depth[ii-1]))/ (depth[ii]-depth[ii-1]);
                  }
		 else if (z1<(E->control.MO_depth_initial/2.)) {
                   age = E->control.plate_age;
                   temp = 0.5*z1/sqrt(age);
                   E->T[node] = erf(temp)*T_file[0];
                  }
		 else {
                   age = E->control.plate_age;
                   temp = 0.5*(E->control.MO_depth_initial-z1)/sqrt(age);
                   E->T[node] = T_file[1]-erf(temp)*T_file[0];
		  }
		}
              else if (E->control.init_age_distr>=150 && E->control.init_age_distr<160)
                {
                if (node==1)
                   { sprintf(filename,"%s.Temp_solidus.input",E->control.data_file);
                     fp = fopen(filename,"r");
                     if (fp==NULL)
                       fprintf(E->fp,"CIND 15X: file reading error\n");
                     for (ii=1001; ii>=1; ii--)
                       { fgets(input_s,200,fp);
                         sscanf(input_s,"%g",&val2);
                         T_file[ii] = val2/E->data.ref_temperature;
                         depth[ii]=((float)(ii-1))/1000.;
                         //fprintf(E->fp, "testoo %d %g %g\n",ii,depth[ii],T_file[ii]);
                       }
                     fclose(fp);
                   }
                 ii = 2 + ((int)(z1*1000));
                 E->T[node] = (T_file[ii-1]*(depth[ii]-z1) + T_file[ii]*(z1-depth[ii-1]))/ (depth[ii]-depth[ii-1]);
                 if (j==1)
                  fprintf(stderr,"CIND 15X: MO_T-final=%g, z=%g, T=%g\n",
                   E->control.MO_temp_final,z1,E->T[node]);
                }
              if (E->control.randomT_init != 0.0)
              {  /* Add noise with ampl randomT to T-inflow field */
                 randomnr = 2*(float)(rand()) / (float)(RAND_MAX) - 1.;
                 /* randomnr = drand48(); */
                 E->T[node] += E->control.randomT_init * randomnr; 
                 /*
                 fprintf(stderr,"Convection: ampl=%g, rand=%g, T=%g\n",
                           E->control.randomT_bc,randomnr, E->T[node]);
                 */
              }
 
 
              if (E->control.init_plume_cyl>0 && z1>E->viscosity.zlm)
              {  /* A hot plume cylinder in the lowermost layer is 
                    imposed as initial cond. Lateral location & T of 
                    cylinder equals bnd.cnd. for plume at bottom.
                 */ 
                 Rp2 = E->segment.plume_radius*E->segment.plume_radius;
                 Xp  = E->segment.plume_coord[1];
                 Yp  = E->segment.plume_coord[3];
                 r2  = (x1-Xp)*(x1-Xp);
                 if (E->mesh.nsd==3)
                    r2+=(y1-Yp)*(y1-Yp);
                 dT  = E->segment.plume_DT*exp(-r2/Rp2);
                 E->T[node] += dT;
              }
 
              E->node[node] = E->node[node] | (INTX | INTZ | INTY);
 	  }
       average_node_f20(E, E->T, E->mesh.levmax);
   }
   else  
   { sprintf(output_file,"%s.temp.%d.%d",E->convection.old_T_file,
                          E->parallel.me,E->monitor.solution_cycles);
     fp = fopen(output_file,"r");
     fgets(input_s,200,fp);
     sscanf(input_s,"%d %d %g",&i,&j,&E->monitor.elapsed_time);
     for(node=1;node<=E->lmesh.nno;node++)
     {  fgets(input_s,200,fp);
        sscanf(input_s,"%g",&E->T[node]);
        E->node[node] = E->node[node] | (INTX | INTZ | INTY);
     }
     fclose(fp);
   }

   temperatures_conform_bcs(E);

   thermal_buoyancy(E);

   return; 
}


void setup_plume_problem(E)
     struct All_variables *E;
{
    int iy, iz;
    float age;
    FILE *fp;
    
    int l,nox,noy,noz;
    char output_file[255];
    double temp,temp1,outflow_age,time_scale,velo_scale;
    char input_s[200];
    int i, j, node;
    double adiab=E->data.adiabatic_gr*E->data.layer_km/E->data.ref_temperature;
    int differentn,nequalone;
    float depth,n1;

    float point_visc();
    int layers();
    void xmin_flux_bc();
    void xmin_flux_bc2();
    void xmin_flux_bc_rst();

    nox=E->lmesh.nox;  
    noy=E->lmesh.noy;  
    noz=E->lmesh.noz;  

    E->control.MO_depth_initial/=E->data.layer_km;
    E->control.MO_depth=E->control.MO_depth_initial;

    if (E->control.BMO_implicit_mode) {
        if (E->control.BMO_z0_init < 1.0e-6f)
            E->control.BMO_z0_init = E->control.MO_depth_initial;
        if (E->control.BMO_TL0 < 1.0e-6f)
            E->control.BMO_TL0 = E->control.MO_temperature;
    }

    E->control.MO_depth_row=0;
    for (i=1; i<=noz; i++)
      { 
        depth=E->X[2][i];
//   fprintf (stderr,"depth %g %d\n", depth,i);
        if (depth>=E->control.MO_depth && E->control.MO_depth_row==0)
         E->control.MO_depth_row=i;
      }
     fprintf (stderr,"MOdepth %g %d\n", E->control.MO_depth_initial, E->control.MO_depth_row);

    init_BMO_rho_nv_prof(E);

    E->control.plate_vel = E->control.plate_vel*0.01/3.155734e7;  
    E->control.old_plate_vel = E->control.old_plate_vel*0.01/3.155734e7;
           /* now  with a unit of m/sec */
    velo_scale = E->data.therm_diff/(E->data.layer_km*1.0e3);
            /* in m/sec */ 
    E->control.plate_vel /= velo_scale;  
    E->control.old_plate_vel /= velo_scale;  
            /* dimensionless */

    time_scale = E->data.layer_km*E->data.layer_km*1.0e6/
                 E->data.therm_diff;
            /* in sec */ 
    E->control.onset_ssc=pow(E->control.Atemp,-0.68)*36.35*51551*1.0e6*3.155734e7;
   
    E->control.plate_age = E->control.plate_age*1.0e6*3.155734e7;
            /* in sec */ 
    E->control.onset_ssc = E->control.onset_ssc/time_scale;
    fprintf(stderr ,"onsetssc: %g \n", E->control.onset_ssc); 
    E->control.plate_age = E->control.plate_age/time_scale;
            /* dimensionless */
    fprintf(stderr ,"plateage: %g \n", E->control.plate_age);

    if (E->control.init_age_distr==4 || E->control.init_age_distr==5)
      E->control.plate_angle = (E->control.plate_angle*1.0e6*3.155734e7)/time_scale;
    else
    E->control.plate_angle *= 3.141592635/180.;
            /* degs -> rads */ 


    /* Essential thermal boundary conditions */
    if (E->control.plate_angle == 0.)
    {  if (E->control.init_age_distr==0) 
       {  for(iz=1;iz<=noz;iz++)
          {  /* Inflow+outflow conditions */
             temp = 0.5*E->X[2][iz]/sqrt(E->control.plate_age);
             E->segment.Tz_in[iz] = erf(temp)+E->X[2][iz]*adiab;
             E->segment.Tz_out[iz] = erf(temp)+E->X[2][iz]*adiab;
          }
       }
       else if (E->control.init_age_distr==1) 
       {  for(iz=1;iz<=noz;iz++)
          {  /* Inflow conditions */
             temp = 0.5*E->X[2][iz]/sqrt(E->control.plate_age);
             E->segment.Tz_in[iz] = erf(temp)+E->X[2][iz]*adiab;
             /* Outflow conditions */
             outflow_age = E->control.plate_age +  
                        E->mesh.layer[1]/E->control.plate_vel;
             temp = 0.5*E->X[2][iz]/sqrt(outflow_age);
             E->segment.Tz_out[iz] = erf(temp)+E->X[2][iz]*adiab;
          }
       }
       else if (E->control.init_age_distr==2) 
       {  /* Read in T-field of former solution */
          sprintf(output_file,"%s.temp.%d.%d",
                          E->convection.old_T_file,
                          E->parallel.me,
                          E->monitor.solution_cycles);
          fp = fopen(output_file,"r");
          fgets(input_s,200,fp);
          sscanf(input_s,"%d %d %g",&i,&j,&E->monitor.elapsed_time);
          for(iz=1;iz<=noz;iz++)
          {  fgets(input_s,200,fp);
             sscanf(input_s,"%g",&E->T[node]);
             E->segment.Tz_in[iz] = E->T[node];
          }
          fclose(fp);
       }
       else if (E->control.init_age_distr==3) 
       {  for(iz=1;iz<=noz;iz++)
          {  /* Inflow conditions */
             temp = 0.5*E->X[2][iz]/sqrt(E->control.plate_age);
             E->segment.Tz_in[iz] = erf(temp)+E->X[2][iz]*adiab;
             /* Outflow conditions */
             E->segment.Tz_out[iz] = 1.0;
          }
       }

       else if (E->control.init_age_distr==4)
       	 printf("ERROR: plate_angle should be > 0.0\n");

       else if (E->control.init_age_distr==8)
	 { /*    ssc-initial-profile, valid for ridge melt   *
	    *    depleted buoyant lid and EA=120 kJ/mol      */
           /* Inflow conditions */
 	   for(iz=1;iz<=noz;iz++) {
             temp = 0.5*E->X[2][iz]/sqrt(E->control.plate_age);
             E->segment.Tz_in[iz] = erf(temp)+E->X[2][iz]*adiab;
	     /* Outflow conditions */
	     outflow_age = min(E->control.onset_ssc,
	                     E->control.plate_age+E->mesh.layer[1]/E->control.plate_vel);
	     temp = 0.5*E->X[2][iz]/sqrt(outflow_age);
	     E->segment.Tz_out[iz] = erf(temp)+E->X[2][iz]*adiab;
	   }
	 }
       else if (E->control.init_age_distr==10) 
       {
          for(iz=1;iz<=noz;iz++)
          {   E->segment.Tz_in[iz] = E->X[2][iz];
              E->segment.Tz_out[iz] = E->X[2][iz];
          }
       }
       else if (E->control.init_age_distr==11 || E->control.init_age_distr>=100) 
       {
          E->segment.Tz_in[1] = 0.;
          E->segment.Tz_out[1] = 0.;
          for(iz=2;iz<=noz;iz++)
          {  E->segment.Tz_in[iz] = 1.;
             E->segment.Tz_out[iz] = 1.;
          }
       }
       else if (E->control.init_age_distr==12)
          {  E->segment.Tz_in[1] = 0.;
             E->segment.Tz_out[1] = 0.;
             for(iz=2;iz<=noz;iz++)
             {  temp = 0.5*E->X[2][iz]/sqrt(E->control.plate_age);
                E->segment.Tz_in[iz] = erf(temp)+E->X[2][iz]*adiab;
             }
          }

    }
    else 
    {  if (E->control.init_age_distr>=10)
       {  fprintf(stderr,"Error in setup_plume_problem: STOP\n");
          fprintf(stderr,"   plate_angle = %g is not 0, while\n", 
                              E->control.plate_angle);
          fprintf(stderr,"   init_age_distr = %d >=10\n", 
                              E->control.init_age_distr);
          exit(10);
       }
       if (E->mesh.nsd<3)
       {  fprintf(stderr,"Error in setup_plume_problem: STOP\n");
          fprintf(stderr,"   Problem=%d-D is not 3-D, while\n", 
                              E->mesh.nsd);
          fprintf(stderr,"   init_age_distr = %d >=10\n", 
                              E->control.init_age_distr);
          exit(10);
       }
       /* Age varies in y-direction along inflow-bnd. (x=0-bnd):
        * age = E->control.plate_angle at (x,y)=(0,1), and 
        * increases with decreasing y
        */
       if (E->control.init_age_distr==1)
       {  /* Functional description of T */
          for(iy=1;iy<=noy;iy++)
          {  age = E->control.plate_age 
                 + sin(E->control.plate_angle)
                   * (E->mesh.layer[3]-E->X[3][iy*nox*noz]) 
                   / E->control.old_plate_vel;
             for(iz=1;iz<=noz;iz++)
             {  temp = 0.5*E->X[2][iz]/sqrt(age);
                E->segment.Tzy_in[iz+noz*(iy-1)] = erf(temp)+E->X[2][iz]*adiab;
             }
             fprintf(E->fp, "age = %g, y = %g, T_iz=6 =%g\n", 
                     age, E->X[3][iy*nox*noz],
                     E->segment.Tzy_in[6+noz*(iy-1)]);
          }
       }
       else if (E->control.init_age_distr==2)
       {  /* Read in T-field of former solution */
          sprintf(output_file,"%s.temp.%d.%d",
                          E->convection.old_T_file,
                          E->parallel.me,E->monitor.solution_cycles);
          fp = fopen(output_file,"r");
          fgets(input_s,200,fp);
          sscanf(input_s,"%d %d %g",&i,&j,&E->monitor.elapsed_time);
          for(node=1;node<=E->lmesh.nno;node++)
          {  fgets(input_s,200,fp);
             sscanf(input_s,"%g",&E->T[node]);
          }
          fclose(fp);
          for(iy=1;iy<=noy;iy++)
             for(iz=1;iz<=noz;iz++)
             {  node=iz+(iy-1)*noz*nox;
                E->segment.Tzy_in[iz+noz*(iy-1)] = E->T[node];
             }
       }
       else if (E->control.init_age_distr==4)
       {  /* Functional description of T */
	 for(iy=1;iy<(noy/2);iy++)
	   {  for(iz=1;iz<=noz;iz++)
              {  temp = 0.5*E->X[2][iz]/sqrt(E->control.plate_age);
                 E->segment.Tzy_in[iz+noz*(iy-1)] =+E->X[2][iz]*adiab;
              }
            fprintf(E->fp, "age = %g, y = %g, T_iz=6 =%g\n", 
                     age, E->X[3][iy*nox*noz],
                     E->segment.Tzy_in[6+noz*(iy-1)]);
            }
	   for(iy=noy/2;iy<=noy;iy++)
	    {  
             for(iz=1;iz<=noz;iz++)
	       {  temp = 0.5*E->X[2][iz]/sqrt(E->control.plate_angle);
		 E->segment.Tzy_in[iz+noz*(iy-1)] = erf(temp)+E->X[2][iz]*adiab;
             }
             fprintf(E->fp, "age = %g, y = %g, T_iz=6 =%g\n", 
                     age, E->X[3][iy*nox*noz],
                     E->segment.Tzy_in[6+noz*(iy-1)]);
            }
       }
     else if (E->control.init_age_distr==5)
       {  /* Functional description of T */
	 for(iy=1;iy<noy/2;iy++)
	   {  for(iz=1;iz<=noz;iz++)
              {  temp = 0.5*E->X[2][iz]/sqrt(E->control.plate_age);
                 E->segment.Tzy_in[iz+noz*(iy-1)] = erf(temp)+E->X[2][iz]*adiab;
              }
            fprintf(E->fp, "age = %g, y = %g, T_iz=6 =%g\n", 
                     age, E->X[3][iy*nox*noz],
                     E->segment.Tzy_in[6+noz*(iy-1)]);
            }
	 for(iy=noy/2;iy<=(noy/2+4);iy++)
	    {  
             for(iz=1;iz<=noz;iz++)
	       {  temp = 0.5*E->X[2][iz]/sqrt(E->control.plate_angle);
                E->segment.Tzy_in[iz+noz*(iy-1)]=erf(temp)+E->X[2][iz]*adiab;
             }
             fprintf(E->fp, "age = %g, y = %g, T_iz=6 =%g\n", 
                     age, E->X[3][iy*nox*noz],
                     E->segment.Tzy_in[6+noz*(iy-1)]);
            }
	 for(iy=noy/2+5;iy<noy;iy++)
	   {  for(iz=1;iz<=noz;iz++)
              {  temp = 0.5*E->X[2][iz]/sqrt(E->control.plate_age);
                 E->segment.Tzy_in[iz+noz*(iy-1)]=erf(temp)+E->X[2][iz]*adiab;;
              }
            fprintf(E->fp, "age = %g, y = %g, T_iz=6 =%g\n", 
                     age, E->X[3][iy*nox*noz],
                     E->segment.Tzy_in[6+noz*(iy-1)]);
	   }
	
       }
    }
    
    /* Material flux boundary condition at xmin: */
    if (E->control.restart)
       xmin_flux_bc_rst(E);
    else 
    {  n1 = E->viscosity.sdepv_expt[0];
       differentn=0;
       nequalone=1;
       for (l=2;l<=E->viscosity.num_mat;l++)
         {if (E->viscosity.sdepv_expt[l-1] != 1)
             nequalone=0;
         if (E->viscosity.sdepv_expt[l-1] != n1)
            differentn=1;}
       if (differentn || nequalone) 
          xmin_flux_bc(E,0); 
       else 
          xmin_flux_bc2(E,0);
	//fprintf(E->fp,"nequalone %i\n",nequalone);
    }
   
    /* non-dimensionalize other parameters */
    E->segment.plume_radius /= E->data.layer_km;
    E->segment.plume_DT /= E->data.ref_temperature;
    E->segment.plume_coord[1] /= E->data.layer_km;
    if (E->mesh.nsd==3)
        E->segment.plume_coord[3] /= E->data.layer_km;

    if(E->control.fluxcontrol==1)
    {  if(E->parallel.me == 0)
          fprintf(stderr,"Dimensional desired plume flux:%g\n",
                   E->control.flux);
       E->control.flux /= (E->data.density*E->data.therm_exp*
                           E->data.ref_temperature*E->data.therm_diff*
                           E->data.layer_km*1000.);
       if(E->parallel.me == 0)
          fprintf(stderr,"Non-dimensional desired plume flux:%g\n",
                   E->control.flux);
    }

  if (E->parallel.me==0)
      fprintf(stderr,"[SPP-exit] setup_plume_problem done. init_age_distr=%d\n",
              E->control.init_age_distr);
  return;
}


void xmin_flux_bc(E, etaknown)
   struct All_variables *E;
   int etaknown;
{
    int iy, iz;

    int l,nox,noy,noz;
    double temp, temp1;
    float compo0,compo1=0.;
    float point_visc();
    float *segment_Tz;
    int layers();

    nox=E->lmesh.nox;
    noy=E->lmesh.noy;
    noz=E->lmesh.noz;

    
    if (E->control.plate_angle == 0.)
    { /* Purely z-dependent inflow viscosity profile determines inflow
       * velocity profile:
       *   eta*dvdz = constant = sigma_zz
       *   Integrating over whole z-range gives:
       *     INT_0^1 sigma_zz/eta dz = INT_0^1 dvdz dz = -vplate
       *   so:
       *     sigma_zz = -vplate / (INT_0^1 1/eta dz)
       *   now dvdz = sigma_zz / eta, and integrating to z gives:
       *     v = vplate + INT_0^z sigma_zz/eta dz
       * Error traps are set above, assuming this is all correct 
       */
       temp = 0.0;
       if (E->control.xminvxtype == 8) segment_Tz=&E->segment.Tz_out[0];
          else segment_Tz=&E->segment.Tz_in[0];
       for(iz=1;iz<=noz;iz++)
       {  /* l = E->mat[iz]; */
          l = layers(E,iz); 
          if (E->tracers.compoviscoused ==1 && E->viscosity.RHEOL == 111)
             {  
		compo0 = E->tracers.init_at_inflow0[iz];
                compo1 = E->tracers.init_at_inflow1[iz];
             }
             else 
                compo0 = 0.;
          if (etaknown==0)
             E->Have.Vi[iz] =
                point_visc(E, E->viscosity.RHEOL, E->viscosity.N0[l-1],
                           E->viscosity.E[l-1], E->viscosity.Z[l-1],
                           E->viscosity.T[l-1], E->viscosity.dN[l-1],
                           E->viscosity.Ttransition, 
                           *(segment_Tz+iz),(double)(E->X[1][iz]), 
                           (double)(.5*(E->X[2][iz]+E->X[2][iz-1])),
			   1.,1.,compo0,compo1);
          if (iz>1) 
             temp += (double)(1/E->Have.Vi[iz])
                     *(E->X[2][iz]-E->X[2][iz-1]);
	fprintf(E->fp,"eta:%g dz:%g temp:%g\n", (double)(1/E->Have.Vi[iz]),
	(E->X[2][iz]-E->X[2][iz-1]),temp);
       }
       temp  = -E->control.plate_vel/temp;
       temp1 = 0.0;
       for(iz=1;iz<=noz;iz++)
       {  if (iz>1)  
            temp1 += (double)(1/E->Have.Vi[iz])
                   * (E->X[2][iz]-E->X[2][iz-1]);
          E->segment.Vx[iz] = (double)(temp*temp1) + E->control.plate_vel;
          fprintf(E->fp,"n:%d eta:%g T:%g v:%g l:%d z:%g",
		  iz,E->Have.Vi[iz], *(segment_Tz+iz),
                 E->segment.Vx[iz], E->mat[iz], E->X[2][iz]);
          if (iz>1)  
             fprintf(E->fp,"edot:%g stress: %g\n",
                (E->segment.Vx[iz]-E->segment.Vx[iz-1])/
                (E->X[2][iz]- E->X[2][iz-1]),
                (E->segment.Vx[iz]-E->segment.Vx[iz-1])/
                (E->X[2][iz]- E->X[2][iz-1])*E->Have.Vi[iz]);

          else 
             fprintf(E->fp,"\n");
       }
    }
    else 
    { /* y- and z-dependent inflow viscosity profile determines inflow
       * velocity profile:
       *   eta*dvdz = constant = sigma_zz
       *   For each y:
       *   | Integrating over whole z-range gives:
       *   |   INT_0^1 sigma_zz/eta dz = INT_0^1 dvdz dz = -vplate
       *   | so:
       *   |   sigma_zz = -vplate / (INT_0^1 1/eta dz)
       *   | now dvdz = sigma_zz / eta, and integrating to z gives:
       *   |   v = vplate + INT_0^z sigma_zz/eta dz
       * Note that each y has a different sigma_zz!
       * Error traps are set above, so I assume this is correct 
       */
       for(iy=1;iy<=noy;iy++)
       {  temp = 0.0;
          for(iz=1;iz<=noz;iz++)
          {  l = layers(E,iz); 
             if (etaknown==0 && iz>1)
             /* 29-6-2004 JvH: last param to point_visc = C0 = 0 added.
                This means that this routine does not yet take 
                composition-dependent viscosity into account */
             /* 9-7-2004 JvH: last param to point_visc = C1 = 0 added.
                This means that this routine does not yet take 
                composition-dependent viscosity into account */
                E->Have.Vi[iz] =
                   point_visc(E, E->viscosity.RHEOL, 
                              E->viscosity.N0[l-1],
                              E->viscosity.E[l-1], E->viscosity.Z[l-1],
                              E->viscosity.T[l-1], E->viscosity.dN[l-1],
                              E->viscosity.Ttransition,
                              E->segment.Tzy_in[iz+noz*(iy-1)],(double)(E->X[1][iz]),
                              (double)(.5*(E->X[2][iz]+E->X[2][iz-1])),
			       1.,1.,compo0,compo1);
             else
             {  fprintf(stderr, 
                   "ERROR xmin_flux_bc: pl.angle/=0 & SDEPV unimpl.\n");
                exit(10);
             }
             if (iz>1) 
                temp += (1/E->Have.Vi[iz])
                        *(E->X[2][iz]-E->X[2][iz-1]);
          }
          temp  = -E->control.plate_vel/temp;
          temp1 = 0.0;
          for(iz=1;iz<=noz;iz++)
          {  if (iz>1)  
               temp1 += (1/E->Have.Vi[iz])
                      * (E->X[2][iz]-E->X[2][iz-1]);
             E->segment.Vx_zy[iz+noz*(iy-1)] = temp*temp1   
                                             + E->control.plate_vel;
             fprintf(E->fp,
                "iz=%d iy=%d visc=%g T=%g v=%g layer=%d z=%g\n",
                iz, iy, E->Have.Vi[iz], 
                E->segment.Tzy_in[iz+noz*(iy-1)],
                E->segment.Vx_zy[iz+noz*(iy-1)], E->mat[iz], 
                E->X[2][iz]);
          }
       }
    }

    return;
}

void xmin_flux_bc2(E)
   struct All_variables *E;
{
    int iy, iz;

    int l,nox,noy,noz;
    double temp, temp1;

    float point_arrh();
    float *segment_Tz;
    int layers();

    nox=E->lmesh.nox;
    noy=E->lmesh.noy;
    noz=E->lmesh.noz;

    
    if (E->control.plate_angle == 0.)
    { /* Purely z-dependent inflow viscosity profile determines inflow
       * velocity profile:
       *   eta*dvdz = constant = sigma_zz
       *   Integrating over whole z-range gives:
       *     INT_0^1 sigma_zz/eta dz = INT_0^1 dvdz dz = -vplate
       *   so:
       *     sigma_zz = -vplate / (INT_0^1 1/eta dz)
       *   now dvdz = sigma_zz / eta, and integrating to z gives:
       *     v = vplate + INT_0^z sigma_zz/eta dz
       * Error traps are set above, assuming this is all correct 
       */
       temp = 0.0;
       if (E->control.xminvxtype == 8) segment_Tz=&E->segment.Tz_out[0];
       else segment_Tz=&E->segment.Tz_in[0];
       for(iz=1;iz<=noz;iz++)
       {  /* l = E->mat[iz]; */
          l = layers(E,iz); 
          E->Have.arrh[iz] =
             point_arrh(E, E->viscosity.RHEOL, E->viscosity.N0[l-1],
                        E->viscosity.E[l-1], E->viscosity.T[l-1], 
                        *(segment_Tz+iz));
          if (iz>1)
             temp += (E->Have.arrh[iz]+E->Have.arrh[iz-1])
                     *(E->X[2][iz]-E->X[2][iz-1])*0.5;
       }
       temp  = -E->control.plate_vel/temp;
       temp1 = 0.0;
       for(iz=1;iz<=noz;iz++)
       {  if (iz>1)  
            temp1 += (E->Have.arrh[iz]+E->Have.arrh[iz-1])
                     *(E->X[2][iz]-E->X[2][iz-1]) * 0.5;
          E->segment.Vx[iz] = temp*temp1 + E->control.plate_vel;
          fprintf(E->fp,"n:%d arrh:%g T:%g v:%g l:%d z:%g",
		  iz,E->Have.arrh[iz], *(segment_Tz+iz),
                 E->segment.Vx[iz], E->mat[iz], E->X[2][iz]);
          if (iz>1)  
             fprintf(E->fp," edot:%g\n",
                (E->segment.Vx[iz]-E->segment.Vx[iz-1])/
                (E->X[2][iz]- E->X[2][iz-1]));
          else 
             fprintf(E->fp,"\n");
       }
    }
    else 
    { /* y- and z-dependent inflow viscosity profile determines inflow
       * velocity profile:
       *   eta*dvdz = constant = sigma_zz
       *   For each y:
       *   | Integrating over whole z-range gives:
       *   |   INT_0^1 sigma_zz/eta dz = INT_0^1 dvdz dz = -vplate
       *   | so:
       *   |   sigma_zz = -vplate / (INT_0^1 1/eta dz)
       *   | now dvdz = sigma_zz / eta, and integrating to z gives:
       *   |   v = vplate + INT_0^z sigma_zz/eta dz
       * Note that each y has a different sigma_zz!
       * Error traps are set above, so I assume this is correct 
       */
       for(iy=1;iy<=noy;iy++)
       {  temp = 0.0;
          for(iz=1;iz<=noz;iz++)
          {  l = layers(E,iz); 
             E->Have.arrh[iz] =
                point_arrh(E, E->viscosity.RHEOL, E->viscosity.N0[l-1],
                           E->viscosity.E[l-1], E->viscosity.T[l-1], 
                           E->segment.Tzy_in[iz+noz*(iy-1)]);
             if (iz>1)
                temp += (E->Have.arrh[iz]+E->Have.arrh[iz-1])
                        *(E->X[2][iz]-E->X[2][iz-1])*0.5;
          }
          temp  = -E->control.plate_vel/temp;
          temp1 = 0.0;
          for(iz=1;iz<=noz;iz++)
          {  if (iz>1)  
               temp1 += (E->Have.arrh[iz]+E->Have.arrh[iz-1])
                        *(E->X[2][iz]-E->X[2][iz-1]) * 0.5;
             E->segment.Vx_zy[iz+noz*(iy-1)] = temp*temp1   
                                             + E->control.plate_vel;
             fprintf(E->fp,
                "iz=%d iy=%d visc=%g T=%g v=%g layer=%d z=%g\n",
                iz, iy, E->Have.Vi[iz], 
                E->segment.Tzy_in[iz+noz*(iy-1)],
                E->segment.Vx_zy[iz+noz*(iy-1)], E->mat[iz], 
                E->X[2][iz]);
          }
       }
    }
    return;
}

void xmin_flux_bc_rst(E)
   struct All_variables *E;
{
   /* In case of restart, the inflow vel. is simply read from the 
    * restart file
    */

   FILE *fp;
   char input_s[200], output_file[255];
   float dum, Vx, i, j;
   int iz, noz;

   noz=E->lmesh.noz;

   sprintf(output_file,"%s.temp.0.%d",E->convection.old_T_file,
                        E->control.first_tstep);
   fprintf(stderr, "Reading restart inflow profile from %s\n",
           output_file);
   fp = fopen(output_file,"r");
   fgets(input_s,200,fp);
   sscanf(input_s,"%g %g %g",&i,&j,&dum);
   for(iz=1;iz<=noz;iz++)
   {  fgets(input_s,200,fp);
      sscanf(input_s,"%g %g",&dum, &Vx);
      E->segment.Vx[iz] = Vx;
      fprintf(stderr, "Vx[%d]=%g\n", iz, E->segment.Vx[iz]);
   }
   fclose(fp);
    
   return;
}

void PG_process(E,ii)
     struct All_variables *E;
     int ii;
{ 
    void remove_horiz_ave();
    void strain_rate_2_inv();
    void generic_data_storage();
    float return_bulk_value();
    
    float *P,*P2;
    float visc[9];

    int i,j,k,p,a1,nint,n,el;
    int node,hnode;        
    int this_node,that_node;

    struct Shape_function GN;
    struct Shape_function_dx GNx;
    struct Shape_function_dA dOmega;
    
    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];

    float x1,x2,int1,int2;

    P=(float *)malloc((1+E->lmesh.nno)*sizeof(double));
    P2=(float *)malloc((1+E->lmesh.nel)*vpts*sizeof(double));
 
    /* This is an ideal point to calculate the energy consistency integral that
       Slava suggested to me.. NOTE, the buoyancy and velocity are most closely in step
       at this point, before advection. */

    strain_rate_2_inv(E,P2,0);  /* strain rate invariant squared */
    v_to_nodes(E,P2,P,E->mesh.levmax);
    
    int1=0.0;
    int2=0.0;

    for(el=1;el<=E->lmesh.nel;el++) {
	
	for(i=1;i<=ends;i++) {
	    n = E->ien[el].node[i];
	    for(j=1;j<=vpts;j++) {
		int1 += E->EVi[(el-1)*vpts+j] * P[n] * E->N.vpt[GNVINDEX(i,j)] * E->gDA[el].vpt[j];
		int2 += E->V[2][n]*E->buoyancy[n] * E->N.vpt[GNVINDEX(i,j)] * E->gDA[el].vpt[j];
	    }
	}
    }

    
    if(int1 != 0.0)
	fprintf(E->fp,"Energy balance integral accurate to  %g%% in %e\n",100.0*fabs((int1-int2)/int1),int1);
  
    
    free((void *)P);
    free((void *)P2);
    return;
}



