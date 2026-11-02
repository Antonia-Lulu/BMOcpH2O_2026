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
/* Functions relating to the determination of viscosity field either
   as a function of the run, as an initial condition or as specified from
   a previous file */

#if 1
#include <mpi.h>
#endif
#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

void viscosity_parameters(E)
     struct All_variables *E;
{
    int i,l;
    float temp;

/*
    int input_int();
    int input_float();
*/
  
  /* default values .... */
    E->viscosity.update_allowed = 0; 
    E->viscosity.SDEPV = E->viscosity.TDEPV = 0;
    E->viscosity.EXPX=0;
  

   for(i=0;i<40;i++) {
       E->viscosity.N0[i]=1.0;
       E->viscosity.T[i] = 0.0;
       E->viscosity.Z[i] = 0.0;
       E->viscosity.E[i] = 0.0;
       E->viscosity.T0[i] = 0.0;
   }

  /* read in information */
    input_int("rheol",&(E->viscosity.RHEOL),"essential");
    input_int("num_mat",&(E->viscosity.num_mat),"1");
 
    input_float_vector("viscT",E->viscosity.num_mat,(E->viscosity.T));  /* redundant */
    input_float_vector("viscT1",E->viscosity.num_mat,(E->viscosity.T));
    input_float_vector("viscZ",E->viscosity.num_mat,(E->viscosity.Z));
    input_float_vector("viscE",E->viscosity.num_mat,(E->viscosity.E));
    input_float_vector("viscT0",E->viscosity.num_mat,(E->viscosity.T0));
    input_float_vector("visc0",E->viscosity.num_mat,(E->viscosity.N0)); /* redundant */
    input_float_vector("viscN0",E->viscosity.num_mat,(E->viscosity.N0));
    input_float_vector("dvisc",E->viscosity.num_mat,(E->viscosity.dN));
    E->viscosity.subdweak_w = 0.1;
    input_float("subdweak_w",&(E->viscosity.subdweak_w),"0.");
    E->viscosity.subdweak_eta = 1.0;
    fprintf(stderr," subdweak_w=%g\n", E->viscosity.subdweak_w);
    input_float("subdweak_eta",&(E->viscosity.subdweak_eta),"1.0");
    fprintf(stderr," subdweak_eta=%g\n", E->viscosity.subdweak_eta);
  
    input_boolean("TDEPV",&(E->viscosity.TDEPV),"on");
    input_boolean("SDEPV",&(E->viscosity.SDEPV),"off");
    E->viscosity.byerlee = 0;
    input_int("byerlee",&(E->viscosity.byerlee),"0");
    fprintf(stderr,"Byerlee's yielding used (0=no, 1=yes):%d\n", 
             E->viscosity.byerlee);
    if (E->viscosity.byerlee == 1) 
    {  E->viscosity.byerleemu=0.6;
       input_float("byerleemu",&(E->viscosity.byerleemu),"0.6");
       fprintf(stderr,"byerleemu=%g\n", E->viscosity.byerleemu);
    }
    
    /* General stress-dependent rheology input */
    E->viscosity.sdepv_misfit=0.001;
    input_float("sdepv_misfit",&(E->viscosity.sdepv_misfit),"0.001");
    E->viscosity.sdepv_relax=0.;
    input_float("sdepv_relax",&(E->viscosity.sdepv_relax),"0.");
    input_float_vector("sdepv_expt",E->viscosity.num_mat,
                       (E->viscosity.sdepv_expt));
    input_float_vector("sdepv_trns",E->viscosity.num_mat,
                       (E->viscosity.sdepv_trns));

    /* Input for background strainrate in 3rd dir in 2-D models 
       used for stress-dependent rheology cases
       (for example in 2-D LR SSC with nonlinear rheology) */
    /* background strainrate available (no if =0, yes if >0) */
    E->viscosity.edot31in2Davail=0;
    input_int("edot31in2Davail",&(E->viscosity.edot31in2Davail),"0");
    fprintf(stderr, "at read-in:edot31in2Davail=%d\n", 
            E->viscosity.edot31in2Davail);
    /* constant backgr. strainrate value (option edot31in2Davail=1) */
    E->viscosity.edot31in2Dval=0.;
    input_float("edot31in2Dval",&(E->viscosity.edot31in2Dval),"0.0");
    /* backgr. strainrate calculated from given vplate in 3rd dir 
       (option edot31in2Davail=2) */
    E->viscosity.vplate3in2Dval=0.;
    input_float("vplate3in2Dval",&(E->viscosity.vplate3in2Dval),"0.0");
    E->viscosity.calcstress=0;

    input_boolean("TDEPV_AVE",&(E->viscosity.TDEPV_AVE),"off");
    input_boolean("VFREEZE",&(E->viscosity.FREEZE),"off");
    input_boolean("VMAX",&(E->viscosity.MAX),"off");
    input_boolean("VMIN",&(E->viscosity.MIN),"off");
    input_boolean("VISC_UPDATE",&(E->viscosity.update_allowed),"on");
  
    input_float("freeze_thresh",&(E->viscosity.freeze_thresh),"0.0");
    input_float("freeze_value",&(E->viscosity.freeze_value),"1.0");
    input_float("visc_max",&(E->viscosity.max_value),"nodefault");
    input_float("visc_min",&(E->viscosity.min_value),"nodefault");

    input_boolean("VISC_GUESS",&(E->viscosity.guess),"off");
    input_string("visc_old_file",E->viscosity.old_file," ");

  if (E->viscosity.RHEOL==1 || E->viscosity.RHEOL==111
     || (E->viscosity.RHEOL>=200 && E->viscosity.RHEOL<=219))
  { for(l=1;l<=E->viscosity.num_mat;l++) 
    {  E->viscosity.E[l-1] = E->viscosity.E[l-1]/(E->data.gas_const*E->data.ref_temperature);
       E->viscosity.T[l-1] = E->viscosity.T[l-1]/E->data.ref_temperature;
       temp = exp(E->viscosity.E[l-1]/(1.0+E->viscosity.T[l-1]));
       E->viscosity.N0[l-1] = E->viscosity.N0[l-1]/temp;
    }
    E->viscosity.byerleemu *= E->control.Atemp/E->data.therm_exp/
                              E->data.ref_temperature;
    fprintf(stderr,"Scaled byerleemu=%g\n", E->viscosity.byerleemu);
  }
  else if (E->viscosity.RHEOL==2)    {
    for(l=1;l<=E->viscosity.num_mat;l++) {
       E->viscosity.E[l-1] = E->viscosity.E[l-1]/(E->data.gas_const*E->data.ref_temperature);
       E->viscosity.Z[l-1] = E->viscosity.Z[l-1]*E->data.density*E->data.grav_acc*E->data.layer_km*1000/(E->data.gas_const*E->data.ref_temperature);
       E->viscosity.T[l-1] = E->viscosity.T[l-1]/E->data.ref_temperature;
       temp = exp((E->viscosity.E[l-1]+200/E->data.layer_km*E->viscosity.Z[l-1])/(1.0+E->viscosity.T[l-1]));
       E->viscosity.N0[l-1] = E->viscosity.N0[l-1]/temp;
       }
    }
  else if (E->viscosity.RHEOL==3 || E->viscosity.RHEOL==33 ||  E->viscosity.RHEOL==34)    {
    for(l=1;l<=E->viscosity.num_mat;l++) {
       E->viscosity.E[l-1] = E->viscosity.E[l-1]/
                            (E->data.gas_const*E->data.ref_temperature);
       E->viscosity.Z[l-1] = E->viscosity.Z[l-1]*E->data.density*
                             E->data.grav_acc*E->data.layer_km*1000/
                            (E->data.gas_const*E->data.ref_temperature);
       E->viscosity.T[l-1] = E->viscosity.T[l-1]/
                             E->data.ref_temperature;
       temp = exp(E->viscosity.E[l-1]+E->viscosity.Z[l-1]/
                  (1.0+E->viscosity.T[l-1])
                 );
       E->viscosity.N0[l-1] = E->viscosity.N0[l-1]/temp;
       }
    }
  else if (E->viscosity.RHEOL==100)    
     { /* Binary viscosity distribution:
             visc = visc_background = N0 for T <  ReferenceT+25 degC
             visc = visc_plume=N0*dN     for T >= ReferenceT+25 degC */ 
       E->viscosity.Ttransition = 1. + 25./E->data.ref_temperature;
     }
  else if (E->viscosity.RHEOL==101)
     { /* Binary viscosity distribution:
             N0_background = N0    for T <  ReferenceT+25 degC
             N0_plume      = N0*dN for T >= ReferenceT+25 degC */
       E->viscosity.Ttransition = 1. + 25./E->data.ref_temperature;
       for(l=1;l<=E->viscosity.num_mat;l++) 
       { E->viscosity.E[l-1] = E->viscosity.E[l-1]/
                     (E->data.gas_const*E->data.ref_temperature);
         E->viscosity.T[l-1] = E->viscosity.T[l-1]/
                     E->data.ref_temperature;
         temp = exp(E->viscosity.E[l-1]/(1.0+E->viscosity.T[l-1]));
         E->viscosity.N0[l-1] = E->viscosity.N0[l-1]/temp;
       }
    }
        
return;
}

void get_viscosity_option(E)
     struct All_variables *E;
{
    void viscosity_for_system();
 
    /* general, essential default */
  
    input_string("Viscosity",E->viscosity.STRUCTURE,NULL);   /* Which form of viscosity */
    
    input_boolean("VISC_EQUIVDD",&(E->viscosity.EQUIVDD),"off");    /* Whether to average it */
    input_int("equivdd_opt",&(E->viscosity.equivddopt),"1");
    input_int("equivdd_x",&(E->viscosity.proflocx),"1");
    input_int("equivdd_y",&(E->viscosity.proflocy),"1");
  
    input_int("update_every_steps",&(E->control.KERNEL),"1");

    input_boolean("VISC_SMOOTH",&(E->viscosity.SMOOTH),"off");
    input_int ("visc_smooth_cycles",&(E->viscosity.smooth_cycles),"0");
    
    if ( strcmp(E->viscosity.STRUCTURE,"system") == 0) /* Interpret */ {
      if(E->parallel.me==0) fprintf(E->fp,"Viscosity derived from system state\n");
      E->viscosity.FROM_SYSTEM = 1;
      viscosity_for_system(E);
    }

  return;

}

/* ============================================ */

void viscosity_for_system(E)
     struct All_variables *E;
{
    void get_system_viscosity();
    void twiddle_thumbs();

    get_system_viscosity(E,1,E->EVI[E->mesh.levmax],
                                E->VI[E->mesh.levmax]);
    return;
}


void apply_viscosity_smoother(E,visc,evisc)
     struct All_variables *E;
     float *visc,*evisc;

{
    void p_to_centres();
    void p_to_nodes();
    
    double  *ViscCentre;
    int i;

    ViscCentre = (double *)malloc((E->lmesh.nno+10)*sizeof(double));
   
    for(i=1;i<=E->viscosity.smooth_cycles;i++)  
    {   p_to_centres(E,visc,ViscCentre,E->mesh.levmax);
	p_to_nodes(E,ViscCentre,visc,E->mesh.levmax);
    }

    free ((void *)ViscCentre);
 
    return;
}

void get_system_viscosity(E,propogate,EEta,Eta)
    struct All_variables *E;
    float *Eta,*EEta;
    int propogate;
{
    int i,j,k,l,z,jj,kk, iy, iz, nox, noy, noz;
    float zero,one,temp,tempa,TT[9], compo0, compo1=0.;
    float *nedot, nexpo, EEta_old, tot, tot1;
    int lnode1, lnode2;
    float seglen;
    double xx, zz,temp2;
    float n1;

    static int visits=0;
    static int differentn=0;
    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];
    const int nel = E->lmesh.nel;
    const int dims=E->mesh.nsd;
    int lv, row;

    float point_visc();
    void strain_rate_2_inv();
    void apply_viscosity_smoother();
    void printarray();
    void visc_from_gint_to_nodes();
    void visc_from_nodes_to_gint();
    void edot_to_nodes();
    void xmin_flux_bc();
    void velocity_boundary_conditions();
    void horizontal_bc_within();

    //eedot = (float *) malloc((2+nel)*sizeof(float));
    nedot = (float *) malloc((2+E->lmesh.nno)*sizeof(float));

    one = 1.0;
    zero = 0.0;

    if (visits==0)
    {  E->viscosity.calcstress=1;
       if (E->viscosity.SDEPV)
       {  n1 = E->viscosity.sdepv_expt[0];
          differentn=0;
          for (l=2;l<=E->viscosity.num_mat;l++)
            if (E->viscosity.sdepv_expt[l-1] != n1)
               differentn=1;
       }
    }
 
    /* Print some info at first call */
    if(propogate && visits==0) 
    {  fprintf(E->fp,"Rheological option: %d\n", E->viscosity.RHEOL);
       for(l=1;l<=E->viscosity.num_mat;l++) 
          fprintf(E->fp,"\tlayer %d/%d: E=%g T1=%g N0=%g Z0=%g\n",
                 l,E->viscosity.num_mat,
                 E->viscosity.E[l-1],E->viscosity.T[l-1],
                 E->viscosity.N0[l-1],E->viscosity.Z[l-1]);
       fflush(E->fp);
    }

    /* In case of nonl.rh., calc. strain rate (per element).
     * First iteration, no old V available, take strainrate=1.
     */
    if (visits==0 || E->viscosity.SDEPV==0 ||
        E->control.eedot_available==0 )
       {if (E->control.shear_heating==0)
	  for(i=1;i<=nel;i++)
            E->eedot[i] = one;
	}
    else if (E->control.shear_heating==0)
    {  strain_rate_2_inv(E,E->eedot,1);
       if (E->viscosity.edot31in2Davail == 3)
          E->viscosity.calcstress=0;
    }
    else
      for(i=1;i<=nel;i++)
        E->eedot[i] = (float)(sqrt((double)(E->eedot[i]))); 

    

    if (visits==0 || E->monitor.solution_cycles%E->control.KERNEL==0) 
    {  for(i=1;i<=nel;i++)   
       {  l = E->mat[i];
          tempa = E->viscosity.N0[l-1];
          j = 0;
          temp2=0.0;
   
          for(kk=1;kk<=ends;kk++)
             TT[kk] = E->T[E->ien[i].node[kk]];
   
          if(E->viscosity.SDEPV)
             nexpo = E->viscosity.sdepv_expt[E->mat[i]-1];
          else 
             nexpo = 1.;
          
      

          for(jj=1;jj<=vpts;jj++) 
          {  temp=1.0e-32;
             xx = 0;
             zz = 0;
             for(kk=1;kk<=ends;kk++) 
             {  temp += max(zero,TT[kk])*E->N.vpt[GNVINDEX(kk,jj)];
                xx   += E->X[1][E->ien[i].node[kk]]*
                        E->N.vpt[GNVINDEX(kk,jj)];
                zz   += E->X[2][E->ien[i].node[kk]]*
                        E->N.vpt[GNVINDEX(kk,jj)];
             }
	     //if (i==1) fprintf(stderr,"%i %i-",E->tracers.compoviscoused,E->viscosity.RHEOL);
             if (E->tracers.compoviscoused == 1)
             {
                if (E->viscosity.RHEOL == 33 && l!=E->viscosity.num_mat) compo0 = 0.;
                else compo0 = E->tracers.felem[0][i];

                if (E->viscosity.RHEOL == 34 || 
                    E->viscosity.RHEOL == 103 ||
                    E->viscosity.RHEOL == 105 || 
                    E->viscosity.RHEOL == 111)
                   compo1 = E->tracers.felem[1][i];
                else 
                   compo1 = 0.;
             }
             else 
                compo0 = 0.;
	     if (E->viscosity.RHEOL==222)
                {
		  if (E->tracers.felem[0][i]<-50 && l==E->viscosity.num_mat) compo0=10.; else compo0=1.;
		}
             if (visits>0 && E->monitor.sdepv_cycles>0) /* Used for relaxation, see below */
                EEta_old = EEta[(i-1)*vpts + jj];
             EEta[(i-1)*vpts + jj] = 
                 point_visc( E, 
                      E->viscosity.RHEOL, E->viscosity.N0[l-1], 
                      E->viscosity.E[l-1], E->viscosity.Z[l-1],
                      E->viscosity.T[l-1], E->viscosity.dN[l-1],
                      E->viscosity.Ttransition, temp, xx, zz, E->eedot[i],
                      nexpo, compo0, compo1);

             
          //  fprintf(stderr, "elem=%d, gaus=%d, visc=%g\n", i, jj, 
            //      EEta[(i-1)*vpts + jj]);
              
                      

             /* Relaxation of viscosity solution: x*old+(1-x)*new 
              * - x is read from input file: sdepv_relax 
              * - Relaxation only performed for nonlin. rheology 
              * - and only AFTER the 1st timestep, otherwise
              *   averaging with the former timestep would be done
              */
             if (visits>0 && E->viscosity.SDEPV &&
                 E->monitor.sdepv_cycles>0)
                EEta[(i-1)*vpts + jj] = 
                   exp( (1.-E->viscosity.sdepv_relax) * 
                                log(EEta[(i-1)*vpts+jj])
                        + E->viscosity.sdepv_relax * 
                                log(EEta_old)
                      );
          }
       }
   
       /* Smoothing of the viscosity to stabilize the solution? 
        * -----------------------------------------------------
        */
       /* Not used
       visc_from_gint_to_nodes(E,EEta,Eta,E->mesh.levmax); 
       visc_from_nodes_to_gint(E,Eta,EEta, E->mesh.levmax);
       */

       /* Material flux boundary condition at xmin: */
       /* update xmin-flux-bc for 1st nonl.rh. timestep:
        * - interpolate eta to nodes (linear/expon.)
        * - calc. eta at inflow bnd. (2-D)
        *         or average eta at inflow plane (3-D)
        * - calc. inflow velocity field
        * - update velocity b.c.'s
        */
       if (E->monitor.solution_cycles==0 && 
           E->viscosity.SDEPV && 
           E->control.eedot_available && 
           differentn)
       {  /* interpolation eta_gauss -> eta_nodes
           * ------------------------------------ */
          visc_from_gint_to_nodes(E,EEta,Eta,E->mesh.levmax);
          nox=E->lmesh.nox;
          noy=E->lmesh.noy;
          noz=E->lmesh.noz;
          for(iz=1;iz<=noz;iz++)
          {  if (dims==2)
                E->Have.Vi[iz] = Eta[iz];
             else if (dims==3)
             {  tot1=0.;
                if (E->parallel.me_loc[1]==0) 
                {  for(iy=1;iy<noy;iy++)
                   {  lnode1 = iz + (iy-1)*noz*nox;
                      lnode2 = iz + iy*noz*nox;
                      seglen = E->X[3][lnode2]-E->X[3][lnode1];
                      tot1 += seglen*Eta[iz+noz*nox*(iy-1)];
                      tot1 += seglen*Eta[iz+noz*nox*(iy-1)];
                   }
                   tot1/=2;
                }
                MPI_Allreduce(&tot1,&tot,1,MPI_FLOAT,MPI_SUM,
                              MPI_COMM_WORLD);
                E->Have.Vi[iz] = tot/E->mesh.layer[3];
             }
             fprintf(E->fp, "inflow eta[%d]=%g\n", iz, E->Have.Vi[iz]);
          }
          xmin_flux_bc(E,1);
          velocity_boundary_conditions(E);
       }
       else if (E->mesh.topvbc >= 2 || E->mesh.botvbc >= 2)  {
         velocity_boundary_conditions(E); /*
        for(lv=E->mesh.levmax;lv>=E->mesh.levmin;lv--)
         {  row = (E->mesh.NOZ[lv] - 1) / 2 + 1;
            horizontal_bc_within(E,E->VB,E->control.MO_depth_row,1,1,0.0,VBX,1,lv);
            horizontal_bc_within(E,E->VB,E->control.MO_depth_row,1,1,0.0,SBX,0,lv);
            horizontal_bc_within(E,E->VB,E->control.MO_depth_row,2,1,0.0,VBZ,1,lv);
            horizontal_bc_within(E,E->VB,E->control.MO_depth_row,2,1,0.0,SBZ,0,lv);
            if(E->mesh.nsd==3)
            {   E->control.VBYtopval = 0.0;
                horizontal_bc_within(E,E->VB,E->control.MO_depth_row,3,1,0.0,VBY,1,lv);
                horizontal_bc_within(E,E->VB,E->control.MO_depth_row,3,1,0.0,SBY,0,lv);
            }
         } */
       }        

       printarray(E, "T", E->T, E->lmesh.nno);
       printarray(E, "V1", E->V[1], E->lmesh.nno);
       printarray(E, "V2", E->V[2], E->lmesh.nno);
       printarray(E, "EE", E->eedot, nel);
       printarray(E, "visc", EEta, nel*vpts);

       visits++;
    }
    

    if (E->viscosity.SMOOTH)
       apply_viscosity_smoother(E,Eta,EEta);

    // free ((void *)eedot);
    free ((void *)nedot);
}

float point_visc(E, RHEOL, N0, Eact, V, T0, dN, 
                 Ttr, T, x, z, eedot, nexpo, compo0, compo1)
    struct All_variables *E;
    int RHEOL;
    float N0, Eact, V, T0, dN, Ttr, T, eedot, nexpo, compo0, compo1;
    double x, z;
{
    int tracnr;
    float visc, exponent1, exponent2;
    float viscold;
    float w,dzz;
    float xmid, xx, etaw, width;

    float fpow();

/* For testing purposes: 
    fprintf(stderr, "RHEOL=%d\n",RHEOL);
    fprintf(stderr, "N0=%g\n",N0);
    fprintf(stderr, "E=%g\n",Eact);
    fprintf(stderr, "V=%g\n",V);
    fprintf(stderr, "T0=%g\n",T0);
    fprintf(stderr, "dN=%g\n",dN);
    fprintf(stderr, "Ttr=%g\n",Ttr);
    fprintf(stderr, "T=%g\n",T);
    fprintf(stderr, "x=%g\n",x);
    fprintf(stderr, "z=%g\n",z);
    fprintf(stderr, "eedot=%g\n",eedot);
    fprintf(stderr, "nexpo=%g\n",nexpo);
    fprintf(stderr, "compo0=%g\n",compo0);
    fprintf(stderr, "compo1=%g\n",compo1);
    fprintf(stderr, "eedot_available=%d\n",E->control.eedot_available);
*/

    if (RHEOL==0)
       if (E->control.eedot_available==0)
          visc=1.;
       else 
          visc = N0;
    else if (RHEOL==222)
       if (E->control.eedot_available==0)
          visc=compo0;
       else
          visc = N0*compo0;
    else if (RHEOL==1)
       /* Suitable for nonlin.rh.: first (lin.rh.-) it., take eta'=1 */
       if (E->viscosity.SDEPV && E->control.eedot_available==0)
          visc = 1.0;
       else 
          visc = N0 * exp(Eact/(T+T0));
    else if (RHEOL==2)
       visc = N0 * exp((Eact+z*V)/(T+T0));
    else if (RHEOL==3)
       visc = N0 * exp((Eact+z*V)/(T+T0));
    else if (RHEOL==4)
       visc = N0 * exp(Eact*(1.-T));
    else if (RHEOL==33)
       visc = N0 * (exp((Eact+z*V)/(T+T0))) * (pow(E->tracers.comporheo,compo0));
    else if (RHEOL==34)
       visc = N0 * (exp((Eact+z*V)/(T+T0))) * (pow(E->tracers.comporheo,compo1));
    else if (RHEOL==100)
    {  if (T < Ttr)
          visc = N0;
       else if (T >= Ttr)
          visc = N0*dN;
    }
    else if (RHEOL==101)
    {  visc = exp(Eact/(T+T0));
       if (T < Ttr)
          visc *= N0;
       else if (T >= Ttr)
          visc *= N0*dN;
    }
    else if (RHEOL==102)
    { /* Compositional dependent rheology from PvK-RT * benchmark:
         Use smooth geometric (or logarithmic) average */
      visc = pow(E->tracers.comporheo, compo0);
    }
    else if (RHEOL==103)
    { visc = pow(1e-4, compo1)*pow(1e2,compo0);
      /*
      fprintf(stderr,"C0=%g, C1=%g, visco=%g\n", compo0, compo1, visc);
      */
    }
    else if (RHEOL==104)
    { if (z < 0.0666) 
         visc=1e-4;
      else
         visc=1e0;
    }
    else if (RHEOL==105) /* sim. Harro Schmeling's viscosity definition 
                            for Subduction benchmark */
    { visc = pow(1e-2, compo1)*pow(1e2,compo0);
      /*
      fprintf(stderr,"C0=%g, C1=%g, visco=%g\n", compo0, compo1, visc);
      */
    }
    else if (RHEOL==111)
       /* MDB IV/2006 MR & H2O dep. rheology */
      /* H2O exhaution an depl=3% stiffens rheology by factor of 10 */
     /* melt lubricates according to Hirth&Kohlstedt and Kohlstedt&Zimmermann*/
      {
       //compo0=E->data.H2O/(0.01+compo0*0.99);
       if (compo0>0.03) compo0=10.;
       else compo0=(compo0*300.)+1.;
 
       if (E->viscosity.SDEPV==1 && E->control.eedot_available==0)
          visc = 1.0;
       else 
	 visc = N0 * exp(Eact/(T+T0)) * exp(-40.0*compo1) * compo0;
      }
    else if (RHEOL==200)
    {  /* Suitable for nonlin.rh.: first (lin.rh.-) it., take eta'=1 */
       if (E->viscosity.SDEPV && E->control.eedot_available==0)
          visc = 1.0;
       else
          visc = N0 * exp(Eact/(T+T0));
       if ((x-z)>5.0 && (x-z) < (5.0 + E->viscosity.subdweak_w) 
                     && z < 0.2) 
                       visc = min (visc, E->viscosity.subdweak_eta);
       //fprintf(stderr, " ---->visc=%g\n", visc);
    }
    else if (RHEOL==201)
    {  /* Suitable for nonlin.rh.: first (lin.rh.-) it., take eta'=1 */
       if (E->viscosity.SDEPV && E->control.eedot_available==0)
          visc = 1.0;
       else
          visc = N0 * exp(Eact/(T+T0));
       if (z<0.2)
       {  w = E->viscosity.subdweak_w;
          if (   ((x-z)>5.0       && (x-z)<(5.0+0.2*w))
              || ((x-z)>5.0+0.8*w && (x-z)<(5.0+w    )) )
                  visc = min (visc, 10*E->viscosity.subdweak_eta);
          else if ( (x-z)>5.0+0.2*w  && (x-z)<(5.0+0.8*w))
                  visc = min (visc, E->viscosity.subdweak_eta);
       }
       //fprintf(stderr, " ---->visc=%g\n", visc);
    }
    else if (RHEOL==210)
    {  /* as 200, but with location of weak zone from tracers */
       if (E->viscosity.SDEPV && E->control.eedot_available==0)
          visc = 1.0;
       else
       {   visc = N0 * exp(Eact/(T+T0));
           dzz = E->tracers.dzz;
           tracnr = floor(z/dzz+dzz/2);
           if (tracnr<E->tracers.lntrac 
              && x>(E->tracers.xtrac[1][tracnr]-E->viscosity.subdweak_w)
              && x< E->tracers.xtrac[1][tracnr]) 
           {  visc = min (visc, E->viscosity.subdweak_eta);
              /* fprintf(stderr, "(x,z=%g,%g, xtr,ztr,trnr=%g,%g,%d\n", 
                   x, z, E->tracers.xtrac[1][tracnr], 
                   E->tracers.xtrac[2][tracnr], tracnr);
               */
           }
       }
    }
    else if (RHEOL==211)
    {  /* as 210, but with smoother eta-variation*/
       if (E->viscosity.SDEPV && E->control.eedot_available==0)
          visc = 1.0;
       else
       {   visc = N0 * exp(Eact/(T+T0));
           dzz = E->tracers.dzz;
           tracnr = floor(z/dzz+dzz/2);
           width=0.5*E->viscosity.subdweak_w;
           xmid=E->tracers.xtrac[1][tracnr]-width;
           xx=max(0.0,min(1.0,2*fabs((x-xmid)/width)-1.0));
           etaw=pow(E->viscosity.subdweak_eta,(1.0-xx))*
                pow(E->viscosity.max_value,xx);
           if (tracnr<E->tracers.lntrac) 
           {  visc = min (visc, etaw);
              /*
              fprintf(stderr, 
                   "(x,z=%g,%g,xtr,ztr,trnr=%g,%g,%d,xx=%g,eta=%g\n", 
                   x, z, E->tracers.xtrac[1][tracnr], 
                   E->tracers.xtrac[2][tracnr], tracnr,
                   xx,visc);
              */
              
           }
       }
    }
    else if (RHEOL==212)
    {  /* as 211, but with constant dip (only top tracer used) */
       if (E->viscosity.SDEPV && E->control.eedot_available==0)
          visc = 1.0;
       else
       {   visc = N0 * exp(Eact/(T+T0));
           width=0.5*E->viscosity.subdweak_w;
           xmid=E->tracers.xtrac[1][0]-width+z;
           xx=max(0.0,min(1.0,2*fabs((x-xmid)/width)-1.0));
           etaw=pow(E->viscosity.subdweak_eta,(1.0-xx))*
                pow(E->viscosity.max_value,xx);
           if (z<0.15) 
           {  visc = min (visc, etaw);
              /*
              fprintf(stderr, 
                   "(x,z=%g,%g,xtr,ztr,trnr=%g,%g,%d,xx=%g,eta=%g\n", 
                   x, z, E->tracers.xtrac[1][tracnr], 
                   E->tracers.xtrac[2][tracnr], tracnr,
                   xx,visc);
              */
              
           }
       }
    }
    else 
    { printf("STOP in point_visc: Wrong RHEOL=%d", RHEOL);
      abort();
    }


    /* Nonlinear rheology part:
     * For strainr.=A tau^n exp[-E/R(1/T-1/T0)], and eta=tau/strainr
     *     eta = A^{-1/n} strainr^{(1-n)/n} exp[E/nR(1/T-1/T0)]
     *         = strainr^{(1-n)/n} {A*exp[-E/R(1/T-1/T0)]}^{-1/n}
     * For n=1, nonl. definition eta = {A*exp[-E/R(1/T-1/T0)]}^-1
     *                               = A^-1 * exp[E/R(1/T-1/T0)]
     * For lin. definition       eta = N0   * exp[E/R(1/T-1/T0)]
     * So A = 1/N0
     * Substituting gives:
     *     eta = strainr^{(1-n)/n} * 
     *              {exp[-E/R(1/T-1/T0)]/N0}^{-1/n}
     *         = strainr^{(1-n)/n} * 
     *              {N0*exp[E/R(1/T-1/T0)]}^{1/n}
     *         = strainr^{(1-n)/n} * {eta_Newt}^{1/n}
     */
    if(E->viscosity.SDEPV && E->control.eedot_available)
    {   exponent1 = 1./nexpo;
        exponent2 = exponent1 - 1.;
        visc = fpow(eedot, exponent2) * fpow(visc, exponent1);
    }   

    if(E->viscosity.SDEPV && E->control.eedot_available &&
       E->viscosity.byerlee)
    {  viscold=visc;
       /* OLD visc = min(visc, (0.1*1.36e7*z)/eedot); */
       visc = min(visc, E->viscosity.byerleemu*z/eedot);
       /*
       fprintf(stderr,"z=%g, C=%g, visc=%g, visc_y=%g, eedot=%g\n", 
              z, compo0, viscold, (0.1*1.36e7*z)/eedot, eedot);
       */
    }
    /* debugging:
    if (E->control.eedot_available)
       fprintf(stderr, "z=%g,N0=%g,T=%g,Eact=%g,e=%g,eta=%g,n=%g\n", 
                     z,N0,T,Eact,eedot,visc,nexpo);
    */

    /* Apply max.&min.-eta-limiters */
    if(E->viscosity.MAX && visc > E->viscosity.max_value)
       visc = E->viscosity.max_value;
    if(E->viscosity.MIN && visc < E->viscosity.min_value)
       visc = E->viscosity.min_value;

    return visc;
}

float point_arrh(E, RHEOL, N0, Eact, T0, T)
    struct All_variables *E;
    int RHEOL;
    float N0, Eact, T0, T;
{
    float arrh;

    /* For testing purposes: 
    fprintf(stderr, "RHEOL=%d\n",RHEOL);
    fprintf(stderr, "N0=%g\n",N0);
    fprintf(stderr, "E=%g\n",Eact);
    fprintf(stderr, "T0=%g\n",T0);
    fprintf(stderr, "T=%g\n",T);
    */

    if (RHEOL==0 || RHEOL==222 || (RHEOL>=102 && RHEOL<=105))
       arrh = 1./N0;
    else if (RHEOL==1 || RHEOL==2 || RHEOL==111)
       arrh = 1./(N0 * exp(Eact/(T+T0)));
    else if (RHEOL>=200 && RHEOL<=219)
       arrh = 1./(N0 * exp(Eact/(T+T0)));
    else 
    {  printf("STOP in point_arrh: Wrong RHEOL=%d", RHEOL);
       abort();
    }

    return arrh;
}

float fpow(base,exponent)
    /* Single precision (or float-) version of pow */
    float base,exponent;
{    
    float result;

    result = (float)(pow((double)(base),(double)(exponent)));

    return result;
}

void printarray(E, text, array, len)
    /* Debugging tools for nonlin.rheology: prints field for every
     * timestep and nonl.rh.interation
     */
    struct All_variables *E;
    char *text; 
    float *array;
    int len;
{
    char filename[200];
    FILE *fp;
    int i;
  
    if (E->control.debug > 0)
    {
       sprintf(filename,"%s.%s.%d.%d",E->control.data_file2,text,
                E->monitor.solution_cycles, E->monitor.sdepv_cycles);
       
       fp = fopen(filename,"w");
       for (i=1;i<=len;i++)
           fprintf(fp, "%d %15.10e\n", i, array[i]);
       fclose(fp);
    }

    return;
}

 int layers(E,node)
    struct All_variables *E;
    int node;
 {
    int llayers = 0;

    if (E->X[2][node] <= E->viscosity.zlith)
      llayers = 1;
    else if (E->X[2][node] >  E->viscosity.zlith &&
             E->X[2][node] <= E->viscosity.zlm)
      llayers = 2;
    else if (E->X[2][node] >  E->viscosity.zlm)
      llayers = 3;

    return (llayers);
    }

float boundary_thickness(E,H)
   struct All_variables *E;
   float *H;
{
   float thickness;
   int i,j;
  
   thickness = 0.0;
  
   for (i=1;i<=E->lmesh.noz;i++) 
      H[i] = -H[i]/E->control.Atemp;
   
   if (E->parallel.me_loc[1]==0)  
   {  for (i=1;i<=E->lmesh.noz;i++) 
      {  if (H[i]>0.45)
         {  j = i;
            break;
         }
      }
      thickness=E->X[2][j];
   }
   
   MPI_Bcast(&thickness,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  
   return (thickness);
}
