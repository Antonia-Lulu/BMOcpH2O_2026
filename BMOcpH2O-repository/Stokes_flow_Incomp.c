/*   Functions which solve for the velocity and pressure fields using Uzawa-type iteration loop.  */

#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

extern int Emergency_stop;

/* Master loop for pressure and (hence) velocity field */


void solve_constrained_flow_iterative(E)
     struct All_variables *E;

{ 
    double *D1;
    double *u;
    double *R,*Bp;
    double residual_ddash;
    double vmag;
    double global_vdot(),global_pdot();

    static int been_here = 0;
   
    float solve_Ahat_p_fhat();
    void assemble_del2_u();
    void assemble_grad_p();
    void assemble_div_u();
    void v_from_vector();
    void p_to_nodes();   
    void strip_bcs_from_residual();
    void velocities_conform_bcs();
   
    int steps,cycles;
    int i,j,k,doff,vel_cycles_previous,vel_calls_previous;
  
    double time,CPU_time0();
 
    const int npno = E->lmesh.npno;
    const int gnpno = E->mesh.npno;
    const int nno = E->lmesh.nno;
    const int dims = E->mesh.nsd;
    const int neq = E->lmesh.neq;
    const int gneq = E->mesh.neq;

    time=CPU_time0();
    
    cycles=E->control.p_iterations;
   
    /* Solve for velocity and pressure, correct for bc's */

    residual_ddash=solve_Ahat_p_fhat(E,E->U,E->P,E->F,E->control.accuracy,&cycles);

    been_here=1;

    v_from_vector(E,E->V,E->U); 
   /* p_to_nodes(E,E->P,E->NP,E->mesh.levmax);   

 */

  return; 
}



/*  ===============================================================  */

float solve_Ahat_p_fhat(E,V,P,F,imp,steps_max)

     struct All_variables *E;
     double *V,*P,*F;
     double imp;
     int *steps_max;
    
{ 
  int i,j,k,ii,count,convergent,valid,problems;
  int lev,lev_low,npno,neq,steps;
  int gnpno,gneq;
  
  static int been_here = 0;
  double *u;
  static double *p1,*r1,*r0,*r2,*z0,*z1,*s1,*s2,*Ah,*u1;
  double *shuffle, *R;
  double alpha,delta,s2dotAhat,r0dotr0,r1dotz1;
  double residual, initial_residual, last_residual,res_magnitude,v_res;
  
  double global_vdot(),global_pdot();
  double *dvector();

  float CPU_time();
  double time0,time,CPU_time0();
  static double timea;
  float dpressure,dvelocity,tole_comp;

  int imax, eqmax, dofmax, nodemax, procmax, elmax;
  double u1max, s2max; 
  float x, y, z;
  int global_vmax(), global_pmax();

  void assemble_div_u();
  void assemble_del2_u();
  void assemble_grad_p();
  void strip_bcs_from_residual();
  int  solve_del2_u();
   
  const int dims=E->mesh.nsd;
  const int n=loc_mat_size[E->mesh.nsd];

  npno=E->lmesh.npno;
  neq=E->lmesh.neq;

  gnpno=E->mesh.npno;
  gneq=E->mesh.neq;

  if (been_here==0)  {
    r0 = (double *)malloc((npno+1)*sizeof(double));
    r1 = (double *)malloc((npno+1)*sizeof(double));
    r2 = (double *)malloc((npno+1)*sizeof(double));
    z0 = (double *)malloc((npno+1)*sizeof(double));
    z1 = (double *)malloc((npno+1)*sizeof(double));
    s1 = (double *)malloc((npno+1)*sizeof(double));
    s2 = (double *)malloc((npno+1)*sizeof(double));
    Ah = (double *)malloc((neq+1)*sizeof(double));
    u1 = (double *)malloc((neq+1)*sizeof(double));
    been_here=1;
    timea=CPU_time0();
    }

  problems=0;
  time0=time=CPU_time0();

  been_here ++;
 
  /* calc. the velocity residual, note there are tricks involved here */

  lev=E->mesh.levmax;

  v_res=sqrt(global_vdot(E,F,F,lev)/gneq);

  if (E->parallel.me==0)
     fprintf(stderr,
        "initial residue of momentum equation %g %d\n",v_res,gneq);

  assemble_grad_p(E,P,Ah,lev);
  assemble_del2_u(E,V,u1,lev,1);

  for(i=0;i<neq;i++) 
      Ah[i] = F[i] - Ah[i] - u1[i]; 

  strip_bcs_from_residual(E,Ah,lev);

  valid=solve_del2_u(E,u1,Ah,imp*v_res,E->mesh.levmax);
  strip_bcs_from_residual(E,u1,lev);

  if(!valid) problems++; 
  
  for(i=0;i<neq;i++) 
      V[i] += u1[i];

  assemble_div_u(E,V,r1,lev);

  residual = initial_residual = sqrt(global_pdot(E,r1,r1,lev)/gnpno);

  E->monitor.vdotv = sqrt(global_vdot(E,V,V,lev)/gneq);

  E->monitor.incompressibility = residual/E->monitor.vdotv;
         
  if(E->control.print_convergence && E->parallel.me==0)
      fprintf(stderr,"Loop to reduce pressure residual %g\n",residual);
     
   count = 0;
   convergent=0;

   dpressure = 1.0;
   dvelocity = 1.0;

   res_magnitude = residual;

   if (E->control.print_convergence && E->parallel.me==0)  {
      fprintf(E->fp,
       "AhatP (%03d) after %g sec %g sec with div/v=%.3e for step %d\n",
       count,CPU_time0()-time0,CPU_time0()-timea,
       E->monitor.incompressibility,E->monitor.solution_cycles); /**/
      fprintf(stderr,
       "AhatP (%03d) after %g sec %g sec with div/v=%.3e for step %d\n",
       count,CPU_time0()-time0,CPU_time0()-timea,
       E->monitor.incompressibility,E->monitor.solution_cycles); /**/
         }         

   /* while( (count < *steps_max) && 
    *       (E->monitor.incompressibility >= E->control.tole_comp || 
    *        dvelocity >= imp) 
    *     )  {     
    */
   while( (count < *steps_max) && 
          (dpressure >= imp || dvelocity >= imp) 
        )  {

        for(j=1;j<=npno;j++)
           z1[j] = E->BPI[lev][j]*r1[j];
     
        r1dotz1 = global_pdot(E,r1,z1,lev);

        if ((count == 0))
 	   for(j=1;j<=npno;j++)
	      s2[j] = z1[j];
        else {
	   r0dotr0=global_pdot(E,r0,z0,lev);
	   assert(r0dotr0 != 0.0  /* Division by zero in head of 
                                     incompressibility iteration */);
	   delta = r1dotz1/r0dotr0;
	   for(j=1;j<=npno;j++)
	      s2[j] = z1[j] + delta * s1[j];
           }
      
        assemble_grad_p(E,s2,Ah,lev); 

        valid=solve_del2_u(E,u1,Ah,imp*v_res,lev);  
        strip_bcs_from_residual(E,u1,lev);

        if(!valid) problems++;
      
        assemble_div_u(E,u1,Ah,lev);

        s2dotAhat=global_pdot(E,s2,Ah,lev);
      
        if(s2dotAhat == 0.0)  {
	  fprintf(E->fp,
            "The pressure residual cannot be improved any further.\n");
	  fprintf(E->fp,
          "To improve this sol.,reduce tolerance of velocity solver\n");
	  alpha = 0.0;
          }
        else {
	  /* alpha defined this way is the same as R&W */
	  alpha = r1dotz1/s2dotAhat; 
	  // fprintf(stderr,"r1dotz1/s2dotAhat: %f/%f\n",r1dotz1,s2dotAhat);
          }   
     
        for(j=1;j<=npno;j++)   {
 	  r2[j] = r1[j] - alpha * Ah[j];
	  P[j] += alpha * s2[j];
          }
     
        for(j=0;j<neq;j++)
	  V[j] -= alpha * u1[j];
      
        assemble_div_u(E,V,Ah,lev);
        E->monitor.vdotv = global_vdot(E,V,V,E->mesh.levmax);
        E->monitor.incompressibility = 
           sqrt( (gneq/gnpno)*
                 (1.0e-32+global_pdot(E,Ah,Ah,lev)/
                         (1.0e-32+E->monitor.vdotv) ) );
        dpressure=alpha*
           sqrt(global_pdot(E,s2,s2,lev)/
                (1.0e-32+global_pdot(E,P,P,lev)) );
        dvelocity=alpha*
           sqrt(global_vdot(E,u1,u1,lev)/(1.0e-32+E->monitor.vdotv));
        
         /*  Locate the location of the worst convergence  */
        if (count > 10) 
        {  imax = global_vmax( E, u1, &u1max, lev );
           neq = E->lmesh.NEQ[lev];
           procmax = imax/neq;
           eqmax = imax%neq;
           dofmax = (eqmax%E->mesh.nsd)+1;
           nodemax = eqmax/E->mesh.nsd;
           if ( E->parallel.me==procmax )
           {  x = E->XX[lev][1][nodemax];
              z = E->XX[lev][2][nodemax];
              if(E->mesh.nsd==3)
              {  y = E->XX[lev][3][nodemax];
                 fprintf(stderr,
                 "  Max.dv=%g,dof=%d,node=%d,proc=%d,(x,y,z)=(%g,%g,%g)\n",
                      u1max,dofmax, nodemax, procmax, x, y, z);
              }
              else 
                 fprintf(stderr,
                 "  Max.dv=%g,dof=%d,node=%d,proc=%d,(x,z)=(%g,%g)\n",
                      u1max,dofmax, nodemax, procmax, x, z);
           }
         
           imax = global_pmax( E, s2, &s2max, lev );
           npno=E->lmesh.NPNO[lev];
           procmax = imax/npno;
           elmax = imax%npno;
           if ( E->parallel.me==procmax )
           {  x = E->ECO[lev][elmax].centre[1];
              z = E->ECO[lev][elmax].centre[2];
              if(E->mesh.nsd==3)
              {  y = E->ECO[lev][elmax].centre[3];
                fprintf(stderr,
                  "  Max.dp=%g,node=%d,proc=%d,(x,y,z)=(%g,%g,%g)\n",
                  s2max, elmax, procmax, x, y, z);
              }
              else 
                fprintf(stderr,
                  "  Max.dp=%g,node=%d,proc=%d,(x,z)=(%g,%g)\n",
                  s2max, elmax, procmax, x, z);
           }
        }
       
         
        count++;
        if (E->control.print_convergence && E->parallel.me==0)  {
         fprintf(E->fp,
             "AhatP (%03d) div/v=%.3e, dv/v=%.3e, dp/p=%.3e, step=%d\n",
             count, E->monitor.incompressibility,
             dvelocity, dpressure, E->monitor.solution_cycles); 
         fprintf(stderr,
             "AhatP (%03d) div/v=%.3e, dv/v=%.3e, dp/p=%.3e, step=%d\n",
             count, E->monitor.incompressibility,
             dvelocity, dpressure, E->monitor.solution_cycles); 
          fflush(E->fp);
          }         

        shuffle=s1;s1=s2;s2=shuffle;
        shuffle=r0;r0=r1;r1=r2;r2=shuffle;
        shuffle=z0;z0=z1;z1=shuffle;
      
        }       /* end loop for conjugate gradient   */

    if(problems) {
      fprintf(E->fp,
         "Convergence of velocity solver may affect continuity\n");
      fprintf(E->fp,
         "Consider running with the `see_convergence=on' option\n");
      fprintf(E->fp,
         "To evaluate performance of current relaxation parameters\n");
      fflush(E->fp);
      }

        if (E->control.print_convergence && E->parallel.me==0)  {
         fprintf(E->fp,
            "after (%03d) pressure loops and %g sec for step %d\n",
            count,CPU_time0()-timea,E->monitor.solution_cycles); /**/
         fprintf(stderr,
            "after (%03d) pressure loops and %g sec for step %d\n",
            count,CPU_time0()-timea,E->monitor.solution_cycles); /**/
         fflush(E->fp);
         }         
   
    
/*
  free((void *) r0);
  free((void *) r1);       
  free((void *) r2);
  free((void *) z0);
  free((void *) z1);
  free((void *) s1);
  free((void *) s2);
  free((void *) Ah);
  free((void *) u1);
*/

    *steps_max=count;

    return(residual);
 }

/*  ===============================================================  */

void v_from_vector(E,V,F)
     struct All_variables *E;
     float **V;
     double *F;
{
  int node,d;
  unsigned int type;

  const int nno = E->lmesh.nno;
  const int dofs = E->mesh.dof;

  for(node=1;node<=nno;node++)     {
      V[1][node] = F[E->id[node].doff[1]]; 
      V[2][node] = F[E->id[node].doff[2]]; 
      if(dofs==3) V[3][node] = F[E->id[node].doff[3]];
      if (E->node[node] & VBX)
             V[1][node] = E->VB[1][node]; 
      if (E->node[node] & VBZ)
             V[2][node] = E->VB[2][node]; 
      if (dofs==3 && E->node[node] & VBY)
             V[3][node] = E->VB[3][node]; 
    }
  return;
}
