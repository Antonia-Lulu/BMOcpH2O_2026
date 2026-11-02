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
#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

void general_stokes_solver(E)
     struct All_variables *E;

{
    void construct_stiffness_B_matrix();
    void velocities_conform_bcs();
    void assemble_forces();
    double global_vdot();
    void parallel_process_termination();

    FILE *fp;
    char filename[100];

    float vmag, radius, dx1, dx2, dx3, angle;
    double *delta_U;
    double *force,  Udot_mag, dUdot_mag;
    double CPU_time0(),time;
    int i,j,k;
    // static double time_saver;

    //static double *p_init_time;
    static double *oldU;
    static int visits=0;

    const int nno = E->lmesh.nno;
    const int neq = E->lmesh.neq;
    const int vpts = vpoints[E->mesh.nsd];
    const int dims = E->mesh.nsd;
    const float PI=3.141592635;

    if(visits==0) {
	//p_init_time=&time_saver;
	oldU = (double *)malloc((neq+2)*sizeof(double));
	for(i=0;i<neq;i++) 
	    oldU[i]=0.0;
    }
    visits ++;

    dUdot_mag = 0.0;

    delta_U = (double *)malloc((neq+2)*sizeof(double));
     
    E->monitor.elapsed_time_vsoln1 =  E->monitor.elapsed_time_vsoln;
    E->monitor.elapsed_time_vsoln = E->monitor.elapsed_time;

    if(E->parallel.me==0) time=CPU_time0();

    E->monitor.sdepv_cycles=0;

    /* Solve Stokes 1st time a little more accurate */
    if (visits == 1)
    {   E->viscosity.sdepv_misfit*=0.2;
        E->control.accuracy*=0.2;
	if (E->control.shear_heating>0) 
	  {E->control.accuracy*=0.01;
           E->viscosity.sdepv_misfit*=0.01;}
    }
    else if (visits == 2)
    {   E->viscosity.sdepv_misfit*=5.;
        E->control.accuracy*=5.;
	if (E->control.shear_heating>0)
	  {E->control.accuracy*=100.;
	   E->viscosity.sdepv_misfit*=100.;}
    }
          
    do  
    {  //time_display(E,p_init_time, 100);
       if(E->viscosity.update_allowed)
          get_system_viscosity(E,1,E->EVI[E->mesh.levmax],
                               E->VI[E->mesh.levmax]);
       //time_display(E,p_init_time, 101);
       if (E->control.testvel == 0)
       { /* Solve Stokes, no analytical testvelocity */
          velocities_conform_bcs(E,E->U);
          //time_display(E,p_init_time, 102);
          assemble_forces(E,0); 
          //time_display(E,p_init_time, 103);
          construct_stiffness_B_matrix(E);
	  //time_display(E,p_init_time, 104);
          solve_constrained_flow_iterative(E);	
          //time_display(E,p_init_time, 105);
          /* Now v, and so also strainrate is available */
          E->control.eedot_available=1;

          if ( E->viscosity.SDEPV )   
          {  for (i=0;i<neq;i++) 
             {  delta_U[i] = E->U[i] - oldU[i]; 
 	        oldU[i] = E->U[i];
 	     }
             Udot_mag  = sqrt(global_vdot(E,E->U,E->U,E->mesh.levmax));
             dUdot_mag = sqrt(global_vdot(E,delta_U,delta_U,
                                          E->mesh.levmax)); 
             /* fprintf(stderr,"Udot_mag=%g, dUdot_mag=%g\n", Udot_mag,
                      dUdot_mag); */
   
             if (Udot_mag !=0.0)
	          dUdot_mag /= Udot_mag;

             if(E->parallel.me==0)
             {  fprintf(stderr,
                        "SDEPV: dUdot/Udot=%.4e (Udot=%.4e) for it. %d\n",
                        dUdot_mag,Udot_mag,E->monitor.sdepv_cycles);
	        fprintf(E->fp,
                        "SDEPV: dUdot/Udot=%.4e (Udot=%.4e) for it. %d\n",
                        dUdot_mag,Udot_mag,E->monitor.sdepv_cycles);
                fflush(E->fp); 
	     }
	     E->monitor.sdepv_cycles++;
          }         /* end for SDEPV   */
       } 
       else if (E->control.testvel == 1)
       { /* Analytical velocity field 1: solid clockwise rotation 
          * around O'=(x,z)=(0.5,0.5) in 2-D
          * and omega'=1000 (i.e. 1 revolution in t'=0.001)
          */
          for (i=1;i<=E->lmesh.nno;i++)
          {  dx1 = E->X[1][i]-0.5;
             dx2 = E->X[2][i]-0.5;
             radius=sqrt(dx1*dx1+dx2*dx2);
             angle = atan2(dx2,dx1); /* Clockwise angle
                                        O'->(x1,x2) with pos. x-axis */
             E->V[1][i]=-2000.*PI*radius*sin(angle);
             E->V[2][i]=2000.*PI*radius*cos(angle);
             if (E->mesh.nsd == 3)
                E->V[3][i]=0.;
             if (visits==0)
                fprintf(E->fp, "At (%g,%g), angle=%g, v=(%g,%g)\n",
                            E->X[1][i], E->X[2][i],
                            angle,
                            E->V[1][i], E->V[2][i]);
             
          }
       }      
       else if (E->control.testvel == 2)
       { /* Analytical velocity field 1: solid clockwise rotation 
          * around O'=(x,z,y)=(0.5,z,0.5) in 3-D
          * and omega'=1000 (i.e. 1 revolution in t'=0.001)
          */
          for (i=1;i<=E->lmesh.nno;i++)
          {  dx1 = E->X[1][i]-0.5;
             dx2 = E->X[3][i]-0.5;
             radius=sqrt(dx1*dx1+dx2*dx2);
             angle = atan2(dx2,dx1); /* Clockwise angle
                                        O'->(x1,x2) with pos. x-axis */
             E->V[1][i]=-2000.*PI*radius*sin(angle);
             E->V[2][i]=0.;
             E->V[3][i]=2000.*PI*radius*cos(angle);
             /* fprintf(stderr, "At (%g,%g), angle=%g, v=(%g,%g)\n",
                            E->X[1][i], E->X[2][i],
                            angle,
                            E->V[1][i], E->V[2][i]);
              */
          }
       }      
       else if (E->control.testvel == 3)
       { /* No flow: vx=vy=vz=0 */
          for (i=1;i<=E->lmesh.nno;i++)
          {  E->V[1][i]=0e0;
             E->V[2][i]=0e0;
             if (E->mesh.nsd == 3)
                E->V[3][i]=0.;
          }
       }      

       else if (E->control.testvel == 4)
       { /* Diagonal flow 1 in 2-D */
          for (i=1;i<=E->lmesh.nno;i++)
          {  E->V[1][i]=E->control.testvelval[1];
             E->V[2][i]=E->control.testvelval[2];
             if (E->mesh.nsd == 3)
                E->V[3][i]=E->control.testvelval[3];
          }
       } 
       else if (E->control.testvel == 5)
       { /* accelerated horizontal v-field */
          for (i=1;i<=E->lmesh.nno;i++)
	    {  E->V[1][i]=260000000*(E->monitor.elapsed_time); // (+ E->advection.timestep) in 1st PC-iteration !!!
             E->V[2][i]=0.;
             if (E->mesh.nsd == 3)
                E->V[3][i]=0.;
          }
       } 
       else if (E->control.testvel == 6)
       { /* horizontal couette flow */
          for (i=1;i<=E->lmesh.nno;i++)
	    {  E->V[1][i]=E->control.testvelval[1]*(1-E->X[1][i]);
             E->V[2][i]=0.;
             if (E->mesh.nsd == 3)
                E->V[3][i]=0.;	 
       }}

       
    } while( (E->monitor.sdepv_cycles < 50)  
             && (dUdot_mag>E->viscosity.sdepv_misfit) 
             && E->viscosity.SDEPV);
    	
    free((void *) delta_U);
      
  return;
}
