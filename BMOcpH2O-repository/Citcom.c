/*
This 2D/3D cartesian version was original from L. Moresi.
Shijie Zhong's contributions are:
1) rewrote the code to include parallel computing with all the MPI
calls.
2) introduced the full multigrid algorithm and consistent projection
scheme.
3) re-ordered stiffness matrix and its matrix operation such that
   only half of the stiffness matrix is stored.

If using this code results in any publication, I suggest that
you reference the following two papers for the method:

Zhong, S., M. T. Zuber, L. Moresi, and M. Gurnis, Role
           of temperature-dependent viscosity and surface plates in
           spherical shell models of mantle convection,
           J. Geophys. Res., 105, 11063-11082, 2000.

and also Moresi's original paper.

Moresi L.N. and M. Gurnis, Constraints on lateral strength of slabs
           from 3-D dynamic flow models, EPSL, 138, 15-28, 1996.
*/

#if 1
#include <mpi.h>
#endif

#include <math.h>
#include <malloc.h>
#include <sys/types.h>

#include "element_definitions.h"
#include "global_defs.h"

extern int Emergency_stop;

main(argc,argv)
     int argc;
     char **argv;
     
{	/* Functions called by main*/
  void general_stokes_solver();
  void read_instructions();
  void solve_constrained_flow();
  void solve_derived_velocities();
  void parallel_process_initilization();
  void parallel_process_termination();
  void process_temp_field();
  void vcopy();
  void tracadv();
  void updatetrac();
  void tracinterp_t2e();
  void tracinterp_e2t();
  void time_display();  

  float dot();

  int k, *temp;
  double CPU_time0(),time,initial_time,start_time;
  double *p_init_time;

  struct All_variables E;

/*
  parallel_process_initilization(&E,argc,argv); 
*/

  E.parallel.me = 0;
  E.parallel.nproc = 1;
  E.parallel.me_loc[1] = 0;
  E.parallel.me_loc[2] = 0;
  E.parallel.me_loc[3] = 0;

  /* flag to tell whether tracers are (already) initialized */
  E.tracers.initialized=0;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &(E.parallel.me) );
  MPI_Comm_size(MPI_COMM_WORLD, &(E.parallel.nproc) );

  gethostname ( E.parallel.machinename, 160); 

 E.monitor.solution_cycles=0;
 E.control.eedot_available=0;

 if (E.parallel.me == 0) {
    start_time = time = CPU_time0();
    }
 
  if (E.parallel.me == 0) fprintf(stderr,"Before read_instructions\n");
  read_instructions(&E,argc,argv);
  if (E.parallel.me == 0) fprintf(stderr,"After read_instructions\n");

  E.control.keep_going=1;

  if (E.parallel.me == 0)  {
     fprintf(stderr,"Input parameters taken from file '%s'\n",argv[1]);
     fprintf(stderr,"Initialization complete after %g seconds\n\n",
                     CPU_time0()-time); fflush(E.fp);
     initial_time = CPU_time0()-time;
     fprintf(E.fp,"Initialization overhead = %f\n",initial_time);
     initial_time = CPU_time0();
     E.control.initial_time = initial_time;
     }
  p_init_time = &initial_time;

  general_stokes_solver(&E);
  process_new_velocity(&E,E.monitor.solution_cycles);

  if (E.control.stokes)  {
    E.control.keep_going = 0;
    E.monitor.solution_cycles++;
    }

 
  while ( E.control.keep_going   &&  (Emergency_stop == 0) )   
  {
     E.monitor.solution_cycles++; 
     if(E.monitor.solution_cycles>E.control.print_convergence)
     E.control.print_convergence=1;
     time_display(&E,p_init_time, 1);
     /**/ report(&E,"Update buoyancy for further `timesteps'");
     (E.next_buoyancy_field)(&E);
     time_display(&E,p_init_time, 2);
     /**/ report(&E,"Process results of buoyancy update");
     process_temp_field(&E,E.monitor.solution_cycles); 
     time_display(&E,p_init_time, 3); 
     /* Tracer advection */
     if (E.tracers.tracers_used!=0)
      { 
	tracadv(&E,1);
        time_display(&E,p_init_time, 4);
        tracinterp_e2t(&E);
        time_display(&E,p_init_time, 5);
        updatetrac(&E,1);
        time_display(&E,p_init_time, 6);
        tracinterp_t2e(&E);
        if (E.tracers.tractest > 0) /* Print tracer coords */
	   printtrac(&E);
        time_display(&E,p_init_time, 7);
      }

     general_stokes_solver(&E);
     time_display(&E,p_init_time, 8);
     /**/ report(&E,"Process results of velocity solver");
     process_new_velocity(&E,E.monitor.solution_cycles);
  }
  
  if (E.parallel.me == 0)  {
     time=CPU_time0()-initial_time;
     fprintf(E.fp,"Average cpu time taken for velocity step = %f\n",
	 time/((float)(E.monitor.solution_cycles)));
     fprintf(stderr,"Average cpu time taken for velocity step = %f\n",
	 time/((float)(E.monitor.solution_cycles)));
     }

  fclose(E.fp);

  parallel_process_termination();

  return;  

} 
