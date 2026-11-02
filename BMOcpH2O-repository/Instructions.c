/* Problem setup: allocate memory and establish mesh, T, viscosity, boundary and initial conditions. */

#include <signal.h>
#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include <string.h>
#include <stdlib.h>
#include "element_definitions.h"
#include "global_defs.h"

int Emergency_stop;

void read_instructions(E,argc,argv)
     struct All_variables *E;
     int argc;
     char **argv;
{
    int get_process_identifier();

    void allocate_common_vars();  
    //void common_initial_fields();
    void initial_pressure();
    void initial_velocity();
    void read_viscosity_option();
    void read_initial_settings();
    void global_default_values();
    void global_derived_values();
    void construct_masks();
    void construct_shape_functions();
    void construct_id();
    void construct_lm();
    void construct_sub_element();
    void mass_matrix();
    void construct_node_ks();
    void construct_node_maps();
    void construct_ien();
    void construct_mat_group();
    void interuption();
    void set_up_nonmg_aliases();
    void check_bc_consistency();
    void node_locations();
    void parallel_domain_decomp1();
    void parallel_shuffle_ele_and_id();
    void parallel_communication_routs();
    void parallel_process_termination();
    void setup_plume_problem(); 
    void init_tracers(); 
    void init_tracers_preliminary();
  
    void setup_parser();
  
    double start_time, CPU_time0();

    int *temp, i;
 
    /* Global interruption handling */

    if (E->parallel.me==0) 
	{start_time=CPU_time0();
    E->control.show_time=1;}
    Emergency_stop = 0;
    signal(SIGINT,interuption);
    signal(SIGTERM,interuption);

    E->control.PID=get_process_identifier(); 

    setup_parser(E,argc,argv);

    global_default_values(E); 
    read_initial_settings(E);  
    (E->problem_derived_values)(E); /* call before global_derived_ */
    global_derived_values(E);

    parallel_domain_decomp1(E);

    allocate_common_vars(E);                 
    (E->problem_allocate_vars)(E);
    (E->solver_allocate_vars)(E);

    construct_ien(E);
    construct_masks(E);		/* order is important here */
    construct_id(E);
    construct_lm(E);
    construct_sub_element(E);

    node_locations (E); 
    construct_mat_group(E);

    init_tracers_preliminary(E); 

    /* Boundary conditions */
    setup_plume_problem (E); 
    (E->problem_boundary_conds)(E);
    check_bc_consistency(E);

    parallel_shuffle_ele_and_id(E);
    parallel_communication_routs(E);

    construct_shape_functions(E);
    mass_matrix(E);

    /* Initial conditions */
    (E->problem_initial_fields)(E);  /* T/C/melting etc */

/* Initialize problem-independent fields; viscosity is done last because it may
   depend on all other fields. */

    report(E,"Initialize pressure field");
    initial_pressure(E);
    report(E,"Initialize velocity field");
    initial_velocity(E);

    init_tracers(E);

    report(E,"Initialize viscosity field");
    get_viscosity_option(E);

    /* ========================================================= */

    shutdown_parser(E);

 
/*if (E->parallel.me==0)
  fprintf (stderr,"done instruction\n");
    parallel_process_termination();

  exit(8);
*/
    return;
}



void allocate_common_vars(E) 
     struct All_variables *E;

{ 
    double **dmatrix();
    float **fmatrix();
    void set_up_nonmg_aliases();
    int nox,noy,noz,i,j,l,nno_l,npno_l,nozl,nnov_l,nxyz;

    E->mesh.fnodal_malloc_size = (E->lmesh.nno+2)*sizeof(float);
    E->mesh.dnodal_malloc_size = (E->lmesh.nno+2)*sizeof(double);
    E->mesh.feqn_malloc_size = (E->mesh.nsd*E->lmesh.nno+2)*sizeof(float);
    E->mesh.deqn_malloc_size = (E->mesh.nsd*E->lmesh.nno+2)*sizeof(double);
  
    E->P = (double *) malloc ((E->lmesh.npno+1)*sizeof(double));
    E->S = (double *) malloc ((E->lmesh.npno+1)*sizeof(double));
    E->eedot = (float *) malloc((E->lmesh.nel+1)*sizeof(float));
    E->visc_diss = (float *) malloc((E->lmesh.nel+1)*sizeof(float));

 
    /* HAVE */
    E->Have.T = (float *) malloc((E->mesh.noz+1)*sizeof(float));
    E->Have.Vi = (float *) malloc((E->mesh.noz+1)*sizeof(float));
    E->Have.Rho = (float *) malloc((E->mesh.noz+1)*sizeof(float));
    E->Have.vrms = (float *) malloc((E->mesh.noz+1)*sizeof(float));

    /* T_y is used to calc. hor. average in y-dir for onsettime ssc 
       The array has a local (l) and global (g) version */
    E->Have.lT_y = (float *) malloc((E->mesh.nox*E->mesh.noz+1)
		    *sizeof(float));
    E->Have.gT_y = (float *) malloc((E->mesh.nox*E->mesh.noz+1)
		    *sizeof(float));

    E->Have.arrh = (float *) malloc((E->mesh.noz+1)*sizeof(float));
    E->Have.Vy = (float *) malloc((E->mesh.noz+1)*sizeof(float));
    E->Have.edot31 = (float *) malloc((E->mesh.noz)*sizeof(float));

    /* T profiles at in and outflow boundaries */
    E->segment.Tz_in=(float *) malloc((E->lmesh.noz+1)*sizeof(float));
    E->segment.Tzy_in=(float *) 
            malloc((E->lmesh.noz+1)*(E->lmesh.noy+1)*sizeof(float));
    E->segment.Tz_out=(float *) malloc((E->lmesh.noz+1)*sizeof(float));
    E->segment.Vx = (float *) malloc((E->lmesh.noz+1)*sizeof(float));
    E->segment.Vx_zy=(float *) 
            malloc((E->lmesh.noz+1)*(E->lmesh.noy+1)*sizeof(float));
    
    E->F = (double *) malloc((E->mesh.nsd*E->lmesh.nnov+1)*sizeof(double));
    E->U = (double *) malloc((E->mesh.nsd*E->lmesh.nnov+1)*sizeof(double));
    E->T = (float *) malloc((E->lmesh.nno+1)*sizeof(float));
    E->Tstar = (float *) malloc((E->lmesh.nno+1)*sizeof(float));
    E->Qmelt = (float *) malloc((E->lmesh.nel+1)*sizeof(float));
    E->buoyancy = (float *) malloc((E->lmesh.nno+1)*sizeof(float));
    E->NP = (float *) malloc((E->lmesh.nno+1)*sizeof(float));
    E->edot = (float *) malloc((E->lmesh.nno+1)*sizeof(float));
 
  
  for(i=1;i<=E->mesh.nsd;i++)  {
      E->TB[i] = (float *)  malloc((E->lmesh.nno+1)*sizeof(float));
      E->V[i]  = (float *)  malloc((E->lmesh.nnov+1)*sizeof(float));
      E->V1[i]  = (float *)  malloc((E->lmesh.nnov+1)*sizeof(float));
      E->Vest[i]  = (float *)  malloc((E->lmesh.nnov+1)*sizeof(float));
      E->VB[i] = (float *)  malloc((E->lmesh.nnov+1)*sizeof(float));
      E->Vstar[i] = (float *)  malloc((E->lmesh.nnov+1)*sizeof(float));
      }

    E->stress = (float *) malloc((E->lmesh.nsf*12+12)*sizeof(float));
    E->slice.tpg = (float *) malloc((E->lmesh.nsf+2)*sizeof(float));
    E->slice.tpgb = (float *) malloc((E->lmesh.nsf+2)*sizeof(float));
    E->slice.vline = (float *) malloc((E->lmesh.nsf+2)*sizeof(float));
    E->slice.vlinek = (float *) malloc((E->lmesh.nsf+2)*sizeof(float));
    E->slice.shflux = (float *) malloc((E->lmesh.nsf+2)*sizeof(float));
    E->slice.bhflux = (float *)malloc((E->lmesh.nsf+2)*sizeof(float));
    E->slice.cen_mflux = (float *)malloc((E->lmesh.nsf+2)*sizeof(float));
    E->slice.vxsurf[1] = (float *)malloc((E->lmesh.nsf+2)*sizeof(float));
    E->slice.vxsurf[2] = (float *)malloc((E->lmesh.nsf+2)*sizeof(float));
   
   E->mat = (int *) malloc((E->lmesh.nel+2)*sizeof(int));

/* set up memory for different grids  */
for(i=E->mesh.levmin;i<=E->mesh.levmax;i++)   {
    for(j=1;j<=E->mesh.nsd;j++)
      E->XX[i][j] = (float *)  malloc((E->lmesh.NNO[i]+1)*sizeof(float));
    E->MASS[i]     = (float *) malloc((E->lmesh.NNO[i]+1)*sizeof(float));

    E->ECO[i] = (struct COORD *) malloc((E->lmesh.NNO[i]+2)*sizeof(struct COORD));
    E->IEN[i] = (struct IEN *) malloc((E->lmesh.NNO[i]+2)*sizeof(struct IEN));
    E->ID[i]  = (struct ID *) malloc((E->lmesh.NNO[i]+2)*sizeof(struct ID));
    E->GNX[i] = (struct Shape_function_dx *) malloc((E->lmesh.NEL[i]+2)*sizeof(struct Shape_function_dx));
    E->GDA[i] = (struct Shape_function_dA *) malloc((E->lmesh.NEL[i]+2)*sizeof(struct Shape_function_dA));
    E->EL[i]  = (struct SUBEL *) malloc((E->lmesh.NEL[i]+2)*sizeof(struct SUBEL));
    E->LMD[i]  = (struct LM *) malloc((E->lmesh.NEL[i]+2)*sizeof(struct LM));

    E->EVI[i]  = (float *) malloc((E->lmesh.NEL[i]+2)*vpoints[E->mesh.nsd]*sizeof(float));

    E->TW[i]  = (float  *) malloc((E->lmesh.NNO[i]+2)*sizeof(float));
    E->VI[i]  = (float *) malloc((E->lmesh.NNO[i]+2)*sizeof(float));
    E->NODE[i] = (unsigned int *) malloc((E->lmesh.NNO[i]+2)*sizeof(unsigned int));
    E->TWW[i] =(struct FNODE *) malloc((E->lmesh.NEL[i]+2)*sizeof(struct FNODE));
    
    E->NEI[i].nels  = (int *) malloc((E->lmesh.NNO[i]+2)*sizeof(int));
    E->NEI[i].lnode  = (int *) malloc((E->lmesh.NNO[i]+2)*enodes[E->mesh.nsd]*sizeof(int));
    E->NEI[i].element  = (int *) malloc((E->lmesh.NNO[i]+2)*enodes[E->mesh.nsd]*sizeof(int));   

    E->elt_del[i]=(struct EG *)malloc((E->lmesh.NEL[i]+1)*sizeof(struct EG));

    E->BI[i] = (double *) malloc((E->lmesh.NEQ[i]+2)*sizeof(double)); 
    E->BPI[i] = (double *) malloc((E->lmesh.NPNO[i]+1)*sizeof(double));
    E->control.B_is_good[i] = 0;
    }
    E->temp = (double *) malloc((E->lmesh.NEQ[E->mesh.levmax]+2)*sizeof(double)); 


for(i=E->mesh.levmin;i<=E->mesh.levmax;i++)  {
      nox = E->lmesh.NOX[i];
      noy = E->lmesh.NOY[i];
      noz = E->lmesh.NOZ[i];
      if (E->mesh.nsd==2)  {
        nxyz = max(nox,noz);
	}
      else if (E->mesh.nsd==3)  {
        nxyz = max(nox*noz,nox*noy);
        nxyz = max(nxyz,noz*noy);
	}

    E->parallel.IDD[i] = (int *) malloc((E->lmesh.NEQ[i]+2)*sizeof(int)); 
    E->parallel.NODE[i]  = (struct BOUND *) malloc((nxyz+2)*sizeof(struct BOUND));
    E->parallel.IDPASS[i]  = (struct BOUND *) malloc((10)*sizeof(struct BOUND));
    E->parallel.EXCHANGE_NODE[i]  = (struct PASS *) malloc((nxyz+2)*sizeof(struct PASS));
    E->parallel.EXCHANGE_ID[i]  = (struct PASS *) malloc((nxyz*E->mesh.nsd+3)*sizeof(struct PASS));

    if (i==E->mesh.levmax)  {
      E->sien         = (struct SIEN *) malloc((nxyz+2)*sizeof(struct SIEN));
      E->surf_element = (int *) malloc((nxyz+2)*sizeof(int));
      E->surf_node    = (int *)malloc((E->lmesh.nsf+2)*sizeof(int));
      }

  }



  for(l=E->mesh.levmin;l<=E->mesh.levmax;l++) {
    for(i=1;i<=E->lmesh.NNO[l];i++)
      { E->NODE[l][i] = (INTX | INTY | INTZ);  /* and any others ... */
	E->VI[l][i] = 1.0;
	E->TW[l][i]= 0.0;
      }

    for(i=0;i<E->lmesh.NEQ[l];i++) {
      E->BI[l][i]=0.0;
      E->parallel.IDD[l][i]=0;
      }
    }
  
  for(i=1;i<=E->lmesh.nnov;i++)
    for(j=1;j<=E->mesh.nsd;j++)
       E->V[j][i] =
           E->VB[j][i] = 0.0;

  for(i=1;i<=E->lmesh.nno;i++)
    {
    E->T[i] = E->buoyancy[i] = 0.0;
    for(j=1;j<=E->mesh.nsd;j++)
      E->TB[j][i] = 0.0;}
 
  for(i=1;i<=E->lmesh.npno;i++)
    E->P[i] = 0.0;

  for(i=1;i<=E->lmesh.nel;i++)  {
      E->mat[i]=1;
      E->Qmelt[i]=0.;
      }

  for(i=0;i<E->lmesh.neq;i++)
    E->U[i] = 0.0;
 
  set_up_nonmg_aliases(E); 
  
  return;   }
  

/*  =========================================================  */


void interuption()
    
{  if (Emergency_stop++) exit(0);
   fprintf(stderr,"Cleaning up before exit\n");
   return;  }


void global_default_values(E)
     struct All_variables *E;
{
    FILE *fp;

  /* FIRST: values which are not changed routinely by the user */
 		
  E->control.v_steps_low = 10;
  E->control.v_steps_upper = 1;
  E->control.max_res_red_each_p_mg = 1.0e-3;
  E->control.accuracy = 1.0e-6;
  E->control.vaccuracy = 1.0e-8;
  E->control.true_vcycle=0;
  E->control.depth_dominated=0;
  E->control.eqn_zigzag=0;
  E->control.verbose=0; /* debugging/profiles */
  E->control.debug=0; /* also for debugging/profiles */
  E->control.stokes=0;
  E->control.restart = 0;
  E->control.testvel=0;
  E->control.testvelval[0]=0.;
  E->control.testvelval[1]=0.;
  E->control.testvelval[2]=0.;
  E->control.testvelval[3]=0.;

  /* SECOND: values for which an obvious default setting is useful */

  E->control.ORTHO = 1; /* for orthogonal meshes by default */
  E->control.ORTHOZ = 1; /* for orthogonal meshes by default */

 
    E->control.KERNEL = 0;
    E->control.CONVECTION = 0;
    E->control.SLAB = 0;
    E->control.CART2D = 0;
    E->control.CART3D = 0;
    E->control.CART2pt5D = 0;
    E->control.AXI = 0;
    E->control.CONJ_GRAD = 0;
    E->control.NMULTIGRID = 0;
    E->control.EMULTIGRID = 0;
    E->control.COMPRESS = 1;
    E->control.augmented_Lagr = 0;
    E->control.augmented = 0.0;

    /* Default: all optional modules set to `off' */
    E->control.MELTING_MODULE = 0;
    E->control.CHEMISTRY_MODULE = 0;

    E->control.GRID_TYPE=1;
    E->mesh.hwidth[1]=E->mesh.hwidth[2]=E->mesh.hwidth[3]=1.0; /* divide by this one ! */
    E->mesh.magnitude[1]=E->mesh.magnitude[2]=E->mesh.magnitude[3]=0.0;
    E->mesh.offset[1]=E->mesh.offset[2]=E->mesh.offset[3]=0.0;

  E->parallel.automa=0; E->parallel.nprocx=1; E->parallel.nprocz=1; E->parallel.nprocy=1;

  E->mesh.levmax=0;
  E->mesh.levmin=0;
  E->mesh.nox = 1;    E->mesh.nxs = 1;  E->lmesh.nox = 1;    E->lmesh.nxs = 1;
  E->mesh.noz = 1;    E->mesh.nzs = 1;  E->lmesh.noz = 1;    E->lmesh.nzs = 1;
  E->mesh.noy = 1;    E->mesh.nys = 1;  E->lmesh.noy = 1;    E->lmesh.nys = 1;

  E->viscosity.guess = 0;
  sprintf(E->viscosity.old_file,"initialize");
 
  E->control.precondition = 0;	/* for larger visc contrasts turn this back on  */
  E->control.vprecondition = 1;	
 
  /* Default bcs for T and v */
  E->mesh.toptbc = 1;   /* T prescribed at top */
  E->mesh.bottbc = 1;   /* T prescribed at bottom */
  E->mesh.xmintbc = 1;  /* T prescribed at xmin */
  E->mesh.xmaxtbc = 0;  /* heat flux prescribed at xmax */
  E->mesh.ymintbc = 0;  /* heat flux prescribed at ymin */
  E->mesh.ymaxtbc = 0;  /* heat flux prescribed at ymax */

  E->mesh.topvbc = 0;   /* stress prescribed at top*/
  E->mesh.botvbc = 0;   /* stress prescribed at bottom*/
  E->mesh.closedbot = 0;   /* default open bottom of the box */

  E->mesh.xminvxbc = 1; /* inflow prescribed at xmin in x- & y-dirs */
  E->mesh.xminvybc = 1; 
  E->mesh.xminvzbc = 0; /* free vertical flow at xmin */

  E->mesh.xmaxvxbc = 1; /* inflow prescribed at xmax in x- & y-dirs */
  E->mesh.xmaxvybc = 1; 
  E->mesh.xmaxvzbc = 0; /* free vertical flow at xmax */

  E->mesh.yminvxbc = 0; /* free-slip ymin bnd. */
  E->mesh.yminvybc = 1; 
  E->mesh.yminvzbc = 0; 

  E->mesh.ymaxvxbc = 0; /* free-slip ymax bnd. */
  E->mesh.ymaxvybc = 1; 
  E->mesh.ymaxvzbc = 0; 

  E->mesh.periodic_x=0; /* reflection is default*/
  E->mesh.periodic_y=0;

  E->control.VBXtopval=0.0;
  E->control.VBYtopval=0.0;
  E->control.VBXbotval=0.0;
  E->control.VBYbotval=0.0;

  E->control.xminvxtype=1;
  E->control.xminvytype=1;
  E->control.xminvztype=1;
  E->control.xmaxvxtype=1;
  E->control.xmaxvytype=1;
  E->control.xmaxvztype=1;

  E->control.xminvxval=0.0;
  E->control.xminvyval=0.0;
  E->control.xminvzval=0.0;
  E->control.xmaxvxval=0.0;
  E->control.xmaxvyval=0.0;
  E->control.xmaxvzval=0.0;

  /* Default model parameters */
  E->data.layer_km = 2800.0; /* Earth, whole mantle defaults */
  E->data.grav_acc = 9.81;
  E->data.therm_exp = 3.28e-5;
  E->data.Cp = 1200.0;
  E->data.therm_diff = 8.0e-7;
  E->data.therm_cond = 3.168;
  E->data.density = 3340.0;
  E->data.res_density = 3295.0;  /* density when X = ... */
  E->data.res_density_X = 0.3;
  E->data.melt_density = 2800.0;
  E->data.permeability = 3.0e-10;
  E->data.density_above = 1030.0;    /* sea water */
  E->data.gas_const = 8.3;
  E->data.surf_heat_flux = 4.4e-2;
  E->data.grav_const = 6.673e-11;
  E->data.surf_temp = 0.0;
  E->data.youngs_mod = 1.0e11;
  E->data.Te = 0.0;
  E->data.T_sol0 = 1373.0;/* Dave's values 1991 (for the earth) */
  E->data.Tsurf = 273.0;
  E->data.dTsol_dz = 3.4e-3 ;
  E->data.dTsol_dF = 440.0;
  E->data.dT_dz = 0.48e-3;
  E->data.delta_S = 250.0;
  E->data.ref_temperature = 2 * 1350.0; /* fixed temperature delta T */

  E->control.onsetssc=0;
  E->control.maxzssc=0.5;
  E->control.minxssc=0.0;
  E->control.dTprint=0;

  E->control.timdepTbc=0;
  E->control.randomT_bc=0;
  E->control.randomT_init=0;
  E->control.init_plume_cyl=0;
  E->control.fluxcontrol=0.;
  E->control.init_age_distr=0;
  E->control.plate_angle=0.;
    
  /* THIRD: fallback defaults */
 
    sprintf(E->control.data_file,"citcom.tmp.%d",getpid());
  
    E->control.NASSEMBLE = 0;
  
    E->mesh.layer[1] =  E->mesh.layer[2] =  E->mesh.layer[3] = 1.0;
    E->monitor.elapsed_time=0.0;
 
  return;  }


void global_derived_values(E)
     struct All_variables *E;

{
    int d,lx,lz,ly,i,nox,noz,noy;
    char logfile[100];
    FILE *fp;
 
    /* Set up the log file as early as possible */

    sprintf(logfile,"%s.log%d",E->control.data_file,E->parallel.me);
    
    if((fp=fopen(logfile,"w")) == NULL)
	E->fp = stdout;
    else
	E->fp = fp;



 if (E->control.NMULTIGRID||E->control.EMULTIGRID)
    { E->mesh.levmax=E->mesh.levels-1;
      E->mesh.nox = E->mesh.mgunitx * (int) pow(2.0,((double)E->mesh.levmax)) + 1;
      E->mesh.noz = E->mesh.mgunitz *(int) pow(2.0,((double)E->mesh.levmax)) + 1; 
      if(E->mesh.nsd == 3 ) 
	  E->mesh.noy = E->mesh.mgunity * (int) pow(2.0,((double)E->mesh.levmax)) + 1; 
       
    }

 if(E->mesh.nsd != 3) 
   E->mesh.noy = 1;

  E->mesh.nnx[1] = E->mesh.nox;	
  E->mesh.nnx[2] = E->mesh.noz;	
  E->mesh.nnx[3] = E->mesh.noy;	
  E->mesh.nex[1] = E->mesh.elx = E->mesh.nox-1;	
  E->mesh.nex[2] = E->mesh.elz = E->mesh.noz-1;
  E->mesh.nex[3] = E->mesh.ely = max(E->mesh.noy-1,1);
  E->mesh.nno = E->mesh.nel = 1;
  E->mesh.nmx = max(E->mesh.nox,max(E->mesh.noz,E->mesh.noy));
  for(d=1;d<=E->mesh.nsd;d++)
    { E->mesh.nno *= E->mesh.nnx[d];
      E->mesh.nel *= E->mesh.nex[d];
    }

  E->mesh.nlno = E->mesh.nno;
  if(E->mesh.periodic_x)
    { E->mesh.nlno *= E->mesh.nox-1;
      E->mesh.nlno /= E->mesh.nox; 
    }

  if(E->mesh.periodic_y && E->mesh.nsd==3)
    { E->mesh.nlno *= E->mesh.noy-1;
      E->mesh.nlno /= E->mesh.noy; 
    }

  E->mesh.nnov = E->mesh.nno;
  E->mesh.neq = E->mesh.nnov*E->mesh.nsd;

  E->mesh.npno = E->mesh.nel;
  E->mesh.nsf = E->mesh.nox*E->mesh.noy;

  for(i=E->mesh.levmax;i>=E->mesh.levmin;i--)  /* set up dimensions for different grids  */
    { if (E->control.NMULTIGRID||E->control.EMULTIGRID)
	{ nox = E->mesh.mgunitx * (int) pow(2.0,(double)i) + 1;
	  noz = E->mesh.mgunitz * (int) pow(2.0,(double)i) + 1;
	  if(E->mesh.nsd==3)
	    noy = E->mesh.mgunity * (int) pow(2.0,(double)i) + 1;
	  else 
	    noy = 1;
	}
    else 
	{ nox = E->mesh.nox;
	  noz = E->mesh.noz;
	  noy = E->mesh.noy;
	}

      E->mesh.ELX[i] = nox-1;
      E->mesh.ELZ[i] = noz-1;
      E->mesh.ELY[i] = max(noy-1,1);
      E->mesh.NNO[i] = nox * noz * noy;
      E->mesh.NLNO[i] = nox * noz * noy;
      E->mesh.NEL[i] = (nox-1) * (noz-1) * max((noy-1),1);
      E->mesh.NPNO[i] = E->mesh.NEL[i] ;
      E->mesh.NOX[i] = nox;
      E->mesh.NOZ[i] = noz;
      E->mesh.NOY[i] = noy;
      E->mesh.NMX[i] = max(nox,max(noz,noy));

      if(E->mesh.periodic_x)
	{ E->mesh.NLNO[i] *= E->mesh.NOX[i]-1;
	  E->mesh.NLNO[i] /= E->mesh.NOX[i]; 
	}
      
      if(E->mesh.periodic_y && E->mesh.nsd==3)
	{ E->mesh.NLNO[i] *= E->mesh.NOY[i]-1;
	  E->mesh.NLNO[i] /= E->mesh.NOY[i]; 
	}

      E->mesh.NNOV[i] = E->mesh.NNO[i];
      E->mesh.NEQ[i] = E->mesh.nsd * E->mesh.NNOV[i] ;  

      lx--; lz--; ly--;
    }

    if(E->control.print_convergence && E->parallel.me==0)
	fprintf(stderr,"Problem has %d x %d x %d nodes\n",E->mesh.nox,E->mesh.noz,E->mesh.noy);
  	


   return; } 


void read_initial_settings(E)
     struct All_variables *E;
  
{
    char tmp_string[100];

    void set_convection_defaults();
    void set_2dc_defaults();
    void set_3dc_defaults();
    void set_cg_defaults();
    void set_mg_defaults();
    void viscosity_parameters();

  /* first the problem type (defines subsequent behaviour) */

    input_string("Problem",E->control.PROBLEM_TYPE,NULL);
    if ( strcmp(E->control.PROBLEM_TYPE,"convection") == 0)  {
	E->control.CONVECTION = 1; 
	set_convection_defaults(E);
    }

    else if ( strcmp(E->control.PROBLEM_TYPE,"convection-chemical") == 0) {
	E->control.CONVECTION = 1;
	E->control.CHEMISTRY_MODULE=1;
	set_convection_defaults(E);
    }
    
    else {
	fprintf(E->fp,"Unable to determine problem type, assuming convection ... \n");
	E->control.CONVECTION = 1;
	set_convection_defaults(E);
    }
      
  input_string("Geometry",E->control.GEOMETRY,NULL); 
  if ( strcmp(E->control.GEOMETRY,"cart2d") == 0)
    { E->control.CART2D = 1; 
      set_2dc_defaults(E);}
  else if ( strcmp(E->control.GEOMETRY,"axi") == 0)
    { E->control.AXI = 1; 
      }
  else if ( strcmp(E->control.GEOMETRY,"cart2pt5d") == 0)
    { E->control.CART2pt5D = 1; 
      set_2pt5dc_defaults(E);}
  else if ( strcmp(E->control.GEOMETRY,"cart3d") == 0)
    { E->control.CART3D = 1;
      set_3dc_defaults(E);}
  else
    { fprintf(E->fp,"Unable to determine geometry, assuming cartesian 2d ... \n");
      E->control.CART2D = 1; 
      set_2dc_defaults(E); }

  input_string("Solver",E->control.SOLVER_TYPE,NULL);
  if ( strcmp(E->control.SOLVER_TYPE,"cgrad") == 0)
    { E->control.CONJ_GRAD = 1;
      set_cg_defaults(E);}
  else if ( strcmp(E->control.SOLVER_TYPE,"multigrid") == 0)
    { E->control.NMULTIGRID = 1;
      set_mg_defaults(E);}
  else if ( strcmp(E->control.SOLVER_TYPE,"multigrid-el") == 0)
    { E->control.EMULTIGRID = 1;
      set_mg_defaults(E);}
  else
    { if (E->parallel.me==0) fprintf(stderr,"Unable to determine how to solve, specify Solver=VALID_OPTION \n");
      exit(0); 
    }

 
  /* admin */

    /* Which files to print and which variables to calculate (default: none) */
    input_string("datatypes",E->control.which_data_files,"");
    input_string("averages",E->control.which_horiz_averages,"");
    input_string("timelog",E->control.which_running_data,"");
    input_string("observables",E->control.which_observable_data,"");

    /* Data location: server (default) or cluster node local scratch disks (use_scratch) */
    input_string("datafile",E->control.data_file,"initialize");
    input_string("use_scratch",E->control.storage,"server");
    if ( strcmp(E->control.storage,"server") == 0 )
          strcpy( E->control.data_file2, E->control.data_file ) ;
    else 
          sprintf( E->control.data_file2, 
             "/scratch_%s/%s/%s", E->parallel.machinename,
             E->control.storage, E->control.data_file ) ;

    input_string("process_command",E->control.output_written_external_command,"");
    input_boolean("AVS",&(E->control.AVS),"off");
    input_boolean("CONMAN",&(E->control.CONMAN),"off");
    
    if (E->control.NMULTIGRID||E->control.EMULTIGRID) {
	input_int("mgunitx",&(E->mesh.mgunitx),"1");
	input_int("mgunitz",&(E->mesh.mgunitz),"1");
	input_int("mgunity",&(E->mesh.mgunity),"1");
	input_int("mgunitxl",&(E->lmesh.mgunitx),"1");
	input_int("mgunitzl",&(E->lmesh.mgunitz),"1");
	input_int("mgunityl",&(E->lmesh.mgunity),"1");
	input_int("levels",&(E->mesh.levels),"0");
    }

    input_int("restart",&(E->control.restart),"0");

    input_int("regular_grid",&(E->control.ORTHOZ),"1");

    input_int("stokes_flow_only",&(E->control.stokes),"0");

    /* testvel > 0: impose an analytical test velocity instead of solving Stokes
     * (testvel=4: constant velocities in x-,z-,y-dir). */
    input_int("testvel",&(E->control.testvel),"0");
    fprintf(stderr,"testvel = %d\n", E->control.testvel);
    input_float_vector("testvelval",3,(E->control.testvelval));
    E->control.testvelval[3] = E->control.testvelval[2];
    E->control.testvelval[2] = E->control.testvelval[1];
    E->control.testvelval[1] = E->control.testvelval[0];
    fprintf(stderr,"testvelval_x = %g\n", E->control.testvelval[1]);
    fprintf(stderr,"testvelval_z = %g\n", E->control.testvelval[2]);
    fprintf(stderr,"testvelval_y = %g\n", E->control.testvelval[3]);

    input_boolean("node_assemble",&(E->control.NASSEMBLE),"off");
                                    /* general mesh structure */

    input_boolean("parallel_auto",&(E->parallel.automa),"off");
    if (E->parallel.automa==0) {
	   input_int("nprocx",&(E->parallel.nprocx),"1");
	   input_int("nprocz",&(E->parallel.nprocz),"1");
	   input_int("nprocy",&(E->parallel.nprocy),"1");
       }
       
    input_boolean("verbose",&(E->control.verbose),"off");
    input_int("debug",&(E->control.debug),"0");
    input_boolean("see_convergence",&(E->control.print_convergence),"off");
    input_boolean("COMPRESS",&(E->control.COMPRESS),"on");
    input_float("sobtol",&(E->control.sob_tolerance),"0.0001");

    input_int("obs_maxlongk",&(E->slice.maxlong),"100,1");
    input_int("obs_minlongk",&(E->slice.minlong),"1,1");

        /* for layers    */
    E->viscosity.zlm = 1.0;
    E->viscosity.zlith = 0.0;
    input_float("z_lmantle",&(E->viscosity.zlm),"1.0");
    input_float("z_lith",&(E->viscosity.zlith),"0.0");
    
    input_int ( "topvbc",     &(E->mesh.topvbc),     "0" );
    input_int ( "botvbc",     &(E->mesh.botvbc),     "0" );
    input_int ( "closedbot",  &(E->mesh.closedbot),  "0" );
    input_int ( "xminvxbc",   &(E->mesh.xminvxbc),   "1" );
    input_int ( "xminvybc",   &(E->mesh.xminvybc),   "1" );
    input_int ( "xminvzbc",   &(E->mesh.xminvzbc),   "0" );
    input_int ( "xmaxvxbc",   &(E->mesh.xmaxvxbc),   "1" );
    input_int ( "xmaxvybc",   &(E->mesh.xmaxvybc),   "1" );
    input_int ( "xmaxvzbc",   &(E->mesh.xmaxvzbc),   "0" );
    input_int ( "yminvxbc",   &(E->mesh.yminvxbc),   "0" );
    input_int ( "yminvybc",   &(E->mesh.yminvybc),   "1" );
    input_int ( "yminvzbc",   &(E->mesh.yminvzbc),   "0" );
    input_int ( "ymaxvxbc",   &(E->mesh.ymaxvxbc),   "0" );
    input_int ( "ymaxvybc",   &(E->mesh.ymaxvybc),   "1" );
    input_int ( "ymaxvzbc",   &(E->mesh.ymaxvzbc),   "0" );

    input_boolean("periodicx",&(E->mesh.periodic_x),"off");
    input_boolean("periodicy",&(E->mesh.periodic_y),"off");
    input_boolean("depthdominated",&(E->control.depth_dominated),"off");
    input_boolean("eqnzigzag",&(E->control.eqn_zigzag),"off");
    input_boolean("eqnviscosity",&(E->control.eqn_viscosity),"off");

    input_float("topvbxval",&(E->control.VBXtopval),"0.0");
    input_float("botvbxval",&(E->control.VBXbotval),"0.0");
    input_float("topvbyval",&(E->control.VBYtopval),"0.0");
    input_float("botvbyval",&(E->control.VBYbotval),"0.0");

    input_int("xminvxtype",&(E->control.xminvxtype),"1");
    input_int("xminvytype",&(E->control.xminvytype),"1");
    input_int("xminvztype",&(E->control.xminvztype),"1");
    input_int("xmaxvxtype",&(E->control.xmaxvxtype),"1");
    input_int("xmaxvytype",&(E->control.xmaxvytype),"1");
    input_int("xmaxvztype",&(E->control.xmaxvztype),"1");

    input_float("xminvxval",&(E->control.xminvxval),"0.0");
    input_float("xminvyval",&(E->control.xminvyval),"0.0");
    input_float("xminvzval",&(E->control.xminvzval),"0.0");
    input_float("xmaxvxval",&(E->control.xmaxvxval),"0.0");
    input_float("xmaxvyval",&(E->control.xmaxvyval),"0.0");
    input_float("xmaxvzval",&(E->control.xmaxvzval),"0.0");
  
    /* Plate b.c.: constant age or age variation at t=0 */
    input_int("init_age_distr",&(E->control.init_age_distr),"0");

    /* Plume b.c.: optional initial hot cylinder in the lowermost layer */
    input_int("init_plume_cyl",&(E->control.init_plume_cyl),"0");

    input_int("toptbc",&(E->mesh.toptbc),"1");
    input_int("bottbc",&(E->mesh.bottbc),"1");
    input_float("toptbcval",&(E->control.TBCtopval),"0.0");
    input_float("bottbcval",&(E->control.TBCbotval),"1.0");
    input_int("xmintbc",&(E->mesh.xmintbc),"1");
    input_int("xmaxtbc",&(E->mesh.xmaxtbc),"0");
    input_int("ymintbc",&(E->mesh.ymintbc),"1");
    input_int("ymaxtbc",&(E->mesh.ymaxtbc),"1");

    /* Optional random T noise (amplitude randomT); requires time-dep bc */
    input_int("timdepTbc",&(E->control.timdepTbc),"0");
    input_float("randomT",&(E->control.randomT_bc),"0.0");
    /* In case randomT is given, T-bc are always time-dep: */
    if (E->control.randomT_bc != 0.0) 
    {  E->control.timdepTbc = 1;
       srand(E->parallel.me+1);
    }
    input_float("randomT_init",&(E->control.randomT_init),"0.0");

    /* Optional plume radius update to obtain the desired plume flux */
    input_int("fluxcontrol",&(E->control.fluxcontrol),"0");
    input_float("flux",&(E->control.flux),"0.");
    input_float("relaxfluxupdate",&(E->control.relaxfluxupdate),"1.");
    input_int("updatefluxevery",&(E->control.updatefluxevery),"100");

    /* Possibility for onset time calc. of small scale conv */
    input_int("onsetssc",&(E->control.onsetssc),"0");
    input_float("maxzssc",&(E->control.maxzssc),"0.5");
    input_float("minxssc",&(E->control.minxssc),"0.0");
    input_int("dTprint",&(E->control.dTprint),"0");

    E->control.plate_vel=0.;
    input_float("plate_velocity",&(E->control.plate_vel),"0.0");
    E->control.old_plate_vel = E->control.plate_vel;
    input_float("old_plate_velocity",&(E->control.old_plate_vel),
                    "E->control.plate_vel");
                    /* Plate velocity in a ssc-restart with rotated
                       plate motion */
    if (E->parallel.me == 0)
    {  fprintf(stderr, "plate-vels.:\n");
       fprintf(stderr, "   old=%g\n", E->control.old_plate_vel);
       fprintf(stderr, "   new=%g\n", E->control.plate_vel);
    }
    input_float("plate_age",&(E->control.plate_age),"0.0");
    input_float("plate_angle",&(E->control.plate_angle),"0.0");
                   /* Input in degs, and NOT in rads !!! */
    input_float("plume_radius",&(E->segment.plume_radius),"0.0");
    input_float("plume_DT",&(E->segment.plume_DT),"0.0");
    input_float("plume_x",&(E->segment.plume_coord[1]),"0.0");
    input_float("plume_z",&(E->segment.plume_coord[2]),"0.0");
    input_float("plume_y",&(E->segment.plume_coord[3]),"0.0");
 
    input_float("blyr_hwx1",&(E->mesh.bl1width[1]),"nodefault");
    input_float("blyr_hwz1",&(E->mesh.bl1width[2]),"nodefault");
    input_float("blyr_hwy1",&(E->mesh.bl1width[3]),"nodefault");
    input_float("blyr_hwx2",&(E->mesh.bl2width[1]),"nodefault");
    input_float("blyr_hwz2",&(E->mesh.bl2width[2]),"nodefault");
    input_float("blyr_hwy2",&(E->mesh.bl2width[3]),"nodefault");
    input_float("blyr_mgx1",&(E->mesh.bl1mag[1]),"nodefault");
    input_float("blyr_mgz1",&(E->mesh.bl1mag[2]),"nodefault");
    input_float("blyr_mgy1",&(E->mesh.bl1mag[3]),"nodefault");
    input_float("blyr_mgx2",&(E->mesh.bl2mag[1]),"nodefault");
    input_float("blyr_mgz2",&(E->mesh.bl2mag[2]),"nodefault");
    input_float("blyr_mgy2",&(E->mesh.bl2mag[3]),"nodefault");
   
   
    input_float("region_wdx",&(E->mesh.width[1]),"nodefault");
    input_float("region_wdz",&(E->mesh.width[2]),"nodefault");
    input_float("region_wdy",&(E->mesh.width[3]),"nodefault");
    input_float("region_hwx",&(E->mesh.hwidth[1]),"nodefault");
    input_float("region_hwz",&(E->mesh.hwidth[2]),"nodefault");
    input_float("region_hwy",&(E->mesh.hwidth[3]),"nodefault");
    input_float("region_mgx",&(E->mesh.magnitude[1]),"nodefault");
    input_float("region_mgz",&(E->mesh.magnitude[2]),"nodefault");
    input_float("region_mgy",&(E->mesh.magnitude[3]),"nodefault");
    input_float("region_ofx",&(E->mesh.offset[1]),"nodefault");
    input_float("region_ofz",&(E->mesh.offset[2]),"nodefault");
    input_float("region_ofy",&(E->mesh.offset[3]),"nodefault");

    input_string("gridxfile",E->mesh.gridfile[1]," ");
    input_string("gridzfile",E->mesh.gridfile[2]," ");
    input_string("gridyfile",E->mesh.gridfile[3]," ");
    
    
    input_float("dimenx",&(E->mesh.layer[1]),"nodefault");
    input_float("dimenz",&(E->mesh.layer[2]),"nodefault");
    input_float("dimeny",&(E->mesh.layer[3]),"nodefault");
    
    input_int("nodex",&(E->mesh.nox),"nodefault,1,nomax");
    input_int("nodez",&(E->mesh.noz),"nodefault,1,nomax");
    input_int("nodey",&(E->mesh.noy),"1,1,nomax");
    input_boolean("aug_lagr",&(E->control.augmented_Lagr),"off");
    input_double("aug_number",&(E->control.augmented),"0.0");

    input_float("tole_compressibility",&(E->control.tole_comp),"0.0");
    input_boolean("orthogonal",&(E->control.ORTHO),"on");
    input_boolean("crust",&(E->control.crust),"off");
    input_float("crust_width",&(E->crust.width),"0.0");

    input_int("storage_spacing",&(E->control.record_every),"10");
    input_int("storage_always_before",&(E->control.record_all_until),"5");
 
    input_boolean("precond",&(E->control.precondition),"off");
    input_boolean("vprecond",&(E->control.vprecondition),"on");
    input_int("mg_cycle",&(E->control.mg_cycle),"2,0,nomax");
    input_int("down_heavy",&(E->control.down_heavy),"1,0,nomax");
    input_int("up_heavy",&(E->control.up_heavy),"1,0,nomax");
    input_double("accuracy",&(E->control.accuracy),"1.0e-4,0.0,1.0");
    input_int("viterations",&(E->control.max_vel_iterations),"250,0,nomax");

 
    input_int("vhighstep",&(E->control.v_steps_high),"1,0,nomax");
    input_int("vlowstep",&(E->control.v_steps_low),"250,0,nomax");
    input_int("vupperstep",&(E->control.v_steps_upper),"1,0,nomax");
    input_int("piterations",&(E->control.p_iterations),"100,0,nomax");
    input_int("maxsamevisc",&(E->control.max_same_visc),"25,0,nomax");

  /* data section */ 

  input_float("ReferenceT",&(E->data.ref_temperature),"2600.0");
  input_float("Q0",&(E->control.Q0),"0.0");
  input_float("layerd",&(E->data.layer_km),"2800.0");
  input_float("gravacc",&(E->data.grav_acc),"9.81");
  input_float("thermexp",&(E->data.therm_exp),"3.28e-5");
  input_float("cp",&(E->data.Cp),"1200.0");
  input_float("H2O",&(E->data.H2O),"0.02");
  input_float("thermdiff",&(E->data.therm_diff),"8.0e-7");
  input_float("thermcond",&(E->data.therm_cond),"3.168");
  input_float("density",&(E->data.density),"3340.0");
  input_float("mdensity",&(E->data.melt_density),"2800.0");
  input_float("wdensity",&(E->data.density_above),"1030.0");
  input_float("rdensity",&(E->data.res_density),"3295.0");
  /* Adiabatic gradient [K/km] = alpha*g*T_ref/Cp; used by Convection.c and BMO core
   * coupling.  Not read from input — derived from thermexp, gravacc, cp, ReferenceT. */
  E->data.adiabatic_gr = E->data.therm_exp * E->data.grav_acc
                       * E->data.ref_temperature / E->data.Cp * 1000.0f;
  if (E->parallel.me == 0)
      fprintf(stderr,
          "[adiabatic_gr] auto = thermexp*g*ReferenceT/cp = %.5g K/km"
          "  (thermexp=%.4g  g=%.4g  T_ref=%.4g  cp=%.4g)\n",
          (double)E->data.adiabatic_gr,
          (double)E->data.therm_exp, (double)E->data.grav_acc,
          (double)E->data.ref_temperature, (double)E->data.Cp);
  if (E->data.adiabatic_gr > 0.)
    if (E->control.TBCbotval==1.0)
      E->control.TBCbotval+=E->data.layer_km*E->data.adiabatic_gr/E->data.ref_temperature;
  input_int("shear_heating",&(E->control.shear_heating),"0");
  input_float("heat_of_melt",&(E->data.heat_of_melt),"400000.");    /*[J/kg]*/
  input_float("heatflux",&(E->data.surf_heat_flux),"4.4e-2");
  input_float("refvisc",&(E->data.ref_viscosity),"nodefault");
  input_float("meltvisc",&(E->data.melt_viscosity),"nodefault");
  input_float("surftemp",&(E->data.surf_temp),"0.0");
  input_float("youngs",&(E->data.youngs_mod),"1.0e11");
  input_float("Te",&(E->data.Te),"0.0");
  input_float("Tsol0",&(E->data.T_sol0),"1373.0");
  input_float("dTsoldz",&(E->data.dTsol_dz),"3.4e-3");
  input_float("dTsoldF",&(E->data.dTsol_dF),"440.0");
  input_float("dTdz",&(E->data.dT_dz),"0.48e-3");
  input_float("deltaS",&(E->data.delta_S),"250.0");
  input_float("gasconst",&(E->data.gas_const),"8.3");     /* not much cause to change these ! */
  input_float("gravconst",&(E->data.grav_const),"6.673e-11");
  input_float("permeability",&(E->data.permeability),"3.0e-10");
  input_float("MO_coolingtime",&(E->control.MO_coolingtime),"1.");
  input_float("MO_decay_tau",&(E->control.MO_decay_tau),"0.33");
  input_float("MO_depth",&(E->control.MO_depth_initial),"2000.");
  input_float("protoM_depth",&(E->control.protoM_depth),"1450."); /* km */
  E->control.protoM_depth/=E->data.layer_km;
  input_float("protoM_density",&(E->control.protoM_density),"2000."); /* kg/m3 */
  input_float("MO_temperature",&(E->control.MO_temperature),"2500.");
  E->control.MO_temperature/=E->data.ref_temperature;
  input_float("MO_temp_final",&(E->control.MO_temp_final),"1500.");
  E->control.MO_temp_final/=E->data.ref_temperature;

  /* ----- Labrosse-2025 implicit BMO coupling (init_age_distr=160) ----- */
  /* BMO_* parameters are dimensionless, except BMO_TL0 (>1: in K) and
     BMO_Qcmb / BMO_HBMO (total power in TW if BMO_heat_in_TW=1). */
  input_int  ("BMO_implicit_mode",&(E->control.BMO_implicit_mode),"0");
  input_int  ("BMO_heat_in_TW"   ,&(E->control.BMO_heat_in_TW)   ,"1");
  /* BMO_xi_l0 / BMO_Cw_l0 / BMO_Chpe_l0: bulk compositions; at t=0 (phi0 = 1 - MO_depth/layer_km)
   * melt X_L = X_bulk/[phi0+(1-phi0)*D], solid X_S = D*X_L; Rayleigh fractionation thereafter. */
  input_float("BMO_xi_l0"        ,&(E->control.BMO_xi_l0)        ,"0.10");
  input_float("BMO_partition_D"  ,&(E->control.BMO_partition_D)  ,"0.5"); /* K_Fe exchange coeff */
  input_float("BMO_TL_dPdim"     ,&(E->control.BMO_TL_dPdim)     ,"0.30");
  input_float("BMO_TL_dxi"       ,&(E->control.BMO_TL_dxi)       ,"-0.50");

  /* Fe# -> FeO mass-fraction mapping (liquidus only; BMO_TL_dxi is per unit w_FeO):
   *   BMO_MgSi_mode = 0 : constant Mg/Si = BMO_Mg_Si (default 1.25, Du et al. 2019).
   *   BMO_MgSi_mode = 1 : Mg/Si(Fe#) table from BMO_MgSi_table_file (columns: Fe#, Mg/Si). */
  input_int  ("BMO_MgSi_mode"    ,&(E->control.BMO_MgSi_mode)    ,"0");
  input_float("BMO_Mg_Si"        ,&(E->control.BMO_Mg_Si)        ,"1.25");
  input_string("BMO_MgSi_table_file", E->control.BMO_MgSi_table_file, "");
  E->control.BMO_MgSi_ntab      = 0;
  E->control.BMO_MgSi_tab_FeNum = NULL;
  E->control.BMO_MgSi_tab_val   = NULL;
  E->control.BMO_w_melt0        = -1.0f;   /* set later from initial Fe#+Mg/Si */
  if (E->control.BMO_MgSi_mode == 1 && E->control.BMO_MgSi_table_file[0] != '\0') {
      FILE *ftab = fopen(E->control.BMO_MgSi_table_file, "r");
      if (ftab == NULL) {
          if (E->parallel.me == 0)
              fprintf(stderr,
                  "[BMO] WARNING: cannot open BMO_MgSi_table_file '%s'; "
                  "falling back to constant Mg/Si=%.4g.\n",
                  E->control.BMO_MgSi_table_file, (double)E->control.BMO_Mg_Si);
          E->control.BMO_MgSi_mode = 0;
      } else {
          int   cap = 64, ntab = 0;
          char  line[256];
          float *fe = (float *) malloc(cap * sizeof(float));
          float *mg = (float *) malloc(cap * sizeof(float));
          while (fgets(line, sizeof(line), ftab) != NULL) {
              float a, b;
              char *p = line;
              while (*p == ' ' || *p == '\t') p++;
              if (*p == '#' || *p == '\n' || *p == '\r' || *p == '\0') continue;
              if (sscanf(line, "%f %f", &a, &b) == 2) {
                  if (ntab >= cap) {
                      cap *= 2;
                      fe = (float *) realloc(fe, cap * sizeof(float));
                      mg = (float *) realloc(mg, cap * sizeof(float));
                  }
                  fe[ntab] = a;
                  mg[ntab] = b;
                  ntab++;
              }
          }
          fclose(ftab);
          if (ntab < 2) {
              if (E->parallel.me == 0)
                  fprintf(stderr,
                      "[BMO] WARNING: BMO_MgSi_table_file '%s' has <2 rows; "
                      "falling back to constant Mg/Si=%.4g.\n",
                      E->control.BMO_MgSi_table_file,
                      (double)E->control.BMO_Mg_Si);
              free(fe); free(mg);
              E->control.BMO_MgSi_mode = 0;
          } else {
              E->control.BMO_MgSi_ntab      = ntab;
              E->control.BMO_MgSi_tab_FeNum = fe;
              E->control.BMO_MgSi_tab_val   = mg;
              if (E->parallel.me == 0)
                  fprintf(stderr,
                      "[BMO] Loaded Mg/Si(Fe#) table '%s' (%d rows): "
                      "Fe#=[%.4g,%.4g], Mg/Si=[%.4g,%.4g]\n",
                      E->control.BMO_MgSi_table_file, ntab,
                      (double)fe[0], (double)fe[ntab-1],
                      (double)mg[0], (double)mg[ntab-1]);
          }
      }
  }

  input_float("BMO_z0_init"      ,&(E->control.BMO_z0_init)      ,"-1.0");
  input_float("BMO_TL0"          ,&(E->control.BMO_TL0)          ,"-1.0");
  input_float("BMO_St"           ,&(E->control.BMO_St)           ,"1.0");
  input_float("BMO_Qcmb"         ,&(E->control.BMO_Qcmb)         ,"0.0");
  input_float("BMO_HBMO"         ,&(E->control.BMO_HBMO)         ,"0.0");

  /* Dynamic Stefan number: if BMO_dS_Cp > 0, St = T_BOMB * BMO_dS_Cp each step (default 0: fixed St). */
  input_float("BMO_dS_Cp"        ,&(E->control.BMO_dS_Cp)        ,"0.0");
  
  /* .BMO_comp.dat sampling: when MO_depth_row changes (BMO_comp_per_row=1), every dz of
   * front motion (BMO_comp_dz>0), or every record_every step (both off). */
  input_int  ("BMO_comp_per_row" ,&(E->control.BMO_comp_per_row) ,"1.0");
  input_float("BMO_comp_dz"      ,&(E->control.BMO_comp_dz)      ,"0.0");
  E->control.BMO_comp_per_row = 1.0;

  /* H2O partitioning in BMO (same bulk -> melt -> solid framework as HPE):
   * BMO_Cw_l0 = bulk H2O mass fraction, BMO_Dw = partition coefficient,
   * BMO_TL_dCw = liquidus slope dT_L/dC_w (<=0; default 0 disables water-liquidus coupling). */
  input_float("BMO_Cw_l0"        ,&(E->control.BMO_Cw_l0)        ,"0.0");
  input_float("BMO_Dw"           ,&(E->control.BMO_Dw)           ,"1.0e-3");
  input_float("BMO_TL_dCw"       ,&(E->control.BMO_TL_dCw)       ,"0.0");

  /* HPE (heat-producing elements): BMO_Chpe_l0 = bulk mass fraction, BMO_Dhpe = partition
   * coefficient, BMO_H0_hpe = reference heating (Q_R = H0_hpe*(Chpe_l/Chpe_l_melt0)*s, s = 1-z),
   * BMO_tau_hpe = dimensionless half-life. */
  input_float("BMO_Chpe_l0"      ,&(E->control.BMO_Chpe_l0)      ,"0.0");
  input_float("BMO_Dhpe"         ,&(E->control.BMO_Dhpe)         ,"0.01");
  input_float("BMO_tau_hpe"      ,&(E->control.BMO_tau_hpe)      ,"0.02");
  input_float("BMO_H0_hpe"       ,&(E->control.BMO_H0_hpe)       ,"0.0");

  /* Solid-mantle HPE heating (BMO_hpe_mantle_on=1, decay-free tracer convention):
   * tracers hold fractionation-only Chpe; decay exp(-t*ln2/tau_hpe) is applied globally in
   * Hbmo and in the mantle source Q_HPE = H_mantle*mult*cHPE(el)*decay.  0 = legacy behaviour. */
  E->control.BMO_hpe_mantle_on = 0;
  E->control.BMO_H_mantle      = 0.0f;
  E->control.BMO_H_mantle_mult = 1.0f;
  input_int  ("BMO_hpe_mantle_on" ,&(E->control.BMO_hpe_mantle_on) ,"0");
  input_float("BMO_H_mantle"      ,&(E->control.BMO_H_mantle)      ,"0.0");
  input_float("BMO_H_mantle_mult" ,&(E->control.BMO_H_mantle_mult) ,"1.0");
  if (E->control.BMO_hpe_mantle_on && E->parallel.me == 0)
      fprintf(stderr,
          "[BMO_hpe_mantle] mantle HPE heating ON: H_mantle=%.4g,"
          " mult=%.4g, tau_hpe=%.4g (decay-free tracer mode)\n",
          (double)E->control.BMO_H_mantle,
          (double)E->control.BMO_H_mantle_mult,
          (double)E->control.BMO_tau_hpe);

   
  /*   BMO_xi_sat   : FeO legacy cap (unused by K_Fe/D_eff path; kept for input compat).
   *   BMO_Cw_sat   : H2O; ~0.10 (10 wt%) recommended; range 0.05-0.20.
   *   BMO_Chpe_sat : HPE; ~0.10 numerical/mass-conservation safeguard. */
  input_float("BMO_xi_sat"       ,&(E->control.BMO_xi_sat)       ,"1.0");
  input_float("BMO_Cw_sat"       ,&(E->control.BMO_Cw_sat)       ,"0.1");
  input_float("BMO_Chpe_sat"     ,&(E->control.BMO_Chpe_sat)     ,"0.2");

  E->control.BMO_xi_sat = 1.0;
  E->control.BMO_Cw_sat = 0.1;
  E->control.BMO_Chpe_sat = 0.2;

  /* If BMO_TL0 is given in K (>1), assume Kelvin and rescale. */
  if (E->control.BMO_TL0 > 1.0)
      E->control.BMO_TL0 /= E->data.ref_temperature;
  /* Sentinel -1 means: take initial values from the legacy MO_* parameters. */
  if (E->control.BMO_TL0 < 0.0)
      E->control.BMO_TL0 = E->control.MO_temperature;
  /* BMO_z0_init must be dimensionless (z in [0,1]); MO_depth_initial is still in km here
   * (converted later in setup_plume_problem), so divide by layer_km. */
  if (E->control.BMO_z0_init < 0.0)
      E->control.BMO_z0_init = E->control.MO_depth_initial / E->data.layer_km;
  else if (E->control.BMO_z0_init > 1.0f)  /* user supplied a km value */
      E->control.BMO_z0_init /= E->data.layer_km;

  /* Depth-layered liquidus pressure slope mP: [BMO_z0_init, 1] is split into
   * BMO_TL_dPdim_nlayers bands (0 = scalar BMO_TL_dPdim everywhere). */
  input_int("BMO_TL_dPdim_nlayers", &(E->control.BMO_TL_dPdim_nlayers), "0");
  if (E->control.BMO_TL_dPdim_nlayers > 8) {
      if (E->parallel.me == 0)
          fprintf(stderr,
              "[BMO] WARNING: BMO_TL_dPdim_nlayers=%d exceeds max (8); capped.\n",
              E->control.BMO_TL_dPdim_nlayers);
      E->control.BMO_TL_dPdim_nlayers = 8;
  }
  if (E->control.BMO_TL_dPdim_nlayers > 0) {
      int il_;
      /* Pre-fill with the scalar fallback so unspecified layers inherit it. */
      for (il_ = 0; il_ < E->control.BMO_TL_dPdim_nlayers; il_++)
          E->control.BMO_TL_dPdim_vals[il_] = E->control.BMO_TL_dPdim;
      input_float_vector("BMO_TL_dPdim_layers",
                         E->control.BMO_TL_dPdim_nlayers,
                         E->control.BMO_TL_dPdim_vals);
      if (E->parallel.me == 0) {
          float dz_ = (1.0f - E->control.BMO_z0_init)
                      / (float)E->control.BMO_TL_dPdim_nlayers;
          fprintf(stderr,
              "[BMO] Depth-layered mP: %d layers in z=[%.4g,1.0], dz=%.4g\n",
              E->control.BMO_TL_dPdim_nlayers,
              (double)E->control.BMO_z0_init, (double)dz_);
          for (il_ = 0; il_ < E->control.BMO_TL_dPdim_nlayers; il_++)
              fprintf(stderr,
                  "[BMO]   layer[%d] z=[%.4g,%.4g]  mP=%.5g\n",
                  il_,
                  (double)(E->control.BMO_z0_init + il_      * dz_),
                  (double)(E->control.BMO_z0_init + (il_+1)  * dz_),
                  (double)E->control.BMO_TL_dPdim_vals[il_]);
      }
  }

  /* Optional TW input for CMB heat flow and BMO radiogenic heating. */
  if (E->control.BMO_heat_in_TW) {
      float BMO_TW_to_Qdimless(), BMO_TW_to_Hdimless();
      float qcmb_tw = E->control.BMO_Qcmb;
      float hbmo_tw = E->control.BMO_HBMO;
      float s0 = 1.0f - E->control.BMO_z0_init;

      E->control.BMO_HBMO_TW = hbmo_tw;
      E->control.BMO_Qcmb = BMO_TW_to_Qdimless(E, qcmb_tw);
      E->control.BMO_HBMO = BMO_TW_to_Hdimless(E, hbmo_tw, s0);
      if (E->parallel.me == 0)
          fprintf(stderr,
              "[BMO_heat] TW input: Q_CMB=%.4g TW -> %.4g (dimless),"
              " H_BMO=%.4g TW -> %.4g (dimless at s0=%.4g)\n",
              (double)qcmb_tw, (double)E->control.BMO_Qcmb,
              (double)hbmo_tw, (double)E->control.BMO_HBMO, (double)s0);
  } else
      E->control.BMO_HBMO_TW = 0.0f;

  E->control.BMO_xi_l_melt0   = 0.0f;
  E->control.BMO_Cw_l_melt0   = 0.0f;
  E->control.BMO_Chpe_l_melt0 = 0.0f;
  E->control.BMO_xi_s0        = 0.0f;
  E->control.BMO_Cw_s0        = 0.0f;
  E->control.BMO_Chpe_s0      = 0.0f;
  {
      void BMO_init_compositions_from_bulk();
      BMO_init_compositions_from_bulk(E);
  }
  E->control.BMO_zdot      = 0.0f;
  E->control.BMO_QbomT     = 0.0f;
  E->control.BMO_XFe_prof  = NULL;
  E->control.BMO_CH2O_prof = NULL;
  E->control.BMO_Chpe_prof = NULL;
  E->control.BMO_prof_size = 0;
  E->control.BMO_rho_nv_prof = NULL;
  E->control.BMO_rho_nv_anchor = 0.0f;
  E->control.BMO_prem_rho_anchor = 0.0f;
  E->control.BMO_rho_nv_prem = 0;
  /* BMO_TL_dPdim_nlayers / _vals are set by the input_* calls above. */

  /* FeO + H2O combined density mode (Du et al. 2019); solid/melt split at MO_depth_initial
   * (z=0 surface, z=1 CMB); BMO_density_solid_above=1 swaps the test. */
  input_int  ("BMO_density_mode"        ,&(E->control.BMO_density_mode)       ,"0");
  input_int  ("BMO_density_solid_above" ,&(E->control.BMO_density_solid_above),"0");
  /* BMO_rho_melt: density of the SOLID zone reference [kg/m^3]         */
  input_float("BMO_rho_melt"            ,&(E->control.BMO_rho_melt)            ,"4000.0");
  /* BMO_rho_solid_ref: anhydrous non-volatile base for Du+2019 MELT rho*/
  input_float("BMO_rho_solid_ref"       ,&(E->control.BMO_rho_solid_ref)       ,"-1.0");
  /* BMO_drho_dXFe kept for backward-compatible input files; not used in formula */
  input_float("BMO_drho_dXFe"           ,&(E->control.BMO_drho_dXFe)           ,"0.0");
  /* BMO_XFe_ref: Fe# fallback used when tracfunc 80 is inactive (0–1)           */
  input_float("BMO_XFe_ref"             ,&(E->control.BMO_XFe_ref)             ,"0.0");
  /* BMO_n_spin: low-spin fraction n in literature formula (0.28 + 0.04*n)       *
   * Du et al. (2019) / Karki et al. (2018) recommend n = 0.8                   */
  input_float("BMO_n_spin"              ,&(E->control.BMO_n_spin)              ,"0.8");
  /* BMO_alpha_water: Du+2019 rho/rho0-1 = alpha*X_H2O[wt%]; tracers store mass fraction */
  input_float("BMO_alpha_water"         ,&(E->control.BMO_alpha_water)         ,"-0.97");
  /* BMO_rho_nv_prem: 1 = rho_nv(z) = BMO_rho_solid_ref * PREM(z)/PREM(MO) */
  input_int  ("BMO_rho_nv_prem"         ,&(E->control.BMO_rho_nv_prem)         ,"0");
  /* BMO_T_override_mode: if 1, reset T=MO_temperature in melt zone     */
  input_int  ("BMO_T_override_mode"     ,&(E->control.BMO_T_override_mode)     ,"0");
  /* Negative sentinel for BMO_rho_solid_ref: use CITCOM reference rho  */
  if (E->control.BMO_rho_solid_ref < 0.0f)
      E->control.BMO_rho_solid_ref = E->data.density;
  if (E->parallel.me == 0 && E->control.BMO_density_mode == 1) {
      float coeff = 0.28f + 0.04f * E->control.BMO_n_spin;
      fprintf(stderr,
          "[BMO_density] mode=1 (literature formula): z_init=%.4f"
          "  rho_0=%.1f  n_spin=%.2f  coeff(0.28+0.04n)=%.4f"
          "  XFe_ref(fallback)=%.4f  alpha_w=%.4f"
          "  rho_nv_prem=%d  solid_above=%d  T_override=%d\n",
          (double)E->control.MO_depth_initial,
          (double)E->control.BMO_rho_solid_ref,
          (double)E->control.BMO_n_spin,
          (double)coeff,
          (double)E->control.BMO_XFe_ref,
          (double)E->control.BMO_alpha_water,
          E->control.BMO_rho_nv_prem,
          E->control.BMO_density_solid_above,
          E->control.BMO_T_override_mode);
  }

  /* Core-BMO thermal coupling (Labrosse 2025, Eqs. 37-40): BMO_core_coupling=1 enables it;
   * BMO_core_beta (<=0: auto-computed from the core parameters rhoN, CpN, Lp, Ap, gamma, rCMB). */
  input_int  ("BMO_core_coupling", &(E->control.BMO_core_coupling), "0");
  input_float("BMO_core_beta",     &(E->control.BMO_core_beta),    "-1.0");
  input_float("BMO_core_rhoN",     &(E->control.BMO_core_rhoN),    "12500.0");
  input_float("BMO_core_CpN",      &(E->control.BMO_core_CpN),     "750.0");
  input_float("BMO_core_Lp_km",    &(E->control.BMO_core_Lp_km),   "7200.0");
  input_float("BMO_core_Ap",       &(E->control.BMO_core_Ap),      "0.484");
  input_float("BMO_core_gamma",    &(E->control.BMO_core_gamma),   "1.5");
  input_float("BMO_core_rcmb_km",  &(E->control.BMO_core_rcmb_km), "3480.0");
  E->control.BMO_T_CMB = 0.0f;

  if (E->control.BMO_core_coupling == 1) {
      float BMO_core_compute_beta();
      /* Auto-compute beta from physical core parameters if not set by user */
      if (E->control.BMO_core_beta <= 0.0f)
          E->control.BMO_core_beta = BMO_core_compute_beta(E);
      if (E->parallel.me == 0)
          fprintf(stderr,
              "[BMO_core] Labrosse-2025 core coupling ON: beta=%.4f"
              "  rhoN=%.4g kg/m3  CpN=%.4g J/kgK  Lp=%.4g km"
              "  Ap=%.4f  gamma=%.4f  rCMB=%.4g km\n",
              (double)E->control.BMO_core_beta,
              (double)E->control.BMO_core_rhoN,
              (double)E->control.BMO_core_CpN,
              (double)E->control.BMO_core_Lp_km,
              (double)E->control.BMO_core_Ap,
              (double)E->control.BMO_core_gamma,
              (double)E->control.BMO_core_rcmb_km);
  }

 (E->problem_settings)(E);

 viscosity_parameters(E);

return; }

void check_bc_consistency(E)
     struct All_variables *E;

{ int i,lev;

  for(i=1;i<=E->lmesh.nno;i++)
    { if ((E->node[i] & VBX) && (E->node[i] & SBX))
	printf("Inconsistent x velocity bc at %d\n",i);
      if ((E->node[i] & VBZ) && (E->node[i] & SBZ))
	printf("Inconsistent z velocity bc at %d\n",i);
      if ((E->node[i] & VBY) && (E->node[i] & SBY))
	printf("Inconsistent y velocity bc at %d\n",i);
      if ((E->node[i] & TBX) && (E->node[i] & FBX))
	printf("Inconsistent x temperature bc at %d\n",i);
      if ((E->node[i] & TBZ) && (E->node[i] & FBZ))
	printf("Inconsistent z temperature bc at %d\n",i);
      if ((E->node[i] & TBY) && (E->node[i] & FBY))
	printf("Inconsistent y temperature bc at %d\n",i); }

  for(lev=E->mesh.levmin;lev<=E->mesh.levmax;lev++)
    for(i=1;i<=E->lmesh.NNO[lev];i++)
      { if ((E->NODE[lev][i] & VBX) && (E->NODE[lev][i]  & SBX))
	  printf("Inconsistent x velocity bc at %d,%d\n",lev,i);
	if ((E->NODE[lev][i]  & VBZ) && (E->NODE[lev][i]  & SBZ))
	  printf("Inconsistent z velocity bc at %d,%d\n",lev,i);
	if ((E->NODE[lev][i]  & VBY) && (E->NODE[lev][i]  & SBY))
	  printf("Inconsistent y velocity bc at %d,%d\n",lev,i);
	/* Tbc's not applicable below top level */ }

  return;

}

void set_up_nonmg_aliases(E)
     struct All_variables *E;
     
{ /* Aliases for functions only interested in the highest mg level */
  int i;

  E->eco = E->ECO[E->mesh.levmax]; 
  E->ien = E->IEN[E->mesh.levmax];
  E->id = E->ID[E->mesh.levmax];
  E->lm = E->LMD[E->mesh.levmax];
  E->Vi = E->VI[E->mesh.levmax];
  E->EVi = E->EVI[E->mesh.levmax];
  E->node = E->NODE[E->mesh.levmax];
  E->tw = E->TW[E->mesh.levmax];
  E->Mass = E->MASS[E->mesh.levmax];
  E->gDA = E->GDA[E->mesh.levmax];
  E->gNX = E->GNX[E->mesh.levmax];

  for (i=1;i<=E->mesh.nsd;i++)
    E->X[i] = E->XX[E->mesh.levmax][i];
 
  return; }

report(E,string)
     struct All_variables *E;
     char * string;
{ if(E->control.verbose && E->parallel.me==0)
    { fprintf(stderr,"%s\n",string);
      fflush(stderr);
    }
  return;
}

record(E,string)
     struct All_variables *E;
     char * string;
{ if(E->control.verbose)
    { fprintf(E->fp,"%s\n",string);
      fflush(E->fp);
    }

  return;
}



/* ========================================== */

void initial_pressure(E)
     struct All_variables *E;
{
    int i,node,ii;

    for(i=1;i<=E->lmesh.npno;i++)
	    E->P[i]=0.0;

  return; 
}

void initial_velocity(E)
     struct All_variables *E;
{
    int i,node,ii;

    for(i=1;i<=E->lmesh.nnov;i++)   {
	E->V[1][i]=0.0;
	E->V[2][i]=0.0;
	}

    if(E->mesh.dof==3)   {
 	for(i=1;i<=E->lmesh.nnov;i++)
	    E->V[3][i]=0.0;
        }

    for(i=1;i<=E->lmesh.nnov;i++) {
       E->Vest[1][i] = E->V[1][i];
       E->Vest[2][i] = E->V[2][i];
       if(E->mesh.dof==3)
          E->Vest[3][i] = E->V[3][i];
	  }

    return; 
}
