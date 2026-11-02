 	/* This file contains the definitions of variables which are passed as arguments */
	/* to functions across the whole filespace of CITCOM. #include this file everywhere !*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#if defined(__osf__) 
void *Malloc1();
#endif

#define Malloc0(a) Malloc1((a),__FILE__,__LINE__)


/* #define Malloc0 malloc */

#define LIDN 0x1
#define VBX 0x2
#define VBZ 0x4
#define VBY 0x8
#define TBX 0x10
#define TBZ 0x20
#define TBY 0x40
#define TZEDGE 0x80
#define TXEDGE 0x100
#define TYEDGE 0x200
#define VXEDGE 0x400
#define VZEDGE 0x800
#define VYEDGE 0x1000
#define INTX 0x2000
#define INTZ 0x4000
#define INTY 0x8000
#define SBX 0x10000
#define SBZ 0x20000
#define SBY 0x40000
#define FBX 0x80000
#define FBZ 0x100000
#define FBY 0x200000

#define OFFSIDE 0x400000

#define LIDE 1

#ifndef COMPRESS_BINARY
#define COMPRESS_BINARY "/usr/bin/compress"
#endif

#define MAX_LEVELS 12
#define MAX_F    10
#define MAX_S    30

/* Macros */

#define max(A,B) (((A) > (B)) ? (A) : (B))
#define min(A,B) (((A) < (B)) ? (A) : (B))
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

typedef float higher_precision;  /* matrix coeffs etc */
typedef double higher_precision1; /* intermediate calculations for finding above coeffs */


/* Common structures */

struct Rect {
    int numb;
    char overlay[40];
    float x1[40];
    float x2[40];
    float z1[40];
    float z2[40];
    float y1[40];
    float y2[40];
    float halfw[40];
    float mag[40];
} ;
 

struct Circ {
    int numb;
    char overlay[40];
    float x[40];
    float z[40];
    float y[40];
    float rad[40];
    float mag[40];
    float halfw[40];
};


struct Harm {
    int numb;
    int harms;
    char overlay[40];
    float off[40];
    float x1[40];
    float x2[40];
    float z1[40];
    float z2[40];
    float y1[40];
    float y2[40];
    float kx[20][40];
    float kz[20][40];
    float ky[20][40];
    float ka[20][40];
    float phx[20][40];
    float phz[20][40];
    float phy[20][40];
};

struct Erfc {
 int numb;


};

struct RectBc {
    int numb;
    char norm[40];
    float intercept[40];
    float x1[40];
    float x2[40];
    float z1[40];
    float z2[40];
    float halfw[40];
    float mag[40];
} ;
 

struct CircBc {
    int numb;
    char norm[40];
    float intercept[40];
    float x[40];
    float z[40];
    float rad[40];
    float mag[40];
    float halfw[40];
};


struct PolyBc {
    int numb;
    int order;
    char norm[40];
    float intercept[40];
    float x1[40];
    float x2[40];
    float z1[40];
    float z2[40];
    float ax[20][40];
    float az[20][40];
};


struct HarmBc {
    int numb;
    int harms;
    char norm[40];
    float off[40];
    float intercept[40];
    float x1[40];
    float x2[40];
    float z1[40];
    float z2[40];
    float kx[20][40];
    float kz[20][40];
    float ka[20][40];
    float phx[20][40];
    float phz[20][40];
 };


struct Shape_function_dA  {
  float vpt[8];
  float ppt[1]; };

struct Shape_function1_dA  {
  double vpt[6*4];
  double ppt[6*1]; };

struct Shape_function1 	{ 
    double vpt[4*4];  /* node & gauss pt */
    double ppt[4*1];  };

struct Shape_function 	{ 
    double vpt[8*8];  /* node & gauss pt */
    double spt[8*4];  /* node & gauss pt */
    double ppt[8*1];  };

struct Shape_function_dx 	{ 
    float vpt[3*8*8]; /* dirn & node & gauss pt */
    float ppt[3*8*1];  };

struct Shape_function1_dx 	{ 
    double vpt[2*4*4]; /* dirn & node & gauss pt */
    double ppt[2*4*1];  };

struct EG { 
    higher_precision g[24][1]; };

struct EK2 { 
    double k[8*8]; };

struct EK { 
    double k[24*24]; };

struct MEK { 
    double nint[9]; };
 
struct NK {
    higher_precision *k;
    int *map;
};

struct COORD {
    float centre[4];
    float size[4];
    float recip_size[4];
    float area;   } ;

struct SUBEL { 
    int sub[9];   };
			
struct ID  { 
    int doff[6];	}; /* can  be 1 or 2 or 3 */
struct IEN {
    int node[9];	};
struct FNODE {
    float node[9];	};
struct SIEN {
    int node[5];	};
struct LM  { 
    struct { int doff[4]; } node[9]; } ;

struct NEI {
    int *nels;
    int *lnode;
    int *element; };

struct Slabs  { float   max_depth    ;
                float   min_depth    ;
		float   *dip  [MAX_F];  
		float   *strike [MAX_F];  
                float   *width[MAX_F];
		float   *depth[MAX_F];
		int     total_node;
		int     *vert;
		float   age_max;
		float   *age;
		int     *node;};

struct Crust  { float   width;
                float   thickness;
                float   *T;
		int     total_node;
		int     *node;};

struct Segment  { 

                int zlayers;
                int nzlayer[40];
                float zzlayer[40];
                int xlayers;
                int nxlayer[40];
                float xxlayer[40];
                int ylayers;
                int nylayer[40];
                float yylayer[40];

		float   plume_coord[4];
		float   plume_radius;
		float   plume_DT;
		float   *Tz_in;   /* For y-independent inflow */
		float   *Tzy_in;   /* For y-dependent inflow */
		float   *Tz_out;   /* For y-independent outflow */
		float   *Vx;   /* For y-independent inflow */
		float   *Vx_zy;   /* For y-dependent inflow */
		int     sym_nx;
		};

struct BOUND  { 
    int bound[8];	}; 

struct PASS  { 
    int pass[27];	}; 

struct Parallel {
		int me;
		int nproc;
		int nprocx;
		int nprocz;
		int nprocy;
		int nprocxy;
		int nproczy;
		int nprocxz;
		int automa;
		int idb;
		int me_loc[4];
		int num_b;
        int *IDD[MAX_LEVELS];
        int *ELE_ORDER[MAX_LEVELS];
        int *NODE_ORDER[MAX_LEVELS];
        struct BOUND *NODE[MAX_LEVELS];
        struct BOUND NUM_NNO[MAX_LEVELS];
        struct BOUND NUM_PASS[MAX_LEVELS];
		struct BOUND NUM_ELE[MAX_LEVELS];
		struct BOUND *IDPASS[MAX_LEVELS];

        int TNUM_PASS[MAX_LEVELS];
        int START_ELE[MAX_LEVELS];
        int START_NODE[MAX_LEVELS];
		struct PASS NUM_NEQ[MAX_LEVELS];
		struct PASS NUM_NODE[MAX_LEVELS];
		struct PASS PROCESSOR[MAX_LEVELS];
		struct PASS *EXCHANGE_ID[MAX_LEVELS];
        struct PASS *EXCHANGE_NODE[MAX_LEVELS];
        char machinename[160];

		};

struct MESH_DATA {/* general information concerning the fe mesh */ 
    int nsd;        /* Spatial extent 1,2,3d*/
    int dof;        /* degrees of freedom per node */
    int levmax;     
    int levmin;
    int levels;
    int mgunitx;
    int mgunitz;
    int mgunity;
    int NEQ[MAX_LEVELS];	/* All other values refer to the biggest mesh (& lid)  */
    int NNO[MAX_LEVELS];
    int NNOV[MAX_LEVELS];
    int NLNO[MAX_LEVELS];
    int NPNO[MAX_LEVELS];
    int NEL[MAX_LEVELS];
    int NOX[MAX_LEVELS];
    int NOZ[MAX_LEVELS];
    int NOY[MAX_LEVELS];
    int NMX[MAX_LEVELS];
    int NXS[MAX_LEVELS];
    int NZS[MAX_LEVELS];
    int NYS[MAX_LEVELS];
    int ELX[MAX_LEVELS];
    int ELZ[MAX_LEVELS];
    int ELY[MAX_LEVELS];
    int LNDS[MAX_LEVELS];
    int LELS[MAX_LEVELS];
    int neqd;
    int neq;
    int nno;
    int nnov;
    int nlno;
    int npno;
    int nel;
    int snel;
    int nex[4]; /* general form of ... */
    int elx;
    int elz;
    int ely;
    int nnx[4]; /* general form of ... */
    int nox;
    int noz;
    int noy;
    int nxs;
    int nzs;
    int nys;
    int nmx;
    int nsf; /* nodes for surface observables */
    int toptbc;
    int bottbc;
    int xmintbc;
    int xmaxtbc;
    int ymintbc;
    int ymaxtbc;
    int topvbc;
    int botvbc;
    int closedbot;
    int xminvxbc;
    int xminvybc;
    int xminvzbc;
    int xmaxvxbc;
    int xmaxvybc;
    int xmaxvzbc;
    int yminvxbc;
    int yminvybc;
    int yminvzbc;
    int ymaxvxbc;
    int ymaxvybc;
    int ymaxvzbc;
    int xminvbc;
    int xmaxvbc;
    int yminvbc;
    int ymaxvbc;

    char topvbc_file[100];
    char botvbc_file[100];
    char sidevbc_file[100];
    char gridfile[4][100];


    int periodic_x;
    int periodic_y;
    float layer[4];			/* dimensionless dimensions */
    float lidz;
    float bl1width[4],bl2width[4],bl1mag[4],bl2mag[4];
    float hwidth[4],magnitude[4],offset[4],width[4]; /* grid compression information */ 
    int fnodal_malloc_size;
    int dnodal_malloc_size;
    int feqn_malloc_size;
    int deqn_malloc_size;
    int bandwidth;
    int null_source;
    int null_sink;
    int matrix_size[MAX_LEVELS];
} ;

struct HAVE {    /* horizontal averages */
    float *T;
    float *Vi;
    float *Rho;
    float *f;
    float *F;
    float *vrms;
    float *V[4];
    float *lT_y; /* local average T in y-dir. for onsettime of SSC */
    float *gT_y; /* global average T in y-dir. for onsettime of SSC */
    float *arrh;
    float *Vy;
    float *edot31;
  
};

struct SLICE {    /* horizontally sliced data, including topography */
    float *tpg;
    float *tpgb;
    float *grv;
    float *geo;
    float *geok;
    float *grvk;
    float *grvb;
    float *geob;
    float *geobk;
    float *grvbk;
    float *tpgk;
    float *tpgbk;
    float *shflux;
    float *bhflux;
    float *cen_mflux;
    float *vxsurf[3];    /* surface velocity vectors */
    float *vline;        /* for kernels, velocity at force term */
    float *vlinek;
    float *tpglong;
    float *tpgblong;
    float *grvlong;
    float *geolong;
    float *grvblong;
    float *geoblong;
    float Nub;
    float Nut;

    int minlong;
    int maxlong;
  };

struct BAVE {
    float T;
    float Vi;
    double V[4]; };


struct TOTAL {
    float melt_prod;  };

struct MONITOR {
    char node_output[100][6];  /* recording the format of the output data */
    char sobs_output[100][6];  /* recording the format of the output data */
    int node_output_cols;
    int sobs_output_cols;

    int solution_cycles;
    int sdepv_cycles;

    float  time_scale;
    float  length_scale;
    float  viscosity_scale;
    float  geoscale;
    float  tpgscale;
    float  grvscale;
  
    float  delta_v_last_soln;
    float  elapsed_time;
    float  elapsed_time_vsoln;
    float  elapsed_time_vsoln1;
    float  reference_stress;
    float  incompressibility;
    float  vdotv;
    float  nond_av_heat_fl;  
    float  nond_av_adv_hfl;  
    float  cpu_time_elapsed;
    float  cpu_time_on_vp_it;
    float  cpu_time_on_forces;
    float  cpu_time_on_mg_maps;
    float  tpgkmag;
    float  grvkmag;
   
    float  Nusselt;
    float  Vmax;
    float  Vsrms;
    float  Vrms, vVrms;
    float  Vrms_surface;
    float  Vrms_base;
    float  F_surface;
    float  F_base;
    float  Frat_surface;
    float  Frat_base;
    float  T_interior;
    float  Sigma_max;
    float  Sigma_interior;
    float  Vi_average;
   
};

struct CONTROL {
    int PID;

    char output_written_external_command[500];   /* a unix command to run when output files have been created */

    int ORTHO,ORTHOZ;   /* indicates levels of mesh symmetry */
    char B_is_good[MAX_LEVELS];  /* general information controlling program flow */
    char Ahat_is_good[MAX_LEVELS];  /* general information controlling program flow */
    char old_P_file[100];
    char data_file[100];
    char data_file2[100];
    char storage[100];

    char which_data_files[1000];
    char which_horiz_averages[1000];
    char which_running_data[1000];
    char which_observable_data[1000];
  
    char PROBLEM_TYPE[20]; /* one of ... */
    int KERNEL;
    int CONVECTION;
    int SLAB;	
    char GEOMETRY[20]; /* one of ... */
    int CART2D;
    int CART2pt5D;
    int CART3D;
    int AXI;	 
    char SOLVER_TYPE[20]; /* one of ... */
    int DIRECT;
    int restart;
    int restart_frame;
    int CONJ_GRAD;
    int NMULTIGRID;
    int EMULTIGRID;
    int DIRECTII;
    char NODE_SPACING[20]; /* turns into ... */
    int GRID_TYPE;
    int COMPRESS;
    int AVS;
    int CONMAN;
    int stokes;
    int eedot_available;
    int shear_heating;
    int first_tstep;
    int testvel;
    float testvelval[4];
    int testweirdmesh;
 
    int timdepTbc; /* 0=notimdepTbc,1=timdepTbc*/
    float randomT_bc;
    float randomT_init;

    int fluxcontrol; /* =1: adjusting plume radius for desired flux */
    float flux; 
    float relaxfluxupdate; 
    int updatefluxevery; 
    float inittreduct;

    int onsetssc; /* Calc. of onsettime of small scale conv */
    float maxzssc; /* Max. (dimless) depth to which onset is 
                      searched for */
    float minxssc; /* Min. (dimless) x-coord from which ssc is 
                      searched for */
    int dTprint; /* print y-av. z-max. T-T_coolinghalfspace(x) */

    int dfact;
    double penalty;
    int augmented_Lagr;
    double augmented;
    int faults;
    int NASSEMBLE;
    int comparison;
    int crust;
    float onset_ssc;
    float plate_vel, old_plate_vel, plate_age;
    float plate_angle;
    float timestepmax;
    float CFL_number;    /* Courant number for advective timestep limiter (default 0.5) */
    float max_dT_frac;   /* max fractional temperature change per step (default 1.0 = disabled) */
    int   n_early_steps; /* number of initial steps with linearly-ramped timestep reduction (default 50) */
    int init_age_distr;
    int init_plume_cyl;

    float tole_comp;
  
    float sob_tolerance;
 
    float Atemp; 
    float VBXtopval;
    float VBXbotval;
    float VBYtopval;
    float VBYbotval;

    float TBCtopval;
    float TBCbotval;
    int MO_depth_row;
    float MO_depth;
    float MO_depth_initial;
    float MO_temperature;
    float MO_temp_final;
    float MO_coolingtime;
    float MO_decay_tau;     /* exponential decay: tau as fraction of MO_coolingtime (default 0.33) */
    float protoM_depth,protoM_density;

    int   BMO_implicit_mode;   /* 0 = legacy prescribed front; 1 = self-consistent */
    float BMO_xi_l;            /* current liquid-side Fe# xi_l (state)               */
    float BMO_xi_l0;           /* bulk Fe# X_Fe (input); molar Fe/(Fe+Mg)             */
    float BMO_xi_l_melt0;      /* melt Fe# at t=0 (liquidus anchor; from bulk+phi0)  */
    float BMO_xi_s0;           /* solid Fe# at t=0: D_eff(K_Fe,xi_m)*xi_m (phi0 eq.) */
    float BMO_partition_D;     /* Fe-Mg exchange coefficient K_Fe (<1); BMO_partition_D name kept */
    float BMO_TL_dPdim;        /* dimensionless m_P = (dT_L/dP)*rho*g*d/dT (scalar fallback) */
    int   BMO_TL_dPdim_nlayers;/* number of depth layers for mP (0 = use scalar above)       */
    float BMO_TL_dPdim_vals[8];/* mP per layer [0..nlayers-1], shallowest-first               */
    float BMO_TL_dxi;          /* dimensionless m_xi = (dT_L/dxi)/dT (typically <0)   */
    float BMO_z0_init;         /* z_0 used as anchor in liquidus expression           */
    float BMO_TL0;             /* T_L(z_0) in dimensionless form (=initial T_BOMB)    */
    float BMO_St;              /* Stefan number L/(Cp*dT)                             */
    int   BMO_heat_in_TW;      /* 1: BMO_Qcmb/BMO_HBMO input is total power [TW]      */
    float BMO_Qcmb;            /* dimensionless heat flux Q_CMB into BMO from below   */
    float BMO_HBMO;            /* dimensionless radiogenic heating in BMO             */
    float BMO_HBMO_TW;         /* total BMO radiogenic power [TW] if BMO_heat_in_TW=1 */
    float BMO_zdot;            /* most recent solved dz_BOMB/dt (diagnostic+CFL)      */
    float BMO_QbomT;           /* most recent thermal heat flux Q_BOMB^T (diagnostic) */

    float BMO_dS_Cp;           /* Ds/Cp ratio for dynamic St; 0 = use fixed BMO_St    */

    float BMO_Cw_l;            /* current liquid-side H2O mass fraction (state var)   */
    float BMO_Cw_l0;           /* bulk H2O mass fraction X_H2O (input)                */
    float BMO_Cw_l_melt0;      /* melt H2O at t=0 (liquidus anchor; from bulk+phi0)  */
    float BMO_Cw_s0;           /* solid H2O at t=0: Dw*Cw_melt0 (phi0 equilibrium)   */
    float BMO_Dw;              /* H2O solid/liquid partition coeff (D_H2O ~ 1e-3)     */
    float BMO_TL_dCw;          /* dT_L/dC_w liquidus slope for water (dim'less, <=0)  */

    float BMO_Chpe_l;          /* current liquid-side HPE mass fraction (state var)   */
    float BMO_Chpe_l0;         /* bulk HPE mass fraction X_HPE (input)              */
    float BMO_Chpe_l_melt0;    /* melt HPE at t=0 (Q_R anchor; from bulk+phi0)        */
    float BMO_Chpe_s0;         /* solid HPE at t=0: Dhpe*Chpe_melt0 (phi0 eq.)        */
    float BMO_Dhpe;            /* HPE solid/liquid partition coeff (K_HPE)            */
    float BMO_tau_hpe;         /* dimensionless HPE radioactive half-life             */
    float BMO_H0_hpe;          /* ref specific heating rate at Chpe_l=Chpe_l0          */

    /* HPE heating of the solid mantle after crystallisation from the BMO.
     * When BMO_hpe_mantle_on=1 the code switches to a "decay-free" tracer
     * convention (Lebec/Labrosse 2023): Chpe stored in tracers/profiles is the
     * fractionation-only concentration (no accumulated decay), and decay is
     * applied as a global factor exp(-t*ln2/tau) in both Hbmo (BMO energy
     * balance) and the mantle source term Q_HPE below.
     *   BMO_H_mantle       : dimensionless specific heating coefficient H such that
     *                        Q_HPE(el,t) = H * mult * cHPE(el) * exp(-t*ln2/tau).
     *   BMO_H_mantle_mult  : dimensionless sensitivity-test multiplier (default 1). */
    int   BMO_hpe_mantle_on;   /* 0=legacy (off); 1=mantle HPE heating active          */
    float BMO_H_mantle;        /* dimensionless specific heating coeff for mantle HPE  */
    float BMO_H_mantle_mult;   /* sensitivity-test multiplier on BMO_H_mantle          */

    /* Melt saturation caps: once the liquid concentration reaches these,    */
    /* the dilute-solution (constant-D) assumption breaks down, the species  */
    /* saturates, and crystallisation switches to equilibrium (D_eff=1):     */
    /* the liquid concentration plateaus (removing the Rayleigh s->0          */
    /* singularity) and the excess crystallises into the solid (mass con*).  */
    float BMO_xi_sat;          /* FeO legacy cap; unused by K_Fe/D_eff FeO path       */
    float BMO_Cw_sat;          /* H2O melt saturation mass fraction                   */
    float BMO_Chpe_sat;        /* HPE melt cap (numerical/mass-conservation safeguard)*/

    /* Ternary Fe# -> FeO mass-fraction mapping for the liquidus (L2 scheme).
     * The Labrosse (2024) liquidus slope BMO_TL_dxi is defined per unit FeO
     * MASS fraction w_FeO in [0,1]; the state variable xi_l is molar Fe#
     * (= Fe/(Fe+Mg)).  We map Fe# -> w_FeO through a FeO-MgO-SiO2 ternary
     * composition of prescribed oxide-molar Mg/Si.  This affects ONLY the
     * liquidus temperature term T_L(w); density/buoyancy remain Fe#-only
     * (Du 2019) and Rayleigh fractionation still evolves Fe#.               */
    int   BMO_MgSi_mode;        /* 0=constant Mg/Si; 1=table Mg/Si(Fe#) (Route 1)      */
    float BMO_Mg_Si;            /* constant oxide-molar Mg/Si (mode 0; table fallback) */
    float BMO_w_melt0;          /* w_FeO at t=0 melt Fe# (fixed liquidus reference w0) */
    char  BMO_MgSi_table_file[100]; /* Fe#  Mg/Si  two-column table (mode 1)          */
    int   BMO_MgSi_ntab;        /* number of rows loaded from the table                */
    float *BMO_MgSi_tab_FeNum;  /* ascending Fe# knots                                 */
    float *BMO_MgSi_tab_val;    /* Mg/Si at each knot                                  */

    float *BMO_XFe_prof;       /* FeO solid profile indexed by grid row [1..NOZ]      */
    float *BMO_CH2O_prof;      /* H2O solid profile indexed by grid row [1..NOZ]      */
    float *BMO_Chpe_prof;      /* HPE solid profile indexed by grid row [1..NOZ]      */
    int   BMO_prof_size;       /* size of profile arrays (0 = not yet allocated)      */
    float *BMO_rho_nv_prof;    /* rho_nv(z) on grid rows [1..NOZ] [kg/m^3]            */
    float BMO_rho_nv_anchor;   /* = BMO_rho_solid_ref at MO_depth [kg/m^3]            */
    float BMO_prem_rho_anchor; /* PREM rho at MO_depth anchor [kg/m^3]                */
    int   BMO_rho_nv_prem;     /* 0 = constant rho_nv; 1 = PREM-scaled rho_nv(z)      */

    int   BMO_density_mode;       /* 0 = legacy comporayl; 1 = literature formula     */
    int   BMO_density_solid_above;/* 0: solid at z<MO_depth_init; 1: solid at z>init  */
    float BMO_rho_melt;           /* kept for backward compat (unused in mode=1)       */
    float BMO_rho_solid_ref;      /* rho_0: anhydrous non-vol. base density [kg/m^3]  */
    float BMO_drho_dXFe;          /* DEPRECATED: no longer used in density formula     */
    float BMO_XFe_ref;            /* Fe# fallback when tracfunc 80 inactive (0–1)      */
    float BMO_n_spin;             /* low-spin fraction n in (0.28+0.04n); default 0.8  */
    float BMO_alpha_water;        /* Du+2019 water coeff alpha (per wt% H2O in rho formula) */
    int   BMO_T_override_mode;    /* 0 = off; 1 = reset T=MO_T in melt zone each step  */

    int   BMO_core_coupling;    /* 0=off; 1=Labrosse-2025 core coupling enabled          */
    float BMO_core_beta;        /* dimensionless core heat capacity coefficient           */
    float BMO_core_rhoN;        /* core central density ρ_N [kg/m³]  (default 12500)     */
    float BMO_core_CpN;         /* core specific heat C_pN [J/kg/K]  (default 750)       */
    float BMO_core_Lp_km;       /* core pressure scale height L_p [km] (default 7200)    */
    float BMO_core_Ap;          /* core structure parameter A_p       (default 0.484)     */
    float BMO_core_gamma;       /* Grüneisen parameter γ              (default 1.5)       */
    float BMO_core_rcmb_km;     /* CMB radius [km]                    (default 3480)      */
    float BMO_T_CMB;            /* current T_CMB [dimensionless] (diagnostic, Eq. 39)    */
    float BMO_comp_dz;          /* min |dMO_depth| between .BMO_comp.dat lines; 0=off   */
    int   BMO_comp_per_row;     /* 1: write when MO_depth_row changes (grid-aligned)    */

    int xminvxtype;
    int xminvytype;
    int xminvztype;
    int xmaxvxtype;
    int xmaxvytype;
    int xmaxvztype;

    float xminvxval;
    float xminvyval;
    float xminvzval;
    float xmaxvxval;
    float xmaxvyval;
    float xmaxvzval;

    float Q0;
    
    int precondition;
    int vprecondition;
    int keep_going;
    int v_steps_low;
    int v_steps_high;
    int v_steps_upper;
    int max_vel_iterations;
    int p_iterations;
    int max_same_visc;
    float max_res_red_each_p_mg;
    float sub_stepping_factor;
    int mg_cycle;
    int true_vcycle;
    int down_heavy;
    int up_heavy;
    int depth_dominated;
    int eqn_viscosity;
    int eqn_zigzag;
    int verbose;
    int debug;
    double accuracy;
    double vaccuracy; 
   
    int total_iteration_cycles;
    int total_v_solver_calls;
    
    int record_every;
    int record_all_until;

    int print_convergence;
    int sdepv_print_convergence;

    int show_time;
    double initial_time;

     /* modules */
    int MELTING_MODULE;
    int CHEMISTRY_MODULE;
};

struct DATA {  
    float  layer_km;
    float   grav_acc;
    float   therm_exp;
    float   Cp;
    float  H2O;
    float  DioverRa;
    float  therm_diff;
    float  therm_cond;
    float   density;
    float  res_density;
    float  res_density_X;
    float   melt_density;
    float   density_above;
    float  adiabatic_gr;
    float  deltarhodepl;
    float  deltarhomelt;
    float  heat_of_melt;
    float   gas_const;
    float   surf_heat_flux;
    float  ref_viscosity;
    float   melt_viscosity;
    float   permeability;
    float   grav_const;
    float   surf_temp;
    float   youngs_mod; 
    float   Te;
    float   ref_temperature;
    float   Tsurf;
    float   T_sol0;
    float   buoy_flux;
    float   delta_S;
    float   dTsol_dz;
    float   dTsol_dF;
    float   dT_dz; };
	      
#include "tracers.h"
struct All_variables {     
#include "viscosity_descriptions.h"
#include "Convection_variables.h"
#include "temperature_descriptions.h"
#include "advection.h"
    struct TRACERS tracers;

    FILE *fp;
    FILE *filed[20];

    struct HAVE Have;
    struct BAVE Bulkave;
    struct TOTAL Total;
    struct MESH_DATA mesh;
    struct MESH_DATA lmesh;
    struct CONTROL control;
    struct MONITOR monitor;
    struct DATA data;
    struct SLICE slice;
    struct COORD *eco;
    struct IEN *ien;  /* global */
    struct SIEN *sien;
    struct ID *id;
    struct COORD *ECO[MAX_LEVELS];
    struct IEN *IEN[MAX_LEVELS]; /* global at each level */
    struct ID *ID[MAX_LEVELS];
    struct NEI NEI[MAX_LEVELS];
    struct SUBEL *EL[MAX_LEVELS];
    struct EG *elt_del[MAX_LEVELS];
    struct EK *elt_k[MAX_LEVELS];
    struct FNODE *TWW[MAX_LEVELS];
    struct Segment segment;
    struct Slabs slabs;
    struct Crust crust;

    struct Parallel parallel;
   
    higher_precision *Eqn_k1[MAX_LEVELS];  
    higher_precision *Eqn_k2[MAX_LEVELS];  
    higher_precision *Eqn_k3[MAX_LEVELS];  
    int *Node_map[MAX_LEVELS];
    int *Node_eqn[MAX_LEVELS];
    int *Node_k_id[MAX_LEVELS];

    double *BI[MAX_LEVELS];      /* inv of  diagonal elements of K matrix */
    double *BPI[MAX_LEVELS];
    float *V[4];                 /* velocity X[dirn][node] can save memory */
    float *V1[4];                /* velocity X[dirn][node] can save memory */
    float *Vest[4];              /* velocity X[dirn][node] can save memory */
    float *Vstar[4];             /* velocity X[dirn][node] can save memory */

    double *P,*F,*H,*S,*U;
    float *stress;
    float *Psi;
    float *NP;
    float *edot;               /* strain rate invariant */
    float *Mass;               /* lumped mass matrix (diagonal elements) for p-g solver etc. */  
    float *tw;
    float *X[4];
    float *XL[4];    /* JvH: 3 1-D arrays for local x-,y-,and 
                             z-distribution seperately */
    float *XX[MAX_LEVELS][4],*MASS[MAX_LEVELS];

    float *Vi,*EVi;
    float *T,*buoyancy;
    float *Tstar;
    float *Qmelt;                /* latent heat of melt */
    float *eedot;
    float *visc_diss;
    float *Tdot;		
    float *VI[MAX_LEVELS];	/* viscosity has to soak down to all levels */
    float *EVI[MAX_LEVELS];	/* element viscosity has to soak down to all levels */
    float *VB[4],*TB[4];	/* boundary conditions for V,T defined everywhere */
    float *TW[MAX_LEVELS];	/* nodal weightings */

    int *surf_node, *surf_element;
    int *mat;	        /* properties of mat */
    unsigned int *node;	        /* properties of node */
    unsigned int *NODE[MAX_LEVELS];
    unsigned int *eqn;
    unsigned int *EQN[MAX_LEVELS];
 		
    double *temp;
    double **global_K;  /* direct solver stuff */
    double **factor_K;
    double *global_F;
    struct LM *lmd;
    struct LM *lm;
    struct LM *LMD[MAX_LEVELS];

    struct Shape_function_dx *GNX[MAX_LEVELS];
    struct Shape_function_dA *GDA[MAX_LEVELS];
    struct Shape_function_dx *gNX;
    struct Shape_function_dA *gDA;

  
    struct Shape_function1 M; /* master-element shape funtions */
    struct Shape_function1_dx Mx; 
    struct Shape_function N;
    struct Shape_function_dx Nx;
    struct Shape_function1 L; /* master-element shape funtions */
    struct Shape_function1_dx Lx; 
 
    void (* build_forcing_term)();
    void (* iterative_solver)();
    void (* next_buoyancy_field)();
    void (* obtain_gravity)();
    void (* problem_settings)();
    void (* problem_derived_values)();
    void (* problem_allocate_vars)();
    void (* problem_boundary_conds)();
    void (* problem_node_positions)();
    void (* problem_update_node_positions)();
    void (* problem_initial_fields)();
    void (* problem_update_bcs)();
    void (* special_process_new_velocity)();
    void (* special_process_new_buoyancy)();
    void (* solve_stokes_problem)(); 
    void (* solver_allocate_vars)(); 
    void (* transform)();

    float (* node_space_function[3])();
 
};


			    







