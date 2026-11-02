/* in this file define the contents of the TRACER data structure */

#define NTRACFUNCMAX 10

static const int trac_sr_code [27] =
{222, /*  0: default: not used */
 122, /*  1: xmin  ------------- 2-D or 3-D ------------- */
 322, /*  2: xmax */
 212, /*  3: zmin */
 232, /*  4: zmax */
 112, /*  5: xmin+zmin */
 312, /*  6: xmax+zmin */
 132, /*  7: xmin+zmax */
 332, /*  8: xmax+zmax */
 221, /*  9: ymin  ------------- 3-D only ------------- */
 223, /* 10: ymax */
 121, /* 11: xmin+ymin */
 321, /* 12: xmax+ymin */
 123, /* 13: xmin+ymax */
 323, /* 14: xmax+ymax */
 211, /* 15: zmin+ymin */
 231, /* 16: zmax+ymin */
 213, /* 17: zmin+ymax */
 233, /* 18: zmax+ymax */
 111, /* 19: xmin+zmin+ymin */
 131, /* 20: xmin+zmax+ymin */
 311, /* 21: xmax+zmin+ymin */
 331, /* 22: xmax+zmax+ymin */
 113, /* 23: xmin+zmin+ymax */
 133, /* 24: xmin+zmax+ymax */
 313, /* 25: xmax+zmin+ymax */
 333  /* 26: xmax+zmax+ymax */
};
/* process-data exchange order: checkerboard shape.
   Order in which neighbors are contacted to exchange tracer info:
   use parallel_order1 if iprocx+iprocy+iprocz = even 
   use parallel_order2 if iprocx+iprocy+iprocz = odd  
   (first proc has iprocx=iprocz=iprocy=0)
   In this way processors are not waiting too much.
 */
static const int parallel_order1[27] = 
{ 0,
  1,2,3,4,5,6,7,8,  /* 2-D and 3-D */
  9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 /* 3-D only */
};
static const int parallel_order2[27] = 
{ 0,
  2,1,4,3,8,7,6,5,  /* 2-D and 3-D */
  10,9,14,13,12,11,18,17,16,15,26,25,24,23,22,21,20,19 /* 3-D only */
};

struct Sendrec   /* parallel: help array to pass info to neighbors */
{   int used;          /* =0 if unused, =1 if used */
    int size;          /* length of the arrays in this struct */
    int procnr;        /* processor to send/recv info to/from */
    int n1;            /* counter / # filled positions in arrays */
    int n2;            /* counter / # filled positions in arrays */
    int *tracnr;       /* storage of tracer nrs. */
    int *elemnr;       /* storage of tracer nrs. */
    float *vec1[4];    /* storage1 of tracer crds/vels (to/from nghb) */
    float *vec2[4];    /* storage2 of tracer crds/vels (to/from nghb) */
    float *func1[NTRACFUNCMAX]; /* storage1 of trac funcs (t/f nghb) */
    float *func2[NTRACFUNCMAX]; /* storage2 of trac funcs (t/f nghb) */
};

struct TRACERS {
    int tracers_used; /* flag: 0=no, 1=yes */
    int initialized; /* flag: 0=no, 1=yes */
    int distr;
    int ntrac;
    int lntrac_usual;
    int lntrac_bottom_usual;
    int lntrac;
    int lntrac1;
    int lntrac_stor;
    int larraylen;
    float array_buf_factor; /* buffer fraction added to lntrac for array allocation (default 0.5) */
    float tracer_density;
    float **xtrac;
    float *xtrac1[4];
    float **xtrac2;
    float *vtrac[4];
    float *vtrac1[4];
    float **ftrac;
    float **ftrac2;
    float *felem[NTRACFUNCMAX];
    int *elemnr;
    int *elemnr1;
    int *ntrac_in_elem;
    int *ntrac_in_elem1;
    int *tracs_in_elem;
    int *tracs_in_elem1;
    int *elemtraccounter;
    int *crossborder; /* parallel: indicates if/which crossing border */
    int ntracfunc;
    int ntracbuoyfunc;
    int *tracfunc;
    int *interp_t2e;
    int *interp_e2t;
    int tractest;
    int Warn;
    float lxmin, lymin, lzmin;
    float lxmax, lymax, lzmax;
    float *comporayl;
    float comporheo;
    int compoviscoused;
    float *compobuoy;
    float mindiffts;
    struct Sendrec sr[27]; /* see above for struct definition */
    int order[27];
    int type;
    float dzz; /* used for tracer.type=2: vertical spacing of tracers*/
    int restart;
    int storage;
    float *trapmelt;
    float melting_at_trac;
    float melting_per_step;
    int extraction;
    float phiC;
    float phiR;
    float Fcorr[102][125];
    int melt_columns;
    int melt_rows;
    int melt_lcolumns;
    int melt_lrows;
    int melt_firstcolumn;
    int melt_lastcolumn;
    float melt_xdistance;
    float melt_ydistance;
    float melttracer_density,thermochemical_layerd;
    float trapmelt_per_step,trapmelt_total;
    int melttracers_used;
    int bot_inflow;
    int melttracer_flag;
    float **melt;
    double *melt_xbound;
    double *melt_ybound;
    float *init_at_inflow0;
    float *init_at_inflow1;
    int Q0compo_tracnr;
    float Q0compo_value;
    float lambda;
    int Q0compo_U_Th_tracnr;
    int hpe_mantle_tracnr;     /* ifunc index of tracfunc 83 (HPE solid profile);
                                * -1 if tracfunc 83 is absent. Set during init. */
} tracers;
