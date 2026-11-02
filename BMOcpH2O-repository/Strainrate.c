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

void strain_rate_2_inv(E,EEDOT,SQRT)
     struct All_variables *E;
     float *EEDOT;
     int SQRT;
{
    double tst;
    double dudx[4][4];
    double edot[4][4];
    float VV[4][9], temp, temp1, age, simvar, Ths, z;
    float velo_scale;
    float Thsold, zold, TT,ZZ,visc,visc1,visc2,stress; 
    int l, iz, ielz;
    char filename[100];
    FILE *fp;
    static int visits=0;
    static float Ataupown=0.;

    void printarray();
    int layers();
    float point_arrh();
    float point_visc();

    int e,i,p,q,n,nel,k;

    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];
    const int lev = E->mesh.levmax;
    const int nno = E->lmesh.nno;
    const int noz = E->lmesh.noz;
    const int vpts = vpoints[dims];

    nel = E->lmesh.nel;

    visits++;

    
    if (E->viscosity.edot31in2Davail == 2 || 
        E->viscosity.edot31in2Davail == 3)
    {  if (visits==1) /* Non-dim. vplate3in2Dval */
       {  /* transform to m/sec */
          E->viscosity.vplate3in2Dval *= (0.01/3.155734e7);
          velo_scale = E->data.therm_diff/(E->data.layer_km*1.0e3);
          /* non-dimensionalize */
          E->viscosity.vplate3in2Dval /= velo_scale; 
       }
       temp = 0.0;
       age    = E->control.plate_age + E->monitor.elapsed_time;
       fprintf(stderr, "age=%g\n", age);
       for(iz=1;iz<=noz;iz++)
       {  l      = layers(E,iz);
          z      = E->X[2][iz];
          if (z == 0.)
             Ths=0.;
          else if (age == 0.) 
             Ths=1.;
          else 
          {  simvar = 0.5*z/sqrt(age);
             Ths    = erf(simvar);
          }
          E->Have.arrh[iz] =
             point_arrh(E, E->viscosity.RHEOL, E->viscosity.N0[l-1],
                E->viscosity.E[l-1], E->viscosity.T[l-1], Ths);
          /*
          if (visits==1)
            fprintf(stderr,"n=%d, z=%g, Ths=%g, arrh=%g\n",iz, 
                   E->X[2][iz], Ths, E->Have.arrh[iz]);
          */
          if (iz>1)
             temp += (E->Have.arrh[iz]+E->Have.arrh[iz-1])
                     *(E->X[2][iz]-E->X[2][iz-1])*0.5;
       }
       temp  = -E->viscosity.vplate3in2Dval/temp;
       if (E->viscosity.calcstress==1) 
          Ataupown = temp;
       fprintf(stderr,"calcstress=%d, t=%g, temp=%g, Ataun=%g\n",
               E->viscosity.calcstress, E->monitor.elapsed_time, 
               temp, Ataupown);
       temp1 = 0.0;
       Ths=0.;
       zold=0.;
       for(iz=1;iz<=noz;iz++) /* iz=vertical node counter */
       {  if (iz>1)
             temp1 = (E->Have.arrh[iz]+E->Have.arrh[iz-1])*0.5;
          /* ielz = vertical element counter */
          ielz = iz-1;
          E->Have.edot31[ielz] = Ataupown*temp1;

          /* Testing if eta*edot31=stress31 is constant over depth: 
             29-6-2004 JvH: last param to point_visc = C = 0 added.
                This means that this routine does not yet take 
                composition-dependent viscosity into account 
          if (iz>1)
          {  l      = layers(E,iz);
             Thsold=Ths;
             zold = z;
             z      = E->X[2][iz];
             if (z == 0.)
                Ths=0.;
             else if (age == 0.) 
                Ths=1.;
             else 
             {  simvar = 0.5*z/sqrt(age);
                Ths    = erf(simvar);
             }
             ZZ=(zold+z)/2.;
             E->control.eedot_available=1;
             TT=(Thsold + Ths)/2.;
             visc = point_visc(E, E->viscosity.RHEOL,
                           E->viscosity.N0[l-1],
                           E->viscosity.E[l-1], E->viscosity.Z[l-1],
                           E->viscosity.T[l-1], E->viscosity.dN[l-1],
                           E->viscosity.Ttransition,
                           TT, dummy, ZZ, -E->Have.edot31[ielz], 
                           E->viscosity.sdepv_expt[l-1], 0.);
             stress=visc*E->Have.edot31[ielz];
             fprintf(stderr, "edot=%g, visc=%g, stress=%g\n", 
                              E->Have.edot31[ielz], visc, stress);
          }
          */

          if (visits==1)
             fprintf(stderr,"1: e=%d edot31=%g\n",
                  ielz,E->Have.edot31[ielz]);
          /* 
          */ 
       }
       /* 
          fprintf(stderr,"e=%d edot31=%g\n",ielz,E->Have.edot31[ielz]);
       */
    }

    for(e=1;e<=nel;e++) 
    {  for(i=1;i<=ends;i++)   
       {  n=E->ien[e].node[i];
             VV[1][i] = E->V[1][n];
             VV[2][i] = E->V[2][n];
             if(dims==3) VV[3][i] = E->V[3][n];
             /*
             if(dims==2)
                fprintf(stderr,"V[%d]=[%g %g]\n",
                     i, VV[1][i], VV[2][i]);
             else
                fprintf(stderr,"V[%d]=[%g %g %g]\n",
                     i, VV[1][i], VV[2][i], VV[3][i]);
             */
       }
 
       for(p=1;p<=dims;p++)
          for(q=1;q<=dims;q++)
             dudx[p][q] = 0.0;
 

       for(i=1;i<=ends;i++)
          for(p=1;p<=dims;p++)
             for(q=1;q<=dims;q++)
                dudx[p][q] += 
                    VV[p][i]*E->gNX[e].ppt[GNPXINDEX(q-1,i,1)];
                   
       for(p=1;p<=dims;p++)
         for(q=1;q<=dims;q++)
             edot[p][q] = dudx[p][q] + dudx[q][p];

 
       if (dims==2)
       {  EEDOT[e] = edot[1][1]*edot[1][1] + edot[2][2]*edot[2][2]
                   + edot[1][2]*edot[1][2]*2.0;
          if (E->viscosity.edot31in2Davail == 1)
             /* Adding constant edot31 component to mimic out-of-plane
                strainrate, which affects nonlinear rheology */
             EEDOT[e] += 2.0*E->viscosity.edot31in2Dval
                            *E->viscosity.edot31in2Dval;
          else if (E->viscosity.edot31in2Davail == 2 ||
                   E->viscosity.edot31in2Davail == 3)
          {  /* Calc. edot31(z) that corresponds to (analytic)
                halfspace cooling profile */
             /* Determine element row in z-dir */
             ielz=(e-1)%E->lmesh.elz+1;
             /*
             fprintf(stderr,"1: e=%d ielz=%d, edot31=%g\n",
                  e, ielz,E->Have.edot31[ielz]);
             */
             EEDOT[e] += 2.0*E->Have.edot31[ielz]*E->Have.edot31[ielz];
          }
             
       }
 
       else if (dims==3)
          EEDOT[e] = edot[1][1]*edot[1][1] + edot[1][2]*edot[1][2]*2.0
                   + edot[2][2]*edot[2][2] + edot[2][3]*edot[2][3]*2.0
                   + edot[3][3]*edot[3][3] + edot[1][3]*edot[1][3]*2.0;
    }

/* Testing 
    fprintf(stderr,"tst2: %g\n",tst);
    fprintf(stderr,"tst2: %p\n",&tst);
    fprintf(stderr,"dudx2: %p\n",&dudx[1][1]);
    fprintf(stderr, "2:tst=%g\n", tst);
    fprintf(stderr,"tst3: %g\n",tst);
    fprintf(stderr,"tst3: %p\n",&tst);
    fprintf(stderr,"dudx3: %p\n",&dudx[1][1]);
    printarray(E, "eedot", EEDOT, nel);
    n=E->ien[40].node[1];
    fprintf(stderr, "2:n=%d, E->Vs=%g %g\n", n, E->V[1][n], E->V[2][n]);
    fprintf(stderr, "2:VVs1=%g %g\n", VV[1][1], VV[2][1]);
    fprintf(stderr, "2:VVs2=%g %g\n", VV[1][2], VV[2][2]);
    fprintf(stderr, "2:VVs3=%g %g\n", VV[1][3], VV[2][3]);
    fprintf(stderr, "2:VVs4=%g %g\n", VV[1][4], VV[2][4]);
    fprintf(stderr, "2:dudx=%g %g %g\n", dudx[1][1], dudx[2][2],
                                            dudx[1][2]);
    fprintf(stderr, "2:edot=%g %g %g\n", edot[1][1], edot[2][2],
                                            edot[1][2]);
    fprintf(stderr, "2:e256=%g\n", EEDOT[256]);
    fprintf(stderr, "2:calc=%g\n", edot[1][1]*edot[1][1] + 
                                   edot[2][2]*edot[2][2] +
                                   edot[1][2]*edot[1][2]*2.0);
    fprintf(stderr, "2:dims,nel,ends=%d %d %d\n", dims, nel, ends);
*/

    if(SQRT)
        for(e=1;e<=nel;e++)
            EEDOT[e] = (float)(sqrt((double)(0.5 *EEDOT[e])));
    else
        for(e=1;e<=nel;e++)
            EEDOT[e] *=  0.5;

    return;
}
