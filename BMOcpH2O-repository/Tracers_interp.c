/* Tracers_interp.c: tracer <-> element/node interpolation, locating tracers in
 * elements, nodal field (velocity/T) to tracer interpolation, and conversion of
 * tracer composition into element compositional buoyancy (compofuncbuoy). */

#include <signal.h>
#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include <string.h>
#include <stdlib.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "prem_density.h"


void tracinterp_t2e(E)
   struct All_variables *E;
{
   void all_trac_to_elem();
   void compofuncbuoy();
   void printtrac();

   /* All tracer functions with E->tracers.interp_t2e - value = 1
    * will be interpolated to the FEM mesh */
   all_trac_to_elem(E);
   /* In case of compos. buoyancy, E->tracers.comporayl is updated */
   compofuncbuoy(E);
}

void tracinterp_e2t(E)
   struct All_variables *E;
{
  int ifunc,itrac;
  float xi[4];
  float func_to_trac();

   for (ifunc=0; ifunc<E->tracers.ntracfunc; ifunc++)
   {  if (E->tracers.interp_e2t[ifunc] == 1)
       for (itrac=0; itrac<E->tracers.lntrac; itrac++) 
         E->tracers.ftrac[ifunc][itrac]=func_to_trac(E,
	           E->tracers.xtrac,E->T,xi,E->tracers.elemnr,itrac);
   }
   return;
}

void all_trac_to_elem(E)
   struct All_variables *E;
{
   int ielem, nelem, ifunc;
   
   float trac_to_elem();

   nelem=E->lmesh.NEL[E->mesh.levmax];

   for (ifunc=0; ifunc<E->tracers.ntracfunc; ifunc++)
   {  if (E->tracers.interp_t2e[ifunc] == 1)
      {   for (ielem=1; ielem<=nelem; ielem++)
              E->tracers.felem[ifunc][ielem] = 
                trac_to_elem(E,ifunc,ielem, E->tracers.ntrac_in_elem,
                             E->tracers.tracs_in_elem);
      }
   }      
   return;
}


float trac_to_elem(E, ifunc, ielem, ntrac_in_elem, tracs_in_elem)
   struct All_variables *E;
   int ifunc, ielem, *ntrac_in_elem, *tracs_in_elem;
{
   int itrac,lev=E->mesh.levmax;
   int nie1,nie2,nie_thiselem;
   float elemval, tracval;

   void tracexit();

   nie1=ntrac_in_elem[ielem-1];
   nie2=ntrac_in_elem[ielem];
   nie_thiselem = nie2 - nie1;
   elemval=0.;
   if (nie_thiselem > 0) 
   {  for (itrac=nie1; itrac<nie2; itrac++)
      {  tracval = E->tracers.ftrac[ifunc][tracs_in_elem[itrac]];
         /* fprintf(stderr, " tracer %d in elem %d: val=%g\n",
                 tracs_in_elem[itrac], ielem, tracval); */
         elemval += tracval;
      }
      elemval /= (nie2-nie1);
      /* fprintf(stderr, " elem %d: val=%g\n", ielem, elemval); */
   }
   else if ((ielem%E->lmesh.ELZ[lev])<E->tracers.bot_inflow)
    { //fprintf(stderr,"tracval undetermined: no tracers found in elem %d\n! ",ielem);

     if (E->tracers.Warn==1)
     { tracexit(E,"no tracers found in elem");
     }
    }

   return elemval;
}


void tracinterp_t2n(E)  /* works just for 2 dim at the moment */
   struct All_variables *E;
{
   void all_trac_to_node();
   void temperatures_conform_bcs();

   /* All tracer functions with E->tracers.interp_t2e - value = 2
    * will be interpolated to the FEM mesh */
   all_trac_to_node(E);
   temperatures_conform_bcs(E);
}


void all_trac_to_node(E)
   struct All_variables *E;
{
   int ifunc, nnode, inode;
   
   float trac_to_node();

   nnode=E->lmesh.nno;

   for (ifunc=0; ifunc<E->tracers.ntracfunc; ifunc++)
   {  if (E->tracers.interp_t2e[ifunc] == 2) 
      {   for (inode=1; inode<=nnode; inode++)
              E->T[inode] = 
                trac_to_node(E,ifunc,inode, E->tracers.ntrac_in_elem,
                             E->tracers.tracs_in_elem,E->tracers.xtrac);
      }
   }      
   return;
}


float trac_to_node(E, ifunc, inode, ntrac_in_elem, tracs_in_elem, xtrac)
   struct All_variables *E;
   int ifunc, inode, *ntrac_in_elem, *tracs_in_elem;
   float **xtrac;
{
   int itrac,level,ielem;
   int nie1, nie2, nie_thiselem, this_element, number_of_elements;
   float tracval;
   float node_value,node_weight,weight,trac_X,trac_Z,node_X,node_Z,Xsize,Zsize;

   void tracexit();

   level=E->mesh.levmax;
   number_of_elements=E->NEI[level].nels[inode];
/*   fprintf(stderr,"inode=%d, noe=%d\n", inode,number_of_elements);
   for (ielem=0; ielem<20; ielem++)
       { fprintf(stderr,"inode=%d, nels=%d \n",ielem,E->NEI[level].nels[ielem]);
       }
   for (ielem=0; ielem<20; ielem++)
       { fprintf(stderr,"index=%d, NEI=%d \n",ielem,E->NEI[level].element[ielem]);
       }*/
   node_X=E->X[1][inode];
   node_Z=E->X[2][inode];
   node_value=0.;
   node_weight=0.;
   for (ielem=1; ielem<=number_of_elements; ielem++) 
      {
      /* ntrac_in_elem[index] is the sum of number of tracers */
      /* in each elements until the index */
      this_element=E->NEI[level].element[4*(inode-1)+ielem-1];
      /*fprintf(stderr,"inode=%d, ielem=%d, level=%d, thisel=%d\n", inode,ielem,level,this_element);*/
      Xsize=E->ECO[level][this_element].size[1];
      Zsize=E->ECO[level][this_element].size[2];
      nie1 = ntrac_in_elem[this_element-1]; 
      nie2 = ntrac_in_elem[this_element];
      nie_thiselem = nie2 - nie1;
      if (nie_thiselem > 0)
         {  for (itrac=nie1; itrac<nie2; itrac++)
	    {  tracval = E->tracers.ftrac[ifunc][tracs_in_elem[itrac]];
	    trac_X=xtrac[1][tracs_in_elem[itrac]];
	    trac_Z=xtrac[2][tracs_in_elem[itrac]];
	    /* fprintf(stderr, " tracer %d in elem %d: val=%g\n",
	       tracs_in_elem[itrac], ielem, tracval); */
	    weight=fabs(Xsize-fabs(node_X-trac_X))*(Zsize-fabs(node_Z-trac_Z));
	    node_weight+=weight;
	    node_value += tracval*weight;
	    }
	 /* fprintf(stderr, " elem %d: val=%g\n", ielem, elemval); */
	 }
      else 
	{ //fprintf(stderr,"tracval undetermined: no tracers found in elem %d\n! ",ielem);
	 if (E->tracers.Warn==1)
	 { tracexit(E,"no tracers found in elem");
	 }
      }

    }
   node_value /= node_weight;
     
   return node_value;
}


void find_trac_in_elem(E, xtrac, elemnr, lntrac, iprintcheck)
   struct All_variables *E;
   float **xtrac; 
   int *elemnr, *lntrac, iprintcheck;
{
   /* Determine for each tracer the element in which it is located.
    * Note that the procedure here assumes a orthogonal mesh */
   int itrac, ix, iy, iz, nox, noy, noz;
   int lev,ifunc;
   //   float xmin,xmax,ymin,ymax,zmin,zmax;

   int bisect();
   float trac_to_elem();
   void tracexit();

   /* Find element number in x- and z-dir 
    * 1st element has elemnr 1.
    */
   lev = E->mesh.levmax;
   nox=E->lmesh.NOX[lev];
   noz=E->lmesh.NOZ[lev];
  
   if (E->mesh.nsd == 2) 
   {  for (itrac=0; itrac<*lntrac; itrac++)
       { 
	 // if (E->segment.xlayers==1) 
	 //   ix = (int)((nox-1)*(xtrac[1][itrac]/E->mesh.layer[1])+1); 
	 // else
	    ix = bisect(xtrac[1][itrac],E->XL[1],nox);
	 // if (E->segment.zlayers==1)
	 //   iz = (int)((noz-1)*(xtrac[2][itrac]/E->mesh.layer[2])+1); 
	 // else
	    iz = bisect(xtrac[2][itrac],E->XL[2],noz); 
         if (ix==0 || iz==0) 
	   {elemnr[itrac] = 0;
            if (iprintcheck==1) /* Print error message to stdin */
             { fprintf(stderr,"Proc %d: Trac %d (%g,%g,%g) is outside domain\n",
                   E->parallel.me, itrac, 
                   xtrac[1][itrac],xtrac[2][itrac],xtrac[3][itrac]);
                /* Now remove this tracer by lowering total nr of tracers */
                (*lntrac)--; 
                /* ... and overprint info wrong tracer with last tracer */
                xtrac[1][itrac]=xtrac[1][*lntrac]; 
                xtrac[2][itrac]=xtrac[2][*lntrac];
                for (ifunc=0; ifunc < E->tracers.ntracfunc;ifunc++)
                E->tracers.ftrac[ifunc][itrac]=E->tracers.ftrac[ifunc][*lntrac];
                itrac--;  /* Test this new tracer info again */
             }
	   }
	 else
         elemnr[itrac] = iz + (ix-1)*(noz-1);
      }
   }
   else if (E->mesh.nsd == 3) 
   {  noy=E->lmesh.NOY[E->mesh.levmax];
      for (itrac=0; itrac<*lntrac; itrac++)
	{
	 //  if (E->segment.xlayers==1) 
	 //   ix = (int)((nox-1)*(xtrac[1][itrac]/E->mesh.layer[1])+1); 
	 // else
	    ix = bisect(xtrac[1][itrac],E->XL[1],nox);
	 // if (E->segment.zlayers==1)
	 //   iz = (int)((noz-1)*(xtrac[2][itrac]/E->mesh.layer[2])+1); 
	 // else
	    iz = bisect(xtrac[2][itrac],E->XL[2],noz); 
	 // if (E->segment.ylayers==1)
	 //   iy = (int)((noy-1)*(xtrac[3][itrac]/E->mesh.layer[3])+1);
	 // else
            iy = bisect(xtrac[3][itrac],E->XL[3],noy); 
         if (ix==0 || iy==0 || iz==0) 
         {  /* Outside proc.bnds.: set elemnr to zero */
            elemnr[itrac] = 0;
            if (iprintcheck==1) /* Print error message to stdin */
            {  fprintf(stderr,"Proc %d: Trac %d (%g,%g,%g) is outside domain\n",
                   E->parallel.me, itrac, 
                   xtrac[1][itrac],xtrac[2][itrac],xtrac[3][itrac]);
                /* Now remove this tracer by lowering total nr of tracers */
                (*lntrac)--; 
                /* ... and overprint info wrong tracer with last tracer */
                xtrac[1][itrac]=xtrac[1][*lntrac]; 
                xtrac[2][itrac]=xtrac[2][*lntrac];
                xtrac[3][itrac]=xtrac[3][*lntrac];
                for (ifunc=0; ifunc < E->tracers.ntracfunc;ifunc++)
                E->tracers.ftrac[ifunc][itrac]=E->tracers.ftrac[ifunc][*lntrac];
                itrac--;  /* Test this new tracer info again */
            }
         }
         else
         {  elemnr[itrac] = iz+(ix-1)*(noz-1)+(iy-1)*(noz-1)*(nox-1);
         }
      }
   }
// for (itrac=0; itrac<lntrac; itrac++)
// { fprintf(stderr,"elemnr[%d]=%d, val2=%g\n",itrac,elemnr[itrac],
//                                    E->tracers.ftrac[2][itrac]);
// }

   /* Orthogonal mesh test: to save CPU perform only once. Note that if
    * remeshing takes place, this test must be redone! *
   if (visits==0)
   {  if (E->mesh.nsd == 2) 
      {  for (itrac=0; itrac<lntrac; itrac++)
         {  node = E->IEN[lev][elemnr[itrac]].node[1];
            xmin= E->X[1][node];
            zmin= E->X[2][node];
            node = E->IEN[lev][elemnr[itrac]].node[3];
            xmax= E->X[1][node];
            zmax= E->X[2][node];
	    if (xmin > xtrac[1][itrac] ||
                xmax < xtrac[1][itrac] ||
                zmin > xtrac[2][itrac] ||
                zmax < xtrac[2][itrac])
            {  fprintf(E->fp, "x-crd:%g<=%g<%g and z-crd:%g<=%g<%g? No   %i\n",
                  xmin, xtrac[1][itrac], xmax,
		       zmin, xtrac[2][itrac], zmax,itrac);
               tracexit(E,"Error find_trac_in_elem: 2-D mesh orthogonal?");
	       }
         }
      }
      else if (E->mesh.nsd == 3) 
      {  for (itrac=0; itrac<lntrac; itrac++)
         {  node = E->IEN[lev][elemnr[itrac]].node[1];
            xmin= E->X[1][node];
            zmin= E->X[2][node];
            ymin= E->X[3][node];
            node = E->IEN[lev][elemnr[itrac]].node[7];
            xmax= E->X[1][node];
            zmax= E->X[2][node];
            ymax= E->X[3][node];
            if (xmin > xtrac[1][itrac] ||
                xmax < xtrac[1][itrac] ||
                zmin > xtrac[2][itrac] ||
                zmax < xtrac[2][itrac] ||
                ymin > xtrac[3][itrac] ||
                ymax < xtrac[3][itrac])  
            {  fprintf(E->fp, "x-crd:%g<=%g<%g, z-crd:%g<=%g<%g",
                  xmin, xtrac[1][itrac], xmax,
                  zmin, xtrac[2][itrac], zmax);
               fprintf(E->fp, " and y-crd:%g<=%g<%g? No  \n",
		       ymin, xtrac[3][itrac], ymax);
               tracexit(E,"Error find_trac_in_elem: 3-D mesh orthogonal?");
            }
         }
      }
   }

   visits ++;
   */
   return;
}

void update_trac_tables(E, elemnr, ntrac_in_elem,
                       tracs_in_elem, lntrac)
   struct All_variables *E;
   int *elemnr, *ntrac_in_elem, *tracs_in_elem, lntrac;
{
   int lev;
   int ielem, nelem, itrac, jindex, iindex;

   lev = E->mesh.levmax;
   nelem=E->lmesh.NEL[lev];
   /* Make table of # tracers per element */
   /* Reset ntrac_in_elem array */
   for (ielem=0; ielem<=nelem; ielem++)
      ntrac_in_elem[ielem]=0;
   /* Count # tracers in each elem */
   for (itrac=0; itrac<lntrac; itrac++)
   {  ntrac_in_elem[elemnr[itrac]]++;
      // fprintf(stderr,"adding tracer %d to element %d\n",itrac,
      //    elemnr[itrac]); 
   }
   /* Make ntrac_in_elem array cumulative to make it useful index for 
    * tracs_in_elem (see below) */
   for (ielem=1; ielem<=nelem; ielem++)
   {  ntrac_in_elem[ielem]+=ntrac_in_elem[ielem-1];
   }
// for (ielem=1; ielem<=nelem; ielem++)
// {  fprintf(stderr,"nie[%d]=%d\n",ielem,ntrac_in_elem[ielem]);
// }

   /* make table of trac.nrs. per element */
   for (ielem=0; ielem<=nelem; ielem++)
      E->tracers.elemtraccounter[ielem]=0;
   for (itrac=0; itrac<lntrac; itrac++)
   {  ielem = elemnr[itrac];
      
      jindex = E->tracers.elemtraccounter[ielem];
      if (ielem!=0) iindex = ntrac_in_elem[ielem-1] + jindex;
      else iindex=jindex;   
      tracs_in_elem[iindex] = itrac;
      E->tracers.elemtraccounter[ielem]++;
   }
// for (itrac=0; itrac<lntrac; itrac++)
//    fprintf(stderr,"tie[%d]=%d\n",itrac,tracs_in_elem[itrac]); 

   return;
}

void vel_to_trac(E, vel, xtrac, vtrac, elemnr, lntrac)
   struct All_variables *E;
   float **vel, **xtrac, **vtrac; 
   int *elemnr, lntrac;
{
   /* Locate each tracer's element, then interpolate the nodal velocity to the
    * tracer (assumes an orthogonal mesh). */

   float xmin, ymin, zmin, xmax, ymax, zmax;
   float xtr[4], xi[4], N[9];
   double  poly;
   int iel, itrac, gnode[9], gn, lnode, lev, d;

   double lpoly();
   void find_trac_in_elem();

   lev = E->mesh.levmax;
 
   /* Determine and store the element nr of each tracer */
   find_trac_in_elem(E, xtrac, elemnr, &(lntrac), 0);

   /* Interpolate v from surrounding nodal points to tracers */
   if (E->mesh.nsd == 3) /* 3-D */
   {  for (itrac=0; itrac<lntrac; itrac++) 
      {  iel = elemnr[itrac];
         if (iel > 0) /* trac inside domain: interp. v */
         {  for(lnode=1;lnode<=enodes[E->mesh.nsd];lnode++)  
            {   gnode[lnode] = E->IEN[lev][iel].node[lnode];
            }
            /* Determine local coordinates of tracer */
            xmin = E->X[1][gnode[1]];
            xmax = E->X[1][gnode[4]];
            xtr[1] = xtrac[1][itrac];
            xi[1] = -1. + 2*(xtr[1]-xmin)/(xmax-xmin);
            zmin = E->X[2][gnode[1]];
            zmax = E->X[2][gnode[2]];
            xtr[2] = xtrac[2][itrac];
            xi[2] = -1. + 2*(xtr[2]-zmin)/(zmax-zmin);
            ymin = E->X[3][gnode[1]];
            ymax = E->X[3][gnode[5]];
            xtr[3] = xtrac[3][itrac];
            xi[3] = -1. + 2*(xtr[3]-ymin)/(ymax-ymin);
            /* fprintf(E->fp, "itrac=%d, (xi, eta, dzeta)=(%g,%g,%g)\n",
             *       itrac, xi[1], xi[2], xi[3]);
             */
            for(lnode=1;lnode<=enodes[E->mesh.nsd];lnode++)  
            {  N[lnode] = 1.;
               for(d=1;d<=E->mesh.nsd;d++)  
               {  poly = lpoly(bb[d-1][lnode],xi[d]);
                  N[lnode] *= (float)(poly);
               }
               /* fprintf(E->fp, " N[%d]=%g\n", lnode, N[lnode]);
                */
            }
            vtrac[1][itrac] = 0.;
            vtrac[2][itrac] = 0.;
            vtrac[3][itrac] = 0.;
            for(lnode=1;lnode<=enodes[E->mesh.nsd];lnode++)  
            {  gn = gnode[lnode];
               vtrac[1][itrac] += N[lnode]*vel[1][gn];
               vtrac[2][itrac] += N[lnode]*vel[2][gn];
               vtrac[3][itrac] += N[lnode]*vel[3][gn];
            }
            /*
            fprintf(E->fp,
               "itr=%d,(x,z,y)=(%g,%g,%g),(vx,vz,vy)=(%g,%g,%g)\n",
               itrac, xtrac[1][itrac], xtrac[2][itrac], xtrac[3][itrac], 
                      vtrac[1][itrac], vtrac[2][itrac], vtrac[3][itrac]);
             */
          }
          else /* trac outside domain: give weird v to indicate that */
          {  vtrac[1][itrac] = 1e40;
             vtrac[2][itrac] = 1e40;
             vtrac[3][itrac] = 1e40;
          }
      }
   }
   else  /* 2-D */
   {  for (itrac=0; itrac<lntrac; itrac++) 
      {  iel = elemnr[itrac];
         if (iel > 0) /* trac inside domain: interp. v */
         {  for(lnode=1;lnode<=enodes[E->mesh.nsd];lnode++)  
            {   gnode[lnode] = E->IEN[lev][iel].node[lnode];
            }
            /* Determine local coordinates of tracer */
            xmin = E->X[1][gnode[1]];
            xmax = E->X[1][gnode[4]];
            xtr[1] = xtrac[1][itrac];
            xi[1] = -1. + 2*(xtr[1]-xmin)/(xmax-xmin);
            zmin = E->X[2][gnode[1]];
            zmax = E->X[2][gnode[2]];
            xtr[2] = xtrac[2][itrac];
            xi[2] = -1. + 2*(xtr[2]-zmin)/(zmax-zmin);
            /* fprintf(E->fp, "itrac=%d, (xi, eta)=(%g,%g)\n",
             *       itrac, xi[1], xi[2]);
             */
            for(lnode=1;lnode<=enodes[E->mesh.nsd];lnode++)  
            {  N[lnode] = 1.;
               for(d=1;d<=E->mesh.nsd;d++)  
               { poly = lpoly(bb[d-1][lnode],xi[d]);
                 N[lnode] *= (float)(poly);
               }
               /* fprintf(E->fp, " N[%d]=%g\n", lnode, N[lnode]);
                */
            }
            vtrac[1][itrac] = 0.;
            vtrac[2][itrac] = 0.;
            for(lnode=1;lnode<=enodes[E->mesh.nsd];lnode++)  
            {  gn = gnode[lnode];
               vtrac[1][itrac] += N[lnode]*vel[1][gn];
               vtrac[2][itrac] += N[lnode]*vel[2][gn];
            }
            /*
            fprintf(E->fp, "itrac=%d,(x,z)=(%g,%g),(vx,vz)=(%g,%g)\n",
              itrac, xtrac[1][itrac], xtrac[2][itrac],
                     vtrac[1][itrac], vtrac[2][itrac]);
            */
         }
         else /* trac outside domain: give weird v to indicate that */
	 {  
	    vtrac[1][itrac] = 1e40;
            vtrac[2][itrac] = 1e40;
         }
      } 
   }
}

float func_to_trac(E, xtrac, field, xi, elemnr, itrac)
   struct All_variables *E;
   float **xtrac, *field, *xi; 
   int *elemnr, itrac;
{
   /* Interpolate a nodal field (e.g. T) to the tracer location
    * (assumes an orthogonal mesh). */

  float xmin, ymin, zmin, xmax, ymax, zmax, ftrac;
   float xtr[4], N[9];
   double  poly;
   int iel, gnode[9], gn, lnode, lev, d;

   double lpoly();
   
   lev = E->mesh.levmax;
 
   /* Interpolate v from surrounding nodal points to tracers */
   if (E->mesh.nsd == 3) /* 3-D */
   {  
         iel = elemnr[itrac];
         if (iel > 0) /* trac inside domain: interp. f */
         {  for(lnode=1;lnode<=enodes[E->mesh.nsd];lnode++)  
            {   gnode[lnode] = E->IEN[lev][iel].node[lnode];
            }
            /* Determine local coordinates of tracer */
            xmin = E->X[1][gnode[1]];
            xmax = E->X[1][gnode[4]];
            xtr[1] = xtrac[1][itrac];
            xi[1] = -1. + 2*(xtr[1]-xmin)/(xmax-xmin);
            zmin = E->X[2][gnode[1]];
            zmax = E->X[2][gnode[2]];
            xtr[2] = xtrac[2][itrac];
            xi[2] = -1. + 2*(xtr[2]-zmin)/(zmax-zmin);
            ymin = E->X[3][gnode[1]];
            ymax = E->X[3][gnode[5]];
            xtr[3] = xtrac[3][itrac];
            xi[3] = -1. + 2*(xtr[3]-ymin)/(ymax-ymin);
            /* fprintf(E->fp, "itrac=%d, (xi, eta, dzeta)=(%g,%g,%g)\n",
             *       itrac, xi[1], xi[2], xi[3]);
             */
            for(lnode=1;lnode<=enodes[E->mesh.nsd];lnode++)  
            {  N[lnode] = 1.;
               for(d=1;d<=E->mesh.nsd;d++)  
               {  poly = lpoly(bb[d-1][lnode],xi[d]);
                  N[lnode] *= (float)(poly);
               }
               /* fprintf(E->fp, " N[%d]=%g\n", lnode, N[lnode]);
                */
            }
            ftrac = 0.;
            for(lnode=1;lnode<=enodes[E->mesh.nsd];lnode++)  
            {  gn = gnode[lnode];
               ftrac += N[lnode]*field[gn];
            }
            /*
            fprintf(E->fp,
               "itr=%d,(x,z,y)=(%g,%g,%g),func=%g\n",
               itrac, xtrac[1][itrac],xtrac[2][itrac],xtrac[3][itrac], 
                      ftrac);
             */
          }
          else /* trac outside domain: give weird ftrac to indicate that */
          {  ftrac = 1e40;
          }
      
   }
   else  /* 2-D */
   {   
         iel = elemnr[itrac];
         if (iel > 0) /* trac inside domain: interp. v */
         {  for(lnode=1;lnode<=enodes[E->mesh.nsd];lnode++)  
            {   gnode[lnode] = E->IEN[lev][iel].node[lnode];
            }
            /* Determine local coordinates of tracer */
            xmin = E->X[1][gnode[1]];
            xmax = E->X[1][gnode[4]];
            xtr[1] = xtrac[1][itrac];
            xi[1] = -1. + 2*(xtr[1]-xmin)/(xmax-xmin);
            zmin = E->X[2][gnode[1]];
            zmax = E->X[2][gnode[2]];
            xtr[2] = xtrac[2][itrac];
            xi[2] = -1. + 2*(xtr[2]-zmin)/(zmax-zmin);
            /* fprintf(E->fp, "itrac=%d, (xi, eta)=(%g,%g)\n",
             *       itrac, xi[1], xi[2]);
             */
            for(lnode=1;lnode<=enodes[E->mesh.nsd];lnode++)  
            {  N[lnode] = 1.;
               for(d=1;d<=E->mesh.nsd;d++)  
               { poly = lpoly(bb[d-1][lnode],xi[d]);
                 N[lnode] *= (float)(poly);
               }
               /* fprintf(E->fp, " N[%d]=%g\n", lnode, N[lnode]);
                */
            }
            ftrac = 0.;
            for(lnode=1;lnode<=enodes[E->mesh.nsd];lnode++)  
            {  gn = gnode[lnode];
               ftrac += N[lnode]*field[gn];
            }
            /*
            fprintf(E->fp, "itrac=%d,(x,z)=(%g,%g),func=%g\n",
              itrac, xtrac[1][itrac], xtrac[2][itrac], ftrac);
            */
         }
         else /* trac outside domain: give weird ftrac to indicate that */
         {  ftrac = 1e40;
         }
       
   }
   return ftrac;
}

int bisect(xx, X, n)
   int n;
   float xx, *X;
{
   /* Bisection: returns m with X[m] <= xx < X[m+1] (X filled 1..n); 0 if out of range. */
   
   int mmin, mmax, mmid;
   
   mmin=1;
   mmax=n;
   mmid=(mmin+mmax)/2;
   if (xx < X[mmin] || xx > X[mmax])
   {  /* out of range */
      return 0;
   }
   while ((mmax-mmin) > 1)
   {  
      /*
      fprintf(stderr,"X[%d]=%g, X[%d]=%g\n",mmin,X[mmin],mmax,X[mmax]);
      */
      if (xx < X[mmid])
         mmax=mmid;
      else
         mmin=mmid;
      mmid=(mmin+mmax)/2;
   }
   return mmin;
}

void compofuncbuoy(E)
   struct All_variables *E;
{
   int iz, ielem, noz, nelem, lev;
   float depth660, depth, depth300, depth410;
   static int been_here=0;
   static int print=0;
   static float SiO2_density[4099];
   void tracexit();

   nelem=E->lmesh.NEL[E->mesh.levmax];

   /* BMO_density_mode == 1: Du et al. (2019) density from FeO (tracfunc 80) and H2O (tracfunc 81):
    *   rho = rho_nv(z) * (1 + X_Fe*(0.28 + 0.04*n)) * (1 + alpha*X_H2O)
    *   compobuoy[el] = (rho - data.density)*Atemp / (data.density*therm_exp*ref_temperature)
    * The melt/solid boundary is the initial front MO_depth_initial (z=0 surface, z=1 CMB):
    * z >= z_init is melt, z < z_init is solid (BMO_density_solid_above=1 swaps this).
    * Falls back to uniform XFe_ref / 0 if tracfunc 80/81 are not active. */

   if (E->control.BMO_density_mode == 1 && E->tracers.ntracbuoyfunc > 0)
   {
       /* Locate felem arrays for FeO (tracfunc 80) and H2O (tracfunc 81) */
       int ifunc_Fe  = -1, ifunc_H2O = -1, ifunc;
       for (ifunc = 0; ifunc < E->tracers.ntracfunc; ifunc++) {
           if (E->tracers.tracfunc[ifunc] == 80 && E->tracers.felem[ifunc] != NULL)
               ifunc_Fe  = ifunc;
           if (E->tracers.tracfunc[ifunc] == 81 && E->tracers.felem[ifunc] != NULL)
               ifunc_H2O = ifunc;
       }
       if (been_here == 0 && E->parallel.me == 0) {
           fprintf(stderr,
               "[BMO_density] Du2019 mode: z_init=%.4f  ifunc_Fe=%d  ifunc_H2O=%d"
               "  rho_nv_prem=%d\n",
               (double)E->control.MO_depth_initial, ifunc_Fe, ifunc_H2O,
               E->control.BMO_rho_nv_prem);
       }

       lev = E->mesh.levmax;
       noz = E->lmesh.NOZ[lev];
       int   elz       = noz - 1;             /* elements in z per column   */
       float rho_ref   = E->data.density;     /* CITCOM reference [kg/m^3]  */
       (void)E->control.BMO_rho_melt;
       float XFe_ref   = E->control.BMO_XFe_ref;
       float n_spin    = E->control.BMO_n_spin;
       float alpha_w   = E->control.BMO_alpha_water;
       float z_init    = E->control.MO_depth_initial;
       int solid_above = E->control.BMO_density_solid_above;
       float factor = E->control.Atemp
                      / (rho_ref * E->data.therm_exp * E->data.ref_temperature);

       for (ielem = 1; ielem <= nelem; ielem++) {
           /* Element z-index and centre coordinate */
           iz = (ielem - 1) % elz + 1;
           float z_ctr = 0.5f * (E->XL[2][iz] + E->XL[2][iz + 1]);

           /* Zone test: melt is z >= z_init (solid_above=0) or z <= z_init (solid_above=1). */
           int in_melt = (solid_above == 0) ? (z_ctr >= z_init)
                                            : (z_ctr <= z_init);

           float drho;
           float rho_nv = lookup_BMO_rho_nv(E, z_ctr);

           if (in_melt) {
               /* BMO melt: density from liquid FeO and H2O (tracfunc 80/81). */
               float XFe = (ifunc_Fe  >= 0)
                           ? E->tracers.felem[ifunc_Fe ][ielem] : XFe_ref;
               float XH2O= (ifunc_H2O >= 0)
                           ? E->tracers.felem[ifunc_H2O][ielem] : 0.0f;

               /* rho/rho_0 = 1 + X_Fe*(0.28 + 0.04*n); X_Fe = molar Fe#, n = BMO_n_spin */
               float rho_anhydrous = rho_nv * (1.0f + XFe * (0.28f + 0.04f * n_spin));

               /* Water correction from Du et al. (2019):
                * rho/rho_0 - 1 = alpha_w * X_H2O (mass fraction) */
               float rho_calc = rho_anhydrous * (1.0f + alpha_w * XH2O);

               drho = rho_calc - rho_ref;
           } else {
               /* Solid zone: same formula with solid-side FeO/H2O tracer values. */
               float XFe_s  = (ifunc_Fe  >= 0)
                              ? E->tracers.felem[ifunc_Fe ][ielem] : XFe_ref;
               float XH2O_s = (ifunc_H2O >= 0)
                              ? E->tracers.felem[ifunc_H2O][ielem] : 0.0f;
               /* Same formula as the melt zone */
               float rho_anhydrous_s = rho_nv * (1.0f + XFe_s * (0.28f + 0.04f * n_spin));
               float rho_calc_s = rho_anhydrous_s * (1.0f + alpha_w * XH2O_s);
               drho = rho_calc_s - rho_ref;
           }

           E->tracers.compobuoy[ielem] = drho * factor;
       }

       been_here++;
       return;
   }

   /* Legacy linear comporayl * C summation (BMO_density_mode == 0) */
   if (E->tracers.ntracbuoyfunc == 1)
   {  for (ielem=1; ielem<=nelem; ielem++)
      {  E->tracers.compobuoy[ielem] = 
            E->tracers.comporayl[0] * E->tracers.felem[0][ielem];
         if (been_here==0 && print==1)
	            { fprintf(stderr, "el=%d, C0=%g, buoy=%g   ", 
                ielem, E->tracers.felem[0][ielem],
                E->tracers.compobuoy[ielem]);
         }
      }
   }
   else if (E->tracers.ntracbuoyfunc == 2 && E->tracers.interp_e2t[0] == 22)
   {if (been_here==0)
      {
        lev = E->mesh.levmax;
        noz=E->lmesh.NOZ[lev];
        depth300=300./E->data.layer_km;
        depth410=410./E->data.layer_km;
        depth660=660./E->data.layer_km;
        for (iz=1; iz<noz; iz++)
         {
	   if (E->XL[2][iz+1]<=depth300) SiO2_density[iz]=500.;
           else if (E->XL[2][iz]<=depth300) SiO2_density[iz]=700.;
           else if (E->XL[2][iz+1]<depth410) SiO2_density[iz]=900.;
           else if (E->XL[2][iz]<depth410) SiO2_density[iz]=700.;
           else if (E->XL[2][iz+1]<depth660) SiO2_density[iz]=500.;
           else if (E->XL[2][iz]<depth660) SiO2_density[iz]=299.1;
           else
             { 
               depth=0.5*(E->XL[2][iz]+E->XL[2][iz+1])*E->data.layer_km;
               SiO2_density[ielem]=E->tracers.comporayl[0]*(176.42-0.12405*depth+1.044e-5*depth*depth);
             }
         }
      }
    for (ielem=1; ielem<=nelem; ielem++)
      {  E->tracers.compobuoy[ielem] = E->tracers.felem[0][ielem] * SiO2_density[ielem];  // depth-dep SiO2 compo-Ra#
         if (been_here==0 && ielem<noz)
         {  fprintf(stderr, "el=%d, C0=%g, SiO2=%g, buoy=%g   ", 
              ielem, E->tracers.felem[0][ielem],
              SiO2_density[ielem],
              E->tracers.compobuoy[ielem]);
         }
      }
   }
   else if (E->tracers.ntracbuoyfunc == 2)
   {  for (ielem=1; ielem<=nelem; ielem++)
      {  E->tracers.compobuoy[ielem] =
             E->tracers.comporayl[0]* E->tracers.felem[0][ielem]
           + E->tracers.comporayl[1]* E->tracers.felem[1][ielem];
         /*if (been_here==0)
         {  fprintf(stderr, "el=%d, C0=%g, C1=%g, buoy=%g   ", 
              ielem, E->tracers.felem[0][ielem],
              E->tracers.felem[1][ielem],
              E->tracers.compobuoy[ielem]);
         }*/
      }
   }      
   else if (E->tracers.ntracbuoyfunc == 3)
   {  for (ielem=1; ielem<=nelem; ielem++)
      {  E->tracers.compobuoy[ielem] =
             E->tracers.comporayl[0]* E->tracers.felem[0][ielem]
           + E->tracers.comporayl[1]* E->tracers.felem[1][ielem]
	   + E->tracers.comporayl[2]* E->tracers.felem[2][ielem];
      }
   }      
   else if (E->tracers.ntracbuoyfunc > 3)
   { fprintf(stderr, "Too many tracer buoyancy functions\n");
     fprintf(stderr, "Not yet implemented\n");
     tracexit(E,"Too many tracer buoyancy functions");
   }

   been_here++;

   return;
}
