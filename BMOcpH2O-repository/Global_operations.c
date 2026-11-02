#if 1
#include <mpi.h>
#endif

#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

/* ===============================================
   strips horizontal average from nodal field X. 
   Assumes orthogonal mesh, otherwise, horizontals
   aren't & another method is required.
   =============================================== */

void remove_horiz_ave(E,X,H,store_or_not)
     struct All_variables *E;
     float *X, *H;
     int store_or_not;

{
    int i,j,k,n,ln,nox,noz,noy;
    void return_horiz_ave();

    const int dims = E->mesh.nsd;

    nox = E->lmesh.nox;
    noy = E->lmesh.noy;
    noz = E->lmesh.noz;

    return_horiz_ave(E,X,H);

    for(i=1;i<=noz;i++)
      for(k=1;k<=noy;k++)
        for(j=1;j<=nox;j++)     {
	   n = i+(j-1)*noz+(k-1)*nox*noz;
	   X[n] -= H[i];  
	   }

   return; 
  }

void return_horiz_sum(E,X,H,nn)
     struct All_variables *E;
     float *X, *H;
     int nn;
{ 
  const int dims = E->mesh.nsd;
  int i,j,k,d,nint,noz,nox,noy,el,elz,elx,ely,j1,j2,i1,i2,k1,k2,nproc;
  int lnode[5], sizeofH, noz2,iroot;
  float *Have,*temp;

  int *processors;

  MPI_Comm world, horizon_p;
  MPI_Group world_g, horizon_g;

  processors = (int *)malloc((E->parallel.nprocxy+2)*sizeof(int));

       /* determine which processors should get the message from me for 
	       computing the layer averages */

  nproc = 0;
  for (j=0;j<E->parallel.nprocy;j++)
    for (i=0;i<E->parallel.nprocx;i++)  {
      d = E->parallel.me_loc[2] + i*E->parallel.nprocz + 
       		                  j*E->parallel.nprocxz;
      processors[nproc] =  d;
      nproc ++;
      }

  if (nproc>0)  {
    world = MPI_COMM_WORLD;
    MPI_Comm_group(world, &world_g);
    MPI_Group_incl(world_g, nproc, processors, &horizon_g);
    MPI_Comm_create(world, horizon_g, &horizon_p);

    MPI_Allreduce(X,H,nn,MPI_FLOAT,MPI_SUM,horizon_p);

    MPI_Comm_free (&horizon_p);
    MPI_Group_free(&horizon_g);
    MPI_Group_free(&world_g);

    }
  else
    for (i=0;i<nn;i++)
      H[i] = X[i];

  free ((void *) processors);
  
  return; 
  }


void return_horiz_ave(E,X,H)
     struct All_variables *E;
     float *X, *H;
{ 
  const int dims = E->mesh.nsd;
  int i,j,k,d,nint,noz,nox,noy,el,elz,elx,ely,j1,j2,i1,i2,k1,k2,nproc;
  int lnode[5], sizeofH, noz2,iroot;
  float *Have,*temp;
  struct Shape_function1 M;
  struct Shape_function1_dA dGamma;
  void get_global_1d_shape_fn();

  int *processors;

  MPI_Comm world, horizon_p;
  MPI_Group world_g, horizon_g;

  sizeofH = (2*E->lmesh.noz+2)*sizeof(float);

  processors = (int *)malloc((E->parallel.nprocxy+2)*sizeof(int));
  Have = (float *)malloc(sizeofH);
  temp = (float *)malloc(sizeofH);

  noz = E->lmesh.noz;
  noy = E->lmesh.noy;
  nox = E->lmesh.nox;
  elz = E->lmesh.elz;
  elx = E->lmesh.elx;
  ely = E->lmesh.ely;
  noz2 = 2*noz;

  for (i=1;i<=elz;i++)  {
    temp[i] = temp[i+noz] = 0.0;
    temp[i+1] = temp[i+1+noz] = 0.0;
    for (j=1;j<=elx;j++)
      for (k=1;k<=ely;k++)  {
 	    el = i + (j-1)*elz + (k-1)*elx*elz;
	    get_global_1d_shape_fn(E,el,&M,&dGamma,0);   

 	    lnode[1] = i+(j-1)*noz+(k-1)*nox*noz;
      	lnode[2] = i+j*noz+(k-1)*nox*noz;
    	lnode[3] = i+j*noz+k*nox*noz;
    	lnode[4] = i+(j-1)*noz+k*nox*noz;

    	for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
    	  for(nint=1;nint<=onedvpoints[E->mesh.nsd];nint++)   {
    	      temp[i] += X[lnode[d]] * E->M.vpt[GMVINDEX(d,nint)] 
			 * dGamma.vpt[GMVGAMMA(1,nint)];
    	      temp[i+noz] += E->M.vpt[GMVINDEX(d,nint)] 
			 * dGamma.vpt[GMVGAMMA(1,nint)]; 
    	      } 

        if (i==elz)  {
    	  lnode[1] = 1+i+(j-1)*noz+(k-1)*nox*noz;
    	  lnode[2] = 1+i+j*noz+(k-1)*nox*noz;
    	  lnode[3] = 1+i+j*noz+k*nox*noz;
    	  lnode[4] = 1+i+(j-1)*noz+k*nox*noz;

    	  for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
    	    for(nint=1;nint<=onedvpoints[E->mesh.nsd];nint++)   {
    	      temp[i+1] += X[lnode[d]] * E->M.vpt[GMVINDEX(d,nint)] 
			 * dGamma.vpt[GMVGAMMA(1,nint)];
    	      temp[i+1+noz] += E->M.vpt[GMVINDEX(d,nint)] 
			 * dGamma.vpt[GMVGAMMA(1,nint)]; 
    	      } 
    	  } 

    	}        /* Done one traverse */ 

     }        /* Done for i */ 


       /* determine which processors should get the message from me for 
	       computing the layer averages */

  nproc = 0;
  for (j=0;j<E->parallel.nprocy;j++)
    for (i=0;i<E->parallel.nprocx;i++)  {
      d = E->parallel.me_loc[2] + i*E->parallel.nprocz + 
       		                  j*E->parallel.nprocxz;
      processors[nproc] =  d;
      nproc ++;
      }

  if (nproc>0)  {
    world = MPI_COMM_WORLD;
    MPI_Comm_group(world, &world_g);
    MPI_Group_incl(world_g, nproc, processors, &horizon_g);
    MPI_Comm_create(world, horizon_g, &horizon_p);

    MPI_Allreduce(temp,Have,noz2+1,MPI_FLOAT,MPI_SUM,horizon_p);

    MPI_Comm_free (&horizon_p);
    MPI_Group_free(&horizon_g);
    MPI_Group_free(&world_g);

    }
  else
    for (i=1;i<=noz2;i++)
      Have[i] = temp[i];

  for (i=1;i<=noz;i++) {
    if(Have[i+noz] != 0.0) 
       H[i] = Have[i]/Have[i+noz];
    }

  free ((void *) Have);
  free ((void *) temp);
  free ((void *) processors);
  
  return; 
  }
	   

float return_bulk_value(E,Z,average) 
     struct All_variables *E;
     float *Z;
     int average;

{  
    void float_global_operation();

    int i,j,k,l,m,n,el,elx,ely,elz,i1,i2,j1,j2,k1,k2;
    float volume,integral,volume1,integral1;
   
    struct Shape_function GN;
    struct Shape_function_dx GNx;
    struct Shape_function_dA dOmega;
    
    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];

    volume1=0.0;
    integral1=0.0;
     
  elz = E->lmesh.elz;
  elx = E->lmesh.elx;
  ely = E->lmesh.ely;


  for (i=1;i<=elz;i++)  
    for (j=1;j<=elx;j++)
      for (k=1;k<=ely;k++)  {
	  el = i+(j-1)*elz+(k-1)*elz*elx;

	
	  for(l=1;l<=vpts;l++)
	    for(m=1;m<=ends;m++) {
		n = E->ien[el].node[m];
		volume1 += E->N.vpt[GNVINDEX(m,l)] * E->gDA[el].vpt[l];
		integral1 += Z[n] * E->N.vpt[GNVINDEX(m,l)] * E->gDA[el].vpt[l];
	    }

      }


    MPI_Allreduce(&volume1  ,&volume  ,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&integral1,&integral,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
     
    if(average && volume != 0.0)
 	   integral /= volume;

    return((float)integral);
}



double global_vdot(E,A,B,lev)
   struct All_variables *E;
   double *A,*B;
   int lev;

{
  int i,neq;
  double prod, temp;

  neq=E->lmesh.NEQ[lev];

  temp = 0.0;
  for (i=0;i<neq;i++)  {
    if (E->parallel.IDD[lev][i]) /* only get the sum from relevant ID */
      temp += A[i]*B[i]; 
    }

  MPI_Allreduce(&temp, &prod,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  return (prod);
}


double global_pdot(E,A,B,lev)
   struct All_variables *E;
   double *A,*B;
   int lev;

{
  int i,npno;
  double prod, temp;

  npno=E->lmesh.NPNO[lev];

  temp = 0.0;
  for (i=1;i<=npno;i++)
      temp += A[i]*B[i]; 

  MPI_Allreduce(&temp, &prod,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  return (prod);
  }


float global_tdot(E,A,B,lev)
   struct All_variables *E;
   float *A,*B;
   int lev;

{
  int i,nno;
  float prod, temp;

  nno=E->lmesh.NNO[lev];

  temp = 0.0;
  for (i=1;i<=nno;i++)
      temp += A[i]*B[i]; 

  MPI_Allreduce(&temp, &prod,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

  return (prod);
  }


float global_fmin(E,a)
   struct All_variables *E;
   float a;
{
  float temp;
  MPI_Allreduce(&a, &temp,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
  return (temp);
  }

float global_fmax(E,a)
   struct All_variables *E;
   float a;
{
  float temp;
  MPI_Allreduce(&a, &temp,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
  return (temp);
  }

int global_vmax( E, A, Amax, lev )
   struct All_variables *E;
   double *A, *Amax;
   int lev;
{
  int i,neq;
  struct  
  {  double value;
     int   index;
  } localmax, globalmax;

  neq=E->lmesh.NEQ[lev];

  /* Find local maximum */
  localmax.value = 0.;
  for (i=0;i<neq;i++) 
  { if (E->parallel.IDD[lev][i]) /* only get the sum from relevant ID */
      if (A[i] > localmax.value) 
      {  localmax.value = A[i];
         localmax.index = i + E->parallel.me*neq;
      }
  }

  /* Find the maximax over all processors */
  MPI_Allreduce(&localmax, &globalmax,1,MPI_DOUBLE_INT,MPI_MAXLOC,
                MPI_COMM_WORLD);
  *Amax = globalmax.value;

  return (globalmax.index);
}

int dummy_global_pmax( E, A, Amax, lev )
   struct All_variables *E;
   double *A, *Amax;
   int lev;
{
  int i,npno;
  struct  
  {  double value;
     int   index;
  } localmax, globalmax;

  *Amax = 12345;
  return (1);
}


int global_pmax( E, A, Amax, lev )
   struct All_variables *E;
   double *A, *Amax;
   int lev;
{
  int i,npno;
  struct  
  {  double value;
     int   index;
  } localmax, globalmax;

  npno=E->lmesh.NPNO[lev];

  /* Find local maximum */
  localmax.value = A[1];
  localmax.index = 1 + E->parallel.me*npno;
  for (i=2;i<=npno;i++) 
     if (A[i] > localmax.value) 
     {  localmax.value = A[i];
        localmax.index = i + E->parallel.me*npno;
     }

  /* Find the maximax over all processors */
  MPI_Allreduce(&localmax, &globalmax,1,MPI_DOUBLE_INT,MPI_MAXLOC,
                MPI_COMM_WORLD);
  *Amax = globalmax.value;

  return (globalmax.index);
}

