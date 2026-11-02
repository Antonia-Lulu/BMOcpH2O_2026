/* Routines here are for intel paragon with MPI */

#if 1
#include <mpi.h>
#endif

#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"


void parallel_process_initilization(E,argc,argv)
  struct All_variables *E;
  int argc;
  char **argv;
  {

  E->parallel.me = 0;
  E->parallel.nproc = 1;
  E->parallel.me_loc[1] = 0;
  E->parallel.me_loc[2] = 0;
  E->parallel.me_loc[3] = 0;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &(E->parallel.me) );
  MPI_Comm_size(MPI_COMM_WORLD, &(E->parallel.nproc) );

  return;
  }


/* ============================================ */
/* ============================================ */

void parallel_process_termination()
{

  MPI_Finalize();
  exit(8);
  return;
  }


/* ============================================ */
/* ============================================ */


void parallel_domain_decomp1(E)
  struct All_variables *E;
  {
  
  int i,j,k,nox,noz,noy,me;

  me = E->parallel.me;

  if (E->parallel.automa)  {
    }

 if (E->control.NMULTIGRID||E->control.EMULTIGRID)
    {
      E->lmesh.elx = E->lmesh.mgunitx * 
                     (int) pow(2.0,((double)E->mesh.levmax));
      E->lmesh.elz = E->lmesh.mgunitz *
                     (int) pow(2.0,((double)E->mesh.levmax));
      if(E->mesh.nsd == 3 )
          E->lmesh.ely = E->lmesh.mgunity * 
                         (int) pow(2.0,((double)E->mesh.levmax));

    }

  E->lmesh.elx = E->mesh.elx / E->parallel.nprocx;
  E->lmesh.elz = E->mesh.elz / E->parallel.nprocz;
  E->lmesh.ely = E->mesh.ely / E->parallel.nprocy;

  E->parallel.nprocxz = E->parallel.nprocx*E->parallel.nprocz;
  E->parallel.nprocxy = E->parallel.nprocx*E->parallel.nprocy;
  E->parallel.nproczy = E->parallel.nprocz*E->parallel.nprocy;

  /* for overlapping domain, good for e by e assemble */

        /* z direction first */
  j = me % E->parallel.nprocz;
  E->parallel.me_loc[2] = j;
  E->lmesh.nzs = j*E->lmesh.elz+1; /* 1st node for proc me in z-dir */

        /* y direction then */
  k = (me+1)/E->parallel.nprocxz-(((me+1)%E->parallel.nprocxz==0)?1:0);
  E->parallel.me_loc[3] = k;
  E->lmesh.nys = k*E->lmesh.ely+1; /* 1st node for proc me in y-dir */

        /* x direction then */
  i = (me+1-k*E->parallel.nprocxz)/E->parallel.nprocz - ( ((me+1-k*E->parallel.nprocxz)%E->parallel.nprocz==0)? 1:0 );
  E->parallel.me_loc[1] = i;
  E->lmesh.nxs = i*E->lmesh.elx+1; /* 1st node for proc me in x-dir */

/* fprintf(stderr,"b %d %d %d %d %d %d
%d\n",E->parallel.me,E->parallel.me_loc[1],E->parallel.me_loc[2],E->parallel.me_loc[3],E->lmesh.nxs,E->lmesh.nzs,E->lmesh.nys); */

  E->lmesh.nox = E->lmesh.nnx[1] = E->lmesh.elx + 1;
  E->lmesh.noz = E->lmesh.nnx[2] = E->lmesh.elz + 1;
  if (E->mesh.nsd == 2)
     E->lmesh.noy = E->lmesh.nnx[3] = 1;
  else if (E->mesh.nsd == 3)
     E->lmesh.noy = E->lmesh.nnx[3] = E->lmesh.ely + 1;

  E->lmesh.nno = E->lmesh.noz*E->lmesh.nox*E->lmesh.noy;
  E->lmesh.nsf = E->lmesh.nox*E->lmesh.noy;
  E->lmesh.npno = E->lmesh.elz*E->lmesh.elx*E->lmesh.ely;
  E->lmesh.nel = E->lmesh.elz*E->lmesh.elx*E->lmesh.ely;
  E->lmesh.nex[1] = E->lmesh.elx;
  E->lmesh.nex[2] = E->lmesh.elz;
  E->lmesh.nex[3] = E->lmesh.ely;
  E->lmesh.nnov = E->lmesh.nno;
  E->lmesh.neq = E->lmesh.nnov*E->mesh.nsd;

  for(i=E->mesh.levmax;i>=E->mesh.levmin;i--)   {
     if (E->control.NMULTIGRID||E->control.EMULTIGRID)  {
	nox = E->lmesh.elx / ((int) pow(2.0,(double)(E->mesh.levmax-i))) + 1;
	noz = E->lmesh.elz / ((int) pow(2.0,(double)(E->mesh.levmax-i))) + 1;
	if(E->mesh.nsd==3)
	  noy = E->lmesh.ely / ((int) pow(2.0,(double)(E->mesh.levmax-i))) + 1;
        else
	  noy = 1;
        }
     else      {
	nox = E->lmesh.nox;
	noz = E->lmesh.noz;
	noy = E->lmesh.noy;
	}

     E->lmesh.ELX[i] = nox-1;
     E->lmesh.ELZ[i] = noz-1;
     E->lmesh.ELY[i] = max(noy-1,1);
     E->lmesh.NNO[i] = nox * noz * noy;
     E->lmesh.NEL[i] = (nox-1) * (noz-1) * max((noy-1),1);
     E->lmesh.NPNO[i] = E->lmesh.NEL[i] ;
     E->lmesh.NOX[i] = nox;
     E->lmesh.NOZ[i] = noz;
     E->lmesh.NOY[i] = noy;
     E->lmesh.NNOV[i] = E->lmesh.NNO[i];
     E->lmesh.NEQ[i] = E->mesh.nsd * E->lmesh.NNOV[i] ;
     E->lmesh.NXS[i] = E->parallel.me_loc[1]*E->lmesh.ELX[i]+1;
     E->lmesh.NZS[i] = E->parallel.me_loc[2]*E->lmesh.ELZ[i]+1;
     E->lmesh.NYS[i] = E->parallel.me_loc[3]*E->lmesh.ELY[i]+1;
/*fprintf(E->fp,"b %d %d %d %d %d %d %d\n",E->parallel.me,E->lmesh.ELX[i],E->lmesh.ELZ[i],E->lmesh.ELY[i],E->lmesh.NNO[i],E->lmesh.NEL[i],E->lmesh.NEQ[i]); */

     }

  return;
  }

/* ============================================ 
 shuffle element order and determine boundary nodes for 
 exchange info across the boundaries
 ============================================ */

void parallel_shuffle_ele_and_id(E)
  struct All_variables *E;
  {

  int i,ii,j,k,l,node,node1,el,elt,lnode,llnode,jj,k1,k2;
  int lev,elx,elz,ely,nel,nno,nox,noz,noy;

  for(lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--){
      nel = E->lmesh.NEL[lev];
      elx = E->lmesh.ELX[lev];
      elz = E->lmesh.ELZ[lev];
      ely = E->lmesh.ELY[lev];
      nox = E->lmesh.NOX[lev];
      noy = E->lmesh.NOY[lev];
      noz = E->lmesh.NOZ[lev];
 
      ii = 0;

      for(i=1;i<=2;i++)       {       /* do the ZOY boundary elements first */
         ii ++;
         if ( (i==1&&E->parallel.me_loc[1]!=0) || (i==2&&E->parallel.me_loc[1]!=E->parallel.nprocx-1))        {
   	     for(k=1;k<=noy;k++)
 	       for(j=1;j<=noz;j++)   {
	         node = j + ( ((i==1)?1:nox)-1 )*noz + (k-1)*noz*nox;
	         llnode = j+(k-1)*noz;
	         E->parallel.NODE[lev][llnode].bound[ii] = node;
                 E->NODE[lev][node] = E->NODE[lev][node] | OFFSIDE;
                 E->NODE[lev][node] = E->NODE[lev][node] | LIDN;
	         node1 = node + ( ((i==1)?1:-1) )*noz;
                 E->NODE[lev][node1] = E->NODE[lev][node1] | LIDN;
	         }
	     E->parallel.NUM_NNO[lev].bound[ii] = noy*noz;
             }
	 }         /* end for i   */

      for(j=1;j<=2;j++)       {   /* do XOY boundary elements */
         ii ++;
         if ( (j==1&&E->parallel.me_loc[2]!=0) || (j==2&&E->parallel.me_loc[2]!=E->parallel.nprocz-1))        {
  	     for(k=1;k<=noy;k++)
 	       for(i=1;i<=nox;i++)   {
	         node = ((j==1)?1:noz) + (i-1)*noz + (k-1)*noz*nox;
	         llnode = i+(k-1)*nox;
	         E->parallel.NODE[lev][llnode].bound[ii] = node;
                 E->NODE[lev][node] = E->NODE[lev][node] | OFFSIDE;
                 E->NODE[lev][node] = E->NODE[lev][node] | LIDN;
	         node1 = node + ( ((j==1)?1:-1) );
                 E->NODE[lev][node1] = E->NODE[lev][node1] | LIDN;
	         }
	     E->parallel.NUM_NNO[lev].bound[ii] = nox*noy;
	     }
	  }         /* end for j   */

       if (E->mesh.nsd==3)  {
 	 for(k=1;k<=2;k++)        {  /* do XOZ boundary elements for 3D */
           ii ++;
           if ( (k==1&&E->parallel.me_loc[3]!=0) || (k==2&&E->parallel.me_loc[3]!=E->parallel.nprocy-1))        {
 	       for(j=1;j<=noz;j++)
 	         for(i=1;i<=nox;i++)   {
	           node = j + (i-1)*noz + (((k==1)?1:noy)-1)*noz*nox;
	           llnode = j+(i-1)*noz;
	           E->parallel.NODE[lev][llnode].bound[ii] = node;
                   E->NODE[lev][node] = E->NODE[lev][node] | OFFSIDE;
                   E->NODE[lev][node] = E->NODE[lev][node] | LIDN;
	           node1 = node + ( ((k==1)?1:-1) )*noz*nox;
                   E->NODE[lev][node1] = E->NODE[lev][node1] | LIDN;
	           }
	       E->parallel.NUM_NNO[lev].bound[ii] = nox*noz;
	       }       /* end for k */
	     }       /* end for k */
	   }       /* end for dims=3  */


      E->parallel.num_b = ii;

      }   /* end for level */
  
  return;
  }


/* ============================================ 
 determine communication routs  and boundary ID for 
 exchange info across the boundaries
 ============================================ */

void parallel_communication_routs(E)
  struct All_variables *E;
  {

  int i,ii,j,k,l,node,el,elt,lnode,jj,doff;
  int lev,elx,elz,ely,nno,nox,noz,noy,p,kkk,kk;
  int me, nprocz,nprocx, nprocy;

  me = E->parallel.me;
  nprocx = E->parallel.nprocx;
  nprocz = E->parallel.nprocz;
  nprocy = E->parallel.nprocy;
  
  for(lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--){
     elx = E->lmesh.ELX[lev];
     elz = E->lmesh.ELZ[lev];
     ely = E->lmesh.ELY[lev];
     nox = E->lmesh.NOX[lev];
     noy = E->lmesh.NOY[lev];
     noz = E->lmesh.NOZ[lev];
     nno = E->lmesh.NNO[lev];
 
     ii = 0;
     kkk = 0;
     if (E->mesh.nsd == 2) {     

 	   for(i=1;i<=2;i++)    {     /* do the OZ boundaries */

         ii ++;

	     E->parallel.NUM_PASS[lev].bound[ii] = 1;
	     if (E->parallel.me_loc[1]==0 && i==1)
	       E->parallel.NUM_PASS[lev].bound[ii] = 0;
	     else if (E->parallel.me_loc[1]==nprocx-1 && i==2)
	       E->parallel.NUM_PASS[lev].bound[ii] = 0;

         for (p=1;p<=E->parallel.NUM_PASS[lev].bound[ii];p++)  {

	       kkk ++;
                 /* determine the pass ID for ii-th boundary and p-th pass */

                 /* for kkk th pass, get the # of data to be passed */
	       E->parallel.NUM_NEQ[lev].pass[kkk] = ((p==1)?E->parallel.NUM_NNO[lev].bound[ii]:noy) *E->mesh.nsd;
	       E->parallel.NUM_NODE[lev].pass[kkk] = ((p==1)?E->parallel.NUM_NNO[lev].bound[ii]:noy);

  /* for kkk th pass, determine the target processor */

	       E->parallel.PROCESSOR[lev].pass[kkk]=me-((i==1)?1:-1)*nprocz;

  /* for kkk th pass, determine the ID of equations to be passed */

		   for (node=1;node<=E->parallel.NUM_NODE[lev].pass[kkk];node++)   {
		      lnode = E->parallel.NODE[lev][node].bound[ii]; 
		      E->parallel.EXCHANGE_NODE[lev][node].pass[kkk] = lnode; 
		      for(doff=1;doff<=E->mesh.nsd;doff++)     {
		         jj = doff + (node-1)*E->mesh.nsd;
	 	         E->parallel.EXCHANGE_ID[lev][jj].pass[kkk] = E->ID[lev][lnode].doff[doff];
		         }
		      }

	       }     /* end for loop p */
         }       /* end for loop i */

       E->parallel.IDPASS[lev][1].bound[1] = kkk;       /* pass in x direction */

       for(j=1;j<=2;j++)       {    /* do OX boundary */

         ii ++;

	     E->parallel.NUM_PASS[lev].bound[ii] = 1;
	     if (E->parallel.me_loc[2]==0 && j==1)
	       E->parallel.NUM_PASS[lev].bound[ii] = 0;
	     else if (E->parallel.me_loc[2]==nprocz-1 && j==2)
	       E->parallel.NUM_PASS[lev].bound[ii] = 0;

         for (p=1;p<=E->parallel.NUM_PASS[lev].bound[ii];p++)  {
	       kkk ++;

                 /* for kkk th pass, get the # of data to be passed */
	       E->parallel.NUM_NEQ[lev].pass[kkk] = E->parallel.NUM_NNO[lev].bound[ii]*E->mesh.nsd;
	       E->parallel.NUM_NODE[lev].pass[kkk] = E->parallel.NUM_NNO[lev].bound[ii];

  /* for kkk th pass, determine the target processor */

	       E->parallel.PROCESSOR[lev].pass[kkk]=me-((j==1)?1:-1);

  /* for kkk th pass, determine the ID of equations to be passed */

	       for (node=1;node<=E->parallel.NUM_NODE[lev].pass[kkk];node++)   {
		     lnode = E->parallel.NODE[lev][node].bound[ii]; 
		     E->parallel.EXCHANGE_NODE[lev][node].pass[kkk] = lnode; 
		     for(doff=1;doff<=E->mesh.nsd;doff++)  {
		       jj = doff + (node-1)*E->mesh.nsd;
	 	       E->parallel.EXCHANGE_ID[lev][jj].pass[kkk] = E->ID[lev][lnode].doff[doff];
		       }
		     }      /* end loop for node */

	       }     /* end for loop p */
         }    /* end for loop j */

       E->parallel.IDPASS[lev][1].bound[2] = kkk;       /* pass in z direction */

	 }         /* end for 2D */



    else if (E->mesh.nsd ==3)  {

 	  for(i=1;i<=2;i++)       {       /* do YOZ boundaries & OY lines */

        ii ++;
        E->parallel.NUM_PASS[lev].bound[ii] = 1;
        if(E->parallel.me_loc[1]==0 && i==1)
          E->parallel.NUM_PASS[lev].bound[ii] = 0;
        else if(E->parallel.me_loc[1]==nprocx-1 && i==2)
          E->parallel.NUM_PASS[lev].bound[ii] = 0;

        for (p=1;p<=E->parallel.NUM_PASS[lev].bound[ii];p++)  {
          kkk ++;
              /* determine the pass ID for ii-th boundary and p-th pass */
          E->parallel.IDPASS[lev][p].bound[ii] = kkk;
          E->parallel.NUM_NEQ[lev].pass[kkk] = ((p==1)?E->parallel.NUM_NNO[lev].bound[ii]:noy) *E->mesh.nsd;
          E->parallel.NUM_NODE[lev].pass[kkk] = ((p==1)?E->parallel.NUM_NNO[lev].bound[ii]:noy);

          E->parallel.PROCESSOR[lev].pass[kkk]=me-((i==1)?1:-1)*nprocz;

          for (k=1;k<=E->parallel.NUM_NODE[lev].pass[kkk];k++)   {
            lnode = k;
            node = E->parallel.NODE[lev][lnode].bound[ii];
            E->parallel.EXCHANGE_NODE[lev][k].pass[kkk] = node;
            for(doff=1;doff<=E->mesh.nsd;doff++)     {
              jj = doff + (k-1)*E->mesh.nsd;
              E->parallel.EXCHANGE_ID[lev][jj].pass[kkk] = E->ID[lev][node].doff[doff];
              }
            }  /* end for node k */

          }   /* end for loop p */
	    }  /* end for i */

 	  for(j=1;j<=2;j++)       {       /* do XOY boundaries & OX lines */

        ii ++;
        E->parallel.NUM_PASS[lev].bound[ii] = 1;
        if(E->parallel.me_loc[2]==0 && j==1)
          E->parallel.NUM_PASS[lev].bound[ii] = 0;
        else if(E->parallel.me_loc[2]==nprocz-1 && j==2)
          E->parallel.NUM_PASS[lev].bound[ii] = 0;

        for (p=1;p<=E->parallel.NUM_PASS[lev].bound[ii];p++)  {
          kkk ++;
              /* determine the pass ID for ii-th boundary and p-th pass */
          E->parallel.IDPASS[lev][p].bound[ii] = kkk;
          E->parallel.NUM_NEQ[lev].pass[kkk] = ((p==1)?E->parallel.NUM_NNO[lev].bound[ii]:nox) *E->mesh.nsd;
          E->parallel.NUM_NODE[lev].pass[kkk] = ((p==1)?E->parallel.NUM_NNO[lev].bound[ii]:nox);

          E->parallel.PROCESSOR[lev].pass[kkk]=me-((j==1)?1:-1);

          for (k=1;k<=E->parallel.NUM_NODE[lev].pass[kkk];k++)   {
            lnode = k;
            node = E->parallel.NODE[lev][lnode].bound[ii];
            E->parallel.EXCHANGE_NODE[lev][k].pass[kkk] = node;
            for(doff=1;doff<=E->mesh.nsd;doff++)     {
              jj = doff + (k-1)*E->mesh.nsd;
              E->parallel.EXCHANGE_ID[lev][jj].pass[kkk] = E->ID[lev][node].doff[doff];
              }
            }  /* end for node k */

          }   /* end for loop p */

	    }     /* end for j */

 	  for(k=1;k<=2;k++)        {      /* do XOZ boundaries & OZ lines */

        ii ++;
        E->parallel.NUM_PASS[lev].bound[ii] = 1;
        if(E->parallel.me_loc[3]==0 && k==1)
          E->parallel.NUM_PASS[lev].bound[ii] = 0;
        else if(E->parallel.me_loc[3]==nprocy-1 && k==2)
          E->parallel.NUM_PASS[lev].bound[ii] = 0;

        for (p=1;p<=E->parallel.NUM_PASS[lev].bound[ii];p++)  {

          kkk ++;
              /* determine the pass ID for ii-th boundary and p-th pass */
          E->parallel.IDPASS[lev][p].bound[ii] = kkk;
          E->parallel.NUM_NEQ[lev].pass[kkk] = ((p==1)?E->parallel.NUM_NNO[lev].bound[ii]:noz) *E->mesh.nsd;
          E->parallel.NUM_NODE[lev].pass[kkk] = ((p==1)?E->parallel.NUM_NNO[lev].bound[ii]:noz);

          E->parallel.PROCESSOR[lev].pass[kkk]=me-((k==1)?1:-1)*nprocx*nprocz;

          for (kk=1;kk<=E->parallel.NUM_NODE[lev].pass[kkk];kk++)   {
            lnode = kk;
            node = E->parallel.NODE[lev][lnode].bound[ii];
            E->parallel.EXCHANGE_NODE[lev][kk].pass[kkk] = node;
            for(doff=1;doff<=E->mesh.nsd;doff++)     {
              jj = doff + (kk-1)*E->mesh.nsd;
              E->parallel.EXCHANGE_ID[lev][jj].pass[kkk] = E->ID[lev][node].doff[doff];
              }
            }  /* end for node kk */

          }   /* end for loop p */

	    }  /* end for k */

	  }       /* end of dims==3 */

      E->parallel.TNUM_PASS[lev] = kkk;

/*
*/
    fprintf(E->fp," me= %d  pass  %d \n",E->parallel.me,E->parallel.TNUM_PASS[lev]);
    for (k=1;k<=E->parallel.TNUM_PASS[lev];k++)  
fprintf(E->fp,"proc %d and pass  %d to proc %d with %d eqn\n",E->parallel.me,k,E->parallel.PROCESSOR[lev].pass[k],E->parallel.NUM_NEQ[lev].pass[k]);
         

    }        /* end for level */
  
  return;
  }


/* ============================================ 
  exchange ID related information across boundaries
 ============================================ */

void exchange_id_d20(E, U, lev)
 struct All_variables *E;
 double *U;
 int lev;
 {
 void  parallel_process_sync(); 
 int e,node,i,ii,j,k,bound,type,idb,msginfo[8];
 static double *S[27],*R[27];

 static int been_here = 0;
 static int mid_recv, sizeofk[27];
 const int levmax = E->mesh.levmax;

 MPI_Status status;

 if (been_here ==0 )   {
    for (k=1;k<=E->parallel.TNUM_PASS[levmax];k++)  {
         sizeofk[k] = (1+E->parallel.NUM_NEQ[levmax].pass[k])*sizeof(double);
         S[k]=(double *)malloc( sizeofk[k] );
         R[k]=(double *)malloc( sizeofk[k] );
	 }
    been_here ++;
    }

    for (k=1;k<=E->parallel.TNUM_PASS[lev];k++)  {
        for (j=1;j<=E->parallel.NUM_NEQ[lev].pass[k];j++)
          S[k][j-1] = U[ E->parallel.EXCHANGE_ID[lev][j].pass[k] ];

        MPI_Sendrecv(S[k],E->parallel.NUM_NEQ[lev].pass[k],MPI_DOUBLE,
             E->parallel.PROCESSOR[lev].pass[k],1,
                     R[k],E->parallel.NUM_NEQ[lev].pass[k],MPI_DOUBLE,
             E->parallel.PROCESSOR[lev].pass[k],1,MPI_COMM_WORLD,&status);

/*        parallel_process_sync();
*/
        for (j=1;j<=E->parallel.NUM_NEQ[lev].pass[k];j++)
          U[ E->parallel.EXCHANGE_ID[lev][j].pass[k] ] += R[k][j-1];
        }

 return;
 }


void exchange_node_f20(E, U, lev)
 struct All_variables *E;
 float *U;
 int lev;
 {
 int e,node,i,ii,j,k,bound,type,idb,msginfo[8];
 static float *S[27],*R[27];

 static int been_here = 0;
 static int mid_recv, sizeofk[27];
 const int levmax = E->mesh.levmax;

 MPI_Status status;

 if (been_here ==0 )   {
    for (k=1;k<=E->parallel.TNUM_PASS[levmax];k++)  {
         sizeofk[k] = (1+E->parallel.NUM_NODE[levmax].pass[k])*sizeof(float);
         S[k]=(float *)malloc( sizeofk[k] );
         R[k]=(float *)malloc( sizeofk[k] );
	 }
    been_here ++;
    }

    for (k=1;k<=E->parallel.TNUM_PASS[lev];k++)  {
        for (j=1;j<=E->parallel.NUM_NODE[lev].pass[k];j++)
          S[k][j-1] = U[ E->parallel.EXCHANGE_NODE[lev][j].pass[k] ];

        MPI_Sendrecv(S[k],E->parallel.NUM_NODE[lev].pass[k],MPI_FLOAT,
             E->parallel.PROCESSOR[lev].pass[k],1,
                     R[k],E->parallel.NUM_NODE[lev].pass[k],MPI_FLOAT,
             E->parallel.PROCESSOR[lev].pass[k],1,MPI_COMM_WORLD,&status);

        for (j=1;j<=E->parallel.NUM_NODE[lev].pass[k];j++)
          U[ E->parallel.EXCHANGE_NODE[lev][j].pass[k] ] += R[k][j-1];
        }


 return;
}

void average_node_f20(E, U, lev)
 struct All_variables *E;
 float *U;
 int lev;
 {
 int e,node,i,ii,j,k,bound,type,idb,msginfo[8];
 static float *S[27],*R[27];

 static int been_here = 0;
 static int mid_recv, sizeofk[27];
 const int levmax = E->mesh.levmax;

 MPI_Status status;

 if (been_here ==0 )   {
    for (k=1;k<=E->parallel.TNUM_PASS[levmax];k++)  {
         sizeofk[k] = 
             (1+E->parallel.NUM_NODE[levmax].pass[k])*sizeof(float);
         S[k]=(float *)malloc( sizeofk[k] );
         R[k]=(float *)malloc( sizeofk[k] );
	 }
    been_here ++;
    }

    for (k=1;k<=E->parallel.TNUM_PASS[lev];k++)  
    {   for (j=1;j<=E->parallel.NUM_NODE[lev].pass[k];j++)
           S[k][j-1] = U[ E->parallel.EXCHANGE_NODE[lev][j].pass[k] ];

        MPI_Sendrecv(S[k],E->parallel.NUM_NODE[lev].pass[k],MPI_FLOAT,
           E->parallel.PROCESSOR[lev].pass[k],1,
                  R[k],E->parallel.NUM_NODE[lev].pass[k],MPI_FLOAT,
           E->parallel.PROCESSOR[lev].pass[k],1,MPI_COMM_WORLD,&status);

        for (j=1;j<=E->parallel.NUM_NODE[lev].pass[k];j++)
        { U[ E->parallel.EXCHANGE_NODE[lev][j].pass[k] ] += R[k][j-1];
          U[ E->parallel.EXCHANGE_NODE[lev][j].pass[k] ] /= 2;
        }
    }


 return;
}


/* ==========================   */

 double CPU_time0()
{

 double time, MPI_Wtime();

 time = MPI_Wtime();

 return (time);
}


void parallel_process_sync()
{

  MPI_Barrier(MPI_COMM_WORLD);
  return;
  }

