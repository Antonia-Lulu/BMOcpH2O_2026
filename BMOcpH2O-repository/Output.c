/* Routine to process the output of the finite element cycles 
   and to turn them into a coherent suite  files  */


#include <fcntl.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>             /* for "system" command */
#ifndef __sunos__               /* string manipulations */
#include <strings.h>
#else
#include <string.h>
#endif

#include "element_definitions.h"
#include "global_defs.h"
#include "prem_density.h"

/* Return 1 if a .BMO_comp.dat line should be written this step (rank 0 only). */
int BMO_comp_should_write(E, step)
  struct All_variables *E;
  int step;
{
  static float last_z_written = -999.0f;
  static int   last_row_written = -1;

  if (E->parallel.me != 0)
      return 0;
  if (E->control.init_age_distr != 160 && E->control.BMO_implicit_mode != 1)
      return 0;

  if (E->control.BMO_comp_per_row == 1) {
      if (last_row_written < 0
          || E->control.MO_depth_row != last_row_written) {
          last_row_written = E->control.MO_depth_row;
          return 1;
      }
      return 0;
  }

  if (E->control.BMO_comp_dz > 0.0f) {
      if (last_z_written < -900.0f
          || fabsf(E->control.MO_depth - last_z_written)
             >= E->control.BMO_comp_dz) {
          last_z_written = E->control.MO_depth;
          return 1;
      }
      return 0;
  }

  return (step % E->control.record_every) == 0;
}

void output_BMO_comp(E, file_number, write_header)
  struct All_variables *E;
  int file_number;
  int write_header;
{
  char output_file[255];

  sprintf(output_file, "%s.BMO_comp.dat", E->control.data_file2);
  FILE *fp_bmo = fopen(output_file, "a");
  if (fp_bmo == NULL)
      return;

  if (write_header) {
      if (E->control.BMO_density_mode == 1
          && E->tracers.ntracbuoyfunc > 0)
          fprintf(fp_bmo,
              "# step  time  z_BOMB  T_BOMB  St_eff  "
              "xi_l(FeO_liq)  Cw_l(H2O_liq)  Chpe_l(HPE_liq)  "
              "XFe_sol(front)  Cw_sol(front)  Chpe_sol(front)  "
              "rho_nv(z_BOMB)  rho_kg_m3(z_BOMB)\n");
      else
          fprintf(fp_bmo,
              "# step  time  z_BOMB  T_BOMB  St_eff  "
              "xi_l(FeO_liq)  Cw_l(H2O_liq)  Chpe_l(HPE_liq)  "
              "XFe_sol(front)  Cw_sol(front)  Chpe_sol(front)\n");
  }

  /* Effective Stefan number for this step */
  float St_eff = (E->control.BMO_dS_Cp > 0.0f)
               ? E->control.MO_temperature * E->control.BMO_dS_Cp
               : E->control.BMO_St;

  /* Solid composition at the row that last crystallised.
   * MO_depth_row is the melt-side front; prof[row_f] is still
   * liquid/uninitialised — use row_f-1 (方案 B).              */
  int   row_f  = E->control.MO_depth_row;
  int   row_s  = row_f - 1;
  float XFe_s  = 0.0f, Cw_s = 0.0f, Chpe_s = 0.0f;
  if (E->control.BMO_XFe_prof  != NULL && row_s >= 1
      && row_s <= E->control.BMO_prof_size)
      XFe_s = E->control.BMO_XFe_prof[row_s];
  if (E->control.BMO_CH2O_prof != NULL && row_s >= 1
      && row_s <= E->control.BMO_prof_size)
      Cw_s  = E->control.BMO_CH2O_prof[row_s];
  if (E->control.BMO_Chpe_prof != NULL && row_s >= 1
      && row_s <= E->control.BMO_prof_size)
      Chpe_s = E->control.BMO_Chpe_prof[row_s];

  if (E->control.BMO_density_mode == 1
      && E->tracers.ntracbuoyfunc > 0)
  {
      float rho_k, rho_nv_k, z_ctr;
      int   row_rho = (row_s >= 1) ? row_s : 1;
      int   lev_rho = E->mesh.levmax;
      int   elz_rho = E->lmesh.NOZ[lev_rho] - 1;

      if (row_rho > elz_rho) row_rho = elz_rho;
      z_ctr    = 0.5f * (E->XL[2][row_rho] + E->XL[2][row_rho + 1]);
      rho_nv_k = lookup_BMO_rho_nv(E, z_ctr);
      rho_k    = BMO_du2019_rho_kg_m3(rho_nv_k, XFe_s, Cw_s,
                                      E->control.BMO_n_spin,
                                      E->control.BMO_alpha_water);
      fprintf(fp_bmo,
          "%6d  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  "
          "%.6e  %.6e  %.6e  %.6e  %.6e\n",
          file_number,
          (double)E->monitor.elapsed_time,
          (double)E->control.MO_depth,
          (double)E->control.MO_temperature,
          (double)St_eff,
          (double)E->control.BMO_xi_l,
          (double)E->control.BMO_Cw_l,
          (double)E->control.BMO_Chpe_l,
          (double)XFe_s,
          (double)Cw_s,
          (double)Chpe_s,
          (double)rho_nv_k,
          (double)rho_k);
  }
  else
  {
      fprintf(fp_bmo,
          "%6d  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  "
          "%.6e  %.6e  %.6e\n",
          file_number,
          (double)E->monitor.elapsed_time,
          (double)E->control.MO_depth,
          (double)E->control.MO_temperature,
          (double)St_eff,
          (double)E->control.BMO_xi_l,
          (double)E->control.BMO_Cw_l,
          (double)E->control.BMO_Chpe_l,
          (double)XFe_s,
          (double)Cw_s,
          (double)Chpe_s);
  }

  fclose(fp_bmo);

  /* --- HPE mantle heating diagnostics (only when BMO_hpe_mantle_on=1) ---
   * for the mantle mean; the actual mean evolves as the BMO crystallises. */
  if (E->control.BMO_hpe_mantle_on) {
      char hpe_file[255];
      FILE *fp_hpe;
      float hpe_decay, cHPE_s0, H_coeff;
      sprintf(hpe_file, "%s.HPE_mantle.dat", E->control.data_file2);
      fp_hpe = fopen(hpe_file, "a");
      if (fp_hpe != NULL) {
          if (write_header)
              fprintf(fp_hpe,
                  "# step  time  H_coeff(t)  cHPE_s0  H_Q0ref(t)\n"
                  "# H_coeff = H_mantle*mult*exp(-t*ln2/tau_hpe)\n"
                  "# H_Q0ref = H_coeff * cHPE_s0  (ref mantle HPE heating rate)\n"
                  "# H_mantle=%.6g  mult=%.6g  tau_hpe=%.6g  Dhpe=%.6g  Chpe_l_melt0=%.6g\n",
                  (double)E->control.BMO_H_mantle,
                  (double)E->control.BMO_H_mantle_mult,
                  (double)E->control.BMO_tau_hpe,
                  (double)E->control.BMO_Dhpe,
                  (double)E->control.BMO_Chpe_l_melt0);
          hpe_decay  = (E->control.BMO_tau_hpe > 1.0e-10f)
                       ? expf(-(float)M_LN2 * (float)E->monitor.elapsed_time
                              / E->control.BMO_tau_hpe)
                       : 1.0f;
          cHPE_s0    = E->control.BMO_Dhpe * E->control.BMO_Chpe_l_melt0;
          H_coeff    = E->control.BMO_H_mantle * E->control.BMO_H_mantle_mult
                       * hpe_decay;
          fprintf(fp_hpe,
              "%6d  %.6e  %.6e  %.6e  %.6e\n",
              file_number,
              (double)E->monitor.elapsed_time,
              (double)H_coeff,
              (double)cHPE_s0,
              (double)(H_coeff * cHPE_s0));
          fclose(fp_hpe);
      }
  }
}

void output_velo_related(E,file_number)
  struct All_variables *E;
  int file_number; 
{
  int el,els,i,j,k,ii,m,node,fd;
  int nox,noz,noy,nfx,nfz,nfy1,nfy2,size1,size2;
  int lev, itrac, nels, iel, nodenr;
  float val;
  char output_file[255];
  float fnode[E->lmesh.nno+1];
  static float *SV,*EV;
  static int been_here=0;
  float vs;

  void get_surface_velo ();
  void edot_from_elem_to_nodes();
  void get_ele_visc ();
  const int nno = E->mesh.nno;

  lev = E->mesh.levmax;

  if (been_here==0 && E->control.restart==0) {
    sprintf(output_file,"%s.coord.%d",E->control.data_file2,E->parallel.me);
    E->filed[13]=fopen(output_file,"w");
    been_here++;
    fprintf(E->filed[13],"%6d\n",E->lmesh.nno);
    if (E->mesh.nsd==3)
       for (i=1;i<=E->lmesh.nno;i++)
	 fprintf(E->filed[13],"%.5e %.5e %.5e\n",
                 E->X[1][i],E->X[2][i],E->X[3][i]);
    else
       for (i=1;i<=E->lmesh.nno;i++)
          fprintf(E->filed[13],"%.5e %.5e\n",E->X[1][i],E->X[2][i]);
    fclose(E->filed[13]);
  }

  if (E->parallel.me<E->parallel.nprocz) {
     sprintf(output_file,"%s.ave.%d.%d",
             E->control.data_file,E->parallel.me,file_number);
     E->filed[10]=fopen(output_file,"w");

     fprintf(E->filed[10],"%6d %6d %.5e %.5e %.5e %.5e %.5e %.5e\n",
	     E->lmesh.noz,E->advection.timesteps,E->monitor.elapsed_time,
             E->slice.Nut,E->data.buoy_flux, E->monitor.Vrms,
	     E->tracers.melting_per_step,E->tracers.melting_at_trac);

     for (i=1;i<=E->lmesh.noz;i++)
        fprintf(E->filed[10],"%.4e %.4e %.4e\n",E->Have.T[i],
                E->Have.vrms[i],E->Have.Vi[i]);
     fclose(E->filed[10]);
     //E->control.MO_depth_row--;
     //if (E->control.MO_depth_row<1) E->control.MO_depth_row=1;
     //fprintf (stderr,"MO_depth %g %d\n", E->control.MO_depth_initial, E->control.MO_depth_row); 
  }


  if (file_number%(1*E->control.record_every)==0)  {
     sprintf(output_file,"%s.temp.%d.%d",E->control.data_file2,
             E->parallel.me,file_number);
     E->filed[10]=fopen(output_file,"w");
     fprintf(E->filed[10],"%6d %6d %.5e\n",E->lmesh.nno,
             E->advection.timesteps,E->monitor.elapsed_time);
     if (file_number%(1*E->control.record_every)==0)  {
        if (E->mesh.nsd==2)
        {  for (i=1;i<=E->lmesh.nno;i++)
           {  fprintf(E->filed[10],"%.5e %.5e %.5e\n",E->T[i],
                         E->V[1][i],E->V[2][i]);
	                }
        }
        else if (E->mesh.nsd==3)
        {  for (i=1;i<=E->lmesh.nno;i++)
              fprintf(E->filed[10],"%.5e %.5e %.5e %.5e\n",E->T[i],
                      E->V[1][i],E->V[2][i],E->V[3][i]);
        }
     }
     else   {
        for (i=1;i<=E->lmesh.nno;i++)
            fprintf(E->filed[10],"%.5e\n",E->T[i]);
     }
     fclose(E->filed[10]);

/* Tracer output */
     if (E->tracers.ntracfunc > 0)
     {  for (itrac=0; itrac< E->tracers.ntracfunc; itrac++)
        {  /* felem exists only when interp_t2e==1 (e.g. tracfunc 82 flag) */
           if (E->tracers.interp_t2e[itrac] != 1
               || E->tracers.felem[itrac] == NULL)
              continue;
           sprintf(output_file,"%s.trac.%d.%d.%d",E->control.data_file2,
                itrac,E->parallel.me,file_number);
           E->filed[10]=fopen(output_file,"w");
           //fprintf(stderr, "file opened\n");
           fprintf(E->filed[10],"%6d %6d %.5e\n",E->lmesh.nno,
                E->advection.timesteps,E->monitor.elapsed_time);
	   /* interp from elems to nodes */
           edot_from_elem_to_nodes(E,E->tracers.felem[itrac],
				   fnode,E->mesh.levmax);
	   for (i=1;i<=E->lmesh.nno;i++)
           { 
              fprintf(E->filed[10],"%g\n", fnode[i]);
           }
           fclose(E->filed[10]);
        }
        for (itrac=0; itrac< E->tracers.ntracfunc; itrac++)
        {  sprintf(output_file,"%s.ftrac.%d.%d.0",E->control.data_file2,
                itrac,E->parallel.me);
        }
     }

/* Shearheat output */
     if (E->control.shear_heating==2 &&
	E->monitor.solution_cycles>0)
        {  
           sprintf(output_file,"%s.shearheat.%d.%d",E->control.data_file2,
                E->parallel.me,file_number);
           E->filed[10]=fopen(output_file,"w");
           //fprintf(stderr, "file opened\n");
           fprintf(E->filed[10],"%6d %6d %.5e\n",E->lmesh.nno,
                E->advection.timesteps,E->monitor.elapsed_time);
	   /* interp from elems to nodes */
           edot_from_elem_to_nodes(E,E->visc_diss,
				   fnode,E->mesh.levmax);
	   for (i=1;i<=E->lmesh.nno;i++)
           {  
	     //if ((i%(145))<75 && (i%(145))>0) 
              fprintf(E->filed[10],"%g\n", fnode[i]);
	     //else
	     // fprintf(E->filed[10],"0.0\n", fnode[i]);
           }
           fclose(E->filed[10]);
     }

     /*
      * Disabled output for files whose name contains ".th_t.".
      * (Requested: comment out instead of deleting.)
      */
     /*
     if (E->parallel.me_loc[2]==0)  {
        sprintf(output_file,"%s.th_t.%d.%d",E->control.data_file2,
                E->parallel.me,file_number);
        E->filed[11]=fopen(output_file,"w");
        fprintf(E->filed[11],"%6d %6d %.5e %.5e %.5e\n",E->lmesh.nsf,
                E->advection.timesteps,E->monitor.elapsed_time,
                E->slice.Nut,E->data.buoy_flux);
        for (i=1;i<=E->lmesh.nsf;i++)  {
           fprintf(E->filed[11],"%.5e %.5e\n",E->slice.tpg[i],
                   E->slice.shflux[i]);
        }
        fclose(E->filed[11]);
     }
     */

  }

/* In comment: no bottom topography
  if (E->parallel.me_loc[2]==E->parallel.nprocz-1)  {
      sprintf(output_file,"%s.th_b.%d.%d",E->control.data_file2,
              E->parallel.me,file_number);
      E->filed[11]=fopen(output_file,"w");
      fprintf(E->filed[11],"%6d %6d %.5e %.5e\n",E->lmesh.nsf,
              E->advection.timesteps,E->monitor.elapsed_time,E->slice.Nub);
      for (i=1;i<=E->lmesh.nsf;i++)  {
         vs = E->V[1][i*E->lmesh.noz];
         fprintf(E->filed[11],"%.5e %.5e\n",E->slice.tpgb[i],
         E->slice.bhflux[i]);
      }
      fclose(E->filed[11]);
  }
*/

return;
}

void output_velo_related_binary(E,file_number)
  struct All_variables *E;
  int file_number; 
{
  int el,els,i,j,k,ii,m,node,fd;
  int nox,noz,noy,nfx,nfz,nfy1,nfy2,size1,size2;
  char output_file[255];
  static float *SV,*EV;
  static int been_here=0;

  void get_surface_velo ();
  void get_ele_visc ();
  const int nno = E->lmesh.nno;
/*
  if (been_here==0 && E->control.restart==0) {
    sprintf(output_file,"%s.velo",E->control.data_file);
    E->filed[10]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.topo_t",E->control.data_file);
    E->filed[11]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.topo_b",E->control.data_file);
    E->filed[12]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.visc",E->control.data_file);
    E->filed[13]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.fas670",E->control.data_file);
    E->filed[14]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.stress",E->control.data_file);
    E->filed[9]=open(output_file,O_RDWR | O_CREAT, 0644);
    }

  if (been_here==0)  {
    ii = E->mesh.nsf;
    SV = (float *) malloc ((2*ii+2)*sizeof(float));

    size2 = (E->mesh.nel+1)*sizeof(float);
    EV = (float *) malloc (size2);
    been_here++;
    }

  ii = E->mesh.nsf;
  size2 = 2*(ii+2)*sizeof(float);
    get_surface_velo (E,SV);
  write(E->filed[10],SV,size2);

  size2 = (E->mesh.nsf+1)*sizeof(float);
  write(E->filed[11],E->slice.tpg,size2);
  write(E->filed[12],E->slice.tpgb,size2);

  size2 = (E->mesh.nel+1)*sizeof(float);
    get_ele_visc (E,EV);
  write(E->filed[13],EV,size2);

  size2 = (E->mesh.nsf+1)*sizeof(float);
  write(E->filed[14],E->Fas670_b,size2);

  size2 = (2*E->mesh.nsf+1)*sizeof(float);
  write(E->filed[9],E->stress,size2);
*/
  return;
  }

/* ====================================================================== */

void output_temp(E,file_number)
  struct All_variables *E;
  int file_number; 
{
  int nno,i,j,fd;
  static int *temp1;
  static int been_here=0;
  static int size2,size1;
  char output_file[255];
/*
  if (been_here==0 && E->control.restart==0) {
    sprintf(output_file,"%s.temp",E->control.data_file);
    E->filed[5]=open(output_file,O_RDWR | O_CREAT, 0644);
    }

  if (been_here==0) {
    temp1 = (int *) malloc ((E->lmesh.noy*6)*sizeof(int));

    sprintf(output_file,"%s.mesh",E->control.data_file);
    E->filed[1]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.x",E->control.data_file);
    E->filed[2]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.z",E->control.data_file);
    E->filed[3]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.y",E->control.data_file);
    E->filed[4]=open(output_file,O_RDWR | O_CREAT, 0644);

    size1 = (E->lmesh.noy*6)*sizeof(int);
    size2= (E->lmesh.nno+1)*sizeof(float);

    temp1[1] = E->lmesh.nno;
    temp1[3] = size2;
    temp1[5] = E->lmesh.nsf;
    temp1[6] = E->lmesh.nel;

        write(E->filed[1],temp1,size1);
        write(E->filed[2],E->X[1],size2);
        write(E->filed[3],E->X[2],size2);
        write(E->filed[4],E->X[3],size2);

    close(E->filed[1]);
    close(E->filed[2]);
    close(E->filed[3]);
    close(E->filed[4]);

    been_here++;
    }


    write(E->filed[5],E->T,size2);
*/
  return; 
}


/* ====================================================================== */

void process_restart(E)
    struct All_variables *E;
{
    int fileid[20];
    int i,j,k,ii,size2;
    char output_file[255],in_file[255];
/*
    sprintf(output_file,"%s.temp",E->control.data_file);
    E->filed[5]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.velo",E->control.data_file);
    E->filed[10]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.topo_t",E->control.data_file);
    E->filed[11]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.topo_b",E->control.data_file);
    E->filed[12]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.visc",E->control.data_file);
    E->filed[13]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.fas670",E->control.data_file);
    E->filed[14]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.stress",E->control.data_file);
    E->filed[9]=open(output_file,O_RDWR | O_CREAT, 0644);

    sprintf(in_file,"%s.temp",E->control.data_file1);
    fileid[5]=open(in_file,O_RDONLY,0);
    sprintf(in_file,"%s.velo",E->control.data_file1);
    fileid[10]=open(in_file,O_RDONLY,0);
    sprintf(in_file,"%s.topo_t",E->control.data_file1);
    fileid[11]=open(in_file,O_RDONLY,0);
    sprintf(in_file,"%s.topo_b",E->control.data_file1);
    fileid[12]=open(in_file,O_RDONLY,0);
    sprintf(in_file,"%s.visc",E->control.data_file1);
    fileid[13]=open(in_file,O_RDONLY,0);
    sprintf(in_file,"%s.fas670",E->control.data_file1);
    fileid[14]=open(in_file,O_RDONLY,0);
    sprintf(in_file,"%s.stress",E->control.data_file1);
    fileid[9]=open(in_file,O_RDONLY,0);


  for (i=1;i<=E->control.restart_frame;i++)  {
    ii = E->mesh.nsf;
    size2 = (2*ii+2)*sizeof(float);
    read(fileid[10],E->T,size2);
    write(E->filed[10],E->T,size2);

    size2 = (E->mesh.nsf+1)*sizeof(float);
    read(fileid[11],E->T,size2);
    write(E->filed[11],E->T,size2);
    read(fileid[12],E->T,size2);
    write(E->filed[12],E->T,size2);

    size2 = (E->mesh.nel+1)*sizeof(float);
    read(fileid[13],E->T,size2);
    write(E->filed[13],E->T,size2);

    size2 = (E->mesh.nsf+1)*sizeof(float);
    read(fileid[14],E->T,size2);
    write(E->filed[14],E->T,size2);

    size2 = (12*E->mesh.nsf+1)*sizeof(float);
    read(fileid[9],E->T,size2);
    write(E->filed[9],E->T,size2);

    size2= (E->mesh.nno+1)*sizeof(float);
    read(fileid[5],E->T,size2);
    write(E->filed[5],E->T,size2);

    }

    close(fileid[5]);
    close(fileid[10]);
    close(fileid[11]);
    close(fileid[12]);
    close(fileid[13]);
    close(fileid[14]);
    close(fileid[9]);
*/
  return;
  }

void time_display(E, initial_time, step)
     /*tad 12-14-05 displays time in stderr, updates initial
 * time*/
     struct All_variables *E;
     double *initial_time;
     int step;

{
  double CPU_time0(),time;

  if (E->control.show_time && E->parallel.me == 0){
       time = CPU_time0() - *initial_time;
       fprintf(stderr, "time to %d  = %f sec\n", step, time);
       *initial_time = CPU_time0();
     }

}

