/* XI/2005 MDB  Updating of tracer depletion values
 *              Handling of melting
 *
 * Tracers_upd.c contains the following functions:
 *   void updatetrac()        - main updating routine
 *   float katz_melting()     - F(T,z) corresponding to Katz et al. 2003
 */ 

#include <signal.h>
#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include <string.h>
#include <stdlib.h>
#include "element_definitions.h"
#include "global_defs.h"


void updatetrac(E, PC_iteration)
   struct All_variables *E; 
   int  PC_iteration;
{
  int ifunc,itrac,ielem, c, columns, r;
  float Ttrac, dFdt, xoffset, yoffset;
  const float factor = E->data.layer_km*E->data.density*E->data.grav_acc*1e-6;
  const float factor2= E->data.heat_of_melt/(E->data.Cp*E->data.ref_temperature);
  const float HperCpA=1.+.002578947*E->data.heat_of_melt/E->data.Cp;
 /* factor2 scales dimensional heat_of_melt (400 kJ/kg=default) nondimensionally as dQ'/dt */
  float func_to_trac();
  int av_visit;
  float H2O;
  float trapmelt[E->lmesh.nel+1];
  int trapmelt_per_step, trapmelt_add;
  float xi[4],MR;
  int NEL_in_a_column,write_dFdt=0;
  int column,lev=E->mesh.levmax;
  char output_file[100];
  const float TZ_thickness=660./E->data.layer_km;
  const float shallower_km250=250./E->data.layer_km;
  float pressure, y,null=0.;
  int p,T,ndim,idim;
  FILE *fp;
  int melting;

  float hirschmann_melting();
  float katz_melting();
  float trac_funcinit();
  float calc_H2O();


  NEL_in_a_column = E->lmesh.ELZ[lev];
  columns = E->tracers.melt_lcolumns;

  /* Locate the optional per-tracer "solidified" flag channel (tracfunc 82). */
  int ifunc_solidflag = -1;
  { int jf;
    for (jf = 0; jf < E->tracers.ntracfunc; jf++)
      if (E->tracers.tracfunc[jf] == 82) ifunc_solidflag = jf;
  }

  for (ifunc=0;ifunc<E->tracers.ntracfunc;ifunc++)
  {
  /* BMO composition tracers (tracfunc 80 = FeO, 81 = H2O, 83 = HPE): "stamp once, then advect".
   *   1. Inside the melt (z >= z_BOMB): set to the liquid composition (xi_l / Cw_l / Chpe_l).
   *   2. Just crossed the front: stamp the solid composition once from BMO_*_prof
   *      (via the flag in tracfunc 82, or only within a thin shell below the front without it).
   *   3. Already solid: untouched (pure advection).
   * Once fully crystallised (z_BOMB -> 1) all resets stop.
   * Use interp_e2t=0 and interp_t2e=1 for these tracfunc codes. */
  if (E->tracers.tracfunc[ifunc] == 80 && E->control.BMO_prof_size > 0)
    {
      float z_bomb = E->control.MO_depth;
      /* Step 2: stop all resets once the BMO is fully crystallised. */
      if (z_bomb < 1.0f - 1.0e-6f)
      {
        float xi_l   = E->control.BMO_xi_l;
        int   noz_p  = E->control.BMO_prof_size;
        int   lev_p  = E->mesh.levmax;
        int   noz_m  = E->lmesh.NOZ[lev_p];
        int   flagf  = ifunc_solidflag;
        /* Thin-shell thickness near the front (fallback when no flag channel). */
        float delta_shell;
        { int rf = E->control.MO_depth_row;
          if (rf < 2)     rf = 2;
          if (rf > noz_m) rf = noz_m;
          delta_shell = 3.0f * (E->X[2][rf] - E->X[2][rf-1]);
        }

        for (itrac = 0; itrac < E->tracers.lntrac; itrac++)
          {
            /* Already-solidified tracers: never rewrite -> pure advection. */
            if (flagf >= 0 && E->tracers.ftrac[flagf][itrac] > 0.5f)
                continue;

            float zt = E->tracers.xtrac[2][itrac];
            if (zt >= z_bomb - 1.0e-5f)
              {
                /* Tracer is inside the BMO melt: liquid FeO composition. */
                E->tracers.ftrac[ifunc][itrac] = xi_l;
              }
            else if (E->control.BMO_XFe_prof != NULL)
              {
                /* Below the front: stamp the recorded solid composition (once per
                 * tracer with the flag; otherwise only within the thin shell). */
                int do_stamp = (flagf >= 0) || (zt >= z_bomb - delta_shell);
                if (do_stamp)
                  {
                    int row = 1, kr;
                    for (kr = 1; kr <= noz_m; kr++) {
                        if (E->X[2][kr] <= zt) row = kr;
                        else break;
                    }
                    if (row < 1)     row = 1;
                    if (row > noz_p) row = noz_p;
                    E->tracers.ftrac[ifunc][itrac] = E->control.BMO_XFe_prof[row];
                  }
                /* else: deep solid, no flag -> leave value (advected). */
              }
          }
      }
    }
  else if (E->tracers.tracfunc[ifunc] == 81 && E->control.BMO_prof_size > 0)
    {
      float z_bomb = E->control.MO_depth;
      if (z_bomb < 1.0f - 1.0e-6f)
      {
        float Cw_l   = E->control.BMO_Cw_l;
        int   noz_p  = E->control.BMO_prof_size;
        int   lev_p  = E->mesh.levmax;
        int   noz_m  = E->lmesh.NOZ[lev_p];
        int   flagf  = ifunc_solidflag;
        float delta_shell;
        { int rf = E->control.MO_depth_row;
          if (rf < 2)     rf = 2;
          if (rf > noz_m) rf = noz_m;
          delta_shell = 3.0f * (E->X[2][rf] - E->X[2][rf-1]);
        }

        for (itrac = 0; itrac < E->tracers.lntrac; itrac++)
          {
            if (flagf >= 0 && E->tracers.ftrac[flagf][itrac] > 0.5f)
                continue;

            float zt = E->tracers.xtrac[2][itrac];
            if (zt >= z_bomb - 1.0e-5f)
              {
                /* Tracer inside BMO melt: liquid H2O composition */
                E->tracers.ftrac[ifunc][itrac] = Cw_l;
              }
            else if (E->control.BMO_CH2O_prof != NULL)
              {
                int do_stamp = (flagf >= 0) || (zt >= z_bomb - delta_shell);
                if (do_stamp)
                  {
                    int row = 1, kr;
                    for (kr = 1; kr <= noz_m; kr++) {
                        if (E->X[2][kr] <= zt) row = kr;
                        else break;
                    }
                    if (row < 1)     row = 1;
                    if (row > noz_p) row = noz_p;
                    E->tracers.ftrac[ifunc][itrac] = E->control.BMO_CH2O_prof[row];
                  }
              }
          }
      }
    }
  else if (E->tracers.tracfunc[ifunc] == 83 && E->control.BMO_prof_size > 0)
    {
      float z_bomb = E->control.MO_depth;
      if (z_bomb < 1.0f - 1.0e-6f)
      {
        float Chpe_l = E->control.BMO_Chpe_l;
        int   noz_p  = E->control.BMO_prof_size;
        int   lev_p  = E->mesh.levmax;
        int   noz_m  = E->lmesh.NOZ[lev_p];
        int   flagf  = ifunc_solidflag;
        float delta_shell;
        { int rf = E->control.MO_depth_row;
          if (rf < 2)     rf = 2;
          if (rf > noz_m) rf = noz_m;
          delta_shell = 3.0f * (E->X[2][rf] - E->X[2][rf-1]);
        }

        for (itrac = 0; itrac < E->tracers.lntrac; itrac++)
          {
            if (flagf >= 0 && E->tracers.ftrac[flagf][itrac] > 0.5f)
                continue;

            float zt = E->tracers.xtrac[2][itrac];
            if (zt >= z_bomb - 1.0e-5f)
              {
                E->tracers.ftrac[ifunc][itrac] = Chpe_l;
              }
            else if (E->control.BMO_Chpe_prof != NULL)
              {
                int do_stamp = (flagf >= 0) || (zt >= z_bomb - delta_shell);
                if (do_stamp)
                  {
                    int row = 1, kr;
                    for (kr = 1; kr <= noz_m; kr++) {
                        if (E->X[2][kr] <= zt) row = kr;
                        else break;
                    }
                    if (row < 1)     row = 1;
                    if (row > noz_p) row = noz_p;
                    E->tracers.ftrac[ifunc][itrac] = E->control.BMO_Chpe_prof[row];
                  }
              }
          }
      }
    }
  else if (E->tracers.tracfunc[ifunc] == 82)
    {
      /* Solidified flag (0 = liquid, 1 = solidified): set when a tracer first lies
       * below the front; runs after the 80/81 branches, so the composition is
       * already stamped. */
      float z_bomb = E->control.MO_depth;
      if (E->control.BMO_prof_size > 0 && z_bomb < 1.0f - 1.0e-6f)
        {
          for (itrac = 0; itrac < E->tracers.lntrac; itrac++)
            {
              if (E->tracers.ftrac[ifunc][itrac] <= 0.5f)
                {
                  float zt = E->tracers.xtrac[2][itrac];
                  if (zt < z_bomb - 1.0e-5f)
                      E->tracers.ftrac[ifunc][itrac] = 1.0f;
                }
            }
        }
    }
  else if (E->tracers.interp_e2t[ifunc] == 20)
    {  trapmelt_per_step = trapmelt_add = 0;
       for (itrac=0; itrac<E->tracers.lntrac; itrac++)
        if (E->tracers.ftrac[ifunc][itrac]==1. && E->tracers.xtrac[2][itrac] < TZ_thickness)
          {
	   if  (E->tracers.ftrac[ifunc][itrac]==0.7)
	     trapmelt_per_step++;
           else
	     {
		trapmelt_per_step++;
		trapmelt_add++;
           	E->tracers.ftrac[ifunc][itrac]=0.7;
	     }
	  }
      E->tracers.trapmelt_total+=(((float) (trapmelt_add))/(E->tracers.tracer_density*E->tracers.tracer_density));
      E->tracers.trapmelt_per_step = (((float) (trapmelt_per_step))/(E->tracers.tracer_density*E->tracers.tracer_density));
    }
  else if (E->tracers.interp_e2t[ifunc] == 22)
    {
      trapmelt_per_step = 0;
      for (itrac=0; itrac<E->tracers.lntrac; itrac++)
        if (E->tracers.ftrac[ifunc][itrac]==1. && E->tracers.xtrac[2][itrac] < shallower_km250)
          {
                E->tracers.ftrac[ifunc][itrac]=0.;
                trapmelt_per_step++;
          }
      E->tracers.trapmelt_total+=(((float) (trapmelt_per_step))/(E->tracers.tracer_density*E->tracers.tracer_density));
      E->tracers.trapmelt_per_step = (((float) (trapmelt_per_step))/(E->tracers.tracer_density*E->tracers.tracer_density));
    }
  else if (E->tracers.interp_e2t[ifunc] == 5)
    {
      for (itrac=0;itrac<=E->tracers.lntrac;itrac++)
    E->tracers.ftrac[ifunc][itrac]= func_to_trac(E,E->tracers.xtrac,E->T,
						 xi,E->tracers.elemnr,itrac); 
    }
  else if ((E->tracers.interp_e2t[ifunc] == 3 || E->tracers.interp_e2t[ifunc] ==4) && E->tracers.interp_e2t[ifunc+1] == 6)
   {
    if (PC_iteration==1)
      {
      for (ielem=1;ielem<=E->lmesh.nel;ielem++)
	  {
	    E->Qmelt[ielem]=0.;
	    trapmelt[ielem]=0.;
	  }
      if (((E->monitor.solution_cycles+1)%200==0 || 
	   (E->monitor.solution_cycles-1)%200==0 || 
	   E->monitor.solution_cycles%200==0) && 
	   E->monitor.solution_cycles>1997)
	{ // open File for dFdt-output
	  ndim = E->mesh.nsd;
	  write_dFdt=E->monitor.solution_cycles+1;
          write_dFdt/=200;
          write_dFdt*=200;
	  sprintf(output_file,"%s.dFdt.%d.%d",E->control.data_file,
                   E->parallel.me,write_dFdt);
	  if ((E->monitor.solution_cycles+1)%200==0)
	    fp=fopen(output_file,"w");
	  else
	    fp=fopen(output_file,"a");
	  fwrite(&E->advection.timestep,1,sizeof(float),fp);
	  fwrite(&ndim,1,sizeof(int),fp);              
	  write_dFdt=1;
	}
     }
    if (E->tracers.interp_e2t[ifunc] == 3) melting=0;
	 else melting=1;
    for (itrac=0;itrac<E->tracers.lntrac;itrac++)
     if (E->tracers.xtrac[1][itrac]>.1)
     {
     ielem=E->tracers.elemnr[itrac];

     Ttrac=func_to_trac(E,E->tracers.xtrac,E->T,
			xi,E->tracers.elemnr,itrac);
     Ttrac*=E->data.ref_temperature;
     pressure=factor*E->tracers.xtrac[2][itrac];
     if (melting == 0)

	{

	 H2O=calc_H2O(E,E->tracers.ftrac[ifunc][itrac]); 
         dFdt=katz_melting(E,Ttrac,pressure,H2O);
        }
     else if (melting == 1) 
       dFdt=hirschmann_melting(E,Ttrac,pressure);

     dFdt-=E->tracers.ftrac[ifunc][itrac];
     

     if (dFdt>0.)
	{
	  if  (E->data.heat_of_melt !=0.)
	   {
	   if (melting == 0)
	    {
	      /* precorrection of dF due to consumption of latent heat */
	      p=(int)(pressure / 0.05);
	      T = ((int)((Ttrac-E->tracers.Fcorr[p][0])/2.))-1;
	      if (T<1) T=1;
	      dFdt*=E->tracers.Fcorr[p][T];
	    }
	   else if (melting == 1)
	     dFdt/=HperCpA;

	   }
	   MR=E->tracers.ftrac[ifunc+1][itrac]; /* melt retention */
	   MR+=dFdt;
	   if (PC_iteration==1)
	     {
	      E->Qmelt[ielem]+=dFdt;

	      if (write_dFdt)
	       {
		for (idim=1; idim <= ndim; idim++)
		  fwrite(&E->tracers.xtrac[idim][itrac],1,sizeof(float),fp);
	        fwrite(&dFdt,1,sizeof(float),fp);
	       }
	     }
	   E->tracers.ftrac[ifunc][itrac]+=dFdt;
	   if (MR>E->tracers.phiC) 
	     {
	     if (PC_iteration==1)
	      {/* melt extraction to static trap*/
	       av_visit=E->tracers.ntrac_in_elem[ielem]
		       -E->tracers.ntrac_in_elem[ielem-1];
	       trapmelt[ielem]+=(MR-E->tracers.phiR)/av_visit;

	       if (E->tracers.melttracers_used !=0)   
	       /* idem melt extraction to dynamic trap (melttracers)*/
		 { //computing column and row ....
		   xoffset = E->tracers.xtrac[1][itrac]- 
		     E->tracers.melt_xbound[E->tracers.melt_firstcolumn];
		   xoffset /= E->tracers.melt_xdistance;
		   c=((int)(xoffset)) + E->tracers.melt_firstcolumn;
		   if (c>columns) c-=columns+1;
		   //if (c>=columns || c<0) fprintf(E->fp,"ERR: c:%d",c);
		   if (E->mesh.nsd==3)
		     { yoffset=E->tracers.xtrac[3][itrac]-E->tracers.melt_ybound[0];
		       yoffset /= E->tracers.melt_ydistance;
		       r=((int)(yoffset));
		     }
		   else r=0;
 		   E->tracers.melt[c][r]+=(MR-E->tracers.phiR)/av_visit;
		 }

	       /* melt extracted */

	      }
	     MR=E->tracers.phiR;
	     }
	     // if (E->monitor.solution_cycles>1)
	     E->tracers.ftrac[ifunc+1][itrac]=MR;
	     // else
	     // E->tracers.ftrac[ifunc+1][itrac]=0.;
        }
     }
     else
       {
         if (E->mesh.nsd == 3) 
          y=E->tracers.xtrac[3][itrac];
         else y=0.;
         E->tracers.ftrac[ifunc][itrac] = trac_funcinit(E, E->tracers.tracfunc[ifunc], E->tracers.xtrac[1][itrac], y, E->tracers.xtrac[2][itrac]);
	 E->tracers.ftrac[ifunc+1][itrac] = trac_funcinit(E, E->tracers.tracfunc[ifunc+1], E->tracers.xtrac[1][itrac], y, E->tracers.xtrac[2][itrac]);
       }


    
    E->tracers.melting_per_step=0.;
    if (PC_iteration==1)
     for (ielem=1;ielem<=E->lmesh.nel;ielem++)
     { 
       av_visit=(E->tracers.ntrac_in_elem[ielem]-E->tracers.ntrac_in_elem[ielem-1]);
       if (E->Qmelt[ielem]!=0.)
	{
	 if (av_visit>0)
	  {
	    E->Qmelt[ielem] /= av_visit;
	    E->tracers.melting_per_step += E->Qmelt[ielem]*
	                                   E->ECO[lev][ielem].area;
            if (E->Qmelt[ielem]!=0.)
	      E->Qmelt[ielem] *= factor2;
	    else E->Qmelt[ielem]=0.;
	  }
	 else E->Qmelt[ielem]=0.;
	}

      if (trapmelt[ielem]>0.0 && E->monitor.solution_cycles>1)
	{  // melt column height dF*z
          trapmelt[ielem]/=av_visit;
	  column=ielem/NEL_in_a_column;
          trapmelt[ielem]*=E->ECO[lev][ielem].size[2];
	  E->tracers.trapmelt[column]+=trapmelt[ielem];
	}  
     }
   if (write_dFdt && (E->monitor.solution_cycles-1)%200==0)
    {
      for (idim=1; idim<=ndim; idim++)
	fwrite(&null,1,sizeof(float),fp);
      fwrite(&null,1,sizeof(float),fp);
      fclose(fp);
    }
  }
  else if (E->tracers.interp_e2t[ifunc] == 3)
    {
       for (itrac=0; itrac<E->tracers.lntrac; itrac++) 
         E->tracers.ftrac[ifunc][itrac]=func_to_trac(E,
	           E->tracers.xtrac,E->T,xi,E->tracers.elemnr,itrac);  
    }
  }
 return;
}


float calc_H2O(E,F)
   struct All_variables *E;
   float F;
{
  float H2O;
  int extraction=E->tracers.extraction;

  if (extraction==1)
    {
      H2O=exp(-54.303*F);
      H2O*=1.8*E->data.H2O/0.02;
    }
  else if (extraction==2)
    {
      H2O=exp(-50.95*F);
      H2O*=1.7477*E->data.H2O/0.02;
    }
  else if (extraction==3)
    {
      H2O=exp(-46.718*F);
      H2O*=1.6721*E->data.H2O/0.02;
    }
  else if (extraction==4)
    {
      H2O=exp(-41.9*F);
      H2O*=1.5884*E->data.H2O/0.02;
    }
  else if (extraction==20)
    {
      H2O=1/(.5+49.5*F);
      H2O*=2.*E->data.H2O/0.02;
    }
  else
    H2O=E->data.H2O;

  return (H2O);
}


float katz_melting(E, Ttrac, pressure, H2O)
   struct All_variables *E;
   float Ttrac;
   float pressure;
   float H2O;
{
   // const float H2O=.02; /***** wt-% *****/
   float F=0., TsolP, TlherzLiq, Tprime;
   float TliqP, TcpxOut, FcpxOut;
   float deltaT;
    
   deltaT=43*pow(H2O,0.75);

   if (pressure < 5.05313)
     {
       TsolP = 1085.7 +(132.9*pressure) - (5.1 * pressure*pressure) - deltaT;
   if (Ttrac > TsolP)
      {
       /*printf ("Solidus at: %f T= °C %f ;\n", TsolP,Ttrac);*/
       TlherzLiq = 1475.0+(80.0*pressure)-(3.2 * pressure*pressure)-deltaT;
       Tprime = (Ttrac - TsolP) / (TlherzLiq - TsolP);
       F = Tprime*sqrt(Tprime);

       FcpxOut = 0.1/(0.5+0.08*pressure);
       
       if (F > FcpxOut)
	F=FcpxOut;
       /* {
	 TcpxOut=pow(FcpxOut,0.666667)*(TlherzLiq-TsolP)+TsolP;
	 TliqP=1780.0+45.*pressure-2.*pressure*pressure-deltaT;

	 F=(Ttrac-TcpxOut)/(TliqP-TcpxOut);
	 F=FcpxOut+(1-FcpxOut)*F*sqrt(F);
	 if (E->data.heat_of_melt !=1.) fprintf (stderr,"Warning: precorrection for Qmelt is not longer valid above Tcpxout\n");
   
	}*/
      }
     }
return (F);
}


float hirschmann_melting(E, Ttrac, pressure)
   struct All_variables *E;
   float Ttrac;
   float pressure;
{
   float F, TsolP;
   F=0.;
   if (pressure<=6.)
    {

     TsolP = 1096.154 +(146.533*pressure)-(6.435*pressure*pressure);
     Ttrac-= TsolP;
     if (Ttrac > 0)
      {
	F = .002578947 * Ttrac;
      }
    }    
    
return (F);
}
