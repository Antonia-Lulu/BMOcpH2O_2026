#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"

float Zj0[1000],Zj1[1000];

void set_axi_defaults(E)
     struct All_variables *E;
{ 
#include  "Zeros_of_J0_and_J1.h"
  E->mesh.nsd = 2;
  E->mesh.dof = 2;
  
  /* add functions for hankel transforms &c */

}

void hankel_transform(E,data0,transform)
     struct All_variables *E;
     double data0[],transform[];     
{
    int i,j,points ;
    double k,integrand,data1,data2,norm,horizontal_1scale;
    float data[2000],Data[2000];
    float return_layer_value();

    points = E->mesh.nnx[1];
    horizontal_1scale = 1.0 / (E->X[1][E->mesh.nno]);
    
    norm=(double)return_layer_value(E,Data,1,1);
    for(j=1;j<=points;j++)
	data[j] = data0[j]-norm;	/*	Improves the accuracy of the integration to 
						use the axial value as the reference level.	*/
    for(i=0;i<points;i++)  {
	k=Zj1[i]*horizontal_1scale;	/*  Wavenumbers for the problem   */
	norm =0.5* j0(Zj1[i]*horizontal_1scale)*j0(Zj1[i]*horizontal_1scale);

	for(j=1;j<=points;j++) 
	    Data[j] = j0(k*E->X[1][1+(j-1)*E->mesh.noz])*data[j];

	integrand=(double)return_layer_value(E,Data,1,0);
	for(j=1;j<=points;j++) 
	    Data[j] = j0(k*E->X[1][1+(j-1)*E->mesh.noz])*j0(k*E->X[1][1+(j-1)*E->mesh.noz]);
	norm=(double)return_layer_value(E,Data,1,0);

	transform[i+1] = integrand/norm;
	}	
    transform[1] = 0.0;	
    return;	
}

void inverse_hankel_transform(E,data,transform)
     struct All_variables *E;
     double data[],transform[];     
{
    int i,j;
    double x,integrand,data1,data2,norm1,norm2;
    double horizontal_1scale;
    int points;

    points = E->mesh.nnx[1];
    horizontal_1scale = 1.0/E->X[1][E->mesh.nno];

    for(i=0;i<points;i++)    {
	x = E->X[1][1+(i-1)*E->mesh.noz];	
	integrand = 0.0;
	for(j=1;j<points;j++) {
	    norm1 =  0.5*j0(Zj1[j-1]*horizontal_1scale)*j0(Zj1[j-1]*horizontal_1scale);
	    norm2 =  0.5*j0(Zj1[j]*horizontal_1scale)*j0(Zj1[j]*horizontal_1scale);
	    data1=norm1*j0(x*Zj1[j-1]*horizontal_1scale)*data[j];
	    data2=norm2*j0(x*Zj1[j]*horizontal_1scale)*data[j+1];
	    integrand += (data1 + data2) * 0.25 * (Zj1[j]*Zj1[j]-Zj1[j-1]*Zj1[j-1])*
		(horizontal_1scale*horizontal_1scale);		}
	transform[i+1] = integrand;	}
    return;	
}

