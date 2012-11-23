#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "em1.h"

double em1_1step(sde_1d *sde,
		 double start,
		 double step_w,
		 double normal_rand){
  return (start+(sde->v0t)(start)*step_w+(sde->v1)(start)*sqrt(step_w)*normal_rand);
}

double em1(sde_1d *sde,
	   double init,
	   double T,
	   int n,
	   double (*rand)(void)
	   )
{
  double point,step_w;
  point=init;
  step_w=T/n;
  printf("%f",step_w);  
  int i;
  for(i=0;i<n;i++){
    point=em1_1step(sde,point,step_w,(*rand)());
    printf("%d:%f\n",i,point);    
  }
  return point;
}

/* sample code */

double mu = 0.0;
double sigma = 0.001;

double bs1_fv0(double y){
  return (mu-0.5*sigma*sigma)*y;
}

double bs1_fv0t(double y){
  return mu*y;
}

double bs1_fv1(double y){
  return sigma*y;
}

double Uniform( void ){
  return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}

double Normal(void){
  srand( (unsigned int)time( NULL ) );
  return sqrt( -log(Uniform())*2.0 ) * sin( 2.0*M_PI*Uniform() );
}

int main(void){
  /*definition of BS equation*/
  sde_1d bs1;
  bs1.v0=&bs1_fv0;
  bs1.v0t=&bs1_fv0t;
  bs1.v1=&bs1_fv1;
  bs1.x=1.0;
  em1(&bs1,100.0,10,100,&Normal);
  return 0;
}
