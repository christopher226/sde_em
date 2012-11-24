#include "sde_mult.h"
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <stdio.h>
     

/*initializer*/
SDE::SDE(
	 int dimention,
	 void (*v0)(gsl_vector *,void *),
	 void (*v1)(gsl_vector *,void *),
	 void *params_setting 
	 )
{
  V[0]=v0;
  V[1]=v1;
  dim = dimention;
  sde_params = params_setting;

  const gsl_rng_type *T;
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set(r,time(NULL));
}

/*finalizer*/
SDE::~SDE()
{
  gsl_rng_free(r);
}

void SDE::simulate_one_step(
			    gsl_vector *tmp_y,
			    double step_w)
{
  gsl_vector *y1 = gsl_vector_alloc(dim);
  gsl_vector *y2 = gsl_vector_alloc(dim);
  gsl_vector *gaussian_noise =gsl_vector_alloc(dim);
  
  gsl_vector_memcpy(y1,tmp_y);
  gsl_vector_memcpy(y2,tmp_y);

  /*drift term*/
  (V[0])(y1,sde_params);
  gsl_vector_scale(y1,step_w);

  /*diffution term*/
  (V[1])(y2,sde_params);
  set_gaussian_noise(gaussian_noise);
  gsl_vector_scale(gaussian_noise,sqrt(step_w));
  gsl_vector_mul(y2,gaussian_noise);  

  /*sum all*/
  gsl_vector_add(tmp_y,y1);
  gsl_vector_add(tmp_y,y2);

  /*free memory*/
  gsl_vector_free(y1);
  gsl_vector_free(y2);
  gsl_vector_free(gaussian_noise);
}

void SDE::simulate_one_chain(
			     gsl_vector *tmp_y,
			     int step_n,
			     double T
			     )
{
  double step_w;
  step_w=T/step_n;

  int i;
  for(i=0;i<step_n;i++)
    {
      simulate_one_step(tmp_y,step_w);
    }
  
}

void SDE::simulate_chains(
			  gsl_vector *init_y,
			  int sample_size,
			  int step_n,
			  double T
			  )
{
  /*vector to sum up*/
  gsl_vector *sum_y = gsl_vector_alloc(dim);
  gsl_vector_set_all(sum_y,0);
  
  int i;
  for(i=0;i<sample_size;i++)
    {
      gsl_vector *tmp_y = gsl_vector_alloc(dim);
      gsl_vector_memcpy(tmp_y,init_y);
      simulate_one_chain(tmp_y,step_n,T);
      /*for test*/
      printf("%dth y\n",i);
      gsl_vector_fprintf(stdout,tmp_y,"%E");
      gsl_vector_add(sum_y,tmp_y);
      gsl_vector_free(tmp_y);
    }
  
  gsl_vector_scale(sum_y,1/sample_size);
  gsl_vector_free(sum_y);
}

void SDE::set_gaussian_noise(gsl_vector *y)
{
  
  
  int i;
  for(i=0;i<dim;i++)
    {
      gsl_vector_set(y,i,gsl_ran_gaussian(r,1));
    }
  
  
}

struct ah_params_s
{
  double mu;
  double alpha;
  double theta;
  double beta;
  double rho;
};


void ah_v0(gsl_vector *y,void *params)
{
  struct ah_params_s *ah_params;
  /*convert the type of void pointer*/
  ah_params=reinterpret_cast<ah_params_s*>(params);

  double y0,y1,y2;
  y0=gsl_vector_get(y,0);
  y1=gsl_vector_get(y,1);
  y2=gsl_vector_get(y,2);
  gsl_vector_set(y,0,ah_params->mu * y0);
  gsl_vector_set(y,1,ah_params->alpha * (ah_params->theta - y1));
  gsl_vector_set(y,2,y0);
}


void ah_v1(gsl_vector *y,void *params)
{
  struct ah_params_s *ah_params;
  /*convert the type of void pointer*/
  ah_params=reinterpret_cast<ah_params_s*>(params);

  double y0,y1,y2;
  y0=gsl_vector_get(y,0);
  y1=gsl_vector_get(y,1);
  y2=gsl_vector_get(y,2);
  gsl_vector_set(y,0,y0*sqrt(y1));
  gsl_vector_set(y,1,ah_params->beta * sqrt(y1) * ah_params->rho);
  gsl_vector_set(y,2,0);
}


int main(void)
{
  struct ah_params_s ah_params;
  ah_params.mu=1;
  ah_params.alpha=1;
  ah_params.theta=1;
  ah_params.beta=1;
  ah_params.rho=1;
  SDE ah(3,&ah_v0,&ah_v1,&ah_params);
  gsl_vector *init_y = gsl_vector_alloc(3);
  gsl_vector_set(init_y,0,1);
  gsl_vector_set(init_y,1,1);
  gsl_vector_set(init_y,2,2);
  ah.simulate_chains(init_y,1000,100000,10);
  gsl_vector_free(init_y);
  
}
