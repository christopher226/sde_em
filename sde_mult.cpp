#include "sde_mult.h"
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
     

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
}

/*finalizer*/
SDE::~SDE()
{
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
      simulate_one_chain(tmp_y,step_n,T);
      gsl_vector_add(sum_y,tmp_y);
      gsl_vector_free(tmp_y);
    }
  
  gsl_vector_scale(sum_y,1/sample_size);
  gsl_vector_free(sum_y);
}

void SDE::set_gaussian_noise(gsl_vector *y)
{
  const gsl_rng_type *T;
  gsl_rng *r;  
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  int i;
  for(i=0;i<dim;i++)
    {
      gsl_vector_set(y,i,gsl_ran_gaussian(r,1));
    }
  gsl_rng_free (r);
}

struct ah_params_s
{
  double alpha;
  double theta;
};


void ah_v0(gsl_vector *y,void *params)
{
  struct ah_params_s *ah_params;
  ah_params=params;
  ah_params->alpha;
}


void ah_v1(gsl_vector *y,void *params)
{
  struct ah_params_s *ah_params;
  ah_params=params;
  ah_params->alpha;

}


int main(void)
{
  struct ah_params_s ah_params;
  ah_params.alpha=1;
  ah_params.theta=1;  
  SDE ah(1,&ah_v0,&ah_v1,&ah_params);
}
