#include "sde_mult.h"
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <stdio.h>
#include <vector>

/*initializer*/
SDE::SDE(
	 int dimention,
	 vector_field v0,
	 vector<vector_field> v1,
	 void *params_setting 
	 )
{
  drift = v0; // including only one term
  diffusion = v1; // including possibly many terms 
  dim = dimention; // dimention of the space is explicitly given
  dim_diffusion = diffusion.size(); // num of diffusion terms is inplicitly identified
  sde_params = params_setting;
  
  const gsl_rng_type *T;
  T = gsl_rng_default;
  r = gsl_rng_alloc (T); // r is a gsl random sequence
  gsl_rng_set(r,time(NULL)); // setting the random seed
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
  //drift term
  gsl_vector *y1 = gsl_vector_alloc(dim);
  //a vector to save the sum of diffusion terms
  gsl_vector *y2 = gsl_vector_alloc(dim);

  /*drift term*/
  gsl_vector_memcpy(y1,tmp_y);
  (drift)(y1,sde_params);
  gsl_vector_scale(y1,step_w);

  int i=0;
  //diffution term
  gsl_vector_memcpy(y2,tmp_y);
  gsl_vector_set_zero(y2);
  
  for(int i=0;i<dim_diffusion;i++)
    {
      // note that i indicates i th term of the diffusion terms
      //alloc memory
      gsl_vector *i_th_y2 = gsl_vector_alloc(dim);
      gsl_vector_memcpy(i_th_y2,tmp_y);      
      gsl_vector *gaussian_noise = gsl_vector_alloc(dim);
      //calulation
      (diffusion[i])(i_th_y2,sde_params);
      set_gaussian_noise(gaussian_noise);
      gsl_vector_scale(gaussian_noise,sqrt(step_w));
      gsl_vector_mul(i_th_y2,gaussian_noise);
      gsl_vector_add(y2,i_th_y2);
      // free memory
      gsl_vector_free(i_th_y2);
      gsl_vector_free(gaussian_noise);   
    }
  

  /*sum all*/
  gsl_vector_add(tmp_y,y1);
  gsl_vector_add(tmp_y,y2);

  /*free memory*/
  gsl_vector_free(y1);
  gsl_vector_free(y2);
}

void SDE::simulate_one_chain(
			     gsl_vector *tmp_y,
			     int step_n,
			     double T
			     )
{
  double step_w;
  step_w=T/step_n;

  int i=0;
  for(i=0;i<step_n;i++)
    {
      simulate_one_step(tmp_y,step_w);
      printf("%dth y\n",i);
      gsl_vector_fprintf(stdout,tmp_y,"%E");
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

  printf("Maturity:%d",T);
  printf("step_width:%d",T/step_n);

  int i =0;
  
  for(i=0;i<sample_size;i++)
    {
      gsl_vector *tmp_y = gsl_vector_alloc(dim);
      gsl_vector_memcpy(tmp_y,init_y);
      simulate_one_chain(tmp_y,step_n,T);
      /*for test*/
      //printf("%dth y\n",i);
      //gsl_vector_fprintf(stdout,tmp_y,"%E");

      gsl_vector_add(sum_y,tmp_y);
      gsl_vector_free(tmp_y);
    }
  gsl_vector_scale(sum_y,1.0/sample_size);
  //gsl_vector_scale(sum_y);
  printf("total");
  gsl_vector_fprintf(stdout,sum_y,"%E");
  
  gsl_vector_free(sum_y);
}

void SDE::set_gaussian_noise(gsl_vector *y)
{
  int i=0;
  
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


// vector[0]=S(t)=StockPrice
// vector[1]=V(t)=Volatility
// vector[2]=OptionPrice

// mu:return
// alpha:tense of regression
// theta:regression level
// rho: correlation between the processes of Price and Volatility

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
  gsl_vector_set(y,0,y0*sqrt(y1<0 ? 0 : y1));
  gsl_vector_set(y,1,ah_params->beta * sqrt(y1<0 ? 0 : y1) * ah_params->rho);
  gsl_vector_set(y,2,0);
}


void ah_v2(gsl_vector *y,void *params)
{
  struct ah_params_s *ah_params;
  /*convert the type of void pointer*/
  ah_params=reinterpret_cast<ah_params_s*>(params);

  double y0,y1,y2;
  y0=gsl_vector_get(y,0);
  y1=gsl_vector_get(y,1);
  y2=gsl_vector_get(y,2);
  gsl_vector_set(y,0,0);
  gsl_vector_set(y,1,ah_params->beta * sqrt(y1<0 ? 0 : y1*(1-pow(ah_params->rho,2))));
  gsl_vector_set(y,2,0);
}

int main(void)
{
  struct ah_params_s ah_params;
  ah_params.mu=1.12;
  ah_params.alpha=1.21;
  ah_params.theta=0.34;
  ah_params.beta=0.944;
  ah_params.rho=0.421;
  // 1-dimentional vector v1
  vector<vector_field> v1;
  v1.push_back(&ah_v1);
  v1.push_back(&ah_v2);
  
  SDE ah(3,&ah_v0,v1,&ah_params);
  gsl_vector *init_y = gsl_vector_alloc(3);
  printf("start");
  
  gsl_vector_set(init_y,0,100);
  gsl_vector_set(init_y,1,0.21);
  gsl_vector_set(init_y,2,0);
  ah.simulate_chains(init_y,100,100,1);
  gsl_vector_free(init_y);
  
}
