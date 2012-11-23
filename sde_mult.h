#include <gsl/gsl_vector.h>

class SDE
{
 public:
  /*constructor*/
  SDE(
      int dimention,
      void (*v0)(gsl_vector *y,void *params),
      void (*v1)(gsl_vector *y,void *params),
      void *
      );
  ~SDE();
  
  /*methods to simulate the SDE*/
  void simulate_one_step(gsl_vector *tmp_y,
			 double step_w);
  void simulate_one_chain(gsl_vector *tmp_y,
			  int step_n,
			  double T);
  void simulate_chains(int sample_size,
		       int step_n,
		       double T);
  
 private:
  
  /*parameters of SDE*/
  void (*V[3])(gsl_vector *,void *);
  void *sde_params;
  int dim;

  void set_gaussian_noise(gsl_vector *y);
  
  /*methods to simulate Wiener process*/
  double MultivariateNormal(int);  
};

