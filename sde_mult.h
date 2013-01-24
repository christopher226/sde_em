#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>
using namespace std;

// what vector field is all about
// basically generating a new vector field by side-effect
// *params contains parameter values used in this function
typedef void (*vector_field)(gsl_vector *y,void *params);

// definition SDE problem
class SDE
{
 public:
  /*constructor*/
  // dimention is not automatically identified,
  // but explicitly indicated as the parameter dimention
  SDE(
      int dimention,
      vector_field, // drift term
      vector<vector_field>, // diffusion term
      void *
      );
  ~SDE();
  
  /*methods to simulate the SDE*/
  // executing simulation by side-effect like other methods in gsl library
  void simulate_one_step(gsl_vector *tmp_y,
			 double step_w);
  void simulate_one_chain(gsl_vector *tmp_y,
			  int step_n,
			  double T);
  void simulate_chains(
		       gsl_vector *init_y,
		       int sample_size,
		       int step_n,
		       double T);
  
 private:
  
  /*parameters of SDE*/
  vector_field drift;
  vector<vector_field> diffusion;  
  void *sde_params;
  int dim; // dimention of the space
  int dim_diffusion; // num of diffusion terms

  /*parameters of gaussian noise*/
  gsl_rng *r;

  void set_gaussian_noise(gsl_vector *y);
  
};

