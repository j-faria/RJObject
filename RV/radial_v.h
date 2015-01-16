#include<cmath>
#include<vector>
namespace Suto
{
  const double D2S=86400.0;
  double ecc_anomaly(double time, double prd,double ecc,double peri_pass);
  double true_anomaly(double time,double prd,double ecc,double peri_pass);
  double radial_velocity(double time,
			       double T,
			       double K,
			       double e,
			       double w,
			       double X);

  double rad_v(double time, 
  	           double prd, double amp, double ecc, double omega, double chi, double vsys, 
  	           int n_planets);

}
