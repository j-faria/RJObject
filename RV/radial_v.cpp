#include<cmath>
#include<iostream>
#include"radial_v.h"

// this function calculates eccentric anomaly
double Suto::ecc_anomaly(double time, double prd, double ecc,double peri_pass)
  {
    time=time*D2S;
    double T=prd*D2S;
    double n=2*M_PI/T;
    double tau=peri_pass*T;
    double M=n*(time+tau);
    double E0;// start value
    double Mstar;// check equation  6.6.9 Danby
    double sigma;
    int  int_M2PI;
    // to get an interger value for M/(2*PI)
    // so that for M we have a value between 0 and 2PI
    int_M2PI=M/(2*M_PI);
    Mstar=M-int_M2PI*2*M_PI;

    // define a SIGN function
    if(sin(Mstar)<0)
      sigma=-1*sin(Mstar);
    else 
      sigma=sin(Mstar);

    E0=Mstar+sigma*0.85*ecc;
  
    // the value for k=0.85 is arbitrary
    // the only condition is 0<= k <=1 check, again Danby   
  
    double TINY=1e-6;
    int count=0;
    for(;;){
    
      double eSinE=ecc*sin(E0);    // a dummy
      double f=E0-eSinE-Mstar;
      if(fabs(f)<TINY || count >100){
	break;
      }
      double eCosE=ecc*cos(E0);
      double f1=1-eCosE;
      double f2=eSinE;
      double f3=eCosE;
      double dE0=-f/f1;
      dE0=-f/(f1+0.5*dE0*f2);
      dE0=-f/(f1+0.5*dE0*f2+dE0*dE0*f3/6);
      E0=E0+dE0;
      ++count;
    }
  
    return E0;
  }

// this function calculates true anomaly
double Suto::true_anomaly(double time,double prd,double ecc,double peri_pass)
  {

    double E=Suto::ecc_anomaly(time,prd,ecc,peri_pass);
    double f=acos( (cos(E)-ecc)/( 1-ecc*cos(E) ) );
    //acos gives the principal values ie [0:PI]
    //when E goes above PI we need another condition
    if(E>M_PI)
      f=2*M_PI-f;
    return f;
  }

// function that calculates the radial velocity
double Suto::radial_velocity(double time,
			     double T,
			     double K,
			     double e,
			     double w,
			     double X)
  {
  
    double f=Suto::true_anomaly(time,T,e,X);
    double v=K*( sin(f+w)+e*sin(w) );
    return v;
  }

double Suto::rad_v(double time,
                   double prd, double amp, double ecc, double omega, double chi, double vsys, 
                   int n_planets)
{
  
  //if(params.size()==5*n_planets+1)
  //{
    //std::vector<double> temp_para;
    //double sys_vel=params[0];
    double sys_vel = vsys;
    double radial_v=sys_vel;

    for(int i=0;i<n_planets;++i)
    {
  //    double prd=params[i*5+1];
  //    double amp=params[i*5+2];
  //    double ecc=params[i*5+3];
  //    double omega=params[i*5+4];
  //    double chi=params[i*5+5];
      radial_v-=Suto::radial_velocity(time,prd,amp,ecc,omega,chi);
      //temp_para.clear();
    }
      
    return radial_v;

  //}
  //else
  //{
  //  std::cout<<"Parameters does not match the number of planets\n";
  //  return 0;
  //}
}
	  
