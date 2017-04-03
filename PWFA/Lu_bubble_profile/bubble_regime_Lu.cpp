/*
 *  Created by F. Mira A. Marocchino C. Gatti
 */

#include "parameters.h"
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <stdlib.h> 

typedef double (*func)(double,double,double);

double density(double);
double density_driver(double);
double density_witness(double);
void compute_beta(double&, double);
void compute_Dbeta(double&,double&, double);
void compute_DDbeta(double&, double&,double&, double);
void compute_coefficients(double&,double&,double&,double,double);
double compute_Ez( double, double, double);

double g(double,double,double);
double f(double,double,double);
void RK4traj(func,func,double&,double&,double&,double);



const double Ebeam=100;//MeV
const double Lchannel=0.03;
const double mu=5;//driver position
double alpha;
double ComputeAlpha(double,double);

int main(int argc, char **argv) {
  
  if (argc<8){
    std::cout<<"Usage: ./bubble.exe n0(10^16 cm-3) Q(driver pC) sigma_z(driver mum) sigma_r(driver mum) Q(witness pC) sigma_z(witness mum)  twd(distance witness from driver mum)  ScanFlag(yes/no)"<<std::endl;    
    return 0;
  }  
  for (int i=0;i<argc;i++){
    std::cout<<"ARG["<<i<<"]="<<argv[i]<<std::endl;        
  } 

  /** perform parameter scan **/
  std::string ScanFlag;
  if (argc==9) {
    ScanFlag=argv[8];
  } else {
    ScanFlag="no";  
  }

  /** Units
      n --> n/n0
      L --> kp L
      E --> eE/(kp me c^2)  
      Q -> Qkp^3/n0e
  **/
  
  n0 = atof(argv[1])*1e16;									//background density input in cm^-3
  
  wp=pow(4.*pi*re*clight*clight*n0*1.e6,0.5);				//plasma frequency in rad/s
  kp=wp/clight;												//kp in m^-1
  lp=2.*pi*clight/wp;										// plasma wavelength in m
  std::cout<<"lambda_p(mum) "<<lp/l0<<std::endl;
  std::cout<<"kp(mum^-1) "<<kp*l0<<std::endl;
  
  /** conversione factors**/

  double Qnorm=(qe*n0*1e6)/pow(kp,3);						//charge in a skin depth volume in C
  double Efield=me*kp;										//MV/m

  Qd     = atof(argv[2])*pC/Qnorm;							//input in pC, Qnorm in C
  Sz_d   = kp*atof(argv[3])*l0;								//input in mum, kp in m^-1
  Sr_d   = kp*atof(argv[4])*l0;								//input in mum, kp in m^-1
  double Qw0     = atof(argv[5])*pC/Qnorm;					//input in pC, Qnorm in C
  double Sz_w0   = kp*atof(argv[6])*l0;						//input in mum, kp in m^-1
  double tm_w0   = kp*atof(argv[7])*l0;						//input in mum, kp in m^-1  


  /** initial values for scan parameters **/

  Qw=Qw0;
  Sz_w=Sz_w0;
  tm_w=tm_w0;
  
  //

  double ndriver=Qd/pow(2*pi,1.5)/Sz_d/Sr_d/Sr_d;			//maximum value of the density of the driver
  double rmax=2.*Sr_d*pow(ndriver,0.5);			 			//bubble dimension

  std::cout<<"Rmax= "<<rmax<<std::endl;
  std::cout<<"ndriver(adimensional)= "<<ndriver<<std::endl;
  std::cout<<"ndriver(10^16 cm^-3)= "<<ndriver*n0*1.e-16<<std::endl;
  std::cout<<"ndriver (radially integrated)= "<<ndriver*2*pi*Sr_d*Sr_d<<std::endl;
  std::cout<<"Q0(pC)= "<<(qe*n0*1e6)/pow(kp,3)/pC<<std::endl;
 
 
  //-apri file di output-//
 
  char fname[1024];  
  sprintf(fname,"bubble_RK4.dat");
  
  std::ofstream outfile(fname);  
  
 
  /** corrispondenza nomi vecchie variabili
     f=function_r
     g=function_u
     y=ub
     h=Dxi
     t=xi
     x=rb
   **/

  double t,x,y,h;  


  
  
  int NScan=20;//point scan
  if (ScanFlag=="no") NScan=1;
  double DeltaQw=3./2.*Qw/double(NScan);
  double DeltaXi=1./10.*tm_w/double(NScan);
for(int ixi=0;ixi<NScan;ixi++) {
    tm_w=tm_w0-(ixi-NScan/2)*DeltaXi;
	
	for (int iq=0;iq<NScan;iq++) {
		Qw=Qw0-(iq-NScan/2)*DeltaQw;

		int nmax=1e9;
		int it=0;      
		t = 0.; x = Sr_d*1e-1; y = 0; h =1e-5;    
		
		double Eav=0,SigmaE=0,Qtot=0,Eznow=0,den=0,tm_wnow=0;
    
		while (it<nmax&&x>0) {
			RK4traj(f,g,t,x,y,h);
			it++;
			Eznow=compute_Ez(x,t,y);
			den=density_witness(t);
			if (int(it)%100==0&&ScanFlag=="no") outfile<<t<<"\t"<<x<<"\t"<<density(t)*2*pi<<"\t"<<Eznow<<"\n";        
			if (Eznow!=Eznow) {
			// for f NAN:  f!=f is true
				std::cout<<"in t="<<t<<" Eznow is "<<Eznow<<std::endl;        
			} 
			else {
				Eav+=Eznow*den*h;
				SigmaE+=Eznow*Eznow*den*h;
				Qtot+=den*h;
			}
		}
		Eav=-Eav/Qtot;
		SigmaE=pow(SigmaE/Qtot-Eav*Eav,0.5)*Efield*Lchannel;
		Eav=Eav*Efield*Lchannel;
		double SEoverE=SigmaE/(Ebeam+Eav);
		Qtot=Qw*Qnorm/pC;
		tm_wnow=tm_w/kp/l0;
		std::cout<<"Witness position (mum)= "<<tm_wnow<<" Qwitness(pC)= "<<Qtot<<" Ebeam(MeV)="<<Ebeam<<" <DEnergy>(MeV)="<<Eav
	     <<" SigmaEnergy(MeV)="<<SigmaE
	     <<" SigmaE/E="<<SEoverE
	     <<std::endl;        
    
	if (ScanFlag=="yes") outfile<<tm_wnow<<"\t"<<"\t"<<Qtot<<"\t"<<"\t"<<Eav<<"\t"<<"\t"<<SEoverE<<std::endl;          
    
		std::cout<<"crossing point reached at iteration n "<<it<<std::endl;
	}
}//Scan
    

  outfile.close();
  
  return 0;
}






void RK4traj(func f,func g,double& t,double &x,double& y,double h) {
  
  double s=x,r=y,q=t;
  double K1 = f(q,s,r),L1=g(q,s,r);
  s= x+(h/2)*K1; r=y+(h/2)*L1; q=t+h/2;
  double K2 = f(q,s,r),L2=g(q,s,r);  
  s= x+(h/2)*K2; r=y+(h/2)*L2; q=t+h/2;
  double K3 = f(q,s,r),L3=g(q,s,r);
  s= x+h*K3; r=y+h*L3; q=t+h;
  double K4 = f(q,s,r),L4=g(q,s,r);
  x += h*(K1+2*(K2+K3)+K4)/6;
  y += h*(L1+2*(L2+L3)+L4)/6;
  t += h;
}

double f(double t,double x,double y) {
  return y;
}

double g(double t,double x,double y) {
  double fu,A,B,C;
  double L = density(t);
  
  compute_coefficients(A,B,C,x,t);
  return (L/x - C*x -B*x*y*y)/A;
}

double density( double t ){
  /** extra 1/(2*pig) since Q= 2 pig \int lambda dxi **/
  return (density_driver(t)+density_witness(t))/(2*pi);
}

double density_driver( double t ){
  
  double sigma, rho;  
  sigma = Sz_d;

  rho = exp(-(t-mu)*(t-mu)/2./sigma/sigma);
  rho = Qd*rho/pow(2.*pi,0.5)/sigma;

  return rho;
}

double density_witness( double t ){
  
  double sigma, rho;
  
  sigma = Sz_w; 

  rho = exp(-(t-mu-tm_w)*(t-mu-tm_w)/2./sigma/sigma);
  rho = Qw*rho/pow(2.*pi,0.5)/sigma;

  return rho;
}

double ComputeAlpha(double r,double t) {
  double DS=0.1;
  double DL=1;
  /**
     if (t>mu+3){
     DL=1-(t-mu-3)*0.8/4.;
  }
  **/
  alpha=DL/r+DS;
  return alpha;
}

void compute_beta(double &beta, double r){
  //alpha=1./r+0.1;
  beta = (1.+alpha)*(1.+alpha) * 2.*log(1.+alpha) / (alpha*(2.+alpha)) - 1.;
}

void compute_Dbeta(double &Dbeta,double& beta, double r){
  compute_beta(beta, r);
  Dbeta = 2.*(1.+alpha)/(r*r*alpha*(2.+alpha))*(beta-log(pow(alpha+1.,2)));
}

void compute_DDbeta(double &DDbeta, double &Dbeta, double &beta, double r){
  compute_Dbeta(Dbeta,beta,r);
  DDbeta = (1+alpha)*(r*r*Dbeta+beta*(1+alpha)/alpha/(2+alpha))/alpha/(2+alpha);
  DDbeta = DDbeta - (r*r*Dbeta+(beta+1)/(1+alpha))/(1+alpha);
  DDbeta = DDbeta + beta/pow(alpha*(2+alpha),2);
  DDbeta = 2*(DDbeta/r/r-Dbeta)/r; 
}

void compute_coefficients(double &A,double &B,double &C,double r,double t) {
  double beta,Dbeta,DDbeta;
  
  //compute_beta(beta, r);
  //compute_Dbeta(Dbeta,beta, r); //first derivative
  ComputeAlpha(r,t);
  compute_DDbeta(DDbeta, Dbeta, beta, r); //second derivative
  
  
  A = 1+(0.25+beta/2.+r*Dbeta/8.) *r*r;
  B = .5 + 0.75*beta+3./4.*r*Dbeta+0.125*r*r*DDbeta;
  C = .25 + 0.25/pow(1.+beta/4.*r*r,2);
  /**
     A = r*r*0.25;
     B=0.5;
     C=0.25;
  **/
}

double compute_Ez( double r, double t, double u){
  double Ez,beta,Dbeta;
    double Du=g(r,t,u);
  
  ComputeAlpha(r,t);
  compute_beta(beta, r);
  compute_Dbeta(Dbeta,beta, r); //first derivative
  
  Ez = 0.5*u*r * (1.+beta+0.5*r*Dbeta);

// 	Ez = 0.5*r*u;	
  return Ez;
}

