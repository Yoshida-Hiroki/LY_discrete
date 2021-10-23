#include <bits/stdc++.h>
#include <fstream>
using namespace std;


double pi = 3.141592;


// independent parameters
int dt = 1;
double omega = 2*pi/50;
double T = 1000;
// int N = 1;

// dependent parameters
int M = pi/(omega*dt);
double T_0 = 2*M*dt;
int N = T/T_0;
// int T = N*dt;

double r(double x){
  return pi/200*sin(x);
}

double phi_a(double x){
  return pi*3/4+pi/10*cos(x);
}

double phi_b(double x){
  return pi/4+pi/10*cos(x);
}

double a_L(double x){
  return 0.5*(1+r(x))*pow(sin(phi_a(x)/2),2.0);
}

double a_R(double x){
  return 0.5*(1+r(x))*pow(cos(phi_a(x)/2),2.0);
}

double b_L(double x){
  return 0.5*(1-r(x))*pow(sin(phi_b(x)/2),2.0);
}

double b_R(double x){
  return 0.5*(1-r(x))*pow(cos(phi_b(x)/2),2.0);
}

double a(double x){
  return a_R(x)+a_L(x);
}

double b(double x){
  return b_R(x)+b_L(x);
}

double z1(double x){
  return (-(pow(b(x)-a(x),2.0)/4+b_L(x)*a_L(x)+b_R(x)*a_R(x))+sqrt(pow((pow((b(x)-a(x)),2.0)/4+b_L(x)*a_L(x)+b_R(x)*a_R(x)),2.0)-4*a_R(x)*a_L(x)*b_R(x)*b_L(x)))/(2*a_R(x)*b_L(x));
}

double z2(double x){
  return (-(pow(b(x)-a(x),2.0)/4+b_L(x)*a_L(x)+b_R(x)*a_R(x))-sqrt(pow((pow((b(x)-a(x)),2.0)/4+b_L(x)*a_L(x)+b_R(x)*a_R(x)),2.0)-4*a_R(x)*a_L(x)*b_R(x)*b_L(x)))/(2*a_R(x)*b_L(x));
}

double sup_func(double x){
  return (double)(a(x)+b(x))*0.5;
}

double ave_sup_func(int i){
  double sum = 0.0;
  for(int i = 0; i < 2*M ; i++){
    double theta = 2*pi/(2*(double)M)*(double)i;
    sum += 1/(2*(double)M)*sup_func(theta);
  }
  return sum;
}

double ave_LHS(double x){
  double sum = 0.0;
  for(int i = 0; i < 2*M ; i++){
    double theta = 2*pi/(2*(double)M)*(double)i;
    sum += 1/(2*(double)M)*sup_func(theta)*sqrt(((x-z1(theta))*(x-z2(theta)))/(x*(1-z1(theta))*(1-z2(theta))));
  }
  return sum;
}

double ave_f(double x){
  double sum =0.0;
  for(int i = 0; i < 2*M ; i++){
    double theta = 2*pi/(2*(double)M)*(double)i;
    sum += 1/(2*(double)M)*sup_func(theta)*sqrt(-((x-z1(theta))*(x-z2(theta)))/(x*(1-z1(theta))*(1-z2(theta))));
  }
  return sum;
}

double ave_g(double x){
  double sum =0.0;
  for(int i = 0; i < 2*M ; i++){
    double theta = 2*pi/(2*(double)M)*(double)i;
    sum += 1/(2*(double)M)*sup_func(theta)*sqrt(-((x-z1(theta))*(x-z2(theta)))/(x*(1-z1(theta))*(1-z2(theta))))*(1/(x-z1(theta))+1/(x-z2(theta))-1/x);
  }
  return sum;
}

double R(int i){
  return ave_sup_func(0)*dt/(1-ave_sup_func(0)*dt);
}

double rho1(double x){
  return R(0)*ave_g(x)/ave_sup_func(0)/(1+pow(R(0),2.0)*pow(ave_f(x)/ave_sup_func(0),2.0));
}

double rho2(double x){
  return -R(0)*ave_g(x)/ave_sup_func(0)/(1+pow(R(0),2.0)*pow(ave_f(x)/ave_sup_func(0),2.0));
}


int main(){
  double zmin = -25;
  double zmax = 0;
  int partnum = 10;

  string path = "C:/Users/hyoshida/Desktop/";
  string file = "den.dat";
  string filename = path + file;
  ofstream writing_file;
  writing_file.open(filename, ios::out);

  for(int i = 0;i< partnum;i++){
    double z = zmin + (zmax-zmin)/partnum*(double)i;

    writing_file << z << " " << rho1(z) << " " << rho2(z) << endl;
  }
}
