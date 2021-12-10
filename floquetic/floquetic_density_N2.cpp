#include <bits/stdc++.h>
#include <fstream>
#include<iostream>
#include <time.h>
#include <Windows.h>
using namespace std;
#define GNUPLOT_PATH "C:/PROGRA~1/gnuplot/bin/gnuplot.exe"

double pi = 3.141592;

double dt = 1;

string date = "211210";
string ver = "_N2_1";

vector<double> z_disc(4);
vector<double> nume(3);
vector<double> deno(3);
vector<double> z_1(2);
vector<double> z_2(2);

double xU_1min,xU_1max,dxU_1;
double xU_2min,xU_2max,dxU_2;
int partnum1 = 1.0e+06;
int partnum2 = 1.0e+06;
vector<double> RhoU(partnum1),RhoU2(partnum2);
vector<double> JU_dat;

double x1_1min,x1_1max,dx1_1;
double x1_2min,x1_2max,dx1_2;
vector<double> Rho1_1(partnum1),Rho1_2(partnum2);
vector<double> J1_dat;

double x2_1min,x2_1max,dx2_1;
double x2_2min,x2_2max,dx2_2;
vector<double> Rho2_1(partnum1),Rho2_2(partnum2);
vector<double> J2_dat;


double Trace(double x){
  return deno[0]*x+deno[1]+deno[2]/x;
}

double R(double x){
  return abs(Trace(1)-2)/Trace(x);
}

double root(double x){
  return sqrt(-(x-z_disc[0])*(x-z_disc[1])*(x-z_disc[2])*(x-z_disc[3])/(pow(x,2.0)*(1-z_disc[0])*(1-z_disc[1])*(1-z_disc[2])*(1-z_disc[3])));
}

double rho(double x){
  return 1/pi*(abs(Trace(1)-2)*root(x))/(pow(Trace(x),2.0)+pow((Trace(1)-2)*root(x),2.0))*(((1/(x-z_disc[0])+1/(x-z_disc[1])+1/(x-z_disc[2])+1/(x-z_disc[3]))-2.0/x)*Trace(x)-2*(nume[0]*pow(x,2.0)+nume[1]*x+nume[2])/(pow(x,2.0)));
}

double J(double z){
  double integ = 0;
  for(int i = 0 ; i< partnum1;i++){
    double xU_1 = xU_1min + (double)dxU_1*i;
    double temp = (double)dxU_1*RhoU[i]*z/(z-xU_1);
    integ += (temp==temp) ? temp : 0;
  }
  for(int i = 0 ; i< partnum2;i++){
    double xU_2 = xU_2min + (double)dxU_2*i;
    double temp = (double)dxU_2*RhoU2[i]*z/(z-xU_2);
    integ += (temp==temp) ? temp : 0;
  }
  return 0.5*integ -1.0;
}

double phi(double z,int j){
  double integ = 0;
  for(int i = 0 ; i< partnum1;i++){
    double xU_1 = xU_1min + (double)dxU_1*i;
    double temp = dxU_1*RhoU[i]*(log((z-xU_1)/(1-xU_1)));
    integ += (temp == temp) ? temp:0;
  }
  for(int i = 0 ; i< partnum2;i++){
    double xU_2 = xU_2min + (double)dxU_2*i;
    double temp = dxU_2*RhoU2[i]*(log((z-xU_2)/(1-xU_2)));
    integ += (temp == temp) ? temp:0;
  }
  return 0.5*integ-JU_dat[j]*log(z)-log(z);
}

///////////// adiabatic approximation ////////////////

double root1(double x){
  return sqrt(-(x-z_1[0])*(x-z_1[1])/(x*(1-z_1[0])*(1-z_1[1])));
}

double rho1(double x){
  return 1/pi*root1(x)/(1+pow(root1(x),2.0))*(1/(x-z_1[0])+1/(x-z_1[1])-1.0/x);
}

double J1(double z){
  double integ = 0;
  for(int i = 0 ; i< partnum1;i++){
    double x = x1_1min + (double)dx1_1*i;
    double temp = (double)dx1_1*Rho1_1[i]*z/(z-x);
    integ += (temp==temp) ? temp : 0;
  }
  for(int i = 0 ; i< partnum2;i++){
    double x = x1_2min + (double)dx1_2*i;
    double temp = (double)dx1_2*Rho1_2[i]*z/(z-x);
    integ += (temp==temp) ? temp : 0;
  }
  return 0.5*integ -0.5;
}

double phi1(double z,int j){
  double integ = 0;
  for(int i = 0 ; i< partnum1;i++){
    double x = x1_1min + (double)dx1_1*i;
    double temp = dx2_1*Rho1_1[i]*(log((z-x)/(1-x)));
    integ += (temp == temp) ? temp:0;
  }
  for(int i = 0 ; i< partnum2;i++){
    double x = x1_2min + (double)dx1_2*i;
    double temp = dx1_2*Rho1_2[i]*(log((z-x)/(1-x)));
    integ += (temp == temp) ? temp:0;
  }
  return 0.5*integ-J1_dat[j]*log(z)-0.5*log(z);
}



double root2(double x){
  return sqrt(-(x-z_2[0])*(x-z_2[1])/(x*(1-z_2[0])*(1-z_2[1])));
}

double rho2(double x){
  return 1/pi*root2(x)/(1+pow(root2(x),2.0))*(1/(x-z_2[0])+1/(x-z_2[1])-1.0/x);
}

double J2(double z){
  double integ = 0;
  for(int i = 0 ; i< partnum1;i++){
    double x = x2_1min + (double)dx2_1*i;
    double temp = (double)dx2_1*Rho2_1[i]*z/(z-x);
    integ += (temp==temp) ? temp : 0;
  }
  for(int i = 0 ; i< partnum2;i++){
    double x = x2_2min + (double)dx2_2*i;
    double temp = (double)dx2_2*Rho2_2[i]*z/(z-x);
    integ += (temp==temp) ? temp : 0;
  }
  return 0.5*integ -0.5;
}

double phi2(double z,int j){
  double integ = 0;
  for(int i = 0 ; i< partnum1;i++){
    double x = x2_1min + (double)dx2_1*i;
    double temp = dx2_1*Rho2_1[i]*(log((z-x)/(1-x)));
    integ += (temp == temp) ? temp:0;
  }
  for(int i = 0 ; i< partnum2;i++){
    double x = x2_2min + (double)dx2_2*i;
    double temp = dx2_2*Rho2_2[i]*(log((z-x)/(1-x)));
    integ += (temp == temp) ? temp:0;
  }
  return 0.5*integ-J2_dat[j]*log(z)-0.5*log(z);
}




int main(){

  ifstream fin("C:/Users/hyoshida/Desktop/floquetic/zero_"+date+ver+".dat");
  fin >> z_disc[0] >> z_disc[1] >> z_disc[2] >> z_disc[3] >> nume[0] >> nume[1] >> nume[2] >> deno[0] >> deno[1] >> deno[2] >> z_1[0] >> z_1[1]>> z_2[0] >> z_2[1];

  xU_1min = z_disc[0];
  xU_1max = z_disc[1];
  dxU_1 = (double)(xU_1max-xU_1min)/partnum1;
  xU_2min = z_disc[2];
  xU_2max = z_disc[3];
  dxU_2 = (double)(xU_2max-xU_2min)/partnum2;

  x1_1min = -100;
  x1_1max = z_1[0];
  dx1_1 = (double)(x1_1max-x1_1min)/partnum1;
  x1_2min = z_1[1];
  x1_2max = 0;
  dx1_2 = (double)(x1_2max-x1_2min)/partnum2;

  x2_1min = -100;
  x2_1max = z_2[0];
  dx2_1 = (double)(x2_1max-x2_1min)/partnum1;
  x2_2min = z_2[1];
  x2_2max = 0;
  dx2_2 = (double)(x2_2max-x2_2min)/partnum2;

  clock_t start = clock();
  for(int j = 0 ;j<partnum1+1;j++){
    double xU_1 = xU_1min + (double)dxU_1*j;
    RhoU[j] = abs(rho(xU_1));
    double x1_1 = x1_1min + (double)dx1_1*j;
    Rho1_1[j] = abs(rho1(x1_1));
    double x2_1 = x2_1min + (double)dx2_1*j;
    Rho2_1[j] = abs(rho2(x2_1));
  }
  for(int j = 0 ;j<partnum2+1;j++){
    double xU_2 = xU_2min + (double)dxU_2*j;
    RhoU2[j] = abs(rho(xU_2));
    double x1_2 = x1_2min + (double)dx1_2*j;
    Rho1_2[j] = abs(rho1(x1_2));
    double x2_2 = x2_2min + (double)dx2_2*j;
    Rho2_2[j] = abs(rho2(x2_2));
  }
  clock_t end = clock();
  cout << "Rho : " << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;

  // long long s1=0;
  // for(int i = 0 ; i < partnum1;i++){
  //   s1 += RhoU[i];
  // }
  // cout << dxU_1 << endl;
  // cout << "sum rho : " << s1*dxU_1 << endl;
  //
  // long long s2=0;
  // for(int i = 0 ; i < partnum2;i++){
  //   s2 += RhoU2[i];
  // }
  // cout << dxU_2 << endl;
  // cout << "sum rho2 : " << dxU_2*s2 << endl;
  // cout << endl;

  string path = "C:/Users/hyoshida/Desktop/floquetic/";
  string ext = ".dat";

  start = clock();
  // //////////// J-phi(J) plot ////////////////////
  string filename = path + "phi_"+date+ver + ext;
  ofstream writing_file;
  writing_file.open(filename, ios::out);

  double chi_min = -4;
  double chi_max = 5;
  int chi_part = 100;

  vector<double> Jad_dat(chi_part);

  for(int j = 0 ; j < chi_part;j++){
    double chi = chi_min+(double)(chi_max-chi_min)/chi_part*j;
    JU_dat.push_back(J(exp(chi)));
    J1_dat.push_back(J1(exp(chi)));
    J2_dat.push_back(J2(exp(chi)));
    Jad_dat[j] = J1_dat[j]+J2_dat[j];
  }
  end = clock();
  cout << "J : " << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;

  start = clock();
  for(int k = 0 ; k < chi_part ; k++){
    double chi = chi_min+(double)(chi_max-chi_min)/chi_part*k;
    writing_file << exp(chi) << " "<< JU_dat[k] << " " << phi(exp(chi),k) << " "<< Jad_dat[k] << " " << phi1(exp(chi),k)+phi2(exp(chi),k) << endl;
  }
  end = clock();
  cout << "phi : " << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;

  // start = clock();
  // ////////////// rho plot //////////////////
  // string filename2 = path + "rho_"+date+ver + ext;
  // ofstream writing_file2;
  // writing_file2.open(filename2, ios::out);
  //
  // for(int l = 0 ; l < partnum1;l++){
  //   double xU_1 = xU_1min + (double)dxU_1*l;
  //   writing_file2 << xU_1 << " "<< RhoU[l] << endl;
  // }
  // for(int l = 0 ; l < partnum2;l++){
  //   double xU_2 = xU_2min + (double)dxU_2*l;
  //   writing_file2 << xU_2 << " "<< RhoU2[l] << endl;
  // }
  //
  // end = clock();
  // cout << "rho write : " << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;

  Beep(750,200);
}
