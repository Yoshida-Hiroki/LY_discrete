#include <bits/stdc++.h>
#include <fstream>
#include<iostream>
#include <time.h>
#include <Windows.h>
using namespace std;
#define GNUPLOT_PATH "C:/PROGRA~1/gnuplot/bin/gnuplot.exe"

double pi = 3.141592;

double dt = 1;

string date = "211225";
string ver = "_N1_2";

vector<double> z_disc(2);

double xU_1min,xU_1max,dxU_1;
double xU_2min,xU_2max,dxU_2;
int partnum1 = 1.0e+07;
int partnum2 = 1.0e+07;
vector<double> RhoU(partnum1),RhoU2(partnum2);
vector<double> JU_dat;

double root(double x){
  return sqrt((-(x-z_disc[0])*(x-z_disc[1]))/(x*(1-z_disc[0])*(1-z_disc[1])));
}

double rho(double x){
  return 1/(2*pi)*root(x)/(1+pow(root(x),2.0))*((1/(x-z_disc[0])+1/(x-z_disc[1]))-1.0/x);
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
  return integ -0.5;
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
  return integ-0.5*log(z)-JU_dat[j]*log(z);
}

int main(){

  ifstream fin("C:/Users/hyoshida/Desktop/floquetic/zero_"+date+ver+".dat");
  fin >> z_disc[0] >> z_disc[1];

  xU_1min = z_disc[0]-100;
  xU_1max = z_disc[0]-1.0e-10;
  dxU_1 = (double)(xU_1max-xU_1min)/partnum1;
  xU_2min = z_disc[1];
  xU_2max = 0;
  dxU_2 = (double)(xU_2max-xU_2min)/partnum2;

  clock_t start = clock();
  for(int j = 0 ;j<partnum1+1;j++){
    double xU_1 = xU_1min + (double)dxU_1*j;
    RhoU[j] = abs(rho(xU_1));
  }
  for(int j = 0 ;j<partnum2+1;j++){
    double xU_2 = xU_2min + (double)dxU_2*j;
    RhoU2[j] = abs(rho(xU_2));
  }
  clock_t end = clock();
  cout << "Rho : " << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;

  string path = "C:/Users/hyoshida/Desktop/floquetic/";
  string ext = ".dat";

  start = clock();
  //////////// J-phi(J) plot ////////////////////
  string filename = path + "phi_"+date+ver + ext;
  ofstream writing_file;
  writing_file.open(filename, ios::out);

  double chi_min = -2.5;
  double chi_max = 2;
  int chi_part = 500;

  vector<double> Jad_dat(chi_part);

  for(int j = 0 ; j < chi_part;j++){
    double chi = chi_min+(double)(chi_max-chi_min)/chi_part*j;
    JU_dat.push_back(J(exp(chi)));
  }
  end = clock();
  cout << "J : " << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;

  start = clock();
  for(int k = 0 ; k < chi_part ; k++){
    double chi = chi_min+(double)(chi_max-chi_min)/chi_part*k;
    writing_file << exp(chi) << " "<< JU_dat[k] << " " << phi(exp(chi),k) << endl;
  }
  end = clock();
  cout << "phi : " << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;

  Beep(750,200);
}
