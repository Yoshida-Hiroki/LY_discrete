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
string ver = "_N1_1";

vector<double> z_disc(2);

double xU_1min,xU_1max,dxU_1;
double xU_2min,xU_2max,dxU_2;
int partnum1 = 1.0e+04;
int partnum2 = 1.0e+04;
vector<double> RhoU(partnum1),RhoU2(partnum2);
vector<double> JU_dat;

double root(double x){
  return sqrt((-(x-z_disc[0])*(x-z_disc[1]))/(x*(1-z_disc[0])*(1-z_disc[1])));
}

double rho(double x){
  return 1/(2*pi)*root(x)/(1+pow(root(x),2.0))*((1/(x-z_disc[0])+1/(x-z_disc[1]))-1.0/x);
}


int main(){

  ifstream fin("C:/Users/hyoshida/Desktop/floquetic/zero_"+date+ver+".dat");
  fin >> z_disc[0] >> z_disc[1];

  xU_1min = z_disc[0]-2;
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
  ////////////// rho plot //////////////////
  string filename2 = path + "rho_"+date+ver + ext;
  ofstream writing_file2;
  writing_file2.open(filename2, ios::out);

  for(int l = 0 ; l < partnum1+1;l++){
    double x = xU_1min + (double)dxU_1*l;
    writing_file2 << x << " "<< RhoU[l] << endl;
  }
  for(int l = 0 ; l < partnum2+1;l++){
    double x = xU_2min + (double)dxU_2*l;
    writing_file2 << x << " "<< RhoU2[l] << endl;
  }

  end = clock();
  cout << "rho write : " << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;

  Beep(750,200);
}
