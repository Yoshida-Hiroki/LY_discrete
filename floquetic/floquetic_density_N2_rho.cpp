#include <bits/stdc++.h>
#include <fstream>
#include<iostream>
#include <time.h>
#include <Windows.h>
using namespace std;
#define GNUPLOT_PATH "C:/PROGRA~1/gnuplot/bin/gnuplot.exe"

double pi = 3.141592;

double dt = 1;

string date = "211213";
string ver = "_N2_2";

vector<double> z_disc(4);
vector<double> nume(3);
vector<double> deno(3);
vector<double> z_1(2);
vector<double> z_2(2);

double x_1_1,x_1_2,x_1_3,x_1_4,dx_1_1,dx_1_2,dx_1_3;
double x_2_1,x_2_2,x_2_3,x_2_4,dx_2_1,dx_2_2,dx_2_3;
int partnum1 = 1.0e+03;
int partnum2 = 1.0e+04;
vector<double> RhoU(partnum1+2*partnum2),RhoU2(partnum1+2*partnum2);

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

int main(){

  ifstream fin("C:/Users/hyoshida/Desktop/floquetic/zero_"+date+ver+".dat");
  fin >> z_disc[0] >> z_disc[1] >> z_disc[2] >> z_disc[3] >> nume[0] >> nume[1] >> nume[2] >> deno[0] >> deno[1] >> deno[2] >> z_1[0] >> z_1[1]>> z_2[0] >> z_2[1];

  x_1_1 = z_disc[0];
  x_1_4 = z_disc[1];
  x_1_2 = z_disc[0]+(x_1_4-x_1_1)*1.0e-03;
  x_1_3 = z_disc[1]-(x_1_4-x_1_1)*1.0e-03;
  dx_1_1 = (double)(x_1_2-x_1_1)/partnum2;
  dx_1_2 = (double)(x_1_3-x_1_2)/partnum1;
  dx_1_3 = (double)(x_1_4-x_1_3)/partnum2;
  x_2_1 = z_disc[2];
  x_2_4 = z_disc[3];
  x_2_2 = z_disc[2]+(x_2_4-x_2_1)*3.0e-02;
  x_2_3 = z_disc[3]-(x_2_4-x_2_1)*4.0e-02;
  dx_2_1 = (double)(x_2_2-x_2_1)/partnum2;
  dx_2_2 = (double)(x_2_3-x_2_2)/partnum1;
  dx_2_3 = (double)(x_2_4-x_2_3)/partnum2;

  clock_t start = clock();
  for(int j = 0 ;j<partnum2;j++){
    double x = x_1_1 + (double)dx_1_1*j;
    RhoU[j] = abs(rho(x));
  }
  for(int j = 0 ;j<partnum1;j++){
    double x = x_1_2 + (double)dx_1_2*j;
    RhoU[partnum2+j] = abs(rho(x));
  }
  for(int j = 0;j<partnum2;j++){
    double x = x_1_3 + (double)dx_1_3*j;
    RhoU[partnum1+partnum2+j] = abs(rho(x));
  }

  for(int j = 0 ;j<partnum2;j++){
    double x = x_2_1 + (double)dx_2_1*j;
    RhoU2[j] = abs(rho(x));
  }
  for(int j = 0 ;j<partnum1;j++){
    double x = x_2_2 + (double)dx_2_2*j;
    RhoU2[partnum2+j] = abs(rho(x));
  }
  for(int j = 0;j<partnum2;j++){
    double x = x_2_3 + (double)dx_2_3*j;
    RhoU2[partnum1+partnum2+j] = abs(rho(x));
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

  for(int l = 0 ; l < partnum2;l++){
    double x = x_1_1 + (double)dx_1_1*l;
    writing_file2 << x << " "<< RhoU[l] << endl;
  }
  for(int l = 0 ; l < partnum1;l++){
    double x = x_1_2 + (double)dx_1_2*l;
    writing_file2 << x << " "<< RhoU[partnum2+l] << endl;
  }
  for(int l = 0 ; l < partnum2;l++){
    double x = x_1_3 + (double)dx_1_3*l;
    writing_file2 << x << " "<< RhoU[partnum1+partnum2+l] << endl;
  }

  for(int l = 0 ; l < partnum2;l++){
    double x = x_2_1 + (double)dx_2_1*l;
    writing_file2 << x << " "<< RhoU2[l] << endl;
  }
  for(int l = 0 ; l < partnum1;l++){
    double x = x_2_2 + (double)dx_2_2*l;
    writing_file2 << x << " "<< RhoU2[partnum2+l] << endl;
  }
  for(int l = 0 ; l < partnum2;l++){
    double x = x_2_3 + (double)dx_2_3*l;
    writing_file2 << x << " "<< RhoU2[partnum1+partnum2+l] << endl;
  }

  end = clock();
  cout << "rho write : " << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;

  Beep(750,200);
}
