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
string ver = "_N4_2";

//////////////////////////////////////
vector<double> z_disc(8);
vector<double> nume(5);
vector<double> deno(5);
vector<double> z_1(2),z_2(2),z_3(2),z_4(2);

double x_1_1,x_1_2,x_1_3,x_1_4,dx_1_1,dx_1_2,dx_1_3;
double x_2_1,x_2_2,x_2_3,x_2_4,dx_2_1,dx_2_2,dx_2_3;
double x_3_1,x_3_2,x_3_3,x_3_4,dx_3_1,dx_3_2,dx_3_3;
double x_4_1,x_4_2,x_4_3,x_4_4,dx_4_1,dx_4_2,dx_4_3;
long long partnum1 = 1.0e+03,partnum2 = 1.0e+04;
vector<double> RhoU_1(partnum1+2*partnum2),RhoU_2(partnum1+2*partnum2),RhoU_3(partnum1+2*partnum2),RhoU_4(partnum1+2*partnum2);

double Trace(double x){
  return deno[0]*pow(x,2.0)+deno[1]*x+deno[2]+deno[3]/x+deno[4]/pow(x,2.0);
}

double R(double x){
  return abs(Trace(1)-2)/Trace(x);
}

double root(double x){
  return sqrt(-(x-z_disc[0])*(x-z_disc[1])*(x-z_disc[2])*(x-z_disc[3])*(x-z_disc[4])*(x-z_disc[5])*(x-z_disc[6])*(x-z_disc[7])/(pow(x,4.0)*(1-z_disc[0])*(1-z_disc[1])*(1-z_disc[2])*(1-z_disc[3])*(1-z_disc[4])*(1-z_disc[5])*(1-z_disc[6])*(1-z_disc[7])));
}

double rho(double x){
  return 1/pi*(R(x)*root(x))/(1+pow(R(x)*root(x),2.0))*(1/(x-z_disc[0])+1/(x-z_disc[1])+1/(x-z_disc[2])+1/(x-z_disc[3])+1/(x-z_disc[4])+1/(x-z_disc[5])+1/(x-z_disc[6])+1/(x-z_disc[7])-4.0/x-2*(nume[0]*pow(x,4.0)+nume[1]*pow(x,3.0)+nume[2]*pow(x,2.0)+nume[3]*x+nume[4])/(Trace(x)*pow(x,3.0)));
}

int main(){
  ifstream fin("C:/Users/hyoshida/Desktop/floquetic/zero_"+date+ver+".dat");
  fin >> z_disc[0] >> z_disc[1] >> z_disc[2] >> z_disc[3] >>z_disc[4] >>z_disc[5] >>z_disc[6] >>z_disc[7] >> nume[0] >> nume[1] >> nume[2] >> nume[3] >>nume[4] >> deno[0] >> deno[1] >> deno[2] >> deno[3] >>deno[4] >> z_1[0] >> z_1[1] >> z_2[0] >> z_2[1] >> z_3[0] >> z_3[1] >> z_4[0] >> z_4[1];

  x_1_1 = z_disc[0];
  x_1_4 = z_disc[1];
  x_1_2 = z_disc[0]+(x_1_4-x_1_1)*1.0e-09;
  x_1_3 = z_disc[1]-(x_1_4-x_1_1)*1.0e-07;
  dx_1_1 = (double)(x_1_2-x_1_1)/partnum2;
  dx_1_2 = (double)(x_1_3-x_1_2)/partnum1;
  dx_1_3 = (double)(x_1_4-x_1_3)/partnum2;
  x_2_1 = z_disc[2];
  x_2_4 = z_disc[3];
  x_2_2 = z_disc[2]+(x_2_4-x_2_1)*1.0e-03;
  x_2_3 = z_disc[3]-(x_2_4-x_2_1)*1.0e-03;
  dx_2_1 = (double)(x_2_2-x_2_1)/partnum2;
  dx_2_2 = (double)(x_2_3-x_2_2)/partnum1;
  dx_2_3 = (double)(x_2_4-x_2_3)/partnum2;
  x_3_1 = z_disc[4];
  x_3_4 = z_disc[5];
  x_3_2 = z_disc[4]+(x_3_4-x_3_1)*1.0e-03;
  x_3_3 = z_disc[5]-(x_3_4-x_3_1)*1.0e-03;
  dx_3_1 = (double)(x_3_2-x_3_1)/partnum2;
  dx_3_2 = (double)(x_3_3-x_3_2)/partnum1;
  dx_3_3 = (double)(x_3_4-x_3_3)/partnum2;
  x_4_1 = z_disc[6];
  x_4_4 = z_disc[7];
  x_4_2 = z_disc[6]+(x_4_4-x_4_1)*1.0e-03;
  x_4_3 = z_disc[7]-(x_4_4-x_4_1)*1.0e-03;
  dx_4_1 = (double)(x_4_2-x_4_1)/partnum2;
  dx_4_2 = (double)(x_4_3-x_4_2)/partnum1;
  dx_4_3 = (double)(x_4_4-x_4_3)/partnum2;

  clock_t start = clock();

  for(int j = 0 ;j<partnum2;j++){
    double x = x_1_1 + (double)dx_1_1*j;
    RhoU_1[j] = abs(rho(x));
  }
  for(int j = 0 ;j<partnum1;j++){
    double x = x_1_2 + (double)dx_1_2*j;
    RhoU_1[partnum2+j] = abs(rho(x));
  }
  for(int j = 0 ;j<partnum2;j++){
    double x = x_1_3 + (double)dx_1_3*j;
    RhoU_1[partnum2+partnum1+j] = abs(rho(x));
  }

  for(int j = 0 ;j<partnum2;j++){
    double x = x_2_1 + (double)dx_2_1*j;
    RhoU_2[j] = abs(rho(x));
  }
  for(int j = 0 ;j<partnum1;j++){
    double x = x_2_2 + (double)dx_2_2*j;
    RhoU_2[partnum2+j] = abs(rho(x));
  }
  for(int j = 0 ;j<partnum2;j++){
    double x = x_2_3 + (double)dx_2_3*j;
    RhoU_2[partnum2+partnum1+j] = abs(rho(x));
  }

  for(int j = 0 ;j<partnum2;j++){
    double x = x_3_1 + (double)dx_3_1*j;
    RhoU_3[j] = abs(rho(x));
  }
  for(int j = 0 ;j<partnum1;j++){
    double x = x_3_2 + (double)dx_3_2*j;
    RhoU_3[partnum2+j] = abs(rho(x));
  }
  for(int j = 0 ;j<partnum2;j++){
    double x = x_3_3 + (double)dx_3_3*j;
    RhoU_3[partnum2+partnum1+j] = abs(rho(x));
  }

  for(int j = 0 ;j<partnum2;j++){
    double x = x_4_1 + (double)dx_4_1*j;
    RhoU_4[j] = abs(rho(x));
  }
  for(int j = 0 ;j<partnum1;j++){
    double x = x_4_2 + (double)dx_4_2*j;
    RhoU_4[partnum2+j] = abs(rho(x));
  }
  for(int j = 0 ;j<partnum2;j++){
    double x = x_4_3 + (double)dx_4_3*j;
    RhoU_4[partnum2+partnum1+j] = abs(rho(x));
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
    writing_file2 << x << " "<< RhoU_1[l] << endl;
  }
  for(int l = 0 ; l < partnum1;l++){
    double x = x_1_2 + (double)dx_1_2*l;
    writing_file2 << x << " "<< RhoU_1[partnum2+l] << endl;
  }
  for(int l = 0 ; l < partnum2;l++){
    double x = x_1_3 + (double)dx_1_3*l;
    writing_file2 << x << " "<< RhoU_1[partnum2+partnum1+l] << endl;
  }

  for(int l = 0 ; l < partnum2;l++){
    double x = x_2_1 + (double)dx_2_1*l;
    writing_file2 << x << " "<< RhoU_2[l] << endl;
  }
  for(int l = 0 ; l < partnum1;l++){
    double x = x_2_2 + (double)dx_2_2*l;
    writing_file2 << x << " "<< RhoU_2[partnum2+l] << endl;
  }
  for(int l = 0 ; l < partnum2;l++){
    double x = x_2_3 + (double)dx_2_3*l;
    writing_file2 << x << " "<< RhoU_2[partnum2+partnum1+l] << endl;
  }

  for(int l = 0 ; l < partnum2;l++){
    double x = x_3_1 + (double)dx_3_1*l;
    writing_file2 << x << " "<< RhoU_3[l] << endl;
  }
  for(int l = 0 ; l < partnum1;l++){
    double x = x_3_2 + (double)dx_3_2*l;
    writing_file2 << x << " "<< RhoU_3[partnum2+l] << endl;
  }
  for(int l = 0 ; l < partnum2;l++){
    double x = x_3_3 + (double)dx_3_3*l;
    writing_file2 << x << " "<< RhoU_3[partnum2+partnum1+l] << endl;
  }

  for(int l = 0 ; l < partnum2;l++){
    double x = x_4_1 + (double)dx_4_1*l;
    writing_file2 << x << " "<< RhoU_4[l] << endl;
  }
  for(int l = 0 ; l < partnum1;l++){
    double x = x_4_2 + (double)dx_4_2*l;
    writing_file2 << x << " "<< RhoU_4[partnum2+l] << endl;
  }
  for(int l = 0 ; l < partnum2;l++){
    double x = x_4_3 + (double)dx_4_3*l;
    writing_file2 << x << " "<< RhoU_4[partnum2+partnum1+l] << endl;
  }

  end = clock();
  cout << "rho write : " << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;

  Beep(750,200);
}
