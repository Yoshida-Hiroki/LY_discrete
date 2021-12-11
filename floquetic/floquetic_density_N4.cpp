#include <bits/stdc++.h>
#include <fstream>
#include<iostream>
#include <time.h>
#include <Windows.h>
using namespace std;
#define GNUPLOT_PATH "C:/PROGRA~1/gnuplot/bin/gnuplot.exe"

double pi = 3.141592;

double dt = 1;

string date = "211211";
string ver = "_N4_2";

//////////////////////////////////////
vector<double> z_disc(8);
vector<double> nume(5);
vector<double> deno(5);
vector<double> z_1(2),z_2(2),z_3(2),z_4(2);

double xU_1min,xU_1max,dxU_1,xU_2min,xU_2max,dxU_2,xU_3min,xU_3max,dxU_3,xU_4min,xU_4max,dxU_4;
long long partnum1 = 1.0e+03,partnum2 = 1.0e+04,partnum3 = 1.0e+03,partnum4 = 1.0e+03;
vector<double> RhoU_1(partnum1),RhoU_2(partnum2),RhoU_3(partnum3),RhoU_4(partnum4);
vector<double> JU_dat;

int Partnum1=1.0e+01,Partnum2 = 1.0e+01;
double x1_1min,x1_1max,dx1_1;
double x1_2min,x1_2max,dx1_2;
vector<double> Rho1_1(Partnum1),Rho1_2(Partnum2);
vector<double> J1_dat;

double x2_1min,x2_1max,dx2_1;
double x2_2min,x2_2max,dx2_2;
vector<double> Rho2_1(Partnum1),Rho2_2(Partnum2);
vector<double> J2_dat;

double x3_1min,x3_1max,dx3_1;
double x3_2min,x3_2max,dx3_2;
vector<double> Rho3_1(Partnum1),Rho3_2(Partnum2);
vector<double> J3_dat;

double x4_1min,x4_1max,dx4_1;
double x4_2min,x4_2max,dx4_2;
vector<double> Rho4_1(Partnum1),Rho4_2(Partnum2);
vector<double> J4_dat;

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


double J(double z){
  double integ = 0;
  for(int i = 0 ; i< partnum1;i++){
    double xU_1 = xU_1min + (double)dxU_1*i;
    double temp = (double)dxU_1*RhoU_1[i]*z/(z-xU_1);
    integ += (temp==temp) ? temp : 0;
  }
  for(int i = 0 ; i< partnum2;i++){
    double xU_2 = xU_2min + (double)dxU_2*i;
    double temp = (double)dxU_2*RhoU_2[i]*z/(z-xU_2);
    integ += (temp==temp) ? temp : 0;
  }
  for(int i = 0 ; i< partnum3;i++){
    double xU_3 = xU_3min + (double)dxU_3*i;
    double temp = (double)dxU_3*RhoU_3[i]*z/(z-xU_3);
    integ += (temp==temp) ? temp : 0;
  }
  for(int i = 0 ; i< partnum4;i++){
    double xU_4 = xU_4min + (double)dxU_4*i;
    double temp = (double)dxU_4*RhoU_4[i]*z/(z-xU_4);
    integ += (temp==temp) ? temp : 0;
  }
  return 0.5*integ -2.0;
}

double phi(double z,int j){
  double integ = 0;
  for(int i = 0 ; i< partnum1;i++){
    double xU_1 = xU_1min + (double)dxU_1*i;
    double temp = dxU_1*RhoU_1[i]*(log((z-xU_1)/(1-xU_1)));
    integ += (temp == temp) ? temp:0;
  }
  for(int i = 0 ; i< partnum2;i++){
    double xU_2 = xU_2min + (double)dxU_2*i;
    double temp = dxU_2*RhoU_2[i]*(log((z-xU_2)/(1-xU_2)));
    integ += (temp == temp) ? temp:0;
  }
  for(int i = 0 ; i< partnum3;i++){
    double xU_3 = xU_3min + (double)dxU_3*i;
    double temp = dxU_3*RhoU_3[i]*(log((z-xU_3)/(1-xU_3)));
    integ += (temp == temp) ? temp:0;
  }
  for(int i = 0 ; i< partnum4;i++){
    double xU_4 = xU_4min + (double)dxU_4*i;
    double temp = dxU_4*RhoU_4[i]*(log((z-xU_4)/(1-xU_4)));
    integ += (temp == temp) ? temp:0;
  }
  return 0.5*integ-JU_dat[j]*log(z)-2.0*log(z);
}

///////////// adiabatic approximation /////////////////////////////////////////////////////////////////////////

double root1(double x){
  return sqrt(-(x-z_1[0])*(x-z_1[1])/(x*(1-z_1[0])*(1-z_1[1])));
}

double rho1(double x){
  return 1/pi*root1(x)/(1+pow(root1(x),2.0))*(1/(x-z_1[0])+1/(x-z_1[1])-1.0/x);
}

double J1(double z){
  double integ = 0;
  for(int i = 0 ; i< Partnum1;i++){
    double x = x1_1min + (double)dx1_1*i;
    double temp = (double)dx1_1*Rho1_1[i]*z/(z-x);
    integ += (temp==temp) ? temp : 0;
  }
  for(int i = 0 ; i< Partnum2;i++){
    double x = x1_2min + (double)dx1_2*i;
    double temp = (double)dx1_2*Rho1_2[i]*z/(z-x);
    integ += (temp==temp) ? temp : 0;
  }
  return 0.5*integ -0.5;
}

double phi1(double z,int j){
  double integ = 0;
  for(int i = 0 ; i< Partnum1;i++){
    double x = x1_1min + (double)dx1_1*i;
    double temp = dx1_1*Rho1_1[i]*(log((z-x)/(1-x)));
    integ += (temp == temp) ? temp:0;
  }
  for(int i = 0 ; i< Partnum2;i++){
    double x = x1_2min + (double)dx1_2*i;
    double temp = dx1_2*Rho1_2[i]*(log((z-x)/(1-x)));
    integ += (temp == temp) ? temp:0;
  }
  return 0.5*integ-J1_dat[j]*log(z)-0.5*log(z);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////

double root2(double x){
  return sqrt(-(x-z_2[0])*(x-z_2[1])/(x*(1-z_2[0])*(1-z_2[1])));
}

double rho2(double x){
  return 1/pi*root2(x)/(1+pow(root2(x),2.0))*(1/(x-z_2[0])+1/(x-z_2[1])-1.0/x);
}

double J2(double z){
  double integ = 0;
  for(int i = 0 ; i< Partnum1;i++){
    double x = x2_1min + (double)dx2_1*i;
    double temp = (double)dx2_1*Rho2_1[i]*z/(z-x);
    integ += (temp==temp) ? temp : 0;
  }
  for(int i = 0 ; i< Partnum2;i++){
    double x = x2_2min + (double)dx2_2*i;
    double temp = (double)dx2_2*Rho2_2[i]*z/(z-x);
    integ += (temp==temp) ? temp : 0;
  }
  return 0.5*integ -0.5;
}

double phi2(double z,int j){
  double integ = 0;
  for(int i = 0 ; i< Partnum1;i++){
    double x = x2_1min + (double)dx2_1*i;
    double temp = dx2_1*Rho2_1[i]*(log((z-x)/(1-x)));
    integ += (temp == temp) ? temp:0;
  }
  for(int i = 0 ; i< Partnum2;i++){
    double x = x2_2min + (double)dx2_2*i;
    double temp = dx2_2*Rho2_2[i]*(log((z-x)/(1-x)));
    integ += (temp == temp) ? temp:0;
  }
  return 0.5*integ-J2_dat[j]*log(z)-0.5*log(z);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////

double root3(double x){
  return sqrt(-(x-z_3[0])*(x-z_3[1])/(x*(1-z_3[0])*(1-z_3[1])));
}

double rho3(double x){
  return 1/pi*root3(x)/(1+pow(root3(x),2.0))*(1/(x-z_3[0])+1/(x-z_3[1])-1.0/x);
}

double J3(double z){
  double integ = 0;
  for(int i = 0 ; i< Partnum1;i++){
    double x = x3_1min + (double)dx3_1*i;
    double temp = (double)dx3_1*Rho3_1[i]*z/(z-x);
    integ += (temp==temp) ? temp : 0;
  }
  for(int i = 0 ; i< Partnum2;i++){
    double x = x3_2min + (double)dx3_2*i;
    double temp = (double)dx3_2*Rho3_2[i]*z/(z-x);
    integ += (temp==temp) ? temp : 0;
  }
  return 0.5*integ -0.5;
}

double phi3(double z,int j){
  double integ = 0;
  for(int i = 0 ; i< Partnum1;i++){
    double x = x3_1min + (double)dx3_1*i;
    double temp = dx3_1*Rho3_1[i]*(log((z-x)/(1-x)));
    integ += (temp == temp) ? temp:0;
  }
  for(int i = 0 ; i< Partnum2;i++){
    double x = x3_2min + (double)dx3_2*i;
    double temp = dx3_2*Rho3_2[i]*(log((z-x)/(1-x)));
    integ += (temp == temp) ? temp:0;
  }
  return 0.5*integ-J3_dat[j]*log(z)-0.5*log(z);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////

double root4(double x){
  return sqrt(-(x-z_4[0])*(x-z_4[1])/(x*(1-z_4[0])*(1-z_4[1])));
}

double rho4(double x){
  return 1/pi*root4(x)/(1+pow(root4(x),2.0))*(1/(x-z_4[0])+1/(x-z_4[1])-1.0/x);
}

double J4(double z){
  double integ = 0;
  for(int i = 0 ; i< Partnum1;i++){
    double x = x4_1min + (double)dx4_1*i;
    double temp = (double)dx4_1*Rho4_1[i]*z/(z-x);
    integ += (temp==temp) ? temp : 0;
  }
  for(int i = 0 ; i< Partnum2;i++){
    double x = x4_2min + (double)dx4_2*i;
    double temp = (double)dx4_2*Rho4_2[i]*z/(z-x);
    integ += (temp==temp) ? temp : 0;
  }
  return 0.5*integ -0.5;
}

double phi4(double z,int j){
  double integ = 0;
  for(int i = 0 ; i< Partnum1;i++){
    double x = x4_1min + (double)dx4_1*i;
    double temp = dx4_1*Rho4_1[i]*(log((z-x)/(1-x)));
    integ += (temp == temp) ? temp:0;
  }
  for(int i = 0 ; i< Partnum2;i++){
    double x = x4_2min + (double)dx4_2*i;
    double temp = dx4_2*Rho4_2[i]*(log((z-x)/(1-x)));
    integ += (temp == temp) ? temp:0;
  }
  return 0.5*integ-J4_dat[j]*log(z)-0.5*log(z);
}



int main(){
  ifstream fin("C:/Users/hyoshida/Desktop/floquetic/zero_"+date+ver+".dat");
  fin >> z_disc[0] >> z_disc[1] >> z_disc[2] >> z_disc[3] >>z_disc[4] >>z_disc[5] >>z_disc[6] >>z_disc[7] >> nume[0] >> nume[1] >> nume[2] >> nume[3] >>nume[4] >> deno[0] >> deno[1] >> deno[2] >> deno[3] >>deno[4] >> z_1[0] >> z_1[1] >> z_2[0] >> z_2[1] >> z_3[0] >> z_3[1] >> z_4[0] >> z_4[1];

  xU_1min = z_disc[0];
  xU_1max = z_disc[1];
  dxU_1 = (double)(xU_1max-xU_1min)/partnum1;
  xU_2min = z_disc[2];
  xU_2max = z_disc[3];
  dxU_2 = (double)(xU_2max-xU_2min)/partnum2;
  xU_3min = z_disc[4];
  xU_3max = z_disc[5];
  dxU_3 = (double)(xU_3max-xU_3min)/partnum3;
  xU_4min = z_disc[6];
  xU_4max = z_disc[7];
  dxU_4 = (double)(xU_4max-xU_4min)/partnum4;

  // x1_1min = z_1[0]-100;
  // x1_1max = z_1[0];
  // dx1_1 = (double)(x1_1max-x1_1min)/Partnum1;
  // x1_2min = z_1[1];
  // x1_2max = 0;
  // dx1_2 = (double)(x1_2max-x1_2min)/Partnum2;
  //
  // x2_1min = z_2[0]-100;
  // x2_1max = z_2[0];
  // dx2_1 = (double)(x2_1max-x2_1min)/Partnum1;
  // x2_2min = z_2[1];
  // x2_2max = 0;
  // dx2_2 = (double)(x2_2max-x2_2min)/Partnum2;
  //
  // x3_1min = z_3[0]-100;
  // x3_1max = z_3[0];
  // dx3_1 = (double)(x3_1max-x3_1min)/Partnum1;
  // x3_2min = z_3[1];
  // x3_2max = 0;
  // dx3_2 = (double)(x3_2max-x3_2min)/Partnum2;
  //
  // x4_1min = z_4[0]-100;
  // x4_1max = z_4[0];
  // dx4_1 = (double)(x4_1max-x4_1min)/Partnum1;
  // x4_2min = z_4[1];
  // x4_2max = 0;
  // dx4_2 = (double)(x4_2max-x4_2min)/Partnum2;

  clock_t start = clock();

  for(int j = 0 ;j<partnum1;j++){
    double xU_1 = xU_1min + (double)dxU_1*j;
    RhoU_1[j] = abs(rho(xU_1));
  }

  for(int j = 0 ;j<partnum2;j++){
    double xU_2 = xU_2min + (double)dxU_2*j;
    RhoU_2[j] = abs(rho(xU_2));
  }
  for(int j = 0 ;j<partnum3;j++){
    double xU_3 = xU_3min + (double)dxU_3*j;
    RhoU_3[j] = abs(rho(xU_3));
  }
  for(int j = 0 ;j<partnum4;j++){
    double xU_4 = xU_4min + (double)dxU_4*j;
    RhoU_4[j] = abs(rho(xU_4));
  }

  // for(int j =0 ; j < Partnum1;j++){
  //   double x1_1 = x1_1min + (double)dx1_1*j;
  //   Rho1_1[j] = abs(rho1(x1_1));
  //   double x2_1 = x2_1min + (double)dx2_1*j;
  //   Rho2_1[j] = abs(rho2(x2_1));
  //   double x3_1 = x3_1min + (double)dx3_1*j;
  //   Rho3_1[j] = abs(rho3(x3_1));
  //   double x4_1 = x4_1min + (double)dx4_1*j;
  //   Rho4_1[j] = abs(rho4(x4_1));
  // }
  // for(int j = 0 ; j < Partnum2; j++){
  //   double x1_2 = x1_2min + (double)dx1_2*j;
  //   Rho1_2[j] = abs(rho1(x1_2));
  //   double x2_2 = x2_2min + (double)dx2_2*j;
  //   Rho2_2[j] = abs(rho2(x2_2));
  //   double x3_2 = x3_2min + (double)dx3_2*j;
  //   Rho3_2[j] = abs(rho3(x3_2));
  //   double x4_2 = x4_2min + (double)dx4_2*j;
  //   Rho4_2[j] = abs(rho4(x4_2));
  // }


  clock_t end = clock();
  cout << "Rho : " << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;

  // long long s1=0;
  // for(int i = 0 ; i < partnum1;i++){
  //   s1 += RhoU_1[i];
  // }
  // cout << dxU_1 << endl;
  // cout << s1*dxU_1 << endl;
  //
  // long long s2=0;
  // for(int i = 0 ; i < partnum2;i++){
  //   s2 += RhoU_2[i];
  // }
  // cout << dxU_2 << endl;
  // cout << dxU_2*s2 << endl;
  //
  // long long s3=0;
  // for(int i = 0 ; i < partnum3;i++){
  //   s3 += RhoU_3[i];
  // }
  // cout << dxU_3 << endl;
  // cout << dxU_3*s3 << endl;
  //
  // long long s4=0;
  // for(int i = 0 ; i < partnum4;i++){
  //   s4 += RhoU_4[i];
  // }
  // cout << dxU_4 << endl;
  // cout << dxU_4*s4 << endl;

  string path = "C:/Users/hyoshida/Desktop/floquetic/";
  string ext = ".dat";

  // start = clock();
  // //////////// J-phi(J) plot ////////////////////
  // string filename = path + "phi_"+date+ver + ext;
  // ofstream writing_file;
  // writing_file.open(filename, ios::out);
  //
  // double chi_min = -4;
  // double chi_max = 5;
  // int chi_part = 100;
  //
  // vector<double> Jad_dat(chi_part);
  //
  // for(int j = 0 ; j < chi_part;j++){
  //   double chi = chi_min+(double)(chi_max-chi_min)/chi_part*j;
  //   JU_dat.push_back(J(exp(chi)));
  //   J1_dat.push_back(J1(exp(chi)));
  //   J2_dat.push_back(J2(exp(chi)));
  //   J3_dat.push_back(J3(exp(chi)));
  //   J4_dat.push_back(J4(exp(chi)));
  //   Jad_dat[j] = J1_dat[j]+J2_dat[j]+J3_dat[j]+J4_dat[j];
  // }
  // end = clock();
  // cout << "J : " << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;
  //
  // start = clock();
  // for(int k = 0 ; k < chi_part ; k++){
  //   double chi = chi_min+(double)(chi_max-chi_min)/chi_part*k;
  //   writing_file << exp(chi) << " "<< JU_dat[k] << " " << phi(exp(chi),k) << " " << Jad_dat[k] << " " << phi1(exp(chi),k)+phi2(exp(chi),k)+phi3(exp(chi),k)+phi4(exp(chi),k) << endl;
  // }
  // end = clock();
  // cout << "phi : " << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;

  start = clock();
  ////////////// rho plot //////////////////
  string filename2 = path + "rho_"+date+ver + ext;
  ofstream writing_file2;
  writing_file2.open(filename2, ios::out);

  for(int l = 0 ; l < partnum1;l++){
    double xU_1 = xU_1min + (double)dxU_1*l;
    writing_file2 << xU_1 << " "<< RhoU_1[l] << endl;
  }
  for(int l = 0 ; l < partnum2;l++){
    double xU_2 = xU_2min + (double)dxU_2*l;
    writing_file2 << xU_2 << " "<< RhoU_2[l] << endl;
  }
  for(int l = 0 ; l < partnum3;l++){
    double xU_3 = xU_3min + (double)dxU_3*l;
    writing_file2 << xU_3 << " "<< RhoU_3[l] << endl;
  }
  for(int l = 0 ; l < partnum4;l++){
    double xU_4 = xU_4min + (double)dxU_4*l;
    writing_file2 << xU_4 << " "<< RhoU_4[l] << endl;
  }
  // for(int l = 0 ; l < partnum1;l++){
  //   double xU_1 = xU_1min + (double)dxU_1*l;
  //   writing_file2 << xU_1 << " "<< RhoU_1[l] << " "<< RhoU_d1[l] << " "<< RhoU_g1[l] << endl;
  // }
  // for(int l = 0 ; l < partnum2;l++){
  //   double xU_2 = xU_2min + (double)dxU_2*l;
  //   writing_file2 << xU_2 << " "<< RhoU_2[l] << " "<< RhoU_d2[l] << " "<< RhoU_g2[l] << endl;
  // }
  // for(int l = 0 ; l < partnum3;l++){
  //   double xU_3 = xU_3min + (double)dxU_3*l;
  //   writing_file2 << xU_3 << " "<< RhoU_3[l] << " "<< RhoU_d3[l] << " "<< RhoU_g3[l] << endl;
  // }
  // for(int l = 0 ; l < partnum4;l++){
  //   double xU_4 = xU_4min + (double)dxU_4*l;
  //   writing_file2 << xU_4 << " "<< RhoU_4[l] << " "<< RhoU_d4[l] << " "<< RhoU_g4[l] << endl;
  // }
  end = clock();
  cout << "rho write : " << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;

  Beep(750,200);
}
