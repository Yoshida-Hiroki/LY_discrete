#include <bits/stdc++.h>
#include <fstream>
#include<iostream>
#include <time.h>
#include <Windows.h>
using namespace std;
#define GNUPLOT_PATH "C:/PROGRA~1/gnuplot/bin/gnuplot.exe"

double pi = 3.141592;

double dt = 1;
int N = 100;

string date = "211121";
string ver = "_1";

//////////////////////////////////////
vector<double> z_disc(8);
vector<double> nume(5);
vector<double> deno(5);

double Trace(double x){
  return deno[0]*pow(x,2.0)+deno[1]*x+deno[2]+deno[3]/x+deno[4]/pow(x,2.0);
}

double R(double x){
  return abs(Trace(1)-2)/Trace(x);
}

double root(double x){
  return sqrt(-(x-z_disc[0])*(x-z_disc[1])*(x-z_disc[2])*(x-z_disc[3])*(x-z_disc[4])*(x-z_disc[5])*(x-z_disc[6])*(x-z_disc[7])/(pow(x,2.0)*(1-z_disc[0])*(1-z_disc[1])*(1-z_disc[2])*(1-z_disc[3])*(1-z_disc[4])*(1-z_disc[5])*(1-z_disc[6])*(1-z_disc[7])));
}

double rho1(double x){
  return 1/pi*(R(x)*root(x))/(1+pow(R(x)*root(x),2.0))*(1/(x-z_disc[0])+1/(x-z_disc[1])+1/(x-z_disc[2])+1/(x-z_disc[3])+1/(x-z_disc[4])+1/(x-z_disc[5])+1/(x-z_disc[6])+1/(x-z_disc[7])-4/x-2*(nume[0]*pow(x,4.0)+nume[1]*pow(x,3.0)+nume[2]*pow(x,2.0)+nume[3]*x+nume[4])/(Trace(x)*pow(x,3.0)));
}


double x1min,x1max,dx1,x2min,x2max,dx2,x3min,x3max,dx3,x4min,x4max,dx4;
long long partnum1 = 100000000,partnum2 = 100000,partnum3 = 10000,partnum4 = 1000;
vector<double> Rho1(partnum1),Rho2(partnum2),Rho3(partnum3),Rho4(partnum4);
vector<double> J_dat;

double J(double z){
  double integ = 0;
  for(int i = 0 ; i< partnum1;i++){
    double x1 = x1min + (double)dx1*i;
    double temp = (double)dx1*Rho1[i]*z/(z-x1);
    integ += (temp==temp) ? temp : 0;
  }
  for(int i = 0 ; i< partnum2;i++){
    double x2 = x2min + (double)dx2*i;
    double temp = (double)dx2*Rho2[i]*z/(z-x2);
    integ += (temp==temp) ? temp : 0;
  }
  for(int i = 0 ; i< partnum3;i++){
    double x3 = x3min + (double)dx3*i;
    double temp = (double)dx3*Rho3[i]*z/(z-x3);
    integ += (temp==temp) ? temp : 0;
  }
  for(int i = 0 ; i< partnum4;i++){
    double x4 = x4min + (double)dx4*i;
    double temp = (double)dx4*Rho4[i]*z/(z-x4);
    integ += (temp==temp) ? temp : 0;
  }
  return 0.5*integ -2.0;
}

double phi(double z,int j){
  double integ = 0;
  for(int i = 0 ; i< partnum1;i++){
    double x1 = x1min + (double)dx1*i;
    double temp = dx1*Rho1[i]*(log((z-x1)/(1-x1)));
    integ += (temp == temp) ? temp:0;
  }
  for(int i = 0 ; i< partnum2;i++){
    double x2 = x2min + (double)dx2*i;
    double temp = dx2*Rho2[i]*(log((z-x2)/(1-x2)));
    integ += (temp == temp) ? temp:0;
  }
  for(int i = 0 ; i< partnum3;i++){
    double x3 = x3min + (double)dx3*i;
    double temp = dx3*Rho3[i]*(log((z-x3)/(1-x3)));
    integ += (temp == temp) ? temp:0;
  }
  for(int i = 0 ; i< partnum4;i++){
    double x4 = x4min + (double)dx4*i;
    double temp = dx4*Rho4[i]*(log((z-x4)/(1-x4)));
    integ += (temp == temp) ? temp:0;
  }
  return 0.5*integ-J_dat[j]*log(z)-2.0*log(z);
}

int main(){
  ifstream fin("C:/Users/hyoshida/Desktop/floquetic/zero_"+date+ver+".dat");
  fin >> z_disc[0] >> z_disc[1] >> z_disc[2] >> z_disc[3] >>z_disc[4] >>z_disc[5] >>z_disc[6] >>z_disc[7] >> nume[0] >> nume[1] >> nume[2] >> nume[3] >>nume[4] >> deno[0] >> deno[1] >> deno[2] >> deno[3] >>deno[4];

  double dx = 0.00000001;
  x1min = z_disc[0]+dx;
  x1max = z_disc[1];
  dx1 = (double)(x1max-x1min)/partnum1;

  x2min = z_disc[2]+dx;
  x2max = z_disc[3];
  dx2 = (double)(x2max-x2min)/partnum2;

  x3min = z_disc[4]+dx;
  x3max = z_disc[5];
  dx3 = (double)(x3max-x3min)/partnum3;

  x4min = z_disc[6]+dx;
  x4max = z_disc[7];
  dx4 = (double)(x4max-x4min)/partnum4;

  clock_t start = clock();

  for(int j = 0 ;j<partnum1;j++){
    double x1 = x1min + (double)dx1*j;
    Rho1[j] = abs(rho1(x1));
  }
  // long long s1=0;
  // for(int i = 0 ; i < partnum1;i++){
  //   s1 += Rho1[i]*dx1;
  // }
  // cout << dx1 << endl;
  // cout << s1 << endl;

  for(int j = 0 ;j<partnum2;j++){
    double x2 = x2min + (double)dx2*j;
    Rho2[j] = abs(rho1(x2));
  }
  // long long s2=0;
  // for(int i = 0 ; i < partnum2;i++){
  //   s2 += Rho2[i];
  // }
  // cout << dx2 << endl;
  // cout << s2 << endl;
  // cout << dx2*s2 << endl;

  for(int j = 0 ;j<partnum3;j++){
    double x3 = x3min + (double)dx3*j;
    Rho3[j] = abs(rho1(x3));
  }
  // long long s3=0;
  // for(int i = 0 ; i < partnum3;i++){
  //   s3 += Rho3[i];
  // }
  // cout << dx3 << endl;
  // cout << s3 << endl;
  // cout << dx3*s3 << endl;

  for(int j = 0 ;j<partnum4;j++){
    double x4 = x4min + (double)dx4*j;
    Rho4[j] = abs(rho1(x4));
  }
  // long long s4=0;
  // for(int i = 0 ; i < partnum4;i++){
  //   s4 += Rho4[i];
  // }
  // cout << dx4 << endl;
  // cout << s4 << endl;
  // cout << dx4*s4 << endl;

  clock_t end = clock();
  cout << "Rho : " << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;



  start = clock();
  //////////// J-phi(J) plot ////////////////////
  string path = "C:/Users/hyoshida/Desktop/floquetic/";
  string ext = ".dat";
  string filename = path + "phi_"+date+ver + ext;
  ofstream writing_file;
  writing_file.open(filename, ios::out);

  double chi_min = -4;
  double chi_max = 5;
  int chi_part = 500;
  for(int j = 0 ; j < chi_part;j++){
    double chi = chi_min+(double)(chi_max-chi_min)/chi_part*j;
    J_dat.push_back(J(exp(chi)));
  }
  end = clock();
  cout << "J : "<<(double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;

  start = clock();
  for(int k = 0 ; k < chi_part ; k++){
    double chi = chi_min+(double)(chi_max-chi_min)/chi_part*k;
    writing_file << exp(chi) << " "<< J_dat[k] << " " << phi(exp(chi),k) << endl;
  }
  end = clock();
  cout << "phi : " << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;

  // ////////////// rho plot //////////////////
  // path = "C:/Users/hyoshida/Desktop/floquetic/";
  // ext = ".dat";
  // filename = path + "rho_"+date+ver + ext;
  // ofstream writing_file2;
  // writing_file2.open(filename, ios::out);
  //
  // for(int l = 0 ; l < partnum1;l++){
  //   double x1 = x1min + (double)dx1*l;
  //   writing_file2 << x1 << " "<< Rho1[l] << endl;
  // }
  // for(int l = 0 ; l < partnum2;l++){
  //   double x2 = x2min + (double)dx2*l;
  //   writing_file2 << x2 << " "<< Rho2[l] << endl;
  // }
  // for(int l = 0 ; l < partnum3;l++){
  //   double x3 = x3min + (double)dx3*l;
  //   writing_file2 << x3 << " "<< Rho3[l] << endl;
  // }
  // for(int l = 0 ; l < partnum4;l++){
  //   double x4 = x4min + (double)dx4*l;
  //   writing_file2 << x4 << " "<< Rho4[l] << endl;
  // }

  Beep(750,200);
}
