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

string date = "211129";
string ver = "_N2_3";

vector<double> z_disc(4);
vector<double> nume(3);
vector<double> deno(3);

double Trace(double x){
  return deno[0]*x+deno[1]+deno[2]/x;
}

double R(double x){
  return abs(Trace(1)-2)/Trace(x);
}

double root(double x){
  return sqrt(-(x-z_disc[0])*(x-z_disc[1])*(x-z_disc[2])*(x-z_disc[3])/(pow(x,2.0)*(1-z_disc[0])*(1-z_disc[1])*(1-z_disc[2])*(1-z_disc[3])));
}

double rho1(double x){
  return 1/pi*(abs(Trace(1)-2)*root(x))/(pow(Trace(x),2.0)+pow((Trace(1)-2)*root(x),2.0))*(((1/(x-z_disc[0])+1/(x-z_disc[1])+1/(x-z_disc[2])+1/(x-z_disc[3]))-2.0/x)*Trace(x)-2*(nume[0]*pow(x,2.0)+nume[1]*x+nume[2])/(pow(x,2.0)));
}

// double rhod(double x){
//   return 1/pi*(abs(Trace(1)-2)*root(x))/(pow(Trace(x),2.0)+pow((Trace(1)-2)*root(x),2.0))*((1/(x-z_disc[0])+1/(x-z_disc[1])+1/(x-z_disc[2])+1/(x-z_disc[3]))*Trace(x)-2.0/x*Trace(x));
// }
//
// double rhog(double x){
//   return 1/pi*(abs(Trace(1)-2)*root(x))/(pow(Trace(x),2.0)+pow((Trace(1)-2)*root(x),2.0))*(-2*(nume[0]*pow(x,2.0)+nume[1]*x+nume[2])/(pow(x,2.0)));
// }


double x1min,x1max,dx1;
double x2min,x2max,dx2;
int partnum1 = 1000000;
int partnum2 = 1000000;
vector<double> Rho1(partnum1),Rho2(partnum2);
// vector<double> Rhod1(partnum1),Rhod2(partnum2);
// vector<double> Rhog1(partnum1),Rhog2(partnum2);
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
  return 0.5*integ -1.0;
}
// double Jd(double z){
//   double integ = 0;
//   for(int i = 0 ; i< partnum1+1;i++){
//     double x1 = x1min + (double)dx1*i;
//     double temp = (double)dx1*Rhod1[i]*z/(z-x1);
//     integ += (temp==temp) ? temp : 0;
//   }
//   for(int i = 0 ; i< partnum2+1;i++){
//     double x2 = x2min + (double)dx2*i;
//     double temp = (double)dx2*Rhod2[i]*z/(z-x2);
//     integ += (temp==temp) ? temp : 0;
//   }
//   return 0.5*integ;
// }
// double Jg(double z){
//   double integ = 0;
//   for(int i = 0 ; i< partnum1+1;i++){
//     double x1 = x1min + (double)dx1*i;
//     double temp = (double)dx1*Rhog1[i]*z/(z-x1);
//     integ += (temp==temp) ? temp : 0;
//   }
//   for(int i = 0 ; i< partnum2+1;i++){
//     double x2 = x2min + (double)dx2*i;
//     double temp = (double)dx2*Rhog2[i]*z/(z-x2);
//     integ += (temp==temp) ? temp : 0;
//   }
//   return 0.5*integ;
// }

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
  return 0.5*integ-J_dat[j]*log(z)-log(z);
}

int main(){

  ifstream fin("C:/Users/hyoshida/Desktop/floquetic/zero_"+date+ver+".dat");
  fin >> z_disc[0] >> z_disc[1] >> z_disc[2] >> z_disc[3] >> nume[0] >> nume[1] >> nume[2] >> deno[0] >> deno[1] >> deno[2];

  double dx = 0;
  x1min = z_disc[0]+dx;
  x1max = z_disc[1];
  dx1 = (double)(x1max-x1min)/partnum1;

  x2min = z_disc[2]+dx;
  x2max = z_disc[3];
  dx2 = (double)(x2max-x2min)/partnum2;

  clock_t start = clock();
  for(int j = 0 ;j<partnum1+1;j++){
    double x1 = x1min + (double)dx1*j;
    Rho1[j] = abs(rho1(x1));
    // Rhod1[j] = abs(rhod(x1));
    // Rhog1[j] = abs(rhog(x1));
  }
  for(int j = 0 ;j<partnum2+1;j++){
    double x2 = x2min + (double)dx2*j;
    Rho2[j] = abs(rho1(x2));
    // Rhod2[j] = abs(rhod(x2));
    // Rhog2[j] = abs(rhog(x2));
  }
  clock_t end = clock();
  cout << "Rho : " << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;

  // long long s1=0;
  // for(int i = 0 ; i < partnum1;i++){
  //   s1 += Rho1[i];
  // }
  // cout << dx1 << endl;
  // cout << "sum rho1 : " << s1*dx1 << endl;
  //
  // long long s2=0;
  // for(int i = 0 ; i < partnum2;i++){
  //   s2 += Rho2[i];
  // }
  // cout << dx2 << endl;
  // cout << "sum rho2 : " << dx2*s2 << endl;
  cout << endl;

  string path = "C:/Users/hyoshida/Desktop/floquetic/";
  string ext = ".dat";

  start = clock();
  // //////////// J-phi(J) plot ////////////////////
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
  cout << "J : " << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;

  start = clock();
  for(int k = 0 ; k < chi_part ; k++){
    double chi = chi_min+(double)(chi_max-chi_min)/chi_part*k;
    writing_file << exp(chi) << " "<< J_dat[k] << " " << phi(exp(chi),k) << endl;
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
  //   double x1 = x1min + (double)dx1*l;
  //   writing_file2 << x1 << " "<< Rho1[l] << " "<< Rhod1[l] << " "<< Rhog1[l] << endl;
  // }
  // for(int l = 0 ; l < partnum2;l++){
  //   double x2 = x2min + (double)dx2*l;
  //   writing_file2 << x2 << " "<< Rho2[l] << " "<< Rhod2[l] << " "<< Rhog2[l] << endl;
  // }
  //
  // // for(int l = 0 ; l < partnum1;l++){
  // //   double x1 = x1min + (double)dx1*l;
  // //   writing_file2 << x1 << " "<< Rho1[l] << endl;
  // // }
  // // for(int l = 0 ; l < partnum2;l++){
  // //   double x2 = x2min + (double)dx2*l;
  // //   writing_file2 << x2 << " "<< Rho2[l] << endl;
  // // }
  //
  // end = clock();
  // cout << "rho write : " << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;

  cout << "J   : " << J(1) << endl;
  // cout << "J_d : " << Jd(1) << endl;
  // cout << "J_g : " << Jg(1) << endl;
  // cout << "J_d+J_d-1 : " << Jd(1)+Jg(1)-1 << endl;

  Beep(750,200);
}
