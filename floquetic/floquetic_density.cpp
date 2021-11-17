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

string date = "211117";
string ver = "_1";

double a1_L(double x){
  // return 0.5*(1+r(x))*pow(sin(phi_a(x)/2),2.0);
  return 0.3;
}

double a1_R(double x){
  // return 0.5*(1+r(x))*pow(cos(phi_a(x)/2),2.0);
  return 0.3;
}

double b1_L(double x){
  // return 0.5*(1-r(x))*pow(sin(phi_b(x)/2),2.0);
  return 0.2;
}

double b1_R(double x){
  // return 0.5*(1-r(x))*pow(cos(phi_b(x)/2),2.0);
  return 0.2;
}

double a1(double x){
  return a1_R(x)+a1_L(x);
}

double b1(double x){
  return b1_R(x)+b1_L(x);
}


double a2_L(double x){
  // return 0.5*(1+r(x))*pow(sin(phi_a(x)/2),2.0);
  return 0.1;
}

double a2_R(double x){
  // return 0.5*(1+r(x))*pow(cos(phi_a(x)/2),2.0);
  return 0.1;
}

double b2_L(double x){
  // return 0.5*(1-r(x))*pow(sin(phi_b(x)/2),2.0);
  return 0.4;
}

double b2_R(double x){
  // return 0.5*(1-r(x))*pow(cos(phi_b(x)/2),2.0);
  return 0.4;
}

double a2(double x){
  return a2_R(x)+a2_L(x);
}

double b2(double x){
  return b2_R(x)+b2_L(x);
}

//////////////////////////////////////
vector<double> z_disc(4);

double Trace(double x){
  return (1-b1(0)*dt)*(1-b2(0)*dt)+(b1_L(0)+b1_R(0)/x)*(a2_L(0)+a2_R(0)*x)*pow(dt,2.0)+(1-a1(0)*dt)*(1-a2(0)*dt)+(a1_L(0)+a1_R(0)*x)*(b2_L(0)+b2_R(0)/x)*pow(dt,2.0);
}

double R(double x){
  return (a1(0)*(1-b2(0)*dt)+a2(0)*(1-a1(0)*dt)+b1(0)*(1-a2(0)*dt)+b2(0)*(1-b1(0)*dt))*dt/Trace(x);
}

double root(double x){
  return sqrt(-(x-z_disc[0])*(x-z_disc[1])*(x-z_disc[2])*(x-z_disc[3])/(pow(x,2.0)*(1-z_disc[0])*(1-z_disc[1])*(1-z_disc[2])*(1-z_disc[3])));
}

double rho1(double x){
  return 1/pi*(R(x)*root(x))/(1+pow(R(x)*root(x),2.0))*(1/(x-z_disc[0])+1/(x-z_disc[1])+1/(x-z_disc[2])+1/(x-z_disc[3])-2/x-2*(a1_R(0)*b2_L(0)+a2_R(0)*b1_L(0)-(a1_L(0)*b2_R(0)+a2_L(0)*b1_R(0))/pow(x,2.0))*pow(dt,2.0)/Trace(x));
}


double x1min,x1max,dx1;
double x2min,x2max,dx2;
int partnum1 = 1000000;
int partnum2 = 1000000;
vector<double> Rho1(partnum1),Rho2(partnum2);
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
  fin >> z_disc[0] >> z_disc[1] >> z_disc[2] >> z_disc[3];

  x1min = z_disc[0];
  x1max = z_disc[1];
  dx1 = (double)(x1max-x1min)/partnum1;

  x2min = z_disc[2];
  x2max = z_disc[3];
  dx2 = (double)(x2max-x2min)/partnum2;

  for(int j = 0 ;j<partnum1;j++){
    double x1 = x1min + (double)dx1*j;
    Rho1[j] = abs(rho1(x1));
  }
  for(int j = 0 ;j<partnum2;j++){
    double x2 = x2min + (double)dx2*j;
    Rho2[j] = abs(rho1(x2));
  }


  // cout << J(0.1) << endl;
  // for(int i =0;i<1000;i++){
  //   cout << Rho1[i] << endl;
  // }
  clock_t start = clock();
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
  clock_t end = clock();
  cout << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;

  start = clock();
  for(int k = 0 ; k < chi_part ; k++){
    double chi = chi_min+(double)(chi_max-chi_min)/chi_part*k;
    writing_file << exp(chi) << " "<< J_dat[k] << " " << phi(exp(chi),k) << endl;
  }

  // ////////////// rho plot //////////////////
  // string path = "C:/Users/hyoshida/Desktop/timedep/";
  // string file2 = "rho_211026_16.dat";
  // string filename2 = path + file2;
  // ofstream writing_file2;
  // writing_file2.open(filename2, ios::out);
  //
  // for(int l = 0 ; l < partnum1;l++){
  //     double x1 = x1min + (double)dx1*l;
  //     writing_file2 << x1 << " "<< abs(rho1(x1)) << endl;
  //   }
  // for(int l = 0 ; l < partnum2;l++){
  //     double x2 = x2min + (double)dx2*l;
  //     writing_file2 << x2 << " "<< abs(rho1(x2)) << endl;
  //   }

  // //////////// f(z) plot ///////////////
  // // string path = "C:/Users/hyoshida/Desktop/timedep/";
  // string filename3 = path + "f_" + date+ ver + to_string(N) + ext;
  // ofstream writing_file3;
  // writing_file3.open(filename3, ios::out);
  //
  // writing_file3 << -100 << " " << ave_f(-100) << " " << x1 << " " << x2 << endl;
  // for(int l = 1 ; l < 10000;l++){
  //   double x = -100 + (double)100/10000*l;
  //   writing_file3 << x << " "<< ave_f(x) << endl;
  // }


  // FILE *gp;
  // gp = _popen(GNUPLOT_PATH, "w");
  // fprintf(gp,"set terminal png\n");
  // fprintf(gp,"set output 'C:/Users/hyoshida/Desktop/timedep/time_density_211026_2.png'\n");
  // fprintf(gp,"plot [-0.45:0][-0.3:0.05]'C:/Users/hyoshida/Desktop/timedep/time_simulation_211026_2.dat'\n");
  // fprintf(gp,"replot 'C:/Users/hyoshida/Desktop/timedep/time_density_211026_2.dat' using 2:3\n");
  // // fprintf(gp,"plot 'C:/Users/hyoshida/Desktop/timedep/time_density_211026.dat' using 1:($2>0 ? $2: 1/0) with line\n");
  // // fprintf(gp,"replot 'C:/Users/hyoshida/Desktop/timedep/time_density_211026.dat' using 1:($3>0 ? $3: 1/0) with line\n");
  // pclose(gp);

  end = clock();
  cout << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;
  // Beep(750,200);
}
