#include <bits/stdc++.h>
#include <fstream>
#include <time.h>
#include <Windows.h>
using namespace std;
#define GNUPLOT_PATH "C:/PROGRA~1/gnuplot/bin/gnuplot.exe"

double pi = 3.141592;

double dt = 1;
int N = 2000;

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

// double sup_func(double x){
//   return (double)(a(x)+b(x))*0.5;
// }

// double ave_sup_func(int i){
//   double sum = 0.0;
//   for(int i = 0; i < N ; i++){
//     double theta = 2*pi/(double)N*(double)i;
//     sum += 1/(double)N*sup_func(theta);
//   }
//   return sum;
// }

vector<double> Z1(N),Z2(N);

double ave_LHS(double x){
  double sum = 0.0;
  for(int i = 0; i < N ; i++){
    sum += 1/(double)N*0.5*sqrt(((x-Z1[i])*(x-Z2[i]))/(x*(1-Z1[i])*(1-Z2[i])));
  }
  return sum;
}

double ave_f(double x){
  double sum =0.0;
  for(int i = 0; i < N ; i++){
    sum += 1/(double)N*0.5*sqrt(-((x-Z1[i])*(x-Z2[i]))/(x*(1-Z1[i])*(1-Z2[i])));
  }
  return sum;
}

double ave_g(double x){
  double sum =0.0;
  for(int i = 0; i < N ; i++){
    sum += 1/(double)N*0.5*sqrt(-((x-Z1[i])*(x-Z2[i]))/(x*(1-Z1[i])*(1-Z2[i])))*(1/(x-Z1[i])+1/(x-Z2[i])-1/x);
  }
  return sum;
}

// double R(int i){
//   return ave_sup_func(0)*dt/(1-ave_sup_func(0)*dt);
// }

double rho1(double x){
  return 1/pi/(1+pow((double)ave_f(x)*2,2.0))*ave_g(x)*2;
}

// double rho2(double x){
//   return -1/pi*ave_g(x)*2/(1+pow(ave_f(x)*2,2.0));
// }

double xmin = -25;
double xmax = -0.0000001;
int partnum = 10000;
double dz = (double)(xmax-xmin)/partnum;

vector<double> Rho(partnum);

double J(double z){
  double integ = 0;
  for(int i = 0 ; i< partnum-1;i++){
    double x = xmin + (double)dz*i;
    double temp = dz*Rho[i]*z/(z-x);
    integ += (temp==temp) ? temp : 0;
  }
  return 0.5*(integ -1.0);
}

double phi(double z){
  double integ = 0;

  for(int i = 0 ; i< partnum-1;i++){
    double x = xmin + (double)dz*i;
    double temp = dz*Rho[i]*(log(z-x)-log(1-x));
    integ += (temp == temp) ? temp:0;
  }
  return 0.5*(integ-2*J(z)*log(z)-log(z));
}

int main(){
  clock_t start = clock();
  for(int i = 0 ; i < N ; i++ ){
    double theta = (double)2*pi/N*i;
    Z1[i] = z1(theta);
    Z2[i] = z2(theta);
  }
  for(int j = 0 ;j<partnum;j++){
    double x = xmin + (double)dz*j;
    Rho[j] = abs(rho1(x));
  }


  // string path = "G:/マイドライブ/research/";
  string path = "C:/Users/hyoshida/Desktop/timedep/";
  string type = "time_density_";
  string date = "211026";
  string file = ".dat";
  string filename = path + type + date + file;
  ofstream writing_file;
  writing_file.open(filename, ios::out);

  double chi_min = -4;
  double chi_max = 5;
  int chi_part = 100;
  for(int k = 0 ; k < chi_part ; k++){
    double chi = chi_min+(double)(chi_max-chi_min)/chi_part*k;
    writing_file << exp(chi) << " "<< J(exp(chi)) << " " << phi(exp(chi)) << endl;
  }

  FILE *gp;
  gp = _popen(GNUPLOT_PATH, "w");
  fprintf(gp,"set terminal png\n");
  fprintf(gp,"set output 'C:/Users/hyoshida/Desktop/timedep/time_density_211026.png'\n");
  fprintf(gp,"plot 'C:/Users/hyoshida/Desktop/timedep/time_density_211026.dat' using 2:3\n");
  // fprintf(gp,"plot 'C:/Users/hyoshida/Desktop/timedep/time_density_211026.dat' using 1:($2>0 ? $2: 1/0) with line\n");
  // fprintf(gp,"replot 'C:/Users/hyoshida/Desktop/timedep/time_density_211026.dat' using 1:($3>0 ? $3: 1/0) with line\n");
  pclose(gp);

  clock_t end = clock();
  cout << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;
  Beep(400,1000);
}
