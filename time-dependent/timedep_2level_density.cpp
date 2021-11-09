#include <bits/stdc++.h>
#include <fstream>
#include <time.h>
#include <Windows.h>
using namespace std;
#define GNUPLOT_PATH "C:/PROGRA~1/gnuplot/bin/gnuplot.exe"

double pi = 3.141592;

double dt = 1;
int N = 100;

string date = "211109";
string ver = "_1_N";

double r(double x){
  return pi/100*sin(x);
}

// double r_prime(double x){
//   return pi/100*cos(x);
// }

double phi_a(double x){
  return pi/4+pi/10*cos(x);
}

double phi_b(double x){
  return pi/4+pi/10*cos(x);
}

double a_L(double x){
  return 0.5*(1+r(x))*pow(sin(phi_a(x)/2),2.0);
  // return 0.3;
}

double a_R(double x){
  return 0.5*(1+r(x))*pow(cos(phi_a(x)/2),2.0);
  // return 0.3;
}

double b_L(double x){
  return 0.5*(1-r(x))*pow(sin(phi_b(x)/2),2.0);
  // return 0.2;
}

double b_R(double x){
  return 0.5*(1-r(x))*pow(cos(phi_b(x)/2),2.0);
  // return 0.2;
}

double a(double x){
  return a_R(x)+a_L(x);
}

double b(double x){
  return b_R(x)+b_L(x);
}

double z1(double x){
  return (-(pow(b(x)-a(x),2.0)*0.25+b_L(x)*a_L(x)+b_R(x)*a_R(x))+sqrt(pow((pow((b(x)-a(x)),2.0)*0.25+b_L(x)*a_L(x)+b_R(x)*a_R(x)),2.0)-4*a_R(x)*a_L(x)*b_R(x)*b_L(x)))/(2*a_R(x)*b_L(x));
}

double z2(double x){
  return (-(pow(b(x)-a(x),2.0)*0.25+b_L(x)*a_L(x)+b_R(x)*a_R(x))-sqrt(pow((pow((b(x)-a(x)),2.0)*0.25+b_L(x)*a_L(x)+b_R(x)*a_R(x)),2.0)-4*a_R(x)*a_L(x)*b_R(x)*b_L(x)))/(2*a_R(x)*b_L(x));
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

double x1min = -1000;
double x1max = -17;
int partnum1 = 100000;
double dx1 = (double)(x1max-x1min)/partnum1;

double x2min = -3;
double x2max = 0;
int partnum2 = 100000;
double dx2 = (double)(x2max-x2min)/partnum2;

vector<double> Rho1(partnum1),Rho2(partnum2);

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
  return 0.5*(integ -1.0);
}

double phi(double z){
  double integ = 0;
  for(int i = 0 ; i< partnum1;i++){
    double x1 = x1min + (double)dx1*i;
    double temp = dx1*Rho1[i]*(log(z-x1)-log(1-x1));
    integ += (temp == temp) ? temp:0;
  }
  for(int i = 0 ; i< partnum2;i++){
    double x2 = x2min + (double)dx2*i;
    double temp = dx2*Rho2[i]*(log(z-x2)-log(1-x2));
    integ += (temp == temp) ? temp:0;
  }
  return 0.5*(integ-2*J(z)*log(z)-log(z));
}

// double J_0(void){
//   double sum = 0 ;
//   for(int i = 0 ; i < 1000 ; i ++){
//     double theta = 2*pi/1000*i;
//     sum += (double)1/1000*(b_L(theta)*a_R(theta)-a_L(theta)*b_R(theta));
//   }
//   return sum;
// }
//
// double J_1(void){
//   double sum = 0 ;
//   double omega = (double)2*pi/(N*dt);
//   for(int i = 0 ; i < 1000 ; i ++){
//     double theta = 2*pi/1000*i;
//     sum += (double)1/1000*(a_R(theta)+b_R(theta))*r_prime(theta)*omega/2;
//   }
//   return sum;
// }

int main(){
  // clock_t start = clock();

  for(int i = 0 ; i < N ; i++ ){
    double theta = (double)2*pi/N*i;
    Z1[i] = z1(theta);
    Z2[i] = z2(theta);
  }

  double f_min1=10000,f_min2=10000;
  for(int j = 0 ;j<partnum1;j++){
    double x1 = x1min + (double)dx1*j;
    Rho1[j] = abs(rho1(x1));
    if(f_min1 > ave_f(x1)) f_min1 = ave_f(x1);
  }
  for(int j = 0 ;j<partnum2;j++){
    double x2 = x2min + (double)dx2*j;
    Rho2[j] = abs(rho1(x2));
    if(f_min2 > ave_f(x2)) f_min2 = ave_f(x2);
  }



  // cout << J_0() << endl;
  // cout << J_1() << endl;
  // cout << J_0()/J_1() << endl;
  //
  double x1,x2;
  x1 = 2/pi*atan(2*f_min1);
  x2 = 2/pi*atan(2*f_min2);

  // cout << x1 << endl;
  // cout << x2 << endl;

  //////////// J-phi(J) plot ////////////////////
  string path = "C:/Users/hyoshida/Desktop/timedep/";
  string ext = ".dat";
  string filename = path + "phi_"+date+ver + to_string(N) + ext;
  ofstream writing_file;
  writing_file.open(filename, ios::out);

  double chi_min = -4;
  double chi_max = 5;
  int chi_part = 100;
  for(int k = 0 ; k < chi_part ; k++){
    double chi = chi_min+(double)(chi_max-chi_min)/chi_part*k;
    writing_file << exp(chi) << " "<< J(exp(chi)) << " " << phi(exp(chi)) << endl;
  }
  //
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

  //////////// f(z) plot ///////////////
  // string path = "C:/Users/hyoshida/Desktop/timedep/";
  string filename3 = path + "f_" + date+ ver + to_string(N) + ext;
  ofstream writing_file3;
  writing_file3.open(filename3, ios::out);

  writing_file3 << -100 << " " << ave_f(-100) << " " << x1 << " " << x2 << endl;
  for(int l = 1 ; l < 10000;l++){
    double x = -100 + (double)100/10000*l;
    writing_file3 << x << " "<< ave_f(x) << endl;
  }


  // FILE *gp;
  // gp = _popen(GNUPLOT_PATH, "w");
  // fprintf(gp,"set terminal png\n");
  // fprintf(gp,"set output 'C:/Users/hyoshida/Desktop/timedep/time_density_211026_2.png'\n");
  // fprintf(gp,"plot [-0.45:0][-0.3:0.05]'C:/Users/hyoshida/Desktop/timedep/time_simulation_211026_2.dat'\n");
  // fprintf(gp,"replot 'C:/Users/hyoshida/Desktop/timedep/time_density_211026_2.dat' using 2:3\n");
  // // fprintf(gp,"plot 'C:/Users/hyoshida/Desktop/timedep/time_density_211026.dat' using 1:($2>0 ? $2: 1/0) with line\n");
  // // fprintf(gp,"replot 'C:/Users/hyoshida/Desktop/timedep/time_density_211026.dat' using 1:($3>0 ? $3: 1/0) with line\n");
  // pclose(gp);

  // clock_t end = clock();
  // cout << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;
  Beep(750,200);
}
