#include <bits/stdc++.h>
#include <fstream>
using namespace std;
#define GNUPLOT_PATH "C:/PROGRA~1/gnuplot/bin/gnuplot.exe"


double pi = 3.141592;


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

double sign(double x){
  return (x>0)-(x<0);
}


int main(){
  // independent parameters
  int dt = 1;
  double omega = 2*pi/50;
  double T = 1000;
  // int N = 1;

  // dependent parameters
  int M = pi/(omega*dt);
  double T_0 = 2*M*dt;
  int N = T/T_0;
  // int T = N*dt;

  int iter = 1000;

  vector<int> J_sim(4*M*N+1);
  vector<double> phi_sim(4*M*N+1);

  srand(time(NULL));
  for(int i = 0 ; i < iter;i++){
    int n_before = 0;
    int n_after = 0;
    int j = 0;

    for(int l = 0 ; l < N ; l++){
      for(int k = 0 ; k < 2*M; k++){
        double rnd = (double)rand()/RAND_MAX;
        double theta = 2*pi/(2*M)*k;

        // n_after = (1-n_before)*0.5*(1+sign(0.4*dt-rnd))+n_before*0.5*(1-sign(0.6*dt-rnd));
        // j += -(1-n_before)*0.5*(1+sign(0.2*dt-rnd))+n_before*0.5*(1+sign(0.3*dt-rnd));
        n_after = (1-n_before)*0.5*(1+sign(b(theta)*dt-rnd))+n_before*0.5*(1-sign(a(theta)*dt-rnd));
        j += -(1-n_before)*0.5*(1+sign(b_R(theta)*dt-rnd))+n_before*0.5*(1+sign(a_R(theta)*dt-rnd));

        n_before = n_after;
      }
    }
    J_sim[(int)(j+2*M*N)] += 1;
  }
  for(int i = 0; i < 4*M*N+1;i++){
    phi_sim[i] = log((double)J_sim[i]/(double)iter)/(2*M*N);
  }

  string path = "C:/Users/hyoshida/Desktop/timedep/";
  string type = "time_simulation_";
  string date = "211024";
  string file = ".dat";
  string filename = path + type + date + file;
  ofstream writing_file;
  writing_file.open(filename, ios::out);

  for(int i = 0; i < 4*M*N+1;i++){
    writing_file << (double)(-2*M*N + i)/(2*M*N) << " " << phi_sim[i] << endl;
  }

  FILE *gp;
  gp = _popen(GNUPLOT_PATH, "w");
  fprintf(gp,"set terminal png\n");
  fprintf(gp,"set output 'C:/Users/hyoshida/Desktop/timedep/time_simulation_211024.png'\n");
  fprintf(gp,"plot 'C:/Users/hyoshida/Desktop/timedep/time_simulation_211024.dat'\n");
  pclose(gp);

}
