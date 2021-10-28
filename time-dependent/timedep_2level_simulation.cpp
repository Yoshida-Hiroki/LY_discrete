#include <bits/stdc++.h>
#include <Windows.h>
#include <time.h>
#include <fstream>
using namespace std;
#define GNUPLOT_PATH "C:/PROGRA~1/gnuplot/bin/gnuplot.exe"


double pi = 3.141592;


double r(double x){
  return pi/50*sin(x);
}

double phi_a(double x){
  return pi*3/4+pi/10*cos(x);
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

double sign(double x){
  return (x>0)-(x<0);
}



int main(){
  clock_t start = clock();

  int dt = 1;
  int N = 100;
  int M = 100;

  int iter = 100000;

  vector<int> J_sim(4*M*N+1);
  vector<double> phi_sim(4*M*N+1);

  vector<double> B(N),A(N),B_R(N),A_R(N);
  for(int i = 0 ; i < N ; i++){
    double theta = (double)2*pi/N*i;
    A[i] = a(theta);
    B[i] = b(theta);
    A_R[i] = a_R(theta);
    B_R[i] = b_R(theta);
  }

  srand(time(NULL));
  for(int i = 0 ; i < iter;i++){
    int n_before = 0;
    int n_after = 0;
    int j = 0;

    for(int l = 0 ; l < 2*M ; l++){
      for(int k = 0 ; k < N; k++){
        double rnd = (double)rand()/RAND_MAX;

        // n_after = (1-n_before)*0.5*(1+sign(0.4*dt-rnd))+n_before*0.5*(1-sign(0.6*dt-rnd));
        // j += -(1-n_before)*0.5*(1+sign(0.2*dt-rnd))+n_before*0.5*(1+sign(0.3*dt-rnd));
        n_after = (1-n_before)*0.5*(1+sign(B[k]*dt-rnd))+n_before*0.5*(1-sign(A[k]*dt-rnd));
        j += -(1-n_before)*0.5*(1+sign(B_R[k]*dt-rnd))+n_before*0.5*(1+sign(A_R[k]*dt-rnd));

        n_before = n_after;
      }
    }
    J_sim[(int)(j+2*M*N)] += 1;
  }
  for(int i = 0; i < 4*M*N+1;i++){
    phi_sim[i] = log((double)J_sim[i]/(double)iter)/(2*(double)M*(double)N);
  }

  string path = "C:/Users/hyoshida/Desktop/timedep/";
  string file = "sim_211028_2.dat";
  string filename = path + file;
  ofstream writing_file;
  writing_file.open(filename, ios::out);

  for(int i = 0; i < 4*M*N+1;i++){
    writing_file << (double)(-2*M*N + i)/(2*M*N) << " " << phi_sim[i] << endl;
  }

  // FILE *gp;
  // gp = _popen(GNUPLOT_PATH, "w");
  // fprintf(gp,"set terminal png\n");
  // fprintf(gp,"set output 'C:/Users/hyoshida/Desktop/timedep/time_simulation_211026.png'\n");
  // fprintf(gp,"plot [-0.45:0][-0.3:0.05]'C:/Users/hyoshida/Desktop/timedep/time_simulation_211026.dat'\n");
  // // fprintf(gp,"replot 'C:/Users/hyoshida/Desktop/timedep/time_density_211026.dat' using 2:3 with line lc 2\n");
  // pclose(gp);

  clock_t end = clock();
  cout << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;
  Beep(660, 500);

}
