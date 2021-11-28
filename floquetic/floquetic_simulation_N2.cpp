#include <bits/stdc++.h>
#include <Windows.h>
#include <time.h>
#include <fstream>
using namespace std;
#define GNUPLOT_PATH "C:/PROGRA~1/gnuplot/bin/gnuplot.exe"


double pi = 3.141592;

int dt = 1;
int N = 1;
int M = 1000;

int iter = 1000000;

string date = "211128";
string ver = "_N2_1";

double r(double x){
  return 0.5+0.4999*sin(x);
}

double phi_a(double x){
  return pi*0.5+pi/10*cos(x);
}

double phi_b(double x){
  return phi_a(x);
}

double a1_L(double x){
  return 0.5*(1+r(x))*pow(sin(phi_a(x)/2),2.0);
  // return 0.4;
}

double a1_R(double x){
  return 0.5*(1+r(x))*pow(cos(phi_a(x)/2),2.0);
  // return 0.1;
}

double b1_L(double x){
  return 0.5*(1-r(x))*pow(sin(phi_b(x)/2),2.0);
  // return 0.2;
}

double b1_R(double x){
  return 0.5*(1-r(x))*pow(cos(phi_b(x)/2),2.0);
  // return 0.3;
}

double a1(double x){
  return a1_R(x)+a1_L(x);
}

double b1(double x){
  return b1_R(x)+b1_L(x);
}


double sign(double x){
  return (x>0)-(x<0);
}

int main(){
  clock_t start = clock();

  vector<int> J_sim(4*M*N+1);
  vector<double> phi_sim(4*M*N+1);

  vector<double> B(2),A(2),B_R(2),A_R(2);
  A[0] = a1(0);
  B[0] = b1(0);
  A_R[0] = a1_R(0);
  B_R[0] = b1_R(0);
  A[1] = a1(pi);
  B[1] = b1(pi);
  A_R[1] = a1_R(pi);
  B_R[1] = b1_R(pi);

  srand(time(NULL));
  for(int i = 0 ; i < iter;i++){
    int n_before = 0;
    int n_after = 0;
    int j = 0;

    for(int l = 0 ; l < 2*M*N ; l++){
      for(int i = 0;i < 2;i++){
        double rnd = (double)rand()/RAND_MAX;

        n_after = (1-n_before)*0.5*(1+sign(B[i]*dt-rnd))+n_before*0.5*(1-sign(A[i]*dt-rnd));
        j += -(1-n_before)*0.5*(1+sign(B_R[i]*dt-rnd))+n_before*0.5*(1+sign(A_R[i]*dt-rnd));

        n_before = n_after;
      }
    }
    J_sim[(int)(j+2*M*N)] += 1;
  }
  for(int i = 0; i < 4*M*N+1;i++){
    phi_sim[i] = log((double)J_sim[i]/(double)iter)/(2*(double)M*(double)N);
  }

  string path = "C:/Users/hyoshida/Desktop/floquetic/";
  // string path = "C:/Users/288B/Desktop/";
  string file = ".dat";
  string filename = path+ "sim_"+date+ver + file;
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
