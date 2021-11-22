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

int iter = 10000;

string date = "211122";
string ver = "_1";

double r(double x){
  return 0.5;
}

double phi_a(double x){
  return pi*0.6666+(double)pi/3.5*cos(x);
}

double phi_b(double x){
  return pi*0.3333+(double)pi/3.5*cos(x);
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
  clock_t start = clock();
  vector<int> J_sim(4*M*N+1);
  vector<double> phi_sim(4*M*N+1);

  vector<double> B(4),A(4),B_R(4),A_R(4);
  for(int i = 0 ; i < 4 ; i++){
    double theta = pi*0.5*i;
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

    for(int l = 0 ; l < 2*M*N ; l++){
      for(int i = 0 ; i < 4; i ++){
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

  clock_t end = clock();
  cout << (double)(end-start) / CLOCKS_PER_SEC<< "sec." << endl;
  Beep(660, 500);

}
