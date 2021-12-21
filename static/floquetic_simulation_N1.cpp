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

string date = "211221";
string ver = "_N1_2";

double a1_L(double x){
  // return 0.3;
  return 0.3;
}

double a1_R(double x){
  // return 0.3;
  return 0.4;
}

double b1_L(double x){
  // return 0.2;
  return 0.2;
}

double b1_R(double x){
  // return 0.2;
  return 0.1;
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

  vector<double> B(1),A(1),B_R(1),A_R(1);
  A[0] = a1(0);
  B[0] = b1(0);
  A_R[0] = a1_R(0);
  B_R[0] = b1_R(0);

  srand(time(NULL));
  for(int i = 0 ; i < iter;i++){
    int n_before = 0;
    int n_after = 0;
    int j = 0;

    for(int l = 0 ; l < 2*M*N ; l++){
      double rnd = (double)rand()/RAND_MAX;

      n_after = (1-n_before)*0.5*(1+sign(B[0]*dt-rnd))+n_before*0.5*(1-sign(A[0]*dt-rnd));
      j += -(1-n_before)*0.5*(1+sign(B_R[0]*dt-rnd))+n_before*0.5*(1+sign(A_R[0]*dt-rnd));

      n_before = n_after;
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
