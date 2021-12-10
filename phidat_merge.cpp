#include <bits/stdc++.h>
#include <Windows.h>
#include <time.h>
#include <fstream>
using namespace std;

int main(){
  vector<double> trash(100),J_N2_1(100),J_N2_1d(100),J_N2_2(100),J_N2_2d(100),J_N4_1(100),J_N4_1d(100),J_N4_2(100),J_N4_2d(100);
  vector<double> phi_N2_1(100),phi_N2_1d(100),phi_N2_2(100),phi_N2_2d(100),phi_N4_1(100),phi_N4_1d(100),phi_N4_2(100),phi_N4_2d(100);

  ifstream fin_N2_1("C:/Users/hyoshida/Desktop/floquetic/phi_211210_N2_1.dat");
  for(int i = 0 ; i <100;i++){
    fin_N2_1 >> trash[i] >> J_N2_1[i] >> phi_N2_1[i] >> J_N2_1d[i] >> phi_N2_1d[i];
  }
  ifstream fin_N2_2("C:/Users/hyoshida/Desktop/floquetic/phi_211210_N2_2.dat");
  for(int i = 0 ; i <100;i++){
    fin_N2_2 >> trash[i] >> J_N2_2[i] >> phi_N2_2[i] >> J_N2_2d[i] >> phi_N2_2d[i];
  }
  ifstream fin_N4_1("C:/Users/hyoshida/Desktop/floquetic/phi_211210_N4_1.dat");
  for(int i = 0 ; i <100;i++){
    fin_N4_1 >> trash[i] >> J_N4_1[i] >> phi_N4_1[i] >> J_N4_1d[i] >> phi_N4_1d[i];
  }
  ifstream fin_N4_2("C:/Users/hyoshida/Desktop/floquetic/phi_211210_N4_2.dat");
  for(int i = 0 ; i <100;i++){
    fin_N4_2 >> trash[i] >> J_N4_2[i] >> phi_N4_2[i] >> J_N4_2d[i] >> phi_N4_2d[i];
  }

  ofstream fout("C:/Users/hyoshida/Desktop/floquetic/phi_211210.dat");
  for(int i = 0 ; i < 100; i++){
    fout << J_N2_1[i] << " " << phi_N2_1[i] << " " << J_N2_1d[i] << " " << phi_N2_1d[i] << " " << J_N2_2[i] << " " << phi_N2_2[i] << " " << J_N2_2d[i] << " " << phi_N2_2d[i] << " " << J_N4_1[i] << " " << phi_N4_1[i] << " " << J_N4_1d[i] << " " << phi_N4_1d[i] << " " << J_N4_2[i] << " " << phi_N4_2[i] << " " << J_N4_2d[i] << " " << phi_N4_2d[i] << endl;
  }

}
