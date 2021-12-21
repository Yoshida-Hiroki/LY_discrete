#include <bits/stdc++.h>
#include <Windows.h>
#include <time.h>
#include <fstream>
using namespace std;

int main(){
  vector<double> trash(100),J_N1_1(100),J_N1_2(100);
  vector<double> phi_N1_1(100),phi_N1_2(100);

  ifstream fin_N1_1("C:/Users/hyoshida/Desktop/floquetic/phi_211221_N1_1.dat");
  for(int i = 0 ; i <100;i++){
    fin_N1_1 >> trash[i] >> J_N1_1[i] >> phi_N1_1[i];
  }

  ifstream fin_N1_2("C:/Users/hyoshida/Desktop/floquetic/phi_211221_N1_2.dat");
  for(int i = 0 ; i <100;i++){
    fin_N1_2 >> trash[i] >> J_N1_2[i] >> phi_N1_2[i];
  }


  ofstream fout("C:/Users/hyoshida/Desktop/floquetic/phiN1.dat");
  for(int i = 0 ; i < 100; i++){
    fout << J_N1_1[i] << " " << phi_N1_1[i] << " " << J_N1_2[i] << " " << phi_N1_2[i] << endl;
  }

}
