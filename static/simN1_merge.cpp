#include <bits/stdc++.h>
#include <Windows.h>
#include <time.h>
#include <fstream>
using namespace std;

int main(){
  vector<string> J(4001),sim_N1_1(4001),sim_N1_2(4001);

  ifstream fin_N2_1("C:/Users/hyoshida/Desktop/floquetic/sim_211221_N1_1.dat");
  for(int i = 0 ; i < 4001 ; i ++){
    fin_N2_1 >> J[i] >> sim_N1_1[i];
  }

  ifstream fin_N2_2("C:/Users/hyoshida/Desktop/floquetic/sim_211221_N1_2.dat");
  for(int i = 0 ; i < 4001 ; i ++){
    fin_N2_2 >> J[i] >> sim_N1_2[i];
  }

  ofstream fout("C:/Users/hyoshida/Desktop/floquetic/simN1.dat");
  for(int i = 0 ; i < 4001 ; i ++){
    fout << J[i] << " " << sim_N1_1[i] << " " << sim_N1_2[i] << endl;
  }

}
