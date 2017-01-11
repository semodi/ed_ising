#include<iostream>
#include<fstream>
#include<cmath>
#include<bitset>
#include<string>
#include<sstream>
#include"ed_shared.h"
#include"ed_bitwise.cpp"

using namespace std;

const char config_path[] = "ed_config.dat";
const char symmetry_path[] = "ed_config.dat";
int N;
int NSym = 1;




int main(){
  ifstream in;
  in.open(config_path);
  ofstream out;
  out.open(".temp.dat");
  
  N = getParameter("$Length");
  int cnt = 0;
  int which = 0;
  
  string line;
  while(!in.eof()){
    in >> line;
    if(line == "$Symcnt"){
      out << line << endl;
      in >> line;
      cnt = stol(line);
      ++cnt;
      out << cnt << endl;
      
      in >> line;
      out << line << endl;
      for(int i = 0; i< cnt ; ++i){
        in >> line;
        out << line << endl;
      }
      which = stol(line);
    }else if(line == "$Sym1"){
      out << line << endl;
      for(int i= 0 ; i < 2*N; ++i){
        in >> line;
        out << line << endl;
      }
      in >> line;
      out << line << endl;
      in >> line;
      out << which << endl;
    }else{
      out << line << endl;
    }
  }
  
  in.close();
  out.close();
  
  in.open(".temp.dat");
  out.open("ed_config.dat",ios::trunc);
  
  cout << endl << endl << "----------------- q = " << which << "-----------------------" <<endl<< endl;
  while(!in.eof()){
    in >> line;
    out << line << endl;
  }
  in.close();
  out.close();
  
}
