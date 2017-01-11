/** ed_shared.cpp
 * 
 * This is a library shared by all exact diagonalization programs.
 * It contains bitwise operations as well as routines to parse the configuration file.
 * Sebastian Dick, Wuerzburg 2016 **/

//TODO: Make parsing safer!

#include<string>
#include<fstream>
#include<bitset>
#include<iostream>
#include"ed_shared.h"
#include"ed_bitwise.cpp"


using namespace std;



// inline void translate(long & field){ //tranlate (move) one bit to the right with PBCs
//   field = field*2;
//   set_to(field,0,test(field,N));
//   unset(field,N);
// }
// 
// inline void mirror(long & field){
//   bool old_field_array[N];
//   bool new_field_array[N];
//   
//   for(int i = 0; i < N; ++i){
//     old_field_array[i] = test(field,i);
//   }
//   int a = (int) floor(N/2.0);
//   if(N%2 == 1){
//     new_field_array[(int)floor(N/2.0)+1] = old_field_array[(int)floor(N/2.0)+1]; 
//   }
//   for(int i = 0; i < a; ++i){
//     new_field_array[i] = old_field_array[N-1-i];
//     new_field_array[N-1-i] = old_field_array[i];
//   }
//   for(int i = 0; i< N; ++i){
//     set_to(field,i,new_field_array[i]);
//   }
//   
// }

/**************** FILE INPUT************************/
  

long stol(string s){ // Turn string into long 
  long l = 0;
  for(unsigned int i = 0;i < s.size();++i){
    l = 10*l+(s[i]-'0');
  }
  return l;
}

    

int getParameter(string parameter_str){
  /**
   * Extract parameters from "config_path"
   * Look for keywords "parameter_str", the according parameter has to be in the line following it 
  **/  
  
  ifstream config;
  config.open(config_path);
  
  int par = 0;
  string line; 
  while(!config.eof()){
    config >> line;
    if(line == parameter_str){
        config >> line;
        par = stol(line);
        break;
    }
  }
  config.close();
  return par;
  
  
}

/*************** SYMMETRIES *************************/


void getSymmetryNames(string names[]){
  /**
   * Get the symmetry names from config_path
  **/  
  
  ifstream config;
  config.open(config_path);
  
  int cnt = 0;
  string line; 
  while(!config.eof()){
    config >> line;
  if(line == "$Symmetries"){
        config >> line;
        while(line != "end" && cnt < 10){
          names[cnt++] = line;
          config >> line;
        }
        break;
    }
  }
  config.close();
  NSym = cnt;
}

void getSymmetry(int old_sites[], int new_sites[],  int & periodicity,int & eigenvalue, string symmetry_string ){
  /**
   * Read symmetry transformation from file.
   * File has to contain transformation properties for all sites + the periodicity (e.g. 2 for mirror sym.) + eigenvalue m
   */
  
  ifstream symmetry;
  symmetry.open(symmetry_path);
  bool successful = false;
  string line; 
  while(!symmetry.eof()){
    symmetry >> line;
    if(line == symmetry_string){
      successful = true;
      for (int i = 0 ; i < N; ++i){
        symmetry >> line;
        old_sites[i] = (int) stol(line);
      }
      for(int i = 0; i< N;++i){
        symmetry >> line;
        new_sites[i] = (int) stol(line);
      }
      symmetry >> line;
      periodicity = (int) stol(line);
      symmetry >> line;
      eigenvalue = (int) stol(line);
      break;
    }
  }
  if(!successful) cerr << "WARNING: SYMMETRY COULD NOT BE FOUND!" << endl;
  symmetry.close();
}

void applySymmetry(long & field,const int new_sites[]){
  /**
   * Apply the symmetry transformation rules obtained in "getSymmetry" routines
   */
  long new_field = 0;
  for(int i = 0; i < N; ++i){
    set_to(new_field,new_sites[i],test(field,i));
  }
  field = new_field;
}

int getMultiplicity(long field, const int periodicity,const int new_sites[]){
  /**
   * Returns mutiplicity (how often is state obtained in one full period of symmetry transformations , -1 if state not representative
   */
  long new_field = field;
  
  for( int i = 0; i < periodicity; ++i){
    applySymmetry(new_field,new_sites);
    if(new_field < field){return -1;} // not representative (representative always is smallest)
    if(new_field == field){return i+1;} // multiplicity
  }
  return periodicity;
}
