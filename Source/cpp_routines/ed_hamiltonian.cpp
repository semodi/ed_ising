
using namespace std;

#include<string>
#include"ed_shared.h"
#include"ed_hamiltonian.h"
#include"ed_bitwise.cpp"


typedef long(*Hamilpointer)(long w,int i,bool & sign);
extern const int hamilton_size = 3;
Hamilpointer Ising[hamilton_size] = {H_sigma_z,H_kin_1,H_kin_2}; //diagonal term has to be first function


//-------------- Hamiltonian 1D transverse Ising --------------
// Diagonal term
long H_sigma_z(long w, int i, bool & sign){
  sign = test(w,i);
  return w;
}

// S-S+
long H_kin_1(long w, int i, bool & sign){
  sign = 1;
  long right;
  
  if(PBC)right = (i+1)%N; //periodic boundary conditions
  else right = i+1;

  if(right == N){ //end of chain
    set(w,N); //test-bit
    return w;
  }else if(test(w,i) && !test(w,right)){
    unset(w,i);
    set(w,right);
  }else set(w,N);
  
  return w;
}

//S+S+
long H_kin_2(long w, int i, bool & sign){
  sign = 1;
  long right;
  
  if(PBC)right = (i+1)%N; //periodic boundary conditions
  else right = i+1;

  if(right == N){ //end of chain
    set(w,N); //test-bit
    return w;
  }else if(!test(w,i) && !test(w,right)){
    set(w,i);
    set(w,right);
  }else set(w,N);
  
  return w;
}


