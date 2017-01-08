
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


