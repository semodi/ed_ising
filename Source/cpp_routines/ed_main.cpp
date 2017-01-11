
#include<iostream>
#include<fstream>
#include<cmath>
#include<bitset>
#include<armadillo>
#include<string>
#include<sstream>
#include<complex>
#include"ed_shared.h"
#include"ed_hamiltonian.h"
#include"ed_bitwise.cpp"

using namespace std;

/* The Parameters are read from a separate file ed_config.dat and stored in global 
 * variables
 */ 

/*** Physical parameters ***/
double par[5]; //array to store parameters 
int N; //number of sites 

/*** Other parameters ***/
int eigenv_div = pow(2,0); // # of computed eigenvalues as fraction of matrix size 
const char config_path[] = "ed_config.dat";
const char dict_path[] = ".ed_dict_gen.dat";
const char symmetry_path[] = "ed_config.dat";
extern int NSym;//Number of symmetries
complex<double>offset; // offset for eigenvalues
bool PBC = false; //Periodic boundary conditions
int m[10]; //Multiplicity for each symmetry

/*** Hamiltonian ***/
typedef long(*Hamilpointer)(long w,int i,bool & sign);
extern const int hamilton_size;
extern Hamilpointer Ising[]; //Hamiltonian defined in ed_hamiltonian.cpp


#define PI 3.14159265
#define SQRTWO = 1.414213562;



/*** Flags ***/
bool WRITEDENSE_FLAG = false;





/**************Building Hamiltonian Matrix********************/

inline void fill(int a,int b,arma::uword & loc0, arma::uword & loc1,arma::uword & loc2, arma::uword & loc3,
          complex<double> & value1,complex<double> & value2, complex<double> fillwith,bool sign = true){

	//TODO: Adding the second matrix entry to make sure the matrix is hermitian is now redundant
	// and is taken care of in main() routine. Adjust arguments in fill();


   /** Adds entries to the Hamiltonian matrix
   Assigned values:
   a: line
   b: column
   fillwith: entry
   Stored in:
   loc0: line
   loc1: column
   value1: entry
   */
    const double s = (-1+2*((int)sign));
    
    loc0 = a;
    loc1 = b;
    value1 = fillwith*s;
	
	
	


}

int findState(long field,long dict [], int dict_size){
  /** finding the pure state in dictionary by bisection (binary search algorithm)
  field : pure state that we want to find 
  dict : dictionary between representative and pure states 
    dict[representative] = pure
  dict_size: size of the dictionary
  */

  int L = 0;
  int R = dict_size-1;
  int m ; 
  bool successful = false;

  while(!successful && L <= R){
    
    m = (L+R) >> 1; //floor of (L+R)/2;
    if(dict[m] < field){
      L = m+1;
      continue;
    }else if(dict[m] > field){
      R = m-1;
      continue;
    }else if(dict[m] == field){
      successful = true;
      break;
    }
  }
  if(successful)return m;
  else return -1; 
  
}



int findRepresentative(long field, int new_sites[], int periodicity, long dict [], int dict_size,int & l){
  /* Inverse operation to findState: takes pure state and returns representative
  field : pure state
  dict : dictionary between rep. and pure states
    dict[rep] = pure
  dict_size: number of entries in dictionary
  */
  long save_field = field;
  int cnt = 0;
  for(int i = 0; i < periodicity;++i){
    applySymmetry(field,new_sites);
    if(field < save_field){ // representative candidate
      save_field = field;
      cnt = i+1;
    }
  }
  l = (periodicity - cnt)%periodicity; // # of inverse translations needed to obtain representative state
  return save_field; 
}

void Build_H(arma::umat & locations,arma::cx_vec & values, Hamilpointer HPointer [],int hamilton_size,long dict[],int dict_size){

  
  //Load Symmetries
  string names[NSym];
  int old_sites[NSym][N];
  int new_sites[NSym][N];
  int periodicity[NSym];
  
  
  long w = 0; // pure state
  long wrep = 0; // representative
  int l = 0;
  int b = 0;
  bool sign = true;
  complex<double> symmetry_factor(1.0,0);// complex phase for sym. transf.
  complex<double> value;
  int cnt = 0;
  

  getSymmetryNames(names);
  
  for(int i = 0; i< NSym; ++i){ // get Symmetry transformations from config file
    getSymmetry(old_sites[i],new_sites[i],periodicity[i],m[i],"$" + names[i]);
  }
  
  
  //-----------------------ISING---------------------------------
  // TODO: generalize to more than 3 Hamiltonian terms

  double coupl[3];
  coupl[0] = .5 * (double) getParameter("$par1"); //Diagonal * 0.5 (2*0.5 = 1 when we make H hermitian later on)
  coupl[1] = (double) getParameter("$par2");
  coupl[2] = (double) getParameter("$par3");
  
  
  
  for(int a = 0; a < dict_size;++a){ // Loop over Hilbert space (rep. states)

    complex<double> entry[dict_size]; //entry[b] is Matrix entry M(b,a);
    for(int i = 0;i < dict_size; ++i){
      entry[i] = 0.0;
    }
    
    //Diagonal entry: M(a,a)
    for(int i = 0; i < N; ++i){
      entry[a] += -1+2*(double)test(dict[a],i);
    }
    entry[a] *= coupl[0];
    entry[a] += offset*.5;
     
    //Off-diagonal entries
    for(int H = 1; H < hamilton_size; ++H){ // Loop over Hamiltonians
      for(int i = 0; i < N;++i){ // Loop over sites, TODO: Simplify diagonal term 
        symmetry_factor = complex<double>(1.0,0);
        w = (*HPointer[H])(dict[a],i,sign); //Apply Hamiltonian
        if (!test(w,N)){ //Test-bit
          wrep = w;
          for(int s = 0; s < NSym;++s){  //Find representative state and keep track of symmetry factor for matrix elements
            wrep = findRepresentative(wrep,new_sites[s],periodicity[s],dict,dict_size,l);
            symmetry_factor*= polar(1.0,(double)(-m[s]*PI*2*l)/(double)N);
          }
          b = findState(wrep,dict,dict_size);
          if(b == -1) continue; // state not compatible with momentum
          
          for(int s = 0; s < NSym;++s){  // Find right normalization
            symmetry_factor*= sqrt((double)getMultiplicity(dict[a],periodicity[s],new_sites[s])/
              (double)getMultiplicity(dict[b],periodicity[s],new_sites[s]));
          }
          entry[b] += symmetry_factor*coupl[H]*(-1+2*(double)sign);
        }
      }
    }
    for(int b = 0 ; b < dict_size; ++b){ //write all matrix entries M(b,a) for fixed a 
      if(abs(entry[b]) != 0){
          fill(a,b,locations(0,cnt),locations(1,cnt),locations(0,cnt+1),locations(1,cnt+1),
                      values(cnt),values(cnt+1),entry[b],1);
          cnt += 1;
      }
    }
  }
  
}



/**************** DIAGONALIZE ***********************/

int Diagonalize(arma::sp_cx_mat & H,arma::cx_vec & eigenvalues,int dict_size,int & n_compute,int & n_write){
  
  cout << "Diagonalizing Hamiltonian..." << flush;
  n_compute = (int)dict_size/eigenv_div - 1;
  n_write = n_compute; // n_compute is dynamical and can be increased when convergence issues; n_write (saved # of EVs) is fixed
  bool compute = true;
  
  if(n_compute < N){
    if(dict_size > N) n_compute = N;
    else n_compute = dict_size-1;
  }
  int out_write = getParameter("$out");
  n_write = n_compute;
  if(out_write < n_write) n_write = out_write;
  
  
  do{
    compute = false;
    try{    
      cout << "Computing " << n_compute << " Eigenvalues... " <<flush;
      eigenvalues = arma::eigs_gen(H,n_compute,"sm");
    }catch(const runtime_error & ex){
     
      cout << "trying again..."; 
      if(n_compute*2 < dict_size){ // If possible increase # of computed eigenvalues to conquer convergence issues
        n_compute *= 2;
        compute = true;
      }else{ //If not possible try to diagonalize dense matrix
        arma::cx_mat H_big(H);
        try{
          cout << "diagonalizing dense matrix...";
          eigenvalues = arma::eig_gen(H_big);
        }catch(const runtime_error & ex2){
          cerr << "Convergence issues not solvable" << endl;
          return -1;
        }
      }    
    }
  }while(compute);    
  cout << "done" << flush;      
  return n_write;
}

/**************** MAIN ***********************/

int main(){

  ofstream out;
  
 

  // Get Parameters from "config_path"
  N = getParameter("$Length");
  NSym = getParameter("$NSym");
  WRITEDENSE_FLAG = (bool)getParameter("$DenseFlag");
  eigenv_div = pow(2,getParameter("$eigenv"));
  offset = getParameter("$offset");
  PBC = (bool) getParameter("$pbc");
  
  
  
  
  
  //--------------------------Load dictionary--------------------------------------------
  
  cout << "Loading Dictionary..." << flush;        
  ifstream in_dict;
  in_dict.open(dict_path);
  
  //Skim through file to get dict_size
  int cnt_elements = 0;
  string element;
  while(!in_dict.eof()){
    in_dict >> element;
    ++cnt_elements;
  }
  
  in_dict.close();
  if((cnt_elements-1)%(NSym+2) != 0)throw runtime_error("Something's wrong with dictionary file!");
  
  const int dict_size = (cnt_elements-1)/(NSym+2); 
  
  //Load elements
  
  // The elements are saved in this way:
  // dictionary_index |  state in binary |  multiplicity for sym1  | multiplcity for sym2 | ...
  cnt_elements = 0;
  in_dict.open(dict_path);
  long dict[NSym+2][dict_size]; //dict[file: Column][file: Row (dictionary entries)]
  
  int i = 0;
  int j = 0;
  
  //Read
  while(!in_dict.eof()){
      in_dict >> dict[j++][i];
      if(j == NSym+2){
        j = 0;
        ++i;
      }
  }
      
  //-------------- Build Hamiltonian -------------------------------------
  
  int factor = 20;
  arma::umat locations(2,factor*dict_size,arma::fill::zeros); //locations of matrix entries (row,column)
  arma::cx_vec values(factor*dict_size,arma::fill::zeros); //values of matrix entries (complex!)
            
  
  
  cout << "Building Hamiltonian..." << flush;
  Build_H(locations, values, Ising,hamilton_size,dict[1],dict_size); //Build Hamiltonian
  cout << "done" << endl;
  
  int stopat = 0;
  //Obtain number of non-zero entries in 'locations' TODO: There has to be a better way... Works pretty well though
  for (int i = 0; i< factor*dict_size-1; ++i){
      if (locations(0,i) == 0.0 and locations(1,i) == 0.0 and abs(values(i)) == 0
        and locations(0,i+1) == 0.0 and locations(1,i+1) == 0.0 and abs(values(i+1)) == 0){
        stopat = i-1;
        break;
      }
  }
  cout << "Density: " << (double)stopat/(double)(dict_size*dict_size) << " : " << (double)stopat/(double)dict_size << endl;
  
  
  arma::sp_cx_mat H(locations.submat(0,0,1,stopat),values.subvec(0,stopat));
  arma::cx_vec eigenvalues;
  
  //Hermiticize
  
  H = H + H.t();  
  
  //Write Dense Hamiltonian in file
  
  if(WRITEDENSE_FLAG){
    ofstream dense;
    dense.open("./Matrix.dat");
    arma::cx_mat H_big(H);
    dense << H_big;
    dense.close();
    cout << "done" << endl;
  }
  
  //-------------------------Diagonalizing Hamiltonian ----------------------------------
  
  //Try to get eigenvalues with Lanczos (for small matrices, switches automatically to dense matrices)
  int n_compute;
  int n_write;
  Diagonalize(H,eigenvalues,dict_size,n_compute,n_write);
  eigenvalues = arma::sort(eigenvalues,"ascend");
  
  //------------------------Write in File------------------------------------------
  out.open("out.dat",ios::app);
  
  string outputstring = "";
  stringstream ss;
 
  if(getParameter("$ParFlag") == 1){
    ss << getParameter("$Parity");
    outputstring += "," + ss.str(); 
  }
  if(getParameter("$NumConFlag") == 1){
    ss << getParameter("$PartNum");
    outputstring += "," + ss.str();
  }
  
  for(int i = 0; i < n_write; ++i){
    out << real(eigenvalues(i));
    for(int j = 0; j < NSym; ++j){
      out << "," << m[j];
    }
    out << outputstring;
    out << endl;
    if(imag(eigenvalues(i)) > 1e-7)cerr << "WARNING: HAMILTONIAN MIGHT NOT BE HERMITIAN!" << endl;
  }
  
  out.close();
  
  return 0;
}

