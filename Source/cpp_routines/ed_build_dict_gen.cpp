/** ed_build_dict_gen.cpp
 * 
 * This file contains routines that are used to implement general symmetries that are specified in ed_config.dat
 * and create the corresponding (Hilbert-space) libraries.
 * Sebastian Dick, Wuerzburg 2016 **/

#include<iostream>
#include<fstream>
#include<cmath>
#include<bitset>
#include<string>
#include<sstream>



using namespace std;

const char config_path[] = "ed_config.dat";
const char symmetry_path[] = "ed_config.dat";
int N;
int NSym; //Number of applied symmetries

#include "ed_shared.cpp"



/***************** NUMBER CONSERVATION ***********************/

int getPartNum(long w){ 
  /**
   * Get the "particle number" (# of 1 bits) in given configuration w
   */
  
  int cnt = 0;
  for (int i = 0; i < N; ++i){
    if (test(w,i)){
      ++cnt;
    }
  }
  return cnt;
}

void usePartNum(string pathin, string pathout, int partnum){ //partnum: number of particles, -1: odd parity, -2: even parity
  /*
   * Create dictionary using either conservation of fermion parity or particle number (particle = 1-bit)
   */
  
  ifstream dictionary_in;
  dictionary_in.open(pathin.c_str());
  
  ofstream dictionary_out;
  dictionary_out.open(pathout.c_str(),ios::trunc);
  
  string element; 
  long entry [NSym +2]; // read entry
  int cnt_element = 0; // counts the number of entries that have been read in one line -> detects line break;
  int cnt_dict = 0; //counts elements written in new dictionary;
  bool use_parity = false;
  bool parity = 0;
  bool check = false;
  
  if(partnum < 0){
    use_parity = true;
    parity = partnum+2;
  }

   
  while(!dictionary_in.eof()){ //read lines in old dict
    dictionary_in >> element;
    entry[cnt_element++] = stol(element);   
    if(cnt_element == NSym +2){ // end of line
      check = false;
      if(use_parity) check = ((getPartNum(entry[1])%2) == parity); //Use parity conservation
      else check = (getPartNum(entry[1]) == partnum); //Use particle number conservation
      if(check){ //right number of particles/ right parity
        //WRITE
        dictionary_out << cnt_dict;
        for(int i = 1 ; i < NSym+2; ++i){
          dictionary_out << "\t" <<  entry[i];
        }
        dictionary_out << endl;
        
        ++cnt_dict; //raise index because entry was added to dictionary
      }
      cnt_element = 0;
    } 
  } 
  
  dictionary_in.close();
  dictionary_out.close();
}

/***************** ABELIAN SYMMETRIES ************************/


/*** General symmetries saved in "ed_symmetry.dat ***/

int useSymmetry(string pathin, string pathout, const int  old_sites[] ,const int  new_sites[] , const int periodicity, const int eigenvalue, 
                    const int symmetryN){
  /**
   * Builds a dictionary using the general symmetry specified by new_sites[] and read from config_file
   */
  
  ifstream dictionary_in;
  dictionary_in.open(pathin.c_str());
  
  ofstream dictionary_out;
  dictionary_out.open(pathout.c_str(),ios::trunc);
  
  string element; 
  long entry [NSym + 2]; // read entry
  int cnt_element = 0; // counts the number of entries that have been read in one line -> detects line break;
  int cnt_dict = 0; //counts elements;Reply
  
  int multiplicity = 0;
  
  
  while(!dictionary_in.eof()){ //read lines in old dict
    dictionary_in >> element;
    entry[cnt_element++] = stol(element);   
    if(cnt_element == NSym + 2){ // end of line
      multiplicity = getMultiplicity(entry[1],periodicity,new_sites);
//       cout << entry[1] << "\t" <<  multiplicity << endl;
      if(multiplicity != -1){ //check if representative
        if((multiplicity*eigenvalue)%periodicity ==0){ //check if eigenvalue compatible
          //WRITE
          dictionary_out << cnt_dict;
          for(int i = 1; i< NSym + 2; ++i){
            if(i == symmetryN+2){
              dictionary_out << "\t" << multiplicity ; //Write multiplicity in file (saves comp. time later on)
            }else{
              dictionary_out << "\t" << entry[i] ;
            }
          }
          dictionary_out << endl;
          ++cnt_dict; //raise index because entry was added to dictionary
        }
      }
      cnt_element = 0;
    } 
  } 
  
  dictionary_in.close();
  dictionary_out.close();
  return 0;
}


/****************** CREATE DICTIONARY *****************/

void createDict(){
  
  N = getParameter("$Length");
  long size = pow(2,N);
  
   //Use number conservation laws?
  bool PART_NUM_FLAG = (bool) getParameter("$NumConFlag"); // use particle number conservation?
  bool FERM_PAR_FLAG = (bool) getParameter("$ParFlag"); // use fermion parity conservation?
  int parnum = getParameter("$PartNum"); // Sector: Number of particles
  bool parity = (bool)getParameter("$Parity"); // Sector: Parity
  if(FERM_PAR_FLAG) {parnum = ((int)parity) - 2;}
  

  
  if(PART_NUM_FLAG && FERM_PAR_FLAG){ 
    cerr << "Warning: Symmetries incompatible" << endl;
  }
  
  //Get names of abelian symmetries from file
  string sym_names[10];
  getSymmetryNames(sym_names);
  
  
  int old_sites[N];
  int new_sites[N];
  int periodicity;
  int eigenvalue;
 
  

  
  //Build standard dictionary (without any symmetries) 
  
  ofstream dictionary;
  dictionary.open(".ed_dict_bare.dat",ios::trunc);
  
  for(long i = 0; i< size; ++i){
    dictionary << i << "\t" << i << "\t";
    for(int j = 2;j < NSym + 2; ++j){ 
      dictionary << "0" << "\t";
    }
    dictionary << endl;
  }
  dictionary.close();
  string pathin = ".ed_dict_bare.dat";
  string pathout;
  
  // Use number conservation laws
  
  if(PART_NUM_FLAG || FERM_PAR_FLAG){ 
    if(NSym == 0){ 
      pathout = ".ed_dict_gen.dat";
    }else{
      pathout = ".ed_dict_pn.dat";
    }
    usePartNum(pathin,pathout,parnum);
    pathin = string(pathout);
  }
  
  
  //Use abelian symmetries
  
  for(int i = 0; i < NSym; ++i){
    if( i == NSym-1){
      pathout = ".ed_dict_gen.dat";
    }else{
      stringstream ss;
      ss << i;
      pathout = ".ed_dict_"+ ss.str() +".dat";
    }
    
    getSymmetry(old_sites,new_sites,periodicity,eigenvalue,"$" + sym_names[i]);
    useSymmetry(pathin,pathout,old_sites, new_sites, periodicity, eigenvalue,i);
    pathin = string(pathout);
  }
  
} 
  
/****************** MAIN *****************************/
int main(){
  
  createDict();
  
  
  return 0;
}
