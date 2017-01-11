#ifndef SHARED_H
#define SHARED_H


extern int N;
extern int NSym;
extern const char symmetry_path[];
extern const char config_path[];

using namespace std;

long stol(string s);
int getParameter(string parameter_str);
void getSymmetryNames(string names[]);
void applySymmetry(long & field,const int new_sites[]);
void getSymmetry(int old_sites[], int new_sites[],  int & periodicity,int & eigenvalue, string symmetry_string );
int getMultiplicity(long field, const int periodicity,const int new_sites[]);

#endif

