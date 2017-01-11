#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

extern bool PBC;
extern int N;

long H_sigma_z(long w, int i, bool & sign);
long H_kin_1(long w, int i, bool & sign);
long H_kin_2(long w, int i, bool & sign);

#endif
