# This is a short script that can be used to generate simple symmetry transformations
# The functions return one line with the un-transformed sites and a second line with the symmetry 
# transformed sites.
#
# Sebastian Dick, Wuerzburg 2016

def trans(N):
  standard(N)
  s = ""
  for i in range(0,N):
    s += str((i+1)%N) + "\t"
  print(s)  
  
def mirror(N):
  standard(N)
  s = ""
  for i in range(0,N):
    s += str(N-i-1) + "\t"
  print(s)
                 
def standard(N):
  s = ""
  for i in range(0,N):
    s += str(i) + "\t"
  print(s)
