
/************** BITWISE OPERATIONS ***************/

inline void  set(long & field,long bit){ //Set bit "bit" to 1 
  field |= (1 << bit);
}
inline void unset(long & field,long bit){ // Unset bit "bit" (Set to 0)
  field &= ~(1 << bit);
}
inline void toggle(long & field,long bit){ // Toggle bit between 1 and 0
  field ^= (1 << bit);
}
inline bool test(long & field,long bit){ // Test whether bit is 1 or 0
 return (field & ( 1 << bit));
}

inline void set_to(long & field, long bit,bool to){ // Set bit to value "to"
  if(to){
    set(field,bit);
  }else{
    unset(field,bit);
  }
}
