CC=g++
CFLAGS= -c -Wall -O2 -larmadillo
LDFLAGS= -larmadillo
SOURCES= ./cpp_routines/ed_main.cpp ./cpp_routines/ed_shared.cpp ./cpp_routines/ed_hamiltonian.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=../ising

main: $(SOURCES) $(EXECUTABLE)
	    
$(EXECUTABLE): $(OBJECTS) 
	    $(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	    $(CC) $(CFLAGS) $< -o $@

build_dict: ./cpp_routines/ed_build_dict_gen.cpp ./cpp_routines/ed_shared.cpp 
	g++ ./cpp_routines/ed_build_dict_gen.cpp ./cpp_routines/ed_shared.cpp -o ../build_dict

ed_update: ./cpp_routines/ed_update_config.cpp ./cpp_routines/ed_shared.cpp 
	g++  ./cpp_routines/ed_update_config.cpp ./cpp_routines/ed_shared.cpp -o ../ed_update
