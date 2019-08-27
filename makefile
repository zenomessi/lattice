CC=g++
CPPFLAGS = -std=gnu++17 -no-pie
INC = -I. -I/usr/local/include -I./Headers
LDLIBS = -lm -lgsl -lgslcblas
CPLUS_INCLUDE_PATH = ./Headers,/usr/local/include

BASFOLD=./
OBJFOLD = $(BASFOLD)Objets/
EXECFOLD = $(BASFOLD)Exec/
SRCFOLD= $(BASFOLD)Source/
HEADFOLD = $(BASFOLD)Headers/
HEADERS = variable.hpp
ANAFOLD = $(SRCFOLD)analysis/

OBJ10 = $(addprefix $(OBJFOLD),lattice_cell_cpp.o parseur_cpp.o system_general_cpp.o system_myo_cpp.o minimizer_cpp_moving.o cell2.o)

$(OBJFOLD)%.o: $(SRCFOLD)%.cpp
	$(CC) $(INC) $(CPPFLAGS) $^ -c -o $@

lattice: $(OBJ10)
	$(CC) $(INC) -O3 -L/usr/local/lib/ $(CPPFLAGS) $^ -static -o $@ $(LDLIBS)

cell2.o: $(HEADFOLD)cell2.hpp $(SRCFOLD)cell2.cpp
	$(CC) $(INC) -O3 $(SRCFOLD)cell2.cpp -c -o $(OBJFOLD)cell2.o

clean:
	rm $(OBJFOLD)*.o
	rm ./displace
