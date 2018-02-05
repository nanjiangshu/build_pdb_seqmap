# Makefile for getdsspfilepath.cpp
CC         = g++
CFLAGS     = -Wall -O3
LIBS       = -lm -lmyfunc
LIB_PATH   = ../lib
INCLUDE    = ../src
SRC        = getdsspfilepath.cpp
OBJ        = $(SRC:.cpp=.o) 
EXE        = getdsspfilepath
RM         = /bin/rm -f
CP         = /bin/cp -f
$(EXE): $(OBJ)
	$(CC) $(CFLAGS)  $(OBJ) -L$(LIB_PATH) -o $(EXE)  $(LIBS)
#compile and assemble C++/C source files into object files
# -c flag tells the compiler to create only OBJ files
$(OBJ): $(SRC)
	$(CC) $(CFLAGS) -c -I$(INCLUDE) $(SRC) 
install:
	$(CP)  $(EXE)  ../bin/
clean:
	$(RM)  $(OBJ) $(EXE)
