# Makefile for txt2bin-shapestring-data.cpp
CC         = g++
CFLAGS     =
LIBS       = -lm -lmyfunc
LIB_PATH   = $(CASIODATA3)/usr/lib
INCLUDE    = $(CASIODATA3)/program/MyInclude
SRC        = txt2bin-shapestring-data.cpp
OBJ        = $(SRC:.cpp=.o) 
EXE        = txt2bin-shapestring-data
RM         = /bin/rm -f
CP         = /bin/cp -f
$(EXE): $(OBJ)
	$(CC) $(CFLAGS)  $(OBJ) -L$(LIB_PATH) -o $(EXE)  $(LIBS)
#compile and assemble C++/C source files into object files
# -c flag tells the compiler to create only OBJ files
$(OBJ): $(SRC)
	$(CC) $(CFLAGS) -c -I$(INCLUDE) $(SRC) 
install:
	$(CP)  $(EXE)  $(HOME)/bin/
clean:
	$(RM)  $(OBJ)
