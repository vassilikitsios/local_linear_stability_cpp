#---------------------------------------------------
BASE_DIR = ../..
include $(BASE_DIR)/Makefile.in

LIBS = -lblitz -lm -L$(LIB_DIR)/ $(STABILTY_LIBS) $(COMPLEX_EIGENSOLVER_LIBS) $(UTILS_LIBS) $(CLAPACK_LIBS)
INCL = -I$(INCLUDE_DIR)/ $(CLAPACK_INC) 

#---------------------------------------------------
.SUFFIXES: .cpp 

.cpp.o:
	$(CPP) $(CPPFLAGS) $(LIBS) $(INCL) -c $*.cpp

.cpp :
	$(CPP) $(CPPFLAGS) $(LIBS) $(INCL) -o $< $@

#---------------------------------------------------
BIN_NAME = local_linear_stability
SRC = Input_Deck_t.cpp main.cpp
OBJS = $(addsuffix .o, $(basename $(SRC))) 

default: $(BIN_NAME)
	(mv $(BIN_NAME) ./bin)

$(BIN_NAME): $(OBJS)
	$(CPP) $(CPPFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -fv ./bin/$(BIN_NAME) *.o *~ ./bin/*~
#---------------------------------------------------
