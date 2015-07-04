#----------------------------------------
BASE_DIR = ../..
include $(BASE_DIR)/Makefile.in

LIBS = -lblitz -lm -L$(LIB_DIR)/ $(COMPLEX_EIGENSOLVER_LIBS) $(UTILS_LIBS) $(CLAPACK_LIBS)
INCL = -I$(INCLUDE_DIR)/ 

#----------------------------------------
.SUFFIXES: .cpp 

.cpp.o:
	$(CPP) $(CPPFLAGS) $(LIBS) $(INCL) -c $*.cpp

.cpp :
	$(CPP) $(CPPFLAGS) $(LIBS) $(INCL) -o $< $@

#----------------------------------------
ARCHIVE_NAME = libLocalStability.a
HEADER_NAME = LocalStability.h
SRC = Chebyshev_t.cpp Base_Flow_t.cpp Linear_Navier_Stokes_t.cpp Dispersion_Relationship_Mapping_t.cpp

OBJS = $(addsuffix .o, $(basename $(SRC))) 
HEADER_FILES = $(addsuffix .h, $(basename $(SRC)))

default: $(ARCHIVE_NAME)

$(ARCHIVE_NAME): $(OBJS)
	mkdir -p $(LIB_DIR) $(INCLUDE_DIR)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@
	cat $(HEADER_FILES) > $(INCLUDE_DIR)/$(HEADER_NAME)
	mv -v $(ARCHIVE_NAME) $(LIB_DIR)

clean:
	rm -fv *.o *~ ._* $(LIB_DIR)/$(ARCHIVE_NAME) $(INCLUDE_DIR)/$(HEADER_NAME) 

#----------------------------------------
