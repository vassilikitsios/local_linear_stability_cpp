# Vassili Kitsios 27/02/2007
# -------------------------
BASE_DIR = ../..
include $(BASE_DIR)/Makefile.in

LIBS = -lblitz -lm -L$(LIB_DIR)/ $(UTILS_LIBS) $(CLAPACK_LIBS)
INCL = -I$(INCLUDE_DIR)/

# --------Suffixes---------
.SUFFIXES: .cpp 

.cpp.o:
	$(CPP) $(CPPFLAGS) $(LIBS) $(INCL) -c $*.cpp

.cpp :
	$(CPP) $(CPPFLAGS) $(LIBS) $(INCL) -o $< $@

# -------Dependencies------
ARCHIVE_NAME = libEigenSolver.a
HEADER_NAME = EigenSolver.h

SRC = lapack_mallocs.cpp Eigen_Solution_t.cpp Eigen_Solver_t.cpp

OBJS = $(addsuffix .o, $(basename $(SRC))) 
HEADER_FILES = $(addsuffix .h, $(basename $(SRC)))

default: $(ARCHIVE_NAME)

$(ARCHIVE_NAME): $(OBJS)
	mkdir -p $(LIB_DIR) $(INCLUDE_DIR)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@
	cat $(HEADER_FILES) > $(INCLUDE_DIR)/$(HEADER_NAME)
	mv -v $(ARCHIVE_NAME) $(LIB_DIR)

# This automatically updates all the dependencies shown below
depend:
	makedepend -Y $(SRC)

clean:
	rm -fv *.o *~ ._* $(LIB_DIR)/$(ARCHIVE_NAME) $(INCLUDE_DIR)/$(HEADER_NAME)
