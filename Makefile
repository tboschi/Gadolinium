.PHONY: clean

INCDIR =	include
APPDIR =	app
BINDIR =	bin
LIBDIR =	lib

ROOTLIBS = $(shell root-config --glibs)
ROOTCXXF = $(shell root-config --cflags)
LDFLAGS  := -Wl,--no-as-needed $(LDFLAGS) -L$(LIBDIR) $(ROOTLIBS)
CXXFLAGS := $(CXXFLAGS) -O3 -mavx $(ROOTCXXF) -I$(INCDIR)

CPP =	LoadTree	\

HPP =	

TGT :=	$(CPP:%=$(APPDIR)/%)
DEP :=	$(HPP:%=$(INCDIR)/%.o)

all: $(TGT)
	mkdir -p $(BINDIR)
	mv $(TGT) $(BINDIR)

$(TGT): $(DEP)

clean:
	find $(INCDIR) -mindepth 1 -name "*.o" -delete
	find $(INCDIR) -mindepth 1 -name "*~"  -delete
	find $(APPDIR) -mindepth 1 -name "*~"  -delete
	find $(BINDIR) -mindepth 1 -name "*"   -delete
old:
	g++ mass.c -lgsl -lgomp -lgslcblas
