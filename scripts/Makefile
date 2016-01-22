# 
# Makefile for programs that include ROOT
# 
ROOTCFLAGS   := $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS     := $(shell $(ROOTSYS)/bin/root-config --libs)

CXXFLAGS      = -O -Wall -fPIC -g -W
CXXFLAGS     += $(ROOTCFLAGS) 

LIBS	      = $(ROOTLIBS) -lTreePlayer -lTMVA -lRooFitCore -lRooFit -lMinuit -lFoam -lXMLIO

CC            = c++
BINDIR        = ${PWD}/bin/
SRCDIR        = ${PWD}/src/

OBJS	      = $(BINDIR)Main.o $(BINDIR)Base.o $(BINDIR)CommonTools.o $(BINDIR)Parameters.o 
OBJS			 += $(BINDIR)DsPhiFitting.o $(BINDIR)DsPhiModel.o

DATAOBJS	  = $(SRCDIR)RooHILLdini.C $(SRCDIR)RooHORNSdini.C 

#-------------------------------------------------------------
all: run
#-------------------------------------------------------------
run : $(OBJS) $(BINDIR)EventDict.o
			${CC} $(LIBS) -o $(BINDIR)run $(OBJS) $(BINDIR)EventDict.o
#-------------------------------------------------------------

$(BINDIR)EventDict.o : ${DATAOBJS}
				@echo "Generating dictionary ..."
				@rm -f ${SRCDIR}EventDict.C ${SRCDIR}EventDict.h
				@rootcint ${SRCDIR}EventDict.C -c ${DATAOBJS}
				${CC} $(CXXFLAGS) -c ${SRCDIR}EventDict.C -o $(BINDIR)EventDict.o

$(BINDIR)Main.o : 		$(SRCDIR)Main.C
				@mkdir -p bin
				${CC} $(CXXFLAGS) -c $(SRCDIR)Main.C -o $(BINDIR)Main.o

$(BINDIR)Base.o : 		$(SRCDIR)Base.C $(SRCDIR)Base.h
				${CC} $(CXXFLAGS) -c $(SRCDIR)Base.C -o $(BINDIR)Base.o

$(BINDIR)CommonTools.o : 	$(SRCDIR)CommonTools.C $(SRCDIR)CommonTools.h
				${CC} $(CXXFLAGS) -c $(SRCDIR)CommonTools.C -o $(BINDIR)CommonTools.o

$(BINDIR)Parameters.o : 	$(SRCDIR)Parameters.C $(SRCDIR)Parameters.h
				${CC} $(CXXFLAGS) -c $(SRCDIR)Parameters.C -o $(BINDIR)Parameters.o

$(BINDIR)DsPhiFitting.o : 	$(SRCDIR)DsPhiFitting.C $(SRCDIR)DsPhiFitting.h
				${CC} $(CXXFLAGS) -c $(SRCDIR)DsPhiFitting.C -o $(BINDIR)DsPhiFitting.o

$(BINDIR)DsPhiModel.o :		$(SRCDIR)DsPhiModel.C $(SRCDIR)DsPhiModel.h
				${CC} $(CXXFLAGS) -c $(SRCDIR)DsPhiModel.C -o $(BINDIR)DsPhiModel.o

$(BINDIR)RooHILLdini.o :        $(SRCDIR)RooHILLdini.C $(SRCDIR)RooHILLdini.h
				${CC} $(CXXFLAGS) -c $(SRCDIR)RooHILLdini.C -o $(BINDIR)RooHILLdini.o

$(BINDIR)RooHORNSdini.o :       $(SRCDIR)RooHORNSdini.C $(SRCDIR)RooHORNSdini.h
				${CC} $(CXXFLAGS) -c $(SRCDIR)RooHORNSdini.C -o $(BINDIR)RooHORNSdini.o


#-------------------------------------------------------
clean:
		@rm -f $(BINDIR)test $(BINDIR)run $(SRCDIR)/*~ core* $(BINDIR)/*.o ./*~ ${SRCDIR}EventDict.C ${SRCDIR}EventDict.h
		@rmdir $(BINDIR)
# DO NOT DELETE
