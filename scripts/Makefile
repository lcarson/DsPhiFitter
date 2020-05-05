# 
# Makefile for programs that include ROOT
# 
ROOTCFLAGS   := $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS     := $(shell $(ROOTSYS)/bin/root-config --libs)

CXXFLAGS      = -O -Wall -fPIC -g -W
CXXFLAGS     += -isystem $(shell $(ROOTSYS)/bin/root-config --incdir)
CXXFLAGS     += $(ROOTCFLAGS) 

LIBS	      = $(ROOTLIBS) -lTreePlayer -lTMVA -lRooFitCore -lRooFit -lRooStats -lMinuit -lFoam -lXMLIO

CC            = c++
BINDIR        = ${PWD}/bin/
SRCDIR        = ${PWD}/src/

OBJS	      = $(BINDIR)Main.o $(BINDIR)Base.o $(BINDIR)CommonTools.o $(BINDIR)Parameters.o 
OBJS			 += $(BINDIR)DsPhiFitting.o $(BINDIR)DsPhiModel.o $(BINDIR)DKst0Fitting.o $(BINDIR)DKst0Model.o $(BINDIR)PhiModel.o $(BINDIR)DsModel.o

DATAOBJS	  = $(SRCDIR)RooHILLdini.C $(SRCDIR)RooHORNSdini.C $(SRCDIR)RooLITTLEHORNSdini.C

#-------------------------------------------------------------
all: run
#-------------------------------------------------------------
run : $(OBJS) $(BINDIR)EventDict.o
			@echo "---------- Doing run..."
			${CC} $(LIBS) -o $(BINDIR)run $(OBJS) $(BINDIR)EventDict.o
#-------------------------------------------------------------

$(BINDIR)EventDict.o : ${DATAOBJS}
				@echo "---------- Generating dictionary ..."
				@rm -f ${SRCDIR}EventDict.C ${SRCDIR}EventDict.h
				@rootcint ${SRCDIR}EventDict.C -c ${DATAOBJS} 
				${CC} $(CXXFLAGS) -I ${PWD} -c ${SRCDIR}EventDict.C -o $(BINDIR)EventDict.o
				@mv ${SRCDIR}EventDict_rdict.pcm ${BINDIR}

$(BINDIR)Main.o : 		$(SRCDIR)Main.C
				@echo "---------- Doing Main.C"
				@mkdir -p bin
				${CC} $(CXXFLAGS) -c $(SRCDIR)Main.C -o $(BINDIR)Main.o

$(BINDIR)Base.o : 		$(SRCDIR)Base.C $(SRCDIR)Base.h
				@echo "---------- Doing Base.C"
				${CC} $(CXXFLAGS) -c $(SRCDIR)Base.C -o $(BINDIR)Base.o

$(BINDIR)CommonTools.o : 	$(SRCDIR)CommonTools.C $(SRCDIR)CommonTools.h
				@echo "---------- Doing CommonTools.C"
				${CC} $(CXXFLAGS) -c $(SRCDIR)CommonTools.C -o $(BINDIR)CommonTools.o

$(BINDIR)Parameters.o : 	$(SRCDIR)Parameters.C $(SRCDIR)Parameters.h
				@echo "---------- Doing Parameters.C"
				${CC} $(CXXFLAGS) -c $(SRCDIR)Parameters.C -o $(BINDIR)Parameters.o

$(BINDIR)DsPhiFitting.o : 	$(SRCDIR)DsPhiFitting.C $(SRCDIR)DsPhiFitting.h
				@echo "---------- Doing DsPhiFitting.C"
				${CC} $(CXXFLAGS) -c $(SRCDIR)DsPhiFitting.C -o $(BINDIR)DsPhiFitting.o

$(BINDIR)DsPhiModel.o :		$(SRCDIR)DsPhiModel.C $(SRCDIR)DsPhiModel.h
				@echo "---------- Doing DsPhiModel.C"
				${CC} $(CXXFLAGS) -c $(SRCDIR)DsPhiModel.C -o $(BINDIR)DsPhiModel.o

$(BINDIR)DKst0Fitting.o :		$(SRCDIR)DKst0Fitting.C $(SRCDIR)DKst0Fitting.h
				@echo "---------- Doing DKst0Fitting.C"
				${CC} $(CXXFLAGS) -c $(SRCDIR)DKst0Fitting.C -o $(BINDIR)DKst0Fitting.o

$(BINDIR)DKst0Model.o :		$(SRCDIR)DKst0Model.C $(SRCDIR)DKst0Model.h
				@echo "---------- Doing DKst0Model.C"
				${CC} $(CXXFLAGS) -c $(SRCDIR)DKst0Model.C -o $(BINDIR)DKst0Model.o

$(BINDIR)PhiModel.o :		$(SRCDIR)PhiModel.C $(SRCDIR)PhiModel.h
				@echo "---------- Doing PhiModel.C"
				${CC} $(CXXFLAGS) -c $(SRCDIR)PhiModel.C -o $(BINDIR)PhiModel.o

$(BINDIR)DsModel.o :		$(SRCDIR)DsModel.C $(SRCDIR)DsModel.h
				@echo "---------- Doing DsModel.C"
				${CC} $(CXXFLAGS) -c $(SRCDIR)DsModel.C -o $(BINDIR)DsModel.o

$(BINDIR)RooHILLdini.o :        $(SRCDIR)RooHILLdini.C $(SRCDIR)RooHILLdini.h
				@echo "---------- Doing RooHILLdini.C"
				${CC} $(CXXFLAGS) -c $(SRCDIR)RooHILLdini.C -o $(BINDIR)RooHILLdini.o

$(BINDIR)RooHORNSdini.o :       $(SRCDIR)RooHORNSdini.C $(SRCDIR)RooHORNSdini.h
				@echo "---------- Doing RooHORNSdini.C"
				${CC} $(CXXFLAGS) -c $(SRCDIR)RooHORNSdini.C -o $(BINDIR)RooHORNSdini.o

$(BINDIR)RooLITTLEHORNSdini.o : $(SRCDIR)RooLITTLEHORNSdini.C $(SRCDIR)RooLITTLEHORNSdini.h
				@echo "---------- Doing RooLITTLEHORNSdini.C"
				${CC} $(CXXFLAGS) -c $(SRCDIR)RooLITTLEHORNSdini.C -o $(BINDIR)RooLITTLEHORNSdini.o


#-------------------------------------------------------
clean:
		@rm -f $(BINDIR)test $(BINDIR)run $(SRCDIR)/*~ core* $(BINDIR)/*.o ./*~ ${SRCDIR}EventDict.C ${SRCDIR}EventDict.h $(BINDIR)*.pcm
		@rmdir $(BINDIR)
# DO NOT DELETE
