
CXX 		=g++
LIBS		=`root-config --libs` -lTMVA 

CXXFLAGS 	=-g3 -std=c++11  `root-config --cflags` 
INCPATH     = include

TARGET		= analysis
all:	$(TARGET) 

main : src/main.cc  src/BMMGAnalysis.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

analysisNtupleMaker : src/analysisNtupleMaker.cc src/BMMGAnalysis.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

mcAnalysis : src/mcAnalysis.cc src/BMMGAnalysis.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

analysis : src/analysis.cc src/BMMGAnalysis.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

mvaDataMaker : src/mvaDataMaker.cc
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

mvaDataMakerMC : src/mvaDataMakerMC.cc
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

runlumi : src/runLumi.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

mvaTrainer : src/makeMVAModels.cc
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

clean :
	@rm *.exe
