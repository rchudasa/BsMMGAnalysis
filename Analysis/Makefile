
CXX 		=g++
LIBS		=`root-config --libs` -lTMVA 

CXXFLAGS 	=-g3 -std=c++11  `root-config --cflags` 
INCPATH     =.

TARGET		= analysis
all:	$(TARGET) 

main : main.cc Util.h BMMGAnalysis.h BMMGAnalysis.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

analysisNtupleMaker : analysisNtupleMaker.cc Util.h BMMGAnalysis.h BMMGAnalysis.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

analysis : analysis.cc Util.h BMMGAnalysis.h BMMGAnalysis.cc RunLumiSelector.h
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

runlumi : runLumi.cc RunLumiSelector.h
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

clean :
	@rm *.exe
