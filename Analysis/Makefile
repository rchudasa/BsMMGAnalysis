
CXX 		=g++
LIBS		=`root-config --libs` -lTMVA 

CXXFLAGS 	=-g3 -std=c++11  `root-config --cflags` 
INCPATH     =.

TARGET		= analysis2018.exe
all:	$(TARGET) 

main.exe : main.cc Util.h BMMGAnalysis.h BMMGAnalysis.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@ $<  $(LIBS)

analysis2018.exe : main.cc Util.h BMMGAnalysis.h BMMGAnalysis.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@ $<  $(LIBS)

clean :
	@rm *.exe
