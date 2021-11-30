
CXX 		=g++
LIBS		=`root-config --libs`
CXXFLAGS 	=-g3 -std=c++11  `root-config --cflags`
INCPATH     =.

TARGET		=main.exe
all:	$(TARGET) 

main.exe : main.cc Util.h BMMGAnalysis.h BMMGAnalysis.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@ $<  $(LIBS)

clean :
	@rm *.exe
