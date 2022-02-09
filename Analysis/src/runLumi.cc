/*
 *      author : aravind t s 
 *      date   : 09 Dec 2021
 *      mail   : aravindsugunan@gmail.com
 */

#include "main.h"
#include "RunLumiSelector.h" 
#include "RunLumiLogger.h"

int  main(int argc,char *argv[])
{

   RunLumiLogger runLogger;
    
  Int_t run,lumi;

/*
  "315265": [[4, 58]],
  "315267": [[1, 244]],
  "315270": [[1, 633]],
  "315322": [[23, 118], [122, 1354]],
  "315339": [[37, 654]],
  "315357": [[44, 732], [736, 770], [780, 831]],
  "315361": [[40, 619]],
  "315363": [[1, 35], [37, 47], [49, 67], [69, 80], [82, 90]],
  "315366": [[10, 61], [67, 750]],
  "315420": [[28, 920], [924, 942], [954, 1748]],
1 , 10 truth :  0  :  checkRunLumi  = 0
315265 , 4 truth :  1  :  checkRunLumi  = 1
315265 , 20 truth :  1  :  checkRunLumi  = 1
315265 , 58 truth :  1  :  checkRunLumi  = 1
315265 , 60 truth :  0  :  checkRunLumi  = 0
315363 , 1 truth :  1  :  checkRunLumi  = 0
315363 , 36 truth :  0  :  checkRunLumi  = 1
315363 , 48 truth :  0  :  checkRunLumi  = 1
315363 , 70 truth :  1  :  checkRunLumi  = 1
315363 , 90 truth :  1  :  checkRunLumi  = 1
*/

  run=1; lumi=10;
  std::cout<<run<< " , "<<lumi << "\n";
  runLogger.addRunLumi(run,lumi);

  run=1; lumi=12;
  std::cout<<run<< " , "<<lumi << "\n";
  runLogger.addRunLumi(run,lumi);

  run=1; lumi=10;
  std::cout<<run<< " , "<<lumi << "\n";
  runLogger.addRunLumi(run,lumi);

  std::cout<<run<< " , "<<lumi << "\n";
  runLogger.addRunLumi(run,lumi);
  
  runLogger.clear();
  run=315265; lumi=4;
  std::cout<<run<< " , "<<lumi << "\n";
  runLogger.addRunLumi(run,lumi);
  run=315265; lumi=5;
  runLogger.addRunLumi(run,lumi);
  run=315265; lumi=6;
  runLogger.addRunLumi(run,lumi);
  run=315265; lumi=7;
  runLogger.addRunLumi(run,lumi);
  run=315265; lumi=8;
  runLogger.addRunLumi(run,lumi);
  run=315265; lumi=11;
  runLogger.addRunLumi(run,lumi);
  
  run=315265; lumi=20;
  std::cout<<run<< " , "<<lumi << "\n";
  runLogger.addRunLumi(run,lumi);

  run=315265; lumi=58;
  std::cout<<run<< " , "<<lumi << "\n";
  runLogger.addRunLumi(run,lumi);

  run=315265; lumi=59;
  std::cout<<run<< " , "<<lumi << "\n";
  runLogger.addRunLumi(run,lumi);

  run=315265; lumi=60;
  std::cout<<run<< " , "<<lumi << "\n";
  runLogger.addRunLumi(run,lumi);
  
  run=315363; lumi=1;
  std::cout<<run<< " , "<<lumi << "\n";
  runLogger.addRunLumi(run,lumi);

  run=315363; lumi=36;
  std::cout<<run<< " , "<<lumi << "\n";
  runLogger.addRunLumi(run,lumi);
    
  runLogger.clear();
  runLogger.addRunLumi(323775,52);
  runLogger.addRunLumi(323775,53);
  runLogger.addRunLumi(323775,54);
  runLogger.addRunLumi(323775,80);
  runLogger.addRunLumi(323775,81);
  runLogger.addRunLumi(323775,150);
  runLogger.addRunLumi(323775,153);
  runLogger.addRunLumi(323775,154);
  runLogger.addRunLumi(323775,171);
  // "315259": [[1, 172]],
  runLogger.addRunLumi(315259,1);
  runLogger.addRunLumi(315259,2);
  runLogger.addRunLumi(315259,3);
  runLogger.addRunLumi(315259,4);
  runLogger.addRunLumi(315259,5);
  runLogger.addRunLumi(315259,100);
  runLogger.addRunLumi(315259,170);
  runLogger.addRunLumi(315259,171);
  

  runLogger.printAllRunLumis();
  runLogger.printJson();
  runLogger.saveAsJson("runLumi.json");
}
