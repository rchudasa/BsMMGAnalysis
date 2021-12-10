/*
 *      author : aravind t s 
 *      date   : 09 Dec 2021
 *      mail   : aravindsugunan@gmail.com
 */

#include "main.h"
#include "RunLumiSelector.h" 

int  main(int argc,char *argv[])
{

   TString fileName("certification/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.root");
   RunLumiSelector goldenJsonSelctor;
   goldenJsonSelctor.loadRunLumiMask(fileName,"RunLumiMask");
    
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
*/

  run=1; lumi=10;
  std::cout<<run<< " , "<<lumi <<"  :  checkRunLumi  = "<<goldenJsonSelctor.checkRunLumi(run,lumi)<<"\n";

  run=315265; lumi=4;
  std::cout<<run<< " , "<<lumi <<"  :  checkRunLumi  = "<<goldenJsonSelctor.checkRunLumi(run,lumi)<<"\n";
  
  run=315265; lumi=20;
  std::cout<<run<< " , "<<lumi <<"  :  checkRunLumi  = "<<goldenJsonSelctor.checkRunLumi(run,lumi)<<"\n";

  run=315265; lumi=58;
  std::cout<<run<< " , "<<lumi <<"  :  checkRunLumi  = "<<goldenJsonSelctor.checkRunLumi(run,lumi)<<"\n";

  run=315265; lumi=60;
  std::cout<<run<< " , "<<lumi <<"  :  checkRunLumi  = "<<goldenJsonSelctor.checkRunLumi(run,lumi)<<"\n";

}
