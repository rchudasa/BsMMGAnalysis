#ifndef _UTIL_H__
#define _UTIL_H__

#include "iostream"

const Double_t TWO_PI(3.141592653589*2) ;
const Double_t PI(3.141592653589)       ;

Double_t getDPHI( Double_t phi1, Double_t phi2) {

  Double_t dphi = phi1 - phi2;
  
  while( dphi > PI)
        dphi-= TWO_PI; 
  while( dphi < -PI)
        dphi += TWO_PI; 
  //std::cout<<"phi1  : "<<phi1<<" , phi2 : "<<phi2<<" dphi : "<<dphi<<"\n";
  if ( TMath::Abs(dphi) > 3.141592653589 ) {
    cout << " commonUtility::getDPHI error!!! dphi is bigger than 3.141592653589 "<< dphi << endl;
  }
  
  return TMath::Abs(dphi);
  //return dphi;
}

Double_t getDETA(Double_t eta1, Double_t eta2)
{

    return TMath::Abs(eta1 - eta2);
}

Double_t getDR( Double_t eta1, Double_t phi1,Double_t eta2 ,Double_t phi2) {

    Double_t de = getDETA(eta1,eta2);
    Double_t dp = getDPHI(phi1,phi2); 
    return sqrt(de*de + dp*dp);

}

void crossProduct(Double_t vect_A[], Double_t vect_B[], Double_t cross_P[])
 
{
    cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
}

Double_t getMag(Double_t x, Double_t y,Double_t z)
{
        return sqrt(x*x+y*y+z*z);
}


Double_t getDCALineAndPoint(
                                Double_t x1[],
                                Double_t x2[],
                                Double_t p[]
                            )
{
    
    Double_t v1[3],v2[3],v3[3];
    v1[0]= x1[0]-p[0]; v1[1] = x1[1]- p[1] ; v1[2] = x1[1] - p[2];
    v2[0]= x2[0]-p[0]; v2[1] = x2[1]- p[1] ; v2[2] = x2[2] - p[2];
    
    crossProduct(v1,v2,v3);

//    std::cout<<"\nV1  : "<<v1[0]<<","<<v1[1]<<" , "<<v1[2]<<"\n";
//    std::cout<<"\nV2  : "<<v2[0]<<","<<v2[1]<<" , "<<v2[2]<<"\n";
//    std::cout<<"\nc p  : "<<v3[0]<<","<<v3[1]<<" , "<<v3[2]<<"\n";
//    std::cout<<" Mag Num : "<<getMag(v3[0],v3[1],v3[2])<<" , Mag Den : "<<getMag(x1[0]-x2[0],x1[1]-x2[1],x1[2]-x2[2])<<"\n";
    return getMag(v3[0],v3[1],v3[2]) / getMag(x1[0]-x2[0],x1[1]-x2[1],x1[2]-x2[2]);

}

void getTStringFromTag(string tag,std::istringstream &strStream, string inputStr  , TString &var)
{   
    std::string field;
    if(inputStr.compare(tag)==0){
           getline(strStream, field);
           var =  field;
           cout<<" SETTING  "<<tag<<" = "<<var<<"\n";
       }
}


void getBoolFromTag(string tag,std::istringstream &strStream, string inputStr  , Bool_t &var)
{   
    std::string field;
    if(inputStr.compare(tag)==0){
           getline(strStream, field);
           var =  std::atoi(field.c_str()) > 0 ? 1 : 0;
           cout<<" SETTING  "<<tag<<" = "<<var<<"\n";
       }
}

void getFloatFromTag(string tag,std::istringstream &strStream, string inputStr  , Float_t &var)
{   
    std::string field;
    if(inputStr.compare(tag)==0){
             getline(strStream, field);
             var=std::atof(field.c_str());
             cout<<" SETTING  "<<tag<<" = "<<var<<"\n";
       }
}

void getIntFromTag(string tag,std::istringstream &strStream, string inputStr  , Int_t &var)
{   
    std::string field;
    if(inputStr.compare(tag)==0){
             getline(strStream, field);
             var=std::atoi(field.c_str());
             cout<<" SETTING  "<<tag<<" = "<<var<<"\n";
       }
}



void getVetorFillledFromConfigFile( fstream &cfgFile , std::vector<string> &vecList, string beginTag,string endTag, bool verbose)
{
	
    bool cfgModeFlag=false;
    cfgModeFlag=false;
    std::istringstream strStream;
    std::string field;
	string line;
    
    // getting flists
    cfgFile.clear();
	cfgFile.seekp(ios::beg);
    cfgModeFlag=false;
	int nItems(0);
    while(std::getline(cfgFile,line))
	{
	   if(line==beginTag) {cfgModeFlag=true;continue;}
	   if(line==endTag) {cfgModeFlag=false;continue;}
	   if(not cfgModeFlag) continue;
	   if(line=="") continue;
	   if(line=="#END") break;
	   vecList.push_back(line);
	   nItems++;
    }

    if(verbose)
    {
       std::cout<<" Added "<<nItems<<" between "<<beginTag<<" and "<< endTag<<"\n";
       for( auto &name : vecList)
       {
           std::cout<<"\t"<<name<<"\n";
       }
    }

}
#endif
