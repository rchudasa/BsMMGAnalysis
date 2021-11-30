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



#endif
