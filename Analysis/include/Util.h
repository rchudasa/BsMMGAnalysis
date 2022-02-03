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


#endif
