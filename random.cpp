#include "random.hpp"
#include <stdlib.h> 

using namespace std;

/***********************************************************************/
/**  ran0 from numerical recipes                                     ***/
/**                                                                  ***/
/***********************************************************************/
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
/***********************************************************************/
float ran0(long *idum)
{
     long k;
     float ans;
          
     *idum ^= MASK;
     k=(*idum)/IQ;
     *idum=IA*(*idum-k*IQ)-IR*k;
     if (*idum <0) *idum += IM;
     ans=AM*(*idum);
     *idum ^= MASK;
     return (ans);
}
/***********************************************************************/
/***********************************************************************/
/**  ran2 from numerical recipes                                     ***/
/**                                                                  ***/
/***********************************************************************/
#define IR 2836	
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
/***********************************************************************/
float ran1(long *idum){
   int j;
   long k;
   static long iy=0;
   static long iv[NTAB];
   float temp;
   
   if (*idum <=0 || !iy){
      if (-(*idum) < 1) *idum=1;
      else *idum = -(*idum);
      for (j=NTAB+7; j>=0; j--){
         k=(*idum)/IQ;
	 *idum=IA*(*idum-k*IQ)-IR*k;
	 if (*idum < 0) *idum += IM;
	 if (j < NTAB) iv[j]=*idum;
      }
      iy=iv[0];
   }
   k=(*idum)/IQ;
   *idum=IA*(*idum-k*IQ)-IR*k;
   if (*idum < 0) *idum +=IM;
   j=iy/NDIV;
   iy=iv[j];
   iv[j]=*idum;
   if ((temp=AM*iy) > RNMX) return RNMX;
   else return temp;
}
/***********************************************************************/
/***********************************************************************/
/**  ran2 from numerical recipes                                     ***/
/**                                                                  ***/
/***********************************************************************/
#define IM1 2147483563
#define IM2 2147483399
#define AM1 (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB2 32
#define NDIV1 (1+IMM1/NTAB2)
#define EPS2 1.2e-7
#define RNMX2 (1.0-EPS2)
/***********************************************************************/
float ran2(long *idum){
   int j;
   long k;
   static long idum2=123456789;   
   static long iy=0;
   static long iv[NTAB];
   float temp;
   
   if (*idum <= 0){
      if (-(*idum) < 1) *idum=1;
      else *idum = -(*idum);
      idum2=(*idum);
      for (j=NTAB2+7; j>=0; j--){
         k=(*idum)/IQ1;
	 *idum=IA1*(*idum-k*IQ1)-k*IR1;
	 if (*idum < 0) *idum += IM1;
	 if (j < NTAB2) iv[j] = *idum;
      }
      iy=iv[0];
   }
   k=(*idum)/IQ1;
   *idum=IA1*(*idum-k*IQ1)-k*IR1;
   if (*idum < 0) *idum += IM1;
   k=idum2/IQ2;
   idum2=IA2*(idum2-k*IQ2)-k*IR2;
   if (idum2 < 0) idum2 += IM2;
   j=iy/NDIV1;
   iy=iv[j]-idum2;
   iv[j] = *idum;
   if (iy < 1) iy += IMM1;
   if ((temp=AM1*iy) > RNMX2) return RNMX2;
   else return temp;
}

/***********************************************************************/
/**  ran3 from numerical recipes                                     ***/
/**                                                                  ***/
/***********************************************************************/
#define MBIG 1000000000
#define MSEED 161803398	
#define MZ 0
#define FAC (1.0/MBIG)
/***********************************************************************/
float ran3(long *idum){
   static int inext, inextp;
   static long ma[56];
   static int iff=0;
   long mj, mk;
   int i, ii, k;
   
   if (*idum < 0 || iff ==0){
      iff=1;
      mj=labs(MSEED-labs(*idum)); //remenant of C
//      mj=abs(MSEED-abs(*idum));
      mj %= MBIG;
      ma[55]=mj;    
      mk=1;
      for (i=1; i<54; i++){
         ii=(21*i) % 55;
	 ma[ii]=mk;
	 mk=mj-mk;
	 if (mk<MZ) mk+= MBIG;
	 mj=ma[ii];
      }
      for (k=1; k<=4; k++)
         for (i=1; i<=55; i++){
	    ma[i] -= ma[1+(i+30) % 55];
	    if (ma[i] < MZ) ma[i] += MBIG;
	 }
       inext = 0;
       inextp = 31;
       *idum = 1;
   }
   if (++inext == 56) inext=1;
   if (++inextp ==56) inextp=1;
   mj=ma[inext]-ma[inextp];
   if (mj < MZ) mj += MBIG;
   ma[inext]=mj;
   return mj*FAC;
}
/***********************************************************************/
/***********************************************************************/
