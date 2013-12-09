//#include <cstdlib>
///#include <cstdio>

#include <math.h>
#include <vector>
#include <errno.h>
#include <iostream> //Trying to avoid confilcts when including Rcpp.h in OSX

//#include <R.h>
//#include <Rdefines.h>
#include "hamming_distance.h"
#include <Rcpp.h>
#define RADIAL_HASH_THRESHOLD 0.90

///////////////////////////////////////////////////////////////////////
// Copied from phash (http://www.phash.org/): src/pHash.cpp to avoid installing all the library and its dependences.
// We would like to thank Evan Klinger & David Starkweather for their library
// (ph_crosscorr was improved a little)

// typedef unsigned long long ulong64; //Moved into separate c file to avoid long-long issues with ISO c++ 1998
// int ph_hamming_distance(const ulong64 hash1,const ulong64 hash2){ //Moved into separate c file to avoid long-long issues with ISO c++ 1998
//   ulong64 x = hash1^hash2;
//   const ulong64 m1  = 0x5555555555555555ULL;
//   const ulong64 m2  = 0x3333333333333333ULL;
//   const ulong64 h01 = 0x0101010101010101ULL;
//   const ulong64 m4  = 0x0f0f0f0f0f0f0f0fULL;
//   x -= (x >> 1) & m1;
//   x = (x & m2) + ((x >> 2) & m2);
//   x = (x + (x >> 4)) & m4;
//   return (x * h01)>>56;
// }

int ph_bitcount8(uint8_t val){
  int num = 0;
  while (val){
    ++num;
    val &= val - 1;
  }
  return num;
}
double ph_hammingdistance2(uint8_t *hashA, int lenA, uint8_t *hashB, int lenB){
  if (lenA != lenB){
    return -1.0;
  }
  if ((hashA == NULL) || (hashB == NULL) || (lenA <= 0)){
    return -1.0;
  }
  double dist = 0;
  uint8_t D = 0;
  for (int i=0;i<lenA;i++){
    D = hashA[i]^hashB[i];
    dist = dist + (double)ph_bitcount8(D);
  }
  double bits = (double)lenA*8;
  return dist/bits;
}
typedef struct ph_digest {
  char *id;                   //hash id
  uint8_t *coeffs;            //the head of the digest integer coefficient array
  int size;                   //the size of the coeff array
} Digest;
int ph_crosscorr(const Digest &x,const Digest &y,double &pcc,double threshold){
  //Changes: Jose Maria Perez Ramos (josem.perez.ramos@gmail.com) (Working for ASASEC)
  // Avoid using vectors, as they are not needed
  // Avoid calculating denx and deny every time, as they don't change
  // More optimizations:
  //  Use FFT: xcorr(X,Y) = ifft(fft(X)*Conj(fft(Y)))
  //    # Warning: FFT method is unnormalized
  //    # For 'cyclic' xcorr, no padding is needed
  //  Use x_sumcuadratic - x_mean * x_sum_subvector = x_sumcuadraticdiff
  //    sum((x-mean(x))^2) = sum(x^2) - mean(x) * sum(x)

  int N = y.size;
  int result = 0;

  uint8_t *x_coeffs = x.coeffs;
  uint8_t *y_coeffs = y.coeffs;

  double sumx = 0.0;
  double sumy = 0.0;
  for (int i=0;i < N;i++){
    sumx += x_coeffs[i];
    sumy += y_coeffs[i];
  }
  double meanx = sumx/N;
  double meany = sumy/N;

  double denx = 0.0;
  double deny = 0.0;
  for (int i=0;i<N;i++){
    double tempx = x_coeffs[i]-meanx;
    double tempy = y_coeffs[i]-meany;
    denx += tempx*tempx;
    deny += tempy*tempy;
  }
  double den = sqrt(denx*deny);

  double max = 0;
  for (int d=0;d<N;d++){
    double num = 0.0;
    for (int i=0;i<N;i++)
      num  += (x_coeffs[i]-meanx)*(y_coeffs[(N+i-d)%N]-meany);

    double val = num/den;
    if (val > max)
      max = val;
  }
  pcc = max;
  if (max > threshold)
    result = 1;

  return result;
}
///////////////////////////////////////////////////////////////////////


/** /brief Distances according to different hashes
 *  /var dct: Discrete cosine transformation
 *  /var mw:  Marr Wavelet
 *  /var rd:  Radial Image Hash
 *  /var str: Hash strada
 **/
struct DISTS {
  int dct;
  double mw;
  double rd;
  int str;
};

/** /brief Hashes
 *  /var dct: Discrete cosine transformation
 *  /var mw:  Marr Wavelet
 *  /var rd:  Radial Image Hash
 *  /var str: Hash strada
 **/
struct HASHSTR {
  //ulong64 dct; //No long long in c++ 1998 :( -> Using separate c (std c99) file
  std::vector<uint8_t> mw;
  std::vector<uint8_t> rd;
  std::vector<uint8_t> str;
};

//////////////////////////////////////////
// ulong64 hatol(char *s) { //No long long in c++ 1998 :( -> Using separate c (std c99) file
//   if (s == 0) return 0;
//   if (*s && *s == '0' && *(s + 1) && (*(s+1) == 'x' || *(s+1) == 'X'))
//     s += 2;

//   ulong64 d;
//   ulong64 l;
//   for(l = 0; *s ; ++s)
//     {
//     if (*s >= '0' && *s <= '9')
//       d = *s - '0';
//     else if (*s >= 'a' && *s <= 'f')
//       d = *s - 'a' + 10;
//     else if (*s >= 'A' && *s <= 'F')
//       d = *s - 'A' + 10;
//     else
//       break;
//     l = (l << 4) + d;
//     }

//   return (l);
// }
//////////////////////////////////////////

extern "C" SEXP VideoDistance(SEXP vm1, SEXP vm2){
  Rcpp::List l1(vm1);
  Rcpp::List l2(vm2);
  std::string cadena, cad;

  Rcpp::IntegerVector mw = l1["mw"];
  Rcpp::IntegerVector hs = l1["hstrada"];
  Rcpp::IntegerVector rd = l1["rd"];
  Rcpp::IntegerVector mws= l2["mw"];
  Rcpp::IntegerVector hss= l2["hstrada"];
  Rcpp::IntegerVector rds= l2["rd"];

  HASHSTR hash1;
  HASHSTR hash2;

  cadena=Rcpp::as<std::string>(l1["dct"]);
  //hash1.dct = std::strtoul(&cadena[0],0,16);

  for (int is=0; is <  mw.size();is++){
    hash1.mw.push_back(mw[is]);
  }
  for (int is=0; is <  rd.size();is++){
    hash1.rd.push_back( rd[is]);
  }
  for (int is=0; is <  hs.size();is++){
    hash1.str.push_back(hs[is]);
  }

  cad = Rcpp::as<std::string>(l2["dct"]);
  //hash2.dct = std::strtoul(&cad[0],0,16);
  for (int is=0; is <  mws.size();is++){
    hash2.mw.push_back( mws[is]);
  }
  for (int is=0; is <  rds.size();is++){
    hash2.rd.push_back( rds[is]);
  }
  for (int is=0; is <  hss.size();is++){
    hash2.str.push_back( hss[is]);
  }

  DISTS distances;  // distances according to several types of hashes

  // dct
  //distances.dct = ph_hamming_distance(hash1.dct, hash2.dct);
  distances.dct =ph_hamming_distance(&cadena[0],&cad[0]);

  // mw
  uint8_t* mw1 = &hash1.mw[0];
  uint8_t* mw2 = &hash2.mw[0];
  distances.mw  = ph_hammingdistance2(mw1,(int)hash1.mw.size(),mw2,(int)hash2.mw.size());

  // rd
  double threshold = RADIAL_HASH_THRESHOLD;
  Digest rd1, rd2;
  rd1.size   = (int)hash1.rd.size();
  rd2.size   = (int)hash2.rd.size();
  rd1.coeffs = &hash1.rd[0];
  rd2.coeffs = &hash2.rd[0];
  ph_crosscorr(rd1,rd2,distances.rd,threshold);

  // str
  distances.str = 0;
  for (unsigned int is=0; is < hash1.str.size();is++){
    distances.str += abs((int)(hash1.str[is] - hash2.str[is]));
  }


/*
  PROTECT(out = NEW_NUMERIC(6));
  NUMERIC_POINTER(out)[0]= distances.dct;
  NUMERIC_POINTER(out)[1]= distances.mw;
  NUMERIC_POINTER(out)[2]= distances.rd;
  NUMERIC_POINTER(out)[3]= distances.str;
  NUMERIC_POINTER(out)[4]= 0;
  UNPROTECT(1);
*/

  //Thanks to Manuel Castejón for providing us with the easier Rcpp interface 
  return (wrap( Rcpp::NumericVector::create(distances.dct, distances.mw, distances.rd, distances.str, 0))); 
}
