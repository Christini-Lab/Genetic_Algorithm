/*
  Christini lab genetic algorithm
  Modified from IlliGAL GA, reference below
*/

/*
  Single & Multi-Objective Real-Coded Genetic Algorithms Code
  Author: Kumara Sastry
  Illinois Genetic Algorithms Laboratory (IlliGAL)
  Deparment of General Engineering
  University of Illinois at Urbana-Champaign
  104 S. Mathews Ave, Urbana, IL 61801
*/

#include "../include/random.hpp"

Random::Random() {
  iset = 0;
  seed = -makeSeed();
}

Random::Random(const long seedValue) {
  iset = 0;
  seed = -seedValue;
}

Random::~Random() {
}

/// Create a new seed for the random number generator using the time of the day.
long Random::makeSeed() {
  time_t tp;
  struct tm *temp;
  unsigned int temp2, temp3;
  static unsigned int temp4 = 0xe0e1;

  time(&tp);
  temp = localtime(&tp);

  temp2 = (temp->tm_sec<<26)+(temp->tm_min<<20)+(temp->tm_hour<<15)+
      (temp->tm_mday<<10)+(temp->tm_mon<<6);
  temp3 = (temp->tm_year<<13)+(temp->tm_wday<<10)+(temp->tm_yday<<1)+
      temp->tm_isdst;
  temp2 ^= clock()^temp3;

  temp4 = (temp4*0xeeee)%0xffff;
  temp2 ^= temp4<<16;
  temp4 = (temp4*0xaeee)%0xffff;
  temp2 ^= temp4;
  temp2 &= 0x7fffffff;
  return (long int)temp2;
}

///Generate a random permutation of the contents of the input array.
void Random::shuffleArray(int *theArray, int arrayLength) {
  static int ii, iTmp, iWhich, D1, D2;

  for(ii = 1, D1 = 1, D2 =(arrayLength-ii+D1)/D1; D2 > 0; D2--, ii += D1) {
    iWhich = (int)((arrayLength-ii)*random01() + ii);
    iTmp = *(theArray+iWhich-1);
    *(theArray+iWhich-1) = *(theArray+ii-1);
    *(theArray+ii-1) = iTmp;
  }
}
// Returns a uniform random deviate between 0.0 and 1.0. This routine is adopted
// from Numerical Recipes in C. The function is an implementation of random
// number generator of L'Ecuyer with Bays-Durham shuffle (this function is ran 2
// in Numerical Recipies in C, pp 282.)
// Reference: Press, W., Flannery, B., Teukolsky, S., & Vettering, W. (1989).
// Numerical Recipes in C. Cambridge: Cambridge University Press.
// http://lib-www.lanl.gov/numerical/bookcpdf/c7-1.pdf. (http://www.nr.com).
double Random::random01() {
  long j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  if (seed <= 0) {
    if (-(seed) < 1) seed=1;
    else seed = -(seed);
    idum2=(seed);
    for (j=NTAB+7;j>=0;j--) {
      k=(seed)/IQ1;
      seed=IA1*(seed-k*IQ1)-k*IR1;
      if (seed < 0) seed += IM1;
      if (j < NTAB) iv[j] = seed;
    }
    iy=iv[0];
  }
  k=(seed)/IQ1;
  seed=IA1*(seed-k*IQ1)-k*IR1;
  if (seed < 0) seed += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = seed;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}

/// Flip a biased coin. Return 1 if the uniform random deviate is less than the
// probability else return 0.
bool Random::flip( double probability ) {
  if (probability == 1.0) return 1;
  return random01() <= probability;
}

///  Generate a uniform random deviate within the range (minimum, maximum)
double Random::boundedRandom(double minimum, double maximum) {
  return minimum + (maximum-minimum)*random01();
}

/// Generate a uniform a random integer within the range [minimum, maximum)
int Random::boundedIntegerRandom(int minimum, int maximum) {
  return (int)(minimum + (maximum-minimum)*random01());
}

// Normally distributed random variable with deviation Sigma.
// Reference: Press, W., Flannery, B., Teukolsky, S., & Vettering, W. (1989).
// Numerical Recipes in C. Cambridge: Cambridge University Press.
// http://lib-www.lanl.gov/numerical/bookcpdf/c7-3.pdf. (http://www.nr.com).
double Random::normalRandom (double standardDeviation) {
  double fac, rsq, v1, v2;
  if(iset == 0) {
    do {
      v1 = 2.0*random01()-1.0;
      v2 = 2.0*random01()-1.0;
      rsq = v1*v1 + v2*v2;
    } while(rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);
    gset = v1*fac*standardDeviation;
    iset = 1;
    return v2*fac*standardDeviation;
  }
  else {
    iset = 0;
    return gset;
  }
}
