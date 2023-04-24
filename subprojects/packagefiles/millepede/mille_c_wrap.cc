#include "Mille.h"
#include "mille_c_wrap.h"


extern "C" Mille *mille_new(const char *outFileName,
                            int asBinary,
                            int writeZero)
{
	return new Mille(outFileName, asBinary, writeZero);
}


extern "C" void mille_free(Mille *m)
{
	delete m;
}


extern "C" void mille_add_measurement(Mille *m,
                                      int NLC, const float *derLc,
                                      int NGL, const  float *derGl,
                                      const int *label, float rMeas,
                                      float sigma)
{
	m->mille(NLC, derLc, NGL, derGl, label, rMeas, sigma);
}


extern "C" void mille_add_special(Mille *m,
                                  int nSpecial,
                                  const float *floatings,
                                  const int *integers)
{
	m->special(nSpecial, floatings, integers);
}


extern "C" void mille_delete_last_record(Mille *m)
{
	m->kill();
}


extern "C" void mille_write_record(Mille *m)
{
	m->end();
}
