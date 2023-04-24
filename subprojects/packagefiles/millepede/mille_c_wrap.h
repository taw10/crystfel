#ifdef __cplusplus
extern "C" {
#else
typedef void *Mille;
#endif


extern Mille *mille_new(const char *outFileName,
                        int asBinary,
                        int writeZero);

extern void mille_add_measurement(Mille *m,
                                  int NLC, const float *derLc,
                                  int NGL, const  float *derGl,
                                  const int *label, float rMeas, float sigma);

extern void mille_add_special(Mille *m,
                              int nSpecial,
                              const float *floatings,
                              const int *integers);

extern void mille_delete_last_record(Mille *m);

extern void mille_write_record(Mille *m);

extern void mille_free(Mille *m);

#ifdef __cplusplus
} /* extern "C" */
#endif
