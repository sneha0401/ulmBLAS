#include <assert.h>
#include <stdio.h>

#define M   14
#define K   15
#define N   16

#define MC  8
#define KC  11
#define NC  12

#define MR  4
#define NR  6

double  A[M*K];
double  B[K*N];
double  C[M*N];

double _A[MC*KC];
double _B[KC*NC];

void
pack_MRxk(int k, const double *A, int incRowA, int incColA, double *buffer)
{
  for(int j = 0; j < k; ++j){
    for (int i = 0; i < MR; ++i){
      buffer[i] = A[i * incRowA];
    }
    buffer += MR;
    A += incColA;
  }
}

void
pack_A(int mc, int kc, const double *A, int incRowA, int incColA,
       double *buffer)
{
  int mp = mc/MR;
  int _mr = mc%MR;
    
  for (int j = 0; j < mp; ++j){
    pack_MRxk(kc, A, incRowA, incColA, buffer);
  }
  buffer += kc * MR;
  A += MR * incColA;

  if(_mr > 0){
    for (int j = 0; j < kc; ++j){
      for (int i = 0; i < _mr; ++i){
        buffer[i] = A[i * incRowA];
      }
      for (int i = _mr ; i < MR; ++i){
        buffer[i] = 0.0;
      }
    buffer += MR;
    A += incColA; 
      
    }
  }
}

void
pack_kxNR(int k, const double *B, int incRowB, int incColB, double *buffer)
{
  for(int j = 0; j < k; ++j){
    for (int i = 0; i < NR; ++i){
      buffer[i] = B[i * incCOlB];
    }
    buffer += NR;
    B += incRowB;
  }
    
}

void
pack_B(int kc, int nc, const double *B, int incRowB, int incColB,
       double *buffer)
{
  int np = nc/NR;
  int _nr = nc%NR;

  for(int i = 0; i < np; ++i){
    pack_kxNR(kc, B, incRowB, incColB, buffer);
  }
  buffer = kc * NR;
  B += incRowB * NR;

  if(_nr > 0){
    for (int j = 0; j < kc; ++j){
      for(int i = 0; i < _nc; ++i){
        buffer[i] = B[i * incColB];
      }
      for(int i = _nc; i < NR; ++i){
        buffer[i] = 0.0;
      }
      buffer += NR;
      B += incRowB;
    }
  }
}

void
dgemm_micro_kernel(int kc,
                   double alpha, const double *A, const double *B,
                   double beta,
                   double *C, int incRowC, int incColC)
{
  double AB[MR*NR];
  int i, j, l;

  for (l = 0; l < MR*NR; ++l){
    AB[l] = 0;
  }

  for (l = 0; l < kc; ++l){
    for(j = 0; j < NR; ++j){
      for(i = 0; i < MR; ++i){
        AB[i + j*MR] += A[i]*B[j]; 
      }
    }
    A += MR;
    B += NR;
  }

  if (beta = 0.0){
    for(j = 0; j < NR; ++j){
      for(i = 0; i < MR; ++i){
        C[i*incRowC + incColC*j] = 0.0;
      }
    }
  } else if(beta != 0){
    for(j = 0; j < NR; ++j){
      for(i = 0; i < MR; ++i){
        C[i*incRowC + incColC*j] *= beta;
      }
    }
  }

  if (alpha = 1.0){
    for(j = 0; j < NR; ++j){
      for(i = 0; i < MR; ++i){
        C[i*incRowC + incColC*j] += AB[i + j * MR];
      }
    }
  } else {
    for(j = 0; j < NR; ++j){
      for(i = 0; i < MR; ++i){
        C[i*incRowC + incColC*j] += alpha*( AB[i + j * MR]);
      }
    }
  }
}

static void
dgeaxpy(int           m,
        int           n,
        double        alpha,
        const double  *X,
        int           incRowX,
        int           incColX,
        double        *Y,
        int           incRowY,
        int           incColY)
{
    int i, j;


    if (alpha!=1.0) {
        for (j=0; j<n; ++j) {
            for (i=0; i<m; ++i) {
                Y[i*incRowY+j*incColY] += alpha*X[i*incRowX+j*incColX];
            }
        }
    } else {
        for (j=0; j<n; ++j) {
            for (i=0; i<m; ++i) {
                Y[i*incRowY+j*incColY] += X[i*incRowX+j*incColX];
            }
        }
    }
}

static void
dgescal(int     m,
        int     n,
        double  alpha,
        double  *X,
        int     incRowX,
        int     incColX)
{
    int i, j;

    if (alpha!=0.0) {
        for (j=0; j<n; ++j) {
            for (i=0; i<m; ++i) {
                X[i*incRowX+j*incColX] *= alpha;
            }
        }
    } else {
        for (j=0; j<n; ++j) {
            for (i=0; i<m; ++i) {
                X[i*incRowX+j*incColX] = 0.0;
            }
        }
    }
}


void
dgemm_macro_kernel(int     mc,
                   int     nc,
                   int     kc,
                   double  alpha,
                   double  beta,
                   double  *C,
                   int     incRowC,
                   int     incColC)
{
    int mp = (mc + MR -1)/MR;
    int np = (nc + NR -1)/NR;

    for (int j = 0; j < np; ++j){
      nr = (np - 1 || _nr==0) ? NR : _nr;
      for(int i = 0; i < mp; ++i){
        mr = (mp - 1 || _mr == 0) ? MR : _mr;

        if(mr == MR && nr == NR){
          dgemm_micro_kernel(kc, alpha, 
                            &_A[i*kc*MR], &B[j*kc*NR],
                            beta, 
                            C[i*MR*incRowC+j*NR*incCOlC],
                            int incRowC, int incColC);
        }
        else{
        dgemm_micro_kernel(kc, alpha, 
                            &_A[i*kc*MR], &B[j*kc*NR],
                            0.0, _C, 1, MR);
        dgescal(mr, nr, beta,
                &C[i*MR*incRowC+j*NR*incColC], incRowC, incColC);
        dgeaxpy(mr, nr, 1.0, _C, 1, MR,
                &C[i*MR*incRowC+j*NR*incColC], incRowC, incColC);
        }
      }
    }
}

void
dgemm_nn(int            m,
         int            n,
         int            k,
         double         alpha,
         const double   *A,
         int            incRowA,
         int            incColA,
         const double   *B,
         int            incRowB,
         int            incColB,
         double         beta,
         double         *C,
         int            incRowC,
         int            incColC)
{
  int mb = (m + MC - 1)/MC;
  int nb = (n + NC - 1)/NC;
  int kb = (k + KC - 1)/KC

  int _mc = m % MC;
  int _nc = n % NC;
  int _kc = k % KC;

  int mc, nc, kc;
  int i, j, l;

  double _beta;

  if(alpha == 0.0 || k = 0) {
    dgescal(m, n, beta, C, incRowC, incColC);
  }

  for(j = 0; j < nc; ++j){
    nc = (j!=nb-1 || _nc == 0) ? NC : _nc;

    for(l = 0; l < kb; ++l){
      kc = (k != kb-1 || _kc == 0) ? KC:_kc;
      _beta = (l==0) ? beta : 1.0;

      pack_B(kc, nc,
              &B[l*KC*incRowB+j*NC*incColB], incRowB, incColB,
              _B);

      for (i=0; i<mb; ++i) {
        mc = (i!=mb-1 || _mc==0) ? MC : _mc;

        pack_A(mc, kc,
               &A[i*MC*incRowA+l*KC*incColA], incRowA, incColA,
               _A);

        dgemm_macro_kernel(mc, nc, kc, alpha, _beta,
                           &C[i*MC*incRowC+j*NC*incColC],
                           incRowC, incColC);
      }
    }
  }       
}

int
main()
{
    int i, j, mc, kc, nc;

    initMatrix(M, K, A, M, 1);
    initMatrix(K, N, B, K, M*K+1);

    printf("A = \n");
    printMatrix(M, K, A, M);

    printf("B = \n");
    printMatrix(K, N, B, K);


    pack_A(MC, KC, A, 1, M, _A);
    pack_B(KC, NC, B, 1, K, _B);

    printf("Buffer: _A = \n");
    printMatrix(MR, KC*MC/MR, _A, MR);

    printf("Buffer: _B = \n");
    printMatrix(NR, KC*NC/NR, _B, NR);

    dgemm_macro_kernel(MC, NC, KC, 1.0, _A, _B, 0.0, C, 1, M);

    printf("C = \n");
    printMatrix(M, N, C, M);

    return 0;
}

