#include <stdio.h>
#include "gmp.h"
#include "point.h"

void double_add(PROJECTIVE_POINT R, PROJECTIVE_POINT P, const mpz_t N)
{
	mpz_t B,C,D,E,F,H,J;

	mpz_inits(B,C,D,E,F,H,J,NULL);

	mpz_add(B,P->X,P->Y); //B = X1+Y1
	mpz_pow_ui(B,B,2); //B = (X1+Y1)^2
	mpz_mod(B,B,N);
	mpz_pow_ui(C,P->X,2); //C = X1^2
	mpz_mod(C,C,N);
	mpz_pow_ui(D,P->Y,2); //D = Y1^2 
	mpz_mod(D,D,N);
	mpz_mul_si(E,C,-1); //E=-C
	mpz_add(F,E,D); //F = E+D
	mpz_pow_ui(H,P->Z,2); //H = Z1^2
	mpz_mod(H,H,N);
	mpz_sub(J,F,H); //J = F-H
	mpz_sub(J,J,H); //J = F-2H
	
	/* calcurate X_3 */
	mpz_sub(R->X,B,C); //X3=B-C
	mpz_sub(R->X,R->X,D); //X3=B-C-D
	mpz_mul_mod(R->X,R->X,J,N); //X3
	
	/* calcurate Y_3 */
	mpz_sub(R->Y,E,D); //Y3 = E-D
	mpz_mul_mod(R->Y,R->Y,F,N); //Y3

	/* calcurate Z_3 */
	mpz_mul_mod(R->Z,F,J,N); //Z3

	mpz_clears(B,C,D,E,F,H,J,NULL);
}


void dedicated_doubling(EXTENDED_POINT R, const EXTENDED_POINT P, const mpz_t N)
{
	mpz_t A,B,C,D,E,G,F,H;

	mpz_inits(A,B,C,D,E,G,F,H,NULL);

	mpz_pow_ui(A, P->X, 2);
	mpz_mod(A, A, N);

	mpz_pow_ui(B, P->Y, 2);
	mpz_mod(B, B, N);

	mpz_pow_ui(C, P->Z, 2);
	mpz_mod(C, C, N);
	mpz_mul_ui(C, C, 2);
	mpz_mod(C, C, N);
	
	mpz_mul_si(D, A, -1);

	mpz_add(E, P->X, P->Y);
	mpz_pow_ui(E, E, 2);
	mpz_mod(E, E, N);
	mpz_sub(E, E, A);
	mpz_sub(E, E, B);

	mpz_add(G, D, B);
	
	mpz_sub(F, G, C);
	
	mpz_sub(H, D, B);

	mpz_mul_mod(R->X, E, F, N);

	mpz_mul_mod(R->Y, G, H, N);
	
	mpz_mul_mod(R->T, E, H, N);

	mpz_mul_mod(R->Z, F, G, N);

	mpz_clears(A,B,C,D,E,G,F,H,NULL);
}

void montgomery_double (MONTGOMERY_POINT R, const MONTGOMERY_POINT P, const mpz_t N,const mpz_t a)
{//montgomery曲線の2倍算の関数
  mpz_t A,B,C,D,E,F,G,four,inv;
  mpz_inits(A,B,C,D,E,F,G,four,inv,NULL);
  mpz_set_ui(four,4);
  //A = P->X + P->Z = x1+z1
  mpz_add(A,P->X,P->Z);
  //A = (A * A) mod N = (x1+z1)^2
  mpz_mul_mod(A,A,A,N);
  //B = P->X - P->Z = x1-z1
  mpz_sub(B,P->X,P->Z);
  //B = (B * B) mod N = (x1-z1)^2
  mpz_mul_mod(B,B,B,N);
  //C = A - B = (x1+z1)^2-(x1-z1)^2
  mpz_sub(C,A,B);
  //R->X = (A * B) mod N = (x1+z1)^2*(x1-z1)^2
  mpz_mul_mod(R->X,A,B,N);
  //D = a + 2 
  mpz_add_ui(D,a,2);
  //E = D / 4 = (a+2)/4
  mpz_invert(inv,four,N);
  mpz_mul(E,D,inv);
  //F = E * C = {(a+2)/4}*{(x1+z1)^2-(x1-z1)^2}
  mpz_mul_mod(F,E,C,N);
  //G = F + B = {(a+2)/4}*{(x1+z1)^2-(x1-z1)^2}+(x1-z1)^2
  mpz_add(G,F,B);
  //R->Z = C * G = {(x1+z1)^2-(x1-z1)^2}*[{(a+2)/4}*{(x1+z1)^2-(x1-z1)^2}+(x1-z1)^2]
  mpz_mul_mod(R->Z,C,G,N);
  mpz_clears(A,B,C,D,E,F,G,four,inv,NULL);
}
