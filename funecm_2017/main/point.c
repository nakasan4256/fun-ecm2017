#include <stdio.h>
#include "gmp.h"
#include <time.h>
#include "point.h"

/* convert affine point to projective point */
void afftopro(PROJECTIVE_POINT R, const AFFINE_POINT P)
{
	mpz_set(R->X, P->x);
	mpz_set(R->Y, P->y);
	mpz_set_ui(R->Z, 1);
}

/* convert projective point to affine point */
void protoaff(AFFINE_POINT R, const PROJECTIVE_POINT P, const mpz_t N)
{
	mpz_t inv;
	mpz_init(inv);

	mpz_invert(inv, P->Z, N);

	mpz_mul(R->x, P->X, inv);
	mpz_mul(R->y, P->Y, inv);

	mpz_mod(R->x, R->x, N);
	mpz_mod(R->y, R->y, N);

	mpz_clear(inv);
}

/* convert projective point to extended point */
void protoext(EXTENDED_POINT R, const PROJECTIVE_POINT P, const mpz_t N)
{
	/* (X:Y:Z) -> (XZ:YZ:XY:Z^2) */
	mpz_mul(R->X, P->X, P->Z);
	mpz_mod(R->X, R->X, N);
	mpz_mul(R->Y, P->Y, P->Z);
	mpz_mod(R->Y, R->Y, N);
	mpz_mul(R->T, P->X, P->Y);
	mpz_mod(R->T, R->T, N);
	mpz_pow_ui(R->Z, P->Z, 2);
	mpz_mod(R->Z, R->Z, N);
}

/* convert extended point to projective point */
void exttopro(PROJECTIVE_POINT R, const EXTENDED_POINT P)
{
	mpz_set(R->X, P->X);
	mpz_set(R->Y, P->Y);
	mpz_set(R->Z, P->Z);
}

void protomon(MONTGOMERY_POINT R, const PROJECTIVE_POINT P, const PROJECTIVE_POINT P1, const PROJECTIVE_POINT P2, const mpz_t a, const mpz_t b, const mpz_t N)
{//projective曲線のkPの座標yをmontgomery曲線のkpの座標yに変更
  mpz_t A,B,C,D,E,F,G,H,I,J,K,L,inv;
  mpz_inits(A,B,C,D,E,F,G,H,I,J,K,L,inv,NULL);
  //A = P1->X * P->X = x1*x
  mpz_mul_mod(A,P1->X,P->X,N);
  //B = A + 1 = x1*x+1
  mpz_add_ui(B,A,1);
  //C = P1->X + P->X = x1+x
  mpz_add(C,P1->X,P->X);
  //D = a * 2 = 2a
  mpz_mul_ui(D,a,2);
  //E = C + D = x1+x+2a
  mpz_add(E,C,D);
  //F = B * E = (x1*x+1)*(x1+x+2a)
  mpz_mul_mod(F,B,E,N);
  //G = F - D = (x1*x+1)*(x1+x+2a)-2a
  mpz_sub(G,F,D);
  //H = P1->X - P->X = (x1-x)
  mpz_sub(H,P1->X,P->X);
  //H = H * H = (x1-x)^2
  mpz_mul_mod(H,H,H,N);
  //I = H * P2->X = (x1-x)^2*x2
  mpz_mul_mod(I,H,P2->X,N);
  //J = G - I = (x1*x+1)*(x1+x+2a)-2a-(x1-x)^2*x2
  mpz_sub(J,G,I);
  //K = b * 2 = 2b
  mpz_mul_ui(K,b,2);
  //L = K * P->Y = 2b*y
  mpz_mul(L,K,P->Y);
  //R->Y = J / L = (x1*x+1)*(x1+x+2a)-2a-(x1-x)^2*x2 / 2b*y
  mpz_invert(inv,L,N);
  mpz_mul_mod(R->Y,J,inv,N);
  mpz_clears(A,B,C,D,E,F,G,H,I,J,K,L,inv,NULL);
}

void exttomon(MONTGOMERY_POINT R, const EXTENDED_POINT P,const mpz_t N)
{//Edwards曲線でのX,Y座標をmontgomery曲線でのX,Z座標に変換
  mpz_t add,sub;
  mpz_inits(add,sub,NULL);
  //add = P->Y + 1 = y+1
  mpz_add_ui(add,P->Y,1);
  //sub = P->Y - 1 = y-1
  mpz_sub_ui(sub,P->Y,1);
  //sub = 1 / sub = 1/(y-1)
  mpz_invert(sub,sub,N);
  //R->X = add * sub = (y+1)*{1/(y-1)} = (y+1)/(y-1)
  mpz_mul_mod(R->X,add,sub,N);
  //R->Z = 1
  mpz_set_ui(R->Z,1);
  mpz_clears(add,sub,NULL);
}

void montgomery_coefficient (mpz_t A, mpz_t B,const mpz_t d, const mpz_t N)
{//montgomery曲線の係数A,Bをedward曲線のdから算出
  mpz_t C,D;
  mpz_inits(C,D,NULL);
  //C = d - 1
  mpz_sub_ui(C,d,1);
  //D = d + 1
  mpz_add_ui(D,d,1);
  //D = D * -1 = -d-1
  mpz_mul_ui(D,D,-1);
  //D = 1 / D = 1/(-d-1)
  mpz_invert(D,D,N);
  //C = C * 2 = 2*(d-1)
  mpz_mul_ui(C,C,2);
  //A = C * D = 2*(d-1)/(-d-1)
  mpz_mul_mod(A,C,D,N);
  //B = D * 4 = 4*(-d-1)
  mpz_mul_ui(B,D,4);
  mpz_clears(C,D,NULL);
}
 
