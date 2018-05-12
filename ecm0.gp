
\\超特異曲線によるWeilペアリング

WeilPairing(E,P,Q,m,f)={
	local(PP,QQ,T);

	PP=P*Mod(1,f);
	QQ=[-Q[1]*Mod(1,f),Q[2]*Mod(x,f)];\\distortion mapの処理

	T=ellweilpairing(E,PP,QQ,m);

	return(T);
}


private_key_create(limit)=
{
  local(private_key);
  private_key=random(limit);
  return(private_key);
}

public_key_create(private_key,Ell,P)=
{
  local(public_key);
  public_key=matrix(2,2);
  public_key[1,]=P;
  public_key[2,]=ellmul(Ell,P,private_key);
  return(public_key);
}

hash1(keyword)=
{
  return(keyword);
}

hash2(a)=
{
  local(k);
  k=123;
  return(k*lift(a));
}

trapdoor_create(private_key,keyword,Ell)=
{
  local(trapdoor);
  trapdoor=ellmul(Ell,hash1(keyword),private_key);
  return(trapdoor);
}
keyword_encrypt(keyword,public_key,Ell,Ell_pair,p,f,limit)=
{
  local(r,A,B,pairing,m,peks);
  r=random(limit);
  A=ellmul(Ell,public_key[1,],r);

  m=(p-1)/2+1;
  pairing=WeilPairing(Ell_pair,hash1(keyword),ellmul(Ell,public_key[2,],r),m,f);
  B=hash2(pairing);
  peks=matrix(2,2);
  peks[1,]=A;
  peks[2,1]=B;
  return(peks);
}

test(peks,trapdoor,Ell,Ell_pair,p,f)=
{
  local(check,m);
  m=(p-1)/2+1;
	\\print("test : ",WeilPairing(Ell_pair,trapdoor,peks[1,],m,f));
  check=hash2(WeilPairing(Ell_pair,trapdoor,peks[1,],m,f));
  if(check-peks[2,1],
    \\print("search_fail");
		return(0),
		\\print("search_success");
		return(1)
  );
}

{
  limit=2^255;

  check=0;
	while(check==0,
		q=random(2^255);
		q=nextprime(q);
		p=2*q+1;
		if(p>2^255 && isprime(p)==1 && p%4==3,check=1);
	);

  E=ellinit([0,0,0,1,0],p);
  f=(x^2+1)*Mod(1,p);
  ff=ffgen(f);
  E_p=ellinit([0,0,0,1,0],ff);
  P=random(E);
	P=ellmul(E,P,2);

  print("zentei");
  print("p : ",p);
  print("E : y^2=x^3+x");
  print("P : ",P);

  private_key=private_key_create(limit);
  public_key=matrix(2,2);
  public_key=public_key_create(private_key,E,P);

  print("private_key : ",private_key);
  print(" public_key : P = ",public_key[1,]);
  print("            a*P = ",public_key[2,]);
	print("-------------------------------------------------");

	n=12;

  keyword=matrix(n,2);
  trapdoor=matrix(n,2);
  for(i=1,n,
		keyword[i,]=random(E);
		keyword[i,]=ellmul(E,keyword[i,],2);
    trapdoor[i,]=trapdoor_create(private_key,keyword[i,],E);
    \\print("keyword  : ",keyword[i,]);
    \\print("trapdoor : ",trapdoor[i,]);
  );

	\\検索のための暗号化、そして検索
	data=matrix(n,n);
  peks=matrix(2,2);
  for(i=1,n,
    peks=keyword_encrypt(keyword[i,],public_key,E,E_p,p,f,limit);
    \\print(peks);
		for(j=1,n,
			data[i,j]=-1;
    	data[i,j]=test(peks,trapdoor[j,],E,E_p,p,f);
		);
  );
	data
}
