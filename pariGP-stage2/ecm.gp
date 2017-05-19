ecm(N, B1,B2,stage2, nb = 100)=
{
  local(k,p,s,U,r,v,v2,u,d,Gx,Hx,G,R,Ecl,mP);
  k=lcm(2,3);
  r=210;
  for(i=4,B1,k=lcm(i,k));

  for(a = 1, nb,


    print("stage1: a= ",a);

    Ecl=ellinit([a,1]*Mod(1,N));

    iferr(mP=ellmul(Ecl, [0,1]*Mod(1,N), k),
    E,return(gcd(lift(component(E,2)),N)),
    errname(E)=="e_INV" && type(component(E,2)) == "t_INTMOD");


    print("stage2: a= ",a);

    U=matrix(r,2);
    for(i=1,r,
      if(gcd(r,i)>1,
        U[i,]=ellmul(Ecl,mP,i-1),U[i,]=[0,0]*Mod(1,N)
      );
    );
    s=nextprime(B1);
    v2=(s-s%r)/r;
    R=ellmul(Ecl,mP,r);
    G=ellmul(Ecl,R,v2);
    d=1;

    while(s<B2,    
      v=(s-s%r)/r;
      if(v!=v2,G=elladd(Ecl,G,R);v2=v);
      u=s%r;
  
      Gx=G[1];

      Hx=U[u+1,1];
    
      d=d*(Gx-Hx);
      if(lift(gcd(d,N))>1,return(lift(gcd(d,N))));

      s=nextprime(s+1);
    );
  );
}