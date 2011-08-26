#include<iostream>
#include<fstream>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZXFactoring.h>
#include <NTL/pair_ZZX_long.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/mat_ZZ.h>

NTL_CLIENT

struct primeIdeal{
  ZZ p;
  ZZ alpha;   //description from (x-alpha)
  ZZX generator;
};

primeIdeal PrincGen(ZZX, long);
ZZX GetIrrePoly(long n);
mat_ZZ setBasis(long n, ZZX gen_p, ZZX gen_q, ZZX f_x);

//Generate random irreducible polynomial of degree n
ZZX GetIrrePoly(long n){
  ZZX irre = ZZX(INIT_SIZE, n+1);
  
  ZZ new_n, sign_s;
  conv(new_n, n+1);
  conv(sign_s, 2); 

  vec_ZZX factors;
  long test = 1;
  while(test){
    ZZ coeff, sn;
    //randomly assigns coef + sign
    for(int i=0; i<n; i++){
      RandomBnd(sn, sign_s);
      RandomBnd(coeff, new_n);
      if(sn==1) coeff = -coeff;  //randomly assign sign
      SetCoeff(irre, i, coeff);
    }
    SetCoeff(irre, n, 1);
    factors = SFFactor(irre,0,0);
    test = (factors(1)!=irre);
  }
  return irre;
}

//outputs prime ideals using PrincGen Algorithm
primeIdeal PrincGen(ZZX irre, long n){
  ZZX s_x = ZZX(INIT_SIZE, n+1);
  ZZX g_x = ZZX(INIT_SIZE, n+1);

  double nu = pow(2.0, sqrt(n))/2;

  ZZ param_nu, sign_s;
  conv(param_nu, nu+1);
  conv(sign_s, 2); 

  ZZ coeff, sn, prime;
  
  while(true){
    for(int i=0; i<n; i++){
      RandomBnd(sn, sign_s);
      RandomBnd(coeff, param_nu);
      if(sn==1) coeff = -coeff;
      SetCoeff(s_x, i, coeff);
    }
    
    //generate g(x) = 2*s(x) +1  
    ZZ co;
    for(long i=0; i<n+1; i++) {
      GetCoeff(co, s_x, i);
      SetCoeff(g_x, i, co*2);
    }
    GetCoeff(co, s_x, 0);
    SetCoeff(g_x, 0, co*2+1);
  
    //gets resultant between F(x) and g(x)
    ZZ res = resultant(irre, g_x, 0);
    if(ProbPrime(res, 100)) {
      prime = res;
      break;
    }
  }

//convert them to be in field F_p
  ZZ_p::init(prime);
  ZZ_pX moded_g_x, moded_irre; // = ZZ_pX(n+1);

//convert to field 
  ZZ co;
  for(int i=0; i<n+1; i++){
    GetCoeff(co, g_x, i);
    SetCoeff(moded_g_x, i, to_ZZ_p(co));
  }
  for(int i=0; i<n+1; i++){
    GetCoeff(co, irre, i);
    SetCoeff(moded_irre, i, to_ZZ_p(co));
  }


//GCD will return one-degree expression: 
  ZZ_pX d_x = GCD(moded_g_x, moded_irre);

  ZZ_p coe;
  primeIdeal pi; 

  GetCoeff(coe, d_x, 0);
  pi.alpha = prime - rep(coe);
  pi.generator = g_x;
  pi.p = prime;

  return pi; 
}

//Forms Basis, which is span of generators.
mat_ZZ setBasis(long n, ZZX gen_p, ZZX gen_q, ZZX f_x){
  mat_ZZ basis;
  basis.SetDims(n-1,n-1);

 //load in each base one at a time
  ZZX mlt = gen_p*gen_q;
  ZZ coe;
  for(int i=0; i<n-1; i++){
    ZZX varBase = ZZX(INIT_SIZE, n-1);
    SetCoeff(varBase, i,1);
    ZZX bi = (mlt*varBase)%f_x;
    for(int j=0; j<n-1; j++){
      GetCoeff(coe, bi, j);
      basis[j][i] = coe;
    }
  }

  return basis; 
}

int main(){
  SetSeed(to_ZZ(time(NULL)));
  long n;
  while(true){
    cout<< "Please enter the security parameter (n): " <<endl;
    if(!(cin>>n)) {
      cout<<"Not a valud number"<<endl;
      cin.clear();
    }
    else if(n<1) cout<<"Not a positive integer"<<endl;
    else break;
  }

  //get monic irreducible poly of degree n
  ZZX irre = GetIrrePoly(n);
  
  //PrincGen Algorithm, for two ideals
  primeIdeal pi_1 = PrincGen(irre, n); 
  primeIdeal pi_2 = PrincGen(irre, n); 
  
  //get TBasis, set as secret key
  mat_ZZ tBasis = setBasis(n+1, pi_1.generator, pi_2.generator, irre);

  //get parameter nu
  long nu = pow(n,4)*log(n)/log(10);
 
  ZZ d;
  while(true){
    cout<< "Please enter d, largest degree of encoded function (d): " <<endl;
    if(!(cin>>d)) {
      cout<<"Not a valud number"<<endl;
      cin.clear();
    }
    else if(n<1) cout<<"Not a positive integer"<<endl;
    else break;
  }
  
  ZZ y;
  while(true){
    cout<< "Please enter y, range for coefficients, [-y/2 , y/2]: " <<endl;
    if(!(cin>>y)) {
      cout<<"Not a valud number"<<endl;
      cin.clear();
    }
    else if(n<1) cout<<"Not a positive integer"<<endl;
    else break;
  }

/*Public key info:
    One-Degree prime ideals as pairs:
      p, a
      q, b
    Parameter:
      nu
    Degree of encoded Fn applied to signatures:
      d
    Coefficient of Fn:
      y
  Private Key:
    Basis made from generators of prime ideals
*/
  string file;
  cout<<"Please name file to store public key."<<endl;
  cin >> file;

  ofstream myfile;
  myfile.open(file.c_str());
  myfile << n <<" "<<pi_1.p<<" "<<pi_2.p<<" "<<pi_1.alpha<<" "<<pi_2.alpha;
  myfile<<" "<<nu<<" "<<y<<" "<<d; 
  myfile.close();

  cout<<"Please name file to store private key."<<endl;
  cin >> file;

  myfile.open(file.c_str());
  myfile << tBasis;
  myfile.close();
 
  return 0;
}
