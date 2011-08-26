#include <string>
#include <vector> 
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <botan/botan.h>
#include <botan/skein_512.h>
#include <botan/buf_comp.h>
#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/matrix.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pE.h>
#include <NTL/LLL.h>

NTL_CLIENT
using namespace std;

ZZX GetHFunction(long n, ZZ a, ZZ b, ZZ p, ZZ q, ZZ_p alph, ZZ_p m);
ZZ SampleInteger(RR c_prime, RR s_prime, long n);
vec_ZZ SamplePre(mat_ZZ tBasis, ZZ nu, ZZX h_x);

int main(){
  try{
     Botan::LibraryInitializer init;
  }catch(exception& e){
    cerr <<e.what() <<"\n";
  }

  ifstream myfile;
  long n;
  mat_ZZ tBasis;
  ZZ p, q, a, b, nu;

  string file;
  cout<<"Please enter the name of the public key file: "<<endl;
  cin >>file;
  myfile.open(file.c_str());
  myfile>> n>> p >> q >>a >>b >> nu;
  myfile.close();
  cout<<"Please enter the name of the private key file: "<<endl;
  cin >>file;
  myfile.open(file.c_str());
  myfile >> tBasis;
  myfile.close();

  cout<<"Please enter the file with data: "<<endl;
  cin >>file;
  
  myfile.open(file.c_str());

  string tag;  
  myfile>>tag;
    
  string hash = "Skein-512";
  Botan::Pipe pipe(new Botan::Hash_Filter(hash), new Botan::Hex_Encoder, new Botan::Hex_Encoder);

  int s;
  myfile>>s;
  
  int i=0;
  vector<vec_ZZ> holder;

  //go through all data + hash tags "grade" + i, where i integers
  while(myfile.good()){
    pipe.start_msg();

    string m = tag;
    ostringstream sin;
    sin << i;
    string index = sin.str();
    m += index;

    pipe.write(m);
    pipe.end_msg();
    string m1 = pipe.read_all_as_string(i);
    char* m2 = (char*)m1.c_str(); 
    ZZ x = to_ZZ(m2);

    ZZ_p::init(q);
    ZZ_p alph = to_ZZ_p(x);
    ZZ_p::init(p);
    ZZ_p message = to_ZZ_p(to_ZZ(s));

    ZZX h_x = GetHFunction(n, a, b, p, q, alph, message);
    vec_ZZ sigma_i = SamplePre(tBasis, nu, h_x);
    holder.push_back(sigma_i);

    while(! (myfile>>s) && myfile.good()){
      cout<<"Warning: a data point is not a valid number."<<endl;
    }
   
    i++;
  }
 
  cout<<"Please enter name of file for signed data."<<endl;   
  cin >>file;
  ofstream outfile;
  outfile.open(file.c_str());
  outfile<<holder.size()<<" ";  //output k value

  //iteratively output each signature
  for(unsigned int i=0; i<holder.size(); i++) outfile<<holder[i]<<" ";

  outfile.close();
  return 1;
}

/*Samples secret key basis and releases signature
  that lies in p*q + h(x) space
 */
vec_ZZ SamplePre(mat_ZZ tBasis, ZZ nu, ZZX h_x){
  mat_RR gs;
  vec_RR c, b_gs, c_conv;
  vec_ZZ v, b;
  ComputeGS(tBasis,gs, c); 

  vec_ZZ c_p;
  c_p.SetLength(tBasis.NumRows());
  v.SetLength(tBasis.NumRows());
  for(int k=0; k<tBasis.NumRows(); k++){
    ZZ coe;
    GetCoeff(coe, -h_x, k);
    c_p[k] = coe;
  }

  b_gs.SetLength(tBasis.NumRows());
  b.SetLength(tBasis.NumRows());
  c_conv.SetLength(tBasis.NumRows());

  //does not deal with last column: zeroed out from gram schmidt
  for(int i=tBasis.NumCols()-2; i>-1; i--){
    RR s_prime, bs, ba, c_prime;
    s_prime=0;
    for(int j=0; j<tBasis.NumRows(); j++){
      b_gs[j] = gs[j][i];
      b[j] = tBasis[j][i];
      s_prime+=b_gs[j]*b_gs[j];
      c_conv[i]=to_RR(c_p[i]);
    }
      
    s_prime = to_RR(nu)/sqrt(s_prime);

    InnerProduct(bs, b_gs, c_conv);
    InnerProduct(ba, b_gs, b_gs);

    c_prime = bs/ba;

    ZZ z_i = SampleInteger(c_prime, s_prime, tBasis.NumRows());

    c_p = c_p -  b*z_i;
    v = v + b*z_i;
  }

  v = v-c_p;
  return v;
}

/*probabilitistcally generate integer from given vectors
 *with gaussian distribution
 */
ZZ SampleInteger(RR c_prime, RR s_prime, long n){
  RR param = s_prime*log10(to_RR(n));
  ZZ lowerbound = CeilToZZ(c_prime-param);
  ZZ upperbound = CeilToZZ(c_prime+param);
  ZZ diff = upperbound-lowerbound;
  while(true){
    ZZ x = FloorToZZ(random_RR()*to_RR(diff)) + lowerbound;
    RR exponent = -(3.1415926)*power((to_RR(x)-c_prime)/s_prime, 2);
    if(random_RR()<exp(exponent)) return x;
  }
}

/*Gets H function that correctly maps to desired signed data
 *Will correctly set h(a) mod p = m, h(b) mod q = alph
 */
ZZX GetHFunction(long n, ZZ a, ZZ b, ZZ p, ZZ q, ZZ_p alph, ZZ_p m){
  mat_ZZ mat;
  mat.SetDims(n,2);

  for(long i=0; i<n; i++) mat[i][0] = power(a,i);  
  for(long i=0; i<n; i++) mat[i][1] = power(b,i);  
   
//basically uses last term to add till proper mod props found
  ZZ fin=rep(alph); 
  while(fin%p!=rep(m)) fin+=q;

  ZZ d;
  vec_ZZ x;
  vec_ZZ z;
  z.SetLength(2);
  z[0]=rep(m);
  z[1]=rep(alph);
  ZZX h = ZZX(INIT_SIZE, n); 

  ZZ sign_s;
  conv(sign_s, 2);

  ZZ coeff;
  ZZ sn;
  SetCoeff(h, 0, fin);
  for(int i=1; i<n; i++){
    RandomBnd(sn, sign_s);
    RandomBnd(coeff, p);
    if(sn==1) coeff = -coeff;
    SetCoeff(h, i, 0);    
  }

  return h;
}
