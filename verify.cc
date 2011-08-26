#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <stdio.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_p.h>
#include <NTL/vec_ZZ.h>
#include <NTL/RR.h>
#include <botan/botan.h>
#include <botan/skein_512.h>
#include <botan/buf_comp.h>

NTL_CLIENT
using namespace std;

ZZ getL(ZZ k, ZZ d);
RR GetMagnitude(vec_ZZ sig);
void VerifyMagnitude(vec_ZZ, ZZ a, ZZ_p m, ZZ p, ZZ b, ZZ q);
ZZ getSum(queue<ZZ>& coe, vector<ZZ> pows, vector<ZZ> args, ZZ deg, ZZ k, ZZ d, ZZ tot, int front, int back);

int main(){
  try{
    Botan::LibraryInitializer init;
  }catch(exception& e){
    cerr <<e.what() <<"\n";
  }
 
  ifstream filestream;
  string file;
  long n;
  vec_ZZ s;
  ZZ p, q, a, b, nu, y, d, k;
  vector<vec_ZZ> sigma;
   
  cout<<"Please enter the public key file"<<endl;
  cin>>file;
  filestream.open(file.c_str());

  filestream>> n >> p >> q >> a >> b >> nu >> y >> d;
  filestream.close();

  vec_ZZ sigmaPrime;
  cout<<"Please enter the file that contains the evaluated signature:"<<endl;
  cin>>file;
  filestream.open(file.c_str());
  filestream>>k;
  filestream>>sigmaPrime;
  filestream.close();


  cout<<"Please type in file with encoded function: "<<endl;
  cin >>file;
  filestream.open(file.c_str(), fstream::in);

  ZZ l;
  ZZ coef;
  l = getL(k, d);
  queue<ZZ> coe;
  for(int i=0; i<l; i++){
    filestream>>coef;
    coe.push(coef);
  }
  filestream.close();
  queue<ZZ> coefCopy = coe;

  //evaulate sigma(b) mod q
  ZZ sigPrimeAns;
  ZZ sigPrime_A;
  ZZ sigPrime_B;
  sigPrimeAns = 0;
  for(int i=0; i<sigmaPrime.length(); i++){
    sigPrimeAns += sigmaPrime[i]*power(b, i); 
  }
  sigPrime_B =  sigPrimeAns%q;

  //evaluate sigma(a) mod p 
  sigPrimeAns = 0;
  for(int i=0; i<sigmaPrime.length(); i++){
    sigPrimeAns += sigmaPrime[i]*power(a, i); 
  }
  sigPrime_A = sigPrimeAns%p;



//need it to get the tag name
  cout<<"Please enter the file with data (to get the tag): "<<endl;
  cin >>file;
  filestream.open(file.c_str(), fstream::in);
  string tag;
  ZZ mes;
  filestream >>tag;


  string hash = "Skein-512";
  Botan::Pipe pipe(new Botan::Hash_Filter(hash), new Botan::Hex_Encoder, new Botan::Hex_Encoder);
  vector<ZZ> alpha;
  vector<ZZ> messages;

//again, retrieves and hashes all "grades" + i
  for(int i=0; i<k; i++){
    pipe.start_msg();
    string m = tag;
    ostringstream sin;
    sin <<i;
    string index = sin.str();
    m+=index;

    pipe.write(m);
    pipe.end_msg();

    string m1 = pipe.read_all_as_string(i);
    char* m2 = (char*) m1.c_str();

    ZZ x = to_ZZ(m2);
    x = x%q;
    alpha.push_back(x);
    filestream>>mes;

    messages.push_back(mes);

  }
  filestream.close();

  ZZ omegaFunc;
  ZZ finalM;
  finalM = 0;
  omegaFunc = 0;
  ZZ deg;
  deg =0;

//evaluates encoded function
  while(deg < d+1){
    vector<ZZ> pows;
    vector<ZZ> pows_2;
    pows.push_back(deg);
    pows_2.push_back(deg);
    for(int i=0; i<k-1; i++){
      pows.push_back(to_ZZ(0));
      pows_2.push_back(to_ZZ(0));
    }
    omegaFunc+=getSum(coe, pows, alpha, deg, k, d, to_ZZ(0), 0, 0);
    finalM+=getSum(coefCopy, pows_2, messages, deg, k, d, to_ZZ(0), 0, 0);
    deg++;
  }
  omegaFunc = omegaFunc%q;
  finalM = finalM%p;

      
  //test one: sigma_prime(a) mod p = m
  if(sigPrime_A != finalM){
    cout<<"NOT AUTHENTICATED: Final Message does not match! "<<endl;
    exit(1);
  }
  //test two: sigma_prime(b) mo q = omega(f)
  else if(sigPrime_B != omegaFunc){
    cout<<"NOT AUTHENTICATED: Omega Function does not match! "<<endl;
    exit(1);
  }
  else cout<<"AUTHENTIC"<<endl;

  return 1;
}

//varifies that the function returns correct result
//return magnitude of signature
RR GetMagnitude(vec_ZZ sig){
  RR tot;
  tot=0;
  for(int i=0; i<sig.length(); i++) tot+=to_RR(sig[i]*sig[i]);

  return sqrt(tot);
}

//Sums up the encoded function, same as in evaluate  
ZZ getSum(queue<ZZ>& coe, vector<ZZ> pows, vector<ZZ> args, ZZ deg, ZZ k, ZZ d, ZZ tot, int front, int back){

  if(back++==k) return to_ZZ(0);
  
  ZZ final;
  final=1;
  ZZ coefficient;
  coefficient = coe.front();
  coe.pop();
  if(coefficient!=0){
    for(int i=0; i<k; i++){
      final*=power(args[i],to_long(pows[i]));
    }
  }
  
  tot=final*coefficient;

  while(back!=k&&pows[front]!=0){
    pows[front]--;
    pows[back]++;
    tot+=getSum(coe, pows, args, deg, k, d, tot, front+1, back);
  }
  return tot;
}

//Returns combination (k+d , d) which is # of terms in encoded fn.
ZZ getL(ZZ k, ZZ d){
  ZZ l; 
  l = 1;
  ZZ top;
  top = k+d;

  while(d<top){
    l*=top; 
    top-=1;
  }

  top = k;
  while(top>1){
    l=l/top;
    top-=1;
  }
  return l;
}
