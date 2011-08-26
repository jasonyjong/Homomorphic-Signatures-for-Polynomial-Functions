#include <string>
#include <iostream>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stdio.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>

using namespace std;
NTL_CLIENT

ZZ getL(ZZ k, ZZ d);
ZZ EvalSig(vec_ZZ sig, ZZ b, ZZ q);
ZZX getSum(queue<ZZ>& coe, vector<ZZ> pows, vector<ZZX> sigma, ZZ deg, ZZ k, ZZ d, ZZX tot, int front, int back);

//Applies quotient mapping to signature
ZZ EvalSig(vec_ZZ sig, ZZ b, ZZ q){
  ZZ tot;
  tot = 0;
  for(int i=0; i<sig.length(); i++) tot+=sig[i]*power(b, i);

  return tot%q;
}

/*Plugs in the signatures to the encoded function to get signture of function.
  The recursive function goes along terms and reduces value of each until zeroed. 
  It will encounter every possible combination of a polynomial given its degree d.
*/
ZZX getSum(queue<ZZ>& coefQ, vector<ZZ> pows, vector<ZZX> sigma, ZZ deg, ZZ k, ZZ d, ZZX tot, int front, int back){
  
  if(back++==k) return to_ZZX(0);
  
  ZZX final;
  ZZ coefficient;
  final=to_ZZX(1);
  coefficient = coefQ.front();
  coefQ.pop();
  
  if(coefficient!=0){
    for(int i=0; i<k; i++){
      for(int j=0; j<pows[i];j++){
        final*=sigma[i]; 
      }
    }
  } 
 
  tot=final*coefficient;

  //advances to next term by carrying one value over to next
  while(back!=k&&pows[front]!=0){
    pows[front]--;
    pows[back]++;
    tot+=getSum(coefQ, pows, sigma, deg, k, d, tot, front+1, back);
  }
  return tot;
} 

//Calculates the number of terms needed in encoded function with combination (k+d, d)
ZZ getL(ZZ k, ZZ d){
  ZZ l, top; 
  l = 1;
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

int main(){
  ifstream filestream;
  string file;

  cout<<"Please enter public key file: "<<endl;
  cin>>file;

  filestream.open(file.c_str());
 
  long n;
  ZZ p, q, a, b, nu, y, d, k;
  vector<ZZX> sigma;

  filestream>> n>> p>> q>> a>> b>> nu>> y>> d;
  filestream.close();

  cout<<"Please enter the signature file: "<<endl;
  cin>>file;
  filestream.open(file.c_str());

  filestream>>k;
  for(int i=0; i<k; i++){
    ZZX s;
    filestream>>s;
    sigma.push_back(s);
  }
  filestream.close();

  //load in function:
  cout<<"Please type in file with encoded function: "<<endl;
  cin >>file;

  filestream.open(file.c_str());
   
  ZZ l, coef;
  l = getL(k, d);    // # of monomials needed for fn
  queue<ZZ> coe;     //store coef in queue
  for(int i=0; i<l; i++){
    filestream>>coef;
    coe.push(coef);
  }

  ZZX omegaT;
  omegaT =0;
  ZZ deg;
  deg =0;

  while(deg < d+1){
    vector<ZZ> pows;
    pows.push_back(deg);

    for(int i=0; i<k-1;i++) pows.push_back(to_ZZ(0));

    omegaT+=getSum(coe, pows, sigma, deg, k, d, to_ZZX(0), 0, 0);
    deg++;
  }
  filestream.close();
  
  ofstream outputFile; 
  cout<<"Where would you like to store evaulated signature?" <<endl;
  cin>>file;
  outputFile.open(file.c_str());
  outputFile<<k<<" ";
  outputFile<<omegaT;
  outputFile.close(); 
  return 1;
}
