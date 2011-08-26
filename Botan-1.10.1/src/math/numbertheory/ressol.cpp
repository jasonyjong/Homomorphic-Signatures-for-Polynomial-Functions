/*
* Shanks-Tonnelli (RESSOL)
* (C) 2007-2008 Falko Strenzke, FlexSecure GmbH
* (C) 2008 Jack Lloyd
*
* Distributed under the terms of the Botan license
*/

#include <botan/numthry.h>
#include <botan/reducer.h>

namespace Botan {

/*
* Shanks-Tonnelli algorithm
*/
BigInt ressol(const BigInt& a, const BigInt& p)
   {
   if(a < 0)
      throw Invalid_Argument("ressol(): a to solve for must be positive");
   if(p <= 1)
      throw Invalid_Argument("ressol(): prime must be > 1");

   if(a == 0)
      return 0;
   if(p == 2)
      return a;

   if(jacobi(a, p) != 1) // not a quadratic residue
      return -BigInt(1);

   if(p % 4 == 3)
      return power_mod(a, ((p+1) >> 2), p);

   size_t s = low_zero_bits(p - 1);
   BigInt q = p >> s;

   q -= 1;
   q >>= 1;

   Modular_Reducer mod_p(p);

   BigInt r = power_mod(a, q, p);
   BigInt n = mod_p.multiply(a, mod_p.square(r));
   r = mod_p.multiply(r, a);

   if(n == 1)
      return r;

   // find random non quadratic residue z
   BigInt z = 2;
   while(jacobi(z, p) == 1) // while z quadratic residue
      ++z;

   BigInt c = power_mod(z, (q << 1) + 1, p);

   while(n > 1)
      {
      q = n;

      size_t i = 0;
      while(q != 1)
         {
         q = mod_p.square(q);
         ++i;
         }

      if(s <= i)
         return -BigInt(1);

      c = power_mod(c, BigInt(BigInt::Power2, s-i-1), p);
      r = mod_p.multiply(r, c);
      c = mod_p.square(c);
      n = mod_p.multiply(n, c);
      s = i;
      }

   return r;
   }

}
