/*
* DSA
* (C) 1999-2010 Jack Lloyd
*
* Distributed under the terms of the Botan license
*/

#include <botan/dsa.h>
#include <botan/numthry.h>
#include <botan/keypair.h>

namespace Botan {

/*
* DSA_PublicKey Constructor
*/
DSA_PublicKey::DSA_PublicKey(const DL_Group& grp, const BigInt& y1)
   {
   group = grp;
   y = y1;
   }

/*
* Create a DSA private key
*/
DSA_PrivateKey::DSA_PrivateKey(RandomNumberGenerator& rng,
                               const DL_Group& grp,
                               const BigInt& x_arg)
   {
   group = grp;
   x = x_arg;

   if(x == 0)
      x = BigInt::random_integer(rng, 2, group_q() - 1);

   y = power_mod(group_g(), x, group_p());

   if(x_arg == 0)
      gen_check(rng);
   else
      load_check(rng);
   }

DSA_PrivateKey::DSA_PrivateKey(const AlgorithmIdentifier& alg_id,
                               const MemoryRegion<byte>& key_bits,
                               RandomNumberGenerator& rng) :
   DL_Scheme_PrivateKey(alg_id, key_bits, DL_Group::ANSI_X9_57)
   {
   y = power_mod(group_g(), x, group_p());

   load_check(rng);
   }

/*
* Check Private DSA Parameters
*/
bool DSA_PrivateKey::check_key(RandomNumberGenerator& rng, bool strong) const
   {
   if(!DL_Scheme_PrivateKey::check_key(rng, strong) || x >= group_q())
      return false;

   if(!strong)
      return true;

   return KeyPair::signature_consistency_check(rng, *this, "EMSA1(SHA-1)");
   }

DSA_Signature_Operation::DSA_Signature_Operation(const DSA_PrivateKey& dsa) :
   q(dsa.group_q()),
   x(dsa.get_x()),
   powermod_g_p(dsa.group_g(), dsa.group_p()),
   mod_q(dsa.group_q())
   {
   }

SecureVector<byte>
DSA_Signature_Operation::sign(const byte msg[], size_t msg_len,
                              RandomNumberGenerator& rng)
   {
   rng.add_entropy(msg, msg_len);

   BigInt i(msg, msg_len);
   BigInt r = 0, s = 0;

   while(r == 0 || s == 0)
      {
      BigInt k;
      do
         k.randomize(rng, q.bits());
      while(k >= q);

      r = mod_q.reduce(powermod_g_p(k));
      s = mod_q.multiply(inverse_mod(k, q), mul_add(x, r, i));
      }

   SecureVector<byte> output(2*q.bytes());
   r.binary_encode(&output[output.size() / 2 - r.bytes()]);
   s.binary_encode(&output[output.size() - s.bytes()]);
   return output;
   }

DSA_Verification_Operation::DSA_Verification_Operation(const DSA_PublicKey& dsa) :
   q(dsa.group_q()), y(dsa.get_y())
   {
   powermod_g_p = Fixed_Base_Power_Mod(dsa.group_g(), dsa.group_p());
   powermod_y_p = Fixed_Base_Power_Mod(y, dsa.group_p());
   mod_p = Modular_Reducer(dsa.group_p());
   mod_q = Modular_Reducer(dsa.group_q());
   }

bool DSA_Verification_Operation::verify(const byte msg[], size_t msg_len,
                                        const byte sig[], size_t sig_len)
   {
   const BigInt& q = mod_q.get_modulus();

   if(sig_len != 2*q.bytes() || msg_len > q.bytes())
      return false;

   BigInt r(sig, q.bytes());
   BigInt s(sig + q.bytes(), q.bytes());
   BigInt i(msg, msg_len);

   if(r <= 0 || r >= q || s <= 0 || s >= q)
      return false;

   s = inverse_mod(s, q);
   s = mod_p.multiply(powermod_g_p(mod_q.multiply(s, i)),
                      powermod_y_p(mod_q.multiply(s, r)));

   return (mod_q.reduce(s) == r);
   }

}
