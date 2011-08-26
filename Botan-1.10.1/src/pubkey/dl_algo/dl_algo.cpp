/*
* DL Scheme
* (C) 1999-2007 Jack Lloyd
*
* Distributed under the terms of the Botan license
*/

#include <botan/dl_algo.h>
#include <botan/numthry.h>
#include <botan/der_enc.h>
#include <botan/ber_dec.h>

namespace Botan {

AlgorithmIdentifier DL_Scheme_PublicKey::algorithm_identifier() const
   {
   return AlgorithmIdentifier(get_oid(),
                              group.DER_encode(group_format()));
   }

MemoryVector<byte> DL_Scheme_PublicKey::x509_subject_public_key() const
   {
   return DER_Encoder().encode(y).get_contents();
   }

DL_Scheme_PublicKey::DL_Scheme_PublicKey(const AlgorithmIdentifier& alg_id,
                                         const MemoryRegion<byte>& key_bits,
                                         DL_Group::Format format)
   {
   DataSource_Memory source(alg_id.parameters);
   group.BER_decode(source, format);

   BER_Decoder(key_bits).decode(y);
   }

MemoryVector<byte> DL_Scheme_PrivateKey::pkcs8_private_key() const
   {
   return DER_Encoder().encode(x).get_contents();
   }

DL_Scheme_PrivateKey::DL_Scheme_PrivateKey(const AlgorithmIdentifier& alg_id,
                                           const MemoryRegion<byte>& key_bits,
                                           DL_Group::Format format)
   {
   DataSource_Memory source(alg_id.parameters);
   group.BER_decode(source, format);

   BER_Decoder(key_bits).decode(x);
   }

/*
* Check Public DL Parameters
*/
bool DL_Scheme_PublicKey::check_key(RandomNumberGenerator& rng,
                                    bool strong) const
   {
   if(y < 2 || y >= group_p())
      return false;
   if(!group.verify_group(rng, strong))
      return false;
   return true;
   }

/*
* Check DL Scheme Private Parameters
*/
bool DL_Scheme_PrivateKey::check_key(RandomNumberGenerator& rng,
                                     bool strong) const
   {
   const BigInt& p = group_p();
   const BigInt& g = group_g();

   if(y < 2 || y >= p || x < 2 || x >= p)
      return false;
   if(!group.verify_group(rng, strong))
      return false;

   if(!strong)
      return true;

   if(y != power_mod(g, x, p))
      return false;

   return true;
   }

}
