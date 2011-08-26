/*
* EAC SIGNED Object
* (C) 1999-2010 Jack Lloyd
*     2007 FlexSecure GmbH
*
* Distributed under the terms of the Botan license
*/

#include <botan/signed_obj.h>
#include <botan/pubkey.h>
#include <botan/oids.h>
#include <memory>

namespace Botan {

/*
* Return a BER encoded X.509 object
*/
SecureVector<byte> EAC_Signed_Object::BER_encode() const
   {
   Pipe ber;
   ber.start_msg();
   encode(ber, RAW_BER);
   ber.end_msg();
   return ber.read_all();
   }

/*
* Return a PEM encoded X.509 object
*/
std::string EAC_Signed_Object::PEM_encode() const
   {
   Pipe pem;
   pem.start_msg();
   encode(pem, PEM);
   pem.end_msg();
   return pem.read_all_as_string();
   }

/*
* Return the algorithm used to sign this object
*/
AlgorithmIdentifier EAC_Signed_Object::signature_algorithm() const
   {
   return sig_algo;
   }

bool EAC_Signed_Object::check_signature(Public_Key& pub_key,
                                        const MemoryRegion<byte>& sig) const
   {
   try
      {
      std::vector<std::string> sig_info =
         split_on(OIDS::lookup(sig_algo.oid), '/');

      if(sig_info.size() != 2 || sig_info[0] != pub_key.algo_name())
         {
         return false;
         }

      std::string padding = sig_info[1];
      Signature_Format format =
         (pub_key.message_parts() >= 2) ? DER_SEQUENCE : IEEE_1363;

      SecureVector<byte> to_sign = tbs_data();

      PK_Verifier verifier(pub_key, padding, format);
      return verifier.verify_message(to_sign, sig);
      }
   catch(...)
      {
      return false;
      }
   }

/*
* Try to decode the actual information
*/
void EAC_Signed_Object::do_decode()
   {
   try {
      force_decode();
   }
   catch(Decoding_Error& e)
      {
      const std::string what = e.what();
      throw Decoding_Error(PEM_label_pref + " decoding failed (" + what + ")");
      }
   catch(Invalid_Argument& e)
      {
      const std::string what = e.what();
      throw Decoding_Error(PEM_label_pref + " decoding failed (" + what + ")");
      }
   }

}
