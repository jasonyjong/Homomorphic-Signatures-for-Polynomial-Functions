/*
* PKCS1 EME
* (C) 1999-2007 Jack Lloyd
*
* Distributed under the terms of the Botan license
*/

#include <botan/eme_pkcs.h>

namespace Botan {

/*
* PKCS1 Pad Operation
*/
SecureVector<byte> EME_PKCS1v15::pad(const byte in[], size_t inlen,
                                     size_t olen,
                                     RandomNumberGenerator& rng) const
   {
   olen /= 8;

   if(olen < 10)
      throw Encoding_Error("PKCS1: Output space too small");
   if(inlen > olen - 10)
      throw Encoding_Error("PKCS1: Input is too large");

   SecureVector<byte> out(olen);

   out[0] = 0x02;
   for(size_t j = 1; j != olen - inlen - 1; ++j)
      while(out[j] == 0)
         out[j] = rng.next_byte();
   out.copy(olen - inlen, in, inlen);

   return out;
   }

/*
* PKCS1 Unpad Operation
*/
SecureVector<byte> EME_PKCS1v15::unpad(const byte in[], size_t inlen,
                                       size_t key_len) const
   {
   if(inlen != key_len / 8 || inlen < 10 || in[0] != 0x02)
      throw Decoding_Error("PKCS1::unpad");

   size_t seperator = 0;
   for(size_t j = 0; j != inlen; ++j)
      if(in[j] == 0)
         {
         seperator = j;
         break;
         }
   if(seperator < 9)
      throw Decoding_Error("PKCS1::unpad");

   return SecureVector<byte>(in + seperator + 1, inlen - seperator - 1);
   }

/*
* Return the max input size for a given key size
*/
size_t EME_PKCS1v15::maximum_input_size(size_t keybits) const
   {
   if(keybits / 8 > 10)
      return ((keybits / 8) - 10);
   else
      return 0;
   }

}
