/*
* DER Encoder
* (C) 1999-2007 Jack Lloyd
*
* Distributed under the terms of the Botan license
*/

#ifndef BOTAN_DER_ENCODER_H__
#define BOTAN_DER_ENCODER_H__

#include <botan/asn1_int.h>
#include <vector>

namespace Botan {

class BigInt;
class ASN1_Object;

/**
* General DER Encoding Object
*/
class BOTAN_DLL DER_Encoder
   {
   public:
      SecureVector<byte> get_contents();

      DER_Encoder& start_cons(ASN1_Tag, ASN1_Tag = UNIVERSAL);
      DER_Encoder& end_cons();

      DER_Encoder& start_explicit(u16bit);
      DER_Encoder& end_explicit();

      DER_Encoder& raw_bytes(const byte[], size_t);
      DER_Encoder& raw_bytes(const MemoryRegion<byte>&);

      DER_Encoder& encode_null();
      DER_Encoder& encode(bool);
      DER_Encoder& encode(size_t);
      DER_Encoder& encode(const BigInt&);
      DER_Encoder& encode(const MemoryRegion<byte>&, ASN1_Tag);
      DER_Encoder& encode(const byte[], size_t, ASN1_Tag);

      DER_Encoder& encode(bool, ASN1_Tag, ASN1_Tag = CONTEXT_SPECIFIC);
      DER_Encoder& encode(size_t, ASN1_Tag, ASN1_Tag = CONTEXT_SPECIFIC);
      DER_Encoder& encode(const BigInt&, ASN1_Tag,
                          ASN1_Tag = CONTEXT_SPECIFIC);
      DER_Encoder& encode(const MemoryRegion<byte>&, ASN1_Tag,
                          ASN1_Tag, ASN1_Tag = CONTEXT_SPECIFIC);
      DER_Encoder& encode(const byte[], size_t, ASN1_Tag,
                          ASN1_Tag, ASN1_Tag = CONTEXT_SPECIFIC);

      template<typename T>
      DER_Encoder& encode_optional(const T& value, const T& default_value)
         {
         if(value != default_value)
            encode(value);
         return (*this);
         }

      template<typename T>
      DER_Encoder& encode_list(const std::vector<T>& values)
         {
         for(size_t i = 0; i != values.size(); ++i)
            encode(values[i]);
         return (*this);
         }

      DER_Encoder& encode(const ASN1_Object&);
      DER_Encoder& encode_if(bool, DER_Encoder&);

      DER_Encoder& add_object(ASN1_Tag, ASN1_Tag, const byte[], size_t);
      DER_Encoder& add_object(ASN1_Tag, ASN1_Tag, const MemoryRegion<byte>&);
      DER_Encoder& add_object(ASN1_Tag, ASN1_Tag, const std::string&);
      DER_Encoder& add_object(ASN1_Tag, ASN1_Tag, byte);
   private:
      class DER_Sequence
         {
         public:
            ASN1_Tag tag_of() const;
            SecureVector<byte> get_contents();
            void add_bytes(const byte[], size_t);
            DER_Sequence(ASN1_Tag, ASN1_Tag);
         private:
            ASN1_Tag type_tag, class_tag;
            SecureVector<byte> contents;
            std::vector< SecureVector<byte> > set_contents;
         };
      SecureVector<byte> contents;
      std::vector<DER_Sequence> subsequences;
   };

}

#endif
