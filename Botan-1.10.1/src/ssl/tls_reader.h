/*
* TLS Data Reader
* (C) 2010 Jack Lloyd
*
* Released under the terms of the Botan license
*/

#ifndef BOTAN_TLS_READER_H__
#define BOTAN_TLS_READER_H__

#include <botan/secmem.h>
#include <botan/loadstor.h>

namespace Botan {

/**
* Helper class for decoding TLS protocol messages
*/
class TLS_Data_Reader
   {
   public:
      TLS_Data_Reader(const MemoryRegion<byte>& buf_in) :
         buf(buf_in), offset(0) {}

      size_t remaining_bytes() const
         {
         return buf.size() - offset;
         }

      bool has_remaining() const
         {
         return (remaining_bytes() > 0);
         }

      void discard_next(size_t bytes)
         {
         assert_at_least(bytes);
         offset += bytes;
         }

      u16bit get_u16bit()
         {
         assert_at_least(2);
         u16bit result = make_u16bit(buf[offset], buf[offset+1]);
         offset += 2;
         return result;
         }

      byte get_byte()
         {
         assert_at_least(1);
         byte result = buf[offset];
         offset += 1;
         return result;
         }

      template<typename T, typename Container>
      Container get_elem(size_t num_elems)
         {
         assert_at_least(num_elems * sizeof(T));

         Container result(num_elems);

         for(size_t i = 0; i != num_elems; ++i)
            result[i] = load_be<T>(&buf[offset], i);

         offset += num_elems * sizeof(T);

         return result;
         }

      template<typename T>
      SecureVector<T> get_range(size_t len_bytes,
                                size_t min_elems,
                                size_t max_elems)
         {
         const size_t num_elems =
            get_num_elems(len_bytes, sizeof(T), min_elems, max_elems);

         return get_elem<T, SecureVector<T> >(num_elems);
         }

      template<typename T>
      std::vector<T> get_range_vector(size_t len_bytes,
                                      size_t min_elems,
                                      size_t max_elems)
         {
         const size_t num_elems =
            get_num_elems(len_bytes, sizeof(T), min_elems, max_elems);

         return get_elem<T, std::vector<T> >(num_elems);
         }

      template<typename T>
      SecureVector<T> get_fixed(size_t size)
         {
         return get_elem<T, SecureVector<T> >(size);
         }

   private:
      size_t get_length_field(size_t len_bytes)
         {
         assert_at_least(len_bytes);

         if(len_bytes == 1)
            return get_byte();
         else if(len_bytes == 2)
            return get_u16bit();

         throw Decoding_Error("TLS_Data_Reader: Bad length size");
         }

      size_t get_num_elems(size_t len_bytes,
                           size_t T_size,
                           size_t min_elems,
                           size_t max_elems)
         {
         const size_t byte_length = get_length_field(len_bytes);

         if(byte_length % T_size != 0)
            throw Decoding_Error("TLS_Data_Reader: Size isn't multiple of T");

         const size_t num_elems = byte_length / T_size;

         if(num_elems < min_elems || num_elems > max_elems)
            throw Decoding_Error("TLS_Data_Reader: Range outside paramaters");

         return num_elems;
         }

      void assert_at_least(size_t n) const
         {
         if(buf.size() - offset < n)
            throw Decoding_Error("TLS_Data_Reader: Corrupt packet");
         }

      const MemoryRegion<byte>& buf;
      size_t offset;
   };

/**
* Helper function for encoding length-tagged vectors
*/
template<typename T>
void append_tls_length_value(MemoryRegion<byte>& buf,
                             const T* vals,
                             size_t vals_size,
                             size_t tag_size)
   {
   const size_t T_size = sizeof(T);
   const size_t val_bytes = T_size * vals_size;

   if(tag_size != 1 && tag_size != 2)
      throw std::invalid_argument("append_tls_length_value: invalid tag size");

   if((tag_size == 1 && val_bytes > 255) ||
      (tag_size == 2 && val_bytes > 65535))
      throw std::invalid_argument("append_tls_length_value: value too large");

   for(size_t i = 0; i != tag_size; ++i)
      buf.push_back(get_byte(sizeof(val_bytes)-tag_size+i, val_bytes));

   for(size_t i = 0; i != vals_size; ++i)
      for(size_t j = 0; j != T_size; ++j)
         buf.push_back(get_byte(j, vals[i]));
   }

template<typename T>
void append_tls_length_value(MemoryRegion<byte>& buf,
                             const MemoryRegion<T>& vals,
                             size_t tag_size)
   {
   append_tls_length_value(buf, &vals[0], vals.size(), tag_size);
   }

template<typename T>
void append_tls_length_value(MemoryRegion<byte>& buf,
                             const std::vector<T>& vals,
                             size_t tag_size)
   {
   append_tls_length_value(buf, &vals[0], vals.size(), tag_size);
   }

}

#endif
