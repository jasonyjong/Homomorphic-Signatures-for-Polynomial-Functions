/*
* Zlib Compressor
* (C) 2001 Peter J Jones
*     2001-2007 Jack Lloyd
*
* Distributed under the terms of the Botan license
*/

#ifndef BOTAN_ZLIB_H__
#define BOTAN_ZLIB_H__

#include <botan/filter.h>

namespace Botan {

/**
* Zlib Compression Filter
*/
class BOTAN_DLL Zlib_Compression : public Filter
   {
   public:
      std::string name() const { return "Zlib_Compression"; }

      void write(const byte input[], size_t length);
      void start_msg();
      void end_msg();

      /**
      * Flush the compressor
      */
      void flush();

      /**
      @param level how much effort to use on compressing (0 to 9);
      higher levels are slower but tend to give better compression
      */
      Zlib_Compression(size_t level = 6);

      ~Zlib_Compression() { clear(); }
   private:
      void clear();
      const size_t level;
      SecureVector<byte> buffer;
      class Zlib_Stream* zlib;
   };

/**
* Zlib Decompression Filter
*/
class BOTAN_DLL Zlib_Decompression : public Filter
   {
   public:
      std::string name() const { return "Zlib_Decompression"; }

      void write(const byte input[], size_t length);
      void start_msg();
      void end_msg();

      Zlib_Decompression();
      ~Zlib_Decompression() { clear(); }
   private:
      void clear();
      SecureVector<byte> buffer;
      class Zlib_Stream* zlib;
      bool no_writes;
   };

}

#endif
