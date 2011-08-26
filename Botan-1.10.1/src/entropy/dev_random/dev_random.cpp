/*
* /dev/random EntropySource
* (C) 1999-2009 Jack Lloyd
*
* Distributed under the terms of the Botan license
*/

#include <botan/internal/dev_random.h>

#include <sys/types.h>
#include <sys/select.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>

namespace Botan {

/**
Close the device, if open
*/
void Device_EntropySource::Device_Reader::close()
   {
   if(fd > 0) { ::close(fd); fd = -1; }
   }

/**
Read bytes from a device file
*/
size_t Device_EntropySource::Device_Reader::get(byte out[], size_t length,
                                                size_t ms_wait_time)
   {
   if(fd < 0)
      return 0;

   if(fd >= FD_SETSIZE)
      return 0;

   fd_set read_set;
   FD_ZERO(&read_set);
   FD_SET(fd, &read_set);

   struct ::timeval timeout;

   timeout.tv_sec = (ms_wait_time / 1000);
   timeout.tv_usec = (ms_wait_time % 1000) * 1000;

   if(::select(fd + 1, &read_set, 0, 0, &timeout) < 0)
      return 0;

   if(!(FD_ISSET(fd, &read_set)))
      return 0;

   const ssize_t got = ::read(fd, out, length);
   if(got <= 0)
      return 0;

   return static_cast<size_t>(got);
   }

/**
Attempt to open a device
*/
Device_EntropySource::Device_Reader::fd_type
Device_EntropySource::Device_Reader::open(const std::string& pathname)
   {
#ifndef O_NONBLOCK
  #define O_NONBLOCK 0
#endif

#ifndef O_NOCTTY
  #define O_NOCTTY 0
#endif

   const int flags = O_RDONLY | O_NONBLOCK | O_NOCTTY;
   return ::open(pathname.c_str(), flags);
   }

/**
Device_EntropySource constructor
Open a file descriptor to each (available) device in fsnames
*/
Device_EntropySource::Device_EntropySource(
   const std::vector<std::string>& fsnames)
   {
   for(size_t i = 0; i != fsnames.size(); ++i)
      {
      Device_Reader::fd_type fd = Device_Reader::open(fsnames[i]);
      if(fd > 0)
         devices.push_back(Device_Reader(fd));
      }
   }

/**
Device_EntropySource destructor: close all open devices
*/
Device_EntropySource::~Device_EntropySource()
   {
   for(size_t i = 0; i != devices.size(); ++i)
      devices[i].close();
   }

/**
* Gather entropy from a RNG device
*/
void Device_EntropySource::poll(Entropy_Accumulator& accum)
   {
   size_t go_get = std::min<size_t>(accum.desired_remaining_bits() / 8, 48);

   size_t read_wait_ms = std::max<size_t>(go_get, 1000);
   MemoryRegion<byte>& io_buffer = accum.get_io_buffer(go_get);

   for(size_t i = 0; i != devices.size(); ++i)
      {
      size_t got = devices[i].get(&io_buffer[0], io_buffer.size(),
                                  read_wait_ms);

      if(got)
         {
         accum.add(&io_buffer[0], got, 8);
         break;
         }
      }
   }

}
