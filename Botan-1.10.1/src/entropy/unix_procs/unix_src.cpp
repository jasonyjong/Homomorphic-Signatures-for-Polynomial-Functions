/*
* Program List for Unix_EntropySource
* (C) 1999-2007 Jack Lloyd
*
* Distributed under the terms of the Botan license
*/

#include <botan/internal/es_unix.h>

namespace Botan {

/**
* Default Commands for Entropy Gathering
*/
std::vector<Unix_Program> Unix_EntropySource::get_default_sources()
   {
   std::vector<Unix_Program> srcs;

   srcs.push_back(Unix_Program("netstat -in",           1));
   srcs.push_back(Unix_Program("pfstat",                1));
   srcs.push_back(Unix_Program("vmstat -s",             1));
   srcs.push_back(Unix_Program("vmstat",                1));

   srcs.push_back(Unix_Program("arp -a -n",             2));
   srcs.push_back(Unix_Program("ifconfig -a",           2));
   srcs.push_back(Unix_Program("iostat",                2));
   srcs.push_back(Unix_Program("ipcs -a",               2));
   srcs.push_back(Unix_Program("mpstat",                2));
   srcs.push_back(Unix_Program("netstat -an",           2));
   srcs.push_back(Unix_Program("netstat -s",            2));
   srcs.push_back(Unix_Program("nfsstat",               2));
   srcs.push_back(Unix_Program("portstat",              2));
   srcs.push_back(Unix_Program("procinfo -a",           2));
   srcs.push_back(Unix_Program("pstat -T",              2));
   srcs.push_back(Unix_Program("pstat -s",              2));
   srcs.push_back(Unix_Program("uname -a",              2));
   srcs.push_back(Unix_Program("uptime",                2));

   srcs.push_back(Unix_Program("listarea",              3));
   srcs.push_back(Unix_Program("listdev",               3));
   srcs.push_back(Unix_Program("ps -A",                 3));
   srcs.push_back(Unix_Program("sysinfo",               3));

   srcs.push_back(Unix_Program("finger",                4));
   srcs.push_back(Unix_Program("mailstats",             4));
   srcs.push_back(Unix_Program("rpcinfo -p localhost",  4));
   srcs.push_back(Unix_Program("who",                   4));

   srcs.push_back(Unix_Program("df -l",                 4));
   srcs.push_back(Unix_Program("dmesg",                 4));
   srcs.push_back(Unix_Program("last -5",               4));
   srcs.push_back(Unix_Program("ls -alni /proc",        4));
   srcs.push_back(Unix_Program("ls -alni /tmp",         4));
   srcs.push_back(Unix_Program("pstat -f",              4));

   srcs.push_back(Unix_Program("ps -elf",               5));
   srcs.push_back(Unix_Program("ps aux",                5));

   srcs.push_back(Unix_Program("lsof -n",               6));
   srcs.push_back(Unix_Program("sar -A",                6));

   return srcs;
   }

}
