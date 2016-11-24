/* ================================================================= *
 *  Util.h: Header file with some utility functions                  *
 *                                                                   *
 *  HISEA: An efficient tool for long read alignment                 *
 *                                                                   *
 *  Copyright (c) 2016, Nilesh Khiste                                *
 *  All rights reserved                                              *
 *                                                                   * 
 *  This program is free software: you can redistribute it and/or    *
 *  modify it under the terms of the GNU General Public License as   *
 *  published by the Free Software Foundation, either version 3 of   *
 *  the License, or (at your option) any later version.              *
 *                                                                   *
 *  This program is distributed in the hope that it will be useful,  *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 *  GNU General Public License for more details.                     *
 *                                                                   *
 *  You should have received a copy of the GNU General Public        *
 *  License along with this program.                                 *
 *                                                                   *
 *  This file is subject to the terms and conditions defined in the  *
 *  file 'LICENSE', which is part of this source code package.       *
 * ================================================================= */

#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <ctime>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <cassert>
#include <new>
#include <list>
#include <iomanip>
#include <unordered_map>
#include <sys/resource.h>
#include <dirent.h>
#include <MemoryTracker.h>
#include <Globals.h>
#include <omp.h>

inline
long int  Max  (long int A, long int B)

//  Return the larger of  A  and  B .

{
   if  (A < B)
       return  B;
     else
       return  A;
}

inline void printMemUsage()
{
    long long vmUsage, rssUsage, vmPeak, rssPeak;
    if(MemoryTracker::getOSMemUsage(vmUsage, rssUsage, vmPeak, rssPeak))
    {       
            std::cerr<< "[MEMORY]  virtual usage: " << vmUsage/1024 << " MB\n";                    
            std::cerr<< "[MEMORY] resident usage: " << rssUsage/1024 << " MB\n";           
            std::cerr<< "[MEMORY]   virtual peak: " << vmPeak/1024 << " MB\n";                     
            std::cerr<< "[MEMORY]  resident peak: " << rssPeak/1024 << " MB\n\n";          
    }
}

/* input starts from 1, while output starts from 0
 * Output is number of uint64_t needed to store the sequence
 */
inline size_t char2byteSize(size_t size)
{
    return (size-1)/CHARACTERS_IN_UINT; //do a +1 for allocation (counting starts with 0)
}

/* input starts from 1, while output starts from 0
 * Output is the offset of the last character in the sequence
 */
inline size_t size2Offset(size_t size)
{
    return (size-1)%CHARACTERS_IN_UINT; 
}

/* compute and store reverse complement of all 8 basepair (16bit) DNA strings */
inline void preComputeReverseComplement()
{
    uint16_t complement=0;
    uint16_t charCount=0;

    for (int i=0; i<65536; i++)
    {
        charCount=0;
        complement=i;
        complement=~complement;

        while(charCount < 8)
        {
            if(charCount < 4)
               preComputedReverseComplement[i] |= ((complement & reverse_mask[charCount]) << (16-4*charCount-2));
           else
               preComputedReverseComplement[i] |= ((complement & reverse_mask[charCount]) >> (4*(charCount%4)+2));
            charCount++;
        }
    }
}

inline std::string trim(const std::string& str,
                 const std::string& whitespace = " ")
{
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

inline void Usage()
{
   std::cerr <<
           "USAGE: hisea  [--self] [--kmerLen  <int>] [--minStoreLen <int>] [--minOlapLen <int>] --ref <file/directory> --query <file/directory>" << std::endl
       <<  std::endl
       <<  "For full list of options type \"hisea --help\"" << std::endl
       <<  std::endl;
}

#endif
