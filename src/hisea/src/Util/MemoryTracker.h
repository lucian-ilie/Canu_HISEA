/* ================================================================= *
 *  MemoryTracker.h: Header file for functions for tracking memory   *
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
#ifndef MEMORY_TRACKER_H
#define MEMORY_TRACKER_H

#include <iostream>
#include <new>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/resource.h>
#include <unistd.h>


class MemoryTracker
{
    public:
        virtual bool dummy()=0; 

        static long currMemUser; //This keeps track of direct allocation from new() calls only
        static long peakMemUser;
        
        static void getUserMemUsage();

        static bool getOSMemUsage(long long &,long long &,long long &,long long &);
};

#endif
