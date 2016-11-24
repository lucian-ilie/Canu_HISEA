/* ================================================================= *
 *  MemoryTracker.cpp: Source file for functions for tracking memory *
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

#include "MemoryTracker.h"

        void MemoryTracker::getUserMemUsage()
        {
            std::cout << "[Peak User] - Memory: " << peakMemUser << std::endl;
            std::cout << "[Current User] - Memory: " << currMemUser << std::endl;
        }

        //Borrowed from SAGE2
        bool MemoryTracker::getOSMemUsage(long long &vm_usage, long long &rss_usage, long long &peak_vm_usage, long long &peak_rss_usage)
        {
            vm_usage = -1;
            rss_usage = -1;
            peak_vm_usage = -1;
            peak_rss_usage = -1;

            std::ifstream stat_stream("/proc/self/status");
            if(stat_stream.is_open()==false)
                return false;
        
            std::string line, name;
            std::istringstream sin;
            while(vm_usage==-1 || rss_usage==-1 || peak_vm_usage==-1 || peak_rss_usage==-1)
            {
                if(!getline(stat_stream, line))
                        return false;
                sin.clear();
                sin.str(line);
                sin >> name;
                if(name == "VmPeak:")
                {
                        sin >> peak_vm_usage;
                }
                else if(name == "VmSize:")
                {
                        sin >> vm_usage;
                }
                else if(name == "VmHWM:")
                {
                        sin >> peak_rss_usage;
                }
                else if(name == "VmRSS:")
                {
                        sin >> rss_usage;
                }
            }
            stat_stream.close();
            return true;
        }

long MemoryTracker::currMemUser=0;
long MemoryTracker::peakMemUser=0;

