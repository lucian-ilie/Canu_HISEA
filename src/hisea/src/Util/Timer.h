/* ================================================================= *
 *  Timer.h: A timer class for keeping track of time during          *
 *           code execution                                          *
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

#ifndef TIMER_H
#define TIMER_H

#include <iostream>
#include <chrono>
#include <ctime>

class Timer
{
    public:
        explicit Timer (const char *str="")
        {
            t_start=std::chrono::system_clock::time_point::min();
            t_string=str;
        }
        void start()
        {
            t_start=std::chrono::system_clock::now();
        }
 
        void reset()
        {
            t_start=std::chrono::system_clock::time_point::min();
        }

        void stop()
        {
             std::cerr << "[TIME] - " << t_string << ": "<< std::chrono::duration<double>(std::chrono::high_resolution_clock::now()-t_start).count() << std::endl;
        }

    private:
        std::string t_string;
        std::chrono::system_clock::time_point t_start;
        
};

#endif
