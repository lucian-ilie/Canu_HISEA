/* ================================================================= *
 *  Options.h: Header file for command line options processing       *
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

#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>
#include <getopt.h>
#include <Util.h>

class Options {

    public:

        Options(); 
        int& getMaxKmerHits();
        void setMaxKmerHits(char *);
        void setMinOlapLength(char *);
        void setMinStoreLength(char *);
        void setMinMatches(char *);
        void setFilterCount(char *);
        int& getFilterCount();
        void setMaxShift(char *);
        float& getMaxShift();
        void setStartSeedIndex(char *);
        int& getSmerSeedIndex();
        void setSmallKmerIndex(char *);
        int& getSmallKmerIndex();
        void setErrorRate(char *);
        float& getErrorRate();
        void setThreads(char *);
        int& getNumThreads();
        void setFileName(char *, FileType);
        std::string& getFileName(int type);
        std::fstream& getOutputStream();
        void checkOptions();
        void checkValidDir(const char *str);
        void help();
        void addFiles(char *str, FileType type);
        std::vector <std::pair <std::string, FileType> > vec_pair;
        inline int& getMinOlapLength()
        {
            return minOlapLength;
        }
        inline int& getMinStoreLength()
        {
            return minStoreLength;
        }
        inline int& getMinMatches()
        {
            return minMatches;
        }

        std::string indexFile;
        std::string queryFile;
    private:

        std::string ignoreFile;
        std::string refFile;
        int minOlapLength;
        int minStoreLength;
        int minMatches;
        int maxKmerHits;
        int startSmerIndex;
        int filterCount;
        int threads;
        float errorRate;
        float maxShift;
        int smallKmerIndex;
};

#endif
