/* ================================================================= *
 *  SequenceReader.h: Header file for read processing classes        *
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

#ifndef SEQ_READER_H
#define SEQ_READER_H

#include <iostream>
#include <cstdint>
#include <string>
#include <cstring>
#include <istream>
#include <ostream>
#include <fstream>
#include <algorithm>
#include <gzstream.h>
#include <Util.h>
#include <Globals.h>
#include <Options.h>

class Sequence{
    public:
    
        Sequence ();
        Sequence (uint32_t, UINT_TYPE *);

        bool getReverseComplement(); 
        bool operator < (const Sequence &seqObj) const;
        bool operator == (const Sequence &seqObj) const;        
        bool operator () (const Sequence &seqObj1, const Sequence &seqObj2) const;
        UINT_TYPE *seq;
        uint32_t size;
        
    private:
};

class SequenceReader
{
    public:
        SequenceReader(); 
        ~SequenceReader(); 
        void Initialize(Options &op); 
        void readSequences(FileType ft);
        void dupSequenceReader(SequenceReader &sr);
        uint64_t getNumberOfReads();
        uint64_t getNumberOfBP();
        uint64_t getReadArraySize();
        UINT_TYPE* getSequence(const uint64_t &, ReadDirection &);
        uint32_t getSequenceSize(const uint64_t &);
        void removeSequences(uint64_t &);
        uint64_t appendNewSequence();
        void checkReadArraySize();

    private:
      
        std::istream* openInStream(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in);
        void readFASTAfile(std::istream &, FileType &);
        void readFASTQfile(std::istream &, FileType &);
        void updateFileType(FileType &ft);
        bool convertToBinary(const std::string &str, FileType &ft);
        void organizeSequences(FileType &ft);

        uint64_t numBP;
        uint64_t numOfReads;
        uint64_t readArraySize;
        Sequence *binReads;
        Options *options;
};

#endif
