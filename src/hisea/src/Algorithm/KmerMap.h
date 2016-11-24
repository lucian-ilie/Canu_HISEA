/* ================================================================= *
 *  KmerMap.h: Source file for k-mer/hash processing algorithms      *
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


#ifndef SMERMAP_H
#define SMERMAP_H

#include <cmath>
#include <array>
#include <Util.h>
#include <Globals.h>
#include <SequenceReader.h>
#include <Options.h>
#include <iostream>

class ReadHashNode {
    public:
        uint64_t readIndex : 43;
        uint64_t pos       : 20;
        uint64_t dir       :  1; // Direction 0=Fwd and 1=Rev
        ReadHashNode(uint64_t &r, uint64_t &p, uint64_t &d);
        ReadHashNode();
};

class SmerNode 
{
    public:
        SmerNode();
        SmerNode(uint64_t &r);
        bool operator == (const SmerNode &node) const;
        uint64_t smerVal;
        std::vector <ReadHashNode> *readMap; 
        
};

class SmerMap 
{
    public:
        bool findSmer(uint64_t &, uint64_t &, uint64_t &);
        void storeAllSmers();
        std::vector<ReadHashNode> *getMapOfReads(uint64_t &key);
        explicit SmerMap (SequenceReader *readerObj, Options *op);
        ~SmerMap ();

    private:
        uint64_t getHashKey(uint64_t &smerVal, int &numLocks, uint64_t &kmserSeed, bool &match_found);
        uint64_t getHashKeyKmerCount(uint64_t &smerVal, int &numLocks, uint64_t &kmserSeed, bool &match_found);
        void destroyLocks();
        void initializeLocks();
        inline int getNumLocks();
        void rehashSmers();
        void readIndexFromFile(std::unordered_map<uint64_t, char> &ignoreKmerHash);
        void writeIndexToFile(std::unordered_map<uint64_t, char> &ignoreKmerHash);
        uint64_t collectKmersParallel(std::unordered_map<uint64_t, char>& ignoreKmerHash);
        void initializeHashForReference(std::unordered_map<uint64_t, char>&ignoreKmerHash, uint64_t &hashSize);
        void populateIgnoreKmerHash(std::unordered_map<uint64_t, char> &ignoreKmerHash);
        void insertReferencePositions(std::unordered_map<uint64_t, char> &ignoreKmerHash);

        uint64_t currHashTabSize;
        uint64_t prevHashTabSize;
        SmerNode *SmerArray; 
        uint64_t *KmerCountArray;
        SequenceReader *reader;
        Options *options;
        omp_lock_t *myLocks;
        UINT_TYPE numElements;
};
#endif
