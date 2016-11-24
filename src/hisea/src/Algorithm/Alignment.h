/* ================================================================= *
 *  Alignment.h: Header file for alignment algorithms                *
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


#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <Util.h>
#include <Globals.h>
#include <SequenceReader.h>
#include <KmerMap.h>
#include <Options.h>

class AlignmentData
{
    public:
        int32_t rStart;
        int32_t qStart;
        AlignmentData ();
        AlignmentData (const int32_t &rs, const int32_t &qs);
        bool operator () (const AlignmentData &Obj1, const AlignmentData &Obj2) const;
};

class Alignment 
{
    public:
        Alignment();
        void findAllAlignments();
        void Initialize(SequenceReader &r, SequenceReader &q, Options &p);

    private:
        void clusterInitialKmers(std::vector<AlignmentData> &v, int32_t &kmerSize, float &maxShift, int32_t &direction, int &leftIndex, int &rightIndex, int &bpCount, int &kmerHits, std::vector<std::vector<int>> &positionVector, int &positionVecIndex);
        uint32_t extendOvlUsingSmallerKmer(int32_t &queryBound, int32_t &refBound, std::vector<AlignmentData> &tmpVec, float &maxShift, const uint32_t &kmerSize, double &probability, bool isLeftExt);
        int32_t estimateExpectedNumOverlaps(uint64_t ovlLength, double &probability, int factor);
        uint64_t computeOverlapLength(uint64_t &qSize, uint64_t &rSize, AlignmentData &leftBound, AlignmentData &rightBound);
        void computeQuerySmers(std::unordered_map<uint64_t, std::vector<uint32_t> > &qHash, UINT_TYPE * &seq, uint64_t &seqSize,  uint64_t &smerSeed, int32_t &smerSize);
        uint32_t addSmallerAlignments(const uint64_t &queryId, const uint64_t &refId, int32_t &qLeft, int32_t &qRight, int32_t &rLeft, int32_t &rRight, std::unordered_map<uint64_t, std::vector<uint32_t> > &qHash, const uint64_t &kmerSeed, const uint32_t &kmerSize, int &refDirection, double &probability, float &maxShift);
         uint32_t computeSmallerKmerHits(const uint64_t &refId, std::vector<AlignmentData> &v, std::vector<int> &positions, std::unordered_map<uint64_t, std::vector<uint32_t> > &qHash, uint64_t &smerSeed, int32_t &smerSize, int &refDirection, int32_t &kmerSize, float &maxShift);


        SequenceReader *ref;
        SequenceReader *query;
        Options *options;
};

#endif
