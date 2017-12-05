/* ================================================================= *
 *  Alignment.cpp: Source file for alignment algorithms              *
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

#include <Alignment.h>
#include <Timer.h>
#include <tuple>
/* Methods for class AlignmentData */

AlignmentData::AlignmentData ()
{
    rStart=0;
    qStart=0;
}

AlignmentData::AlignmentData (const int32_t &rs, const int32_t &qs)
{
    rStart=rs;
    qStart=qs;
}

bool AlignmentData::operator () (const AlignmentData &Obj1, const AlignmentData &Obj2) const
{      
    if (Obj1.rStart < Obj2.rStart)
        return true;
    else if (Obj1.rStart == Obj2.rStart)
    {
        if (Obj1.qStart < Obj2.qStart)
            return true;
        else
            return false;
    }else
        return false;
}    

/* Methods for class Alignment */
Alignment::Alignment()
{
}

void Alignment::Initialize(SequenceReader &r, SequenceReader &q, Options &p)
{
    ref=&r;
    query=&q;
    options=&p;
}

uint64_t Alignment::computeOverlapLength(uint64_t &qSize, uint64_t &rSize, AlignmentData &leftBound, AlignmentData &rightBound)
{
    uint64_t ovlLength=0;
    int kmerSize = kmerSizeArray[options->getSmerSeedIndex()];
    uint32_t rDiff = rSize - rightBound.rStart-kmerSize;
    uint32_t qDiff = qSize - rightBound.qStart-kmerSize;
    if (leftBound.rStart >= leftBound.qStart)
    {
        /*  ------------|----|--- r
         *           ---|----|-------- q
         *            or
         *  ---------------|----|--- r
         *  ---------------|----|---------- q
         *            or
         *  ---------------|----|---------- r
         *  ---------------|----|---------- q
         */ 
        if (rDiff <= qDiff)
            ovlLength = (rightBound.qStart + rDiff); 
        else {
            /*  ------------|----|----------- r
             *           ---|----|-------- q
             *            or
             *  ---------------|----|--- q
             *  ---------------|----|---------- r
             */
            ovlLength = (rightBound.qStart + qDiff); 
        }
    }else if (leftBound.rStart < leftBound.qStart)
    {
        if (rDiff <= qDiff){
           /*  ------------|----|---------------- q
            *           ---|----|-------- r
            */ 
           ovlLength = (rightBound.rStart + rDiff); 
        }else {
           /*  ------------|----|--- q
            *           ---|----|-------- r
            */ 
           ovlLength = (rightBound.rStart + qDiff); 
        }
    }
    return ovlLength;
}


int32_t Alignment::estimateExpectedNumOverlaps(uint64_t ovlLength, double &probability, int factor)
{
        float result = probability*ovlLength - factor*std::sqrt(probability*(1 - probability)*ovlLength); 
        if (result <= 0)
            return 0;
        else
            return static_cast<int32_t>(result);
}


void Alignment::clusterInitialKmers(std::vector<AlignmentData> &v, int32_t &kmerSize, float &maxShift, int32_t &direction, int &leftIndex, int &rightIndex, int &bpCount, int &kmerHits, std::vector<std::vector<int>> &positionVector, int &positionVecIndex)
{
    int32_t rDiff=0, qDiff=0, size=v.size();
    bool found=false;
    std::vector<std::tuple<int, int, int, int> > matchPositions; /* leftIndex, rightIndex, bpCount, kmerHitCount */
    positionVector.resize(1);
    positionVector[0].resize(1);

    if (direction == FORWARD)
    {
        matchPositions.emplace_back(0,0, kmerSize, 1);
        positionVector[0][0]=0;

        for (int32_t k=1; k < size; k++)
        {
            found=false;
            for (uint32_t j=0; j < matchPositions.size(); j++)
            {
                rDiff = v[k].rStart - v[std::get<1>(matchPositions[j])].rStart;
                if (rDiff < 0)
                    continue;
              
                qDiff = v[k].qStart - v[std::get<1>(matchPositions[j])].qStart;
                if (qDiff < 0)
                    continue;
            
                if (rDiff==qDiff) // Same distance apart
                {
                    std::get<2>(matchPositions[j])+=(qDiff>kmerSize)?kmerSize:qDiff;
                    std::get<3>(matchPositions[j])+=1;
                    std::get<1>(matchPositions[j])=k;
                    positionVector[j].emplace_back(k);
                    found=true;
                }else if (rDiff > qDiff)
                {
                    if (qDiff >= (rDiff*(1-maxShift))) 
                    {
                        std::get<2>(matchPositions[j])+=(qDiff>kmerSize)?kmerSize:qDiff;
                        std::get<3>(matchPositions[j])+=1;
                        std::get<1>(matchPositions[j])=k;
                        positionVector[j].emplace_back(k);
                        found=true;
                    } 
                }else 
                {
                    if (rDiff >= (qDiff*(1-maxShift))) 
                    {
                        std::get<2>(matchPositions[j])+=(qDiff>kmerSize)?kmerSize:qDiff;
                        std::get<3>(matchPositions[j])+=1;
                        std::get<1>(matchPositions[j])=k;
                        positionVector[j].emplace_back(k);
                        found=true;
                    } 
                }
            }     
            if (!found) {
                matchPositions.emplace_back(k,k, kmerSize,1);
                positionVector.resize(positionVector.size()+1);
                positionVector[positionVector.size()-1].emplace_back(k);
            }
        }
    }else { //REVERSE
        matchPositions.emplace_back(size-1, size-1, kmerSize, 1);
        positionVector[0][0]=size-1;

        for (int32_t k=size-2; k >=0; k--)
        {
            found=false;
            for (uint32_t j=0; j < matchPositions.size(); j++)
            {
                rDiff = v[std::get<0>(matchPositions[j])].rStart - v[k].rStart;
                if (rDiff < 0)
                    continue;
              
                qDiff = v[std::get<0>(matchPositions[j])].qStart - v[k].qStart;
                if (qDiff < 0)
                    continue;
            
                if (rDiff==qDiff) // Same distance apart
                {
                    std::get<2>(matchPositions[j])+=(qDiff>kmerSize)?kmerSize:qDiff;
                    std::get<3>(matchPositions[j])+=1;
                    std::get<0>(matchPositions[j])=k;
                    positionVector[j].emplace_back(k);
                    found=true;
                }else if (rDiff > qDiff)
                {
                    if (qDiff >= (rDiff*(1-maxShift))) 
                    {
                        std::get<2>(matchPositions[j])+=(qDiff>kmerSize)?kmerSize:qDiff;
                        std::get<3>(matchPositions[j])+=1;
                        std::get<0>(matchPositions[j])=k;
                        positionVector[j].emplace_back(k);
                        found=true;
                    } 
                }else 
                {
                    if (rDiff >= (qDiff*(1-maxShift))) 
                    {
                        std::get<2>(matchPositions[j])+=(qDiff>kmerSize)?kmerSize:qDiff;
                        std::get<3>(matchPositions[j])+=1;
                        std::get<0>(matchPositions[j])=k;
                        positionVector[j].emplace_back(k);
                        found=true;
                    } 
                }
            }     
            if (!found) {
                matchPositions.emplace_back(k,k, kmerSize,1);
                positionVector.resize(positionVector.size()+1);
                positionVector[positionVector.size()-1].emplace_back(k);
            }
        }
    }

    positionVecIndex=0;
    leftIndex =  std::get<0>(matchPositions[0]);
    rightIndex =  std::get<1>(matchPositions[0]);
    bpCount =  std::get<2>(matchPositions[0]);
    kmerHits =  std::get<3>(matchPositions[0]);
    for (uint32_t j=1; j < matchPositions.size(); j++){
        if (bpCount < std::get<2>(matchPositions[j]))
        {
            leftIndex = std::get<0>(matchPositions[j]);
            rightIndex = std::get<1>(matchPositions[j]);
            bpCount = std::get<2>(matchPositions[j]);
            kmerHits =  std::get<3>(matchPositions[j]);
            positionVecIndex = j;
        }
    }

    if (direction==REVERSE)
    {
        std::sort(positionVector[positionVecIndex].begin(), positionVector[positionVecIndex].end());
    }
}

void Alignment::findAllAlignments()
{
    float maxShift=options->getMaxShift();
    int32_t kmerSize = kmerSizeArray[options->getSmerSeedIndex()], smallKmerSize = kmerSizeArray[options->getSmallKmerIndex()];
    UINT_TYPE kmerSeed= kmerSeedArray[options->getSmerSeedIndex()], smallKmerSeed =  kmerSeedArray[options->getSmallKmerIndex()];
    int32_t shiftSize = CHARACTERS_IN_UINT-kmerSize;
    ReadDirection direction = FORWARD;
    uint64_t num_Ref_Seq =  ref->getNumberOfReads(), minStoreLength = static_cast<size_t> (options->getMinStoreLength());
    double probability = std::pow((1-options->getErrorRate()), 2*kmerSize); // (1-e)^2k     
    double probability_12 = std::pow((1-options->getErrorRate()), 2*smallKmerSize); // (1-e)^2k     
    int32_t minBPMatches = options->getMinMatches()*kmerSize;

    SmerMap mapRef(ref, options);
    mapRef.storeAllSmers();
    printMemUsage();
    if (index_write_flag) {
        exit(EXIT_SUCCESS);
    }
    Timer t("Done - Finding Alignments: ");
    t.start();
    #pragma omp parallel  num_threads(options->getNumThreads())
    {
      std::unordered_map<uint64_t, std::vector<AlignmentData>[2] > alignData;
      std::unordered_map<uint64_t, std::vector<uint32_t> > qHash;
      uint32_t count=0, m=0, positionR=0;
      UINT_TYPE *seq=NULL, readLoopMax=0, index=0, offset=0, key=0, currKmer=0,qSize=0, read=0;
      int32_t qLeft=0, qRight=0, rLeft=0, rRight=0; 
      std::vector<ReadHashNode> *referenceReadsMap;
      std::vector<std::vector<int>> positionVector;

      #pragma omp for schedule (dynamic, 1) 
      for (uint64_t i=0; i < query->getNumberOfReads(); i++)
      {
        qSize = query->getSequenceSize(i);
        if (qSize < minStoreLength)
            continue;
        if (self_flag && qSize < minStoreLength && i < num_Ref_Seq)
            continue;
             
        readLoopMax = qSize - kmerSize;
        alignData.clear();
        seq = query->getSequence(i, direction);
        count=0;
        for (uint32_t currKmerPos=0; currKmerPos<=readLoopMax; currKmerPos++)
        {
            if (count)
            {
                currKmer <<= 2;
                count--;
            }else {
                count=shiftSize;
                offset = size2Offset(currKmerPos+1);
                index=char2byteSize(currKmerPos+1);// current location in binReads 
                currKmer = seq[index];
            
                if (offset) {
                    currKmer <<= 2*offset;
                    currKmer |= (seq[index+1] >> 2*(CHARACTERS_IN_UINT - offset));
                }
            }

            if (mapRef.findSmer(currKmer,  kmerSeed, key)) 
            {
                referenceReadsMap=mapRef.getMapOfReads(key);
                if (referenceReadsMap->size()) 
                {
                    // We have a match - iterate over all entries in vector
                    for (std::vector<ReadHashNode>::iterator itMap = referenceReadsMap->begin(); itMap != referenceReadsMap->end(); ++itMap)
                    {
                        read = itMap->readIndex;
                        uint64_t refSize = ref->getSequenceSize(read);

                        if (refSize < minStoreLength && qSize < minStoreLength) 
                            continue;

                        if (self_flag && read <= i && i < num_Ref_Seq && refSize >= minStoreLength && qSize >= minStoreLength) //Self alignment is not needed
                            continue;

                        m=itMap->dir;
                        positionR=itMap->pos;
                        std::vector <AlignmentData> *v = alignData[read];
                        v[m].emplace_back(positionR, currKmerPos);
                    }
                }
            }
        }

        if (!alignData.size())
            continue;    // No matches to be processed here.

        computeQuerySmers(qHash, seq, qSize, smallKmerSeed, smallKmerSize);

        for (auto it = alignData.begin(); it != alignData.end(); ++it)
        {
            int32_t m=FORWARD, left=0, right=0, bpCount=0, kmerHits=0, est=0, positionVecIndex;
            uint64_t refSize=0, ovlLength=0, size=0;
            std::vector <AlignmentData> *v = NULL;

            read = it->first;
            refSize=ref->getSequenceSize(read), ovlLength=0;
            v = alignData[read];
            if (v[REVERSE].size() > v[m].size())
                m=REVERSE;

            if (v[m].size()<3)
                continue;

            clusterInitialKmers(v[m], kmerSize, maxShift, m, left, right, bpCount, kmerHits, positionVector, positionVecIndex);

            if (bpCount < minBPMatches)
            {
                std::sort(v[m].begin(), v[m].end(), AlignmentData());
                clusterInitialKmers(v[m], kmerSize, maxShift, m, left, right, bpCount, kmerHits, positionVector, positionVecIndex);
                if (bpCount < minBPMatches)
                    continue;
            }

            ovlLength = computeOverlapLength(qSize, refSize, v[m][left], v[m][right]);

            est = estimateExpectedNumOverlaps(ovlLength, probability, 2);
            if (!est) {
                continue;
            }else if (kmerHits < est){ 
                continue;
            }

            int score1 = computeSmallerKmerHits(read, v[m], positionVector[positionVecIndex],  qHash, smallKmerSeed, smallKmerSize, m, kmerSize, maxShift);

            qLeft =  v[m][left].qStart;
            qRight = v[m][right].qStart;
            rLeft =  v[m][left].rStart;
            rRight = v[m][right].rStart;
       

            int score = addSmallerAlignments(i, read, qLeft, qRight, rLeft, rRight, qHash, smallKmerSeed, smallKmerSize, m, probability_12, maxShift);

            if (rRight ==  v[m][right].rStart)
                size = kmerSize-1;
            else
                size = smallKmerSize-1;

            est = estimateExpectedNumOverlaps((qRight-qLeft > rRight-rLeft) ? rRight+size-rLeft : qRight+size-qLeft, probability_12, 3);
            if (est < options->getMinMatches()) {
                continue;
            }else if ((score + score1) < est){ 
                continue;
            }

            #pragma omp critical
            {
                if (self_flag) {
                    if (m == REVERSE){
                        std::cout << std::setw(5) << i+1 << " " << std::setw(5) << read+1 << " " << std::setw(5) << "0" << " " << "0" 
                            << " " << "0" << " " << qLeft+1 <<" " << (qRight+size) << " " << qSize << " " << m << " " 
                            << refSize - (rRight +  size) + 1 << " " <<  refSize - rLeft << " " << refSize <<  "\n";
                    }else {
                        std::cout << std::setw(5) << i+1 << " " << std::setw(5) << read+1 << " " << std::setw(5) << "0" << " " << "0" 
                            << " " << "0" << " " << qLeft+1 <<" " << (qRight+size) << " " << qSize << " " << m << " " 
                            << rLeft+1 << " " <<  (rRight+size) << " " << refSize <<  "\n";
                    }
                }else {
                    if (m == REVERSE){
                        std::cout << std::setw(5) << num_Ref_Seq+i+1 << " " << std::setw(5) << read+1 << " " << std::setw(5) << "0" << " " << "0" 
                            << " " << "0" << " " << qLeft+1 <<" " << (qRight+size) << " " << qSize << " " << m << " " 
                            << refSize - (rRight +  size) + 1 << " " <<  refSize - rLeft << " " << refSize <<  "\n";
                    }else {
                        std::cout << std::setw(5) << num_Ref_Seq+i+1 << " " << std::setw(5) << read+1 << " " << std::setw(5) << "0" << " " << "0" 
                            << " " << "0" << " " << qLeft+1 <<" " << (qRight+size) << " " << qSize << " " << m << " " 
                            << rLeft+1 << " " <<  (rRight+size) << " " << refSize <<  "\n";
                    }
                }
            }
        }
        qHash.clear();
      }
    }
    t.stop();
}

uint32_t Alignment::computeSmallerKmerHits(const uint64_t &refId, std::vector<AlignmentData> &v, std::vector<int> &positions, std::unordered_map<uint64_t, std::vector<uint32_t> > &qHash, uint64_t &smallKmerSeed, int32_t &smallKmerSize, int &refDirection, int32_t &kmerSize, float &maxShift)
{
    ReadDirection direction = static_cast<ReadDirection> (refDirection);
    UINT_TYPE *refSeq = ref->getSequence(refId, direction);
    UINT_TYPE refSeqSize = ref->getSequenceSize(refId);
    uint64_t index=0, offset=0, currSmer=0, readLoopMax = refSeqSize - smallKmerSize;
    int64_t rDiff, qDiff;
    std::vector<AlignmentData> tmpVecLeft, tmpVecRight;
    std::unordered_map<uint64_t, std::vector<uint32_t> >::const_iterator it;
    AlignmentData lastMatch;
    int delta = kmerSize-smallKmerSize+1, currSmerPos=0, score=0;
 
    score += 1;  //For the last 16mer converted to 12 mers
    for(uint32_t i=0; i<positions.size()-1; i++)
    {
        score += 1;  
        currSmerPos=v[positions[i]].rStart+delta;
        readLoopMax=v[positions[i+1]].rStart-1;
        lastMatch.rStart=v[positions[i]].rStart;
        lastMatch.qStart=v[positions[i]].qStart;
        while(currSmerPos <= static_cast<int64_t>(readLoopMax))
        {
            offset = size2Offset(currSmerPos+1);
            index=char2byteSize(currSmerPos+1);
            currSmer = refSeq[index];
            if (offset)
            {
                if (CHARACTERS_IN_UINT-offset < static_cast<uint32_t>(smallKmerSize))
                    currSmer = ((currSmer << 2*offset)| (refSeq[index+1] >> 2*(CHARACTERS_IN_UINT - offset)));
                else
                    currSmer <<= 2*offset;
            }

            it=qHash.find(currSmer & smallKmerSeed);
            if (it!=qHash.end())
            {
                std::vector<uint32_t>::const_iterator low,up;
                low=std::lower_bound (it->second.begin(), it->second.end(), lastMatch.qStart+1);

                if (low==it->second.end()){
                    currSmerPos++;
                    continue;
                }

                if (*low > static_cast<uint32_t>(v[positions[i+1]].qStart-1)){
                    currSmerPos++;
                    continue;
                }
                up=std::upper_bound (low, it->second.end(),  v[positions[i+1]].qStart-1); 
                if (up==it->second.end() && *(up-1) >  static_cast<uint32_t>(v[positions[i+1]].qStart-1)){
                    currSmerPos++;
                    continue;
                }
                while (low != it->second.end() && *low <= *(up-1))
                {
                    rDiff = currSmerPos - lastMatch.rStart;
                    qDiff = *low - lastMatch.qStart;
                    if (rDiff == qDiff) {
                        score+=1; //Kmer-count
                        lastMatch.rStart=currSmerPos;
                        lastMatch.qStart=*low;
                        break; // one match per positions is good enough
                    }else if (rDiff > qDiff) {
                        if (qDiff >= (rDiff*(1-maxShift))) { 
                            score+=1;
                            lastMatch.rStart=currSmerPos;
                            lastMatch.qStart=*low;
                            break; // one match per positions is good enough
                        }
                    }else {
                        if (rDiff >= ((qDiff*(1-maxShift)))) { 
                            lastMatch.rStart=currSmerPos;
                            lastMatch.qStart=*low;
                            score+=1;
                            break; // one match per positions is good enough
                        }
                    }
                    low++;
                }
            }
            currSmerPos++;
        }
    }

    return score;
}

uint32_t Alignment::extendOvlUsingSmallerKmer(int32_t &queryBound, int32_t &refBound, std::vector<AlignmentData> &tmpVec, float &maxShift, const uint32_t &kmerSize, double &probability, bool isLeftExt)
{
    int32_t useIndex=-1, rDiff=0, qDiff=0;
    int32_t est=0, score=1;

    for (uint32_t i=0; i < tmpVec.size(); i++)
    {
        if (isLeftExt)
        {
            if (useIndex==-1)
            {
                rDiff = refBound-tmpVec[i].rStart;
                qDiff = queryBound-tmpVec[i].qStart;
            }else {
                if (tmpVec[i].rStart >=  tmpVec[useIndex].rStart)
                    continue;
                if (tmpVec[i].qStart >=  tmpVec[useIndex].qStart)
                    continue;
                rDiff = tmpVec[useIndex].rStart-tmpVec[i].rStart;
                qDiff = tmpVec[useIndex].qStart-tmpVec[i].qStart;
            }
        }else { //Right extension
            if (useIndex==-1)
            {
                rDiff = tmpVec[i].rStart - refBound;
                qDiff = tmpVec[i].qStart - queryBound;
            }else {
                if (tmpVec[i].rStart <=  tmpVec[useIndex].rStart)
                    continue;
                if (tmpVec[i].qStart <=  tmpVec[useIndex].qStart)
                    continue;
                rDiff = tmpVec[i].rStart-tmpVec[useIndex].rStart;
                qDiff = tmpVec[i].qStart-tmpVec[useIndex].qStart;
            }
        }

        if (rDiff == qDiff) {
            est = estimateExpectedNumOverlaps((isLeftExt ? refBound-tmpVec[i].rStart : tmpVec[i].rStart-refBound+kmerSize), probability, 3);
            if (score >= est) {
                useIndex=i;
                score++;
            }else if (useIndex != -1)
                break;
        }else if (rDiff > qDiff) {
            if (qDiff >= (rDiff*(1-maxShift))) { 
                est = estimateExpectedNumOverlaps(isLeftExt ? refBound-tmpVec[i].rStart : tmpVec[i].rStart-refBound+kmerSize, probability, 3);
                if (score >= est) {
                    useIndex=i;
                    score++;
                }else if (useIndex != -1)
                    break;
            }
        }else {
            if (rDiff >= (qDiff*(1-maxShift))) { 
                est = estimateExpectedNumOverlaps(isLeftExt ? refBound-tmpVec[i].rStart : tmpVec[i].rStart-refBound+kmerSize, probability, 3);
                if (score >= est) {
                    useIndex=i;
                    score++;
                }else if (useIndex != -1)
                    break;
            }
        }
    }
    
    if (useIndex != -1) 
    {
        if (isLeftExt)
        {
            queryBound=tmpVec[useIndex].qStart;
            refBound=tmpVec[useIndex].rStart;
        }else {
            refBound=tmpVec[useIndex].rStart;
            queryBound=tmpVec[useIndex].qStart;
        }
        return score;
    }else
        return 0;
}

uint32_t Alignment::addSmallerAlignments(const uint64_t &queryId, const uint64_t &refId, int32_t &qLeft, int32_t &qRight, int32_t &rLeft, int32_t &rRight, std::unordered_map<uint64_t, std::vector<uint32_t> > &qHash, const uint64_t &kmerSeed, const uint32_t &kmerSize, int &refDirection, double &probability, float &maxShift)
{
    ReadDirection direction = static_cast<ReadDirection> (refDirection);
    UINT_TYPE *refSeq = ref->getSequence(refId, direction);
    int64_t refSeqSize = ref->getSequenceSize(refId), querySeqSize = query->getSequenceSize(queryId);
    int64_t readLoopMax = refSeqSize - kmerSize;
    int32_t currSmerPos=0, score=0, lengthRight=0, leftRef=0, leftQuery=0, rightRef=0, rightQuery=0, i=0;
    uint64_t index=0, offset=0, currSmer=0;
    bool last=false;
    std::vector<AlignmentData> tmpVecLeft, tmpVecRight;
    std::unordered_map<uint64_t, std::vector<uint32_t> >::const_iterator it;

    if (rLeft > qLeft)
    {
        /*  ------------|----|---
         *           ---|----|--------
         */ 
        leftQuery=0;
        leftRef = rLeft - qLeft; 
        if (leftRef > (maxShift*qLeft))
            leftRef -= (maxShift*qLeft);
        else  
            leftRef=0;
    }else if (rLeft < qLeft)
    {
        leftRef = 0;
        leftQuery =  qLeft - rLeft;
        if (leftQuery > (maxShift*rLeft))
            leftQuery -= (maxShift*rLeft);
        else  
            leftQuery=0;
    }else //Equal
    {
        leftRef=0;
        leftQuery=0;
    }

    if (refSeqSize - rRight > querySeqSize - qRight)
    {
        lengthRight = querySeqSize - qRight;
        rightQuery=querySeqSize-kmerSize;
        rightRef=rRight+ ((maxShift + 1.0) * lengthRight);
        if (rightRef > readLoopMax)
            rightRef=readLoopMax;
    }else if (refSeqSize - rRight < querySeqSize - qRight)
    {
        lengthRight = refSeqSize - rRight;
        rightRef=readLoopMax;
        rightQuery=qRight + ((maxShift + 1.0) * lengthRight);
        if (rightQuery > querySeqSize)
            rightQuery=querySeqSize-kmerSize;
    }else {
        rightRef=readLoopMax;
        rightQuery=querySeqSize-kmerSize;
    }

    i=0;
    currSmerPos=rLeft-1;
    while(true)
    {
        while(currSmerPos <= static_cast<int64_t>(readLoopMax))
        {
            if (!last)
            {
                if (!i && currSmerPos < static_cast<int32_t>(leftRef))
                    break;
                else if (i == 1)
                {
                    /*                    ---|----|--------------
                     *          -------------|----|----
                     */ 
                    last=true;
                }
            }else if (currSmerPos > rightRef)
                break;

            offset = size2Offset(currSmerPos+1);
            index=char2byteSize(currSmerPos+1);
            currSmer = refSeq[index];
            if (offset)
            {
                if (CHARACTERS_IN_UINT-offset < kmerSize)
                    currSmer = ((currSmer << 2*offset)| (refSeq[index+1] >> 2*(CHARACTERS_IN_UINT - offset)));
                else
                    currSmer <<= 2*offset;
            }

            it=qHash.find(currSmer & kmerSeed);
            if (it!=qHash.end())
            {
                std::vector<uint32_t>::const_iterator low,up;
                if (last)
                    low=std::lower_bound (it->second.begin(), it->second.end(), static_cast<uint32_t>(qRight)+1);
                else 
                    low=std::lower_bound (it->second.begin(), it->second.end(), leftQuery); 

                if (low==it->second.end()){
                    if (last)
                        currSmerPos++;
                    else
                        currSmerPos--;
                    continue;
                }

                if (last){
                    if (*low > static_cast<uint32_t>(rightQuery)){
                        currSmerPos++;
                        continue;
                    }
                    up=std::upper_bound (low, it->second.end(), rightQuery); 
                    if (up==it->second.end() && *(up-1) > static_cast<uint32_t>(rightQuery)){
                        currSmerPos++;
                        continue;
                    }
                    while (low != it->second.end() && *low <= *(up-1)){
                        tmpVecRight.emplace_back(currSmerPos, *low);
                        low++;
                    }
                }else{
                    while (low != it->second.end() && *low < static_cast<uint32_t>(qLeft)){
                        tmpVecLeft.emplace_back(currSmerPos, *low);
                        low++;
                    }
                }
            }
            if (last)
                currSmerPos++;
            else
                currSmerPos--;
        }
        if (!i){
            i=1;
            currSmerPos=rRight+1;
            if (currSmerPos > rightRef)
                break;
        }
        if (last)
            break;
    }

    score += extendOvlUsingSmallerKmer(qLeft, rLeft, tmpVecLeft, maxShift, kmerSize, probability, true);
    score += extendOvlUsingSmallerKmer(qRight, rRight, tmpVecRight, maxShift,  kmerSize, probability, false);

    return score;
}


void Alignment::computeQuerySmers(std::unordered_map<uint64_t, std::vector<uint32_t> > &qHash, UINT_TYPE * &seq, uint64_t &seqSize,  uint64_t &smallKmerSeed, int32_t &smallKmerSize)
{
    UINT_TYPE currSmer;
    uint64_t index=0, offset=0; 
    uint64_t  readLoopMax = seqSize - smallKmerSize;
    int shiftSize=CHARACTERS_IN_UINT-smallKmerSize;
    int count=0;

    qHash.reserve(readLoopMax+1);

    for (uint64_t currSmerPos=0; currSmerPos <= readLoopMax; currSmerPos++)
    {
        if (count)
        {
            currSmer <<= 2;
            count--;
        }else
        {
            count=shiftSize;
            offset = size2Offset(currSmerPos+1);
            index=char2byteSize(currSmerPos+1);// current location in binReads 
            currSmer = seq[index];
    
            if (offset) {
                currSmer <<= 2*offset;
                currSmer |= (seq[index+1] >> 2*(CHARACTERS_IN_UINT - offset));
            }
        }
       
        qHash[currSmer & smallKmerSeed].emplace_back(currSmerPos);
   }
}
