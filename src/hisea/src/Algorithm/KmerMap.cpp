/* ================================================================= *
 *  KmerMap.cpp: Source file for k-mer/hash processing algorithms    *
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


#include <KmerMap.h>
#include <Timer.h>
#include <Globals.h>
#include <iostream>
#include <bitset>
#include <iterator>
#include <fcntl.h>

/* Methods for ReadMap class */
ReadHashNode::ReadHashNode(uint64_t &r, uint64_t &p, uint64_t &d)
{
    readIndex=r;
    pos=p;
    dir=d;
}

ReadHashNode::ReadHashNode()
{
    readIndex=0;
    pos=0;
    dir=0;
}

/* Methods for SmerNode class */
SmerNode::SmerNode(uint64_t& r)
{
    this->smerVal=r;
    this->readMap=NULL;
}

SmerNode::SmerNode()
{
    this->smerVal=HASH_INITIAL_VALUE;
    this->readMap=NULL;
}

bool SmerNode::operator == (const SmerNode &node) const
{       
    if (this->smerVal != node.smerVal)
        return false;
    else
        return true;
} 

/* Methods for SmerMap class */
bool SmerMap::findSmer(uint64_t &smerVal, uint64_t& kmerSeed, uint64_t &key)
{
    bool match=false;
    int numLocks = 0;
    uint64_t smerToFind = smerVal & kmerSeed;
    key=this->getHashKey(smerToFind, numLocks, kmerSeed, match);
    if (match) {
       return true;
    }else
       return false;
}

uint64_t SmerMap::getHashKeyKmerCount(uint64_t &smerVal, int &numLocks, uint64_t& kmerSeed, bool &match_found)
{
    uint64_t key=smerVal%currHashTabSize;
    uint64_t step = 1 + (smerVal%(currHashTabSize-1));

    match_found=false;
    if(numLocks)
        omp_set_lock(myLocks + key%numLocks);

    while (KmerCountArray[key] != HASH_INITIAL_VALUE) {
        if ((KmerCountArray[key] & kmerSeed)==smerVal){
            match_found=true;
            return key;
        }
        else {
            if(numLocks)
                omp_unset_lock(myLocks + key%numLocks);

            key = (key + step) % currHashTabSize;

            if(numLocks)
                omp_set_lock(myLocks + key%numLocks);
        }
    }
    return key; //no collision
}

uint64_t SmerMap::getHashKey(uint64_t &smerVal, int &numLocks, uint64_t& kmerSeed, bool &match_found)
{
    uint64_t key=smerVal%currHashTabSize;
    uint64_t step = 1 + (smerVal%(currHashTabSize-1));

    match_found=false;
    if(numLocks)
        omp_set_lock(myLocks + key%numLocks);

    while (SmerArray[key].smerVal != HASH_INITIAL_VALUE) {
        if ((SmerArray[key].smerVal & kmerSeed)==smerVal){
            match_found=true;
            return key;
        }
        else {
            if(numLocks)
                omp_unset_lock(myLocks + key%numLocks);

            key = (key + step) % currHashTabSize;

            if(numLocks)
                omp_set_lock(myLocks + key%numLocks);
        }
    }
    return key; //no collision
}


uint64_t SmerMap::collectKmersParallel(std::unordered_map<uint64_t, char> &ignoreKmerHash)
{
    UINT_TYPE kmerSeed= kmerSeedArray[options->getSmerSeedIndex()];
    int kmerSize = kmerSizeArray[options->getSmerSeedIndex()];
    int filterCount=options->getFilterCount();
    UINT_TYPE goodKmer=0, hashFull=0;
    int shiftSize=CHARACTERS_IN_UINT-kmerSize, numLocks=0, rehash=1;
    ReadDirection direction;
    bool ignoreFile = options->getFileName(IGNORE_FILE).size() ? true : false;
    uint64_t maxKmerHits = options->getMaxKmerHits();

    while(rehash)
    {
        rehash=0;
        numElements=0;
        hashFull = (currHashTabSize*50)/100;
        for (int j=FORWARD; j<=REVERSE; j++)
        {
          direction = static_cast<ReadDirection>(j);

          #pragma omp parallel  num_threads(options->getNumThreads())
          {
            UINT_TYPE seqIndex=0, currKmerPos=0, currKmer=0, seqOffset=0, key=0;
            UINT_TYPE *sequenceObj = NULL, readSizeLoopMax=0, kmerToFind=0, count=0, refSize=0;
            bool match=false;
            #pragma omp for  
            for (uint64_t i=0; i<reader->getNumberOfReads(); i++)
            {
                if(!omp_get_thread_num() && !rehash) // Possibly master
                {
                    if (hashFull < numElements)
                        rehash++;
                }

                if (!rehash) 
                {
                    refSize=reader->getSequenceSize(i);
                    if (options->getMinOlapLength() > static_cast<int64_t>(refSize)) 
                        continue;
                    readSizeLoopMax=refSize - kmerSize;
                    sequenceObj = reader->getSequence(i, direction);
                    count=0;
                    for (currKmerPos=0; currKmerPos<=readSizeLoopMax; currKmerPos++)
                    {
                        if (count)
                        {
                            currKmer <<= 2;
                            count--;
                        }else
                        { 
                            count=shiftSize;
                            seqOffset = size2Offset(currKmerPos+1);
                            seqIndex=char2byteSize(currKmerPos+1);
                            currKmer = sequenceObj[seqIndex];
        
                            if (seqOffset) {
                                currKmer <<= 2*seqOffset;
                                currKmer |= (sequenceObj[seqIndex+1] >> 2*(CHARACTERS_IN_UINT - seqOffset));
                            }
                        }

                        kmerToFind=currKmer & kmerSeed;
                        key=this->getHashKeyKmerCount(kmerToFind, numLocks, kmerSeed, match);

                        if (self_flag && !query_flag && !index_write_flag && j==REVERSE && !match) {
                           if(numLocks)
                               omp_unset_lock(myLocks + key%numLocks);
                            continue;
                        }
                        if (!match) { 
                            kmerToFind++;
                            KmerCountArray[key]=kmerToFind; 
                            #pragma omp atomic
                            numElements++;
                        }else if ((KmerCountArray[key] & KMER_COUNT_MASK) < KMER_COUNT_MASK)
                            KmerCountArray[key]++; 
                        
                        if(numLocks)
                            omp_unset_lock(myLocks + key%numLocks);
                    }
                }
            }
          }
          if (rehash)
              break;
        }

        if (rehash)
            rehashSmers();
    }

    std::cerr << "Number of K-mers before filtering: " << numElements << std::endl; 

    #pragma omp parallel num_threads(options->getNumThreads())
    {
        #pragma omp for 
        for (uint64_t i=0; i < currHashTabSize; i++)
        {
            if (KmerCountArray[i] != HASH_INITIAL_VALUE) 
            {
                double count = KmerCountArray[i] & KMER_COUNT_MASK;
                if (ignoreFile) 
                { 
                    if (ignoreKmerHash.find(KmerCountArray[i] & kmerSeed) != ignoreKmerHash.end())
                        continue;
                }else {
                    if(count < filterCount || count > maxKmerHits) 
                    {
                        if (index_write_flag)
                        {
                            #pragma omp critical
                            ignoreKmerHash[KmerCountArray[i] & kmerSeed]='1'; 
                        }
                        continue;
                    }
                }

                #pragma omp atomic
                goodKmer++;
            }
        }
    }

    std::cerr << "Number of K-mers after filtering: " << goodKmer << std::endl; 
    return goodKmer;
}

void SmerMap::writeIndexToFile(std::unordered_map<uint64_t, char> &ignoreKmerHash)
{
    uint64_t num_reads = reader->getNumberOfReads();
    uint64_t kmerIgnore=0;
    std::string filename = options->getFileName(INDEX_FILE) + ".dat";
    FILE *output_file = fopen(filename.c_str(), "wb");
    fwrite(&self_flag, sizeof(char), sizeof(int), output_file);
    fwrite(&num_reads, sizeof(char), sizeof(uint64_t), output_file);
    fwrite(&currHashTabSize, sizeof(char), sizeof(uint64_t), output_file);

    if (!options->getFileName(IGNORE_FILE).size()) // No ignore file - kmers computed
    {
        for (auto it=ignoreKmerHash.begin(); it != ignoreKmerHash.end(); ++it)
        {
            kmerIgnore=it->first;
            fwrite(&kmerIgnore, sizeof(char), sizeof(uint64_t), output_file);
        } 
    } 

    fclose(output_file);
}

void SmerMap::readIndexFromFile(std::unordered_map<uint64_t, char> &ignoreKmerHash)
{
    uint64_t num_reads=0, kmerIgnore=0; 
    int self=0;
    std::string filename = options->getFileName(INDEX_FILE) + ".dat";
    FILE *read_file = fopen(filename.c_str(), "rb");
    fread(&self, sizeof(char), sizeof(int), read_file);
    if (self && query_flag && !self_flag) {
        std::cerr << "[ERROR] - Index must be created without self flag when using a query file and no self alignment with reference is needed" << std::endl;
        fclose(read_file);
        exit(EXIT_FAILURE); 
    }
    fread(&num_reads, sizeof(char), sizeof(uint64_t), read_file);
    if (num_reads != reader->getNumberOfReads())
    {
        std::cerr << "[ERROR] - It seems current referece file and file used for index creation are not same" << std::endl;
        std::cerr << "[ERROR] - Index was written using a file with " << num_reads << " reads. The current file has " << reader->getNumberOfReads() << " reads" << std::endl;
        std::cerr << "[ERROR] - Plese use same file for index creation and reading" << std::endl;
        fclose(read_file);
        exit(EXIT_FAILURE); 
    }
    fread(&currHashTabSize, sizeof(char), sizeof(uint64_t), read_file);
    std::cerr << "Hash size from file: " << currHashTabSize << std::endl;
    if (!options->getFileName(IGNORE_FILE).size()) // No ignore file - kmers computed
    {
        while(fread(&kmerIgnore, sizeof(char), sizeof(uint64_t), read_file) == sizeof(uint64_t))
            ignoreKmerHash[kmerIgnore]='1';
    } 
    fclose(read_file);
}

void SmerMap::insertReferencePositions(std::unordered_map<uint64_t, char> &ignoreKmerHash)
{
    UINT_TYPE kmerSeed= kmerSeedArray[options->getSmerSeedIndex()];
    int kmerSize = kmerSizeArray[options->getSmerSeedIndex()];
    int shiftSize=CHARACTERS_IN_UINT-kmerSize;
    ReadDirection direction;
    int numLocks = getNumLocks();
   
    if (index_read_flag) 
        SmerArray = new SmerNode[currHashTabSize];

    for (uint64_t j=FORWARD; j <=REVERSE; j++)
    {
        direction = static_cast<ReadDirection>(j);
        #pragma omp parallel  num_threads(options->getNumThreads())
        {
          UINT_TYPE seqIndex=0, currKmer=0, seqOffset=0, key=0, smerToFind=0;
          uint64_t currKmerPos=0, readSizeLoopMax=0;
          UINT_TYPE *sequenceObj = NULL;
          bool match=false;
          int count;
          #pragma omp for
          for (uint64_t i=0; i<reader->getNumberOfReads(); i++)
          {
            readSizeLoopMax=reader->getSequenceSize(i) - kmerSize;
            sequenceObj = reader->getSequence(i, direction);
            count=0;
            for (currKmerPos=0; currKmerPos<=readSizeLoopMax; currKmerPos++)
            {
                if (count)
                {
                    currKmer <<= 2;
                    count--;
                }else
                { 
                    count=shiftSize;
                    seqOffset = size2Offset(currKmerPos+1);
                    seqIndex=char2byteSize(currKmerPos+1);
                    currKmer = sequenceObj[seqIndex];
        
                    if (seqOffset) {
                        currKmer <<= 2*seqOffset;
                        currKmer |= (sequenceObj[seqIndex+1] >> 2*(CHARACTERS_IN_UINT - seqOffset));
                    }
                }
                smerToFind=currKmer & kmerSeed;

                if (index_read_flag && ignoreKmerHash.find(smerToFind) != ignoreKmerHash.end())
                    continue;
       
                key=this->getHashKey(smerToFind, numLocks, kmerSeed, match);
                if (match){
                    if (!index_read_flag) 
                    { 
                        if (SmerArray[key].readMap == NULL)
                            SmerArray[key].readMap=new std::vector<ReadHashNode>; 
                    }
                       
                    SmerArray[key].readMap->emplace_back(i, currKmerPos, j); 
                }else {
                    if (index_read_flag) { 
                        SmerArray[key].smerVal=smerToFind; 
                        SmerArray[key].readMap=new std::vector<ReadHashNode>; 
                        SmerArray[key].readMap->emplace_back(i, currKmerPos, j); 
                    }
                }
                if(numLocks)
                    omp_unset_lock(myLocks + key%numLocks);
            }
          }
        }
    }

    #pragma omp for
    for (uint64_t i=0; i < currHashTabSize; i++)
    {
        if (SmerArray[i].smerVal != HASH_INITIAL_VALUE)
        {
            SmerArray[i].readMap->shrink_to_fit();
        }
    }
}

void SmerMap::initializeHashForReference(std::unordered_map<uint64_t, char> &ignoreKmerHash, uint64_t &hashSize)
{
    int filterCount=options->getFilterCount(), n=0;
    UINT_TYPE kmerSeed= kmerSeedArray[options->getSmerSeedIndex()];
    int numLocks = getNumLocks();
    uint64_t oldHashSize=currHashTabSize;
    bool ignoreFile = options->getFileName(IGNORE_FILE).size() ? true : false;
    int64_t maxKmerHits = options->getMaxKmerHits();

    for (n=0; n<450; ++n)
    {
        if (hashTableSize[n] > hashSize)
           break;
    }

    currHashTabSize = hashTableSize[n];

    std::cerr << "readPos Current Hash Table Size: " << currHashTabSize << std::endl;
    if (n)
        prevHashTabSize = hashTableSize[n-1];
    else
        prevHashTabSize = 1;

    std::cerr << "Size of SMER node: " << sizeof(SmerNode) << std::endl;
    SmerArray=  new SmerNode[currHashTabSize];

    #pragma omp parallel num_threads(options->getNumThreads())
    {
      uint64_t key=0, kmerToFind=0;
      bool match=false;
      #pragma omp for
      for (uint64_t i=0; i < oldHashSize; i++)
      {
        if (KmerCountArray[i] != HASH_INITIAL_VALUE) {
            int count = KmerCountArray[i] & KMER_COUNT_MASK;
            if (ignoreFile) 
            { 
                if (ignoreKmerHash.find(KmerCountArray[i] & kmerSeed) != ignoreKmerHash.end())
                    continue;
            }else {
                if(count < filterCount || count > maxKmerHits)
                    continue;
            }
            kmerToFind=KmerCountArray[i] & kmerSeed;
            key=this->getHashKey(kmerToFind, numLocks, kmerSeed, match);
            SmerArray[key].smerVal=kmerToFind; 
            if(numLocks)
                omp_unset_lock(myLocks + key%numLocks);
        }
      }
    }
}

void SmerMap::populateIgnoreKmerHash(std::unordered_map<uint64_t, char> &ignoreKmerHash)
{
    std::string line;
    char *pch=NULL, *pch1=NULL;
    char str[1000];

    if (options->getFileName(IGNORE_FILE).size())
    {
        ignoreKmerHash.reserve(10);
        int kmerSize = kmerSizeArray[options->getSmerSeedIndex()];
        std::string filename = options->getFileName(IGNORE_FILE);
        std::istream *fs;
        if (!std::strcmp(filename.substr(filename.length()-std::strlen(".gz"), strlen(".gz")).c_str(), ".gz")) {
            igzstream* gfp = new igzstream(filename.c_str(), std::ios_base::in);
            if(!gfp->good()) {
                std::cerr << "[ERROR] - file open failed - ignore file" << filename << std::endl;
                exit(EXIT_FAILURE);
            }
            fs=gfp;
        }else {
            std::ifstream* fp = new std::ifstream(filename.c_str(), std::ios_base::in);
            if (!fp->is_open()) {
                std::cerr << "[ERROR] - file open failed " << filename << std::endl;
                exit(EXIT_FAILURE);
            }
            fs=fp;
        }

        getline(*fs, line); // line with numebr of kmers
        while(getline(*fs, line).good())
        {
            if(line.empty()) //blank line
                continue;

            strcpy(str, line.c_str());
            pch = strtok (str," \t");
            pch1 = strtok (NULL, " \t");
            if(pch == NULL || pch1 == NULL)
            {
                std::cerr << "[ERROR] - Ignore File format error: " << options->getFileName(IGNORE_FILE) << std::endl;
                exit(EXIT_FAILURE);
            }

            uint64_t kmerToIgnore=0, kmerToIgnoreRev=0;
            int count=0;
            while(count < kmerSize)
            {
                if (pch[count] == 'A' || pch[count] == 'a') {
                    kmerToIgnore <<= 2;
                    kmerToIgnoreRev |= 0xC000000000000000;
                    kmerToIgnoreRev >>= 2;
                }else if (pch[count] == 'C' || pch[count] == 'c') {
                    kmerToIgnore |= 0x1;
                    kmerToIgnore <<= 2;
                    kmerToIgnoreRev |= 0x8000000000000000;
                    kmerToIgnoreRev >>= 2;
                }else if (pch[count] == 'G' || pch[count] == 'g') {
                    kmerToIgnore |= 0x2;
                    kmerToIgnore <<= 2;
                    kmerToIgnoreRev |= 0x4000000000000000;
                    kmerToIgnoreRev >>= 2;
                }else if (pch[count] == 'T' || pch[count] == 't') {
                    kmerToIgnore |= 0x3;
                    kmerToIgnore <<= 2;
                    kmerToIgnoreRev >>= 2;
                }
                count++;
            }
            kmerToIgnore <<= 2*(CHARACTERS_IN_UINT-kmerSize-1);
            kmerToIgnoreRev <<= 2;

            ignoreKmerHash[kmerToIgnore]='1'; 
            ignoreKmerHash[kmerToIgnoreRev]='1'; 
        }
        delete fs;
    }else {
        std::cerr << "[WARNING] - No ignore file passed on the command line" << std::endl;
    }
}

void SmerMap::storeAllSmers()
{
    std::unordered_map<uint64_t, char> ignoreKmerHash;
    uint64_t hashSize;
    std::cerr << "Num of Reference Reads " << reader->getNumberOfReads() <<std::endl;
    std::cerr << "Num of Base Pairs " << reader->getNumberOfBP() <<std::endl;

    if (index_read_flag){
        if (options->getFileName(IGNORE_FILE).size())
            populateIgnoreKmerHash(ignoreKmerHash);
        readIndexFromFile(ignoreKmerHash);
        std::cerr << "Index succesfully loaded from the disk" << std::endl;
        printMemUsage();
    }else {
        if (options->getFileName(IGNORE_FILE).size())
            populateIgnoreKmerHash(ignoreKmerHash);

        Timer t("Kmer Initial Counting: ");
        t.start();
        hashSize = 2*collectKmersParallel(ignoreKmerHash); //50%
        printMemUsage();
        t.stop();

        if (!index_write_flag) {
            Timer t2("populate hash: ");
            t2.start();
            initializeHashForReference(ignoreKmerHash, hashSize);
            t2.stop();
        }
        delete [] KmerCountArray;
        printMemUsage();
        if (index_write_flag) {
            writeIndexToFile(ignoreKmerHash);
            std::cerr << "Index succesfully written to the disk" << std::endl;
            return;
        }
    }
    Timer t1("Store all Kmer with read positions: ");
    t1.start();
    insertReferencePositions(ignoreKmerHash);
    t1.stop();
    printMemUsage();
}

int SmerMap::getNumLocks()
{
    int numThreads = options->getNumThreads();

    if (numThreads==1)
        return 0;
    else
        return 1000000;
}

void SmerMap::rehashSmers()
{
    int n=0;
    uint64_t numberOfSmers=2*currHashTabSize;
    /* Set the size of the hash table to the numberofSmers. */
    for (n=0; n<450; ++n)
    {
        if (hashTableSize[n] > numberOfSmers)
           break;
    }

    currHashTabSize = hashTableSize[n];
    if (n)
        prevHashTabSize = hashTableSize[n-1];
    else
        prevHashTabSize = 1;

    delete [] KmerCountArray;
    KmerCountArray = new uint64_t[currHashTabSize];
    #pragma omp parallel for num_threads(options->getNumThreads())
    for (uint64_t i=0; i < currHashTabSize; i++)
        KmerCountArray[i] = HASH_INITIAL_VALUE;

    std::cerr << "Rehashed, New Hash Table Size: " << currHashTabSize << std::endl;
}

SmerMap::SmerMap (SequenceReader *readerObj, Options *op)
{
    int n=0;
    reader=readerObj;
    options=op;
    numElements=0;

    initializeLocks();
    if (index_read_flag) 
        prevHashTabSize = 1;
    else {
        uint64_t numberOfReads=1000*reader->getNumberOfReads();
        for (n=0; n<450; ++n)
        {
            if (hashTableSize[n] > numberOfReads)
               break;
        }

        currHashTabSize = hashTableSize[n];
        if (n)
            prevHashTabSize = hashTableSize[n-1];
        else
            prevHashTabSize = 1;

        KmerCountArray = new uint64_t[currHashTabSize];
        #pragma omp parallel for num_threads(options->getNumThreads())
        for (uint64_t i=0; i < currHashTabSize; i++)
            KmerCountArray[i] = HASH_INITIAL_VALUE;
        std::cerr << "Current Hash Table Size: " << currHashTabSize << std::endl;
    }
}

std::vector<ReadHashNode> * SmerMap::getMapOfReads(uint64_t &key)
{
    return SmerArray[key].readMap;
}

void SmerMap::initializeLocks()
{
    int numLocks = getNumLocks();
    std::cerr << "Number of Locks:" << numLocks <<std::endl;
           
    if (numLocks) 
    { 
        myLocks = new omp_lock_t [numLocks];
        #pragma omp parallel for  
        for (int i=0; i<numLocks; i++)
            omp_init_lock(myLocks+i);
    }
}

void SmerMap::destroyLocks()
{
    int numLocks = getNumLocks();
    if (numLocks) 
    { 
        #pragma omp parallel for  
        for (int i=0; i<numLocks; i++)
            omp_destroy_lock(myLocks+i);

        delete [] myLocks;
    }
}

SmerMap::~SmerMap()
{
    destroyLocks();
    delete [] SmerArray;
}

