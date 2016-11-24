/* ================================================================= *
 *  SequenceReader.cpp: Source file with methods for read processing *
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

#include <SequenceReader.h>
#include <vector>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <Timer.h>

/*
 * Member functions for Sequence Class
 */

Sequence::Sequence() 
{ //Constructor
        size=0;
        seq=NULL;
}

Sequence::Sequence(uint32_t s, UINT_TYPE *se) 
{ //Constructor
        size=s;
        seq=se;
}

bool Sequence::operator () (const Sequence &seqObj1, const Sequence &seqObj2) const
{      
    if (seqObj1.size < seqObj2.size){
        return true;
    }else if (seqObj1.size == seqObj2.size){
        uint32_t seqMaxIndex = char2byteSize(seqObj1.size);
        for (uint32_t i=0; i<=seqMaxIndex; i++){
            if (seqObj1.seq[i] < seqObj2.seq[i])
                return true;
            else if (seqObj1.seq[i] > seqObj2.seq[i])
                return false; 
        }
        return false;
    }else
        return false;
}    
 
bool Sequence::operator < (const Sequence &seqObj) const
{       
    if (this->size < seqObj.size){
        return true;
    }else if (this->size == seqObj.size){
        uint32_t seqMaxIndex = char2byteSize(seqObj.size);
        for (uint32_t i=0; i<=seqMaxIndex; i++){
            if (this->seq[i] < seqObj.seq[i])
                return true;
            else if (this->seq[i] > seqObj.seq[i])
                return false; 
        }
        return false;
    }else
        return false;
}    

bool Sequence::operator == (const Sequence &seqObj) const
{       
    if (this->size != seqObj.size){
        return false;
    }else{
        uint32_t seqMaxIndex = char2byteSize(seqObj.size);
        for (uint32_t i=0; i<=seqMaxIndex; i++){
            if (this->seq[i] != seqObj.seq[i])
                return false; 
        }
        return true;
    }
}    

/* 
 *  Returns true, if reverse complement has been swaped.
 */
bool Sequence::getReverseComplement()
{
    uint32_t size=char2byteSize(this->size);
    UINT_TYPE currObj=0;
    uint32_t offset=size2Offset(this->size), last=0;
    uint32_t reverseIndex=0, forwardIndex=size;    
   
    UINT_TYPE *seqReverse = (seq+size+1);

    while(reverseIndex <= size)
    {
        seqReverse[reverseIndex]=0;
        currObj = this->seq[forwardIndex];
        currObj >>= 2*(CHARACTERS_IN_UINT-offset-1);
        if (forwardIndex) {
            forwardIndex--;
            if (offset != CHARACTERS_IN_UINT-1) // Offset==31, means full object
               currObj |= this->seq[forwardIndex]<<2*(offset+1);
        }else {
            last=1;
        }

#ifdef UINT16

        if (last){
                seqReverse[reverseIndex] |= (preComputedReverseComplement[currObj] >> 2*(8-offset-1)) << 2*(8-offset-1);
        }else {
            seqReverse[reverseIndex] |= preComputedReverseComplement[currObj];  
        }

#elif UINT32

        if (last){
            if (offset < 8) {
                seqReverse[reverseIndex] |= static_cast<UINT_TYPE>(preComputedReverseComplement[currObj] >> 2*(8-offset-1)) << (16+2*(8-offset-1));
            }else {
                seqReverse[reverseIndex] |= (static_cast<UINT_TYPE>(preComputedReverseComplement[currObj & 0xffff]) << 16 |
                                            static_cast<UINT_TYPE>(preComputedReverseComplement[(currObj >> 16) & 0xffff] >> 2*(16-offset-1)) << 2*(16-offset-1));
            }
        }else {
            seqReverse[reverseIndex] |= (static_cast<UINT_TYPE>(preComputedReverseComplement[currObj & 0xffff]) << 16 |
                                        (static_cast<UINT_TYPE>(preComputedReverseComplement[(currObj >> 16) & 0xffff]))); 
        }

#elif UINT64

        if (last){
            if (offset < 8) {
                seqReverse[reverseIndex] |= static_cast<UINT_TYPE>(preComputedReverseComplement[currObj] >> 2*(8-offset-1)) << (48+2*(8-offset-1));
            }else if (offset < 16) {
                seqReverse[reverseIndex] |= (static_cast<UINT_TYPE>(preComputedReverseComplement[currObj & 0xffff]) << 48 |
                                            static_cast<UINT_TYPE>(preComputedReverseComplement[(currObj >> 16) & 0xffff] >> 2*(16-offset-1)) << (32+2*(16-offset-1)));
            }else if (offset < 24) {
                seqReverse[reverseIndex] |= (static_cast<UINT_TYPE>(preComputedReverseComplement[currObj & 0xffff]) << 48 |
                                            static_cast<UINT_TYPE>(preComputedReverseComplement[(currObj >> 16) & 0xffff]) << 32 |
                                            static_cast<UINT_TYPE>(preComputedReverseComplement[(currObj >> 32) & 0xffff] >> 2*(24-offset-1)) << (16+2*(24-offset-1)));
            }else{
                seqReverse[reverseIndex] |= (static_cast<UINT_TYPE>(preComputedReverseComplement[currObj & 0xffff]) << 48 |
                                            static_cast<UINT_TYPE>(preComputedReverseComplement[(currObj >> 16) & 0xffff]) << 32 |
                                            static_cast<UINT_TYPE>(preComputedReverseComplement[(currObj >> 32) & 0xffff]) << 16 |
                                            static_cast<UINT_TYPE>(preComputedReverseComplement[(currObj >> 48) & 0xffff] >> 2*(32-offset-1)) << 2*(32-offset-1));
            }
        }else {
            seqReverse[reverseIndex] |= (static_cast<UINT_TYPE>(preComputedReverseComplement[currObj & 0xffff]) << 48 |
                                        (static_cast<UINT_TYPE>(preComputedReverseComplement[(currObj >> 16) & 0xffff]) << 32) |
                                        (static_cast<UINT_TYPE>(preComputedReverseComplement[(currObj >> 32) & 0xffff]) << 16) |
                                        (static_cast<UINT_TYPE>(preComputedReverseComplement[(currObj >> 48) & 0xffff])));
        }
#endif
        reverseIndex++;
    }

    return false;
}

/*
 * Member functions for SequenceReader Class
 */
SequenceReader::SequenceReader()
{
}

void SequenceReader::Initialize(Options &op) 
{

    options = &op;
    readArraySize = 10000000;
    binReads = new Sequence[readArraySize];
    numOfReads=0;
    numBP=0;

}

SequenceReader::~SequenceReader() 
{
}

std::istream* SequenceReader::openInStream(const std::string &filename, std::ios_base::openmode mode)
{
    std::string c = filename.substr(filename.length()-std::strlen(".gz"), strlen(".gz"));
    if (!std::strcmp(filename.substr(filename.length()-std::strlen(".gz"), strlen(".gz")).c_str(), ".gz")) {
        igzstream* gfp = new igzstream(filename.c_str(), mode);
        if(!gfp->good()) {
            std::cerr << "[ERROR] - file open failed " << filename << std::endl;
            exit(EXIT_FAILURE);
        }
        return gfp;
    }else {
        std::ifstream* fp = new std::ifstream(filename.c_str(), mode);
        if (!fp->is_open()) {
            std::cerr << "[ERROR] - file open failed " << filename << std::endl;
            exit(EXIT_FAILURE);
        }
        return fp;
    }

}

void SequenceReader::checkReadArraySize()
{
     if (numOfReads+10 > readArraySize)
     {
         Sequence *tmpBinReads = new Sequence[readArraySize+10000000];
         #pragma omp parallel for
         for (uint64_t i=0; i < numOfReads; i++)
             tmpBinReads[i]=binReads[i];
         delete [] binReads;
         readArraySize += 10000000;
         binReads = tmpBinReads;
     }
}

uint64_t SequenceReader::getNumberOfReads()
{
    return numOfReads;
}

uint64_t SequenceReader::getNumberOfBP()
{
    return numBP;
}

uint64_t SequenceReader::getReadArraySize()
{
    return readArraySize;
}

void SequenceReader::removeSequences(uint64_t &index)
{
     delete [] binReads[index].seq;
     binReads[index].seq=NULL;
}

//Returns the index of the last element
uint64_t SequenceReader::appendNewSequence()
{
     checkReadArraySize();
     numOfReads++;     
     return numOfReads-1;
}

uint32_t SequenceReader::getSequenceSize(const uint64_t &index)
{
     return binReads[index].size;
}

UINT_TYPE* SequenceReader::getSequence(const uint64_t &index, ReadDirection &direction)
{
     uint32_t binSize = char2byteSize(binReads[index].size);
     if (direction==FORWARD)
         return binReads[index].seq;
     else if (direction == REVERSE){
         return &(binReads[index].seq[binSize+1]);
     }else{
        std::cerr << "[ERROR] - Invalid type passed" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void SequenceReader::updateFileType(FileType &ft)
{

    if (ft != FASTA_REF)
        return;
    std::istream *fp;
    std::string line;
    for(uint32_t i=0; i < options->vec_pair.size();i++)
    {
        fp = openInStream(options->vec_pair[i].first);
        while(getline(*fp, line).good()){
            if(line.empty())
                continue;

            if(line[0] == '>'){
              break;
            }else if(line[0] == '@'){
              if (options->vec_pair[i].second == FASTA_REF)
                  options->vec_pair[i].second = FASTQ_REF; 
              else
                  options->vec_pair[i].second = FASTQ_QUERY; 
              getline(*fp, line); // line 2 with sequence
              getline(*fp, line); // line 3 with a +
              getline(*fp, line); // line 4 with Quality score
              break;
            }else {
              std::cerr << "[INFO] - Unknown file format for file: " << options->vec_pair[i].first << std::endl;
              std::cerr << "          This file will be ignored from further processing" << options->vec_pair[i].first << std::endl;
                  options->vec_pair[i].second = UNKNOWN; 
              break;
            }
        }
        delete fp;
    }
}

void SequenceReader::dupSequenceReader(SequenceReader &sr)
{
        ReadDirection direction = FORWARD;
        numBP = sr.getNumberOfBP();
        numOfReads = sr.getNumberOfReads();
        readArraySize = sr.getReadArraySize();
        delete [] binReads;
        binReads = new Sequence[readArraySize];
        #pragma omp parallel for
        for (uint64_t i=0; i < numOfReads; i++)
        {
            binReads[i].size = sr.getSequenceSize(i); 
            binReads[i].seq = sr.getSequence(i, direction);
        }
}
        
void SequenceReader::readSequences(FileType ft)
{
    std::string line, variable,fileName;
    std::istream *fp;

    Timer t1("Reading and storing sequences");
    t1.start();
    updateFileType(ft);
    for(uint32_t i=0; i < options->vec_pair.size();i++)
    {
        if (options->vec_pair[i].second == UNKNOWN)
            continue;
        if (ft == FASTA_REF)
        {
            if (options->vec_pair[i].second != FASTA_REF && options->vec_pair[i].second != FASTQ_REF)
                continue;
        }else{
            if (options->vec_pair[i].second != FASTA_QUERY && options->vec_pair[i].second != FASTQ_QUERY)
                continue;
        }

        fp = openInStream(options->vec_pair[i].first);

        std::cerr << "[INFO] - Reading file: " << options->vec_pair[i].first << std::endl;
        if (options->vec_pair[i].second == FASTA_REF || options->vec_pair[i].second == FASTA_QUERY) { 
            readFASTAfile(*fp, ft);
        }else {   
            readFASTQfile(*fp, ft);
        }
        
        delete fp;
    }
    t1.stop();

    /* Remove duplicates and create reverse complement */
    Timer t2("Organizing sequences");
    t2.start();
    organizeSequences(ft);
    t2.stop();
}

void SequenceReader::readFASTAfile(std::istream &fp, FileType &ft)
{
    std::string strName, line, content;
    uint64_t goodSeq=0, badSeq=0;
    
    while(getline(fp, line).good() ){
        if(line[0] == '>'){
            if(!strName.empty()){ // Process what we read from the last entry
                if(convertToBinary(content, ft))
                    goodSeq++;
                else
                    badSeq++;
                strName.clear();
                content.clear();
            }
            if(!line.empty()){
                strName = line.substr(1);
            }
        } else if(!strName.empty()){
            content+=line;
        }
    }

    if(!strName.empty()){ // Process what we read from the last entry
        if(convertToBinary(content, ft))
            goodSeq++;
        else
            badSeq++;
        content.clear();
        strName.clear();
    }
    std::cerr << "\t[INFO] - Number of good sequences:" << goodSeq << std::endl;
    std::cerr << "\t[INFO] - Number of bad sequences (Sequence length smaller than minLength) :" << badSeq << std::endl;
}


void SequenceReader::readFASTQfile(std::istream &fp, FileType &ft)
{
    std::string content;
    uint64_t goodSeq=0, badSeq=0;
    
    while(getline(fp, content).good() ){
        if(content[0] == '@'){
            getline(fp, content); // line 2 with sequence
            if(convertToBinary(content, ft))
                goodSeq++;
            else
                badSeq++;
            getline(fp, content); // line 3 with a +
            getline(fp, content); // line 4 with Quality score
        }
    }
    std::cerr << "\t[INFO] - Number of good sequences:" << goodSeq << std::endl;
    std::cerr << "\t[INFO] - Number of bad sequences (Sequence length smaller than minLength) :" << badSeq << std::endl;
}

/*
 * Function converts a character sequence into an array of integers.
 * Input: character string
 * Output: array of integers, total number of bases
 */
bool SequenceReader::convertToBinary(const std::string &content, FileType &ft)
{
    uint64_t binReadsIndex=0, totalBases=0;
    uint32_t binSize=0; 
    uint64_t lastIndex=0; 
    int chooseLetter=0;
 
    binSize = char2byteSize(content.length());
    lastIndex = appendNewSequence();
    binReads[lastIndex].size = content.length();
    if (ft == FASTA_REF)
        binReads[lastIndex].seq = new UINT_TYPE[2*(binSize+1)]; // +1 because of integer division
    else
        binReads[lastIndex].seq = new UINT_TYPE[binSize+1]; // +1 because of integer division
    binReads[lastIndex].seq[binSize]=0; //only last index needs to be initialized

    /* Processing the sequences by encoding the base pairs into 2 bits. */
    for (std::string::const_iterator it=content.begin(); it!=content.end(); ++it)
    {
        switch(*it)
        {
            case 'A':
            case 'a':
                /* Left shift is gauranteed to fill vacant bits
                 * with zero's. The initialization will not be a 
                 * problem. Just to note - the right shift for our
                 * case (unsigned) is also gauratiees zero fill.
                 */
                binReads[lastIndex].seq[binReadsIndex] <<= 2;
                break;
            case 'C':
            case 'c':
                binReads[lastIndex].seq[binReadsIndex] <<= 2;
                binReads[lastIndex].seq[binReadsIndex] |= 1;
                break;
            case 'G':
            case 'g':
                binReads[lastIndex].seq[binReadsIndex] <<= 2;
                binReads[lastIndex].seq[binReadsIndex] |= 2;
                break;
            case 'T':
            case 't':
                binReads[lastIndex].seq[binReadsIndex] <<= 2;
                binReads[lastIndex].seq[binReadsIndex] |= 3;
                break;
            default:
                /* Reads with N are replaced with random bases */
                chooseLetter = rand() % 4;
                if (chooseLetter == 0)
                    binReads[lastIndex].seq[binReadsIndex] <<= 2;
                else if (chooseLetter == 1)
                {
                    binReads[lastIndex].seq[binReadsIndex] <<= 2;
                    binReads[lastIndex].seq[binReadsIndex] |= 1;
                }
                else if (chooseLetter == 2)
                {
                    binReads[lastIndex].seq[binReadsIndex] <<= 2;
                    binReads[lastIndex].seq[binReadsIndex] |= 2;
                }
                else
                {
                    binReads[lastIndex].seq[binReadsIndex] <<= 2;
                    binReads[lastIndex].seq[binReadsIndex] |= 3;
                }
        }
        totalBases++;
        if (!(totalBases%CHARACTERS_IN_UINT)){
            if (binReadsIndex!=binSize)
                binReadsIndex++;
        }
    }

    /* Align bases to the left for last uint64
     * Varibale binSize reused in dfferent context
     */
    if ((binSize = (totalBases%CHARACTERS_IN_UINT))) {
        binReads[lastIndex].seq[binReadsIndex] <<= 2*(CHARACTERS_IN_UINT-binSize);
    }
    numBP+=totalBases;
    return true; //Indicates a good sequence
}

void SequenceReader::organizeSequences(FileType &ft)
{
    uint64_t numOfReads=getNumberOfReads();
    
    if (!numOfReads){
        std::cerr << "[ERROR] - No valid reads found for alignment" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (ft == FASTA_REF)
    {
        #pragma omp parallel for 
        for (uint64_t i=0; i < numOfReads; i++)
        {
            binReads[i].getReverseComplement();
        }
    }
}


