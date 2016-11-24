/* ================================================================= *
 *  Options.cpp: Source file for command line options processing     *
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

#include <Options.h>
#include <iostream>
#include <sys/stat.h>
#include <boost/algorithm/string.hpp>

bool vec_pairCompare(const std::pair<std::string, FileType>& firstElem, const std::pair<std::string, FileType>& secondElem) {
  return firstElem.first < secondElem.first;
}

Options::Options() 
{
    minStoreLength = 100;
    minOlapLength = 100;
    minMatches = 3;
    threads = 1;
    startSmerIndex=4;
    filterCount=2;
    errorRate=0.15;
    maxKmerHits=10000;
    maxShift=0.20;
    smallKmerIndex=8;
}

void Options::setMaxKmerHits(char *str)
{
    maxKmerHits=atoi(str);
}

int &Options::getMaxKmerHits()
{
    return maxKmerHits;
}

void Options::setFilterCount(char *str)
{
    filterCount=atoi(str);
}

void Options::setMaxShift(char *str)
{
    maxShift=std::stof(str, NULL);
}

float& Options::getMaxShift()
{
    return maxShift;
}

void Options::setStartSeedIndex(char *str)
{
    int seed = atoi(str);
    if (seed == 20)
        startSmerIndex=KMER_20;
    else if (seed == 19)
        startSmerIndex=KMER_19;
    else if (seed == 18)
        startSmerIndex=KMER_18;
    else if (seed == 17)
        startSmerIndex=KMER_17;
    else if (seed == 16)
        startSmerIndex=KMER_16;
    else if (seed == 15)
        startSmerIndex=KMER_15;
    else if (seed == 14)
        startSmerIndex=KMER_14;
    else if (seed == 13)
        startSmerIndex=KMER_13;
    else if (seed == 12)
        startSmerIndex=KMER_12;
    else if (seed == 11)
        startSmerIndex=KMER_11;
    else if (seed == 10)
        startSmerIndex=KMER_10;
    else if (seed == 9)
        startSmerIndex=KMER_10;
    else if (seed == 8)
        startSmerIndex=KMER_10;
    else
        startSmerIndex=100000;
}

void Options::setSmallKmerIndex(char *str)
{
    int seed = atoi(str);
    if (seed == 20)
        smallKmerIndex=KMER_20;
    else if (seed == 19)
        smallKmerIndex=KMER_19;
    else if (seed == 18)
        smallKmerIndex=KMER_18;
    else if (seed == 17)
        smallKmerIndex=KMER_17;
    else if (seed == 16)
        smallKmerIndex=KMER_16;
    else if (seed == 15)
        smallKmerIndex=KMER_15;
    else if (seed == 14)
        smallKmerIndex=KMER_14;
    else if (seed == 13)
        smallKmerIndex=KMER_13;
    else if (seed == 12)
        smallKmerIndex=KMER_12;
    else if (seed == 11)
        smallKmerIndex=KMER_11;
    else if (seed == 10)
        smallKmerIndex=KMER_10;
    else
        smallKmerIndex=100000;
}
        
void Options::setMinMatches(char *str)
{
    minMatches=atoi(str);
}
        
void Options::setMinOlapLength(char *str)
{
    minOlapLength=atoi(str);
}

void Options::setMinStoreLength(char *str)
{
    minStoreLength=atoi(str);
}

void Options::setErrorRate(char *str)
{
    errorRate=std::stof(str);
}

void Options::setThreads(char *str)
{
    threads=atoi(str);
}

void Options::setFileName(char *str, FileType type)
{
    if (type==FASTA_REF) {
        refFile.assign(str);
        addFiles(str, type);
    }
    else if (type==FASTA_QUERY) {
        queryFile.assign(str);
        addFiles(str,type);
    }
    else if (type==IGNORE_FILE) {
        ignoreFile.assign(str);
    }
    else if (type==INDEX_FILE) {
        indexFile.assign(str);
    }
}

int& Options::getNumThreads()
{
    return threads;
}

int& Options::getFilterCount()
{
    return filterCount;
}

int& Options::getSmerSeedIndex()
{
    return startSmerIndex;
}

int& Options::getSmallKmerIndex()
{
    return smallKmerIndex;
}

float& Options::getErrorRate()
{
    return errorRate;
}

std::string& Options::getFileName(int type)
{
    if (type==INDEX_FILE)
        return indexFile;
    else if(type==IGNORE_FILE)
        return ignoreFile;
    else if (type==1)
        return refFile;
    else if (type==2)
        return queryFile;
    else {
        std::cerr<< "Invalid type passed" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void Options::checkOptions()
{
    if (errorRate < 0 && errorRate > 1) {
        std::cerr<< "Invalid Error rate value. Valid values are between 0.0 and 1.0" << std::endl;
        Usage();
        exit(EXIT_FAILURE);
    }

    if (maxShift < 0 && maxShift > 1) {
        std::cerr<< "Invalid maxShift value. Valid values are between 0.0 and 1.0" << std::endl;
        Usage();
        exit(EXIT_FAILURE);
    }

    if (maxKmerHits < 50) {
        std::cerr<< "MaxKmerHits cannot be less than 50" << std::endl;
        Usage();
        exit(EXIT_FAILURE);
    }

    if (minOlapLength < 1) {
        std::cerr<< "Minimum Read length cannot be less than 1" << std::endl;
        Usage();
        exit(EXIT_FAILURE);
    }

    if (minStoreLength < 1) {
        std::cerr<< "Minimum Store length cannot be less than 1" << std::endl;
        Usage();
        exit(EXIT_FAILURE);
    }
    
    if (minStoreLength < minOlapLength) {
        minOlapLength += kmerSizeArray[getSmerSeedIndex()];
        minStoreLength = minOlapLength;
        std::cerr<< "Minimum Store length cannot be smaller than minOlapLen" << std::endl;
        std::cerr<< "Setting Minimum Store length equal to minOlapLen" << std::endl;
    }

    if (minMatches < 1) {
        std::cerr<< "Minimum Matches cannot be less than 1" << std::endl;
        Usage();
        exit(EXIT_FAILURE);
    }

    if (threads < 1) {
        std::cerr<< "Number of threads cannot be less than 1" << std::endl;
        Usage();
        exit(EXIT_FAILURE);
    }

    if (startSmerIndex > 1000) {
        std::cerr<< "Invalid kmerLen value, allowed values are between 8 and 20 both inclusive" << std::endl;
        Usage();
        exit(EXIT_FAILURE);
    }

    if (smallKmerIndex > 1000) {
        std::cerr<< "Invalid smallKmerLen value, allowed values are between 88888888 and 20 both inclusive" << std::endl;
        Usage();
        exit(EXIT_FAILURE);
    }

    if (kmerSizeArray[smallKmerIndex] >= kmerSizeArray[startSmerIndex]) {
        std::cerr<< "The smallKmerLen cannot be more than kmerLen" << std::endl;
        Usage();
        exit(EXIT_FAILURE);
    }


    if (filterCount < 1 ) {
        std::cerr<< "Invalid value for filterCount, allowed values are >=1" << std::endl;
        Usage();
        exit(EXIT_FAILURE);
    }


    if (index_write_flag && queryFile.size())
    {
        std::cerr<< "[ERROR] Queryfile is not required when writing index" << std::endl;
        Usage();
        exit(EXIT_FAILURE);
    }

    if (index_read_flag && index_write_flag)
    {
        std::cerr<< "[ERROR] Index read/write flags are mutually exclusive" << std::endl;
        Usage();
        exit(EXIT_FAILURE);
    }
    
    if (query_flag || index_write_flag || index_read_flag) 
        filterCount=1;

    if (!refFile.size()) {
        std::cerr<< "Reference file/directory is a mandatory argument" << std::endl;
        Usage();
        exit(EXIT_FAILURE);
    }

    if (!self_flag && !index_write_flag && !queryFile.size()) {
        std::cerr<< "When not is self mode, query file/directory is a mandatory argument" << std::endl;
        Usage();
        exit(EXIT_FAILURE);
    }
    

    if (index_read_flag || index_write_flag) 
    {
        char *pch=NULL, *pch1=NULL;
        char str[1000];
        checkValidDir(indexFile.c_str());
        strcpy(str, refFile.c_str());
        pch = strtok (str,"/");
        while (pch != NULL)
        {
          pch1=pch;
          pch = strtok (NULL, "/");
        }
        strcpy(str, pch1);
        
        pch = strtok (str,".");
        
        indexFile += "/";
        indexFile += pch;
    }

    if (ignoreFile.size())
    {
        std::fstream fs;
        fs.open(ignoreFile, std::fstream::in);
        if (!fs.is_open()) {
            std::cerr << "[ERROR] - Unable to open ignore File: " << ignoreFile << std::endl;
            exit(EXIT_FAILURE);
        }
        fs.close();
    }
}

void Options::addFiles(char *str, FileType type)
{
    struct stat sb={0};

    stat(str, &sb);
    if (S_ISREG(sb.st_mode)) {
        vec_pair.push_back(std::make_pair(str, type));
        return;
    }else {  // It is a directory
        char path[20000];
        DIR *dir;
        struct dirent *ent;
        if ((dir = opendir (str)) != NULL) 
        {
          while ((ent = readdir (dir)) != NULL) 
          {
              if(!strcmp(ent->d_name, ".") || !strcmp(ent->d_name, ".."))
                  continue;

              strcpy(path, str);
              strcat(path,"/");
              strcat(path,ent->d_name);
              stat(path, &sb);
              if (S_ISREG(sb.st_mode)) {
                  vec_pair.push_back(std::make_pair(path, type));      
              }
          }
          closedir (dir);
        } else {
          /* could not open directory */
          std::cerr<< "[ERROR] -  Could not access file/directory: " << str << std::endl;
          exit(EXIT_FAILURE);
        }
    }
    std::sort(vec_pair.begin(), vec_pair.end(), vec_pairCompare); 
}


void Options::checkValidDir(const char *str)
{
    struct stat sb={0};

    stat(str, &sb);
    if (!S_ISREG(sb.st_mode)) {
        DIR *dir;
        if ((dir = opendir (str)) != NULL) 
        {
          closedir (dir);
        } else {
          /* could not open directory */
          std::cerr<< "[ERROR] -  Could not access directory path: " << str << std::endl;
          exit(EXIT_FAILURE);
        }
    }else {  // It is not a valid path
        std::cerr<< "[ERROR] -  Seems to be file, need a path/directory: " << str << std::endl;
        exit(EXIT_FAILURE);
    }
}

void Options::help()
{
   std::cerr <<
           "USAGE: hisea  [--self] [--kmerLen  <int>] [--minStoreLen <int>] [--minOlapLen <int>] --ref <file/directory> --query <file/directory>" << std::endl
       <<  std::endl
       <<  "HISEA is an efficient all-vs-all long read alignment program." << std::endl
       <<  "It is designed and tested for SMRT sequecning technology, but" << std::endl
       <<  "should work well for Oxford Nanopore as well." << std::endl
       <<  std::endl
       <<  "Options:" << std::endl
       <<  "--help         no_arg     Print this help message" << std::endl
       <<  "--self         no_arg     Align set of reads with itself. Use only --ref option" << std::endl
       <<  "--ref          <file/dir> The name of the reference fasta/fastq file or directory" << std::endl
       <<  "                          containing these files" << std::endl
       <<  "--query        <file/dir> The name of the query fasta/fastq file or directory" << std::endl
       <<  "                          containing these files" << std::endl
       <<  "--ignore       <file>     The file containing kmers to be ignored" << std::endl
       <<  "--index_write  <dir>      The directory where index is stored" << std::endl
       <<  "--index_read   <dir>      The directory from whcih index is read" << std::endl
       <<  "--kmerLen      <int>      This is the kmer length used for initial hashing." << std::endl
       <<  "                          The possible values are 10-20, both inclusive" << std::endl
       <<  "                          default=16" << std::endl
       <<  "--smallkmerLen <int>      This is the kmer length used during alignment extension." << std::endl
       <<  "                          The possible values are 10-20, both inclusive" << std::endl
       <<  "                          default=12" << std::endl
       <<  "--filterCount  <int>      This is used for initial filtering. This must be set to 1," << std::endl
       <<  "                          if index is created with split set of reads." << std::endl
       <<  "                          default=2" << std::endl
       <<  "--threads      <int>      Number of threads for parallel run" << std::endl
       <<  "                          default=1" << std::endl
       <<  "--minOlapLen   <int>      Minimum Overlap length" << std::endl
       <<  "                          default=100" << std::endl
       <<  "--minStoreLen  <int>      Minimum read length. Overlap of two smaller reads is ignored" << std::endl
       <<  "                          default=100" << std::endl
       <<  "--minMatches   <int>      Minimum number of matches to be considered for alignment" << std::endl
       <<  "                          default=3" << std::endl
       <<  "--maxKmerHits  <int>      The maximum number of repeat kmers" << std::endl
       <<  "                          default=10000" << std::endl
       <<  "--errorRate    <float>    Error rate. Valid values are 0-1 percent" << std::endl
       <<  "                          default=0.15" << std::endl
       <<  "--maxShift     <float>    The value of shift to accomodate indels. Valid values" << std::endl
       <<  "                          are 0-1 percent" << std::endl
       <<  "                          default=0.20" << std::endl;
}

