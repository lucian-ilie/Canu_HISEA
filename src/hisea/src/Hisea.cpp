/* ================================================================= *
 *  hisea.cpp: Main source file for hisea                            *
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

#include <iostream>
#include <Globals.h>
#include <Util.h>
#include <KmerMap.h>
#include <Alignment.h>
#include <Options.h>
#ifdef DEBUG
#include <signal.h>
#include <execinfo.h>
#endif

using namespace std;
#ifdef DEBUG
void handler(int sig) 
{
    void *array[30];
    size_t size;
    size = backtrace(array, 30);
    cerr << "Error: signal: " << sig << std::endl;
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(EXIT_SUCCESS);
}
#endif

static struct option long_options[] =
        {
          /* These options set a flag. */
          {"self",           no_argument,       &self_flag, 1},
          {"help",           no_argument,       &help_flag, 1},
          /* These options don’t set a flag.
           * We distinguish them by their indices. */
          {"minStoreLen",               required_argument, 0, 'a'},
          {"minOlapLen",                required_argument, 0, 'o'},
          {"minMatches",                required_argument, 0, 'b'},
          {"kmerLen",                   required_argument, 0, 'd'},
          {"smallKmerLen",              required_argument, 0, 'r'},
          {"filterCount",  		required_argument, 0, 'e'},
          {"threads",  			required_argument, 0, 'f'},
          {"errorRate",                 required_argument, 0, 'g'},
          {"ref",                       required_argument, 0, 'h'},
          {"query",                     required_argument, 0, 'i'},
          {"index_write",               required_argument, 0, 'j'},
          {"index_read",                required_argument, 0, 'k'},
          {"ignore",                    required_argument, 0, 'l'},
          {"maxKmerHits",               required_argument, 0, 'p'},
          {"maxShift",                  required_argument, 0, 'q'},
          {0, 0, 0, 0}
        };

int main(int argc, char** argv)
{
    SequenceReader refReader, queryReader;
    Alignment alignmentObj;
    int opt=0;
   
    Options options;

    if (argc==1){
        Usage();
        exit(EXIT_SUCCESS);
    }

    // Read Options
    while ((opt = getopt_long(argc, argv, "", long_options, NULL)) != EOF)
    {
 
        switch (opt) {
            case 0:
                break;
            case 'a':
                options.setMinStoreLength(optarg);
                break;
            case 'b':
                options.setMinMatches(optarg);
                break;
            case 'd':
                options.setStartSeedIndex(optarg);
                break;
            case 'e':
                options.setFilterCount(optarg);
                break;
            case 'f':
                options.setThreads(optarg);
                break;
            case 'g':
                options.setErrorRate(optarg);
                break;
            case 'h':
                options.setFileName(optarg, FASTA_REF);
                break;
            case 'i':
                query_flag=1;
                options.setFileName(optarg, FASTA_QUERY);
                break;
            case 'j':
                options.setFileName(optarg, INDEX_FILE);
                index_write_flag=1;
                break;
            case 'k':
                options.setFileName(optarg, INDEX_FILE);
                index_read_flag=1;
                break;
            case 'l':
                options.setFileName(optarg, IGNORE_FILE);
                break;
            case 'o':
                options.setMinOlapLength(optarg);
                break;
            case 'p':
                options.setMaxKmerHits(optarg);
                break;
            case 'q':
                options.setMaxShift(optarg);
                break;
            case 'r':
                options.setSmallKmerIndex(optarg);
                break;
            case ':':
                std::cerr << "Missing argument for option." << std::endl;
                Usage();
                exit(EXIT_FAILURE);
                break;
            case '?':
            default: 
                std::cerr << "Invalid Option." << std::endl;
                Usage();
                exit(EXIT_FAILURE);
        }
    }

    if (help_flag){
        options.help();
        exit(EXIT_SUCCESS);
    }
#ifdef DEBUG   
    signal(SIGSEGV, handler); 
#endif

    // Ensure options intigrity
    options.checkOptions();

    // PreCompute the reverse complement for all combinations of uint16_t
    preComputeReverseComplement();

    omp_set_num_threads(options.getNumThreads());
    cerr << "***************   Reading sequences   ***************" << endl;

    refReader.Initialize(options);
    refReader.readSequences(FASTA_REF);
    if (query_flag)
    {
        queryReader.Initialize(options);
        if (self_flag)
        {
            cerr << "[WARNING] Self mode option is passed with a query file" << endl;
            cerr << "[WARNING] Self alignment will be generated with reference sequences only" << endl;
            queryReader.dupSequenceReader(refReader);
        }
        queryReader.readSequences(FASTA_QUERY);
    }
    printMemUsage();
 
    cerr << "***************   Done Reading sequences   ***************" << endl;
    cerr << "***************   Finding Alignments   ***************" << endl;
    if (!query_flag)
    {
        alignmentObj.Initialize(refReader, refReader, options);
    }else
        alignmentObj.Initialize(refReader, queryReader, options);
    alignmentObj.findAllAlignments();
    cerr << "***************   Done Finding Alignments   ***************" << endl;
}

