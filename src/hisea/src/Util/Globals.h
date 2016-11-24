/* ================================================================= *
 *  Global.h: Header file with extern definitions                    *
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

#ifndef GLOBALS_H
#define GLOBALS_H

#include <cstdint>

#ifdef UINT64 
  typedef uint64_t UINT_TYPE;
  #define CHARACTERS_IN_UINT    32
#elif UINT32 
  typedef uint32_t UINT_TYPE;
  #define CHARACTERS_IN_UINT    16
#else //UINT16
  typedef uint16_t UINT_TYPE;
  #define CHARACTERS_IN_UINT    8
#endif

#define KMER_COUNT_MASK		0xFFFFFF /* This is the max it can go - 16,777, 214 (0xFFFFFE) kmers, this must change if max kmerSize goes beyond 20 */
#define HASH_INITIAL_VALUE	0xFFFFFFFFFFFFFFFF 
#define MAX_ANCHOR_SIZE_FOR_SCORE	14

typedef enum {FORWARD=0, REVERSE=1} ReadDirection;

typedef enum {FASTA_REF, FASTA_QUERY, FASTQ_REF, FASTQ_QUERY, OUTPUT, INDEX_FILE, IGNORE_FILE, UNKNOWN} FileType;

typedef enum {KMER_20=0, KMER_19, KMER_18, KMER_17, KMER_16, KMER_15, KMER_14, KMER_13, KMER_12, KMER_11, KMER_10, KMER_9, KMER_8}KmerType;


/* Global mask for biti reverse manipulation */
extern int self_flag;
extern int query_flag;
extern int help_flag;
extern int index_read_flag;
extern int index_write_flag;
extern uint16_t reverse_mask[8];
extern  uint16_t preComputedReverseComplement[65536];

extern UINT_TYPE global_mask_left[CHARACTERS_IN_UINT];
extern UINT_TYPE global_mask_right[CHARACTERS_IN_UINT];
extern UINT_TYPE clearBits_mask[CHARACTERS_IN_UINT];
extern UINT_TYPE kmerSeedArray[13]; 
extern uint32_t kmerSizeArray[13]; 
extern uint64_t hashTableSize[450];
#endif
