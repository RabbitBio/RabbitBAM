//
// Created by 赵展 on 2021/4/7.
//

#ifndef BAMSTATUS_CONFIG_H
#define BAMSTATUS_CONFIG_H
#define MAXLEN 512
typedef char int8;
typedef unsigned char uint8;
typedef short int16;
typedef unsigned short uint16;
typedef int int32;
typedef unsigned int uint32;
typedef long long int64;
typedef unsigned long long uint64;
// A 1 C 3 T 4 N 6 G 7
static uint8 valAGCT[]={0,0,0,1,2,0,0xC0,3};
static int8 DupvalAGCT[]={-1,0,-1,1,2,-1,-1,3};

static uint8 reverseAGCTN[]={2,3,0,1,0xc0,0xc0,0xc0,0xc0};
static uint8 lenAGCTN[]={2,2,2,2,2,2,0,2};
static int8 charAGCTN[] = {'A','C','T','G','N','N','N','N'};
static uint8 intAGCTN[] = {0,1,2,3,0xc0,0xc0,0xc0,0xc0};
const int8 seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };
const uint8 StatusBase[16] = {0, 65, 67, 0, 71, 0, 0, 0, 84, 0, 0, 0, 0, 0, 0, 78};
const uint8 StatusBaseRever[16] = {0, 84, 71, 0, 67, 0, 0, 0, 65, 0, 0, 0, 0, 0, 0, 78};

#endif //BAMSTATUS_CONFIG_H
