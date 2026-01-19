#ifndef SWLU_SWDEF_H
#define SWLU_SWDEF_H

#include <stdint.h>

struct stack_frame {
    uint64_t :64;
    uint64_t :64;
    uint64_t :64;
    uint64_t :64;
    uint64_t :64;
    uint64_t :64;
    uint64_t :64;
    uint64_t :64;
    uint64_t :64;
    uint64_t :64;
    uint64_t :64;
    uint64_t :64;
    uint64_t :64;
    uint64_t :64;
    uint64_t :64;
    uint64_t fp;
    uint64_t :64;
    uint64_t :64;
    uint64_t :64;
    uint64_t :64;
    uint64_t :64;
    uint64_t :64;
    uint64_t t8;
    uint64_t t9;
    uint64_t :64;
    uint64_t :64;
    uint64_t ra;
    uint64_t t12;
    uint64_t :64;
    uint64_t gp;
    uint64_t sp;
    uint64_t zero;
    uint64_t pc;
};

#define REG(x) #x
#define REG_FP 15
#define REG_T8 22
#define REG_T9 23
#define REG_RA 26
#define REG_T12 27
#define REG_GP 29
#define REG_SP 30
#define REG_ZERO 31
#define SREG_FP "15"
#define SREG_T8 "22"
#define SREG_T9 "23"
#define SREG_RA "26"
#define SREG_T12 "27"
#define SREG_GP "29"
#define SREG_SP "30"
#define SREG_ZERO "31"

// LSB
#define __BITFIELD_FIELD(field, more)   \
    more                                \
    field;

struct instr_s_format { // storage
    __BITFIELD_FIELD(unsigned int opcode : 6,
    __BITFIELD_FIELD(unsigned int ra : 5,
    __BITFIELD_FIELD(unsigned int rb : 5,
    __BITFIELD_FIELD(signed int disp : 16,
    ))));
};

struct instr_b_format { // storage
    __BITFIELD_FIELD(unsigned int opcode : 6,
    __BITFIELD_FIELD(unsigned int ra : 5,
    __BITFIELD_FIELD(signed int disp : 21,
    )));
};

struct instr_c2_r_format { // calculation, 2 src, register
    __BITFIELD_FIELD(unsigned int opcode : 6,
    __BITFIELD_FIELD(unsigned int ra : 5,
    __BITFIELD_FIELD(unsigned int rb : 5,
    __BITFIELD_FIELD(unsigned int sbz : 3,
    __BITFIELD_FIELD(unsigned int func : 8,
    __BITFIELD_FIELD(unsigned int rc : 5,
    ))))));
};

struct instr_c2_i_format { // calculation, 2 src, immediate
    __BITFIELD_FIELD(unsigned int opcode : 6,
    __BITFIELD_FIELD(unsigned int ra : 5,
    __BITFIELD_FIELD(unsigned int imm : 8,
    __BITFIELD_FIELD(unsigned int func : 8,
    __BITFIELD_FIELD(unsigned int rc : 5,
    )))));
};

union sw_instruction {
    uint32_t word;
    struct instr_s_format s_format;
    struct instr_b_format b_format;
    struct instr_c2_r_format c2_r_format;
    struct instr_c2_i_format c2_i_format;
};

#endif
