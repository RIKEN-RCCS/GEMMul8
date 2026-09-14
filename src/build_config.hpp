#pragma once
// Defaults for builds outside the top-level Makefile. 
// Selection is supplied with -D flags; the public enum and public header declarations are unchanged.
#ifndef GEMMUL8_BUILD_INT8
#define GEMMUL8_BUILD_INT8 1
#endif
#ifndef GEMMUL8_BUILD_FP8
#define GEMMUL8_BUILD_FP8 1
#endif
#ifndef GEMMUL8_BUILD_OP_gemm
#define GEMMUL8_BUILD_OP_gemm 1
#endif
#ifndef GEMMUL8_BUILD_OP_symm
#define GEMMUL8_BUILD_OP_symm 1
#endif
#ifndef GEMMUL8_BUILD_OP_syrk
#define GEMMUL8_BUILD_OP_syrk 1
#endif
#ifndef GEMMUL8_BUILD_OP_syr2k
#define GEMMUL8_BUILD_OP_syr2k 1
#endif
#ifndef GEMMUL8_BUILD_OP_syrkx
#define GEMMUL8_BUILD_OP_syrkx 1
#endif
#ifndef GEMMUL8_BUILD_OP_hemm
#define GEMMUL8_BUILD_OP_hemm 1
#endif
#ifndef GEMMUL8_BUILD_OP_herk
#define GEMMUL8_BUILD_OP_herk 1
#endif
#ifndef GEMMUL8_BUILD_OP_her2k
#define GEMMUL8_BUILD_OP_her2k 1
#endif
#ifndef GEMMUL8_BUILD_OP_herkx
#define GEMMUL8_BUILD_OP_herkx 1
#endif
#ifndef GEMMUL8_BUILD_OP_trmm
#define GEMMUL8_BUILD_OP_trmm 1
#endif
#ifndef GEMMUL8_BUILD_OP_trsm
#define GEMMUL8_BUILD_OP_trsm 1
#endif
#ifndef GEMMUL8_BUILD_OP_trtrmm
#define GEMMUL8_BUILD_OP_trtrmm 1
#endif
#if !GEMMUL8_BUILD_INT8 && !GEMMUL8_BUILD_FP8
#error "At least one Ozaki-II arithmetic backend must be enabled"
#endif
