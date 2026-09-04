/**
 * GEMMul8 (GEMMulate)
 * -------------------
 * GEMMul8 emulates high-precision BLAS-like matrix operations using INT8/FP8
 * matrix engines based on the Ozaki Scheme II.
 *
 * The public API is designed to be close to the corresponding cuBLAS/hipBLAS
 * routines.  For example, gemm() follows the cuBLAS/hipBLAS Xgemm argument
 * convention, with additional GEMMul8-specific arguments.
 *
 * cuBLASLt/hipBLASLt-handle variants are also provided for supported routines.
 * Their argument lists are similar to the corresponding cuBLAS/hipBLAS variants,
 * except that the Lt handle is used and a stream argument is appended at the end.
 *
 * See the corresponding header file for the detailed specification of each function:
 *
 *   gemm.hpp     : general matrix-matrix multiplication
 *   symm.hpp     : symmetric matrix-matrix multiplication
 *   syrk.hpp     : symmetric rank-k update
 *   syr2k.hpp    : symmetric rank-2k update
 *   syrkx.hpp    : symmetric rank-k update with two input matrices
 *   trmm.hpp     : triangular matrix-matrix multiplication
 *   hemm.hpp     : Hermitian matrix-matrix multiplication
 *   herk.hpp     : Hermitian rank-k update
 *   her2k.hpp    : Hermitian rank-2k update
 *   herkx.hpp    : Hermitian rank-k update with two input matrices
 *   trsm.hpp     : triangular solve with multiple right-hand sides
 *   trtrmm.hpp   : triangular-by-triangular matrix multiplication
 *   config.hpp   : handle-local execution configuration
 *   worksize.hpp : required workspace size for the specified GEMMul8 operation
 *
 * Performance note:
 *
 *   For better performance of the cuBLAS/hipBLAS-handle routines,
 *   set a sufficiently large workspace using cublasSetWorkspace / hipblasSetWorkspace,
 *   e.g., 32 MiB = (32 << 20) Bytes, for the input BLAS handle before invoking the routines.
 *
 *   This note applies to the cuBLAS/hipBLAS-handle routines,
 *   not to the cuBLASLt/hipBLASLt-handle variants.
 * 
 * GEMMul8 version information:
 * 
 *   See version.hpp.
 * 
 * Developer:
 *   Yuki Uchino
 *   RIKEN Center for Computational Science, Japan
 *
 * License:
 *   MIT License
 */

#pragma once

#include "types.hpp"
#include "version.hpp"
#include "config.hpp"
#include "worksize.hpp"
#include "gemm.hpp"
#include "symm.hpp"
#include "syrk.hpp"
#include "syr2k.hpp"
#include "syrkx.hpp"
#include "trmm.hpp"
#include "hemm.hpp"
#include "herk.hpp"
#include "her2k.hpp"
#include "herkx.hpp"
#include "trsm.hpp"
#include "trtrmm.hpp"
