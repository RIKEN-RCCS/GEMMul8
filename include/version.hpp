/**
 * Version information
 * -------------------
 * This header defines compile-time GEMMul8 version information.
 *
 * Public macros:
 *
 *   GEMMUL8_VERSION_MAJOR
 *   GEMMUL8_VERSION_MINOR
 *   GEMMUL8_VERSION_PATCH
 *   GEMMUL8_VERSION_CODE
 *   GEMMUL8_VERSION_STRING
 *
 * Public C++ constants:
 *
 *   gemmul8::version_major
 *   gemmul8::version_minor
 *   gemmul8::version_patch
 *   gemmul8::version_code
 *   gemmul8::version_string
 *
 * GEMMUL8_VERSION_CODE is encoded as:
 *
 *   major * 10000 + minor * 100 + patch
 */
#pragma once

#define GEMMUL8_VERSION_MAJOR  3
#define GEMMUL8_VERSION_MINOR  0
#define GEMMUL8_VERSION_PATCH  5
#define GEMMUL8_VERSION_CODE   (GEMMUL8_VERSION_MAJOR * 10000 + GEMMUL8_VERSION_MINOR * 100 + GEMMUL8_VERSION_PATCH)
#define GEMMUL8_VERSION_STRING "3.0.5"

// C++ constants corresponding to the public GEMMUL8_VERSION_* macros.
namespace gemmul8 {

inline constexpr int version_major          = GEMMUL8_VERSION_MAJOR;
inline constexpr int version_minor          = GEMMUL8_VERSION_MINOR;
inline constexpr int version_patch          = GEMMUL8_VERSION_PATCH;
inline constexpr int version_code           = GEMMUL8_VERSION_CODE;
inline constexpr const char *version_string = GEMMUL8_VERSION_STRING;

} // namespace gemmul8
