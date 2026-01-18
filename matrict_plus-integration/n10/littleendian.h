#ifndef _LITTLEENDIAN_H
#define _LITTLEENDIAN_H

#include <stdint.h>

static inline uint64_t load_48(const unsigned char *x) {
  return ((uint64_t)(*x)) | (((uint64_t)(*(x + 1))) << 8) |
         (((uint64_t)(*(x + 2))) << 16) | (((uint64_t)(*(x + 3))) << 24) |
         (((uint64_t)(*(x + 4))) << 32) | (((uint64_t)(*(x + 5))) << 40);
}

static inline uint64_t load_40(const unsigned char *x) {
  return ((uint64_t)(*x)) | (((uint64_t)(*(x + 1))) << 8) |
         (((uint64_t)(*(x + 2))) << 16) | (((uint64_t)(*(x + 3))) << 24) |
         (((uint64_t)(*(x + 4))) << 32);
}

static inline uint64_t load_32(const unsigned char *x) {
  return ((uint64_t)(*x)) | (((uint64_t)(*(x + 1))) << 8) |
         (((uint64_t)(*(x + 2))) << 16) | (((uint64_t)(*(x + 3))) << 24);
}

static inline uint64_t load_24(const unsigned char *x) {
  return ((uint64_t)(*x)) | (((uint64_t)(*(x + 1))) << 8) |
         (((uint64_t)(*(x + 2))) << 16);
}

static inline uint64_t load_16(const unsigned char *x) {
  return ((uint64_t)(*x)) | (((uint64_t)(*(x + 1))) << 8);
}

static inline void store_32(unsigned char *x, const uint64_t c) {
  *x = c;
  *(x + 1) = c >> 8;
  *(x + 2) = c >> 16;
  *(x + 3) = c >> 24;
}

static inline void store_16(unsigned char *x, const uint64_t c) {
  *x = c;
  *(x + 1) = c >> 8;
}

#endif
