#ifndef _COMP_H
#define _COMP_H

#include <stdint.h>

/* constant time a < b ? 1 : 0 */
static inline uint64_t ct_lt(uint64_t x, uint64_t y) { return (x - y) >> 63; }

/* return (x >= y) */
static inline uint64_t ct_ge(uint64_t x, uint64_t y) {
  return 1 ^ ((x - y) >> 63);
}

/* return (x == y) */
static inline uint64_t ct_eq(uint64_t x, uint64_t y) {
  return 1 ^ (((x - y) | (y - x)) >> 63);
}
#endif
