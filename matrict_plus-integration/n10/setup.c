#include "setup.h"
#include "littleendian.h"
#include "param.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <libXKCP.a.headers/SimpleFIPS202.h>

#include "params.h"
#include "utils.h"
#include "params-GR19.h"
#include "base_gadgets.h"

#define PP_HASH_LEN (CRYPTO_BYTES + 9 + Q_BYTE + QHAT_BYTE)

void setup(unsigned char *out, unsigned char *seed) {
  unsigned char r[PP_HASH_LEN];

  memcpy(r, seed, CRYPTO_BYTES);
  r[CRYPTO_BYTES] = R;
  r[CRYPTO_BYTES + 1] = N;
  r[CRYPTO_BYTES + 2] = K;
  r[CRYPTO_BYTES + 3] = N_BAR;
  r[CRYPTO_BYTES + 4] = N_HAT;
  r[CRYPTO_BYTES + 5] = K_HAT;
  store_16(r + CRYPTO_BYTES + 6, D);
  STORE_Q(r + CRYPTO_BYTES + 8, Q);
  STORE_QHAT(r + CRYPTO_BYTES + 8 + Q_BYTE, Q);
  r[CRYPTO_BYTES + 8 + Q_BYTE + QHAT_BYTE] = W;

  SHAKE256(out, CRYPTO_BYTES, r, PP_HASH_LEN);
  
}
