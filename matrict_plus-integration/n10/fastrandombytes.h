#ifndef FASTRANDOMBYTES_H
#define FASTRANDOMBYTES_H

void fastrandombytes_setseed_pub(const unsigned char *randomness);
void fastrandombytes_setseed_prv(const unsigned char *randomness);
void fastrandombytes_setseed_tmp(const unsigned char *randomness);
void fastrandombytes_setseed_ch(const unsigned char *randomness);

void fastrandombytes_pub(unsigned char *r, unsigned long long rlen);
void fastrandombytes_prv(unsigned char *r, unsigned long long rlen);
void fastrandombytes_tmp(unsigned char *r, unsigned long long rlen);
void fastrandombytes_ch(unsigned char *r, unsigned long long rlen);

void fastrandombytes_setiv_mat1();
void fastrandombytes_setiv_mat2();
void fastrandombytes_setiv_mat3();

#endif
