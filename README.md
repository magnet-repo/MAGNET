# MAGNET: MAsked Gaussian Now Efficient and Table-less 

This repository provides reference implementations of **MAGNET**, a table-less and efficient masked discrete Gaussian sampler, along with comparative implementation. It includes C and ARM Cortex-M4 implementations, integrations of the masked samplers into [MatRiCT+](https://gitlab.com/raykzhao/matrict_plus) private payment protocol, and leakage assessment artifacts for first-order side-channel evaluation.

## C Implementations

**MAGNET** can be compiled and executed with:

``make && ./MAGNET``

To see the runtime of unmasked **DP-DGS**, under the DP_DGS folder, execute:

``make && ./DP-DGS``

To see the runtime of **GR19-DGS**, under the GR19 folder, execute:

``make && ./GR19-DGS``

## MatRiCT+ Integration

The integration of **MAGNET** and **GR19-DGS** into the MatRiCT+ payment protocol can be found in the _matrict_plus-integration_ folder. To run the implementation, navigate to the _n10_ folder and execute:

``make && ./ringct``

By default, **MAGNET** is used as the discrete Gaussian sampler. The active sampler can be changed by modifying the corresponding #define directive in the _n10/param.h_ file.

## ARM Cortex-M4 Implementation

MAGNET ARM Cortex-M4 (Thumb-2 ASM) implementation can be found in _MAGNET_m4_ folder. The procedure for running the implementation using **nrfutil** is described in the README under this folder.

## Test Vector Leakage Assessment

The leakage assessment tests and the corresponding hardened implementations for each gadget and masked sampling algorithm are located in the _TVLA_ directory.
