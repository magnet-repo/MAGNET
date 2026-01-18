# MAGNET: MAsked Gaussian Now Efficient and Table-less 

**MAGNET** can be compiled and executed with:

``make && ./MAGNET``

To see the runtime of unmasked **DP-DGS**, under the DP_DGS folder, execute:

``make && ./DP-DGS``

To see the runtime of **GR19-DGS**, under the GR19 folder, execute:

``make && ./GR19-DGS``

The integration of **MAGNET** and **GR19-DGS** into the MatRiCT+ payment protocol can be found in the _matrict_plus-integration_ folder. To run the implementation, navigate to the _n10_ folder and execute:

``make && ./ringct``

By default, **MAGNET** is used as the discrete Gaussian sampler. The active sampler can be changed by modifying the corresponding #define directive in the _n10/param.h_ file.
