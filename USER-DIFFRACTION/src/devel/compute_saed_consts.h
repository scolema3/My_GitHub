/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_COMPUTE_SAED_CONSTS_H
#define LMP_COMPUTE_SAED_CONSTS_H

/*
The parameters for analytical approximation of the atomic scattering factors
used for electron diffraction are gathered from the resources: 

Colliex C et al 2004 Electron diffraction International Tables for 
Crystallography Volume C: Mathematical, Physical, and Chemical Tables 
ed E Prince (Norwell, MA: Kluwer) pp 259429

Peng L-M, Ren G, Dudarev S L and Whelan MJ 1996 Robust parameterization of 
elastic and absorptive electron atomic scattering factors 
Acta Crystallogr.A 52 25776
*/

#define SAEDmaxType 99

// list of element types associated with atomic scattering factor constants 
const char *SAEDtypeList[SAEDmaxType] = {
                "H",      "He",      "Li",      "Be",       "B",
                "C",       "N",       "O",       "F",      "Ne",
               "Na",      "Mg",      "Al",      "Si",       "P",
                "S",      "Cl",      "Ar",       "K",      "Ca",
               "Sc",      "Ti",       "V",      "Cr",      "Mn",
               "Fe",      "Co",      "Ni",      "Cu",      "Zn",
               "Ga",      "Ge",      "As",      "Se",      "Br",
               "Kr",      "Rb",      "Sr",       "Y",      "Zr",
               "Nb",      "Mo",      "Tc",      "Ru",      "Rh",
               "Pd",      "Ag",      "Cd",      "In",      "Sn",
               "Sb",      "Te",       "I",      "Xe",      "Cs",
               "Ba",      "La",      "Ce",      "Pr",      "Nd",
               "Pm",      "Sm",      "Eu",      "Gd",      "Tb",
               "Dy",      "Ho",      "Er",      "Tm",      "Yb",
               "Lu",      "Hf",      "Ta",       "W",      "Re",
               "Os",      "Ir",      "Pt",      "Au",      "Hg",
               "Tl",      "Pb",      "Bi",      "Po",      "At",
               "Rn",      "Fr",      "Ra",      "Ac",      "Th",
               "Pa",       "U",      "Np",      "Pu",      "Am",
               "Cm",      "Bk",      "Cf",      "NULL"};

// list of atomic scattering factor constants for electron diffraction
#ifdef _LMP_INTEL_OFFLOAD 
__declspec(target(mic)) const double ASFSAED[SAEDmaxType][20] = {
#else
const double ASFSAED[SAEDmaxType][20] = {
#endif
  /*  Each set of four rows in this file represents a single row in the matrix
      First two rows are constants for 0 < sin(theta)/lambda < 2
      Z-number     A1         A2         A3         A4         A5
                   B1         B2         B3         B4         B5
      Second two rows are constants for 2 < sin(theta)/lambda < 6
                   A1         A2         A3         A4         A5
                   B1         B2         B3         B4         B5*/
  /*  1*/ {    0.0349,    0.1201,    0.1970,    0.0573,    0.1195,
               0.5347,    3.5867,   12.3471,   18.9525,   38.6269,
               0.0088,    0.0449,    0.1481,    0.2356,    0.0914,
               0.1152,    1.0867,    4.9755,   16.5591,   43.2743},
  /*  2*/ {    0.0317,    0.0838,    0.1526,    0.1334,    0.0164,
               0.2507,    1.4751,    4.4938,   12.6646,   31.1653,
               0.0084,    0.0443,    0.1314,    0.1671,    0.0666,
               0.0596,    0.5360,    2.4274,    7.7852,   20.3126},
  /*  3*/ {    0.0750,    0.2249,    0.5548,    1.4954,    0.9354,
               0.3864,    2.9383,   15.3829,   53.5545,  138.7337,
               0.0478,    0.2048,    0.5253,    1.5225,    0.9853,
               0.2258,    2.1032,   12.9349,   50.7501,  136.6280},
  /*  4*/ {    0.0780,    0.2210,    0.6740,    1.3867,    0.6925,
               0.3131,    2.2381,   10.1517,   30.9061,   78.3273,
               0.0423,    0.1874,    0.6019,    1.4311,    0.7891,
               0.1445,    1.4180,    8.1165,   27.9705,   74.8684},
  /*  5*/ {    0.0909,    0.2551,    0.7738,    1.2136,    0.4606,
               0.2995,    2.1155,    8.3816,   24.1292,   63.1314,
               0.0436,    0.1898,    0.6788,    1.3273,    0.5544,
               0.1207,    1.1595,    6.2474,   21.0460,   59.3619},
  /*  6*/ {    0.0893,    0.2563,    0.7570,    1.0487,    0.3575,
               0.2465,    1.7100,    6.4094,   18.6113,   50.2523,
               0.0489,    0.2091,    0.7537,    1.1420,    0.3555,
               0.1140,    1.0825,    5.4281,   17.8811,   51.1341},
  /*  7*/ {    0.1022,    0.3219,    0.7982,    0.8197,    0.1715,
               0.2451,    1.7481,    6.1925,   17.3894,   48.1431,
               0.0267,    0.1328,    0.5301,    1.1020,    0.4215,
               0.0541,    0.5165,    2.8207,   10.6297,   34.3764},
  /*  8*/ {    0.0974,    0.2921,    0.6910,    0.6990,    0.2039,
               0.2067,    1.3815,    4.6943,   12.7105,   32.4726,
               0.0365,    0.1729,    0.5805,    0.8814,    0.3121,
               0.0652,    0.6184,    2.9449,    9.6298,   28.2194},
  /*  9*/ {    0.1083,    0.3175,    0.6487,    0.5846,    0.1421,
               0.2057,    1.3439,    4.2788,   11.3932,   28.7881,
               0.0382,    0.1822,    0.5972,    0.7707,    0.2130,
               0.0613,    0.5753,    2.6858,    8.8214,   25.6668},
  /* 10*/ {    0.1269,    0.3535,    0.5582,    0.4674,    0.1460,
               0.2200,    1.3779,    4.0203,    9.4934,   23.1278,
               0.0380,    0.1785,    0.5494,    0.6942,    0.1918,
               0.0554,    0.5087,    2.2639,    7.3316,   21.6912},
  /* 11*/ {    0.2142,    0.6853,    0.7692,    1.6589,    1.4482,
               0.3334,    2.3446,   10.0830,   48.3037,  138.2700,
               0.1260,    0.6442,    0.8893,    1.8197,    1.2988,
               0.1684,    1.7150,    8.8386,   50.8265,  147.2073},
  /* 12*/ {    0.2314,    0.6866,    0.9677,    2.1882,    1.1339,
               0.3278,    2.2720,   10.9241,   39.2898,  101.9748,
               0.1130,    0.5575,    0.9046,    2.1580,    1.4735,
               0.1356,    1.3579,    6.9255,   32.3165,   92.1138},
  /* 13*/ {    0.2390,    0.6573,    1.2011,    2.5586,    1.2312,
               0.3138,    2.1063,   10.4163,   34.4552,   98.5344,
               0.1165,    0.5504,    1.0179,    2.6295,    1.5711,
               0.1295,    1.2619,    6.8242,   28.4577,   88.4750},
  /* 14*/ {    0.2519,    0.6372,    1.3795,    2.5082,    1.0500,
               0.3075,    2.0174,    9.6746,   29.3744,   80.4732,
               0.0567,    0.3365,    0.8104,    2.4960,    2.1186,
               0.0582,    0.6155,    3.2522,   16.7929,   57.6767},
  /* 15*/ {    0.2548,    0.6106,    1.4541,    2.3204,    0.8477,
               0.2908,    1.8740,    8.5176,   24.3434,   63.2996,
               0.1005,    0.4615,    1.0663,    2.5854,    1.2725,
               0.0977,    0.9084,    4.9654,   18.5471,   54.3648},
  /* 16*/ {    0.2497,    0.5628,    1.3899,    2.1865,    0.7715,
               0.2681,    1.6711,    7.0267,   19.5377,   50.3888,
               0.0915,    0.4312,    1.0847,    2.4671,    1.0852,
               0.0838,    0.7788,    4.3462,   15.5846,   44.6365},
  /* 17*/ {    0.2443,    0.5397,    1.3919,    2.0197,    0.6621,
               0.2468,    1.5242,    6.1537,   16.6687,   42.3086,
               0.0799,    0.3891,    1.0037,    2.3332,    1.0507,
               0.0694,    0.6443,    3.5351,   12.5058,   35.8633},
  /* 18*/ {    0.2385,    0.5017,    1.3428,    1.8899,    0.6079,
               0.2289,    1.3694,    5.2561,   14.0928,   35.5361,
               0.1044,    0.4551,    1.4232,    2.1533,    0.4459,
               0.0853,    0.7701,    4.4684,   14.5864,   41.2474},
  /* 19*/ {    0.4115,    1.4031,    2.2784,    2.6742,    2.2162,
               0.3703,    3.3874,   13.1029,   68.9592,  194.4329,
               0.2149,    0.8703,    2.4999,    2.3591,    3.0318,
               0.1660,    1.6906,    8.7447,   46.7825,  165.6923},
  /* 20*/ {    0.4054,    1.3880,    2.1602,    3.7532,    2.2063,
               0.3499,    3.0991,   11.9608,   53.9353,  142.3892,
               0.2355,    0.9916,    2.3959,    3.7252,    2.5647,
               0.1742,    1.8329,    8.8407,   47.4583,  134.9613},
  /* 21*/ {    0.3787,    1.2181,    2.0594,    3.2618,    2.3870,
               0.3133,    2.5856,    9.5813,   41.7688,  116.7282,
               0.4636,    2.0802,    2.9003,    1.4193,    2.4323,
               0.3682,    4.0312,   22.6493,   71.8200,  103.3691},
  /* 22*/ {    0.3825,    1.2598,    2.0008,    3.0617,    2.0694,
               0.3040,    2.4863,    9.2783,   39.0751,  109.4583,
               0.2123,    0.8960,    2.1765,    3.0436,    2.4439,
               0.1399,    1.4568,    6.7534,   33.1168,  101.8238},
  /* 23*/ {    0.3876,    1.2750,    1.9109,    2.8314,    1.8979,
               0.2967,    2.3780,    8.7981,   35.9528,  101.7201,
               0.2369,    1.0774,    2.1894,    3.0825,    1.7190,
               0.1505,    1.6392,    7.5691,   36.8741,  107.8517},
  /* 24*/ {    0.4046,    1.3696,    1.8941,    2.0800,    1.2196,
               0.2986,    2.3958,    9.1406,   37.4701,  113.7121,
               0.1970,    0.8228,    2.0200,    2.1717,    1.7516,
               0.1197,    1.1985,    5.4097,   25.2361,   94.4290},
  /* 25*/ {    0.3796,    1.2094,    1.7815,    2.5420,    1.5937,
               0.2699,    2.0455,    7.4726,   31.0604,   91.5622,
               0.1943,    0.8190,    1.9296,    2.4968,    2.0625,
               0.1135,    1.1313,    5.0341,   24.1798,   80.5598},
  /* 26*/ {    0.3946,    1.2725,    1.7031,    2.3140,    1.4795,
               0.2717,    2.0443,    7.6007,   29.9714,   86.2265,
               0.1929,    0.8239,    1.8689,    2.3694,    1.9060,
               0.1087,    1.0806,    4.7637,   22.8500,   76.7309},
  /* 27*/ {    0.4118,    1.3161,    1.6493,    2.1930,    1.2830,
               0.2742,    2.0372,    7.7205,   29.9680,   84.9383,
               0.2186,    0.9861,    1.8540,    2.3258,    1.4685,
               0.1182,    1.2300,    5.4177,   25.7602,   80.8542},
  /* 28*/ {    0.3860,    1.1765,    1.5451,    2.0730,    1.3814,
               0.2478,    1.7660,    6.3107,   25.2204,   74.3146,
               0.2313,    1.0657,    1.8229,    2.2609,    1.1883,
               0.1210,    1.2691,    5.6870,   27.0917,   83.0285},
  /* 29*/ {    0.4314,    1.3208,    1.5236,    1.4671,    0.8562,
               0.2694,    1.9223,    7.3474,   28.9892,   90.6246,
               0.3501,    1.6558,    1.9582,    0.2134,    1.4109,
               0.1867,    1.9917,   11.3396,   53.2619,   63.2520},
  /* 30*/ {    0.4288,    1.2646,    1.4472,    1.8294,    1.0934,
               0.2593,    1.7998,    6.7500,   25.5860,   73.5284,
               0.1780,    0.8096,    1.6744,    1.9499,    1.4495,
               0.0876,    0.8650,    3.8612,   18.8726,   64.7016},
  /* 31*/ {    0.4818,    1.4032,    1.6561,    2.4605,    1.1054,
               0.2825,    1.9785,    8.7546,   32.5238,   98.5523,
               0.2135,    0.9768,    1.6669,    2.5662,    1.6790,
               0.1020,    1.0219,    4.6275,   22.8742,   80.1535},
  /* 32*/ {    0.4655,    1.3014,    1.6088,    2.6998,    1.3003,
               0.2647,    1.7926,    7.6071,   26.5541,   77.5238,
               0.2135,    0.9761,    1.6555,    2.8938,    1.6356,
               0.0989,    0.9845,    4.5527,   21.5563,   70.3903},
  /* 33*/ {    0.4517,    1.2229,    1.5852,    2.7958,    1.2638,
               0.2493,    1.6436,    6.8154,   22.3681,   62.0390,
               0.2059,    0.9518,    1.6372,    3.0490,    1.4756,
               0.0926,    0.9182,    4.3291,   19.2996,   58.9329},
  /* 34*/ {    0.4477,    1.1678,    1.5843,    2.8087,    1.1956,
               0.2405,    1.5442,    6.3231,   19.4610,   52.0233,
               0.1574,    0.7614,    1.4834,    3.0016,    1.7978,
               0.0686,    0.6808,    3.1163,   14.3458,   44.0455},
  /* 35*/ {    0.4798,    1.1948,    1.8695,    2.6953,    0.8203,
               0.2504,    1.5963,    6.9653,   19.8492,   50.3233,
               0.1899,    0.8983,    1.6358,    3.1845,    1.1518,
               0.0810,    0.7957,    3.9054,   15.7701,   45.6124},
  /* 36*/ {    0.4546,    1.0993,    1.7696,    2.7068,    0.8672,
               0.2309,    1.4279,    5.9449,   16.6752,   42.2243,
               0.1742,    0.8447,    1.5944,    3.1507,    1.1338,
               0.0723,    0.7123,    3.5192,   13.7724,   39.1148},
  /* 37*/ {    1.0160,    2.8528,    3.5466,   -7.7804,   12.1148,
               0.4853,    5.0925,   25.7851,  130.4510,  138.6775,
               0.3781,    1.4904,    3.5753,    3.0031,    3.3272,
               0.1557,    1.5347,    9.9947,   51.4251,  185.9828},
  /* 38*/ {    0.6703,    1.4926,    3.3368,    4.4600,    3.1501,
               0.3190,    2.2287,   10.3504,   52.3291,  151.2216,
               0.3723,    1.4598,    3.5124,    4.4612,    3.3031,
               0.1480,    1.4643,    9.2320,   49.8807,  148.0937},
  /* 39*/ {    0.6894,    1.5474,    3.2450,    4.2126,    2.9764,
               0.3189,    2.2904,   10.0062,   44.0771,  125.0120,
               0.3234,    1.2737,    3.2115,    4.0563,    3.7962,
               0.1244,    1.1948,    7.2756,   34.1430,  111.2079},
  /* 40*/ {    0.6719,    1.4684,    3.1668,    3.9557,    2.8920,
               0.3036,    2.1249,    8.9236,   36.8458,  108.2049,
               0.2997,    1.1879,    3.1075,    3.9740,    3.5769,
               0.1121,    1.0638,    6.3891,   28.7081,   97.4289},
  /* 41*/ {    0.6123,    1.2677,    3.0348,    3.3841,    2.3683,
               0.2709,    1.7683,    7.2489,   27.9465,   98.5624,
               0.1680,    0.9370,    2.7300,    3.8150,    3.0053,
               0.0597,    0.6524,    4.4317,   19.5540,   85.5011},
  /* 42*/ {    0.6773,    1.4798,    3.1788,    3.0824,    1.8384,
               0.2920,    2.0606,    8.1129,   30.5336,  100.0658,
               0.3069,    1.1714,    3.2293,    3.4254,    2.1224,
               0.1101,    1.0222,    5.9613,   25.1965,   93.5831},
  /* 43*/ {    0.7082,    1.6392,    3.1993,    3.4327,    1.8711,
               0.2976,    2.2106,    8.5246,   33.1456,   96.6377,
               0.2928,    1.1267,    3.1675,    3.6619,    2.5942,
               0.1020,    0.9481,    5.4713,   23.8153,   82.8991},
  /* 44*/ {    0.6735,    1.4934,    3.0966,    2.7254,    1.5597,
               0.2773,    1.9716,    7.3249,   26.6891,   90.5581,
               0.2604,    1.0442,    3.0761,    3.2175,    1.9448,
               0.0887,    0.8240,    4.8278,   19.8977,   80.4566},
  /* 45*/ {    0.6413,    1.3690,    2.9854,    2.6952,    1.5433,
               0.2580,    1.7721,    6.3854,   23.2549,   85.1517,
               0.2713,    1.0556,    3.1416,    3.0451,    1.7179,
               0.0907,    0.8324,    4.7702,   19.7862,   80.2540},
  /* 46*/ {    0.5904,    1.1775,    2.6519,    2.2875,    0.8689,
               0.2324,    1.5019,    5.1591,   15.5428,   46.8213,
               0.2003,    0.8779,    2.6135,    2.8594,    1.0258,
               0.0659,    0.6111,    3.5563,   12.7638,   44.4283},
  /* 47*/ {    0.6377,    1.3790,    2.8294,    2.3631,    1.4553,
               0.2466,    1.6974,    5.7656,   20.0943,   76.7372,
               0.2739,    1.0503,    3.1564,    2.7543,    1.4328,
               0.0881,    0.8028,    4.4451,   18.7011,   79.2633},
  /* 48*/ {    0.6364,    1.4247,    2.7802,    2.5973,    1.7886,
               0.2407,    1.6823,    5.6588,   20.7219,   69.1109,
               0.3072,    1.1303,    3.2046,    2.9329,    1.6560,
               0.0966,    0.8856,    4.6273,   20.6789,   73.4723},
  /* 49*/ {    0.6768,    1.6589,    2.7740,    3.1835,    2.1326,
               0.2522,    1.8545,    6.2936,   25.1457,   84.5448,
               0.3564,    1.3011,    3.2424,    3.4839,    2.0459,
               0.1091,    1.0452,    5.0900,   24.6578,   88.0513},
  /* 50*/ {    0.7224,    1.9610,    2.7161,    3.5603,    1.8972,
               0.2651,    2.0604,    7.3011,   27.5493,   81.3349,
               0.2966,    1.1157,    3.0973,    3.8156,    2.5281,
               0.0896,    0.8268,    4.2242,   20.6900,   71.3399},
  /* 51*/ {    0.7106,    1.9247,    2.6149,    3.8322,    1.8899,
               0.2562,    1.9646,    6.8852,   24.7648,   68.9168,
               0.2725,    1.0651,    2.9940,    4.0697,    2.5682,
               0.0809,    0.7488,    3.8710,   18.8800,   60.6499},
  /* 52*/ {    0.6947,    1.8690,    2.5356,    4.0013,    1.8955,
               0.2459,    1.8542,    6.4411,   22.1730,   59.2206,
               0.2422,    0.9692,    2.8114,    4.1509,    2.8161,
               0.0708,    0.6472,    3.3609,   16.0752,   50.1724},
  /* 53*/ {    0.7047,    1.9484,    2.5940,    4.1526,    1.5057,
               0.2455,    1.8638,    6.7639,   21.8007,   56.4395,
               0.2617,    1.0325,    2.8097,    4.4809,    2.3190,
               0.0749,    0.6914,    3.4634,   16.3603,   48.2522},
  /* 54*/ {    0.6737,    1.7908,    2.4129,    4.2100,    1.7058,
               0.2305,    1.6890,    5.8218,   18.3928,   47.2496,
               0.2334,    0.9496,    2.6381,    4.4680,    2.5020,
               0.0655,    0.6050,    3.0389,   14.0809,   41.0005},
  /* 55*/ {    1.2704,    3.8018,    5.6618,    0.9205,    4.8105,
               0.4356,    4.2058,   23.4342,  136.7780,  171.7561,
               0.5713,    2.4866,    4.9795,    4.0198,    4.4403,
               0.1626,    1.8213,   11.1049,   49.0568,  202.9987},
  /* 56*/ {    0.9049,    2.6076,    4.8498,    5.1603,    4.7388,
               0.3066,    2.4363,   12.1821,   54.6135,  161.9978,
               0.5229,    2.2874,    4.7243,    5.0807,    5.6389,
               0.1434,    1.6019,    9.4511,   42.7685,  148.4969},
  /* 57*/ {    0.8405,    2.3863,    4.6139,    5.1514,    4.7949,
               0.2791,    2.1410,   10.3400,   41.9148,  132.0204,
               0.5461,    2.3856,    5.0653,    5.7601,    4.0463,
               0.1479,    1.6552,   10.0059,   47.3245,  145.8464},
  /* 58*/ {    0.8551,    2.3915,    4.5772,    5.0278,    4.5118,
               0.2805,    2.1200,   10.1808,   42.0633,  130.9893,
               0.2227,    1.0760,    2.9482,    5.8496,    7.1834,
               0.0571,    0.5946,    3.2022,   16.4253,   95.7030},
  /* 59*/ {    0.9096,    2.5313,    4.5266,    4.6376,    4.3690,
               0.2939,    2.2471,   10.8266,   48.8842,  147.6020,
               0.5237,    2.2913,    4.6161,    4.7233,    4.8173,
               0.1360,    1.5068,    8.8213,   41.9536,  141.2424},
  /* 60*/ {    0.8807,    2.4183,    4.4448,    4.6858,    4.1725,
               0.2802,    2.0836,   10.0357,   47.4506,  146.9976,
               0.5368,    2.3301,    4.6058,    4.6621,    4.4622,
               0.1378,    1.5140,    8.8719,   43.5967,  141.8065},
  /* 61*/ {    0.9471,    2.5463,    4.3523,    4.4789,    3.9080,
               0.2977,    2.2276,   10.5762,   49.3619,  145.3580,
               0.5232,    2.2627,    4.4552,    4.4787,    4.5073,
               0.1317,    1.4336,    8.3087,   40.6010,  135.9196},
  /* 62*/ {    0.9699,    2.5837,    4.2778,    4.4575,    3.5985,
               0.3003,    2.2447,   10.6487,   50.7994,  146.4179,
               0.5162,    2.2302,    4.3449,    4.3598,    4.4292,
               0.1279,    1.3811,    7.9629,   39.1213,  132.7846},
  /* 63*/ {    0.8694,    2.2413,    3.9196,    3.9694,    4.5498,
               0.2653,    1.8590,    8.3998,   36.7397,  125.7089,
               0.5272,    2.2844,    4.3361,    4.3178,    4.0908,
               0.1285,    1.3943,    8.1081,   40.9631,  134.1233},
  /* 64*/ {    0.9673,    2.4702,    4.1148,    4.4972,    3.2099,
               0.2909,    2.1014,    9.7067,   43.4270,  125.9474,
               0.9664,    3.4052,    5.0803,    1.4991,    4.2528,
               0.2641,    2.6586,   16.2213,   80.2060,   92.5359},
  /* 65*/ {    0.9325,    2.3673,    3.8791,    3.9674,    3.7996,
               0.2761,    1.9511,    8.9296,   41.5937,  131.0122,
               0.5110,    2.1570,    4.0308,    3.9936,    4.2466,
               0.1210,    1.2704,    7.1368,   35.0354,  123.5062},
  /* 66*/ {    0.9505,    2.3705,    3.8218,    4.0471,    3.4451,
               0.2773,    1.9469,    8.8862,   43.0938,  133.1396,
               0.4974,    2.1097,    3.8906,    3.8100,    4.3084,
               0.1157,    1.2108,    6.7377,   32.4150,  116.9225},
  /* 67*/ {    0.9248,    2.2428,    3.6182,    3.7910,    3.7912,
               0.2660,    1.8183,    7.9655,   33.1129,  101.8139,
               0.4679,    1.9693,    3.7191,    3.9632,    4.2432,
               0.1069,    1.0994,    5.9769,   27.1491,   96.3119},
  /* 68*/ {    1.0373,    2.4824,    3.6558,    3.8925,    3.0056,
               0.2944,    2.0797,    9.4156,   45.8056,  132.7720,
               0.5034,    2.1088,    3.8232,    3.7299,    3.8963,
               0.1141,    1.1769,    6.6087,   33.4332,  116.4913},
  /* 69*/ {    1.0075,    2.3787,    3.5440,    3.6932,    3.1759,
               0.2816,    1.9486,    8.7162,   41.8420,  125.0320,
               0.4839,    2.0262,    3.6851,    3.5874,    4.0037,
               0.1081,    1.1012,    6.1114,   30.3728,  110.5988},
  /* 70*/ {    1.0347,    2.3911,    3.4619,    3.6556,    3.0052,
               0.2855,    1.9679,    8.7619,   42.3304,  125.6499,
               0.5221,    2.1695,    3.7567,    3.6685,    3.4274,
               0.1148,    1.1860,    6.7520,   35.6807,  118.0692},
  /* 71*/ {    0.9927,    2.2436,    3.3554,    3.7813,    3.0994,
               0.2701,    1.8073,    7.8112,   34.4849,  103.3526,
               0.4680,    1.9466,    3.5428,    3.8490,    3.6594,
               0.1015,    1.0195,    5.6058,   27.4899,   95.2846},
  /* 72*/ {    1.0295,    2.2911,    3.4110,    3.9497,    2.4925,
               0.2761,    1.8625,    8.0961,   34.2712,   98.5295,
               0.4048,    1.7370,    3.3399,    3.9448,    3.7293,
               0.0868,    0.8585,    4.6378,   21.6900,   80.2408},
  /* 73*/ {    1.0190,    2.2291,    3.4097,    3.9252,    2.2679,
               0.2694,    1.7962,    7.6944,   31.0942,   91.1089,
               0.3835,    1.6747,    3.2986,    4.0462,    3.4303,
               0.0810,    0.8020,    4.3545,   19.9644,   73.6337},
  /* 74*/ {    0.9853,    2.1167,    3.3570,    3.7981,    2.2798,
               0.2569,    1.6745,    7.0098,   26.9234,   81.3910,
               0.3661,    1.6191,    3.2455,    4.0856,    3.2064,
               0.0761,    0.7543,    4.0952,   18.2886,   68.0967},
  /* 75*/ {    0.9914,    2.0858,    3.4531,    3.8812,    1.8526,
               0.2548,    1.6518,    6.8845,   26.7234,   81.7215,
               0.3933,    1.6973,    3.4202,    4.1274,    2.6158,
               0.0806,    0.7972,    4.4237,   19.5692,   68.7477},
  /* 76*/ {    0.9813,    2.0322,    3.3665,    3.6235,    1.9741,
               0.2487,    1.5973,    6.4737,   23.2817,   70.9254,
               0.3854,    1.6555,    3.4129,    4.1111,    2.4106,
               0.0787,    0.7638,    4.2441,   18.3700,   65.1071},
  /* 77*/ {    1.0194,    2.0645,    3.4425,    3.4914,    1.6976,
               0.2554,    1.6475,    6.5966,   23.2269,   70.0272,
               0.3510,    1.5620,    3.2946,    4.0615,    2.4382,
               0.0706,    0.6904,    3.8266,   16.0812,   58.7638},
  /* 78*/ {    0.9148,    1.8096,    3.2134,    3.2953,    1.5754,
               0.2263,    1.3813,    5.3243,   17.5987,   60.0171,
               0.3083,    1.4158,    2.9662,    3.9349,    2.1709,
               0.0609,    0.5993,    3.1921,   12.5285,   49.7675},
  /* 79*/ {    0.9674,    1.8916,    3.3993,    3.0524,    1.2607,
               0.2358,    1.4712,    5.6758,   18.7119,   61.5286,
               0.3055,    1.3945,    2.9617,    3.8990,    2.0026,
               0.0596,    0.5827,    3.1035,   11.9693,   47.9106},
  /* 80*/ {    1.0033,    1.9469,    3.4396,    3.1548,    1.4180,
               0.2413,    1.5298,    5.8009,   19.4520,   60.5753,
               0.3593,    1.5736,    3.5237,    3.8109,    1.6953,
               0.0694,    0.6758,    3.8457,   15.6203,   56.6614},
  /* 81*/ {    1.0689,    2.1038,    3.6039,    3.4927,    1.8283,
               0.2540,    1.6715,    6.3509,   23.1531,   78.7099,
               0.3511,    1.5489,    3.5676,    4.0900,    2.5251,
               0.0672,    0.6522,    3.7420,   15.9791,   65.1354},
  /* 82*/ {    1.0891,    2.1867,    3.6160,    3.8031,    1.8994,
               0.2552,    1.7174,    6.5131,   23.9170,   74.7039,
               0.3540,    1.5453,    3.5975,    4.3152,    2.7743,
               0.0668,    0.6465,    3.6968,   16.2056,   61.4909},
  /* 83*/ {    1.1007,    2.2306,    3.5689,    4.1549,    2.0382,
               0.2546,    1.7351,    6.4948,   23.6464,   70.3780,
               0.3530,    1.5258,    3.5815,    4.5532,    3.0714,
               0.0661,    0.6324,    3.5906,   15.9962,   57.5760},
  /* 84*/ {    1.1568,    2.4353,    3.6459,    4.4064,    1.7179,
               0.2648,    1.8786,    7.1749,   25.1766,   69.2821,
               0.3673,    1.5772,    3.7079,    4.8582,    2.8440,
               0.0678,    0.6527,    3.7396,   17.0668,   55.9789},
  /* 85*/ {    1.0909,    2.1976,    3.3831,    4.6700,    2.1277,
               0.2466,    1.6707,    6.0197,   20.7657,   57.2663,
               0.3547,    1.5206,    3.5621,    5.0184,    3.0075,
               0.0649,    0.6188,    3.4696,   15.6090,   49.4818},
  /* 86*/ {    1.0756,    2.1630,    3.3178,    4.8852,    2.0489,
               0.2402,    1.6169,    5.7644,   19.4568,   52.5009,
               0.4586,    1.7781,    3.9877,    5.7273,    1.5460,
               0.0831,    0.7840,    4.3599,   20.0128,   62.1535},
  /* 87*/ {    1.4282,    3.5081,    5.6767,    4.1964,    3.8946,
               0.3183,    2.6889,   13.4816,   54.3866,  200.8321,
               0.8282,    2.9941,    5.6597,    4.9292,    4.2889,
               0.1515,    1.6163,    9.7752,   42.8480,  190.7366},
  /* 88*/ {    1.3127,    3.1243,    5.2988,    5.3891,    5.4133,
               0.2887,    2.2897,   10.8276,   43.5389,  145.6109,
               1.4129,    4.4269,    7.0460,   -1.0573,    8.6430,
               0.2921,    3.1381,   19.6767,  102.0430,  113.9798},
  /* 89*/ {    1.3128,    3.1021,    5.3385,    5.9611,    4.7562,
               0.2861,    2.2509,   10.5287,   41.7796,  128.2973,
               0.7169,    2.5710,    5.1791,    6.3484,    5.6474,
               0.1263,    1.2900,    7.3686,   32.4490,  118.0558},
  /* 90*/ {    1.2553,    2.9178,    5.0862,    6.1206,    4.7122,
               0.2701,    2.0636,    9.3051,   34.5977,  107.9200,
               0.6958,    2.4936,    5.1269,    6.6988,    5.0799,
               0.1211,    1.2247,    6.9398,   30.0991,  105.1960},
  /* 91*/ {    1.3218,    3.1444,    5.4371,    5.6444,    4.0107,
               0.2827,    2.2250,   10.2454,   41.1162,  124.4449,
               1.2502,    4.2284,    7.0489,    1.1390,    5.8222,
               0.2415,    2.6442,   16.3313,   73.5757,   91.9401},
  /* 92*/ {    1.3382,    3.2043,    5.4558,    5.4839,    3.6342,
               0.2838,    2.2452,   10.2519,   41.7251,  124.9023,
               0.6410,    2.2643,    4.8713,    5.9287,    5.3935,
               0.1097,    1.0644,    5.7907,   25.0261,  101.3899},
  /* 93*/ {    1.5193,    4.0053,    6.5327,   -0.1402,    6.7489,
               0.3213,    2.8206,   14.8878,   68.9103,   81.7257,
               0.6938,    2.4652,    5.1227,    5.5965,    4.8543,
               0.1171,    1.1757,    6.4053,   27.5217,  103.0482},
  /* 94*/ {    1.3517,    3.2937,    5.3213,    4.6466,    3.5714,
               0.2813,    2.2418,    9.9952,   42.7939,  132.1739,
               0.6902,    2.4509,    5.1284,    5.0339,    4.8575,
               0.1153,    1.1545,    6.2291,   27.0741,  111.3150},
  /* 95*/ {    1.2135,    2.7962,    4.7545,    4.5731,    4.4786,
               0.2483,    1.8437,    7.5421,   29.3841,  112.4579,
               0.7577,    2.7264,    5.4184,    4.8198,    4.1013,
               0.1257,    1.3044,    7.1035,   32.4649,  118.8647},
  /* 96*/ {    1.2937,    3.1100,    5.0393,    4.7546,    3.5031,
               0.2638,    2.0341,    8.7101,   35.2992,  109.4972,
               0.7567,    2.7565,    5.4364,    5.1918,    3.5643,
               0.1239,    1.2979,    7.0798,   32.7871,  110.1512},
  /* 97*/ {    1.2915,    3.1023,    4.9309,    4.6009,    3.4661,
               0.2611,    2.0023,    8.4377,   34.1559,  105.8911,
               0.7492,    2.7267,    5.3521,    5.0369,    3.5321,
               0.1217,    1.2651,    6.8101,   31.6088,  106.4853},
  /* 98*/ {    1.2089,    2.7391,    4.3482,    4.0047,    4.6497,
               0.2421,    1.7487,    6.7262,   23.2153,   80.3108,
               0.8100,    3.0001,    5.4635,    4.1756,    3.5066,
               0.1310,    1.4038,    7.6057,   34.0186,   90.5226},
  /* 99*/ {    0.2000,    0.2000,    0.2000,    0.2000,    0.2000,
               0.0000,    0.0000,    0.0000,    0.0000,    0.0000,
               0.2000,    0.2000,    0.2000,    0.2000,    0.2000,
               0.0000,    0.0000,    0.0000,    0.0000,    0.0000},               
  };

#endif
