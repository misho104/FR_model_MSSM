# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 13.0.0 for Mac OS X x86 (64-bit) (December 3, 2021)
# Date: Sat 12 Feb 2022 18:33:51



from object_library import all_parameters, Parameter


from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot

# This is a default parameter object representing 0.
ZERO = Parameter(name = 'ZERO',
                 nature = 'internal',
                 type = 'real',
                 value = '0.0',
                 texname = '0')

# User-defined parameters.
RRd1x1 = Parameter(name = 'RRd1x1',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRd1x1}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 1, 1 ])

RRd2x2 = Parameter(name = 'RRd2x2',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRd2x2}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 2, 2 ])

RRd3x3 = Parameter(name = 'RRd3x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.938737896,
                   texname = '\\text{RRd3x3}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 3, 3 ])

RRd3x6 = Parameter(name = 'RRd3x6',
                   nature = 'external',
                   type = 'real',
                   value = 0.344631925,
                   texname = '\\text{RRd3x6}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 3, 6 ])

RRd4x4 = Parameter(name = 'RRd4x4',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRd4x4}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 4, 4 ])

RRd5x5 = Parameter(name = 'RRd5x5',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRd5x5}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 5, 5 ])

RRd6x3 = Parameter(name = 'RRd6x3',
                   nature = 'external',
                   type = 'real',
                   value = -0.344631925,
                   texname = '\\text{RRd6x3}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 6, 3 ])

RRd6x6 = Parameter(name = 'RRd6x6',
                   nature = 'external',
                   type = 'real',
                   value = 0.938737896,
                   texname = '\\text{RRd6x6}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 6, 6 ])

alp = Parameter(name = 'alp',
                nature = 'external',
                type = 'real',
                value = -0.11382521,
                texname = '\\alpha',
                lhablock = 'FRALPHA',
                lhacode = [ 1 ])

RMUH = Parameter(name = 'RMUH',
                 nature = 'external',
                 type = 'real',
                 value = 357.680977,
                 texname = '\\text{RMUH}',
                 lhablock = 'HMIX',
                 lhacode = [ 1 ])

tb = Parameter(name = 'tb',
               nature = 'external',
               type = 'real',
               value = 9.74862403,
               texname = 't_b',
               lhablock = 'HMIX',
               lhacode = [ 2 ])

MA2 = Parameter(name = 'MA2',
                nature = 'external',
                type = 'real',
                value = 166439.065,
                texname = 'm_A^2',
                lhablock = 'HMIX',
                lhacode = [ 4 ])

INN1x1 = Parameter(name = 'INN1x1',
                   nature = 'external',
                   type = 'real',
                   value = 1.e-8,
                   texname = '\\text{INN1x1}',
                   lhablock = 'IMNMIX',
                   lhacode = [ 1, 1 ])

INN1x2 = Parameter(name = 'INN1x2',
                   nature = 'external',
                   type = 'real',
                   value = 1.e-8,
                   texname = '\\text{INN1x2}',
                   lhablock = 'IMNMIX',
                   lhacode = [ 1, 2 ])

INN1x3 = Parameter(name = 'INN1x3',
                   nature = 'external',
                   type = 'real',
                   value = 1.e-8,
                   texname = '\\text{INN1x3}',
                   lhablock = 'IMNMIX',
                   lhacode = [ 1, 3 ])

INN1x4 = Parameter(name = 'INN1x4',
                   nature = 'external',
                   type = 'real',
                   value = 1.e-8,
                   texname = '\\text{INN1x4}',
                   lhablock = 'IMNMIX',
                   lhacode = [ 1, 4 ])

INN2x1 = Parameter(name = 'INN2x1',
                   nature = 'external',
                   type = 'real',
                   value = 1.e-8,
                   texname = '\\text{INN2x1}',
                   lhablock = 'IMNMIX',
                   lhacode = [ 2, 1 ])

INN2x2 = Parameter(name = 'INN2x2',
                   nature = 'external',
                   type = 'real',
                   value = 1.e-8,
                   texname = '\\text{INN2x2}',
                   lhablock = 'IMNMIX',
                   lhacode = [ 2, 2 ])

INN2x3 = Parameter(name = 'INN2x3',
                   nature = 'external',
                   type = 'real',
                   value = 1.e-8,
                   texname = '\\text{INN2x3}',
                   lhablock = 'IMNMIX',
                   lhacode = [ 2, 3 ])

INN2x4 = Parameter(name = 'INN2x4',
                   nature = 'external',
                   type = 'real',
                   value = 1.e-8,
                   texname = '\\text{INN2x4}',
                   lhablock = 'IMNMIX',
                   lhacode = [ 2, 4 ])

INN3x1 = Parameter(name = 'INN3x1',
                   nature = 'external',
                   type = 'real',
                   value = 1.e-8,
                   texname = '\\text{INN3x1}',
                   lhablock = 'IMNMIX',
                   lhacode = [ 3, 1 ])

INN3x2 = Parameter(name = 'INN3x2',
                   nature = 'external',
                   type = 'real',
                   value = 1.e-8,
                   texname = '\\text{INN3x2}',
                   lhablock = 'IMNMIX',
                   lhacode = [ 3, 2 ])

INN3x3 = Parameter(name = 'INN3x3',
                   nature = 'external',
                   type = 'real',
                   value = 1.e-8,
                   texname = '\\text{INN3x3}',
                   lhablock = 'IMNMIX',
                   lhacode = [ 3, 3 ])

INN3x4 = Parameter(name = 'INN3x4',
                   nature = 'external',
                   type = 'real',
                   value = 1.e-8,
                   texname = '\\text{INN3x4}',
                   lhablock = 'IMNMIX',
                   lhacode = [ 3, 4 ])

INN4x1 = Parameter(name = 'INN4x1',
                   nature = 'external',
                   type = 'real',
                   value = 1.e-8,
                   texname = '\\text{INN4x1}',
                   lhablock = 'IMNMIX',
                   lhacode = [ 4, 1 ])

INN4x2 = Parameter(name = 'INN4x2',
                   nature = 'external',
                   type = 'real',
                   value = 1.e-8,
                   texname = '\\text{INN4x2}',
                   lhablock = 'IMNMIX',
                   lhacode = [ 4, 2 ])

INN4x3 = Parameter(name = 'INN4x3',
                   nature = 'external',
                   type = 'real',
                   value = 1.e-8,
                   texname = '\\text{INN4x3}',
                   lhablock = 'IMNMIX',
                   lhacode = [ 4, 3 ])

INN4x4 = Parameter(name = 'INN4x4',
                   nature = 'external',
                   type = 'real',
                   value = 1.e-8,
                   texname = '\\text{INN4x4}',
                   lhablock = 'IMNMIX',
                   lhacode = [ 4, 4 ])

RmD21x1 = Parameter(name = 'RmD21x1',
                    nature = 'external',
                    type = 'real',
                    value = 273684.674,
                    texname = '\\text{RmD21x1}',
                    lhablock = 'MSD2',
                    lhacode = [ 1, 1 ])

RmD22x2 = Parameter(name = 'RmD22x2',
                    nature = 'external',
                    type = 'real',
                    value = 273684.674,
                    texname = '\\text{RmD22x2}',
                    lhablock = 'MSD2',
                    lhacode = [ 2, 2 ])

RmD23x3 = Parameter(name = 'RmD23x3',
                    nature = 'external',
                    type = 'real',
                    value = 270261.969,
                    texname = '\\text{RmD23x3}',
                    lhablock = 'MSD2',
                    lhacode = [ 3, 3 ])

RmE21x1 = Parameter(name = 'RmE21x1',
                    nature = 'external',
                    type = 'real',
                    value = 18630.6287,
                    texname = '\\text{RmE21x1}',
                    lhablock = 'MSE2',
                    lhacode = [ 1, 1 ])

RmE22x2 = Parameter(name = 'RmE22x2',
                    nature = 'external',
                    type = 'real',
                    value = 18630.6287,
                    texname = '\\text{RmE22x2}',
                    lhablock = 'MSE2',
                    lhacode = [ 2, 2 ])

RmE23x3 = Parameter(name = 'RmE23x3',
                    nature = 'external',
                    type = 'real',
                    value = 17967.6406,
                    texname = '\\text{RmE23x3}',
                    lhablock = 'MSE2',
                    lhacode = [ 3, 3 ])

RmL21x1 = Parameter(name = 'RmL21x1',
                    nature = 'external',
                    type = 'real',
                    value = 38155.67,
                    texname = '\\text{RmL21x1}',
                    lhablock = 'MSL2',
                    lhacode = [ 1, 1 ])

RmL22x2 = Parameter(name = 'RmL22x2',
                    nature = 'external',
                    type = 'real',
                    value = 38155.67,
                    texname = '\\text{RmL22x2}',
                    lhablock = 'MSL2',
                    lhacode = [ 2, 2 ])

RmL23x3 = Parameter(name = 'RmL23x3',
                    nature = 'external',
                    type = 'real',
                    value = 37828.6769,
                    texname = '\\text{RmL23x3}',
                    lhablock = 'MSL2',
                    lhacode = [ 3, 3 ])

RMx1 = Parameter(name = 'RMx1',
                 nature = 'external',
                 type = 'real',
                 value = 101.396534,
                 texname = '\\text{RMx1}',
                 lhablock = 'MSOFT',
                 lhacode = [ 1 ])

RMx2 = Parameter(name = 'RMx2',
                 nature = 'external',
                 type = 'real',
                 value = 191.504241,
                 texname = '\\text{RMx2}',
                 lhablock = 'MSOFT',
                 lhacode = [ 2 ])

RMx3 = Parameter(name = 'RMx3',
                 nature = 'external',
                 type = 'real',
                 value = 588.263031,
                 texname = '\\text{RMx3}',
                 lhablock = 'MSOFT',
                 lhacode = [ 3 ])

mHd2 = Parameter(name = 'mHd2',
                 nature = 'external',
                 type = 'real',
                 value = 32337.4943,
                 texname = 'm_{H_d}^2',
                 lhablock = 'MSOFT',
                 lhacode = [ 21 ])

mHu2 = Parameter(name = 'mHu2',
                 nature = 'external',
                 type = 'real',
                 value = -128800.134,
                 texname = 'm_{H_u}^2',
                 lhablock = 'MSOFT',
                 lhacode = [ 22 ])

RmQ21x1 = Parameter(name = 'RmQ21x1',
                    nature = 'external',
                    type = 'real',
                    value = 299836.701,
                    texname = '\\text{RmQ21x1}',
                    lhablock = 'MSQ2',
                    lhacode = [ 1, 1 ])

RmQ22x2 = Parameter(name = 'RmQ22x2',
                    nature = 'external',
                    type = 'real',
                    value = 299836.701,
                    texname = '\\text{RmQ22x2}',
                    lhablock = 'MSQ2',
                    lhacode = [ 2, 2 ])

RmQ23x3 = Parameter(name = 'RmQ23x3',
                    nature = 'external',
                    type = 'real',
                    value = 248765.367,
                    texname = '\\text{RmQ23x3}',
                    lhablock = 'MSQ2',
                    lhacode = [ 3, 3 ])

RmU21x1 = Parameter(name = 'RmU21x1',
                    nature = 'external',
                    type = 'real',
                    value = 280382.106,
                    texname = '\\text{RmU21x1}',
                    lhablock = 'MSU2',
                    lhacode = [ 1, 1 ])

RmU22x2 = Parameter(name = 'RmU22x2',
                    nature = 'external',
                    type = 'real',
                    value = 280382.106,
                    texname = '\\text{RmU22x2}',
                    lhablock = 'MSU2',
                    lhacode = [ 2, 2 ])

RmU23x3 = Parameter(name = 'RmU23x3',
                    nature = 'external',
                    type = 'real',
                    value = 179137.072,
                    texname = '\\text{RmU23x3}',
                    lhablock = 'MSU2',
                    lhacode = [ 3, 3 ])

RNN1x1 = Parameter(name = 'RNN1x1',
                   nature = 'external',
                   type = 'real',
                   value = 0.98636443,
                   texname = '\\text{RNN1x1}',
                   lhablock = 'NMIX',
                   lhacode = [ 1, 1 ])

RNN1x2 = Parameter(name = 'RNN1x2',
                   nature = 'external',
                   type = 'real',
                   value = -0.0531103553,
                   texname = '\\text{RNN1x2}',
                   lhablock = 'NMIX',
                   lhacode = [ 1, 2 ])

RNN1x3 = Parameter(name = 'RNN1x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.146433995,
                   texname = '\\text{RNN1x3}',
                   lhablock = 'NMIX',
                   lhacode = [ 1, 3 ])

RNN1x4 = Parameter(name = 'RNN1x4',
                   nature = 'external',
                   type = 'real',
                   value = -0.0531186117,
                   texname = '\\text{RNN1x4}',
                   lhablock = 'NMIX',
                   lhacode = [ 1, 4 ])

RNN2x1 = Parameter(name = 'RNN2x1',
                   nature = 'external',
                   type = 'real',
                   value = 0.0993505358,
                   texname = '\\text{RNN2x1}',
                   lhablock = 'NMIX',
                   lhacode = [ 2, 1 ])

RNN2x2 = Parameter(name = 'RNN2x2',
                   nature = 'external',
                   type = 'real',
                   value = 0.944949299,
                   texname = '\\text{RNN2x2}',
                   lhablock = 'NMIX',
                   lhacode = [ 2, 2 ])

RNN2x3 = Parameter(name = 'RNN2x3',
                   nature = 'external',
                   type = 'real',
                   value = -0.26984672,
                   texname = '\\text{RNN2x3}',
                   lhablock = 'NMIX',
                   lhacode = [ 2, 3 ])

RNN2x4 = Parameter(name = 'RNN2x4',
                   nature = 'external',
                   type = 'real',
                   value = 0.156150698,
                   texname = '\\text{RNN2x4}',
                   lhablock = 'NMIX',
                   lhacode = [ 2, 4 ])

RNN3x1 = Parameter(name = 'RNN3x1',
                   nature = 'external',
                   type = 'real',
                   value = -0.0603388002,
                   texname = '\\text{RNN3x1}',
                   lhablock = 'NMIX',
                   lhacode = [ 3, 1 ])

RNN3x2 = Parameter(name = 'RNN3x2',
                   nature = 'external',
                   type = 'real',
                   value = 0.0877004854,
                   texname = '\\text{RNN3x2}',
                   lhablock = 'NMIX',
                   lhacode = [ 3, 2 ])

RNN3x3 = Parameter(name = 'RNN3x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.695877493,
                   texname = '\\text{RNN3x3}',
                   lhablock = 'NMIX',
                   lhacode = [ 3, 3 ])

RNN3x4 = Parameter(name = 'RNN3x4',
                   nature = 'external',
                   type = 'real',
                   value = 0.710226984,
                   texname = '\\text{RNN3x4}',
                   lhablock = 'NMIX',
                   lhacode = [ 3, 4 ])

RNN4x1 = Parameter(name = 'RNN4x1',
                   nature = 'external',
                   type = 'real',
                   value = -0.116507132,
                   texname = '\\text{RNN4x1}',
                   lhablock = 'NMIX',
                   lhacode = [ 4, 1 ])

RNN4x2 = Parameter(name = 'RNN4x2',
                   nature = 'external',
                   type = 'real',
                   value = 0.310739017,
                   texname = '\\text{RNN4x2}',
                   lhablock = 'NMIX',
                   lhacode = [ 4, 2 ])

RNN4x3 = Parameter(name = 'RNN4x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.64922596,
                   texname = '\\text{RNN4x3}',
                   lhablock = 'NMIX',
                   lhacode = [ 4, 3 ])

RNN4x4 = Parameter(name = 'RNN4x4',
                   nature = 'external',
                   type = 'real',
                   value = -0.684377823,
                   texname = '\\text{RNN4x4}',
                   lhablock = 'NMIX',
                   lhacode = [ 4, 4 ])

RRl1x1 = Parameter(name = 'RRl1x1',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRl1x1}',
                   lhablock = 'SELMIX',
                   lhacode = [ 1, 1 ])

RRl2x2 = Parameter(name = 'RRl2x2',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRl2x2}',
                   lhablock = 'SELMIX',
                   lhacode = [ 2, 2 ])

RRl3x3 = Parameter(name = 'RRl3x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.28248719,
                   texname = '\\text{RRl3x3}',
                   lhablock = 'SELMIX',
                   lhacode = [ 3, 3 ])

RRl3x6 = Parameter(name = 'RRl3x6',
                   nature = 'external',
                   type = 'real',
                   value = 0.959271071,
                   texname = '\\text{RRl3x6}',
                   lhablock = 'SELMIX',
                   lhacode = [ 3, 6 ])

RRl4x4 = Parameter(name = 'RRl4x4',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRl4x4}',
                   lhablock = 'SELMIX',
                   lhacode = [ 4, 4 ])

RRl5x5 = Parameter(name = 'RRl5x5',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRl5x5}',
                   lhablock = 'SELMIX',
                   lhacode = [ 5, 5 ])

RRl6x3 = Parameter(name = 'RRl6x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.959271071,
                   texname = '\\text{RRl6x3}',
                   lhablock = 'SELMIX',
                   lhacode = [ 6, 3 ])

RRl6x6 = Parameter(name = 'RRl6x6',
                   nature = 'external',
                   type = 'real',
                   value = -0.28248719,
                   texname = '\\text{RRl6x6}',
                   lhablock = 'SELMIX',
                   lhacode = [ 6, 6 ])

aEWM1 = Parameter(name = 'aEWM1',
                  nature = 'external',
                  type = 'real',
                  value = 127.934,
                  texname = '\\alpha _w^{-1}',
                  lhablock = 'SMINPUTS',
                  lhacode = [ 1 ])

aS = Parameter(name = 'aS',
               nature = 'external',
               type = 'real',
               value = 0.118,
               texname = '\\alpha _s',
               lhablock = 'SMINPUTS',
               lhacode = [ 3 ])

RRn1x1 = Parameter(name = 'RRn1x1',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRn1x1}',
                   lhablock = 'SNUMIX',
                   lhacode = [ 1, 1 ])

RRn2x2 = Parameter(name = 'RRn2x2',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRn2x2}',
                   lhablock = 'SNUMIX',
                   lhacode = [ 2, 2 ])

RRn3x3 = Parameter(name = 'RRn3x3',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRn3x3}',
                   lhablock = 'SNUMIX',
                   lhacode = [ 3, 3 ])

Rtd3x3 = Parameter(name = 'Rtd3x3',
                   nature = 'external',
                   type = 'real',
                   value = -110.693742,
                   texname = '\\text{Rtd3x3}',
                   lhablock = 'TD',
                   lhacode = [ 3, 3 ])

Rte3x3 = Parameter(name = 'Rte3x3',
                   nature = 'external',
                   type = 'real',
                   value = -25.4019727,
                   texname = '\\text{Rte3x3}',
                   lhablock = 'TE',
                   lhacode = [ 3, 3 ])

Rtu3x3 = Parameter(name = 'Rtu3x3',
                   nature = 'external',
                   type = 'real',
                   value = -444.752457,
                   texname = '\\text{Rtu3x3}',
                   lhablock = 'TU',
                   lhacode = [ 3, 3 ])

RUU1x1 = Parameter(name = 'RUU1x1',
                   nature = 'external',
                   type = 'real',
                   value = 0.916834859,
                   texname = '\\text{RUU1x1}',
                   lhablock = 'UMIX',
                   lhacode = [ 1, 1 ])

RUU1x2 = Parameter(name = 'RUU1x2',
                   nature = 'external',
                   type = 'real',
                   value = -0.399266629,
                   texname = '\\text{RUU1x2}',
                   lhablock = 'UMIX',
                   lhacode = [ 1, 2 ])

RUU2x1 = Parameter(name = 'RUU2x1',
                   nature = 'external',
                   type = 'real',
                   value = 0.399266629,
                   texname = '\\text{RUU2x1}',
                   lhablock = 'UMIX',
                   lhacode = [ 2, 1 ])

RUU2x2 = Parameter(name = 'RUU2x2',
                   nature = 'external',
                   type = 'real',
                   value = 0.916834859,
                   texname = '\\text{RUU2x2}',
                   lhablock = 'UMIX',
                   lhacode = [ 2, 2 ])

RMNS1x1 = Parameter(name = 'RMNS1x1',
                    nature = 'external',
                    type = 'real',
                    value = 1.,
                    texname = '\\text{RMNS1x1}',
                    lhablock = 'UPMNS',
                    lhacode = [ 1, 1 ])

RMNS2x2 = Parameter(name = 'RMNS2x2',
                    nature = 'external',
                    type = 'real',
                    value = 1.,
                    texname = '\\text{RMNS2x2}',
                    lhablock = 'UPMNS',
                    lhacode = [ 2, 2 ])

RMNS3x3 = Parameter(name = 'RMNS3x3',
                    nature = 'external',
                    type = 'real',
                    value = 1.,
                    texname = '\\text{RMNS3x3}',
                    lhablock = 'UPMNS',
                    lhacode = [ 3, 3 ])

RRu1x1 = Parameter(name = 'RRu1x1',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRu1x1}',
                   lhablock = 'USQMIX',
                   lhacode = [ 1, 1 ])

RRu2x2 = Parameter(name = 'RRu2x2',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRu2x2}',
                   lhablock = 'USQMIX',
                   lhacode = [ 2, 2 ])

RRu3x3 = Parameter(name = 'RRu3x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.55364496,
                   texname = '\\text{RRu3x3}',
                   lhablock = 'USQMIX',
                   lhacode = [ 3, 3 ])

RRu3x6 = Parameter(name = 'RRu3x6',
                   nature = 'external',
                   type = 'real',
                   value = 0.83275282,
                   texname = '\\text{RRu3x6}',
                   lhablock = 'USQMIX',
                   lhacode = [ 3, 6 ])

RRu4x4 = Parameter(name = 'RRu4x4',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRu4x4}',
                   lhablock = 'USQMIX',
                   lhacode = [ 4, 4 ])

RRu5x5 = Parameter(name = 'RRu5x5',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRu5x5}',
                   lhablock = 'USQMIX',
                   lhacode = [ 5, 5 ])

RRu6x3 = Parameter(name = 'RRu6x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.83275282,
                   texname = '\\text{RRu6x3}',
                   lhablock = 'USQMIX',
                   lhacode = [ 6, 3 ])

RRu6x6 = Parameter(name = 'RRu6x6',
                   nature = 'external',
                   type = 'real',
                   value = -0.55364496,
                   texname = '\\text{RRu6x6}',
                   lhablock = 'USQMIX',
                   lhacode = [ 6, 6 ])

RCKM1x1 = Parameter(name = 'RCKM1x1',
                    nature = 'external',
                    type = 'real',
                    value = 1.,
                    texname = '\\text{RCKM1x1}',
                    lhablock = 'VCKM',
                    lhacode = [ 1, 1 ])

RCKM2x2 = Parameter(name = 'RCKM2x2',
                    nature = 'external',
                    type = 'real',
                    value = 1.,
                    texname = '\\text{RCKM2x2}',
                    lhablock = 'VCKM',
                    lhacode = [ 2, 2 ])

RCKM3x3 = Parameter(name = 'RCKM3x3',
                    nature = 'external',
                    type = 'real',
                    value = 1.,
                    texname = '\\text{RCKM3x3}',
                    lhablock = 'VCKM',
                    lhacode = [ 3, 3 ])

RVV1x1 = Parameter(name = 'RVV1x1',
                   nature = 'external',
                   type = 'real',
                   value = 0.972557835,
                   texname = '\\text{RVV1x1}',
                   lhablock = 'VMIX',
                   lhacode = [ 1, 1 ])

RVV1x2 = Parameter(name = 'RVV1x2',
                   nature = 'external',
                   type = 'real',
                   value = -0.232661249,
                   texname = '\\text{RVV1x2}',
                   lhablock = 'VMIX',
                   lhacode = [ 1, 2 ])

RVV2x1 = Parameter(name = 'RVV2x1',
                   nature = 'external',
                   type = 'real',
                   value = 0.232661249,
                   texname = '\\text{RVV2x1}',
                   lhablock = 'VMIX',
                   lhacode = [ 2, 1 ])

RVV2x2 = Parameter(name = 'RVV2x2',
                   nature = 'external',
                   type = 'real',
                   value = 0.972557835,
                   texname = '\\text{RVV2x2}',
                   lhablock = 'VMIX',
                   lhacode = [ 2, 2 ])

Ryd3x3 = Parameter(name = 'Ryd3x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.138840206,
                   texname = '\\text{Ryd3x3}',
                   lhablock = 'YD',
                   lhacode = [ 3, 3 ])

Rye3x3 = Parameter(name = 'Rye3x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.10089081,
                   texname = '\\text{Rye3x3}',
                   lhablock = 'YE',
                   lhacode = [ 3, 3 ])

Ryu3x3 = Parameter(name = 'Ryu3x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.89284455,
                   texname = '\\text{Ryu3x3}',
                   lhablock = 'YU',
                   lhacode = [ 3, 3 ])

MZ = Parameter(name = 'MZ',
               nature = 'external',
               type = 'real',
               value = 91.1876,
               texname = '\\text{MZ}',
               lhablock = 'MASS',
               lhacode = [ 23 ])

MW = Parameter(name = 'MW',
               nature = 'external',
               type = 'real',
               value = 79.8290131,
               texname = '\\text{MW}',
               lhablock = 'MASS',
               lhacode = [ 24 ])

Mneu1 = Parameter(name = 'Mneu1',
                  nature = 'external',
                  type = 'real',
                  value = 96.6880686,
                  texname = '\\text{Mneu1}',
                  lhablock = 'MASS',
                  lhacode = [ 1000022 ])

Mneu2 = Parameter(name = 'Mneu2',
                  nature = 'external',
                  type = 'real',
                  value = 181.088157,
                  texname = '\\text{Mneu2}',
                  lhablock = 'MASS',
                  lhacode = [ 1000023 ])

Mneu3 = Parameter(name = 'Mneu3',
                  nature = 'external',
                  type = 'real',
                  value = -363.756027,
                  texname = '\\text{Mneu3}',
                  lhablock = 'MASS',
                  lhacode = [ 1000025 ])

Mneu4 = Parameter(name = 'Mneu4',
                  nature = 'external',
                  type = 'real',
                  value = 381.729382,
                  texname = '\\text{Mneu4}',
                  lhablock = 'MASS',
                  lhacode = [ 1000035 ])

Mch1 = Parameter(name = 'Mch1',
                 nature = 'external',
                 type = 'real',
                 value = 181.696474,
                 texname = '\\text{Mch1}',
                 lhablock = 'MASS',
                 lhacode = [ 1000024 ])

Mch2 = Parameter(name = 'Mch2',
                 nature = 'external',
                 type = 'real',
                 value = 379.93932,
                 texname = '\\text{Mch2}',
                 lhablock = 'MASS',
                 lhacode = [ 1000037 ])

Mgo = Parameter(name = 'Mgo',
                nature = 'external',
                type = 'real',
                value = 607.713704,
                texname = '\\text{Mgo}',
                lhablock = 'MASS',
                lhacode = [ 1000021 ])

MH01 = Parameter(name = 'MH01',
                 nature = 'external',
                 type = 'real',
                 value = 110.899057,
                 texname = '\\text{MH01}',
                 lhablock = 'MASS',
                 lhacode = [ 25 ])

MH02 = Parameter(name = 'MH02',
                 nature = 'external',
                 type = 'real',
                 value = 399.960116,
                 texname = '\\text{MH02}',
                 lhablock = 'MASS',
                 lhacode = [ 35 ])

MA0 = Parameter(name = 'MA0',
                nature = 'external',
                type = 'real',
                value = 399.583917,
                texname = '\\text{MA0}',
                lhablock = 'MASS',
                lhacode = [ 36 ])

MH = Parameter(name = 'MH',
               nature = 'external',
               type = 'real',
               value = 407.879012,
               texname = '\\text{MH}',
               lhablock = 'MASS',
               lhacode = [ 37 ])

Mta = Parameter(name = 'Mta',
                nature = 'external',
                type = 'real',
                value = 1.777,
                texname = '\\text{Mta}',
                lhablock = 'MASS',
                lhacode = [ 15 ])

MT = Parameter(name = 'MT',
               nature = 'external',
               type = 'real',
               value = 175.,
               texname = '\\text{MT}',
               lhablock = 'MASS',
               lhacode = [ 6 ])

MB = Parameter(name = 'MB',
               nature = 'external',
               type = 'real',
               value = 4.88991651,
               texname = '\\text{MB}',
               lhablock = 'MASS',
               lhacode = [ 5 ])

Msn1 = Parameter(name = 'Msn1',
                 nature = 'external',
                 type = 'real',
                 value = 185.258326,
                 texname = '\\text{Msn1}',
                 lhablock = 'MASS',
                 lhacode = [ 1000012 ])

Msn2 = Parameter(name = 'Msn2',
                 nature = 'external',
                 type = 'real',
                 value = 185.258326,
                 texname = '\\text{Msn2}',
                 lhablock = 'MASS',
                 lhacode = [ 1000014 ])

Msn3 = Parameter(name = 'Msn3',
                 nature = 'external',
                 type = 'real',
                 value = 184.708464,
                 texname = '\\text{Msn3}',
                 lhablock = 'MASS',
                 lhacode = [ 1000016 ])

Msl1 = Parameter(name = 'Msl1',
                 nature = 'external',
                 type = 'real',
                 value = 202.91569,
                 texname = '\\text{Msl1}',
                 lhablock = 'MASS',
                 lhacode = [ 1000011 ])

Msl2 = Parameter(name = 'Msl2',
                 nature = 'external',
                 type = 'real',
                 value = 202.91569,
                 texname = '\\text{Msl2}',
                 lhablock = 'MASS',
                 lhacode = [ 1000013 ])

Msl3 = Parameter(name = 'Msl3',
                 nature = 'external',
                 type = 'real',
                 value = 134.490864,
                 texname = '\\text{Msl3}',
                 lhablock = 'MASS',
                 lhacode = [ 1000015 ])

Msl4 = Parameter(name = 'Msl4',
                 nature = 'external',
                 type = 'real',
                 value = 144.102799,
                 texname = '\\text{Msl4}',
                 lhablock = 'MASS',
                 lhacode = [ 2000011 ])

Msl5 = Parameter(name = 'Msl5',
                 nature = 'external',
                 type = 'real',
                 value = 144.102799,
                 texname = '\\text{Msl5}',
                 lhablock = 'MASS',
                 lhacode = [ 2000013 ])

Msl6 = Parameter(name = 'Msl6',
                 nature = 'external',
                 type = 'real',
                 value = 206.867805,
                 texname = '\\text{Msl6}',
                 lhablock = 'MASS',
                 lhacode = [ 2000015 ])

Msu1 = Parameter(name = 'Msu1',
                 nature = 'external',
                 type = 'real',
                 value = 561.119014,
                 texname = '\\text{Msu1}',
                 lhablock = 'MASS',
                 lhacode = [ 1000002 ])

Msu2 = Parameter(name = 'Msu2',
                 nature = 'external',
                 type = 'real',
                 value = 561.119014,
                 texname = '\\text{Msu2}',
                 lhablock = 'MASS',
                 lhacode = [ 1000004 ])

Msu3 = Parameter(name = 'Msu3',
                 nature = 'external',
                 type = 'real',
                 value = 399.668493,
                 texname = '\\text{Msu3}',
                 lhablock = 'MASS',
                 lhacode = [ 1000006 ])

Msu4 = Parameter(name = 'Msu4',
                 nature = 'external',
                 type = 'real',
                 value = 549.259265,
                 texname = '\\text{Msu4}',
                 lhablock = 'MASS',
                 lhacode = [ 2000002 ])

Msu5 = Parameter(name = 'Msu5',
                 nature = 'external',
                 type = 'real',
                 value = 549.259265,
                 texname = '\\text{Msu5}',
                 lhablock = 'MASS',
                 lhacode = [ 2000004 ])

Msu6 = Parameter(name = 'Msu6',
                 nature = 'external',
                 type = 'real',
                 value = 585.785818,
                 texname = '\\text{Msu6}',
                 lhablock = 'MASS',
                 lhacode = [ 2000006 ])

Msd1 = Parameter(name = 'Msd1',
                 nature = 'external',
                 type = 'real',
                 value = 568.441109,
                 texname = '\\text{Msd1}',
                 lhablock = 'MASS',
                 lhacode = [ 1000001 ])

Msd2 = Parameter(name = 'Msd2',
                 nature = 'external',
                 type = 'real',
                 value = 568.441109,
                 texname = '\\text{Msd2}',
                 lhablock = 'MASS',
                 lhacode = [ 1000003 ])

Msd3 = Parameter(name = 'Msd3',
                 nature = 'external',
                 type = 'real',
                 value = 513.065179,
                 texname = '\\text{Msd3}',
                 lhablock = 'MASS',
                 lhacode = [ 1000005 ])

Msd4 = Parameter(name = 'Msd4',
                 nature = 'external',
                 type = 'real',
                 value = 545.228462,
                 texname = '\\text{Msd4}',
                 lhablock = 'MASS',
                 lhacode = [ 2000001 ])

Msd5 = Parameter(name = 'Msd5',
                 nature = 'external',
                 type = 'real',
                 value = 545.228462,
                 texname = '\\text{Msd5}',
                 lhablock = 'MASS',
                 lhacode = [ 2000003 ])

Msd6 = Parameter(name = 'Msd6',
                 nature = 'external',
                 type = 'real',
                 value = 543.726676,
                 texname = '\\text{Msd6}',
                 lhablock = 'MASS',
                 lhacode = [ 2000005 ])

WZ = Parameter(name = 'WZ',
               nature = 'external',
               type = 'real',
               value = 2.41143316,
               texname = '\\text{WZ}',
               lhablock = 'DECAY',
               lhacode = [ 23 ])

WW = Parameter(name = 'WW',
               nature = 'external',
               type = 'real',
               value = 2.00282196,
               texname = '\\text{WW}',
               lhablock = 'DECAY',
               lhacode = [ 24 ])

Wneu2 = Parameter(name = 'Wneu2',
                  nature = 'external',
                  type = 'real',
                  value = 0.0207770048,
                  texname = '\\text{Wneu2}',
                  lhablock = 'DECAY',
                  lhacode = [ 1000023 ])

Wneu3 = Parameter(name = 'Wneu3',
                  nature = 'external',
                  type = 'real',
                  value = 1.91598495,
                  texname = '\\text{Wneu3}',
                  lhablock = 'DECAY',
                  lhacode = [ 1000025 ])

Wneu4 = Parameter(name = 'Wneu4',
                  nature = 'external',
                  type = 'real',
                  value = 2.58585079,
                  texname = '\\text{Wneu4}',
                  lhablock = 'DECAY',
                  lhacode = [ 1000035 ])

Wch1 = Parameter(name = 'Wch1',
                 nature = 'external',
                 type = 'real',
                 value = 0.0170414503,
                 texname = '\\text{Wch1}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000024 ])

Wch2 = Parameter(name = 'Wch2',
                 nature = 'external',
                 type = 'real',
                 value = 2.4868951,
                 texname = '\\text{Wch2}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000037 ])

Wgo = Parameter(name = 'Wgo',
                nature = 'external',
                type = 'real',
                value = 5.50675438,
                texname = '\\text{Wgo}',
                lhablock = 'DECAY',
                lhacode = [ 1000021 ])

WH01 = Parameter(name = 'WH01',
                 nature = 'external',
                 type = 'real',
                 value = 0.00198610799,
                 texname = '\\text{WH01}',
                 lhablock = 'DECAY',
                 lhacode = [ 25 ])

WH02 = Parameter(name = 'WH02',
                 nature = 'external',
                 type = 'real',
                 value = 0.574801389,
                 texname = '\\text{WH02}',
                 lhablock = 'DECAY',
                 lhacode = [ 35 ])

WA0 = Parameter(name = 'WA0',
                nature = 'external',
                type = 'real',
                value = 0.632178488,
                texname = '\\text{WA0}',
                lhablock = 'DECAY',
                lhacode = [ 36 ])

WH = Parameter(name = 'WH',
               nature = 'external',
               type = 'real',
               value = 0.546962813,
               texname = '\\text{WH}',
               lhablock = 'DECAY',
               lhacode = [ 37 ])

WT = Parameter(name = 'WT',
               nature = 'external',
               type = 'real',
               value = 1.56194983,
               texname = '\\text{WT}',
               lhablock = 'DECAY',
               lhacode = [ 6 ])

Wsn1 = Parameter(name = 'Wsn1',
                 nature = 'external',
                 type = 'real',
                 value = 0.149881634,
                 texname = '\\text{Wsn1}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000012 ])

Wsn2 = Parameter(name = 'Wsn2',
                 nature = 'external',
                 type = 'real',
                 value = 0.149881634,
                 texname = '\\text{Wsn2}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000014 ])

Wsn3 = Parameter(name = 'Wsn3',
                 nature = 'external',
                 type = 'real',
                 value = 0.147518977,
                 texname = '\\text{Wsn3}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000016 ])

Wsl1 = Parameter(name = 'Wsl1',
                 nature = 'external',
                 type = 'real',
                 value = 0.213682161,
                 texname = '\\text{Wsl1}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000011 ])

Wsl2 = Parameter(name = 'Wsl2',
                 nature = 'external',
                 type = 'real',
                 value = 0.213682161,
                 texname = '\\text{Wsl2}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000013 ])

Wsl3 = Parameter(name = 'Wsl3',
                 nature = 'external',
                 type = 'real',
                 value = 0.148327268,
                 texname = '\\text{Wsl3}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000015 ])

Wsl4 = Parameter(name = 'Wsl4',
                 nature = 'external',
                 type = 'real',
                 value = 0.216121626,
                 texname = '\\text{Wsl4}',
                 lhablock = 'DECAY',
                 lhacode = [ 2000011 ])

Wsl5 = Parameter(name = 'Wsl5',
                 nature = 'external',
                 type = 'real',
                 value = 0.216121626,
                 texname = '\\text{Wsl5}',
                 lhablock = 'DECAY',
                 lhacode = [ 2000013 ])

Wsl6 = Parameter(name = 'Wsl6',
                 nature = 'external',
                 type = 'real',
                 value = 0.269906096,
                 texname = '\\text{Wsl6}',
                 lhablock = 'DECAY',
                 lhacode = [ 2000015 ])

Wsu1 = Parameter(name = 'Wsu1',
                 nature = 'external',
                 type = 'real',
                 value = 5.47719539,
                 texname = '\\text{Wsu1}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000002 ])

Wsu2 = Parameter(name = 'Wsu2',
                 nature = 'external',
                 type = 'real',
                 value = 5.47719539,
                 texname = '\\text{Wsu2}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000004 ])

Wsu3 = Parameter(name = 'Wsu3',
                 nature = 'external',
                 type = 'real',
                 value = 2.02159578,
                 texname = '\\text{Wsu3}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000006 ])

Wsu4 = Parameter(name = 'Wsu4',
                 nature = 'external',
                 type = 'real',
                 value = 1.15297292,
                 texname = '\\text{Wsu4}',
                 lhablock = 'DECAY',
                 lhacode = [ 2000002 ])

Wsu5 = Parameter(name = 'Wsu5',
                 nature = 'external',
                 type = 'real',
                 value = 1.15297292,
                 texname = '\\text{Wsu5}',
                 lhablock = 'DECAY',
                 lhacode = [ 2000004 ])

Wsu6 = Parameter(name = 'Wsu6',
                 nature = 'external',
                 type = 'real',
                 value = 7.37313275,
                 texname = '\\text{Wsu6}',
                 lhablock = 'DECAY',
                 lhacode = [ 2000006 ])

Wsd1 = Parameter(name = 'Wsd1',
                 nature = 'external',
                 type = 'real',
                 value = 5.31278772,
                 texname = '\\text{Wsd1}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000001 ])

Wsd2 = Parameter(name = 'Wsd2',
                 nature = 'external',
                 type = 'real',
                 value = 5.31278772,
                 texname = '\\text{Wsd2}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000003 ])

Wsd3 = Parameter(name = 'Wsd3',
                 nature = 'external',
                 type = 'real',
                 value = 3.73627601,
                 texname = '\\text{Wsd3}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000005 ])

Wsd4 = Parameter(name = 'Wsd4',
                 nature = 'external',
                 type = 'real',
                 value = 0.285812308,
                 texname = '\\text{Wsd4}',
                 lhablock = 'DECAY',
                 lhacode = [ 2000001 ])

Wsd5 = Parameter(name = 'Wsd5',
                 nature = 'external',
                 type = 'real',
                 value = 0.285812308,
                 texname = '\\text{Wsd5}',
                 lhablock = 'DECAY',
                 lhacode = [ 2000003 ])

Wsd6 = Parameter(name = 'Wsd6',
                 nature = 'external',
                 type = 'real',
                 value = 0.801566294,
                 texname = '\\text{Wsd6}',
                 lhablock = 'DECAY',
                 lhacode = [ 2000005 ])

beta = Parameter(name = 'beta',
                 nature = 'internal',
                 type = 'real',
                 value = 'cmath.atan(tb)',
                 texname = '\\beta')

CKM1x1 = Parameter(name = 'CKM1x1',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RCKM1x1',
                   texname = '\\text{CKM1x1}')

CKM2x2 = Parameter(name = 'CKM2x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RCKM2x2',
                   texname = '\\text{CKM2x2}')

CKM3x3 = Parameter(name = 'CKM3x3',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RCKM3x3',
                   texname = '\\text{CKM3x3}')

cw = Parameter(name = 'cw',
               nature = 'internal',
               type = 'real',
               value = 'MW/MZ',
               texname = 'c_w')

mD21x1 = Parameter(name = 'mD21x1',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmD21x1',
                   texname = '\\text{mD21x1}')

mD22x2 = Parameter(name = 'mD22x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmD22x2',
                   texname = '\\text{mD22x2}')

mD23x3 = Parameter(name = 'mD23x3',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmD23x3',
                   texname = '\\text{mD23x3}')

mE21x1 = Parameter(name = 'mE21x1',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmE21x1',
                   texname = '\\text{mE21x1}')

mE22x2 = Parameter(name = 'mE22x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmE22x2',
                   texname = '\\text{mE22x2}')

mE23x3 = Parameter(name = 'mE23x3',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmE23x3',
                   texname = '\\text{mE23x3}')

mL21x1 = Parameter(name = 'mL21x1',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmL21x1',
                   texname = '\\text{mL21x1}')

mL22x2 = Parameter(name = 'mL22x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmL22x2',
                   texname = '\\text{mL22x2}')

mL23x3 = Parameter(name = 'mL23x3',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmL23x3',
                   texname = '\\text{mL23x3}')

mQ21x1 = Parameter(name = 'mQ21x1',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmQ21x1',
                   texname = '\\text{mQ21x1}')

mQ22x2 = Parameter(name = 'mQ22x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmQ22x2',
                   texname = '\\text{mQ22x2}')

mQ23x3 = Parameter(name = 'mQ23x3',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmQ23x3',
                   texname = '\\text{mQ23x3}')

mU21x1 = Parameter(name = 'mU21x1',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmU21x1',
                   texname = '\\text{mU21x1}')

mU22x2 = Parameter(name = 'mU22x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmU22x2',
                   texname = '\\text{mU22x2}')

mU23x3 = Parameter(name = 'mU23x3',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmU23x3',
                   texname = '\\text{mU23x3}')

MUH = Parameter(name = 'MUH',
                nature = 'internal',
                type = 'complex',
                value = 'RMUH',
                texname = '\\mu')

Mx1 = Parameter(name = 'Mx1',
                nature = 'internal',
                type = 'complex',
                value = 'RMx1',
                texname = 'M_1')

Mx2 = Parameter(name = 'Mx2',
                nature = 'internal',
                type = 'complex',
                value = 'RMx2',
                texname = 'M_2')

Mx3 = Parameter(name = 'Mx3',
                nature = 'internal',
                type = 'complex',
                value = 'RMx3',
                texname = 'M_3')

NN1x1 = Parameter(name = 'NN1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complex(0,1)*INN1x1 + RNN1x1',
                  texname = '\\text{NN1x1}')

NN1x2 = Parameter(name = 'NN1x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complex(0,1)*INN1x2 + RNN1x2',
                  texname = '\\text{NN1x2}')

NN1x3 = Parameter(name = 'NN1x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complex(0,1)*INN1x3 + RNN1x3',
                  texname = '\\text{NN1x3}')

NN1x4 = Parameter(name = 'NN1x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complex(0,1)*INN1x4 + RNN1x4',
                  texname = '\\text{NN1x4}')

NN2x1 = Parameter(name = 'NN2x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complex(0,1)*INN2x1 + RNN2x1',
                  texname = '\\text{NN2x1}')

NN2x2 = Parameter(name = 'NN2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complex(0,1)*INN2x2 + RNN2x2',
                  texname = '\\text{NN2x2}')

NN2x3 = Parameter(name = 'NN2x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complex(0,1)*INN2x3 + RNN2x3',
                  texname = '\\text{NN2x3}')

NN2x4 = Parameter(name = 'NN2x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complex(0,1)*INN2x4 + RNN2x4',
                  texname = '\\text{NN2x4}')

NN3x1 = Parameter(name = 'NN3x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complex(0,1)*INN3x1 + RNN3x1',
                  texname = '\\text{NN3x1}')

NN3x2 = Parameter(name = 'NN3x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complex(0,1)*INN3x2 + RNN3x2',
                  texname = '\\text{NN3x2}')

NN3x3 = Parameter(name = 'NN3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complex(0,1)*INN3x3 + RNN3x3',
                  texname = '\\text{NN3x3}')

NN3x4 = Parameter(name = 'NN3x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complex(0,1)*INN3x4 + RNN3x4',
                  texname = '\\text{NN3x4}')

NN4x1 = Parameter(name = 'NN4x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complex(0,1)*INN4x1 + RNN4x1',
                  texname = '\\text{NN4x1}')

NN4x2 = Parameter(name = 'NN4x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complex(0,1)*INN4x2 + RNN4x2',
                  texname = '\\text{NN4x2}')

NN4x3 = Parameter(name = 'NN4x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complex(0,1)*INN4x3 + RNN4x3',
                  texname = '\\text{NN4x3}')

NN4x4 = Parameter(name = 'NN4x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complex(0,1)*INN4x4 + RNN4x4',
                  texname = '\\text{NN4x4}')

Rd1x1 = Parameter(name = 'Rd1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd1x1',
                  texname = '\\text{Rd1x1}')

Rd2x2 = Parameter(name = 'Rd2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd2x2',
                  texname = '\\text{Rd2x2}')

Rd3x3 = Parameter(name = 'Rd3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd3x3',
                  texname = '\\text{Rd3x3}')

Rd3x6 = Parameter(name = 'Rd3x6',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd3x6',
                  texname = '\\text{Rd3x6}')

Rd4x4 = Parameter(name = 'Rd4x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd4x4',
                  texname = '\\text{Rd4x4}')

Rd5x5 = Parameter(name = 'Rd5x5',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd5x5',
                  texname = '\\text{Rd5x5}')

Rd6x3 = Parameter(name = 'Rd6x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd6x3',
                  texname = '\\text{Rd6x3}')

Rd6x6 = Parameter(name = 'Rd6x6',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd6x6',
                  texname = '\\text{Rd6x6}')

Rl1x1 = Parameter(name = 'Rl1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl1x1',
                  texname = '\\text{Rl1x1}')

Rl2x2 = Parameter(name = 'Rl2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl2x2',
                  texname = '\\text{Rl2x2}')

Rl3x3 = Parameter(name = 'Rl3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl3x3',
                  texname = '\\text{Rl3x3}')

Rl3x6 = Parameter(name = 'Rl3x6',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl3x6',
                  texname = '\\text{Rl3x6}')

Rl4x4 = Parameter(name = 'Rl4x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl4x4',
                  texname = '\\text{Rl4x4}')

Rl5x5 = Parameter(name = 'Rl5x5',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl5x5',
                  texname = '\\text{Rl5x5}')

Rl6x3 = Parameter(name = 'Rl6x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl6x3',
                  texname = '\\text{Rl6x3}')

Rl6x6 = Parameter(name = 'Rl6x6',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl6x6',
                  texname = '\\text{Rl6x6}')

Rn1x1 = Parameter(name = 'Rn1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRn1x1',
                  texname = '\\text{Rn1x1}')

Rn2x2 = Parameter(name = 'Rn2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRn2x2',
                  texname = '\\text{Rn2x2}')

Rn3x3 = Parameter(name = 'Rn3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRn3x3',
                  texname = '\\text{Rn3x3}')

Ru1x1 = Parameter(name = 'Ru1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu1x1',
                  texname = '\\text{Ru1x1}')

Ru2x2 = Parameter(name = 'Ru2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu2x2',
                  texname = '\\text{Ru2x2}')

Ru3x3 = Parameter(name = 'Ru3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu3x3',
                  texname = '\\text{Ru3x3}')

Ru3x6 = Parameter(name = 'Ru3x6',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu3x6',
                  texname = '\\text{Ru3x6}')

Ru4x4 = Parameter(name = 'Ru4x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu4x4',
                  texname = '\\text{Ru4x4}')

Ru5x5 = Parameter(name = 'Ru5x5',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu5x5',
                  texname = '\\text{Ru5x5}')

Ru6x3 = Parameter(name = 'Ru6x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu6x3',
                  texname = '\\text{Ru6x3}')

Ru6x6 = Parameter(name = 'Ru6x6',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu6x6',
                  texname = '\\text{Ru6x6}')

UU1x1 = Parameter(name = 'UU1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RUU1x1',
                  texname = '\\text{UU1x1}')

UU1x2 = Parameter(name = 'UU1x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RUU1x2',
                  texname = '\\text{UU1x2}')

UU2x1 = Parameter(name = 'UU2x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RUU2x1',
                  texname = '\\text{UU2x1}')

UU2x2 = Parameter(name = 'UU2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RUU2x2',
                  texname = '\\text{UU2x2}')

VV1x1 = Parameter(name = 'VV1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RVV1x1',
                  texname = '\\text{VV1x1}')

VV1x2 = Parameter(name = 'VV1x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RVV1x2',
                  texname = '\\text{VV1x2}')

VV2x1 = Parameter(name = 'VV2x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RVV2x1',
                  texname = '\\text{VV2x1}')

VV2x2 = Parameter(name = 'VV2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RVV2x2',
                  texname = '\\text{VV2x2}')

ee = Parameter(name = 'ee',
               nature = 'internal',
               type = 'real',
               value = '2*cmath.sqrt(1/aEWM1)*cmath.sqrt(cmath.pi)',
               texname = 'e')

G = Parameter(name = 'G',
              nature = 'internal',
              type = 'real',
              value = '2*cmath.sqrt(aS)*cmath.sqrt(cmath.pi)',
              texname = 'G')

td3x3 = Parameter(name = 'td3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rtd3x3',
                  texname = '\\text{td3x3}')

te3x3 = Parameter(name = 'te3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rte3x3',
                  texname = '\\text{te3x3}')

tu3x3 = Parameter(name = 'tu3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rtu3x3',
                  texname = '\\text{tu3x3}')

yd3x3 = Parameter(name = 'yd3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Ryd3x3',
                  texname = '\\text{yd3x3}')

ye3x3 = Parameter(name = 'ye3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rye3x3',
                  texname = '\\text{ye3x3}')

yu3x3 = Parameter(name = 'yu3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Ryu3x3',
                  texname = '\\text{yu3x3}')

bb = Parameter(name = 'bb',
               nature = 'internal',
               type = 'complex',
               value = '((-mHd2 + mHu2)*cmath.tan(2*alp))/2. - MZ**2*(cmath.sin(2*beta)/2. + cmath.cos(2*beta)*cmath.tan(2*alp))',
               texname = 'b')

sw = Parameter(name = 'sw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(1 - cw**2)',
               texname = 's_w')

gp = Parameter(name = 'gp',
               nature = 'internal',
               type = 'real',
               value = 'ee/cw',
               texname = 'g\'')

gw = Parameter(name = 'gw',
               nature = 'internal',
               type = 'real',
               value = 'ee/sw',
               texname = 'g_w')

vev = Parameter(name = 'vev',
                nature = 'internal',
                type = 'real',
                value = '(2*cw*MZ*sw)/ee',
                texname = 'v')

vd = Parameter(name = 'vd',
               nature = 'internal',
               type = 'real',
               value = 'vev*cmath.cos(beta)',
               texname = 'v_d')

vu = Parameter(name = 'vu',
               nature = 'internal',
               type = 'real',
               value = 'vev*cmath.sin(beta)',
               texname = 'v_u')

I1a33 = Parameter(name = 'I1a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complexconjugate(CKM3x3)*complexconjugate(yu3x3)',
                  texname = '\\text{I1a33}')

I10a33 = Parameter(name = 'I10a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(yd3x3)',
                   texname = '\\text{I10a33}')

I10a36 = Parameter(name = 'I10a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(yd3x3)',
                   texname = '\\text{I10a36}')

I100a11 = Parameter(name = 'I100a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*complexconjugate(Rl1x1)',
                    texname = '\\text{I100a11}')

I100a22 = Parameter(name = 'I100a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl2x2)',
                    texname = '\\text{I100a22}')

I100a33 = Parameter(name = 'I100a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*complexconjugate(Rl3x3)',
                    texname = '\\text{I100a33}')

I100a36 = Parameter(name = 'I100a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*complexconjugate(Rl6x3)',
                    texname = '\\text{I100a36}')

I101a33 = Parameter(name = 'I101a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                    texname = '\\text{I101a33}')

I101a36 = Parameter(name = 'I101a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                    texname = '\\text{I101a36}')

I102a33 = Parameter(name = 'I102a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rl3x6*ye3x3*complexconjugate(Rn3x3)',
                    texname = '\\text{I102a33}')

I102a36 = Parameter(name = 'I102a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rl6x6*ye3x3*complexconjugate(Rn3x3)',
                    texname = '\\text{I102a36}')

I103a11 = Parameter(name = 'I103a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1',
                    texname = '\\text{I103a11}')

I103a22 = Parameter(name = 'I103a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2',
                    texname = '\\text{I103a22}')

I103a33 = Parameter(name = 'I103a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3',
                    texname = '\\text{I103a33}')

I104a33 = Parameter(name = 'I104a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*complexconjugate(ye3x3)',
                    texname = '\\text{I104a33}')

I105a11 = Parameter(name = 'I105a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*complexconjugate(Rl1x1)',
                    texname = '\\text{I105a11}')

I105a22 = Parameter(name = 'I105a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl2x2)',
                    texname = '\\text{I105a22}')

I105a33 = Parameter(name = 'I105a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*complexconjugate(Rl3x3)',
                    texname = '\\text{I105a33}')

I105a36 = Parameter(name = 'I105a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*complexconjugate(Rl6x3)',
                    texname = '\\text{I105a36}')

I106a33 = Parameter(name = 'I106a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*complexconjugate(Rl3x6)*complexconjugate(te3x3)',
                    texname = '\\text{I106a33}')

I106a36 = Parameter(name = 'I106a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*complexconjugate(Rl6x6)*complexconjugate(te3x3)',
                    texname = '\\text{I106a36}')

I107a33 = Parameter(name = 'I107a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                    texname = '\\text{I107a33}')

I107a36 = Parameter(name = 'I107a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                    texname = '\\text{I107a36}')

I108a33 = Parameter(name = 'I108a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*ye3x3*complexconjugate(Rl3x3)*complexconjugate(ye3x3)',
                    texname = '\\text{I108a33}')

I108a36 = Parameter(name = 'I108a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*ye3x3*complexconjugate(Rl6x3)*complexconjugate(ye3x3)',
                    texname = '\\text{I108a36}')

I109a11 = Parameter(name = 'I109a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*complexconjugate(Rl1x1)',
                    texname = '\\text{I109a11}')

I109a22 = Parameter(name = 'I109a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl2x2)',
                    texname = '\\text{I109a22}')

I109a33 = Parameter(name = 'I109a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*complexconjugate(Rl3x3)',
                    texname = '\\text{I109a33}')

I109a36 = Parameter(name = 'I109a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*complexconjugate(Rl6x3)',
                    texname = '\\text{I109a36}')

I11a33 = Parameter(name = 'I11a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3',
                   texname = '\\text{I11a33}')

I11a36 = Parameter(name = 'I11a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3',
                   texname = '\\text{I11a36}')

I110a33 = Parameter(name = 'I110a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*ye3x3*complexconjugate(Rl3x3)*complexconjugate(ye3x3)',
                    texname = '\\text{I110a33}')

I110a36 = Parameter(name = 'I110a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*ye3x3*complexconjugate(Rl6x3)*complexconjugate(ye3x3)',
                    texname = '\\text{I110a36}')

I111a33 = Parameter(name = 'I111a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I111a33}')

I111a36 = Parameter(name = 'I111a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I111a36}')

I112a33 = Parameter(name = 'I112a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'yu3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I112a33}')

I112a36 = Parameter(name = 'I112a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'yu3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I112a36}')

I113a11 = Parameter(name = 'I113a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I113a11}')

I113a22 = Parameter(name = 'I113a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I113a22}')

I113a33 = Parameter(name = 'I113a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I113a33}')

I113a36 = Parameter(name = 'I113a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I113a36}')

I113a63 = Parameter(name = 'I113a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I113a63}')

I113a66 = Parameter(name = 'I113a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I113a66}')

I114a33 = Parameter(name = 'I114a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I114a33}')

I114a36 = Parameter(name = 'I114a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I114a36}')

I114a44 = Parameter(name = 'I114a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I114a44}')

I114a55 = Parameter(name = 'I114a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I114a55}')

I114a63 = Parameter(name = 'I114a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I114a63}')

I114a66 = Parameter(name = 'I114a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I114a66}')

I115a11 = Parameter(name = 'I115a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                    texname = '\\text{I115a11}')

I115a22 = Parameter(name = 'I115a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I115a22}')

I115a33 = Parameter(name = 'I115a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                    texname = '\\text{I115a33}')

I115a36 = Parameter(name = 'I115a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                    texname = '\\text{I115a36}')

I115a63 = Parameter(name = 'I115a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                    texname = '\\text{I115a63}')

I115a66 = Parameter(name = 'I115a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                    texname = '\\text{I115a66}')

I116a33 = Parameter(name = 'I116a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I116a33}')

I116a36 = Parameter(name = 'I116a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I116a36}')

I116a63 = Parameter(name = 'I116a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I116a63}')

I116a66 = Parameter(name = 'I116a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I116a66}')

I117a33 = Parameter(name = 'I117a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x6)*complexconjugate(tu3x3)',
                    texname = '\\text{I117a33}')

I117a36 = Parameter(name = 'I117a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x6)*complexconjugate(tu3x3)',
                    texname = '\\text{I117a36}')

I117a63 = Parameter(name = 'I117a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x6)*complexconjugate(tu3x3)',
                    texname = '\\text{I117a63}')

I117a66 = Parameter(name = 'I117a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x6)*complexconjugate(tu3x3)',
                    texname = '\\text{I117a66}')

I118a33 = Parameter(name = 'I118a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x6*td3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                    texname = '\\text{I118a33}')

I118a36 = Parameter(name = 'I118a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x6*td3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                    texname = '\\text{I118a36}')

I118a63 = Parameter(name = 'I118a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x6*td3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                    texname = '\\text{I118a63}')

I118a66 = Parameter(name = 'I118a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x6*td3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                    texname = '\\text{I118a66}')

I119a33 = Parameter(name = 'I119a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)*complexconjugate(yd3x3)',
                    texname = '\\text{I119a33}')

I119a36 = Parameter(name = 'I119a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)*complexconjugate(yd3x3)',
                    texname = '\\text{I119a36}')

I119a63 = Parameter(name = 'I119a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)*complexconjugate(yd3x3)',
                    texname = '\\text{I119a63}')

I119a66 = Parameter(name = 'I119a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)*complexconjugate(yd3x3)',
                    texname = '\\text{I119a66}')

I12a11 = Parameter(name = 'I12a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(CKM1x1)',
                   texname = '\\text{I12a11}')

I12a22 = Parameter(name = 'I12a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(CKM2x2)',
                   texname = '\\text{I12a22}')

I12a33 = Parameter(name = 'I12a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(CKM3x3)',
                   texname = '\\text{I12a33}')

I12a36 = Parameter(name = 'I12a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(CKM3x3)',
                   texname = '\\text{I12a36}')

I120a33 = Parameter(name = 'I120a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                    texname = '\\text{I120a33}')

I120a36 = Parameter(name = 'I120a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                    texname = '\\text{I120a36}')

I120a63 = Parameter(name = 'I120a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                    texname = '\\text{I120a63}')

I120a66 = Parameter(name = 'I120a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                    texname = '\\text{I120a66}')

I121a33 = Parameter(name = 'I121a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I121a33}')

I121a36 = Parameter(name = 'I121a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I121a36}')

I121a63 = Parameter(name = 'I121a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I121a63}')

I121a66 = Parameter(name = 'I121a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I121a66}')

I122a33 = Parameter(name = 'I122a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*yu3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I122a33}')

I122a36 = Parameter(name = 'I122a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*yu3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I122a36}')

I122a63 = Parameter(name = 'I122a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*yu3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I122a63}')

I122a66 = Parameter(name = 'I122a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*yu3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I122a66}')

I123a11 = Parameter(name = 'I123a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                    texname = '\\text{I123a11}')

I123a22 = Parameter(name = 'I123a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I123a22}')

I123a33 = Parameter(name = 'I123a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                    texname = '\\text{I123a33}')

I123a36 = Parameter(name = 'I123a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                    texname = '\\text{I123a36}')

I123a63 = Parameter(name = 'I123a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                    texname = '\\text{I123a63}')

I123a66 = Parameter(name = 'I123a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                    texname = '\\text{I123a66}')

I124a33 = Parameter(name = 'I124a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)*complexconjugate(yd3x3)',
                    texname = '\\text{I124a33}')

I124a36 = Parameter(name = 'I124a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)*complexconjugate(yd3x3)',
                    texname = '\\text{I124a36}')

I124a63 = Parameter(name = 'I124a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)*complexconjugate(yd3x3)',
                    texname = '\\text{I124a63}')

I124a66 = Parameter(name = 'I124a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)*complexconjugate(yd3x3)',
                    texname = '\\text{I124a66}')

I125a33 = Parameter(name = 'I125a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I125a33}')

I125a36 = Parameter(name = 'I125a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I125a36}')

I125a63 = Parameter(name = 'I125a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I125a63}')

I125a66 = Parameter(name = 'I125a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I125a66}')

I126a33 = Parameter(name = 'I126a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*yu3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I126a33}')

I126a36 = Parameter(name = 'I126a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*yu3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I126a36}')

I126a63 = Parameter(name = 'I126a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*yu3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I126a63}')

I126a66 = Parameter(name = 'I126a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*yu3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I126a66}')

I127a11 = Parameter(name = 'I127a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                    texname = '\\text{I127a11}')

I127a22 = Parameter(name = 'I127a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I127a22}')

I127a33 = Parameter(name = 'I127a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                    texname = '\\text{I127a33}')

I127a36 = Parameter(name = 'I127a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                    texname = '\\text{I127a36}')

I127a63 = Parameter(name = 'I127a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                    texname = '\\text{I127a63}')

I127a66 = Parameter(name = 'I127a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                    texname = '\\text{I127a66}')

I128a11 = Parameter(name = 'I128a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*complexconjugate(Rl1x1)',
                    texname = '\\text{I128a11}')

I128a22 = Parameter(name = 'I128a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl2x2)',
                    texname = '\\text{I128a22}')

I128a33 = Parameter(name = 'I128a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*complexconjugate(Rl3x3)',
                    texname = '\\text{I128a33}')

I128a36 = Parameter(name = 'I128a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*complexconjugate(Rl6x3)',
                    texname = '\\text{I128a36}')

I129a33 = Parameter(name = 'I129a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                    texname = '\\text{I129a33}')

I129a36 = Parameter(name = 'I129a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                    texname = '\\text{I129a36}')

I13a33 = Parameter(name = 'I13a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(yu3x3)',
                   texname = '\\text{I13a33}')

I13a36 = Parameter(name = 'I13a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(yu3x3)',
                   texname = '\\text{I13a36}')

I130a33 = Parameter(name = 'I130a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                    texname = '\\text{I130a33}')

I130a36 = Parameter(name = 'I130a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                    texname = '\\text{I130a36}')

I130a63 = Parameter(name = 'I130a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                    texname = '\\text{I130a63}')

I130a66 = Parameter(name = 'I130a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                    texname = '\\text{I130a66}')

I131a33 = Parameter(name = 'I131a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(yu3x3)',
                    texname = '\\text{I131a33}')

I131a36 = Parameter(name = 'I131a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(yu3x3)',
                    texname = '\\text{I131a36}')

I132a33 = Parameter(name = 'I132a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*yu3x3',
                    texname = '\\text{I132a33}')

I132a36 = Parameter(name = 'I132a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*yu3x3',
                    texname = '\\text{I132a36}')

I133a11 = Parameter(name = 'I133a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1',
                    texname = '\\text{I133a11}')

I133a22 = Parameter(name = 'I133a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2',
                    texname = '\\text{I133a22}')

I133a33 = Parameter(name = 'I133a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3',
                    texname = '\\text{I133a33}')

I133a36 = Parameter(name = 'I133a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3',
                    texname = '\\text{I133a36}')

I134a33 = Parameter(name = 'I134a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*complexconjugate(yd3x3)',
                    texname = '\\text{I134a33}')

I134a36 = Parameter(name = 'I134a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*complexconjugate(yd3x3)',
                    texname = '\\text{I134a36}')

I135a33 = Parameter(name = 'I135a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x6*yu3x3',
                    texname = '\\text{I135a33}')

I135a36 = Parameter(name = 'I135a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x6*yu3x3',
                    texname = '\\text{I135a36}')

I136a11 = Parameter(name = 'I136a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*complexconjugate(Rd1x1)',
                    texname = '\\text{I136a11}')

I136a22 = Parameter(name = 'I136a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I136a22}')

I136a33 = Parameter(name = 'I136a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*complexconjugate(Rd3x3)',
                    texname = '\\text{I136a33}')

I136a36 = Parameter(name = 'I136a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*complexconjugate(Rd3x3)',
                    texname = '\\text{I136a36}')

I136a63 = Parameter(name = 'I136a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*complexconjugate(Rd6x3)',
                    texname = '\\text{I136a63}')

I136a66 = Parameter(name = 'I136a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*complexconjugate(Rd6x3)',
                    texname = '\\text{I136a66}')

I137a33 = Parameter(name = 'I137a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*complexconjugate(Rd3x6)*complexconjugate(td3x3)',
                    texname = '\\text{I137a33}')

I137a36 = Parameter(name = 'I137a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*complexconjugate(Rd3x6)*complexconjugate(td3x3)',
                    texname = '\\text{I137a36}')

I137a63 = Parameter(name = 'I137a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*complexconjugate(Rd6x6)*complexconjugate(td3x3)',
                    texname = '\\text{I137a63}')

I137a66 = Parameter(name = 'I137a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*complexconjugate(Rd6x6)*complexconjugate(td3x3)',
                    texname = '\\text{I137a66}')

I138a33 = Parameter(name = 'I138a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I138a33}')

I138a36 = Parameter(name = 'I138a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I138a36}')

I138a63 = Parameter(name = 'I138a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I138a63}')

I138a66 = Parameter(name = 'I138a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I138a66}')

I139a33 = Parameter(name = 'I139a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x6*tu3x3*complexconjugate(Rd3x3)',
                    texname = '\\text{I139a33}')

I139a36 = Parameter(name = 'I139a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x6*tu3x3*complexconjugate(Rd3x3)',
                    texname = '\\text{I139a36}')

I139a63 = Parameter(name = 'I139a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x6*tu3x3*complexconjugate(Rd6x3)',
                    texname = '\\text{I139a63}')

I139a66 = Parameter(name = 'I139a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x6*tu3x3*complexconjugate(Rd6x3)',
                    texname = '\\text{I139a66}')

I14a33 = Parameter(name = 'I14a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(CKM3x3)',
                   texname = '\\text{I14a33}')

I14a36 = Parameter(name = 'I14a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(CKM3x3)',
                   texname = '\\text{I14a36}')

I140a33 = Parameter(name = 'I140a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*yd3x3*complexconjugate(Rd3x3)*complexconjugate(yd3x3)',
                    texname = '\\text{I140a33}')

I140a36 = Parameter(name = 'I140a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*yd3x3*complexconjugate(Rd3x3)*complexconjugate(yd3x3)',
                    texname = '\\text{I140a36}')

I140a63 = Parameter(name = 'I140a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*yd3x3*complexconjugate(Rd6x3)*complexconjugate(yd3x3)',
                    texname = '\\text{I140a63}')

I140a66 = Parameter(name = 'I140a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*yd3x3*complexconjugate(Rd6x3)*complexconjugate(yd3x3)',
                    texname = '\\text{I140a66}')

I141a33 = Parameter(name = 'I141a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x6*yu3x3*complexconjugate(Rd3x3)',
                    texname = '\\text{I141a33}')

I141a36 = Parameter(name = 'I141a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x6*yu3x3*complexconjugate(Rd3x3)',
                    texname = '\\text{I141a36}')

I141a63 = Parameter(name = 'I141a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x6*yu3x3*complexconjugate(Rd6x3)',
                    texname = '\\text{I141a63}')

I141a66 = Parameter(name = 'I141a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x6*yu3x3*complexconjugate(Rd6x3)',
                    texname = '\\text{I141a66}')

I142a33 = Parameter(name = 'I142a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*yu3x3*complexconjugate(Rd3x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I142a33}')

I142a36 = Parameter(name = 'I142a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*yu3x3*complexconjugate(Rd3x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I142a36}')

I142a63 = Parameter(name = 'I142a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*yu3x3*complexconjugate(Rd6x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I142a63}')

I142a66 = Parameter(name = 'I142a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*yu3x3*complexconjugate(Rd6x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I142a66}')

I143a33 = Parameter(name = 'I143a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x6*yu3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I143a33}')

I143a36 = Parameter(name = 'I143a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x6*yu3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I143a36}')

I143a63 = Parameter(name = 'I143a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x6*yu3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I143a63}')

I143a66 = Parameter(name = 'I143a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x6*yu3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I143a66}')

I144a11 = Parameter(name = 'I144a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*complexconjugate(Rd1x1)',
                    texname = '\\text{I144a11}')

I144a22 = Parameter(name = 'I144a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I144a22}')

I144a33 = Parameter(name = 'I144a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*complexconjugate(Rd3x3)',
                    texname = '\\text{I144a33}')

I144a36 = Parameter(name = 'I144a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*complexconjugate(Rd3x3)',
                    texname = '\\text{I144a36}')

I144a63 = Parameter(name = 'I144a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*complexconjugate(Rd6x3)',
                    texname = '\\text{I144a63}')

I144a66 = Parameter(name = 'I144a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*complexconjugate(Rd6x3)',
                    texname = '\\text{I144a66}')

I145a33 = Parameter(name = 'I145a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*yd3x3*complexconjugate(Rd3x3)*complexconjugate(yd3x3)',
                    texname = '\\text{I145a33}')

I145a36 = Parameter(name = 'I145a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*yd3x3*complexconjugate(Rd3x3)*complexconjugate(yd3x3)',
                    texname = '\\text{I145a36}')

I145a63 = Parameter(name = 'I145a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*yd3x3*complexconjugate(Rd6x3)*complexconjugate(yd3x3)',
                    texname = '\\text{I145a63}')

I145a66 = Parameter(name = 'I145a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*yd3x3*complexconjugate(Rd6x3)*complexconjugate(yd3x3)',
                    texname = '\\text{I145a66}')

I146a33 = Parameter(name = 'I146a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*yu3x3*complexconjugate(Rd3x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I146a33}')

I146a36 = Parameter(name = 'I146a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*yu3x3*complexconjugate(Rd3x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I146a36}')

I146a63 = Parameter(name = 'I146a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*yu3x3*complexconjugate(Rd6x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I146a63}')

I146a66 = Parameter(name = 'I146a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*yu3x3*complexconjugate(Rd6x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I146a66}')

I147a33 = Parameter(name = 'I147a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x6*yu3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I147a33}')

I147a36 = Parameter(name = 'I147a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x6*yu3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I147a36}')

I147a63 = Parameter(name = 'I147a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x6*yu3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I147a63}')

I147a66 = Parameter(name = 'I147a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x6*yu3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I147a66}')

I148a11 = Parameter(name = 'I148a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I148a11}')

I148a22 = Parameter(name = 'I148a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I148a22}')

I148a33 = Parameter(name = 'I148a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I148a33}')

I148a36 = Parameter(name = 'I148a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I148a36}')

I148a63 = Parameter(name = 'I148a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I148a63}')

I148a66 = Parameter(name = 'I148a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I148a66}')

I149a33 = Parameter(name = 'I149a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I149a33}')

I149a36 = Parameter(name = 'I149a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I149a36}')

I149a44 = Parameter(name = 'I149a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I149a44}')

I149a55 = Parameter(name = 'I149a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I149a55}')

I149a63 = Parameter(name = 'I149a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I149a63}')

I149a66 = Parameter(name = 'I149a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I149a66}')

I15a11 = Parameter(name = 'I15a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I15a11}')

I15a22 = Parameter(name = 'I15a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I15a22}')

I15a33 = Parameter(name = 'I15a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I15a33}')

I15a36 = Parameter(name = 'I15a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I15a36}')

I15a63 = Parameter(name = 'I15a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I15a63}')

I15a66 = Parameter(name = 'I15a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I15a66}')

I150a33 = Parameter(name = 'I150a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I150a33}')

I150a36 = Parameter(name = 'I150a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I150a36}')

I150a63 = Parameter(name = 'I150a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I150a63}')

I150a66 = Parameter(name = 'I150a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I150a66}')

I151a33 = Parameter(name = 'I151a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru3x6)*complexconjugate(tu3x3)',
                    texname = '\\text{I151a33}')

I151a36 = Parameter(name = 'I151a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru3x6)*complexconjugate(tu3x3)',
                    texname = '\\text{I151a36}')

I151a63 = Parameter(name = 'I151a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru6x6)*complexconjugate(tu3x3)',
                    texname = '\\text{I151a63}')

I151a66 = Parameter(name = 'I151a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru6x6)*complexconjugate(tu3x3)',
                    texname = '\\text{I151a66}')

I152a33 = Parameter(name = 'I152a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*tu3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I152a33}')

I152a36 = Parameter(name = 'I152a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*tu3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I152a36}')

I152a63 = Parameter(name = 'I152a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*tu3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I152a63}')

I152a66 = Parameter(name = 'I152a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*tu3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I152a66}')

I153a33 = Parameter(name = 'I153a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*yu3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I153a33}')

I153a36 = Parameter(name = 'I153a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*yu3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I153a36}')

I153a63 = Parameter(name = 'I153a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*yu3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I153a63}')

I153a66 = Parameter(name = 'I153a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*yu3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I153a66}')

I154a33 = Parameter(name = 'I154a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*yu3x3*complexconjugate(Ru3x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I154a33}')

I154a36 = Parameter(name = 'I154a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*yu3x3*complexconjugate(Ru3x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I154a36}')

I154a63 = Parameter(name = 'I154a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*yu3x3*complexconjugate(Ru6x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I154a63}')

I154a66 = Parameter(name = 'I154a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*yu3x3*complexconjugate(Ru6x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I154a66}')

I155a33 = Parameter(name = 'I155a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*yu3x3*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I155a33}')

I155a36 = Parameter(name = 'I155a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*yu3x3*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I155a36}')

I155a63 = Parameter(name = 'I155a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*yu3x3*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I155a63}')

I155a66 = Parameter(name = 'I155a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*yu3x3*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I155a66}')

I156a33 = Parameter(name = 'I156a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*yu3x3*complexconjugate(Ru3x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I156a33}')

I156a36 = Parameter(name = 'I156a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*yu3x3*complexconjugate(Ru3x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I156a36}')

I156a63 = Parameter(name = 'I156a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*yu3x3*complexconjugate(Ru6x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I156a63}')

I156a66 = Parameter(name = 'I156a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*yu3x3*complexconjugate(Ru6x3)*complexconjugate(yu3x3)',
                    texname = '\\text{I156a66}')

I157a33 = Parameter(name = 'I157a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*yu3x3*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I157a33}')

I157a36 = Parameter(name = 'I157a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*yu3x3*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I157a36}')

I157a63 = Parameter(name = 'I157a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*yu3x3*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I157a63}')

I157a66 = Parameter(name = 'I157a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*yu3x3*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I157a66}')

I158a33 = Parameter(name = 'I158a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)*complexconjugate(yd3x3)',
                    texname = '\\text{I158a33}')

I158a36 = Parameter(name = 'I158a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)*complexconjugate(yd3x3)',
                    texname = '\\text{I158a36}')

I158a63 = Parameter(name = 'I158a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)*complexconjugate(yd3x3)',
                    texname = '\\text{I158a63}')

I158a66 = Parameter(name = 'I158a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)*complexconjugate(yd3x3)',
                    texname = '\\text{I158a66}')

I159a11 = Parameter(name = 'I159a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*complexconjugate(Rd1x1)',
                    texname = '\\text{I159a11}')

I159a22 = Parameter(name = 'I159a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I159a22}')

I159a33 = Parameter(name = 'I159a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*complexconjugate(Rd3x3)',
                    texname = '\\text{I159a33}')

I159a36 = Parameter(name = 'I159a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*complexconjugate(Rd3x3)',
                    texname = '\\text{I159a36}')

I159a63 = Parameter(name = 'I159a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*complexconjugate(Rd6x3)',
                    texname = '\\text{I159a63}')

I159a66 = Parameter(name = 'I159a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*complexconjugate(Rd6x3)',
                    texname = '\\text{I159a66}')

I16a33 = Parameter(name = 'I16a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*complexconjugate(Rd3x6)',
                   texname = '\\text{I16a33}')

I16a36 = Parameter(name = 'I16a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*complexconjugate(Rd3x6)',
                   texname = '\\text{I16a36}')

I16a44 = Parameter(name = 'I16a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I16a44}')

I16a55 = Parameter(name = 'I16a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*complexconjugate(Rd5x5)',
                   texname = '\\text{I16a55}')

I16a63 = Parameter(name = 'I16a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*complexconjugate(Rd6x6)',
                   texname = '\\text{I16a63}')

I16a66 = Parameter(name = 'I16a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*complexconjugate(Rd6x6)',
                   texname = '\\text{I16a66}')

I160a33 = Parameter(name = 'I160a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I160a33}')

I160a36 = Parameter(name = 'I160a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I160a36}')

I160a63 = Parameter(name = 'I160a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I160a63}')

I160a66 = Parameter(name = 'I160a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I160a66}')

I161a11 = Parameter(name = 'I161a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I161a11}')

I161a22 = Parameter(name = 'I161a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I161a22}')

I161a33 = Parameter(name = 'I161a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I161a33}')

I161a36 = Parameter(name = 'I161a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I161a36}')

I161a63 = Parameter(name = 'I161a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I161a63}')

I161a66 = Parameter(name = 'I161a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I161a66}')

I162a33 = Parameter(name = 'I162a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I162a33}')

I162a36 = Parameter(name = 'I162a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I162a36}')

I162a44 = Parameter(name = 'I162a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I162a44}')

I162a55 = Parameter(name = 'I162a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I162a55}')

I162a63 = Parameter(name = 'I162a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I162a63}')

I162a66 = Parameter(name = 'I162a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I162a66}')

I163a11 = Parameter(name = 'I163a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                    texname = '\\text{I163a11}')

I163a22 = Parameter(name = 'I163a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I163a22}')

I163a33 = Parameter(name = 'I163a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                    texname = '\\text{I163a33}')

I163a36 = Parameter(name = 'I163a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                    texname = '\\text{I163a36}')

I163a63 = Parameter(name = 'I163a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                    texname = '\\text{I163a63}')

I163a66 = Parameter(name = 'I163a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                    texname = '\\text{I163a66}')

I164a33 = Parameter(name = 'I164a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I164a33}')

I164a36 = Parameter(name = 'I164a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I164a36}')

I164a63 = Parameter(name = 'I164a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I164a63}')

I164a66 = Parameter(name = 'I164a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I164a66}')

I165a33 = Parameter(name = 'I165a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                    texname = '\\text{I165a33}')

I165a36 = Parameter(name = 'I165a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                    texname = '\\text{I165a36}')

I165a63 = Parameter(name = 'I165a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                    texname = '\\text{I165a63}')

I165a66 = Parameter(name = 'I165a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                    texname = '\\text{I165a66}')

I166a33 = Parameter(name = 'I166a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I166a33}')

I166a36 = Parameter(name = 'I166a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I166a36}')

I166a63 = Parameter(name = 'I166a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I166a63}')

I166a66 = Parameter(name = 'I166a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I166a66}')

I167a33 = Parameter(name = 'I167a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x6*yu3x3*complexconjugate(Rd3x3)',
                    texname = '\\text{I167a33}')

I167a36 = Parameter(name = 'I167a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x6*yu3x3*complexconjugate(Rd3x3)',
                    texname = '\\text{I167a36}')

I167a63 = Parameter(name = 'I167a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x6*yu3x3*complexconjugate(Rd6x3)',
                    texname = '\\text{I167a63}')

I167a66 = Parameter(name = 'I167a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x6*yu3x3*complexconjugate(Rd6x3)',
                    texname = '\\text{I167a66}')

I168a11 = Parameter(name = 'I168a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I168a11}')

I168a22 = Parameter(name = 'I168a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I168a22}')

I168a33 = Parameter(name = 'I168a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I168a33}')

I168a36 = Parameter(name = 'I168a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I168a36}')

I168a63 = Parameter(name = 'I168a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I168a63}')

I168a66 = Parameter(name = 'I168a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I168a66}')

I169a11 = Parameter(name = 'I169a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I169a11}')

I169a22 = Parameter(name = 'I169a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I169a22}')

I169a33 = Parameter(name = 'I169a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I169a33}')

I169a36 = Parameter(name = 'I169a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I169a36}')

I169a63 = Parameter(name = 'I169a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I169a63}')

I169a66 = Parameter(name = 'I169a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I169a66}')

I17a33 = Parameter(name = 'I17a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd3x6)*complexconjugate(td3x3)',
                   texname = '\\text{I17a33}')

I17a36 = Parameter(name = 'I17a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd3x6)*complexconjugate(td3x3)',
                   texname = '\\text{I17a36}')

I17a63 = Parameter(name = 'I17a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd6x6)*complexconjugate(td3x3)',
                   texname = '\\text{I17a63}')

I17a66 = Parameter(name = 'I17a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd6x6)*complexconjugate(td3x3)',
                   texname = '\\text{I17a66}')

I170a11 = Parameter(name = 'I170a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I170a11}')

I170a22 = Parameter(name = 'I170a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I170a22}')

I170a33 = Parameter(name = 'I170a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I170a33}')

I170a36 = Parameter(name = 'I170a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I170a36}')

I170a63 = Parameter(name = 'I170a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I170a63}')

I170a66 = Parameter(name = 'I170a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I170a66}')

I171a11 = Parameter(name = 'I171a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I171a11}')

I171a22 = Parameter(name = 'I171a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I171a22}')

I171a33 = Parameter(name = 'I171a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I171a33}')

I171a36 = Parameter(name = 'I171a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I171a36}')

I171a63 = Parameter(name = 'I171a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I171a63}')

I171a66 = Parameter(name = 'I171a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I171a66}')

I172a33 = Parameter(name = 'I172a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I172a33}')

I172a36 = Parameter(name = 'I172a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I172a36}')

I172a44 = Parameter(name = 'I172a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I172a44}')

I172a55 = Parameter(name = 'I172a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I172a55}')

I172a63 = Parameter(name = 'I172a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I172a63}')

I172a66 = Parameter(name = 'I172a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I172a66}')

I173a11 = Parameter(name = 'I173a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I173a11}')

I173a22 = Parameter(name = 'I173a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I173a22}')

I173a33 = Parameter(name = 'I173a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I173a33}')

I173a36 = Parameter(name = 'I173a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I173a36}')

I173a63 = Parameter(name = 'I173a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I173a63}')

I173a66 = Parameter(name = 'I173a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I173a66}')

I174a33 = Parameter(name = 'I174a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I174a33}')

I174a36 = Parameter(name = 'I174a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I174a36}')

I174a44 = Parameter(name = 'I174a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I174a44}')

I174a55 = Parameter(name = 'I174a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I174a55}')

I174a63 = Parameter(name = 'I174a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I174a63}')

I174a66 = Parameter(name = 'I174a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I174a66}')

I175a33 = Parameter(name = 'I175a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I175a33}')

I175a36 = Parameter(name = 'I175a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I175a36}')

I175a44 = Parameter(name = 'I175a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I175a44}')

I175a55 = Parameter(name = 'I175a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I175a55}')

I175a63 = Parameter(name = 'I175a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I175a63}')

I175a66 = Parameter(name = 'I175a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I175a66}')

I176a33 = Parameter(name = 'I176a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I176a33}')

I176a36 = Parameter(name = 'I176a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I176a36}')

I176a44 = Parameter(name = 'I176a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I176a44}')

I176a55 = Parameter(name = 'I176a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I176a55}')

I176a63 = Parameter(name = 'I176a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I176a63}')

I176a66 = Parameter(name = 'I176a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I176a66}')

I177a33 = Parameter(name = 'I177a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I177a33}')

I177a36 = Parameter(name = 'I177a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I177a36}')

I177a44 = Parameter(name = 'I177a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I177a44}')

I177a55 = Parameter(name = 'I177a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I177a55}')

I177a63 = Parameter(name = 'I177a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I177a63}')

I177a66 = Parameter(name = 'I177a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I177a66}')

I178a33 = Parameter(name = 'I178a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I178a33}')

I178a36 = Parameter(name = 'I178a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I178a36}')

I178a44 = Parameter(name = 'I178a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I178a44}')

I178a55 = Parameter(name = 'I178a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I178a55}')

I178a63 = Parameter(name = 'I178a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I178a63}')

I178a66 = Parameter(name = 'I178a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I178a66}')

I179a33 = Parameter(name = 'I179a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I179a33}')

I179a36 = Parameter(name = 'I179a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I179a36}')

I179a63 = Parameter(name = 'I179a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I179a63}')

I179a66 = Parameter(name = 'I179a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I179a66}')

I18a33 = Parameter(name = 'I18a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I18a33}')

I18a36 = Parameter(name = 'I18a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I18a36}')

I18a63 = Parameter(name = 'I18a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I18a63}')

I18a66 = Parameter(name = 'I18a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I18a66}')

I180a33 = Parameter(name = 'I180a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*yu3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I180a33}')

I180a36 = Parameter(name = 'I180a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*yu3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I180a36}')

I180a63 = Parameter(name = 'I180a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*yu3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I180a63}')

I180a66 = Parameter(name = 'I180a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*yu3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I180a66}')

I181a33 = Parameter(name = 'I181a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I181a33}')

I181a36 = Parameter(name = 'I181a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I181a36}')

I181a63 = Parameter(name = 'I181a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I181a63}')

I181a66 = Parameter(name = 'I181a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I181a66}')

I182a33 = Parameter(name = 'I182a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*yu3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I182a33}')

I182a36 = Parameter(name = 'I182a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*yu3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I182a36}')

I182a63 = Parameter(name = 'I182a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*yu3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I182a63}')

I182a66 = Parameter(name = 'I182a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*yu3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I182a66}')

I183a33 = Parameter(name = 'I183a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I183a33}')

I183a36 = Parameter(name = 'I183a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I183a36}')

I183a63 = Parameter(name = 'I183a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I183a63}')

I183a66 = Parameter(name = 'I183a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I183a66}')

I184a33 = Parameter(name = 'I184a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*yu3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I184a33}')

I184a36 = Parameter(name = 'I184a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*yu3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I184a36}')

I184a63 = Parameter(name = 'I184a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*yu3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I184a63}')

I184a66 = Parameter(name = 'I184a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*yu3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I184a66}')

I185a33 = Parameter(name = 'I185a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I185a33}')

I185a36 = Parameter(name = 'I185a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I185a36}')

I185a63 = Parameter(name = 'I185a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I185a63}')

I185a66 = Parameter(name = 'I185a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I185a66}')

I186a33 = Parameter(name = 'I186a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*yu3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I186a33}')

I186a36 = Parameter(name = 'I186a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*yu3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I186a36}')

I186a63 = Parameter(name = 'I186a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*yu3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I186a63}')

I186a66 = Parameter(name = 'I186a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*yu3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I186a66}')

I187a11 = Parameter(name = 'I187a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*complexconjugate(Rd1x1)',
                    texname = '\\text{I187a11}')

I187a22 = Parameter(name = 'I187a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I187a22}')

I187a33 = Parameter(name = 'I187a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*complexconjugate(Rd3x3)',
                    texname = '\\text{I187a33}')

I187a36 = Parameter(name = 'I187a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*complexconjugate(Rd6x3)',
                    texname = '\\text{I187a36}')

I188a33 = Parameter(name = 'I188a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I188a33}')

I188a36 = Parameter(name = 'I188a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                    texname = '\\text{I188a36}')

I189a33 = Parameter(name = 'I189a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*yu3x3*complexconjugate(Rd3x3)',
                    texname = '\\text{I189a33}')

I189a36 = Parameter(name = 'I189a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*yu3x3*complexconjugate(Rd6x3)',
                    texname = '\\text{I189a36}')

I19a33 = Parameter(name = 'I19a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*td3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I19a33}')

I19a36 = Parameter(name = 'I19a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*td3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I19a36}')

I19a63 = Parameter(name = 'I19a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*td3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I19a63}')

I19a66 = Parameter(name = 'I19a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*td3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I19a66}')

I190a11 = Parameter(name = 'I190a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Rl1x1)',
                    texname = '\\text{I190a11}')

I190a22 = Parameter(name = 'I190a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Rl2x2)',
                    texname = '\\text{I190a22}')

I190a33 = Parameter(name = 'I190a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Rl3x3)',
                    texname = '\\text{I190a33}')

I190a36 = Parameter(name = 'I190a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Rl6x3)',
                    texname = '\\text{I190a36}')

I191a33 = Parameter(name = 'I191a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                    texname = '\\text{I191a33}')

I191a36 = Parameter(name = 'I191a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                    texname = '\\text{I191a36}')

I192a11 = Parameter(name = 'I192a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Rn1x1)',
                    texname = '\\text{I192a11}')

I192a22 = Parameter(name = 'I192a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Rn2x2)',
                    texname = '\\text{I192a22}')

I192a33 = Parameter(name = 'I192a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Rn3x3)',
                    texname = '\\text{I192a33}')

I193a33 = Parameter(name = 'I193a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'ye3x3*complexconjugate(Rn3x3)',
                    texname = '\\text{I193a33}')

I194a11 = Parameter(name = 'I194a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                    texname = '\\text{I194a11}')

I194a22 = Parameter(name = 'I194a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I194a22}')

I194a33 = Parameter(name = 'I194a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                    texname = '\\text{I194a33}')

I194a36 = Parameter(name = 'I194a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                    texname = '\\text{I194a36}')

I195a33 = Parameter(name = 'I195a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(CKM3x3)*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I195a33}')

I195a36 = Parameter(name = 'I195a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(CKM3x3)*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                    texname = '\\text{I195a36}')

I196a33 = Parameter(name = 'I196a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                    texname = '\\text{I196a33}')

I196a36 = Parameter(name = 'I196a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                    texname = '\\text{I196a36}')

I197a11 = Parameter(name = 'I197a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*complexconjugate(Rd1x1)',
                    texname = '\\text{I197a11}')

I197a22 = Parameter(name = 'I197a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I197a22}')

I197a33 = Parameter(name = 'I197a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*complexconjugate(Rd3x3)',
                    texname = '\\text{I197a33}')

I197a36 = Parameter(name = 'I197a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*complexconjugate(Rd3x3)',
                    texname = '\\text{I197a36}')

I197a63 = Parameter(name = 'I197a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru3x3*complexconjugate(Rd6x3)',
                    texname = '\\text{I197a63}')

I197a66 = Parameter(name = 'I197a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM3x3*Ru6x3*complexconjugate(Rd6x3)',
                    texname = '\\text{I197a66}')

I198a11 = Parameter(name = 'I198a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*complexconjugate(Rl1x1)',
                    texname = '\\text{I198a11}')

I198a22 = Parameter(name = 'I198a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl2x2)',
                    texname = '\\text{I198a22}')

I198a33 = Parameter(name = 'I198a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*complexconjugate(Rl3x3)',
                    texname = '\\text{I198a33}')

I198a36 = Parameter(name = 'I198a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn3x3*complexconjugate(Rl6x3)',
                    texname = '\\text{I198a36}')

I199a11 = Parameter(name = 'I199a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                    texname = '\\text{I199a11}')

I199a22 = Parameter(name = 'I199a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I199a22}')

I199a33 = Parameter(name = 'I199a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                    texname = '\\text{I199a33}')

I199a36 = Parameter(name = 'I199a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                    texname = '\\text{I199a36}')

I199a63 = Parameter(name = 'I199a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                    texname = '\\text{I199a63}')

I199a66 = Parameter(name = 'I199a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                    texname = '\\text{I199a66}')

I2a33 = Parameter(name = 'I2a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yd3x3*complexconjugate(CKM3x3)',
                  texname = '\\text{I2a33}')

I20a33 = Parameter(name = 'I20a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*yd3x3*complexconjugate(Rd3x3)*complexconjugate(yd3x3)',
                   texname = '\\text{I20a33}')

I20a36 = Parameter(name = 'I20a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*yd3x3*complexconjugate(Rd3x3)*complexconjugate(yd3x3)',
                   texname = '\\text{I20a36}')

I20a63 = Parameter(name = 'I20a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*yd3x3*complexconjugate(Rd6x3)*complexconjugate(yd3x3)',
                   texname = '\\text{I20a63}')

I20a66 = Parameter(name = 'I20a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*yd3x3*complexconjugate(Rd6x3)*complexconjugate(yd3x3)',
                   texname = '\\text{I20a66}')

I200a11 = Parameter(name = 'I200a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rl1x1*complexconjugate(Rn1x1)',
                    texname = '\\text{I200a11}')

I200a22 = Parameter(name = 'I200a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rl2x2*complexconjugate(Rn2x2)',
                    texname = '\\text{I200a22}')

I200a33 = Parameter(name = 'I200a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rl3x3*complexconjugate(Rn3x3)',
                    texname = '\\text{I200a33}')

I200a36 = Parameter(name = 'I200a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rl6x3*complexconjugate(Rn3x3)',
                    texname = '\\text{I200a36}')

I201a11 = Parameter(name = 'I201a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I201a11}')

I201a22 = Parameter(name = 'I201a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I201a22}')

I201a33 = Parameter(name = 'I201a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I201a33}')

I201a36 = Parameter(name = 'I201a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru3x3)',
                    texname = '\\text{I201a36}')

I201a63 = Parameter(name = 'I201a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I201a63}')

I201a66 = Parameter(name = 'I201a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x3*complexconjugate(Ru6x3)',
                    texname = '\\text{I201a66}')

I202a33 = Parameter(name = 'I202a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'ye3x3',
                    texname = '\\text{I202a33}')

I203a33 = Parameter(name = 'I203a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I203a33}')

I203a36 = Parameter(name = 'I203a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I203a36}')

I203a44 = Parameter(name = 'I203a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I203a44}')

I203a55 = Parameter(name = 'I203a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I203a55}')

I203a63 = Parameter(name = 'I203a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I203a63}')

I203a66 = Parameter(name = 'I203a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I203a66}')

I21a33 = Parameter(name = 'I21a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I21a33}')

I21a36 = Parameter(name = 'I21a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I21a36}')

I21a63 = Parameter(name = 'I21a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I21a63}')

I21a66 = Parameter(name = 'I21a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I21a66}')

I22a33 = Parameter(name = 'I22a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I22a33}')

I22a36 = Parameter(name = 'I22a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I22a36}')

I22a63 = Parameter(name = 'I22a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I22a63}')

I22a66 = Parameter(name = 'I22a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I22a66}')

I23a33 = Parameter(name = 'I23a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*yd3x3*complexconjugate(Rd3x3)*complexconjugate(yd3x3)',
                   texname = '\\text{I23a33}')

I23a36 = Parameter(name = 'I23a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*yd3x3*complexconjugate(Rd3x3)*complexconjugate(yd3x3)',
                   texname = '\\text{I23a36}')

I23a63 = Parameter(name = 'I23a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*yd3x3*complexconjugate(Rd6x3)*complexconjugate(yd3x3)',
                   texname = '\\text{I23a63}')

I23a66 = Parameter(name = 'I23a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*yd3x3*complexconjugate(Rd6x3)*complexconjugate(yd3x3)',
                   texname = '\\text{I23a66}')

I24a33 = Parameter(name = 'I24a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I24a33}')

I24a36 = Parameter(name = 'I24a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I24a36}')

I24a63 = Parameter(name = 'I24a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I24a63}')

I24a66 = Parameter(name = 'I24a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I24a66}')

I25a33 = Parameter(name = 'I25a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Rd3x3*yu3x3*complexconjugate(CKM3x3)*complexconjugate(Rd3x3)*complexconjugate(yu3x3)',
                   texname = '\\text{I25a33}')

I25a36 = Parameter(name = 'I25a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Rd6x3*yu3x3*complexconjugate(CKM3x3)*complexconjugate(Rd3x3)*complexconjugate(yu3x3)',
                   texname = '\\text{I25a36}')

I25a63 = Parameter(name = 'I25a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Rd3x3*yu3x3*complexconjugate(CKM3x3)*complexconjugate(Rd6x3)*complexconjugate(yu3x3)',
                   texname = '\\text{I25a63}')

I25a66 = Parameter(name = 'I25a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Rd6x3*yu3x3*complexconjugate(CKM3x3)*complexconjugate(Rd6x3)*complexconjugate(yu3x3)',
                   texname = '\\text{I25a66}')

I26a11 = Parameter(name = 'I26a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I26a11}')

I26a22 = Parameter(name = 'I26a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I26a22}')

I26a33 = Parameter(name = 'I26a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I26a33}')

I26a36 = Parameter(name = 'I26a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I26a36}')

I26a63 = Parameter(name = 'I26a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I26a63}')

I26a66 = Parameter(name = 'I26a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I26a66}')

I27a11 = Parameter(name = 'I27a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I27a11}')

I27a22 = Parameter(name = 'I27a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I27a22}')

I27a33 = Parameter(name = 'I27a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I27a33}')

I27a36 = Parameter(name = 'I27a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I27a36}')

I27a63 = Parameter(name = 'I27a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I27a63}')

I27a66 = Parameter(name = 'I27a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I27a66}')

I28a11 = Parameter(name = 'I28a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I28a11}')

I28a22 = Parameter(name = 'I28a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I28a22}')

I28a33 = Parameter(name = 'I28a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I28a33}')

I28a36 = Parameter(name = 'I28a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I28a36}')

I28a63 = Parameter(name = 'I28a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I28a63}')

I28a66 = Parameter(name = 'I28a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I28a66}')

I29a11 = Parameter(name = 'I29a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I29a11}')

I29a22 = Parameter(name = 'I29a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I29a22}')

I29a33 = Parameter(name = 'I29a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I29a33}')

I29a36 = Parameter(name = 'I29a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I29a36}')

I29a63 = Parameter(name = 'I29a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I29a63}')

I29a66 = Parameter(name = 'I29a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I29a66}')

I3a33 = Parameter(name = 'I3a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM3x3*complexconjugate(yd3x3)',
                  texname = '\\text{I3a33}')

I30a33 = Parameter(name = 'I30a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*complexconjugate(Rd3x6)',
                   texname = '\\text{I30a33}')

I30a36 = Parameter(name = 'I30a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*complexconjugate(Rd3x6)',
                   texname = '\\text{I30a36}')

I30a44 = Parameter(name = 'I30a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I30a44}')

I30a55 = Parameter(name = 'I30a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*complexconjugate(Rd5x5)',
                   texname = '\\text{I30a55}')

I30a63 = Parameter(name = 'I30a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*complexconjugate(Rd6x6)',
                   texname = '\\text{I30a63}')

I30a66 = Parameter(name = 'I30a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*complexconjugate(Rd6x6)',
                   texname = '\\text{I30a66}')

I31a11 = Parameter(name = 'I31a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I31a11}')

I31a22 = Parameter(name = 'I31a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I31a22}')

I31a33 = Parameter(name = 'I31a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I31a33}')

I31a36 = Parameter(name = 'I31a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I31a36}')

I31a63 = Parameter(name = 'I31a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I31a63}')

I31a66 = Parameter(name = 'I31a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I31a66}')

I32a33 = Parameter(name = 'I32a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*complexconjugate(Rd3x6)',
                   texname = '\\text{I32a33}')

I32a36 = Parameter(name = 'I32a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*complexconjugate(Rd3x6)',
                   texname = '\\text{I32a36}')

I32a44 = Parameter(name = 'I32a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I32a44}')

I32a55 = Parameter(name = 'I32a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*complexconjugate(Rd5x5)',
                   texname = '\\text{I32a55}')

I32a63 = Parameter(name = 'I32a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*complexconjugate(Rd6x6)',
                   texname = '\\text{I32a63}')

I32a66 = Parameter(name = 'I32a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*complexconjugate(Rd6x6)',
                   texname = '\\text{I32a66}')

I33a33 = Parameter(name = 'I33a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*complexconjugate(Rd3x6)',
                   texname = '\\text{I33a33}')

I33a36 = Parameter(name = 'I33a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*complexconjugate(Rd3x6)',
                   texname = '\\text{I33a36}')

I33a44 = Parameter(name = 'I33a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I33a44}')

I33a55 = Parameter(name = 'I33a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*complexconjugate(Rd5x5)',
                   texname = '\\text{I33a55}')

I33a63 = Parameter(name = 'I33a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*complexconjugate(Rd6x6)',
                   texname = '\\text{I33a63}')

I33a66 = Parameter(name = 'I33a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*complexconjugate(Rd6x6)',
                   texname = '\\text{I33a66}')

I34a33 = Parameter(name = 'I34a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*complexconjugate(Rd3x6)',
                   texname = '\\text{I34a33}')

I34a36 = Parameter(name = 'I34a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*complexconjugate(Rd3x6)',
                   texname = '\\text{I34a36}')

I34a44 = Parameter(name = 'I34a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I34a44}')

I34a55 = Parameter(name = 'I34a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*complexconjugate(Rd5x5)',
                   texname = '\\text{I34a55}')

I34a63 = Parameter(name = 'I34a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*complexconjugate(Rd6x6)',
                   texname = '\\text{I34a63}')

I34a66 = Parameter(name = 'I34a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*complexconjugate(Rd6x6)',
                   texname = '\\text{I34a66}')

I35a33 = Parameter(name = 'I35a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*complexconjugate(Rd3x6)',
                   texname = '\\text{I35a33}')

I35a36 = Parameter(name = 'I35a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*complexconjugate(Rd3x6)',
                   texname = '\\text{I35a36}')

I35a44 = Parameter(name = 'I35a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I35a44}')

I35a55 = Parameter(name = 'I35a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*complexconjugate(Rd5x5)',
                   texname = '\\text{I35a55}')

I35a63 = Parameter(name = 'I35a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*complexconjugate(Rd6x6)',
                   texname = '\\text{I35a63}')

I35a66 = Parameter(name = 'I35a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*complexconjugate(Rd6x6)',
                   texname = '\\text{I35a66}')

I36a33 = Parameter(name = 'I36a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*complexconjugate(Rd3x6)',
                   texname = '\\text{I36a33}')

I36a36 = Parameter(name = 'I36a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*complexconjugate(Rd3x6)',
                   texname = '\\text{I36a36}')

I36a44 = Parameter(name = 'I36a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I36a44}')

I36a55 = Parameter(name = 'I36a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*complexconjugate(Rd5x5)',
                   texname = '\\text{I36a55}')

I36a63 = Parameter(name = 'I36a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*complexconjugate(Rd6x6)',
                   texname = '\\text{I36a63}')

I36a66 = Parameter(name = 'I36a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*complexconjugate(Rd6x6)',
                   texname = '\\text{I36a66}')

I37a33 = Parameter(name = 'I37a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I37a33}')

I37a36 = Parameter(name = 'I37a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I37a36}')

I37a63 = Parameter(name = 'I37a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I37a63}')

I37a66 = Parameter(name = 'I37a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I37a66}')

I38a33 = Parameter(name = 'I38a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I38a33}')

I38a36 = Parameter(name = 'I38a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I38a36}')

I38a63 = Parameter(name = 'I38a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I38a63}')

I38a66 = Parameter(name = 'I38a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I38a66}')

I39a33 = Parameter(name = 'I39a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I39a33}')

I39a36 = Parameter(name = 'I39a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I39a36}')

I39a63 = Parameter(name = 'I39a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I39a63}')

I39a66 = Parameter(name = 'I39a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I39a66}')

I4a33 = Parameter(name = 'I4a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM3x3*yu3x3',
                  texname = '\\text{I4a33}')

I40a33 = Parameter(name = 'I40a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I40a33}')

I40a36 = Parameter(name = 'I40a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I40a36}')

I40a63 = Parameter(name = 'I40a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I40a63}')

I40a66 = Parameter(name = 'I40a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I40a66}')

I41a33 = Parameter(name = 'I41a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I41a33}')

I41a36 = Parameter(name = 'I41a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I41a36}')

I41a63 = Parameter(name = 'I41a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I41a63}')

I41a66 = Parameter(name = 'I41a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I41a66}')

I42a33 = Parameter(name = 'I42a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I42a33}')

I42a36 = Parameter(name = 'I42a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I42a36}')

I42a63 = Parameter(name = 'I42a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I42a63}')

I42a66 = Parameter(name = 'I42a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I42a66}')

I43a33 = Parameter(name = 'I43a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I43a33}')

I43a36 = Parameter(name = 'I43a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I43a36}')

I43a63 = Parameter(name = 'I43a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I43a63}')

I43a66 = Parameter(name = 'I43a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I43a66}')

I44a33 = Parameter(name = 'I44a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I44a33}')

I44a36 = Parameter(name = 'I44a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I44a36}')

I44a63 = Parameter(name = 'I44a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I44a63}')

I44a66 = Parameter(name = 'I44a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I44a66}')

I45a33 = Parameter(name = 'I45a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I45a33}')

I45a36 = Parameter(name = 'I45a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I45a36}')

I46a33 = Parameter(name = 'I46a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'ye3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I46a33}')

I46a36 = Parameter(name = 'I46a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'ye3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I46a36}')

I47a11 = Parameter(name = 'I47a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I47a11}')

I47a22 = Parameter(name = 'I47a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I47a22}')

I47a33 = Parameter(name = 'I47a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I47a33}')

I47a36 = Parameter(name = 'I47a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I47a36}')

I47a63 = Parameter(name = 'I47a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I47a63}')

I47a66 = Parameter(name = 'I47a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I47a66}')

I48a33 = Parameter(name = 'I48a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I48a33}')

I48a36 = Parameter(name = 'I48a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I48a36}')

I48a44 = Parameter(name = 'I48a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I48a44}')

I48a55 = Parameter(name = 'I48a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I48a55}')

I48a63 = Parameter(name = 'I48a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I48a63}')

I48a66 = Parameter(name = 'I48a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I48a66}')

I49a33 = Parameter(name = 'I49a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(ye3x3)',
                   texname = '\\text{I49a33}')

I49a36 = Parameter(name = 'I49a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(ye3x3)',
                   texname = '\\text{I49a36}')

I5a33 = Parameter(name = 'I5a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complexconjugate(ye3x3)',
                  texname = '\\text{I5a33}')

I50a33 = Parameter(name = 'I50a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3',
                   texname = '\\text{I50a33}')

I50a36 = Parameter(name = 'I50a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3',
                   texname = '\\text{I50a36}')

I51a11 = Parameter(name = 'I51a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1',
                   texname = '\\text{I51a11}')

I51a22 = Parameter(name = 'I51a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2',
                   texname = '\\text{I51a22}')

I51a33 = Parameter(name = 'I51a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3',
                   texname = '\\text{I51a33}')

I51a36 = Parameter(name = 'I51a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3',
                   texname = '\\text{I51a36}')

I52a33 = Parameter(name = 'I52a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3',
                   texname = '\\text{I52a33}')

I52a36 = Parameter(name = 'I52a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3',
                   texname = '\\text{I52a36}')

I53a11 = Parameter(name = 'I53a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I53a11}')

I53a22 = Parameter(name = 'I53a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I53a22}')

I53a33 = Parameter(name = 'I53a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I53a33}')

I53a36 = Parameter(name = 'I53a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I53a36}')

I53a63 = Parameter(name = 'I53a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I53a63}')

I53a66 = Parameter(name = 'I53a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I53a66}')

I54a33 = Parameter(name = 'I54a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I54a33}')

I54a36 = Parameter(name = 'I54a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I54a36}')

I54a44 = Parameter(name = 'I54a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I54a44}')

I54a55 = Parameter(name = 'I54a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I54a55}')

I54a63 = Parameter(name = 'I54a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I54a63}')

I54a66 = Parameter(name = 'I54a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I54a66}')

I55a33 = Parameter(name = 'I55a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x6)*complexconjugate(te3x3)',
                   texname = '\\text{I55a33}')

I55a36 = Parameter(name = 'I55a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x6)*complexconjugate(te3x3)',
                   texname = '\\text{I55a36}')

I55a63 = Parameter(name = 'I55a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x6)*complexconjugate(te3x3)',
                   texname = '\\text{I55a63}')

I55a66 = Parameter(name = 'I55a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x6)*complexconjugate(te3x3)',
                   texname = '\\text{I55a66}')

I56a33 = Parameter(name = 'I56a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I56a33}')

I56a36 = Parameter(name = 'I56a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I56a36}')

I56a63 = Parameter(name = 'I56a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I56a63}')

I56a66 = Parameter(name = 'I56a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I56a66}')

I57a33 = Parameter(name = 'I57a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*te3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I57a33}')

I57a36 = Parameter(name = 'I57a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*te3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I57a36}')

I57a63 = Parameter(name = 'I57a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*te3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I57a63}')

I57a66 = Parameter(name = 'I57a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*te3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I57a66}')

I58a33 = Parameter(name = 'I58a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*ye3x3*complexconjugate(Rl3x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I58a33}')

I58a36 = Parameter(name = 'I58a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*ye3x3*complexconjugate(Rl3x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I58a36}')

I58a63 = Parameter(name = 'I58a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*ye3x3*complexconjugate(Rl6x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I58a63}')

I58a66 = Parameter(name = 'I58a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*ye3x3*complexconjugate(Rl6x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I58a66}')

I59a33 = Parameter(name = 'I59a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I59a33}')

I59a36 = Parameter(name = 'I59a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I59a36}')

I59a63 = Parameter(name = 'I59a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I59a63}')

I59a66 = Parameter(name = 'I59a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I59a66}')

I6a33 = Parameter(name = 'I6a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                  texname = '\\text{I6a33}')

I6a36 = Parameter(name = 'I6a36',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                  texname = '\\text{I6a36}')

I60a33 = Parameter(name = 'I60a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I60a33}')

I60a36 = Parameter(name = 'I60a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I60a36}')

I60a63 = Parameter(name = 'I60a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I60a63}')

I60a66 = Parameter(name = 'I60a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I60a66}')

I61a33 = Parameter(name = 'I61a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*ye3x3*complexconjugate(Rl3x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I61a33}')

I61a36 = Parameter(name = 'I61a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*ye3x3*complexconjugate(Rl3x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I61a36}')

I61a63 = Parameter(name = 'I61a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*ye3x3*complexconjugate(Rl6x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I61a63}')

I61a66 = Parameter(name = 'I61a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*ye3x3*complexconjugate(Rl6x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I61a66}')

I62a33 = Parameter(name = 'I62a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I62a33}')

I62a36 = Parameter(name = 'I62a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I62a36}')

I62a63 = Parameter(name = 'I62a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I62a63}')

I62a66 = Parameter(name = 'I62a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I62a66}')

I63a11 = Parameter(name = 'I63a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I63a11}')

I63a22 = Parameter(name = 'I63a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I63a22}')

I63a33 = Parameter(name = 'I63a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I63a33}')

I63a36 = Parameter(name = 'I63a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I63a36}')

I63a63 = Parameter(name = 'I63a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I63a63}')

I63a66 = Parameter(name = 'I63a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I63a66}')

I64a11 = Parameter(name = 'I64a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I64a11}')

I64a22 = Parameter(name = 'I64a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I64a22}')

I64a33 = Parameter(name = 'I64a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I64a33}')

I64a36 = Parameter(name = 'I64a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I64a36}')

I64a63 = Parameter(name = 'I64a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I64a63}')

I64a66 = Parameter(name = 'I64a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I64a66}')

I65a33 = Parameter(name = 'I65a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*complexconjugate(Rd3x6)',
                   texname = '\\text{I65a33}')

I65a36 = Parameter(name = 'I65a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*complexconjugate(Rd3x6)',
                   texname = '\\text{I65a36}')

I65a44 = Parameter(name = 'I65a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I65a44}')

I65a55 = Parameter(name = 'I65a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*complexconjugate(Rd5x5)',
                   texname = '\\text{I65a55}')

I65a63 = Parameter(name = 'I65a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*complexconjugate(Rd6x6)',
                   texname = '\\text{I65a63}')

I65a66 = Parameter(name = 'I65a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*complexconjugate(Rd6x6)',
                   texname = '\\text{I65a66}')

I66a33 = Parameter(name = 'I66a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I66a33}')

I66a36 = Parameter(name = 'I66a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I66a36}')

I66a44 = Parameter(name = 'I66a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I66a44}')

I66a55 = Parameter(name = 'I66a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I66a55}')

I66a63 = Parameter(name = 'I66a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I66a63}')

I66a66 = Parameter(name = 'I66a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I66a66}')

I67a33 = Parameter(name = 'I67a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I67a33}')

I67a36 = Parameter(name = 'I67a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I67a36}')

I67a63 = Parameter(name = 'I67a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I67a63}')

I67a66 = Parameter(name = 'I67a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I67a66}')

I68a33 = Parameter(name = 'I68a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I68a33}')

I68a36 = Parameter(name = 'I68a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I68a36}')

I68a63 = Parameter(name = 'I68a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I68a63}')

I68a66 = Parameter(name = 'I68a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I68a66}')

I69a33 = Parameter(name = 'I69a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I69a33}')

I69a36 = Parameter(name = 'I69a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I69a36}')

I69a63 = Parameter(name = 'I69a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I69a63}')

I69a66 = Parameter(name = 'I69a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I69a66}')

I7a33 = Parameter(name = 'I7a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yd3x3*complexconjugate(Rd3x3)',
                  texname = '\\text{I7a33}')

I7a36 = Parameter(name = 'I7a36',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yd3x3*complexconjugate(Rd6x3)',
                  texname = '\\text{I7a36}')

I70a33 = Parameter(name = 'I70a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I70a33}')

I70a36 = Parameter(name = 'I70a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I70a36}')

I70a63 = Parameter(name = 'I70a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I70a63}')

I70a66 = Parameter(name = 'I70a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I70a66}')

I71a11 = Parameter(name = 'I71a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I71a11}')

I71a22 = Parameter(name = 'I71a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I71a22}')

I71a33 = Parameter(name = 'I71a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I71a33}')

I71a36 = Parameter(name = 'I71a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I71a36}')

I71a63 = Parameter(name = 'I71a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I71a63}')

I71a66 = Parameter(name = 'I71a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I71a66}')

I72a11 = Parameter(name = 'I72a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I72a11}')

I72a22 = Parameter(name = 'I72a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I72a22}')

I72a33 = Parameter(name = 'I72a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I72a33}')

I72a36 = Parameter(name = 'I72a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I72a36}')

I72a63 = Parameter(name = 'I72a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I72a63}')

I72a66 = Parameter(name = 'I72a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I72a66}')

I73a11 = Parameter(name = 'I73a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I73a11}')

I73a22 = Parameter(name = 'I73a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I73a22}')

I73a33 = Parameter(name = 'I73a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I73a33}')

I73a36 = Parameter(name = 'I73a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I73a36}')

I73a63 = Parameter(name = 'I73a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I73a63}')

I73a66 = Parameter(name = 'I73a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I73a66}')

I74a11 = Parameter(name = 'I74a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I74a11}')

I74a22 = Parameter(name = 'I74a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I74a22}')

I74a33 = Parameter(name = 'I74a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I74a33}')

I74a36 = Parameter(name = 'I74a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I74a36}')

I74a63 = Parameter(name = 'I74a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I74a63}')

I74a66 = Parameter(name = 'I74a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I74a66}')

I75a11 = Parameter(name = 'I75a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I75a11}')

I75a22 = Parameter(name = 'I75a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I75a22}')

I75a33 = Parameter(name = 'I75a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I75a33}')

I75a36 = Parameter(name = 'I75a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I75a36}')

I75a63 = Parameter(name = 'I75a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I75a63}')

I75a66 = Parameter(name = 'I75a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I75a66}')

I76a33 = Parameter(name = 'I76a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I76a33}')

I76a36 = Parameter(name = 'I76a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I76a36}')

I76a44 = Parameter(name = 'I76a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I76a44}')

I76a55 = Parameter(name = 'I76a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I76a55}')

I76a63 = Parameter(name = 'I76a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I76a63}')

I76a66 = Parameter(name = 'I76a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I76a66}')

I77a33 = Parameter(name = 'I77a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I77a33}')

I77a36 = Parameter(name = 'I77a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I77a36}')

I77a44 = Parameter(name = 'I77a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I77a44}')

I77a55 = Parameter(name = 'I77a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I77a55}')

I77a63 = Parameter(name = 'I77a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I77a63}')

I77a66 = Parameter(name = 'I77a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I77a66}')

I78a33 = Parameter(name = 'I78a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I78a33}')

I78a36 = Parameter(name = 'I78a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I78a36}')

I78a44 = Parameter(name = 'I78a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I78a44}')

I78a55 = Parameter(name = 'I78a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I78a55}')

I78a63 = Parameter(name = 'I78a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I78a63}')

I78a66 = Parameter(name = 'I78a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I78a66}')

I79a33 = Parameter(name = 'I79a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I79a33}')

I79a36 = Parameter(name = 'I79a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I79a36}')

I79a44 = Parameter(name = 'I79a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I79a44}')

I79a55 = Parameter(name = 'I79a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I79a55}')

I79a63 = Parameter(name = 'I79a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I79a63}')

I79a66 = Parameter(name = 'I79a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I79a66}')

I8a11 = Parameter(name = 'I8a11',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd1x1*complexconjugate(Rd1x1)',
                  texname = '\\text{I8a11}')

I8a22 = Parameter(name = 'I8a22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd2x2*complexconjugate(Rd2x2)',
                  texname = '\\text{I8a22}')

I8a33 = Parameter(name = 'I8a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd3x3*complexconjugate(Rd3x3)',
                  texname = '\\text{I8a33}')

I8a36 = Parameter(name = 'I8a36',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd6x3*complexconjugate(Rd3x3)',
                  texname = '\\text{I8a36}')

I8a63 = Parameter(name = 'I8a63',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd3x3*complexconjugate(Rd6x3)',
                  texname = '\\text{I8a63}')

I8a66 = Parameter(name = 'I8a66',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd6x3*complexconjugate(Rd6x3)',
                  texname = '\\text{I8a66}')

I80a33 = Parameter(name = 'I80a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I80a33}')

I80a36 = Parameter(name = 'I80a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I80a36}')

I80a44 = Parameter(name = 'I80a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I80a44}')

I80a55 = Parameter(name = 'I80a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I80a55}')

I80a63 = Parameter(name = 'I80a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I80a63}')

I80a66 = Parameter(name = 'I80a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I80a66}')

I81a33 = Parameter(name = 'I81a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I81a33}')

I81a36 = Parameter(name = 'I81a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I81a36}')

I81a44 = Parameter(name = 'I81a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I81a44}')

I81a55 = Parameter(name = 'I81a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I81a55}')

I81a63 = Parameter(name = 'I81a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I81a63}')

I81a66 = Parameter(name = 'I81a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I81a66}')

I82a33 = Parameter(name = 'I82a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I82a33}')

I82a36 = Parameter(name = 'I82a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I82a36}')

I82a63 = Parameter(name = 'I82a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I82a63}')

I82a66 = Parameter(name = 'I82a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I82a66}')

I83a33 = Parameter(name = 'I83a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I83a33}')

I83a36 = Parameter(name = 'I83a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I83a36}')

I83a63 = Parameter(name = 'I83a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I83a63}')

I83a66 = Parameter(name = 'I83a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I83a66}')

I84a33 = Parameter(name = 'I84a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I84a33}')

I84a36 = Parameter(name = 'I84a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I84a36}')

I84a63 = Parameter(name = 'I84a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I84a63}')

I84a66 = Parameter(name = 'I84a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I84a66}')

I85a33 = Parameter(name = 'I85a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I85a33}')

I85a36 = Parameter(name = 'I85a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I85a36}')

I85a63 = Parameter(name = 'I85a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I85a63}')

I85a66 = Parameter(name = 'I85a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I85a66}')

I86a33 = Parameter(name = 'I86a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I86a33}')

I86a36 = Parameter(name = 'I86a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I86a36}')

I86a63 = Parameter(name = 'I86a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I86a63}')

I86a66 = Parameter(name = 'I86a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I86a66}')

I87a33 = Parameter(name = 'I87a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I87a33}')

I87a36 = Parameter(name = 'I87a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I87a36}')

I87a63 = Parameter(name = 'I87a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I87a63}')

I87a66 = Parameter(name = 'I87a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I87a66}')

I88a33 = Parameter(name = 'I88a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I88a33}')

I88a36 = Parameter(name = 'I88a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I88a36}')

I88a63 = Parameter(name = 'I88a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I88a63}')

I88a66 = Parameter(name = 'I88a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I88a66}')

I89a33 = Parameter(name = 'I89a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I89a33}')

I89a36 = Parameter(name = 'I89a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I89a36}')

I89a63 = Parameter(name = 'I89a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I89a63}')

I89a66 = Parameter(name = 'I89a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I89a66}')

I9a33 = Parameter(name = 'I9a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd3x6*complexconjugate(Rd3x6)',
                  texname = '\\text{I9a33}')

I9a36 = Parameter(name = 'I9a36',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd6x6*complexconjugate(Rd3x6)',
                  texname = '\\text{I9a36}')

I9a44 = Parameter(name = 'I9a44',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd4x4*complexconjugate(Rd4x4)',
                  texname = '\\text{I9a44}')

I9a55 = Parameter(name = 'I9a55',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd5x5*complexconjugate(Rd5x5)',
                  texname = '\\text{I9a55}')

I9a63 = Parameter(name = 'I9a63',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd3x6*complexconjugate(Rd6x6)',
                  texname = '\\text{I9a63}')

I9a66 = Parameter(name = 'I9a66',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd6x6*complexconjugate(Rd6x6)',
                  texname = '\\text{I9a66}')

I90a11 = Parameter(name = 'I90a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rn1x1)',
                   texname = '\\text{I90a11}')

I90a22 = Parameter(name = 'I90a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rn2x2)',
                   texname = '\\text{I90a22}')

I90a33 = Parameter(name = 'I90a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rn3x3)',
                   texname = '\\text{I90a33}')

I90a36 = Parameter(name = 'I90a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rn3x3)',
                   texname = '\\text{I90a36}')

I91a33 = Parameter(name = 'I91a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*te3x3*complexconjugate(Rn3x3)',
                   texname = '\\text{I91a33}')

I91a36 = Parameter(name = 'I91a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*te3x3*complexconjugate(Rn3x3)',
                   texname = '\\text{I91a36}')

I92a33 = Parameter(name = 'I92a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*ye3x3*complexconjugate(Rn3x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I92a33}')

I92a36 = Parameter(name = 'I92a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*ye3x3*complexconjugate(Rn3x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I92a36}')

I93a33 = Parameter(name = 'I93a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rn3x3)',
                   texname = '\\text{I93a33}')

I93a36 = Parameter(name = 'I93a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rn3x3)',
                   texname = '\\text{I93a36}')

I94a11 = Parameter(name = 'I94a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rn1x1)',
                   texname = '\\text{I94a11}')

I94a22 = Parameter(name = 'I94a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rn2x2)',
                   texname = '\\text{I94a22}')

I94a33 = Parameter(name = 'I94a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rn3x3)',
                   texname = '\\text{I94a33}')

I94a36 = Parameter(name = 'I94a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rn3x3)',
                   texname = '\\text{I94a36}')

I95a33 = Parameter(name = 'I95a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*ye3x3*complexconjugate(Rn3x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I95a33}')

I95a36 = Parameter(name = 'I95a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*ye3x3*complexconjugate(Rn3x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I95a36}')

I96a33 = Parameter(name = 'I96a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn3x3*ye3x3*complexconjugate(Rn3x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I96a33}')

I97a11 = Parameter(name = 'I97a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I97a11}')

I97a22 = Parameter(name = 'I97a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I97a22}')

I97a33 = Parameter(name = 'I97a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I97a33}')

I97a36 = Parameter(name = 'I97a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I97a36}')

I97a63 = Parameter(name = 'I97a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I97a63}')

I97a66 = Parameter(name = 'I97a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I97a66}')

I98a33 = Parameter(name = 'I98a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I98a33}')

I98a36 = Parameter(name = 'I98a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I98a36}')

I98a44 = Parameter(name = 'I98a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I98a44}')

I98a55 = Parameter(name = 'I98a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I98a55}')

I98a63 = Parameter(name = 'I98a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I98a63}')

I98a66 = Parameter(name = 'I98a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I98a66}')

I99a11 = Parameter(name = 'I99a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rn1x1)',
                   texname = '\\text{I99a11}')

I99a22 = Parameter(name = 'I99a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rn2x2)',
                   texname = '\\text{I99a22}')

I99a33 = Parameter(name = 'I99a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rn3x3)',
                   texname = '\\text{I99a33}')

I99a36 = Parameter(name = 'I99a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rn3x3)',
                   texname = '\\text{I99a36}')

