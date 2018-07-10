# This file was automatically created by FeynRules 2.3.21
# Mathematica version: 10.3.0 for Mac OS X x86 (64-bit) (October 9, 2015)
# Date: Sun 1 May 2016 10:59:45



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
                   value = 0.707106781,
                   texname = '\\text{RRd1x1}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 1, 1 ])

RRd1x4 = Parameter(name = 'RRd1x4',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRd1x4}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 1, 4 ])

RRd2x2 = Parameter(name = 'RRd2x2',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRd2x2}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 2, 2 ])

RRd2x5 = Parameter(name = 'RRd2x5',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRd2x5}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 2, 5 ])

RRd3x3 = Parameter(name = 'RRd3x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRd3x3}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 3, 3 ])

RRd3x6 = Parameter(name = 'RRd3x6',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRd3x6}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 3, 6 ])

RRd4x1 = Parameter(name = 'RRd4x1',
                   nature = 'external',
                   type = 'real',
                   value = -0.707106781,
                   texname = '\\text{RRd4x1}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 4, 1 ])

RRd4x4 = Parameter(name = 'RRd4x4',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRd4x4}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 4, 4 ])

RRd5x2 = Parameter(name = 'RRd5x2',
                   nature = 'external',
                   type = 'real',
                   value = -0.707106781,
                   texname = '\\text{RRd5x2}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 5, 2 ])

RRd5x5 = Parameter(name = 'RRd5x5',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRd5x5}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 5, 5 ])

RRd6x3 = Parameter(name = 'RRd6x3',
                   nature = 'external',
                   type = 'real',
                   value = -0.707106781,
                   texname = '\\text{RRd6x3}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 6, 3 ])

RRd6x6 = Parameter(name = 'RRd6x6',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
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
                texname = '\\text{Subsuperscript}[m,A,2]',
                lhablock = 'HMIX',
                lhacode = [ 4 ])

RmD21x1 = Parameter(name = 'RmD21x1',
                    nature = 'external',
                    type = 'real',
                    value = 276684.674,
                    texname = '\\text{RmD21x1}',
                    lhablock = 'MSD2',
                    lhacode = [ 1, 1 ])

RmD22x2 = Parameter(name = 'RmD22x2',
                    nature = 'external',
                    type = 'real',
                    value = 275684.674,
                    texname = '\\text{RmD22x2}',
                    lhablock = 'MSD2',
                    lhacode = [ 2, 2 ])

RmD23x3 = Parameter(name = 'RmD23x3',
                    nature = 'external',
                    type = 'real',
                    value = 274261.969,
                    texname = '\\text{RmD23x3}',
                    lhablock = 'MSD2',
                    lhacode = [ 3, 3 ])

RmE21x1 = Parameter(name = 'RmE21x1',
                    nature = 'external',
                    type = 'real',
                    value = 18930.6287,
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
                    value = 38455.67,
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
                 texname = '\\text{Subsuperscript}\\left[m,H_d,2\\right]',
                 lhablock = 'MSOFT',
                 lhacode = [ 21 ])

mHu2 = Parameter(name = 'mHu2',
                 nature = 'external',
                 type = 'real',
                 value = -128800.134,
                 texname = '\\text{Subsuperscript}\\left[m,H_u,2\\right]',
                 lhablock = 'MSOFT',
                 lhacode = [ 22 ])

RmQ21x1 = Parameter(name = 'RmQ21x1',
                    nature = 'external',
                    type = 'real',
                    value = 309836.701,
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
                    value = 290382.106,
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
                   value = 0.707106781,
                   texname = '\\text{RRl1x1}',
                   lhablock = 'SELMIX',
                   lhacode = [ 1, 1 ])

RRl1x4 = Parameter(name = 'RRl1x4',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRl1x4}',
                   lhablock = 'SELMIX',
                   lhacode = [ 1, 4 ])

RRl2x2 = Parameter(name = 'RRl2x2',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRl2x2}',
                   lhablock = 'SELMIX',
                   lhacode = [ 2, 2 ])

RRl2x5 = Parameter(name = 'RRl2x5',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRl2x5}',
                   lhablock = 'SELMIX',
                   lhacode = [ 2, 5 ])

RRl3x3 = Parameter(name = 'RRl3x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRl3x3}',
                   lhablock = 'SELMIX',
                   lhacode = [ 3, 3 ])

RRl3x6 = Parameter(name = 'RRl3x6',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRl3x6}',
                   lhablock = 'SELMIX',
                   lhacode = [ 3, 6 ])

RRl4x1 = Parameter(name = 'RRl4x1',
                   nature = 'external',
                   type = 'real',
                   value = -0.707106781,
                   texname = '\\text{RRl4x1}',
                   lhablock = 'SELMIX',
                   lhacode = [ 4, 1 ])

RRl4x4 = Parameter(name = 'RRl4x4',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRl4x4}',
                   lhablock = 'SELMIX',
                   lhacode = [ 4, 4 ])

RRl5x2 = Parameter(name = 'RRl5x2',
                   nature = 'external',
                   type = 'real',
                   value = -0.707106781,
                   texname = '\\text{RRl5x2}',
                   lhablock = 'SELMIX',
                   lhacode = [ 5, 2 ])

RRl5x5 = Parameter(name = 'RRl5x5',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRl5x5}',
                   lhablock = 'SELMIX',
                   lhacode = [ 5, 5 ])

RRl6x3 = Parameter(name = 'RRl6x3',
                   nature = 'external',
                   type = 'real',
                   value = -0.707106781,
                   texname = '\\text{RRl6x3}',
                   lhablock = 'SELMIX',
                   lhacode = [ 6, 3 ])

RRl6x6 = Parameter(name = 'RRl6x6',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRl6x6}',
                   lhablock = 'SELMIX',
                   lhacode = [ 6, 6 ])

aEWM1 = Parameter(name = 'aEWM1',
                  nature = 'external',
                  type = 'real',
                  value = 127.934,
                  texname = '\\text{Subsuperscript}[\\alpha ,w,-1]',
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

Rtd1x1 = Parameter(name = 'Rtd1x1',
                   nature = 'external',
                   type = 'real',
                   value = -108.693742,
                   texname = '\\text{Rtd1x1}',
                   lhablock = 'TD',
                   lhacode = [ 1, 1 ])

Rtd2x2 = Parameter(name = 'Rtd2x2',
                   nature = 'external',
                   type = 'real',
                   value = -109.693742,
                   texname = '\\text{Rtd2x2}',
                   lhablock = 'TD',
                   lhacode = [ 2, 2 ])

Rtd3x3 = Parameter(name = 'Rtd3x3',
                   nature = 'external',
                   type = 'real',
                   value = -110.693742,
                   texname = '\\text{Rtd3x3}',
                   lhablock = 'TD',
                   lhacode = [ 3, 3 ])

Rte1x1 = Parameter(name = 'Rte1x1',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{Rte1x1}',
                   lhablock = 'TE',
                   lhacode = [ 1, 1 ])

Rte2x2 = Parameter(name = 'Rte2x2',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{Rte2x2}',
                   lhablock = 'TE',
                   lhacode = [ 2, 2 ])

Rte3x3 = Parameter(name = 'Rte3x3',
                   nature = 'external',
                   type = 'real',
                   value = -25.4019727,
                   texname = '\\text{Rte3x3}',
                   lhablock = 'TE',
                   lhacode = [ 3, 3 ])

Rtu1x1 = Parameter(name = 'Rtu1x1',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{Rtu1x1}',
                   lhablock = 'TU',
                   lhacode = [ 1, 1 ])

Rtu2x2 = Parameter(name = 'Rtu2x2',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{Rtu2x2}',
                   lhablock = 'TU',
                   lhacode = [ 2, 2 ])

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
                   value = 0.707106781,
                   texname = '\\text{RRu1x1}',
                   lhablock = 'USQMIX',
                   lhacode = [ 1, 1 ])

RRu1x4 = Parameter(name = 'RRu1x4',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRu1x4}',
                   lhablock = 'USQMIX',
                   lhacode = [ 1, 4 ])

RRu2x2 = Parameter(name = 'RRu2x2',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRu2x2}',
                   lhablock = 'USQMIX',
                   lhacode = [ 2, 2 ])

RRu2x5 = Parameter(name = 'RRu2x5',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRu2x5}',
                   lhablock = 'USQMIX',
                   lhacode = [ 2, 5 ])

RRu3x3 = Parameter(name = 'RRu3x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRu3x3}',
                   lhablock = 'USQMIX',
                   lhacode = [ 3, 3 ])

RRu3x6 = Parameter(name = 'RRu3x6',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRu3x6}',
                   lhablock = 'USQMIX',
                   lhacode = [ 3, 6 ])

RRu4x1 = Parameter(name = 'RRu4x1',
                   nature = 'external',
                   type = 'real',
                   value = -0.707106781,
                   texname = '\\text{RRu4x1}',
                   lhablock = 'USQMIX',
                   lhacode = [ 4, 1 ])

RRu4x4 = Parameter(name = 'RRu4x4',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRu4x4}',
                   lhablock = 'USQMIX',
                   lhacode = [ 4, 4 ])

RRu5x2 = Parameter(name = 'RRu5x2',
                   nature = 'external',
                   type = 'real',
                   value = -0.707106781,
                   texname = '\\text{RRu5x2}',
                   lhablock = 'USQMIX',
                   lhacode = [ 5, 2 ])

RRu5x5 = Parameter(name = 'RRu5x5',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
                   texname = '\\text{RRu5x5}',
                   lhablock = 'USQMIX',
                   lhacode = [ 5, 5 ])

RRu6x3 = Parameter(name = 'RRu6x3',
                   nature = 'external',
                   type = 'real',
                   value = -0.707106781,
                   texname = '\\text{RRu6x3}',
                   lhablock = 'USQMIX',
                   lhacode = [ 6, 3 ])

RRu6x6 = Parameter(name = 'RRu6x6',
                   nature = 'external',
                   type = 'real',
                   value = 0.707106781,
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

Ryd1x1 = Parameter(name = 'Ryd1x1',
                   nature = 'external',
                   type = 'real',
                   value = 0.00001,
                   texname = '\\text{Ryd1x1}',
                   lhablock = 'YD',
                   lhacode = [ 1, 1 ])

Ryd2x2 = Parameter(name = 'Ryd2x2',
                   nature = 'external',
                   type = 'real',
                   value = 0.001,
                   texname = '\\text{Ryd2x2}',
                   lhablock = 'YD',
                   lhacode = [ 2, 2 ])

Ryd3x3 = Parameter(name = 'Ryd3x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.138840206,
                   texname = '\\text{Ryd3x3}',
                   lhablock = 'YD',
                   lhacode = [ 3, 3 ])

Rye1x1 = Parameter(name = 'Rye1x1',
                   nature = 'external',
                   type = 'real',
                   value = 0.00001,
                   texname = '\\text{Rye1x1}',
                   lhablock = 'YE',
                   lhacode = [ 1, 1 ])

Rye2x2 = Parameter(name = 'Rye2x2',
                   nature = 'external',
                   type = 'real',
                   value = 0.001,
                   texname = '\\text{Rye2x2}',
                   lhablock = 'YE',
                   lhacode = [ 2, 2 ])

Rye3x3 = Parameter(name = 'Rye3x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.10089081,
                   texname = '\\text{Rye3x3}',
                   lhablock = 'YE',
                   lhacode = [ 3, 3 ])

Ryu1x1 = Parameter(name = 'Ryu1x1',
                   nature = 'external',
                   type = 'real',
                   value = 0.00001,
                   texname = '\\text{Ryu1x1}',
                   lhablock = 'YU',
                   lhacode = [ 1, 1 ])

Ryu2x2 = Parameter(name = 'Ryu2x2',
                   nature = 'external',
                   type = 'real',
                   value = 0.001,
                   texname = '\\text{Ryu2x2}',
                   lhablock = 'YU',
                   lhacode = [ 2, 2 ])

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

Me = Parameter(name = 'Me',
               nature = 'external',
               type = 'real',
               value = 0.0001,
               texname = '\\text{Me}',
               lhablock = 'MASS',
               lhacode = [ 11 ])

Mm = Parameter(name = 'Mm',
               nature = 'external',
               type = 'real',
               value = 0.1,
               texname = '\\text{Mm}',
               lhablock = 'MASS',
               lhacode = [ 13 ])

Mta = Parameter(name = 'Mta',
                nature = 'external',
                type = 'real',
                value = 1.777,
                texname = '\\text{Mta}',
                lhablock = 'MASS',
                lhacode = [ 15 ])

MU = Parameter(name = 'MU',
               nature = 'external',
               type = 'real',
               value = 0.001,
               texname = 'M',
               lhablock = 'MASS',
               lhacode = [ 2 ])

MC = Parameter(name = 'MC',
               nature = 'external',
               type = 'real',
               value = 1.,
               texname = '\\text{MC}',
               lhablock = 'MASS',
               lhacode = [ 4 ])

MT = Parameter(name = 'MT',
               nature = 'external',
               type = 'real',
               value = 175.,
               texname = '\\text{MT}',
               lhablock = 'MASS',
               lhacode = [ 6 ])

MD = Parameter(name = 'MD',
               nature = 'external',
               type = 'real',
               value = 0.001,
               texname = '\\text{MD}',
               lhablock = 'MASS',
               lhacode = [ 1 ])

MS = Parameter(name = 'MS',
               nature = 'external',
               type = 'real',
               value = 0.1,
               texname = '\\text{MS}',
               lhablock = 'MASS',
               lhacode = [ 3 ])

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
                 value = 182.258326,
                 texname = '\\text{Msn1}',
                 lhablock = 'MASS',
                 lhacode = [ 1000012 ])

Msn2 = Parameter(name = 'Msn2',
                 nature = 'external',
                 type = 'real',
                 value = 184.258326,
                 texname = '\\text{Msn2}',
                 lhablock = 'MASS',
                 lhacode = [ 1000014 ])

Msn3 = Parameter(name = 'Msn3',
                 nature = 'external',
                 type = 'real',
                 value = 186.708464,
                 texname = '\\text{Msn3}',
                 lhablock = 'MASS',
                 lhacode = [ 1000016 ])

Msl1 = Parameter(name = 'Msl1',
                 nature = 'external',
                 type = 'real',
                 value = 201.91569,
                 texname = '\\text{Msl1}',
                 lhablock = 'MASS',
                 lhacode = [ 1000011 ])

Msl2 = Parameter(name = 'Msl2',
                 nature = 'external',
                 type = 'real',
                 value = 203.91569,
                 texname = '\\text{Msl2}',
                 lhablock = 'MASS',
                 lhacode = [ 1000013 ])

Msl3 = Parameter(name = 'Msl3',
                 nature = 'external',
                 type = 'real',
                 value = 135.490864,
                 texname = '\\text{Msl3}',
                 lhablock = 'MASS',
                 lhacode = [ 1000015 ])

Msl4 = Parameter(name = 'Msl4',
                 nature = 'external',
                 type = 'real',
                 value = 141.102799,
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
                 value = 562.119014,
                 texname = '\\text{Msu1}',
                 lhablock = 'MASS',
                 lhacode = [ 1000002 ])

Msu2 = Parameter(name = 'Msu2',
                 nature = 'external',
                 type = 'real',
                 value = 564.119014,
                 texname = '\\text{Msu2}',
                 lhablock = 'MASS',
                 lhacode = [ 1000004 ])

Msu3 = Parameter(name = 'Msu3',
                 nature = 'external',
                 type = 'real',
                 value = 396.668493,
                 texname = '\\text{Msu3}',
                 lhablock = 'MASS',
                 lhacode = [ 1000006 ])

Msu4 = Parameter(name = 'Msu4',
                 nature = 'external',
                 type = 'real',
                 value = 542.259265,
                 texname = '\\text{Msu4}',
                 lhablock = 'MASS',
                 lhacode = [ 2000002 ])

Msu5 = Parameter(name = 'Msu5',
                 nature = 'external',
                 type = 'real',
                 value = 544.259265,
                 texname = '\\text{Msu5}',
                 lhablock = 'MASS',
                 lhacode = [ 2000004 ])

Msu6 = Parameter(name = 'Msu6',
                 nature = 'external',
                 type = 'real',
                 value = 586.785818,
                 texname = '\\text{Msu6}',
                 lhablock = 'MASS',
                 lhacode = [ 2000006 ])

Msd1 = Parameter(name = 'Msd1',
                 nature = 'external',
                 type = 'real',
                 value = 561.441109,
                 texname = '\\text{Msd1}',
                 lhablock = 'MASS',
                 lhacode = [ 1000001 ])

Msd2 = Parameter(name = 'Msd2',
                 nature = 'external',
                 type = 'real',
                 value = 563.441109,
                 texname = '\\text{Msd2}',
                 lhablock = 'MASS',
                 lhacode = [ 1000003 ])

Msd3 = Parameter(name = 'Msd3',
                 nature = 'external',
                 type = 'real',
                 value = 515.065179,
                 texname = '\\text{Msd3}',
                 lhablock = 'MASS',
                 lhacode = [ 1000005 ])

Msd4 = Parameter(name = 'Msd4',
                 nature = 'external',
                 type = 'real',
                 value = 541.228462,
                 texname = '\\text{Msd4}',
                 lhablock = 'MASS',
                 lhacode = [ 2000001 ])

Msd5 = Parameter(name = 'Msd5',
                 nature = 'external',
                 type = 'real',
                 value = 543.228462,
                 texname = '\\text{Msd5}',
                 lhablock = 'MASS',
                 lhacode = [ 2000003 ])

Msd6 = Parameter(name = 'Msd6',
                 nature = 'external',
                 type = 'real',
                 value = 545.726676,
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

Wneu1 = Parameter(name = 'Wneu1',
                  nature = 'external',
                  type = 'real',
                  value = 0.14,
                  texname = '\\text{Wneu1}',
                  lhablock = 'DECAY',
                  lhacode = [ 1000022 ])

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
                  value = 'RNN1x1',
                  texname = '\\text{NN1x1}')

NN1x2 = Parameter(name = 'NN1x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN1x2',
                  texname = '\\text{NN1x2}')

NN1x3 = Parameter(name = 'NN1x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN1x3',
                  texname = '\\text{NN1x3}')

NN1x4 = Parameter(name = 'NN1x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN1x4',
                  texname = '\\text{NN1x4}')

NN2x1 = Parameter(name = 'NN2x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN2x1',
                  texname = '\\text{NN2x1}')

NN2x2 = Parameter(name = 'NN2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN2x2',
                  texname = '\\text{NN2x2}')

NN2x3 = Parameter(name = 'NN2x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN2x3',
                  texname = '\\text{NN2x3}')

NN2x4 = Parameter(name = 'NN2x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN2x4',
                  texname = '\\text{NN2x4}')

NN3x1 = Parameter(name = 'NN3x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN3x1',
                  texname = '\\text{NN3x1}')

NN3x2 = Parameter(name = 'NN3x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN3x2',
                  texname = '\\text{NN3x2}')

NN3x3 = Parameter(name = 'NN3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN3x3',
                  texname = '\\text{NN3x3}')

NN3x4 = Parameter(name = 'NN3x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN3x4',
                  texname = '\\text{NN3x4}')

NN4x1 = Parameter(name = 'NN4x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN4x1',
                  texname = '\\text{NN4x1}')

NN4x2 = Parameter(name = 'NN4x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN4x2',
                  texname = '\\text{NN4x2}')

NN4x3 = Parameter(name = 'NN4x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN4x3',
                  texname = '\\text{NN4x3}')

NN4x4 = Parameter(name = 'NN4x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN4x4',
                  texname = '\\text{NN4x4}')

Rd1x1 = Parameter(name = 'Rd1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd1x1',
                  texname = '\\text{Rd1x1}')

Rd1x4 = Parameter(name = 'Rd1x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd1x4',
                  texname = '\\text{Rd1x4}')

Rd2x2 = Parameter(name = 'Rd2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd2x2',
                  texname = '\\text{Rd2x2}')

Rd2x5 = Parameter(name = 'Rd2x5',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd2x5',
                  texname = '\\text{Rd2x5}')

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

Rd4x1 = Parameter(name = 'Rd4x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd4x1',
                  texname = '\\text{Rd4x1}')

Rd4x4 = Parameter(name = 'Rd4x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd4x4',
                  texname = '\\text{Rd4x4}')

Rd5x2 = Parameter(name = 'Rd5x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd5x2',
                  texname = '\\text{Rd5x2}')

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

Rl1x4 = Parameter(name = 'Rl1x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl1x4',
                  texname = '\\text{Rl1x4}')

Rl2x2 = Parameter(name = 'Rl2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl2x2',
                  texname = '\\text{Rl2x2}')

Rl2x5 = Parameter(name = 'Rl2x5',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl2x5',
                  texname = '\\text{Rl2x5}')

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

Rl4x1 = Parameter(name = 'Rl4x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl4x1',
                  texname = '\\text{Rl4x1}')

Rl4x4 = Parameter(name = 'Rl4x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl4x4',
                  texname = '\\text{Rl4x4}')

Rl5x2 = Parameter(name = 'Rl5x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl5x2',
                  texname = '\\text{Rl5x2}')

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

Ru1x4 = Parameter(name = 'Ru1x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu1x4',
                  texname = '\\text{Ru1x4}')

Ru2x2 = Parameter(name = 'Ru2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu2x2',
                  texname = '\\text{Ru2x2}')

Ru2x5 = Parameter(name = 'Ru2x5',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu2x5',
                  texname = '\\text{Ru2x5}')

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

Ru4x1 = Parameter(name = 'Ru4x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu4x1',
                  texname = '\\text{Ru4x1}')

Ru4x4 = Parameter(name = 'Ru4x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu4x4',
                  texname = '\\text{Ru4x4}')

Ru5x2 = Parameter(name = 'Ru5x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu5x2',
                  texname = '\\text{Ru5x2}')

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

td1x1 = Parameter(name = 'td1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rtd1x1',
                  texname = '\\text{td1x1}')

td2x2 = Parameter(name = 'td2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rtd2x2',
                  texname = '\\text{td2x2}')

td3x3 = Parameter(name = 'td3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rtd3x3',
                  texname = '\\text{td3x3}')

te1x1 = Parameter(name = 'te1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rte1x1',
                  texname = '\\text{te1x1}')

te2x2 = Parameter(name = 'te2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rte2x2',
                  texname = '\\text{te2x2}')

te3x3 = Parameter(name = 'te3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rte3x3',
                  texname = '\\text{te3x3}')

tu1x1 = Parameter(name = 'tu1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rtu1x1',
                  texname = '\\text{tu1x1}')

tu2x2 = Parameter(name = 'tu2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rtu2x2',
                  texname = '\\text{tu2x2}')

tu3x3 = Parameter(name = 'tu3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rtu3x3',
                  texname = '\\text{tu3x3}')

yd1x1 = Parameter(name = 'yd1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Ryd1x1',
                  texname = '\\text{yd1x1}')

yd2x2 = Parameter(name = 'yd2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Ryd2x2',
                  texname = '\\text{yd2x2}')

yd3x3 = Parameter(name = 'yd3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Ryd3x3',
                  texname = '\\text{yd3x3}')

ye1x1 = Parameter(name = 'ye1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rye1x1',
                  texname = '\\text{ye1x1}')

ye2x2 = Parameter(name = 'ye2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rye2x2',
                  texname = '\\text{ye2x2}')

ye3x3 = Parameter(name = 'ye3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rye3x3',
                  texname = '\\text{ye3x3}')

yu1x1 = Parameter(name = 'yu1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Ryu1x1',
                  texname = '\\text{yu1x1}')

yu2x2 = Parameter(name = 'yu2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Ryu2x2',
                  texname = '\\text{yu2x2}')

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

I1a11 = Parameter(name = 'I1a11',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complexconjugate(CKM1x1)*complexconjugate(yu1x1)',
                  texname = '\\text{I1a11}')

I1a22 = Parameter(name = 'I1a22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complexconjugate(CKM2x2)*complexconjugate(yu2x2)',
                  texname = '\\text{I1a22}')

I1a33 = Parameter(name = 'I1a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complexconjugate(CKM3x3)*complexconjugate(yu3x3)',
                  texname = '\\text{I1a33}')

I10a11 = Parameter(name = 'I10a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(yd1x1)',
                   texname = '\\text{I10a11}')

I10a14 = Parameter(name = 'I10a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(yd1x1)',
                   texname = '\\text{I10a14}')

I10a22 = Parameter(name = 'I10a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(yd2x2)',
                   texname = '\\text{I10a22}')

I10a25 = Parameter(name = 'I10a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(yd2x2)',
                   texname = '\\text{I10a25}')

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

I100a14 = Parameter(name = 'I100a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*complexconjugate(Rl4x1)',
                    texname = '\\text{I100a14}')

I100a22 = Parameter(name = 'I100a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl2x2)',
                    texname = '\\text{I100a22}')

I100a25 = Parameter(name = 'I100a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl5x2)',
                    texname = '\\text{I100a25}')

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

I101a11 = Parameter(name = 'I101a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                    texname = '\\text{I101a11}')

I101a14 = Parameter(name = 'I101a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                    texname = '\\text{I101a14}')

I101a22 = Parameter(name = 'I101a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                    texname = '\\text{I101a22}')

I101a25 = Parameter(name = 'I101a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                    texname = '\\text{I101a25}')

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

I102a11 = Parameter(name = 'I102a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rl1x4*ye1x1*complexconjugate(Rn1x1)',
                    texname = '\\text{I102a11}')

I102a14 = Parameter(name = 'I102a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rl4x4*ye1x1*complexconjugate(Rn1x1)',
                    texname = '\\text{I102a14}')

I102a22 = Parameter(name = 'I102a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rl2x5*ye2x2*complexconjugate(Rn2x2)',
                    texname = '\\text{I102a22}')

I102a25 = Parameter(name = 'I102a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rl5x5*ye2x2*complexconjugate(Rn2x2)',
                    texname = '\\text{I102a25}')

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

I104a11 = Parameter(name = 'I104a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*complexconjugate(ye1x1)',
                    texname = '\\text{I104a11}')

I104a22 = Parameter(name = 'I104a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(ye2x2)',
                    texname = '\\text{I104a22}')

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

I105a14 = Parameter(name = 'I105a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*complexconjugate(Rl4x1)',
                    texname = '\\text{I105a14}')

I105a22 = Parameter(name = 'I105a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl2x2)',
                    texname = '\\text{I105a22}')

I105a25 = Parameter(name = 'I105a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl5x2)',
                    texname = '\\text{I105a25}')

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

I106a11 = Parameter(name = 'I106a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*complexconjugate(Rl1x4)*complexconjugate(te1x1)',
                    texname = '\\text{I106a11}')

I106a14 = Parameter(name = 'I106a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*complexconjugate(Rl4x4)*complexconjugate(te1x1)',
                    texname = '\\text{I106a14}')

I106a22 = Parameter(name = 'I106a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl2x5)*complexconjugate(te2x2)',
                    texname = '\\text{I106a22}')

I106a25 = Parameter(name = 'I106a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl5x5)*complexconjugate(te2x2)',
                    texname = '\\text{I106a25}')

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

I107a11 = Parameter(name = 'I107a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                    texname = '\\text{I107a11}')

I107a14 = Parameter(name = 'I107a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                    texname = '\\text{I107a14}')

I107a22 = Parameter(name = 'I107a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                    texname = '\\text{I107a22}')

I107a25 = Parameter(name = 'I107a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                    texname = '\\text{I107a25}')

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

I108a11 = Parameter(name = 'I108a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*ye1x1*complexconjugate(Rl1x1)*complexconjugate(ye1x1)',
                    texname = '\\text{I108a11}')

I108a14 = Parameter(name = 'I108a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*ye1x1*complexconjugate(Rl4x1)*complexconjugate(ye1x1)',
                    texname = '\\text{I108a14}')

I108a22 = Parameter(name = 'I108a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*ye2x2*complexconjugate(Rl2x2)*complexconjugate(ye2x2)',
                    texname = '\\text{I108a22}')

I108a25 = Parameter(name = 'I108a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*ye2x2*complexconjugate(Rl5x2)*complexconjugate(ye2x2)',
                    texname = '\\text{I108a25}')

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

I109a14 = Parameter(name = 'I109a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*complexconjugate(Rl4x1)',
                    texname = '\\text{I109a14}')

I109a22 = Parameter(name = 'I109a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl2x2)',
                    texname = '\\text{I109a22}')

I109a25 = Parameter(name = 'I109a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl5x2)',
                    texname = '\\text{I109a25}')

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

I11a11 = Parameter(name = 'I11a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*yd1x1',
                   texname = '\\text{I11a11}')

I11a14 = Parameter(name = 'I11a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*yd1x1',
                   texname = '\\text{I11a14}')

I11a22 = Parameter(name = 'I11a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*yd2x2',
                   texname = '\\text{I11a22}')

I11a25 = Parameter(name = 'I11a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*yd2x2',
                   texname = '\\text{I11a25}')

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

I110a11 = Parameter(name = 'I110a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*ye1x1*complexconjugate(Rl1x1)*complexconjugate(ye1x1)',
                    texname = '\\text{I110a11}')

I110a14 = Parameter(name = 'I110a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*ye1x1*complexconjugate(Rl4x1)*complexconjugate(ye1x1)',
                    texname = '\\text{I110a14}')

I110a22 = Parameter(name = 'I110a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*ye2x2*complexconjugate(Rl2x2)*complexconjugate(ye2x2)',
                    texname = '\\text{I110a22}')

I110a25 = Parameter(name = 'I110a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*ye2x2*complexconjugate(Rl5x2)*complexconjugate(ye2x2)',
                    texname = '\\text{I110a25}')

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

I111a11 = Parameter(name = 'I111a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I111a11}')

I111a14 = Parameter(name = 'I111a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I111a14}')

I111a22 = Parameter(name = 'I111a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I111a22}')

I111a25 = Parameter(name = 'I111a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I111a25}')

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

I112a11 = Parameter(name = 'I112a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'yu1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I112a11}')

I112a14 = Parameter(name = 'I112a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'yu1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I112a14}')

I112a22 = Parameter(name = 'I112a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'yu2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I112a22}')

I112a25 = Parameter(name = 'I112a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'yu2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I112a25}')

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

I113a14 = Parameter(name = 'I113a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I113a14}')

I113a22 = Parameter(name = 'I113a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I113a22}')

I113a25 = Parameter(name = 'I113a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I113a25}')

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

I113a41 = Parameter(name = 'I113a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I113a41}')

I113a44 = Parameter(name = 'I113a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I113a44}')

I113a52 = Parameter(name = 'I113a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I113a52}')

I113a55 = Parameter(name = 'I113a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I113a55}')

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

I114a11 = Parameter(name = 'I114a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*complexconjugate(Ru1x4)',
                    texname = '\\text{I114a11}')

I114a14 = Parameter(name = 'I114a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru1x4)',
                    texname = '\\text{I114a14}')

I114a22 = Parameter(name = 'I114a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*complexconjugate(Ru2x5)',
                    texname = '\\text{I114a22}')

I114a25 = Parameter(name = 'I114a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru2x5)',
                    texname = '\\text{I114a25}')

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

I114a41 = Parameter(name = 'I114a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I114a41}')

I114a44 = Parameter(name = 'I114a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I114a44}')

I114a52 = Parameter(name = 'I114a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I114a52}')

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

I115a14 = Parameter(name = 'I115a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)',
                    texname = '\\text{I115a14}')

I115a22 = Parameter(name = 'I115a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I115a22}')

I115a25 = Parameter(name = 'I115a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)',
                    texname = '\\text{I115a25}')

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

I115a41 = Parameter(name = 'I115a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                    texname = '\\text{I115a41}')

I115a44 = Parameter(name = 'I115a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)',
                    texname = '\\text{I115a44}')

I115a52 = Parameter(name = 'I115a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I115a52}')

I115a55 = Parameter(name = 'I115a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)',
                    texname = '\\text{I115a55}')

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

I116a11 = Parameter(name = 'I116a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I116a11}')

I116a14 = Parameter(name = 'I116a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I116a14}')

I116a22 = Parameter(name = 'I116a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I116a22}')

I116a25 = Parameter(name = 'I116a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I116a25}')

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

I116a41 = Parameter(name = 'I116a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I116a41}')

I116a44 = Parameter(name = 'I116a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I116a44}')

I116a52 = Parameter(name = 'I116a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I116a52}')

I116a55 = Parameter(name = 'I116a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I116a55}')

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

I117a11 = Parameter(name = 'I117a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x4)*complexconjugate(tu1x1)',
                    texname = '\\text{I117a11}')

I117a14 = Parameter(name = 'I117a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x4)*complexconjugate(tu1x1)',
                    texname = '\\text{I117a14}')

I117a22 = Parameter(name = 'I117a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x5)*complexconjugate(tu2x2)',
                    texname = '\\text{I117a22}')

I117a25 = Parameter(name = 'I117a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x5)*complexconjugate(tu2x2)',
                    texname = '\\text{I117a25}')

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

I117a41 = Parameter(name = 'I117a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x4)*complexconjugate(tu1x1)',
                    texname = '\\text{I117a41}')

I117a44 = Parameter(name = 'I117a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x4)*complexconjugate(tu1x1)',
                    texname = '\\text{I117a44}')

I117a52 = Parameter(name = 'I117a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x5)*complexconjugate(tu2x2)',
                    texname = '\\text{I117a52}')

I117a55 = Parameter(name = 'I117a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x5)*complexconjugate(tu2x2)',
                    texname = '\\text{I117a55}')

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

I118a11 = Parameter(name = 'I118a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x4*td1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                    texname = '\\text{I118a11}')

I118a14 = Parameter(name = 'I118a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x4*td1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)',
                    texname = '\\text{I118a14}')

I118a22 = Parameter(name = 'I118a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x5*td2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I118a22}')

I118a25 = Parameter(name = 'I118a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x5*td2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)',
                    texname = '\\text{I118a25}')

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

I118a41 = Parameter(name = 'I118a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x4*td1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                    texname = '\\text{I118a41}')

I118a44 = Parameter(name = 'I118a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x4*td1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)',
                    texname = '\\text{I118a44}')

I118a52 = Parameter(name = 'I118a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x5*td2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I118a52}')

I118a55 = Parameter(name = 'I118a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x5*td2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)',
                    texname = '\\text{I118a55}')

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

I119a11 = Parameter(name = 'I119a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)*complexconjugate(yd1x1)',
                    texname = '\\text{I119a11}')

I119a14 = Parameter(name = 'I119a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)*complexconjugate(yd1x1)',
                    texname = '\\text{I119a14}')

I119a22 = Parameter(name = 'I119a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)*complexconjugate(yd2x2)',
                    texname = '\\text{I119a22}')

I119a25 = Parameter(name = 'I119a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)*complexconjugate(yd2x2)',
                    texname = '\\text{I119a25}')

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

I119a41 = Parameter(name = 'I119a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)*complexconjugate(yd1x1)',
                    texname = '\\text{I119a41}')

I119a44 = Parameter(name = 'I119a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)*complexconjugate(yd1x1)',
                    texname = '\\text{I119a44}')

I119a52 = Parameter(name = 'I119a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)*complexconjugate(yd2x2)',
                    texname = '\\text{I119a52}')

I119a55 = Parameter(name = 'I119a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)*complexconjugate(yd2x2)',
                    texname = '\\text{I119a55}')

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

I12a14 = Parameter(name = 'I12a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(CKM1x1)',
                   texname = '\\text{I12a14}')

I12a22 = Parameter(name = 'I12a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(CKM2x2)',
                   texname = '\\text{I12a22}')

I12a25 = Parameter(name = 'I12a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(CKM2x2)',
                   texname = '\\text{I12a25}')

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

I120a11 = Parameter(name = 'I120a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x4*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                    texname = '\\text{I120a11}')

I120a14 = Parameter(name = 'I120a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x4*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)',
                    texname = '\\text{I120a14}')

I120a22 = Parameter(name = 'I120a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x5*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I120a22}')

I120a25 = Parameter(name = 'I120a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x5*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)',
                    texname = '\\text{I120a25}')

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

I120a41 = Parameter(name = 'I120a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x4*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                    texname = '\\text{I120a41}')

I120a44 = Parameter(name = 'I120a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x4*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)',
                    texname = '\\text{I120a44}')

I120a52 = Parameter(name = 'I120a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x5*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I120a52}')

I120a55 = Parameter(name = 'I120a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x5*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)',
                    texname = '\\text{I120a55}')

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

I121a11 = Parameter(name = 'I121a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x4*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I121a11}')

I121a14 = Parameter(name = 'I121a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x4*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I121a14}')

I121a22 = Parameter(name = 'I121a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x5*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I121a22}')

I121a25 = Parameter(name = 'I121a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x5*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I121a25}')

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

I121a41 = Parameter(name = 'I121a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x4*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I121a41}')

I121a44 = Parameter(name = 'I121a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x4*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I121a44}')

I121a52 = Parameter(name = 'I121a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x5*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I121a52}')

I121a55 = Parameter(name = 'I121a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x5*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I121a55}')

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

I122a11 = Parameter(name = 'I122a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*yu1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I122a11}')

I122a14 = Parameter(name = 'I122a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*yu1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I122a14}')

I122a22 = Parameter(name = 'I122a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*yu2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I122a22}')

I122a25 = Parameter(name = 'I122a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*yu2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I122a25}')

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

I122a41 = Parameter(name = 'I122a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*yu1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I122a41}')

I122a44 = Parameter(name = 'I122a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*yu1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I122a44}')

I122a52 = Parameter(name = 'I122a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*yu2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I122a52}')

I122a55 = Parameter(name = 'I122a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*yu2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I122a55}')

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

I123a14 = Parameter(name = 'I123a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)',
                    texname = '\\text{I123a14}')

I123a22 = Parameter(name = 'I123a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I123a22}')

I123a25 = Parameter(name = 'I123a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)',
                    texname = '\\text{I123a25}')

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

I123a41 = Parameter(name = 'I123a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                    texname = '\\text{I123a41}')

I123a44 = Parameter(name = 'I123a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)',
                    texname = '\\text{I123a44}')

I123a52 = Parameter(name = 'I123a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I123a52}')

I123a55 = Parameter(name = 'I123a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)',
                    texname = '\\text{I123a55}')

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

I124a11 = Parameter(name = 'I124a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)*complexconjugate(yd1x1)',
                    texname = '\\text{I124a11}')

I124a14 = Parameter(name = 'I124a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)*complexconjugate(yd1x1)',
                    texname = '\\text{I124a14}')

I124a22 = Parameter(name = 'I124a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)*complexconjugate(yd2x2)',
                    texname = '\\text{I124a22}')

I124a25 = Parameter(name = 'I124a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)*complexconjugate(yd2x2)',
                    texname = '\\text{I124a25}')

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

I124a41 = Parameter(name = 'I124a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)*complexconjugate(yd1x1)',
                    texname = '\\text{I124a41}')

I124a44 = Parameter(name = 'I124a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)*complexconjugate(yd1x1)',
                    texname = '\\text{I124a44}')

I124a52 = Parameter(name = 'I124a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)*complexconjugate(yd2x2)',
                    texname = '\\text{I124a52}')

I124a55 = Parameter(name = 'I124a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)*complexconjugate(yd2x2)',
                    texname = '\\text{I124a55}')

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

I125a11 = Parameter(name = 'I125a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x4*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I125a11}')

I125a14 = Parameter(name = 'I125a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x4*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I125a14}')

I125a22 = Parameter(name = 'I125a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x5*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I125a22}')

I125a25 = Parameter(name = 'I125a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x5*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I125a25}')

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

I125a41 = Parameter(name = 'I125a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x4*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I125a41}')

I125a44 = Parameter(name = 'I125a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x4*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I125a44}')

I125a52 = Parameter(name = 'I125a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x5*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I125a52}')

I125a55 = Parameter(name = 'I125a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x5*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I125a55}')

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

I126a11 = Parameter(name = 'I126a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*yu1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I126a11}')

I126a14 = Parameter(name = 'I126a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*yu1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I126a14}')

I126a22 = Parameter(name = 'I126a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*yu2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I126a22}')

I126a25 = Parameter(name = 'I126a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*yu2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I126a25}')

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

I126a41 = Parameter(name = 'I126a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*yu1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I126a41}')

I126a44 = Parameter(name = 'I126a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*yu1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I126a44}')

I126a52 = Parameter(name = 'I126a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*yu2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I126a52}')

I126a55 = Parameter(name = 'I126a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*yu2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I126a55}')

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

I127a14 = Parameter(name = 'I127a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)',
                    texname = '\\text{I127a14}')

I127a22 = Parameter(name = 'I127a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I127a22}')

I127a25 = Parameter(name = 'I127a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)',
                    texname = '\\text{I127a25}')

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

I127a41 = Parameter(name = 'I127a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                    texname = '\\text{I127a41}')

I127a44 = Parameter(name = 'I127a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)',
                    texname = '\\text{I127a44}')

I127a52 = Parameter(name = 'I127a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I127a52}')

I127a55 = Parameter(name = 'I127a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)',
                    texname = '\\text{I127a55}')

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

I128a14 = Parameter(name = 'I128a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*complexconjugate(Rl4x1)',
                    texname = '\\text{I128a14}')

I128a22 = Parameter(name = 'I128a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl2x2)',
                    texname = '\\text{I128a22}')

I128a25 = Parameter(name = 'I128a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl5x2)',
                    texname = '\\text{I128a25}')

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

I129a11 = Parameter(name = 'I129a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                    texname = '\\text{I129a11}')

I129a14 = Parameter(name = 'I129a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                    texname = '\\text{I129a14}')

I129a22 = Parameter(name = 'I129a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                    texname = '\\text{I129a22}')

I129a25 = Parameter(name = 'I129a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                    texname = '\\text{I129a25}')

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

I13a11 = Parameter(name = 'I13a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(CKM1x1)*complexconjugate(yu1x1)',
                   texname = '\\text{I13a11}')

I13a14 = Parameter(name = 'I13a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(CKM1x1)*complexconjugate(yu1x1)',
                   texname = '\\text{I13a14}')

I13a22 = Parameter(name = 'I13a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(yu2x2)',
                   texname = '\\text{I13a22}')

I13a25 = Parameter(name = 'I13a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(CKM2x2)*complexconjugate(yu2x2)',
                   texname = '\\text{I13a25}')

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

I130a11 = Parameter(name = 'I130a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x4*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                    texname = '\\text{I130a11}')

I130a14 = Parameter(name = 'I130a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x4*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)',
                    texname = '\\text{I130a14}')

I130a22 = Parameter(name = 'I130a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x5*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I130a22}')

I130a25 = Parameter(name = 'I130a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x5*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)',
                    texname = '\\text{I130a25}')

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

I130a41 = Parameter(name = 'I130a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x4*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                    texname = '\\text{I130a41}')

I130a44 = Parameter(name = 'I130a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x4*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)',
                    texname = '\\text{I130a44}')

I130a52 = Parameter(name = 'I130a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x5*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I130a52}')

I130a55 = Parameter(name = 'I130a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x5*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)',
                    texname = '\\text{I130a55}')

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

I131a11 = Parameter(name = 'I131a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(yu1x1)',
                    texname = '\\text{I131a11}')

I131a14 = Parameter(name = 'I131a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(yu1x1)',
                    texname = '\\text{I131a14}')

I131a22 = Parameter(name = 'I131a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(yu2x2)',
                    texname = '\\text{I131a22}')

I131a25 = Parameter(name = 'I131a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(yu2x2)',
                    texname = '\\text{I131a25}')

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

I132a11 = Parameter(name = 'I132a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*yu1x1',
                    texname = '\\text{I132a11}')

I132a14 = Parameter(name = 'I132a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*yu1x1',
                    texname = '\\text{I132a14}')

I132a22 = Parameter(name = 'I132a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*yu2x2',
                    texname = '\\text{I132a22}')

I132a25 = Parameter(name = 'I132a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*yu2x2',
                    texname = '\\text{I132a25}')

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

I133a14 = Parameter(name = 'I133a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1',
                    texname = '\\text{I133a14}')

I133a22 = Parameter(name = 'I133a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2',
                    texname = '\\text{I133a22}')

I133a25 = Parameter(name = 'I133a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2',
                    texname = '\\text{I133a25}')

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

I134a11 = Parameter(name = 'I134a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*complexconjugate(yd1x1)',
                    texname = '\\text{I134a11}')

I134a14 = Parameter(name = 'I134a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*complexconjugate(yd1x1)',
                    texname = '\\text{I134a14}')

I134a22 = Parameter(name = 'I134a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(yd2x2)',
                    texname = '\\text{I134a22}')

I134a25 = Parameter(name = 'I134a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*complexconjugate(yd2x2)',
                    texname = '\\text{I134a25}')

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

I135a11 = Parameter(name = 'I135a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x4*yu1x1',
                    texname = '\\text{I135a11}')

I135a14 = Parameter(name = 'I135a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x4*yu1x1',
                    texname = '\\text{I135a14}')

I135a22 = Parameter(name = 'I135a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x5*yu2x2',
                    texname = '\\text{I135a22}')

I135a25 = Parameter(name = 'I135a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x5*yu2x2',
                    texname = '\\text{I135a25}')

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

I136a14 = Parameter(name = 'I136a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*complexconjugate(Rd1x1)',
                    texname = '\\text{I136a14}')

I136a22 = Parameter(name = 'I136a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I136a22}')

I136a25 = Parameter(name = 'I136a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I136a25}')

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

I136a41 = Parameter(name = 'I136a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*complexconjugate(Rd4x1)',
                    texname = '\\text{I136a41}')

I136a44 = Parameter(name = 'I136a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*complexconjugate(Rd4x1)',
                    texname = '\\text{I136a44}')

I136a52 = Parameter(name = 'I136a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(Rd5x2)',
                    texname = '\\text{I136a52}')

I136a55 = Parameter(name = 'I136a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*complexconjugate(Rd5x2)',
                    texname = '\\text{I136a55}')

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

I137a11 = Parameter(name = 'I137a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*complexconjugate(Rd1x4)*complexconjugate(td1x1)',
                    texname = '\\text{I137a11}')

I137a14 = Parameter(name = 'I137a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*complexconjugate(Rd1x4)*complexconjugate(td1x1)',
                    texname = '\\text{I137a14}')

I137a22 = Parameter(name = 'I137a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(Rd2x5)*complexconjugate(td2x2)',
                    texname = '\\text{I137a22}')

I137a25 = Parameter(name = 'I137a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*complexconjugate(Rd2x5)*complexconjugate(td2x2)',
                    texname = '\\text{I137a25}')

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

I137a41 = Parameter(name = 'I137a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*complexconjugate(Rd4x4)*complexconjugate(td1x1)',
                    texname = '\\text{I137a41}')

I137a44 = Parameter(name = 'I137a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*complexconjugate(Rd4x4)*complexconjugate(td1x1)',
                    texname = '\\text{I137a44}')

I137a52 = Parameter(name = 'I137a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(Rd5x5)*complexconjugate(td2x2)',
                    texname = '\\text{I137a52}')

I137a55 = Parameter(name = 'I137a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*complexconjugate(Rd5x5)*complexconjugate(td2x2)',
                    texname = '\\text{I137a55}')

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

I138a11 = Parameter(name = 'I138a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I138a11}')

I138a14 = Parameter(name = 'I138a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I138a14}')

I138a22 = Parameter(name = 'I138a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I138a22}')

I138a25 = Parameter(name = 'I138a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I138a25}')

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

I138a41 = Parameter(name = 'I138a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I138a41}')

I138a44 = Parameter(name = 'I138a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I138a44}')

I138a52 = Parameter(name = 'I138a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I138a52}')

I138a55 = Parameter(name = 'I138a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I138a55}')

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

I139a11 = Parameter(name = 'I139a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x4*tu1x1*complexconjugate(Rd1x1)',
                    texname = '\\text{I139a11}')

I139a14 = Parameter(name = 'I139a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x4*tu1x1*complexconjugate(Rd1x1)',
                    texname = '\\text{I139a14}')

I139a22 = Parameter(name = 'I139a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x5*tu2x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I139a22}')

I139a25 = Parameter(name = 'I139a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x5*tu2x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I139a25}')

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

I139a41 = Parameter(name = 'I139a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x4*tu1x1*complexconjugate(Rd4x1)',
                    texname = '\\text{I139a41}')

I139a44 = Parameter(name = 'I139a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x4*tu1x1*complexconjugate(Rd4x1)',
                    texname = '\\text{I139a44}')

I139a52 = Parameter(name = 'I139a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x5*tu2x2*complexconjugate(Rd5x2)',
                    texname = '\\text{I139a52}')

I139a55 = Parameter(name = 'I139a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x5*tu2x2*complexconjugate(Rd5x2)',
                    texname = '\\text{I139a55}')

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

I14a11 = Parameter(name = 'I14a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*yd1x1*complexconjugate(CKM1x1)',
                   texname = '\\text{I14a11}')

I14a14 = Parameter(name = 'I14a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*yd1x1*complexconjugate(CKM1x1)',
                   texname = '\\text{I14a14}')

I14a22 = Parameter(name = 'I14a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*yd2x2*complexconjugate(CKM2x2)',
                   texname = '\\text{I14a22}')

I14a25 = Parameter(name = 'I14a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*yd2x2*complexconjugate(CKM2x2)',
                   texname = '\\text{I14a25}')

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

I140a11 = Parameter(name = 'I140a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*yd1x1*complexconjugate(Rd1x1)*complexconjugate(yd1x1)',
                    texname = '\\text{I140a11}')

I140a14 = Parameter(name = 'I140a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*yd1x1*complexconjugate(Rd1x1)*complexconjugate(yd1x1)',
                    texname = '\\text{I140a14}')

I140a22 = Parameter(name = 'I140a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*yd2x2*complexconjugate(Rd2x2)*complexconjugate(yd2x2)',
                    texname = '\\text{I140a22}')

I140a25 = Parameter(name = 'I140a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*yd2x2*complexconjugate(Rd2x2)*complexconjugate(yd2x2)',
                    texname = '\\text{I140a25}')

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

I140a41 = Parameter(name = 'I140a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*yd1x1*complexconjugate(Rd4x1)*complexconjugate(yd1x1)',
                    texname = '\\text{I140a41}')

I140a44 = Parameter(name = 'I140a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*yd1x1*complexconjugate(Rd4x1)*complexconjugate(yd1x1)',
                    texname = '\\text{I140a44}')

I140a52 = Parameter(name = 'I140a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*yd2x2*complexconjugate(Rd5x2)*complexconjugate(yd2x2)',
                    texname = '\\text{I140a52}')

I140a55 = Parameter(name = 'I140a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*yd2x2*complexconjugate(Rd5x2)*complexconjugate(yd2x2)',
                    texname = '\\text{I140a55}')

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

I141a11 = Parameter(name = 'I141a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x4*yu1x1*complexconjugate(Rd1x1)',
                    texname = '\\text{I141a11}')

I141a14 = Parameter(name = 'I141a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x4*yu1x1*complexconjugate(Rd1x1)',
                    texname = '\\text{I141a14}')

I141a22 = Parameter(name = 'I141a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x5*yu2x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I141a22}')

I141a25 = Parameter(name = 'I141a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x5*yu2x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I141a25}')

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

I141a41 = Parameter(name = 'I141a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x4*yu1x1*complexconjugate(Rd4x1)',
                    texname = '\\text{I141a41}')

I141a44 = Parameter(name = 'I141a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x4*yu1x1*complexconjugate(Rd4x1)',
                    texname = '\\text{I141a44}')

I141a52 = Parameter(name = 'I141a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x5*yu2x2*complexconjugate(Rd5x2)',
                    texname = '\\text{I141a52}')

I141a55 = Parameter(name = 'I141a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x5*yu2x2*complexconjugate(Rd5x2)',
                    texname = '\\text{I141a55}')

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

I142a11 = Parameter(name = 'I142a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*yu1x1*complexconjugate(Rd1x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I142a11}')

I142a14 = Parameter(name = 'I142a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*yu1x1*complexconjugate(Rd1x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I142a14}')

I142a22 = Parameter(name = 'I142a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*yu2x2*complexconjugate(Rd2x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I142a22}')

I142a25 = Parameter(name = 'I142a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*yu2x2*complexconjugate(Rd2x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I142a25}')

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

I142a41 = Parameter(name = 'I142a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*yu1x1*complexconjugate(Rd4x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I142a41}')

I142a44 = Parameter(name = 'I142a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*yu1x1*complexconjugate(Rd4x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I142a44}')

I142a52 = Parameter(name = 'I142a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*yu2x2*complexconjugate(Rd5x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I142a52}')

I142a55 = Parameter(name = 'I142a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*yu2x2*complexconjugate(Rd5x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I142a55}')

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

I143a11 = Parameter(name = 'I143a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x4*yu1x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I143a11}')

I143a14 = Parameter(name = 'I143a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x4*yu1x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I143a14}')

I143a22 = Parameter(name = 'I143a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x5*yu2x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I143a22}')

I143a25 = Parameter(name = 'I143a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x5*yu2x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I143a25}')

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

I143a41 = Parameter(name = 'I143a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x4*yu1x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I143a41}')

I143a44 = Parameter(name = 'I143a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x4*yu1x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I143a44}')

I143a52 = Parameter(name = 'I143a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x5*yu2x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I143a52}')

I143a55 = Parameter(name = 'I143a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x5*yu2x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I143a55}')

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

I144a14 = Parameter(name = 'I144a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*complexconjugate(Rd1x1)',
                    texname = '\\text{I144a14}')

I144a22 = Parameter(name = 'I144a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I144a22}')

I144a25 = Parameter(name = 'I144a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I144a25}')

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

I144a41 = Parameter(name = 'I144a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*complexconjugate(Rd4x1)',
                    texname = '\\text{I144a41}')

I144a44 = Parameter(name = 'I144a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*complexconjugate(Rd4x1)',
                    texname = '\\text{I144a44}')

I144a52 = Parameter(name = 'I144a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(Rd5x2)',
                    texname = '\\text{I144a52}')

I144a55 = Parameter(name = 'I144a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*complexconjugate(Rd5x2)',
                    texname = '\\text{I144a55}')

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

I145a11 = Parameter(name = 'I145a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*yd1x1*complexconjugate(Rd1x1)*complexconjugate(yd1x1)',
                    texname = '\\text{I145a11}')

I145a14 = Parameter(name = 'I145a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*yd1x1*complexconjugate(Rd1x1)*complexconjugate(yd1x1)',
                    texname = '\\text{I145a14}')

I145a22 = Parameter(name = 'I145a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*yd2x2*complexconjugate(Rd2x2)*complexconjugate(yd2x2)',
                    texname = '\\text{I145a22}')

I145a25 = Parameter(name = 'I145a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*yd2x2*complexconjugate(Rd2x2)*complexconjugate(yd2x2)',
                    texname = '\\text{I145a25}')

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

I145a41 = Parameter(name = 'I145a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*yd1x1*complexconjugate(Rd4x1)*complexconjugate(yd1x1)',
                    texname = '\\text{I145a41}')

I145a44 = Parameter(name = 'I145a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*yd1x1*complexconjugate(Rd4x1)*complexconjugate(yd1x1)',
                    texname = '\\text{I145a44}')

I145a52 = Parameter(name = 'I145a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*yd2x2*complexconjugate(Rd5x2)*complexconjugate(yd2x2)',
                    texname = '\\text{I145a52}')

I145a55 = Parameter(name = 'I145a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*yd2x2*complexconjugate(Rd5x2)*complexconjugate(yd2x2)',
                    texname = '\\text{I145a55}')

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

I146a11 = Parameter(name = 'I146a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*yu1x1*complexconjugate(Rd1x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I146a11}')

I146a14 = Parameter(name = 'I146a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*yu1x1*complexconjugate(Rd1x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I146a14}')

I146a22 = Parameter(name = 'I146a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*yu2x2*complexconjugate(Rd2x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I146a22}')

I146a25 = Parameter(name = 'I146a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*yu2x2*complexconjugate(Rd2x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I146a25}')

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

I146a41 = Parameter(name = 'I146a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*yu1x1*complexconjugate(Rd4x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I146a41}')

I146a44 = Parameter(name = 'I146a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*yu1x1*complexconjugate(Rd4x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I146a44}')

I146a52 = Parameter(name = 'I146a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*yu2x2*complexconjugate(Rd5x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I146a52}')

I146a55 = Parameter(name = 'I146a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*yu2x2*complexconjugate(Rd5x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I146a55}')

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

I147a11 = Parameter(name = 'I147a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x4*yu1x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I147a11}')

I147a14 = Parameter(name = 'I147a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x4*yu1x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I147a14}')

I147a22 = Parameter(name = 'I147a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x5*yu2x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I147a22}')

I147a25 = Parameter(name = 'I147a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x5*yu2x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I147a25}')

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

I147a41 = Parameter(name = 'I147a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x4*yu1x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I147a41}')

I147a44 = Parameter(name = 'I147a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x4*yu1x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I147a44}')

I147a52 = Parameter(name = 'I147a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x5*yu2x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I147a52}')

I147a55 = Parameter(name = 'I147a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x5*yu2x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I147a55}')

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

I148a14 = Parameter(name = 'I148a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I148a14}')

I148a22 = Parameter(name = 'I148a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I148a22}')

I148a25 = Parameter(name = 'I148a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I148a25}')

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

I148a41 = Parameter(name = 'I148a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I148a41}')

I148a44 = Parameter(name = 'I148a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I148a44}')

I148a52 = Parameter(name = 'I148a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I148a52}')

I148a55 = Parameter(name = 'I148a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I148a55}')

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

I149a11 = Parameter(name = 'I149a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*complexconjugate(Ru1x4)',
                    texname = '\\text{I149a11}')

I149a14 = Parameter(name = 'I149a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru1x4)',
                    texname = '\\text{I149a14}')

I149a22 = Parameter(name = 'I149a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*complexconjugate(Ru2x5)',
                    texname = '\\text{I149a22}')

I149a25 = Parameter(name = 'I149a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru2x5)',
                    texname = '\\text{I149a25}')

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

I149a41 = Parameter(name = 'I149a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I149a41}')

I149a44 = Parameter(name = 'I149a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I149a44}')

I149a52 = Parameter(name = 'I149a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I149a52}')

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

I15a14 = Parameter(name = 'I15a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I15a14}')

I15a22 = Parameter(name = 'I15a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I15a22}')

I15a25 = Parameter(name = 'I15a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I15a25}')

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

I15a41 = Parameter(name = 'I15a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I15a41}')

I15a44 = Parameter(name = 'I15a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I15a44}')

I15a52 = Parameter(name = 'I15a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I15a52}')

I15a55 = Parameter(name = 'I15a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I15a55}')

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

I150a11 = Parameter(name = 'I150a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I150a11}')

I150a14 = Parameter(name = 'I150a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I150a14}')

I150a22 = Parameter(name = 'I150a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I150a22}')

I150a25 = Parameter(name = 'I150a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I150a25}')

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

I150a41 = Parameter(name = 'I150a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I150a41}')

I150a44 = Parameter(name = 'I150a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I150a44}')

I150a52 = Parameter(name = 'I150a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I150a52}')

I150a55 = Parameter(name = 'I150a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I150a55}')

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

I151a11 = Parameter(name = 'I151a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru1x4)*complexconjugate(tu1x1)',
                    texname = '\\text{I151a11}')

I151a14 = Parameter(name = 'I151a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru1x4)*complexconjugate(tu1x1)',
                    texname = '\\text{I151a14}')

I151a22 = Parameter(name = 'I151a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x5)*complexconjugate(tu2x2)',
                    texname = '\\text{I151a22}')

I151a25 = Parameter(name = 'I151a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru2x5)*complexconjugate(tu2x2)',
                    texname = '\\text{I151a25}')

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

I151a41 = Parameter(name = 'I151a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru4x4)*complexconjugate(tu1x1)',
                    texname = '\\text{I151a41}')

I151a44 = Parameter(name = 'I151a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru4x4)*complexconjugate(tu1x1)',
                    texname = '\\text{I151a44}')

I151a52 = Parameter(name = 'I151a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru5x5)*complexconjugate(tu2x2)',
                    texname = '\\text{I151a52}')

I151a55 = Parameter(name = 'I151a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru5x5)*complexconjugate(tu2x2)',
                    texname = '\\text{I151a55}')

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

I152a11 = Parameter(name = 'I152a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*tu1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I152a11}')

I152a14 = Parameter(name = 'I152a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*tu1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I152a14}')

I152a22 = Parameter(name = 'I152a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*tu2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I152a22}')

I152a25 = Parameter(name = 'I152a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*tu2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I152a25}')

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

I152a41 = Parameter(name = 'I152a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*tu1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I152a41}')

I152a44 = Parameter(name = 'I152a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*tu1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I152a44}')

I152a52 = Parameter(name = 'I152a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*tu2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I152a52}')

I152a55 = Parameter(name = 'I152a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*tu2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I152a55}')

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

I153a11 = Parameter(name = 'I153a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*yu1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I153a11}')

I153a14 = Parameter(name = 'I153a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*yu1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I153a14}')

I153a22 = Parameter(name = 'I153a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*yu2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I153a22}')

I153a25 = Parameter(name = 'I153a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*yu2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I153a25}')

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

I153a41 = Parameter(name = 'I153a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*yu1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I153a41}')

I153a44 = Parameter(name = 'I153a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*yu1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I153a44}')

I153a52 = Parameter(name = 'I153a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*yu2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I153a52}')

I153a55 = Parameter(name = 'I153a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*yu2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I153a55}')

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

I154a11 = Parameter(name = 'I154a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*yu1x1*complexconjugate(Ru1x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I154a11}')

I154a14 = Parameter(name = 'I154a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*yu1x1*complexconjugate(Ru1x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I154a14}')

I154a22 = Parameter(name = 'I154a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*yu2x2*complexconjugate(Ru2x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I154a22}')

I154a25 = Parameter(name = 'I154a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*yu2x2*complexconjugate(Ru2x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I154a25}')

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

I154a41 = Parameter(name = 'I154a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*yu1x1*complexconjugate(Ru4x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I154a41}')

I154a44 = Parameter(name = 'I154a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*yu1x1*complexconjugate(Ru4x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I154a44}')

I154a52 = Parameter(name = 'I154a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*yu2x2*complexconjugate(Ru5x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I154a52}')

I154a55 = Parameter(name = 'I154a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*yu2x2*complexconjugate(Ru5x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I154a55}')

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

I155a11 = Parameter(name = 'I155a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*yu1x1*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I155a11}')

I155a14 = Parameter(name = 'I155a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*yu1x1*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I155a14}')

I155a22 = Parameter(name = 'I155a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*yu2x2*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I155a22}')

I155a25 = Parameter(name = 'I155a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*yu2x2*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I155a25}')

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

I155a41 = Parameter(name = 'I155a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*yu1x1*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I155a41}')

I155a44 = Parameter(name = 'I155a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*yu1x1*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I155a44}')

I155a52 = Parameter(name = 'I155a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*yu2x2*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I155a52}')

I155a55 = Parameter(name = 'I155a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*yu2x2*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I155a55}')

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

I156a11 = Parameter(name = 'I156a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*yu1x1*complexconjugate(Ru1x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I156a11}')

I156a14 = Parameter(name = 'I156a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*yu1x1*complexconjugate(Ru1x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I156a14}')

I156a22 = Parameter(name = 'I156a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*yu2x2*complexconjugate(Ru2x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I156a22}')

I156a25 = Parameter(name = 'I156a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*yu2x2*complexconjugate(Ru2x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I156a25}')

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

I156a41 = Parameter(name = 'I156a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*yu1x1*complexconjugate(Ru4x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I156a41}')

I156a44 = Parameter(name = 'I156a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*yu1x1*complexconjugate(Ru4x1)*complexconjugate(yu1x1)',
                    texname = '\\text{I156a44}')

I156a52 = Parameter(name = 'I156a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*yu2x2*complexconjugate(Ru5x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I156a52}')

I156a55 = Parameter(name = 'I156a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*yu2x2*complexconjugate(Ru5x2)*complexconjugate(yu2x2)',
                    texname = '\\text{I156a55}')

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

I157a11 = Parameter(name = 'I157a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*yu1x1*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I157a11}')

I157a14 = Parameter(name = 'I157a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*yu1x1*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I157a14}')

I157a22 = Parameter(name = 'I157a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*yu2x2*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I157a22}')

I157a25 = Parameter(name = 'I157a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*yu2x2*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I157a25}')

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

I157a41 = Parameter(name = 'I157a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*yu1x1*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I157a41}')

I157a44 = Parameter(name = 'I157a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*yu1x1*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I157a44}')

I157a52 = Parameter(name = 'I157a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*yu2x2*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I157a52}')

I157a55 = Parameter(name = 'I157a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*yu2x2*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I157a55}')

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

I158a11 = Parameter(name = 'I158a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)*complexconjugate(yd1x1)',
                    texname = '\\text{I158a11}')

I158a14 = Parameter(name = 'I158a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)*complexconjugate(yd1x1)',
                    texname = '\\text{I158a14}')

I158a22 = Parameter(name = 'I158a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)*complexconjugate(yd2x2)',
                    texname = '\\text{I158a22}')

I158a25 = Parameter(name = 'I158a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)*complexconjugate(yd2x2)',
                    texname = '\\text{I158a25}')

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

I158a41 = Parameter(name = 'I158a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)*complexconjugate(yd1x1)',
                    texname = '\\text{I158a41}')

I158a44 = Parameter(name = 'I158a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)*complexconjugate(yd1x1)',
                    texname = '\\text{I158a44}')

I158a52 = Parameter(name = 'I158a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)*complexconjugate(yd2x2)',
                    texname = '\\text{I158a52}')

I158a55 = Parameter(name = 'I158a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)*complexconjugate(yd2x2)',
                    texname = '\\text{I158a55}')

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

I159a14 = Parameter(name = 'I159a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*complexconjugate(Rd1x1)',
                    texname = '\\text{I159a14}')

I159a22 = Parameter(name = 'I159a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I159a22}')

I159a25 = Parameter(name = 'I159a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I159a25}')

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

I159a41 = Parameter(name = 'I159a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*complexconjugate(Rd4x1)',
                    texname = '\\text{I159a41}')

I159a44 = Parameter(name = 'I159a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*complexconjugate(Rd4x1)',
                    texname = '\\text{I159a44}')

I159a52 = Parameter(name = 'I159a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(Rd5x2)',
                    texname = '\\text{I159a52}')

I159a55 = Parameter(name = 'I159a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*complexconjugate(Rd5x2)',
                    texname = '\\text{I159a55}')

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

I16a11 = Parameter(name = 'I16a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*complexconjugate(Rd1x4)',
                   texname = '\\text{I16a11}')

I16a14 = Parameter(name = 'I16a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd1x4)',
                   texname = '\\text{I16a14}')

I16a22 = Parameter(name = 'I16a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*complexconjugate(Rd2x5)',
                   texname = '\\text{I16a22}')

I16a25 = Parameter(name = 'I16a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*complexconjugate(Rd2x5)',
                   texname = '\\text{I16a25}')

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

I16a41 = Parameter(name = 'I16a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I16a41}')

I16a44 = Parameter(name = 'I16a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I16a44}')

I16a52 = Parameter(name = 'I16a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*complexconjugate(Rd5x5)',
                   texname = '\\text{I16a52}')

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

I160a11 = Parameter(name = 'I160a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I160a11}')

I160a14 = Parameter(name = 'I160a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I160a14}')

I160a22 = Parameter(name = 'I160a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I160a22}')

I160a25 = Parameter(name = 'I160a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I160a25}')

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

I160a41 = Parameter(name = 'I160a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I160a41}')

I160a44 = Parameter(name = 'I160a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I160a44}')

I160a52 = Parameter(name = 'I160a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I160a52}')

I160a55 = Parameter(name = 'I160a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I160a55}')

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

I161a14 = Parameter(name = 'I161a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I161a14}')

I161a22 = Parameter(name = 'I161a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I161a22}')

I161a25 = Parameter(name = 'I161a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I161a25}')

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

I161a41 = Parameter(name = 'I161a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I161a41}')

I161a44 = Parameter(name = 'I161a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I161a44}')

I161a52 = Parameter(name = 'I161a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I161a52}')

I161a55 = Parameter(name = 'I161a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I161a55}')

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

I162a11 = Parameter(name = 'I162a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*complexconjugate(Ru1x4)',
                    texname = '\\text{I162a11}')

I162a14 = Parameter(name = 'I162a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru1x4)',
                    texname = '\\text{I162a14}')

I162a22 = Parameter(name = 'I162a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*complexconjugate(Ru2x5)',
                    texname = '\\text{I162a22}')

I162a25 = Parameter(name = 'I162a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru2x5)',
                    texname = '\\text{I162a25}')

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

I162a41 = Parameter(name = 'I162a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I162a41}')

I162a44 = Parameter(name = 'I162a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I162a44}')

I162a52 = Parameter(name = 'I162a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I162a52}')

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

I163a14 = Parameter(name = 'I163a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)',
                    texname = '\\text{I163a14}')

I163a22 = Parameter(name = 'I163a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I163a22}')

I163a25 = Parameter(name = 'I163a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)',
                    texname = '\\text{I163a25}')

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

I163a41 = Parameter(name = 'I163a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                    texname = '\\text{I163a41}')

I163a44 = Parameter(name = 'I163a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)',
                    texname = '\\text{I163a44}')

I163a52 = Parameter(name = 'I163a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I163a52}')

I163a55 = Parameter(name = 'I163a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)',
                    texname = '\\text{I163a55}')

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

I164a11 = Parameter(name = 'I164a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I164a11}')

I164a14 = Parameter(name = 'I164a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I164a14}')

I164a22 = Parameter(name = 'I164a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I164a22}')

I164a25 = Parameter(name = 'I164a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I164a25}')

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

I164a41 = Parameter(name = 'I164a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I164a41}')

I164a44 = Parameter(name = 'I164a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I164a44}')

I164a52 = Parameter(name = 'I164a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I164a52}')

I164a55 = Parameter(name = 'I164a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I164a55}')

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

I165a11 = Parameter(name = 'I165a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x4*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                    texname = '\\text{I165a11}')

I165a14 = Parameter(name = 'I165a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x4*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)',
                    texname = '\\text{I165a14}')

I165a22 = Parameter(name = 'I165a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x5*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I165a22}')

I165a25 = Parameter(name = 'I165a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x5*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)',
                    texname = '\\text{I165a25}')

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

I165a41 = Parameter(name = 'I165a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x4*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                    texname = '\\text{I165a41}')

I165a44 = Parameter(name = 'I165a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x4*yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)',
                    texname = '\\text{I165a44}')

I165a52 = Parameter(name = 'I165a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x5*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I165a52}')

I165a55 = Parameter(name = 'I165a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x5*yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)',
                    texname = '\\text{I165a55}')

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

I166a11 = Parameter(name = 'I166a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I166a11}')

I166a14 = Parameter(name = 'I166a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I166a14}')

I166a22 = Parameter(name = 'I166a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I166a22}')

I166a25 = Parameter(name = 'I166a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I166a25}')

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

I166a41 = Parameter(name = 'I166a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I166a41}')

I166a44 = Parameter(name = 'I166a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I166a44}')

I166a52 = Parameter(name = 'I166a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I166a52}')

I166a55 = Parameter(name = 'I166a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I166a55}')

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

I167a11 = Parameter(name = 'I167a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x4*yu1x1*complexconjugate(Rd1x1)',
                    texname = '\\text{I167a11}')

I167a14 = Parameter(name = 'I167a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x4*yu1x1*complexconjugate(Rd1x1)',
                    texname = '\\text{I167a14}')

I167a22 = Parameter(name = 'I167a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x5*yu2x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I167a22}')

I167a25 = Parameter(name = 'I167a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x5*yu2x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I167a25}')

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

I167a41 = Parameter(name = 'I167a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x4*yu1x1*complexconjugate(Rd4x1)',
                    texname = '\\text{I167a41}')

I167a44 = Parameter(name = 'I167a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x4*yu1x1*complexconjugate(Rd4x1)',
                    texname = '\\text{I167a44}')

I167a52 = Parameter(name = 'I167a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x5*yu2x2*complexconjugate(Rd5x2)',
                    texname = '\\text{I167a52}')

I167a55 = Parameter(name = 'I167a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x5*yu2x2*complexconjugate(Rd5x2)',
                    texname = '\\text{I167a55}')

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

I168a14 = Parameter(name = 'I168a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I168a14}')

I168a22 = Parameter(name = 'I168a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I168a22}')

I168a25 = Parameter(name = 'I168a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I168a25}')

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

I168a41 = Parameter(name = 'I168a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I168a41}')

I168a44 = Parameter(name = 'I168a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I168a44}')

I168a52 = Parameter(name = 'I168a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I168a52}')

I168a55 = Parameter(name = 'I168a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I168a55}')

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

I169a14 = Parameter(name = 'I169a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I169a14}')

I169a22 = Parameter(name = 'I169a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I169a22}')

I169a25 = Parameter(name = 'I169a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I169a25}')

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

I169a41 = Parameter(name = 'I169a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I169a41}')

I169a44 = Parameter(name = 'I169a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I169a44}')

I169a52 = Parameter(name = 'I169a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I169a52}')

I169a55 = Parameter(name = 'I169a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I169a55}')

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

I17a11 = Parameter(name = 'I17a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd1x4)*complexconjugate(td1x1)',
                   texname = '\\text{I17a11}')

I17a14 = Parameter(name = 'I17a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd1x4)*complexconjugate(td1x1)',
                   texname = '\\text{I17a14}')

I17a22 = Parameter(name = 'I17a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x5)*complexconjugate(td2x2)',
                   texname = '\\text{I17a22}')

I17a25 = Parameter(name = 'I17a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd2x5)*complexconjugate(td2x2)',
                   texname = '\\text{I17a25}')

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

I17a41 = Parameter(name = 'I17a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd4x4)*complexconjugate(td1x1)',
                   texname = '\\text{I17a41}')

I17a44 = Parameter(name = 'I17a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd4x4)*complexconjugate(td1x1)',
                   texname = '\\text{I17a44}')

I17a52 = Parameter(name = 'I17a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd5x5)*complexconjugate(td2x2)',
                   texname = '\\text{I17a52}')

I17a55 = Parameter(name = 'I17a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd5x5)*complexconjugate(td2x2)',
                   texname = '\\text{I17a55}')

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

I170a14 = Parameter(name = 'I170a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I170a14}')

I170a22 = Parameter(name = 'I170a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I170a22}')

I170a25 = Parameter(name = 'I170a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I170a25}')

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

I170a41 = Parameter(name = 'I170a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I170a41}')

I170a44 = Parameter(name = 'I170a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I170a44}')

I170a52 = Parameter(name = 'I170a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I170a52}')

I170a55 = Parameter(name = 'I170a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I170a55}')

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

I171a14 = Parameter(name = 'I171a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I171a14}')

I171a22 = Parameter(name = 'I171a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I171a22}')

I171a25 = Parameter(name = 'I171a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I171a25}')

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

I171a41 = Parameter(name = 'I171a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I171a41}')

I171a44 = Parameter(name = 'I171a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I171a44}')

I171a52 = Parameter(name = 'I171a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I171a52}')

I171a55 = Parameter(name = 'I171a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I171a55}')

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

I172a11 = Parameter(name = 'I172a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*complexconjugate(Ru1x4)',
                    texname = '\\text{I172a11}')

I172a14 = Parameter(name = 'I172a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru1x4)',
                    texname = '\\text{I172a14}')

I172a22 = Parameter(name = 'I172a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*complexconjugate(Ru2x5)',
                    texname = '\\text{I172a22}')

I172a25 = Parameter(name = 'I172a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru2x5)',
                    texname = '\\text{I172a25}')

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

I172a41 = Parameter(name = 'I172a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I172a41}')

I172a44 = Parameter(name = 'I172a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I172a44}')

I172a52 = Parameter(name = 'I172a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I172a52}')

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

I173a14 = Parameter(name = 'I173a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I173a14}')

I173a22 = Parameter(name = 'I173a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I173a22}')

I173a25 = Parameter(name = 'I173a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I173a25}')

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

I173a41 = Parameter(name = 'I173a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I173a41}')

I173a44 = Parameter(name = 'I173a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I173a44}')

I173a52 = Parameter(name = 'I173a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I173a52}')

I173a55 = Parameter(name = 'I173a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I173a55}')

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

I174a11 = Parameter(name = 'I174a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*complexconjugate(Ru1x4)',
                    texname = '\\text{I174a11}')

I174a14 = Parameter(name = 'I174a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru1x4)',
                    texname = '\\text{I174a14}')

I174a22 = Parameter(name = 'I174a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*complexconjugate(Ru2x5)',
                    texname = '\\text{I174a22}')

I174a25 = Parameter(name = 'I174a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru2x5)',
                    texname = '\\text{I174a25}')

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

I174a41 = Parameter(name = 'I174a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I174a41}')

I174a44 = Parameter(name = 'I174a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I174a44}')

I174a52 = Parameter(name = 'I174a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I174a52}')

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

I175a11 = Parameter(name = 'I175a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*complexconjugate(Ru1x4)',
                    texname = '\\text{I175a11}')

I175a14 = Parameter(name = 'I175a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru1x4)',
                    texname = '\\text{I175a14}')

I175a22 = Parameter(name = 'I175a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*complexconjugate(Ru2x5)',
                    texname = '\\text{I175a22}')

I175a25 = Parameter(name = 'I175a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru2x5)',
                    texname = '\\text{I175a25}')

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

I175a41 = Parameter(name = 'I175a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I175a41}')

I175a44 = Parameter(name = 'I175a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I175a44}')

I175a52 = Parameter(name = 'I175a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I175a52}')

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

I176a11 = Parameter(name = 'I176a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*complexconjugate(Ru1x4)',
                    texname = '\\text{I176a11}')

I176a14 = Parameter(name = 'I176a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru1x4)',
                    texname = '\\text{I176a14}')

I176a22 = Parameter(name = 'I176a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*complexconjugate(Ru2x5)',
                    texname = '\\text{I176a22}')

I176a25 = Parameter(name = 'I176a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru2x5)',
                    texname = '\\text{I176a25}')

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

I176a41 = Parameter(name = 'I176a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I176a41}')

I176a44 = Parameter(name = 'I176a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I176a44}')

I176a52 = Parameter(name = 'I176a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I176a52}')

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

I177a11 = Parameter(name = 'I177a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*complexconjugate(Ru1x4)',
                    texname = '\\text{I177a11}')

I177a14 = Parameter(name = 'I177a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru1x4)',
                    texname = '\\text{I177a14}')

I177a22 = Parameter(name = 'I177a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*complexconjugate(Ru2x5)',
                    texname = '\\text{I177a22}')

I177a25 = Parameter(name = 'I177a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru2x5)',
                    texname = '\\text{I177a25}')

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

I177a41 = Parameter(name = 'I177a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I177a41}')

I177a44 = Parameter(name = 'I177a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I177a44}')

I177a52 = Parameter(name = 'I177a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I177a52}')

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

I178a11 = Parameter(name = 'I178a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*complexconjugate(Ru1x4)',
                    texname = '\\text{I178a11}')

I178a14 = Parameter(name = 'I178a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru1x4)',
                    texname = '\\text{I178a14}')

I178a22 = Parameter(name = 'I178a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*complexconjugate(Ru2x5)',
                    texname = '\\text{I178a22}')

I178a25 = Parameter(name = 'I178a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru2x5)',
                    texname = '\\text{I178a25}')

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

I178a41 = Parameter(name = 'I178a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I178a41}')

I178a44 = Parameter(name = 'I178a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I178a44}')

I178a52 = Parameter(name = 'I178a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I178a52}')

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

I179a11 = Parameter(name = 'I179a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I179a11}')

I179a14 = Parameter(name = 'I179a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I179a14}')

I179a22 = Parameter(name = 'I179a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I179a22}')

I179a25 = Parameter(name = 'I179a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I179a25}')

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

I179a41 = Parameter(name = 'I179a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I179a41}')

I179a44 = Parameter(name = 'I179a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I179a44}')

I179a52 = Parameter(name = 'I179a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I179a52}')

I179a55 = Parameter(name = 'I179a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I179a55}')

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

I18a11 = Parameter(name = 'I18a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I18a11}')

I18a14 = Parameter(name = 'I18a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I18a14}')

I18a22 = Parameter(name = 'I18a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I18a22}')

I18a25 = Parameter(name = 'I18a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I18a25}')

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

I18a41 = Parameter(name = 'I18a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I18a41}')

I18a44 = Parameter(name = 'I18a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I18a44}')

I18a52 = Parameter(name = 'I18a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I18a52}')

I18a55 = Parameter(name = 'I18a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I18a55}')

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

I180a11 = Parameter(name = 'I180a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*yu1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I180a11}')

I180a14 = Parameter(name = 'I180a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*yu1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I180a14}')

I180a22 = Parameter(name = 'I180a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*yu2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I180a22}')

I180a25 = Parameter(name = 'I180a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*yu2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I180a25}')

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

I180a41 = Parameter(name = 'I180a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*yu1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I180a41}')

I180a44 = Parameter(name = 'I180a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*yu1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I180a44}')

I180a52 = Parameter(name = 'I180a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*yu2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I180a52}')

I180a55 = Parameter(name = 'I180a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*yu2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I180a55}')

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

I181a11 = Parameter(name = 'I181a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I181a11}')

I181a14 = Parameter(name = 'I181a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I181a14}')

I181a22 = Parameter(name = 'I181a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I181a22}')

I181a25 = Parameter(name = 'I181a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I181a25}')

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

I181a41 = Parameter(name = 'I181a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I181a41}')

I181a44 = Parameter(name = 'I181a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I181a44}')

I181a52 = Parameter(name = 'I181a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I181a52}')

I181a55 = Parameter(name = 'I181a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I181a55}')

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

I182a11 = Parameter(name = 'I182a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*yu1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I182a11}')

I182a14 = Parameter(name = 'I182a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*yu1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I182a14}')

I182a22 = Parameter(name = 'I182a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*yu2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I182a22}')

I182a25 = Parameter(name = 'I182a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*yu2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I182a25}')

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

I182a41 = Parameter(name = 'I182a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*yu1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I182a41}')

I182a44 = Parameter(name = 'I182a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*yu1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I182a44}')

I182a52 = Parameter(name = 'I182a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*yu2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I182a52}')

I182a55 = Parameter(name = 'I182a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*yu2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I182a55}')

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

I183a11 = Parameter(name = 'I183a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I183a11}')

I183a14 = Parameter(name = 'I183a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I183a14}')

I183a22 = Parameter(name = 'I183a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I183a22}')

I183a25 = Parameter(name = 'I183a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I183a25}')

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

I183a41 = Parameter(name = 'I183a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I183a41}')

I183a44 = Parameter(name = 'I183a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I183a44}')

I183a52 = Parameter(name = 'I183a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I183a52}')

I183a55 = Parameter(name = 'I183a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I183a55}')

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

I184a11 = Parameter(name = 'I184a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*yu1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I184a11}')

I184a14 = Parameter(name = 'I184a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*yu1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I184a14}')

I184a22 = Parameter(name = 'I184a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*yu2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I184a22}')

I184a25 = Parameter(name = 'I184a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*yu2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I184a25}')

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

I184a41 = Parameter(name = 'I184a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*yu1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I184a41}')

I184a44 = Parameter(name = 'I184a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*yu1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I184a44}')

I184a52 = Parameter(name = 'I184a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*yu2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I184a52}')

I184a55 = Parameter(name = 'I184a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*yu2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I184a55}')

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

I185a11 = Parameter(name = 'I185a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I185a11}')

I185a14 = Parameter(name = 'I185a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I185a14}')

I185a22 = Parameter(name = 'I185a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I185a22}')

I185a25 = Parameter(name = 'I185a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I185a25}')

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

I185a41 = Parameter(name = 'I185a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I185a41}')

I185a44 = Parameter(name = 'I185a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I185a44}')

I185a52 = Parameter(name = 'I185a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I185a52}')

I185a55 = Parameter(name = 'I185a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I185a55}')

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

I186a11 = Parameter(name = 'I186a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*yu1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I186a11}')

I186a14 = Parameter(name = 'I186a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*yu1x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I186a14}')

I186a22 = Parameter(name = 'I186a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*yu2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I186a22}')

I186a25 = Parameter(name = 'I186a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*yu2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I186a25}')

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

I186a41 = Parameter(name = 'I186a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*yu1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I186a41}')

I186a44 = Parameter(name = 'I186a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*yu1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I186a44}')

I186a52 = Parameter(name = 'I186a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*yu2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I186a52}')

I186a55 = Parameter(name = 'I186a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*yu2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I186a55}')

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

I187a14 = Parameter(name = 'I187a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*complexconjugate(Rd4x1)',
                    texname = '\\text{I187a14}')

I187a22 = Parameter(name = 'I187a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I187a22}')

I187a25 = Parameter(name = 'I187a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*complexconjugate(Rd5x2)',
                    texname = '\\text{I187a25}')

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

I188a11 = Parameter(name = 'I188a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I188a11}')

I188a14 = Parameter(name = 'I188a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                    texname = '\\text{I188a14}')

I188a22 = Parameter(name = 'I188a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I188a22}')

I188a25 = Parameter(name = 'I188a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                    texname = '\\text{I188a25}')

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

I189a11 = Parameter(name = 'I189a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*yu1x1*complexconjugate(Rd1x1)',
                    texname = '\\text{I189a11}')

I189a14 = Parameter(name = 'I189a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*yu1x1*complexconjugate(Rd4x1)',
                    texname = '\\text{I189a14}')

I189a22 = Parameter(name = 'I189a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*yu2x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I189a22}')

I189a25 = Parameter(name = 'I189a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*yu2x2*complexconjugate(Rd5x2)',
                    texname = '\\text{I189a25}')

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

I19a11 = Parameter(name = 'I19a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*td1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I19a11}')

I19a14 = Parameter(name = 'I19a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*td1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I19a14}')

I19a22 = Parameter(name = 'I19a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*td2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I19a22}')

I19a25 = Parameter(name = 'I19a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*td2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I19a25}')

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

I19a41 = Parameter(name = 'I19a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*td1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I19a41}')

I19a44 = Parameter(name = 'I19a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*td1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I19a44}')

I19a52 = Parameter(name = 'I19a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*td2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I19a52}')

I19a55 = Parameter(name = 'I19a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*td2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I19a55}')

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

I190a14 = Parameter(name = 'I190a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Rl4x1)',
                    texname = '\\text{I190a14}')

I190a22 = Parameter(name = 'I190a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Rl2x2)',
                    texname = '\\text{I190a22}')

I190a25 = Parameter(name = 'I190a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Rl5x2)',
                    texname = '\\text{I190a25}')

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

I191a11 = Parameter(name = 'I191a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                    texname = '\\text{I191a11}')

I191a14 = Parameter(name = 'I191a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                    texname = '\\text{I191a14}')

I191a22 = Parameter(name = 'I191a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                    texname = '\\text{I191a22}')

I191a25 = Parameter(name = 'I191a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                    texname = '\\text{I191a25}')

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

I193a11 = Parameter(name = 'I193a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'ye1x1*complexconjugate(Rn1x1)',
                    texname = '\\text{I193a11}')

I193a22 = Parameter(name = 'I193a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'ye2x2*complexconjugate(Rn2x2)',
                    texname = '\\text{I193a22}')

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

I194a14 = Parameter(name = 'I194a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(CKM1x1)*complexconjugate(Ru4x1)',
                    texname = '\\text{I194a14}')

I194a22 = Parameter(name = 'I194a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I194a22}')

I194a25 = Parameter(name = 'I194a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(CKM2x2)*complexconjugate(Ru5x2)',
                    texname = '\\text{I194a25}')

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

I195a11 = Parameter(name = 'I195a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(CKM1x1)*complexconjugate(Ru1x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I195a11}')

I195a14 = Parameter(name = 'I195a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(CKM1x1)*complexconjugate(Ru4x4)*complexconjugate(yu1x1)',
                    texname = '\\text{I195a14}')

I195a22 = Parameter(name = 'I195a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(CKM2x2)*complexconjugate(Ru2x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I195a22}')

I195a25 = Parameter(name = 'I195a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'complexconjugate(CKM2x2)*complexconjugate(Ru5x5)*complexconjugate(yu2x2)',
                    texname = '\\text{I195a25}')

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

I196a11 = Parameter(name = 'I196a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                    texname = '\\text{I196a11}')

I196a14 = Parameter(name = 'I196a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'yd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)',
                    texname = '\\text{I196a14}')

I196a22 = Parameter(name = 'I196a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I196a22}')

I196a25 = Parameter(name = 'I196a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'yd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)',
                    texname = '\\text{I196a25}')

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

I197a14 = Parameter(name = 'I197a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*complexconjugate(Rd1x1)',
                    texname = '\\text{I197a14}')

I197a22 = Parameter(name = 'I197a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I197a22}')

I197a25 = Parameter(name = 'I197a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*complexconjugate(Rd2x2)',
                    texname = '\\text{I197a25}')

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

I197a41 = Parameter(name = 'I197a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru1x1*complexconjugate(Rd4x1)',
                    texname = '\\text{I197a41}')

I197a44 = Parameter(name = 'I197a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM1x1*Ru4x1*complexconjugate(Rd4x1)',
                    texname = '\\text{I197a44}')

I197a52 = Parameter(name = 'I197a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru2x2*complexconjugate(Rd5x2)',
                    texname = '\\text{I197a52}')

I197a55 = Parameter(name = 'I197a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'CKM2x2*Ru5x2*complexconjugate(Rd5x2)',
                    texname = '\\text{I197a55}')

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

I198a14 = Parameter(name = 'I198a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn1x1*complexconjugate(Rl4x1)',
                    texname = '\\text{I198a14}')

I198a22 = Parameter(name = 'I198a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl2x2)',
                    texname = '\\text{I198a22}')

I198a25 = Parameter(name = 'I198a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rn2x2*complexconjugate(Rl5x2)',
                    texname = '\\text{I198a25}')

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

I199a14 = Parameter(name = 'I199a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)',
                    texname = '\\text{I199a14}')

I199a22 = Parameter(name = 'I199a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I199a22}')

I199a25 = Parameter(name = 'I199a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)',
                    texname = '\\text{I199a25}')

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

I199a41 = Parameter(name = 'I199a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                    texname = '\\text{I199a41}')

I199a44 = Parameter(name = 'I199a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x1*complexconjugate(CKM1x1)*complexconjugate(Ru4x1)',
                    texname = '\\text{I199a44}')

I199a52 = Parameter(name = 'I199a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                    texname = '\\text{I199a52}')

I199a55 = Parameter(name = 'I199a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x2*complexconjugate(CKM2x2)*complexconjugate(Ru5x2)',
                    texname = '\\text{I199a55}')

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

I2a11 = Parameter(name = 'I2a11',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yd1x1*complexconjugate(CKM1x1)',
                  texname = '\\text{I2a11}')

I2a22 = Parameter(name = 'I2a22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yd2x2*complexconjugate(CKM2x2)',
                  texname = '\\text{I2a22}')

I2a33 = Parameter(name = 'I2a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yd3x3*complexconjugate(CKM3x3)',
                  texname = '\\text{I2a33}')

I20a11 = Parameter(name = 'I20a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*yd1x1*complexconjugate(Rd1x1)*complexconjugate(yd1x1)',
                   texname = '\\text{I20a11}')

I20a14 = Parameter(name = 'I20a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*yd1x1*complexconjugate(Rd1x1)*complexconjugate(yd1x1)',
                   texname = '\\text{I20a14}')

I20a22 = Parameter(name = 'I20a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*yd2x2*complexconjugate(Rd2x2)*complexconjugate(yd2x2)',
                   texname = '\\text{I20a22}')

I20a25 = Parameter(name = 'I20a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*yd2x2*complexconjugate(Rd2x2)*complexconjugate(yd2x2)',
                   texname = '\\text{I20a25}')

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

I20a41 = Parameter(name = 'I20a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*yd1x1*complexconjugate(Rd4x1)*complexconjugate(yd1x1)',
                   texname = '\\text{I20a41}')

I20a44 = Parameter(name = 'I20a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*yd1x1*complexconjugate(Rd4x1)*complexconjugate(yd1x1)',
                   texname = '\\text{I20a44}')

I20a52 = Parameter(name = 'I20a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*yd2x2*complexconjugate(Rd5x2)*complexconjugate(yd2x2)',
                   texname = '\\text{I20a52}')

I20a55 = Parameter(name = 'I20a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*yd2x2*complexconjugate(Rd5x2)*complexconjugate(yd2x2)',
                   texname = '\\text{I20a55}')

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

I200a14 = Parameter(name = 'I200a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rl4x1*complexconjugate(Rn1x1)',
                    texname = '\\text{I200a14}')

I200a22 = Parameter(name = 'I200a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rl2x2*complexconjugate(Rn2x2)',
                    texname = '\\text{I200a22}')

I200a25 = Parameter(name = 'I200a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rl5x2*complexconjugate(Rn2x2)',
                    texname = '\\text{I200a25}')

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

I201a14 = Parameter(name = 'I201a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru1x1)',
                    texname = '\\text{I201a14}')

I201a22 = Parameter(name = 'I201a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I201a22}')

I201a25 = Parameter(name = 'I201a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru2x2)',
                    texname = '\\text{I201a25}')

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

I201a41 = Parameter(name = 'I201a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I201a41}')

I201a44 = Parameter(name = 'I201a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x1*complexconjugate(Ru4x1)',
                    texname = '\\text{I201a44}')

I201a52 = Parameter(name = 'I201a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I201a52}')

I201a55 = Parameter(name = 'I201a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x2*complexconjugate(Ru5x2)',
                    texname = '\\text{I201a55}')

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

I202a11 = Parameter(name = 'I202a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'ye1x1',
                    texname = '\\text{I202a11}')

I202a22 = Parameter(name = 'I202a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'ye2x2',
                    texname = '\\text{I202a22}')

I202a33 = Parameter(name = 'I202a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'ye3x3',
                    texname = '\\text{I202a33}')

I203a11 = Parameter(name = 'I203a11',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*complexconjugate(Ru1x4)',
                    texname = '\\text{I203a11}')

I203a14 = Parameter(name = 'I203a14',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru1x4)',
                    texname = '\\text{I203a14}')

I203a22 = Parameter(name = 'I203a22',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*complexconjugate(Ru2x5)',
                    texname = '\\text{I203a22}')

I203a25 = Parameter(name = 'I203a25',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru2x5)',
                    texname = '\\text{I203a25}')

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

I203a41 = Parameter(name = 'I203a41',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru1x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I203a41}')

I203a44 = Parameter(name = 'I203a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I203a44}')

I203a52 = Parameter(name = 'I203a52',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru2x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I203a52}')

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

I21a11 = Parameter(name = 'I21a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*yd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I21a11}')

I21a14 = Parameter(name = 'I21a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*yd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I21a14}')

I21a22 = Parameter(name = 'I21a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*yd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I21a22}')

I21a25 = Parameter(name = 'I21a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*yd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I21a25}')

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

I21a41 = Parameter(name = 'I21a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*yd1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I21a41}')

I21a44 = Parameter(name = 'I21a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*yd1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I21a44}')

I21a52 = Parameter(name = 'I21a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*yd2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I21a52}')

I21a55 = Parameter(name = 'I21a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*yd2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I21a55}')

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

I22a11 = Parameter(name = 'I22a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*yd1x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I22a11}')

I22a14 = Parameter(name = 'I22a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*yd1x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I22a14}')

I22a22 = Parameter(name = 'I22a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*yd2x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I22a22}')

I22a25 = Parameter(name = 'I22a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*yd2x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I22a25}')

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

I22a41 = Parameter(name = 'I22a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*yd1x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I22a41}')

I22a44 = Parameter(name = 'I22a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*yd1x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I22a44}')

I22a52 = Parameter(name = 'I22a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*yd2x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I22a52}')

I22a55 = Parameter(name = 'I22a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*yd2x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I22a55}')

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

I23a11 = Parameter(name = 'I23a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*yd1x1*complexconjugate(Rd1x1)*complexconjugate(yd1x1)',
                   texname = '\\text{I23a11}')

I23a14 = Parameter(name = 'I23a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*yd1x1*complexconjugate(Rd1x1)*complexconjugate(yd1x1)',
                   texname = '\\text{I23a14}')

I23a22 = Parameter(name = 'I23a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*yd2x2*complexconjugate(Rd2x2)*complexconjugate(yd2x2)',
                   texname = '\\text{I23a22}')

I23a25 = Parameter(name = 'I23a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*yd2x2*complexconjugate(Rd2x2)*complexconjugate(yd2x2)',
                   texname = '\\text{I23a25}')

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

I23a41 = Parameter(name = 'I23a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*yd1x1*complexconjugate(Rd4x1)*complexconjugate(yd1x1)',
                   texname = '\\text{I23a41}')

I23a44 = Parameter(name = 'I23a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*yd1x1*complexconjugate(Rd4x1)*complexconjugate(yd1x1)',
                   texname = '\\text{I23a44}')

I23a52 = Parameter(name = 'I23a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*yd2x2*complexconjugate(Rd5x2)*complexconjugate(yd2x2)',
                   texname = '\\text{I23a52}')

I23a55 = Parameter(name = 'I23a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*yd2x2*complexconjugate(Rd5x2)*complexconjugate(yd2x2)',
                   texname = '\\text{I23a55}')

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

I24a11 = Parameter(name = 'I24a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*yd1x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I24a11}')

I24a14 = Parameter(name = 'I24a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*yd1x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I24a14}')

I24a22 = Parameter(name = 'I24a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*yd2x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I24a22}')

I24a25 = Parameter(name = 'I24a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*yd2x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I24a25}')

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

I24a41 = Parameter(name = 'I24a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*yd1x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I24a41}')

I24a44 = Parameter(name = 'I24a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*yd1x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I24a44}')

I24a52 = Parameter(name = 'I24a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*yd2x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I24a52}')

I24a55 = Parameter(name = 'I24a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*yd2x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I24a55}')

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

I25a11 = Parameter(name = 'I25a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM1x1*Rd1x1*yu1x1*complexconjugate(CKM1x1)*complexconjugate(Rd1x1)*complexconjugate(yu1x1)',
                   texname = '\\text{I25a11}')

I25a14 = Parameter(name = 'I25a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM1x1*Rd4x1*yu1x1*complexconjugate(CKM1x1)*complexconjugate(Rd1x1)*complexconjugate(yu1x1)',
                   texname = '\\text{I25a14}')

I25a22 = Parameter(name = 'I25a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM2x2*Rd2x2*yu2x2*complexconjugate(CKM2x2)*complexconjugate(Rd2x2)*complexconjugate(yu2x2)',
                   texname = '\\text{I25a22}')

I25a25 = Parameter(name = 'I25a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM2x2*Rd5x2*yu2x2*complexconjugate(CKM2x2)*complexconjugate(Rd2x2)*complexconjugate(yu2x2)',
                   texname = '\\text{I25a25}')

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

I25a41 = Parameter(name = 'I25a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM1x1*Rd1x1*yu1x1*complexconjugate(CKM1x1)*complexconjugate(Rd4x1)*complexconjugate(yu1x1)',
                   texname = '\\text{I25a41}')

I25a44 = Parameter(name = 'I25a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM1x1*Rd4x1*yu1x1*complexconjugate(CKM1x1)*complexconjugate(Rd4x1)*complexconjugate(yu1x1)',
                   texname = '\\text{I25a44}')

I25a52 = Parameter(name = 'I25a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM2x2*Rd2x2*yu2x2*complexconjugate(CKM2x2)*complexconjugate(Rd5x2)*complexconjugate(yu2x2)',
                   texname = '\\text{I25a52}')

I25a55 = Parameter(name = 'I25a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM2x2*Rd5x2*yu2x2*complexconjugate(CKM2x2)*complexconjugate(Rd5x2)*complexconjugate(yu2x2)',
                   texname = '\\text{I25a55}')

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

I26a14 = Parameter(name = 'I26a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I26a14}')

I26a22 = Parameter(name = 'I26a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I26a22}')

I26a25 = Parameter(name = 'I26a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I26a25}')

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

I26a41 = Parameter(name = 'I26a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I26a41}')

I26a44 = Parameter(name = 'I26a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I26a44}')

I26a52 = Parameter(name = 'I26a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I26a52}')

I26a55 = Parameter(name = 'I26a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I26a55}')

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

I27a14 = Parameter(name = 'I27a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I27a14}')

I27a22 = Parameter(name = 'I27a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I27a22}')

I27a25 = Parameter(name = 'I27a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I27a25}')

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

I27a41 = Parameter(name = 'I27a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I27a41}')

I27a44 = Parameter(name = 'I27a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I27a44}')

I27a52 = Parameter(name = 'I27a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I27a52}')

I27a55 = Parameter(name = 'I27a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I27a55}')

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

I28a14 = Parameter(name = 'I28a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I28a14}')

I28a22 = Parameter(name = 'I28a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I28a22}')

I28a25 = Parameter(name = 'I28a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I28a25}')

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

I28a41 = Parameter(name = 'I28a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I28a41}')

I28a44 = Parameter(name = 'I28a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I28a44}')

I28a52 = Parameter(name = 'I28a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I28a52}')

I28a55 = Parameter(name = 'I28a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I28a55}')

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

I29a14 = Parameter(name = 'I29a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I29a14}')

I29a22 = Parameter(name = 'I29a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I29a22}')

I29a25 = Parameter(name = 'I29a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I29a25}')

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

I29a41 = Parameter(name = 'I29a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I29a41}')

I29a44 = Parameter(name = 'I29a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I29a44}')

I29a52 = Parameter(name = 'I29a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I29a52}')

I29a55 = Parameter(name = 'I29a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I29a55}')

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

I3a11 = Parameter(name = 'I3a11',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM1x1*complexconjugate(yd1x1)',
                  texname = '\\text{I3a11}')

I3a22 = Parameter(name = 'I3a22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM2x2*complexconjugate(yd2x2)',
                  texname = '\\text{I3a22}')

I3a33 = Parameter(name = 'I3a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM3x3*complexconjugate(yd3x3)',
                  texname = '\\text{I3a33}')

I30a11 = Parameter(name = 'I30a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*complexconjugate(Rd1x4)',
                   texname = '\\text{I30a11}')

I30a14 = Parameter(name = 'I30a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd1x4)',
                   texname = '\\text{I30a14}')

I30a22 = Parameter(name = 'I30a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*complexconjugate(Rd2x5)',
                   texname = '\\text{I30a22}')

I30a25 = Parameter(name = 'I30a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*complexconjugate(Rd2x5)',
                   texname = '\\text{I30a25}')

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

I30a41 = Parameter(name = 'I30a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I30a41}')

I30a44 = Parameter(name = 'I30a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I30a44}')

I30a52 = Parameter(name = 'I30a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*complexconjugate(Rd5x5)',
                   texname = '\\text{I30a52}')

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

I31a14 = Parameter(name = 'I31a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I31a14}')

I31a22 = Parameter(name = 'I31a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I31a22}')

I31a25 = Parameter(name = 'I31a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I31a25}')

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

I31a41 = Parameter(name = 'I31a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I31a41}')

I31a44 = Parameter(name = 'I31a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I31a44}')

I31a52 = Parameter(name = 'I31a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I31a52}')

I31a55 = Parameter(name = 'I31a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I31a55}')

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

I32a11 = Parameter(name = 'I32a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*complexconjugate(Rd1x4)',
                   texname = '\\text{I32a11}')

I32a14 = Parameter(name = 'I32a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd1x4)',
                   texname = '\\text{I32a14}')

I32a22 = Parameter(name = 'I32a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*complexconjugate(Rd2x5)',
                   texname = '\\text{I32a22}')

I32a25 = Parameter(name = 'I32a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*complexconjugate(Rd2x5)',
                   texname = '\\text{I32a25}')

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

I32a41 = Parameter(name = 'I32a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I32a41}')

I32a44 = Parameter(name = 'I32a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I32a44}')

I32a52 = Parameter(name = 'I32a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*complexconjugate(Rd5x5)',
                   texname = '\\text{I32a52}')

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

I33a11 = Parameter(name = 'I33a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*complexconjugate(Rd1x4)',
                   texname = '\\text{I33a11}')

I33a14 = Parameter(name = 'I33a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd1x4)',
                   texname = '\\text{I33a14}')

I33a22 = Parameter(name = 'I33a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*complexconjugate(Rd2x5)',
                   texname = '\\text{I33a22}')

I33a25 = Parameter(name = 'I33a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*complexconjugate(Rd2x5)',
                   texname = '\\text{I33a25}')

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

I33a41 = Parameter(name = 'I33a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I33a41}')

I33a44 = Parameter(name = 'I33a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I33a44}')

I33a52 = Parameter(name = 'I33a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*complexconjugate(Rd5x5)',
                   texname = '\\text{I33a52}')

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

I34a11 = Parameter(name = 'I34a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*complexconjugate(Rd1x4)',
                   texname = '\\text{I34a11}')

I34a14 = Parameter(name = 'I34a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd1x4)',
                   texname = '\\text{I34a14}')

I34a22 = Parameter(name = 'I34a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*complexconjugate(Rd2x5)',
                   texname = '\\text{I34a22}')

I34a25 = Parameter(name = 'I34a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*complexconjugate(Rd2x5)',
                   texname = '\\text{I34a25}')

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

I34a41 = Parameter(name = 'I34a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I34a41}')

I34a44 = Parameter(name = 'I34a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I34a44}')

I34a52 = Parameter(name = 'I34a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*complexconjugate(Rd5x5)',
                   texname = '\\text{I34a52}')

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

I35a11 = Parameter(name = 'I35a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*complexconjugate(Rd1x4)',
                   texname = '\\text{I35a11}')

I35a14 = Parameter(name = 'I35a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd1x4)',
                   texname = '\\text{I35a14}')

I35a22 = Parameter(name = 'I35a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*complexconjugate(Rd2x5)',
                   texname = '\\text{I35a22}')

I35a25 = Parameter(name = 'I35a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*complexconjugate(Rd2x5)',
                   texname = '\\text{I35a25}')

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

I35a41 = Parameter(name = 'I35a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I35a41}')

I35a44 = Parameter(name = 'I35a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I35a44}')

I35a52 = Parameter(name = 'I35a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*complexconjugate(Rd5x5)',
                   texname = '\\text{I35a52}')

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

I36a11 = Parameter(name = 'I36a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*complexconjugate(Rd1x4)',
                   texname = '\\text{I36a11}')

I36a14 = Parameter(name = 'I36a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd1x4)',
                   texname = '\\text{I36a14}')

I36a22 = Parameter(name = 'I36a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*complexconjugate(Rd2x5)',
                   texname = '\\text{I36a22}')

I36a25 = Parameter(name = 'I36a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*complexconjugate(Rd2x5)',
                   texname = '\\text{I36a25}')

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

I36a41 = Parameter(name = 'I36a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I36a41}')

I36a44 = Parameter(name = 'I36a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I36a44}')

I36a52 = Parameter(name = 'I36a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*complexconjugate(Rd5x5)',
                   texname = '\\text{I36a52}')

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

I37a11 = Parameter(name = 'I37a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I37a11}')

I37a14 = Parameter(name = 'I37a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I37a14}')

I37a22 = Parameter(name = 'I37a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I37a22}')

I37a25 = Parameter(name = 'I37a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I37a25}')

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

I37a41 = Parameter(name = 'I37a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I37a41}')

I37a44 = Parameter(name = 'I37a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I37a44}')

I37a52 = Parameter(name = 'I37a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I37a52}')

I37a55 = Parameter(name = 'I37a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I37a55}')

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

I38a11 = Parameter(name = 'I38a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*yd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I38a11}')

I38a14 = Parameter(name = 'I38a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*yd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I38a14}')

I38a22 = Parameter(name = 'I38a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*yd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I38a22}')

I38a25 = Parameter(name = 'I38a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*yd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I38a25}')

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

I38a41 = Parameter(name = 'I38a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*yd1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I38a41}')

I38a44 = Parameter(name = 'I38a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*yd1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I38a44}')

I38a52 = Parameter(name = 'I38a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*yd2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I38a52}')

I38a55 = Parameter(name = 'I38a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*yd2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I38a55}')

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

I39a11 = Parameter(name = 'I39a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I39a11}')

I39a14 = Parameter(name = 'I39a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I39a14}')

I39a22 = Parameter(name = 'I39a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I39a22}')

I39a25 = Parameter(name = 'I39a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I39a25}')

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

I39a41 = Parameter(name = 'I39a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I39a41}')

I39a44 = Parameter(name = 'I39a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I39a44}')

I39a52 = Parameter(name = 'I39a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I39a52}')

I39a55 = Parameter(name = 'I39a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I39a55}')

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

I4a11 = Parameter(name = 'I4a11',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM1x1*yu1x1',
                  texname = '\\text{I4a11}')

I4a22 = Parameter(name = 'I4a22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM2x2*yu2x2',
                  texname = '\\text{I4a22}')

I4a33 = Parameter(name = 'I4a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM3x3*yu3x3',
                  texname = '\\text{I4a33}')

I40a11 = Parameter(name = 'I40a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*yd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I40a11}')

I40a14 = Parameter(name = 'I40a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*yd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I40a14}')

I40a22 = Parameter(name = 'I40a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*yd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I40a22}')

I40a25 = Parameter(name = 'I40a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*yd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I40a25}')

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

I40a41 = Parameter(name = 'I40a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*yd1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I40a41}')

I40a44 = Parameter(name = 'I40a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*yd1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I40a44}')

I40a52 = Parameter(name = 'I40a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*yd2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I40a52}')

I40a55 = Parameter(name = 'I40a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*yd2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I40a55}')

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

I41a11 = Parameter(name = 'I41a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I41a11}')

I41a14 = Parameter(name = 'I41a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I41a14}')

I41a22 = Parameter(name = 'I41a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I41a22}')

I41a25 = Parameter(name = 'I41a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I41a25}')

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

I41a41 = Parameter(name = 'I41a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I41a41}')

I41a44 = Parameter(name = 'I41a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I41a44}')

I41a52 = Parameter(name = 'I41a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I41a52}')

I41a55 = Parameter(name = 'I41a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I41a55}')

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

I42a11 = Parameter(name = 'I42a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*yd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I42a11}')

I42a14 = Parameter(name = 'I42a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*yd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I42a14}')

I42a22 = Parameter(name = 'I42a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*yd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I42a22}')

I42a25 = Parameter(name = 'I42a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*yd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I42a25}')

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

I42a41 = Parameter(name = 'I42a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*yd1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I42a41}')

I42a44 = Parameter(name = 'I42a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*yd1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I42a44}')

I42a52 = Parameter(name = 'I42a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*yd2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I42a52}')

I42a55 = Parameter(name = 'I42a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*yd2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I42a55}')

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

I43a11 = Parameter(name = 'I43a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I43a11}')

I43a14 = Parameter(name = 'I43a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I43a14}')

I43a22 = Parameter(name = 'I43a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I43a22}')

I43a25 = Parameter(name = 'I43a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I43a25}')

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

I43a41 = Parameter(name = 'I43a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I43a41}')

I43a44 = Parameter(name = 'I43a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I43a44}')

I43a52 = Parameter(name = 'I43a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I43a52}')

I43a55 = Parameter(name = 'I43a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I43a55}')

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

I44a11 = Parameter(name = 'I44a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*yd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I44a11}')

I44a14 = Parameter(name = 'I44a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*yd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I44a14}')

I44a22 = Parameter(name = 'I44a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*yd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I44a22}')

I44a25 = Parameter(name = 'I44a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*yd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I44a25}')

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

I44a41 = Parameter(name = 'I44a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*yd1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I44a41}')

I44a44 = Parameter(name = 'I44a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*yd1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I44a44}')

I44a52 = Parameter(name = 'I44a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*yd2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I44a52}')

I44a55 = Parameter(name = 'I44a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*yd2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I44a55}')

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

I45a11 = Parameter(name = 'I45a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I45a11}')

I45a14 = Parameter(name = 'I45a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I45a14}')

I45a22 = Parameter(name = 'I45a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I45a22}')

I45a25 = Parameter(name = 'I45a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I45a25}')

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

I46a11 = Parameter(name = 'I46a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'ye1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I46a11}')

I46a14 = Parameter(name = 'I46a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'ye1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I46a14}')

I46a22 = Parameter(name = 'I46a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'ye2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I46a22}')

I46a25 = Parameter(name = 'I46a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'ye2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I46a25}')

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

I47a14 = Parameter(name = 'I47a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I47a14}')

I47a22 = Parameter(name = 'I47a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I47a22}')

I47a25 = Parameter(name = 'I47a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I47a25}')

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

I47a41 = Parameter(name = 'I47a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I47a41}')

I47a44 = Parameter(name = 'I47a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I47a44}')

I47a52 = Parameter(name = 'I47a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I47a52}')

I47a55 = Parameter(name = 'I47a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I47a55}')

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

I48a11 = Parameter(name = 'I48a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*complexconjugate(Rl1x4)',
                   texname = '\\text{I48a11}')

I48a14 = Parameter(name = 'I48a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl1x4)',
                   texname = '\\text{I48a14}')

I48a22 = Parameter(name = 'I48a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*complexconjugate(Rl2x5)',
                   texname = '\\text{I48a22}')

I48a25 = Parameter(name = 'I48a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl2x5)',
                   texname = '\\text{I48a25}')

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

I48a41 = Parameter(name = 'I48a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I48a41}')

I48a44 = Parameter(name = 'I48a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I48a44}')

I48a52 = Parameter(name = 'I48a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I48a52}')

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

I49a11 = Parameter(name = 'I49a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(ye1x1)',
                   texname = '\\text{I49a11}')

I49a14 = Parameter(name = 'I49a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(ye1x1)',
                   texname = '\\text{I49a14}')

I49a22 = Parameter(name = 'I49a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(ye2x2)',
                   texname = '\\text{I49a22}')

I49a25 = Parameter(name = 'I49a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(ye2x2)',
                   texname = '\\text{I49a25}')

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

I5a11 = Parameter(name = 'I5a11',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complexconjugate(ye1x1)',
                  texname = '\\text{I5a11}')

I5a22 = Parameter(name = 'I5a22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complexconjugate(ye2x2)',
                  texname = '\\text{I5a22}')

I5a33 = Parameter(name = 'I5a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complexconjugate(ye3x3)',
                  texname = '\\text{I5a33}')

I50a11 = Parameter(name = 'I50a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*ye1x1',
                   texname = '\\text{I50a11}')

I50a14 = Parameter(name = 'I50a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*ye1x1',
                   texname = '\\text{I50a14}')

I50a22 = Parameter(name = 'I50a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*ye2x2',
                   texname = '\\text{I50a22}')

I50a25 = Parameter(name = 'I50a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*ye2x2',
                   texname = '\\text{I50a25}')

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

I51a14 = Parameter(name = 'I51a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1',
                   texname = '\\text{I51a14}')

I51a22 = Parameter(name = 'I51a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2',
                   texname = '\\text{I51a22}')

I51a25 = Parameter(name = 'I51a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2',
                   texname = '\\text{I51a25}')

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

I52a11 = Parameter(name = 'I52a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*ye1x1',
                   texname = '\\text{I52a11}')

I52a14 = Parameter(name = 'I52a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*ye1x1',
                   texname = '\\text{I52a14}')

I52a22 = Parameter(name = 'I52a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*ye2x2',
                   texname = '\\text{I52a22}')

I52a25 = Parameter(name = 'I52a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*ye2x2',
                   texname = '\\text{I52a25}')

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

I53a14 = Parameter(name = 'I53a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I53a14}')

I53a22 = Parameter(name = 'I53a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I53a22}')

I53a25 = Parameter(name = 'I53a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I53a25}')

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

I53a41 = Parameter(name = 'I53a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I53a41}')

I53a44 = Parameter(name = 'I53a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I53a44}')

I53a52 = Parameter(name = 'I53a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I53a52}')

I53a55 = Parameter(name = 'I53a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I53a55}')

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

I54a11 = Parameter(name = 'I54a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*complexconjugate(Rl1x4)',
                   texname = '\\text{I54a11}')

I54a14 = Parameter(name = 'I54a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl1x4)',
                   texname = '\\text{I54a14}')

I54a22 = Parameter(name = 'I54a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*complexconjugate(Rl2x5)',
                   texname = '\\text{I54a22}')

I54a25 = Parameter(name = 'I54a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl2x5)',
                   texname = '\\text{I54a25}')

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

I54a41 = Parameter(name = 'I54a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I54a41}')

I54a44 = Parameter(name = 'I54a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I54a44}')

I54a52 = Parameter(name = 'I54a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I54a52}')

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

I55a11 = Parameter(name = 'I55a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl1x4)*complexconjugate(te1x1)',
                   texname = '\\text{I55a11}')

I55a14 = Parameter(name = 'I55a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl1x4)*complexconjugate(te1x1)',
                   texname = '\\text{I55a14}')

I55a22 = Parameter(name = 'I55a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x5)*complexconjugate(te2x2)',
                   texname = '\\text{I55a22}')

I55a25 = Parameter(name = 'I55a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl2x5)*complexconjugate(te2x2)',
                   texname = '\\text{I55a25}')

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

I55a41 = Parameter(name = 'I55a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl4x4)*complexconjugate(te1x1)',
                   texname = '\\text{I55a41}')

I55a44 = Parameter(name = 'I55a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl4x4)*complexconjugate(te1x1)',
                   texname = '\\text{I55a44}')

I55a52 = Parameter(name = 'I55a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl5x5)*complexconjugate(te2x2)',
                   texname = '\\text{I55a52}')

I55a55 = Parameter(name = 'I55a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl5x5)*complexconjugate(te2x2)',
                   texname = '\\text{I55a55}')

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

I56a11 = Parameter(name = 'I56a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I56a11}')

I56a14 = Parameter(name = 'I56a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I56a14}')

I56a22 = Parameter(name = 'I56a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I56a22}')

I56a25 = Parameter(name = 'I56a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I56a25}')

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

I56a41 = Parameter(name = 'I56a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I56a41}')

I56a44 = Parameter(name = 'I56a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I56a44}')

I56a52 = Parameter(name = 'I56a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I56a52}')

I56a55 = Parameter(name = 'I56a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I56a55}')

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

I57a11 = Parameter(name = 'I57a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*te1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I57a11}')

I57a14 = Parameter(name = 'I57a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*te1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I57a14}')

I57a22 = Parameter(name = 'I57a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*te2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I57a22}')

I57a25 = Parameter(name = 'I57a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*te2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I57a25}')

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

I57a41 = Parameter(name = 'I57a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*te1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I57a41}')

I57a44 = Parameter(name = 'I57a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*te1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I57a44}')

I57a52 = Parameter(name = 'I57a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*te2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I57a52}')

I57a55 = Parameter(name = 'I57a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*te2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I57a55}')

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

I58a11 = Parameter(name = 'I58a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*ye1x1*complexconjugate(Rl1x1)*complexconjugate(ye1x1)',
                   texname = '\\text{I58a11}')

I58a14 = Parameter(name = 'I58a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*ye1x1*complexconjugate(Rl1x1)*complexconjugate(ye1x1)',
                   texname = '\\text{I58a14}')

I58a22 = Parameter(name = 'I58a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*ye2x2*complexconjugate(Rl2x2)*complexconjugate(ye2x2)',
                   texname = '\\text{I58a22}')

I58a25 = Parameter(name = 'I58a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*ye2x2*complexconjugate(Rl2x2)*complexconjugate(ye2x2)',
                   texname = '\\text{I58a25}')

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

I58a41 = Parameter(name = 'I58a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*ye1x1*complexconjugate(Rl4x1)*complexconjugate(ye1x1)',
                   texname = '\\text{I58a41}')

I58a44 = Parameter(name = 'I58a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*ye1x1*complexconjugate(Rl4x1)*complexconjugate(ye1x1)',
                   texname = '\\text{I58a44}')

I58a52 = Parameter(name = 'I58a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*ye2x2*complexconjugate(Rl5x2)*complexconjugate(ye2x2)',
                   texname = '\\text{I58a52}')

I58a55 = Parameter(name = 'I58a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*ye2x2*complexconjugate(Rl5x2)*complexconjugate(ye2x2)',
                   texname = '\\text{I58a55}')

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

I59a11 = Parameter(name = 'I59a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*ye1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I59a11}')

I59a14 = Parameter(name = 'I59a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*ye1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I59a14}')

I59a22 = Parameter(name = 'I59a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*ye2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I59a22}')

I59a25 = Parameter(name = 'I59a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*ye2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I59a25}')

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

I59a41 = Parameter(name = 'I59a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*ye1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I59a41}')

I59a44 = Parameter(name = 'I59a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*ye1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I59a44}')

I59a52 = Parameter(name = 'I59a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*ye2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I59a52}')

I59a55 = Parameter(name = 'I59a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*ye2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I59a55}')

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

I6a11 = Parameter(name = 'I6a11',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                  texname = '\\text{I6a11}')

I6a14 = Parameter(name = 'I6a14',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                  texname = '\\text{I6a14}')

I6a22 = Parameter(name = 'I6a22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                  texname = '\\text{I6a22}')

I6a25 = Parameter(name = 'I6a25',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                  texname = '\\text{I6a25}')

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

I60a11 = Parameter(name = 'I60a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*ye1x1*complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I60a11}')

I60a14 = Parameter(name = 'I60a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*ye1x1*complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I60a14}')

I60a22 = Parameter(name = 'I60a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*ye2x2*complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I60a22}')

I60a25 = Parameter(name = 'I60a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*ye2x2*complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I60a25}')

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

I60a41 = Parameter(name = 'I60a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*ye1x1*complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I60a41}')

I60a44 = Parameter(name = 'I60a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*ye1x1*complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I60a44}')

I60a52 = Parameter(name = 'I60a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*ye2x2*complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I60a52}')

I60a55 = Parameter(name = 'I60a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*ye2x2*complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I60a55}')

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

I61a11 = Parameter(name = 'I61a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*ye1x1*complexconjugate(Rl1x1)*complexconjugate(ye1x1)',
                   texname = '\\text{I61a11}')

I61a14 = Parameter(name = 'I61a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*ye1x1*complexconjugate(Rl1x1)*complexconjugate(ye1x1)',
                   texname = '\\text{I61a14}')

I61a22 = Parameter(name = 'I61a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*ye2x2*complexconjugate(Rl2x2)*complexconjugate(ye2x2)',
                   texname = '\\text{I61a22}')

I61a25 = Parameter(name = 'I61a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*ye2x2*complexconjugate(Rl2x2)*complexconjugate(ye2x2)',
                   texname = '\\text{I61a25}')

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

I61a41 = Parameter(name = 'I61a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*ye1x1*complexconjugate(Rl4x1)*complexconjugate(ye1x1)',
                   texname = '\\text{I61a41}')

I61a44 = Parameter(name = 'I61a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*ye1x1*complexconjugate(Rl4x1)*complexconjugate(ye1x1)',
                   texname = '\\text{I61a44}')

I61a52 = Parameter(name = 'I61a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*ye2x2*complexconjugate(Rl5x2)*complexconjugate(ye2x2)',
                   texname = '\\text{I61a52}')

I61a55 = Parameter(name = 'I61a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*ye2x2*complexconjugate(Rl5x2)*complexconjugate(ye2x2)',
                   texname = '\\text{I61a55}')

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

I62a11 = Parameter(name = 'I62a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*ye1x1*complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I62a11}')

I62a14 = Parameter(name = 'I62a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*ye1x1*complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I62a14}')

I62a22 = Parameter(name = 'I62a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*ye2x2*complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I62a22}')

I62a25 = Parameter(name = 'I62a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*ye2x2*complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I62a25}')

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

I62a41 = Parameter(name = 'I62a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*ye1x1*complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I62a41}')

I62a44 = Parameter(name = 'I62a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*ye1x1*complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I62a44}')

I62a52 = Parameter(name = 'I62a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*ye2x2*complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I62a52}')

I62a55 = Parameter(name = 'I62a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*ye2x2*complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I62a55}')

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

I63a14 = Parameter(name = 'I63a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I63a14}')

I63a22 = Parameter(name = 'I63a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I63a22}')

I63a25 = Parameter(name = 'I63a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I63a25}')

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

I63a41 = Parameter(name = 'I63a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I63a41}')

I63a44 = Parameter(name = 'I63a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I63a44}')

I63a52 = Parameter(name = 'I63a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I63a52}')

I63a55 = Parameter(name = 'I63a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I63a55}')

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

I64a14 = Parameter(name = 'I64a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I64a14}')

I64a22 = Parameter(name = 'I64a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I64a22}')

I64a25 = Parameter(name = 'I64a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I64a25}')

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

I64a41 = Parameter(name = 'I64a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I64a41}')

I64a44 = Parameter(name = 'I64a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I64a44}')

I64a52 = Parameter(name = 'I64a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I64a52}')

I64a55 = Parameter(name = 'I64a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I64a55}')

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

I65a11 = Parameter(name = 'I65a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*complexconjugate(Rd1x4)',
                   texname = '\\text{I65a11}')

I65a14 = Parameter(name = 'I65a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd1x4)',
                   texname = '\\text{I65a14}')

I65a22 = Parameter(name = 'I65a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*complexconjugate(Rd2x5)',
                   texname = '\\text{I65a22}')

I65a25 = Parameter(name = 'I65a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*complexconjugate(Rd2x5)',
                   texname = '\\text{I65a25}')

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

I65a41 = Parameter(name = 'I65a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I65a41}')

I65a44 = Parameter(name = 'I65a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I65a44}')

I65a52 = Parameter(name = 'I65a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*complexconjugate(Rd5x5)',
                   texname = '\\text{I65a52}')

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

I66a11 = Parameter(name = 'I66a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*complexconjugate(Rl1x4)',
                   texname = '\\text{I66a11}')

I66a14 = Parameter(name = 'I66a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl1x4)',
                   texname = '\\text{I66a14}')

I66a22 = Parameter(name = 'I66a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*complexconjugate(Rl2x5)',
                   texname = '\\text{I66a22}')

I66a25 = Parameter(name = 'I66a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl2x5)',
                   texname = '\\text{I66a25}')

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

I66a41 = Parameter(name = 'I66a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I66a41}')

I66a44 = Parameter(name = 'I66a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I66a44}')

I66a52 = Parameter(name = 'I66a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I66a52}')

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

I67a11 = Parameter(name = 'I67a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I67a11}')

I67a14 = Parameter(name = 'I67a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I67a14}')

I67a22 = Parameter(name = 'I67a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I67a22}')

I67a25 = Parameter(name = 'I67a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I67a25}')

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

I67a41 = Parameter(name = 'I67a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I67a41}')

I67a44 = Parameter(name = 'I67a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I67a44}')

I67a52 = Parameter(name = 'I67a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I67a52}')

I67a55 = Parameter(name = 'I67a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I67a55}')

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

I68a11 = Parameter(name = 'I68a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*yd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I68a11}')

I68a14 = Parameter(name = 'I68a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*yd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I68a14}')

I68a22 = Parameter(name = 'I68a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*yd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I68a22}')

I68a25 = Parameter(name = 'I68a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*yd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I68a25}')

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

I68a41 = Parameter(name = 'I68a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x4*yd1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I68a41}')

I68a44 = Parameter(name = 'I68a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*yd1x1*complexconjugate(Rd4x1)',
                   texname = '\\text{I68a44}')

I68a52 = Parameter(name = 'I68a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x5*yd2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I68a52}')

I68a55 = Parameter(name = 'I68a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*yd2x2*complexconjugate(Rd5x2)',
                   texname = '\\text{I68a55}')

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

I69a11 = Parameter(name = 'I69a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I69a11}')

I69a14 = Parameter(name = 'I69a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd1x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I69a14}')

I69a22 = Parameter(name = 'I69a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I69a22}')

I69a25 = Parameter(name = 'I69a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd2x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I69a25}')

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

I69a41 = Parameter(name = 'I69a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I69a41}')

I69a44 = Parameter(name = 'I69a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x1*complexconjugate(Rd4x4)*complexconjugate(yd1x1)',
                   texname = '\\text{I69a44}')

I69a52 = Parameter(name = 'I69a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I69a52}')

I69a55 = Parameter(name = 'I69a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x2*complexconjugate(Rd5x5)*complexconjugate(yd2x2)',
                   texname = '\\text{I69a55}')

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

I7a11 = Parameter(name = 'I7a11',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yd1x1*complexconjugate(Rd1x1)',
                  texname = '\\text{I7a11}')

I7a14 = Parameter(name = 'I7a14',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yd1x1*complexconjugate(Rd4x1)',
                  texname = '\\text{I7a14}')

I7a22 = Parameter(name = 'I7a22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yd2x2*complexconjugate(Rd2x2)',
                  texname = '\\text{I7a22}')

I7a25 = Parameter(name = 'I7a25',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yd2x2*complexconjugate(Rd5x2)',
                  texname = '\\text{I7a25}')

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

I70a11 = Parameter(name = 'I70a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*ye1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I70a11}')

I70a14 = Parameter(name = 'I70a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*ye1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I70a14}')

I70a22 = Parameter(name = 'I70a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*ye2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I70a22}')

I70a25 = Parameter(name = 'I70a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*ye2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I70a25}')

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

I70a41 = Parameter(name = 'I70a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*ye1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I70a41}')

I70a44 = Parameter(name = 'I70a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*ye1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I70a44}')

I70a52 = Parameter(name = 'I70a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*ye2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I70a52}')

I70a55 = Parameter(name = 'I70a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*ye2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I70a55}')

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

I71a14 = Parameter(name = 'I71a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I71a14}')

I71a22 = Parameter(name = 'I71a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I71a22}')

I71a25 = Parameter(name = 'I71a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I71a25}')

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

I71a41 = Parameter(name = 'I71a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I71a41}')

I71a44 = Parameter(name = 'I71a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I71a44}')

I71a52 = Parameter(name = 'I71a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I71a52}')

I71a55 = Parameter(name = 'I71a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I71a55}')

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

I72a14 = Parameter(name = 'I72a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I72a14}')

I72a22 = Parameter(name = 'I72a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I72a22}')

I72a25 = Parameter(name = 'I72a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I72a25}')

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

I72a41 = Parameter(name = 'I72a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I72a41}')

I72a44 = Parameter(name = 'I72a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I72a44}')

I72a52 = Parameter(name = 'I72a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I72a52}')

I72a55 = Parameter(name = 'I72a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I72a55}')

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

I73a14 = Parameter(name = 'I73a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I73a14}')

I73a22 = Parameter(name = 'I73a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I73a22}')

I73a25 = Parameter(name = 'I73a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I73a25}')

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

I73a41 = Parameter(name = 'I73a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I73a41}')

I73a44 = Parameter(name = 'I73a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I73a44}')

I73a52 = Parameter(name = 'I73a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I73a52}')

I73a55 = Parameter(name = 'I73a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I73a55}')

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

I74a14 = Parameter(name = 'I74a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I74a14}')

I74a22 = Parameter(name = 'I74a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I74a22}')

I74a25 = Parameter(name = 'I74a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I74a25}')

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

I74a41 = Parameter(name = 'I74a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I74a41}')

I74a44 = Parameter(name = 'I74a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I74a44}')

I74a52 = Parameter(name = 'I74a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I74a52}')

I74a55 = Parameter(name = 'I74a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I74a55}')

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

I75a14 = Parameter(name = 'I75a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I75a14}')

I75a22 = Parameter(name = 'I75a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I75a22}')

I75a25 = Parameter(name = 'I75a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I75a25}')

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

I75a41 = Parameter(name = 'I75a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I75a41}')

I75a44 = Parameter(name = 'I75a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I75a44}')

I75a52 = Parameter(name = 'I75a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I75a52}')

I75a55 = Parameter(name = 'I75a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I75a55}')

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

I76a11 = Parameter(name = 'I76a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*complexconjugate(Rl1x4)',
                   texname = '\\text{I76a11}')

I76a14 = Parameter(name = 'I76a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl1x4)',
                   texname = '\\text{I76a14}')

I76a22 = Parameter(name = 'I76a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*complexconjugate(Rl2x5)',
                   texname = '\\text{I76a22}')

I76a25 = Parameter(name = 'I76a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl2x5)',
                   texname = '\\text{I76a25}')

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

I76a41 = Parameter(name = 'I76a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I76a41}')

I76a44 = Parameter(name = 'I76a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I76a44}')

I76a52 = Parameter(name = 'I76a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I76a52}')

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

I77a11 = Parameter(name = 'I77a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*complexconjugate(Rl1x4)',
                   texname = '\\text{I77a11}')

I77a14 = Parameter(name = 'I77a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl1x4)',
                   texname = '\\text{I77a14}')

I77a22 = Parameter(name = 'I77a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*complexconjugate(Rl2x5)',
                   texname = '\\text{I77a22}')

I77a25 = Parameter(name = 'I77a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl2x5)',
                   texname = '\\text{I77a25}')

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

I77a41 = Parameter(name = 'I77a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I77a41}')

I77a44 = Parameter(name = 'I77a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I77a44}')

I77a52 = Parameter(name = 'I77a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I77a52}')

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

I78a11 = Parameter(name = 'I78a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*complexconjugate(Rl1x4)',
                   texname = '\\text{I78a11}')

I78a14 = Parameter(name = 'I78a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl1x4)',
                   texname = '\\text{I78a14}')

I78a22 = Parameter(name = 'I78a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*complexconjugate(Rl2x5)',
                   texname = '\\text{I78a22}')

I78a25 = Parameter(name = 'I78a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl2x5)',
                   texname = '\\text{I78a25}')

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

I78a41 = Parameter(name = 'I78a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I78a41}')

I78a44 = Parameter(name = 'I78a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I78a44}')

I78a52 = Parameter(name = 'I78a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I78a52}')

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

I79a11 = Parameter(name = 'I79a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*complexconjugate(Rl1x4)',
                   texname = '\\text{I79a11}')

I79a14 = Parameter(name = 'I79a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl1x4)',
                   texname = '\\text{I79a14}')

I79a22 = Parameter(name = 'I79a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*complexconjugate(Rl2x5)',
                   texname = '\\text{I79a22}')

I79a25 = Parameter(name = 'I79a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl2x5)',
                   texname = '\\text{I79a25}')

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

I79a41 = Parameter(name = 'I79a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I79a41}')

I79a44 = Parameter(name = 'I79a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I79a44}')

I79a52 = Parameter(name = 'I79a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I79a52}')

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

I8a14 = Parameter(name = 'I8a14',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd4x1*complexconjugate(Rd1x1)',
                  texname = '\\text{I8a14}')

I8a22 = Parameter(name = 'I8a22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd2x2*complexconjugate(Rd2x2)',
                  texname = '\\text{I8a22}')

I8a25 = Parameter(name = 'I8a25',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd5x2*complexconjugate(Rd2x2)',
                  texname = '\\text{I8a25}')

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

I8a41 = Parameter(name = 'I8a41',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd1x1*complexconjugate(Rd4x1)',
                  texname = '\\text{I8a41}')

I8a44 = Parameter(name = 'I8a44',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd4x1*complexconjugate(Rd4x1)',
                  texname = '\\text{I8a44}')

I8a52 = Parameter(name = 'I8a52',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd2x2*complexconjugate(Rd5x2)',
                  texname = '\\text{I8a52}')

I8a55 = Parameter(name = 'I8a55',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd5x2*complexconjugate(Rd5x2)',
                  texname = '\\text{I8a55}')

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

I80a11 = Parameter(name = 'I80a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*complexconjugate(Rl1x4)',
                   texname = '\\text{I80a11}')

I80a14 = Parameter(name = 'I80a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl1x4)',
                   texname = '\\text{I80a14}')

I80a22 = Parameter(name = 'I80a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*complexconjugate(Rl2x5)',
                   texname = '\\text{I80a22}')

I80a25 = Parameter(name = 'I80a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl2x5)',
                   texname = '\\text{I80a25}')

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

I80a41 = Parameter(name = 'I80a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I80a41}')

I80a44 = Parameter(name = 'I80a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I80a44}')

I80a52 = Parameter(name = 'I80a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I80a52}')

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

I81a11 = Parameter(name = 'I81a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*complexconjugate(Rl1x4)',
                   texname = '\\text{I81a11}')

I81a14 = Parameter(name = 'I81a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl1x4)',
                   texname = '\\text{I81a14}')

I81a22 = Parameter(name = 'I81a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*complexconjugate(Rl2x5)',
                   texname = '\\text{I81a22}')

I81a25 = Parameter(name = 'I81a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl2x5)',
                   texname = '\\text{I81a25}')

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

I81a41 = Parameter(name = 'I81a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I81a41}')

I81a44 = Parameter(name = 'I81a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I81a44}')

I81a52 = Parameter(name = 'I81a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I81a52}')

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

I82a11 = Parameter(name = 'I82a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I82a11}')

I82a14 = Parameter(name = 'I82a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I82a14}')

I82a22 = Parameter(name = 'I82a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I82a22}')

I82a25 = Parameter(name = 'I82a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I82a25}')

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

I82a41 = Parameter(name = 'I82a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I82a41}')

I82a44 = Parameter(name = 'I82a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I82a44}')

I82a52 = Parameter(name = 'I82a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I82a52}')

I82a55 = Parameter(name = 'I82a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I82a55}')

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

I83a11 = Parameter(name = 'I83a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*ye1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I83a11}')

I83a14 = Parameter(name = 'I83a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*ye1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I83a14}')

I83a22 = Parameter(name = 'I83a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*ye2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I83a22}')

I83a25 = Parameter(name = 'I83a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*ye2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I83a25}')

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

I83a41 = Parameter(name = 'I83a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*ye1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I83a41}')

I83a44 = Parameter(name = 'I83a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*ye1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I83a44}')

I83a52 = Parameter(name = 'I83a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*ye2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I83a52}')

I83a55 = Parameter(name = 'I83a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*ye2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I83a55}')

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

I84a11 = Parameter(name = 'I84a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I84a11}')

I84a14 = Parameter(name = 'I84a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I84a14}')

I84a22 = Parameter(name = 'I84a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I84a22}')

I84a25 = Parameter(name = 'I84a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I84a25}')

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

I84a41 = Parameter(name = 'I84a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I84a41}')

I84a44 = Parameter(name = 'I84a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I84a44}')

I84a52 = Parameter(name = 'I84a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I84a52}')

I84a55 = Parameter(name = 'I84a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I84a55}')

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

I85a11 = Parameter(name = 'I85a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*ye1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I85a11}')

I85a14 = Parameter(name = 'I85a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*ye1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I85a14}')

I85a22 = Parameter(name = 'I85a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*ye2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I85a22}')

I85a25 = Parameter(name = 'I85a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*ye2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I85a25}')

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

I85a41 = Parameter(name = 'I85a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*ye1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I85a41}')

I85a44 = Parameter(name = 'I85a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*ye1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I85a44}')

I85a52 = Parameter(name = 'I85a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*ye2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I85a52}')

I85a55 = Parameter(name = 'I85a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*ye2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I85a55}')

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

I86a11 = Parameter(name = 'I86a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I86a11}')

I86a14 = Parameter(name = 'I86a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I86a14}')

I86a22 = Parameter(name = 'I86a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I86a22}')

I86a25 = Parameter(name = 'I86a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I86a25}')

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

I86a41 = Parameter(name = 'I86a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I86a41}')

I86a44 = Parameter(name = 'I86a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I86a44}')

I86a52 = Parameter(name = 'I86a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I86a52}')

I86a55 = Parameter(name = 'I86a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I86a55}')

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

I87a11 = Parameter(name = 'I87a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*ye1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I87a11}')

I87a14 = Parameter(name = 'I87a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*ye1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I87a14}')

I87a22 = Parameter(name = 'I87a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*ye2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I87a22}')

I87a25 = Parameter(name = 'I87a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*ye2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I87a25}')

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

I87a41 = Parameter(name = 'I87a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*ye1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I87a41}')

I87a44 = Parameter(name = 'I87a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*ye1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I87a44}')

I87a52 = Parameter(name = 'I87a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*ye2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I87a52}')

I87a55 = Parameter(name = 'I87a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*ye2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I87a55}')

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

I88a11 = Parameter(name = 'I88a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I88a11}')

I88a14 = Parameter(name = 'I88a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl1x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I88a14}')

I88a22 = Parameter(name = 'I88a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I88a22}')

I88a25 = Parameter(name = 'I88a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl2x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I88a25}')

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

I88a41 = Parameter(name = 'I88a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I88a41}')

I88a44 = Parameter(name = 'I88a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl4x4)*complexconjugate(ye1x1)',
                   texname = '\\text{I88a44}')

I88a52 = Parameter(name = 'I88a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I88a52}')

I88a55 = Parameter(name = 'I88a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl5x5)*complexconjugate(ye2x2)',
                   texname = '\\text{I88a55}')

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

I89a11 = Parameter(name = 'I89a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*ye1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I89a11}')

I89a14 = Parameter(name = 'I89a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*ye1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I89a14}')

I89a22 = Parameter(name = 'I89a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*ye2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I89a22}')

I89a25 = Parameter(name = 'I89a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*ye2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I89a25}')

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

I89a41 = Parameter(name = 'I89a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*ye1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I89a41}')

I89a44 = Parameter(name = 'I89a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*ye1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I89a44}')

I89a52 = Parameter(name = 'I89a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*ye2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I89a52}')

I89a55 = Parameter(name = 'I89a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*ye2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I89a55}')

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

I9a11 = Parameter(name = 'I9a11',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd1x4*complexconjugate(Rd1x4)',
                  texname = '\\text{I9a11}')

I9a14 = Parameter(name = 'I9a14',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd4x4*complexconjugate(Rd1x4)',
                  texname = '\\text{I9a14}')

I9a22 = Parameter(name = 'I9a22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd2x5*complexconjugate(Rd2x5)',
                  texname = '\\text{I9a22}')

I9a25 = Parameter(name = 'I9a25',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd5x5*complexconjugate(Rd2x5)',
                  texname = '\\text{I9a25}')

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

I9a41 = Parameter(name = 'I9a41',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd1x4*complexconjugate(Rd4x4)',
                  texname = '\\text{I9a41}')

I9a44 = Parameter(name = 'I9a44',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd4x4*complexconjugate(Rd4x4)',
                  texname = '\\text{I9a44}')

I9a52 = Parameter(name = 'I9a52',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd2x5*complexconjugate(Rd5x5)',
                  texname = '\\text{I9a52}')

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

I90a14 = Parameter(name = 'I90a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rn1x1)',
                   texname = '\\text{I90a14}')

I90a22 = Parameter(name = 'I90a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rn2x2)',
                   texname = '\\text{I90a22}')

I90a25 = Parameter(name = 'I90a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rn2x2)',
                   texname = '\\text{I90a25}')

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

I91a11 = Parameter(name = 'I91a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*te1x1*complexconjugate(Rn1x1)',
                   texname = '\\text{I91a11}')

I91a14 = Parameter(name = 'I91a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*te1x1*complexconjugate(Rn1x1)',
                   texname = '\\text{I91a14}')

I91a22 = Parameter(name = 'I91a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*te2x2*complexconjugate(Rn2x2)',
                   texname = '\\text{I91a22}')

I91a25 = Parameter(name = 'I91a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*te2x2*complexconjugate(Rn2x2)',
                   texname = '\\text{I91a25}')

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

I92a11 = Parameter(name = 'I92a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*ye1x1*complexconjugate(Rn1x1)*complexconjugate(ye1x1)',
                   texname = '\\text{I92a11}')

I92a14 = Parameter(name = 'I92a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*ye1x1*complexconjugate(Rn1x1)*complexconjugate(ye1x1)',
                   texname = '\\text{I92a14}')

I92a22 = Parameter(name = 'I92a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*ye2x2*complexconjugate(Rn2x2)*complexconjugate(ye2x2)',
                   texname = '\\text{I92a22}')

I92a25 = Parameter(name = 'I92a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*ye2x2*complexconjugate(Rn2x2)*complexconjugate(ye2x2)',
                   texname = '\\text{I92a25}')

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

I93a11 = Parameter(name = 'I93a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*ye1x1*complexconjugate(Rn1x1)',
                   texname = '\\text{I93a11}')

I93a14 = Parameter(name = 'I93a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*ye1x1*complexconjugate(Rn1x1)',
                   texname = '\\text{I93a14}')

I93a22 = Parameter(name = 'I93a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*ye2x2*complexconjugate(Rn2x2)',
                   texname = '\\text{I93a22}')

I93a25 = Parameter(name = 'I93a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*ye2x2*complexconjugate(Rn2x2)',
                   texname = '\\text{I93a25}')

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

I94a14 = Parameter(name = 'I94a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rn1x1)',
                   texname = '\\text{I94a14}')

I94a22 = Parameter(name = 'I94a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rn2x2)',
                   texname = '\\text{I94a22}')

I94a25 = Parameter(name = 'I94a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rn2x2)',
                   texname = '\\text{I94a25}')

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

I95a11 = Parameter(name = 'I95a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*ye1x1*complexconjugate(Rn1x1)*complexconjugate(ye1x1)',
                   texname = '\\text{I95a11}')

I95a14 = Parameter(name = 'I95a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*ye1x1*complexconjugate(Rn1x1)*complexconjugate(ye1x1)',
                   texname = '\\text{I95a14}')

I95a22 = Parameter(name = 'I95a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*ye2x2*complexconjugate(Rn2x2)*complexconjugate(ye2x2)',
                   texname = '\\text{I95a22}')

I95a25 = Parameter(name = 'I95a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*ye2x2*complexconjugate(Rn2x2)*complexconjugate(ye2x2)',
                   texname = '\\text{I95a25}')

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

I96a11 = Parameter(name = 'I96a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn1x1*ye1x1*complexconjugate(Rn1x1)*complexconjugate(ye1x1)',
                   texname = '\\text{I96a11}')

I96a22 = Parameter(name = 'I96a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn2x2*ye2x2*complexconjugate(Rn2x2)*complexconjugate(ye2x2)',
                   texname = '\\text{I96a22}')

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

I97a14 = Parameter(name = 'I97a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I97a14}')

I97a22 = Parameter(name = 'I97a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I97a22}')

I97a25 = Parameter(name = 'I97a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I97a25}')

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

I97a41 = Parameter(name = 'I97a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I97a41}')

I97a44 = Parameter(name = 'I97a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rl4x1)',
                   texname = '\\text{I97a44}')

I97a52 = Parameter(name = 'I97a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I97a52}')

I97a55 = Parameter(name = 'I97a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rl5x2)',
                   texname = '\\text{I97a55}')

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

I98a11 = Parameter(name = 'I98a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*complexconjugate(Rl1x4)',
                   texname = '\\text{I98a11}')

I98a14 = Parameter(name = 'I98a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl1x4)',
                   texname = '\\text{I98a14}')

I98a22 = Parameter(name = 'I98a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*complexconjugate(Rl2x5)',
                   texname = '\\text{I98a22}')

I98a25 = Parameter(name = 'I98a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl2x5)',
                   texname = '\\text{I98a25}')

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

I98a41 = Parameter(name = 'I98a41',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I98a41}')

I98a44 = Parameter(name = 'I98a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I98a44}')

I98a52 = Parameter(name = 'I98a52',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I98a52}')

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

I99a14 = Parameter(name = 'I99a14',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x1*complexconjugate(Rn1x1)',
                   texname = '\\text{I99a14}')

I99a22 = Parameter(name = 'I99a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rn2x2)',
                   texname = '\\text{I99a22}')

I99a25 = Parameter(name = 'I99a25',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x2*complexconjugate(Rn2x2)',
                   texname = '\\text{I99a25}')

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

