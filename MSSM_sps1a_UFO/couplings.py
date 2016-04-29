# This file was automatically created by FeynRules 2.3.21
# Mathematica version: 10.3.0 for Mac OS X x86 (64-bit) (October 9, 2015)
# Date: Sat 30 Apr 2016 19:04:51


from object_library import all_couplings, Coupling

from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot



GC_1 = Coupling(name = 'GC_1',
                value = '-(ee*complex(0,1))/3.',
                order = {'QED':1})

GC_2 = Coupling(name = 'GC_2',
                value = '(2*ee*complex(0,1))/3.',
                order = {'QED':1})

GC_3 = Coupling(name = 'GC_3',
                value = '-(ee*complex(0,1))',
                order = {'QED':1})

GC_4 = Coupling(name = 'GC_4',
                value = 'ee*complex(0,1)',
                order = {'QED':1})

GC_5 = Coupling(name = 'GC_5',
                value = 'ee**2*complex(0,1)',
                order = {'QED':2})

GC_6 = Coupling(name = 'GC_6',
                value = '-G',
                order = {'QCD':1})

GC_7 = Coupling(name = 'GC_7',
                value = 'complex(0,1)*G',
                order = {'QCD':1})

GC_8 = Coupling(name = 'GC_8',
                value = 'G/2.',
                order = {'QCD':1})

GC_9 = Coupling(name = 'GC_9',
                value = 'complex(0,1)*G**2',
                order = {'QCD':2})

GC_10 = Coupling(name = 'GC_10',
                 value = '(-2*ee*complex(0,1)*I113a11)/3.',
                 order = {'QED':1})

GC_11 = Coupling(name = 'GC_11',
                 value = '-(complex(0,1)*G*I113a11)',
                 order = {'QCD':1})

GC_12 = Coupling(name = 'GC_12',
                 value = '(-2*ee*complex(0,1)*I113a22)/3.',
                 order = {'QED':1})

GC_13 = Coupling(name = 'GC_13',
                 value = '-(complex(0,1)*G*I113a22)',
                 order = {'QCD':1})

GC_14 = Coupling(name = 'GC_14',
                 value = '(-2*ee*complex(0,1)*I113a33)/3. - (2*ee*complex(0,1)*I114a33)/3.',
                 order = {'QED':1})

GC_15 = Coupling(name = 'GC_15',
                 value = '-(complex(0,1)*G*I113a33) - complex(0,1)*G*I114a33',
                 order = {'QCD':1})

GC_16 = Coupling(name = 'GC_16',
                 value = '(-2*ee*complex(0,1)*I113a36)/3. - (2*ee*complex(0,1)*I114a36)/3.',
                 order = {'QED':1})

GC_17 = Coupling(name = 'GC_17',
                 value = '-(complex(0,1)*G*I113a36) - complex(0,1)*G*I114a36',
                 order = {'QCD':1})

GC_18 = Coupling(name = 'GC_18',
                 value = '(-2*ee*complex(0,1)*I114a44)/3.',
                 order = {'QED':1})

GC_19 = Coupling(name = 'GC_19',
                 value = '-(complex(0,1)*G*I114a44)',
                 order = {'QCD':1})

GC_20 = Coupling(name = 'GC_20',
                 value = '(-2*ee*complex(0,1)*I114a55)/3.',
                 order = {'QED':1})

GC_21 = Coupling(name = 'GC_21',
                 value = '-(complex(0,1)*G*I114a55)',
                 order = {'QCD':1})

GC_22 = Coupling(name = 'GC_22',
                 value = '(2*ee*complex(0,1)*I113a63)/3. + (2*ee*complex(0,1)*I114a63)/3.',
                 order = {'QED':1})

GC_23 = Coupling(name = 'GC_23',
                 value = 'complex(0,1)*G*I113a63 + complex(0,1)*G*I114a63',
                 order = {'QCD':1})

GC_24 = Coupling(name = 'GC_24',
                 value = '(-2*ee*complex(0,1)*I113a66)/3. - (2*ee*complex(0,1)*I114a66)/3.',
                 order = {'QED':1})

GC_25 = Coupling(name = 'GC_25',
                 value = '-(complex(0,1)*G*I113a66) - complex(0,1)*G*I114a66',
                 order = {'QCD':1})

GC_26 = Coupling(name = 'GC_26',
                 value = '(8*ee**2*complex(0,1)*I148a11)/9.',
                 order = {'QED':2})

GC_27 = Coupling(name = 'GC_27',
                 value = '(4*ee*complex(0,1)*G*I148a11)/3.',
                 order = {'QCD':1,'QED':1})

GC_28 = Coupling(name = 'GC_28',
                 value = 'complex(0,1)*G**2*I148a11',
                 order = {'QCD':2})

GC_29 = Coupling(name = 'GC_29',
                 value = '(8*ee**2*complex(0,1)*I148a22)/9.',
                 order = {'QED':2})

GC_30 = Coupling(name = 'GC_30',
                 value = '(4*ee*complex(0,1)*G*I148a22)/3.',
                 order = {'QCD':1,'QED':1})

GC_31 = Coupling(name = 'GC_31',
                 value = 'complex(0,1)*G**2*I148a22',
                 order = {'QCD':2})

GC_32 = Coupling(name = 'GC_32',
                 value = '(8*ee**2*complex(0,1)*I148a33)/9. + (8*ee**2*complex(0,1)*I149a33)/9.',
                 order = {'QED':2})

GC_33 = Coupling(name = 'GC_33',
                 value = '(4*ee*complex(0,1)*G*I148a33)/3. + (4*ee*complex(0,1)*G*I149a33)/3.',
                 order = {'QCD':1,'QED':1})

GC_34 = Coupling(name = 'GC_34',
                 value = 'complex(0,1)*G**2*I148a33 + complex(0,1)*G**2*I149a33',
                 order = {'QCD':2})

GC_35 = Coupling(name = 'GC_35',
                 value = '(8*ee**2*complex(0,1)*I148a36)/9. + (8*ee**2*complex(0,1)*I149a36)/9.',
                 order = {'QED':2})

GC_36 = Coupling(name = 'GC_36',
                 value = '(4*ee*complex(0,1)*G*I148a36)/3. + (4*ee*complex(0,1)*G*I149a36)/3.',
                 order = {'QCD':1,'QED':1})

GC_37 = Coupling(name = 'GC_37',
                 value = 'complex(0,1)*G**2*I148a36 + complex(0,1)*G**2*I149a36',
                 order = {'QCD':2})

GC_38 = Coupling(name = 'GC_38',
                 value = '(8*ee**2*complex(0,1)*I149a44)/9.',
                 order = {'QED':2})

GC_39 = Coupling(name = 'GC_39',
                 value = '(4*ee*complex(0,1)*G*I149a44)/3.',
                 order = {'QCD':1,'QED':1})

GC_40 = Coupling(name = 'GC_40',
                 value = 'complex(0,1)*G**2*I149a44',
                 order = {'QCD':2})

GC_41 = Coupling(name = 'GC_41',
                 value = '(8*ee**2*complex(0,1)*I149a55)/9.',
                 order = {'QED':2})

GC_42 = Coupling(name = 'GC_42',
                 value = '(4*ee*complex(0,1)*G*I149a55)/3.',
                 order = {'QCD':1,'QED':1})

GC_43 = Coupling(name = 'GC_43',
                 value = 'complex(0,1)*G**2*I149a55',
                 order = {'QCD':2})

GC_44 = Coupling(name = 'GC_44',
                 value = '(8*ee**2*complex(0,1)*I148a63)/9. + (8*ee**2*complex(0,1)*I149a63)/9.',
                 order = {'QED':2})

GC_45 = Coupling(name = 'GC_45',
                 value = '(4*ee*complex(0,1)*G*I148a63)/3. + (4*ee*complex(0,1)*G*I149a63)/3.',
                 order = {'QCD':1,'QED':1})

GC_46 = Coupling(name = 'GC_46',
                 value = 'complex(0,1)*G**2*I148a63 + complex(0,1)*G**2*I149a63',
                 order = {'QCD':2})

GC_47 = Coupling(name = 'GC_47',
                 value = '(8*ee**2*complex(0,1)*I148a66)/9. + (8*ee**2*complex(0,1)*I149a66)/9.',
                 order = {'QED':2})

GC_48 = Coupling(name = 'GC_48',
                 value = '(4*ee*complex(0,1)*G*I148a66)/3. + (4*ee*complex(0,1)*G*I149a66)/3.',
                 order = {'QCD':1,'QED':1})

GC_49 = Coupling(name = 'GC_49',
                 value = 'complex(0,1)*G**2*I148a66 + complex(0,1)*G**2*I149a66',
                 order = {'QCD':2})

GC_50 = Coupling(name = 'GC_50',
                 value = '(2*ee**2*complex(0,1)*I15a11)/9.',
                 order = {'QED':2})

GC_51 = Coupling(name = 'GC_51',
                 value = '(-2*ee*complex(0,1)*G*I15a11)/3.',
                 order = {'QCD':1,'QED':1})

GC_52 = Coupling(name = 'GC_52',
                 value = 'complex(0,1)*G**2*I15a11',
                 order = {'QCD':2})

GC_53 = Coupling(name = 'GC_53',
                 value = '(2*ee**2*complex(0,1)*I15a22)/9.',
                 order = {'QED':2})

GC_54 = Coupling(name = 'GC_54',
                 value = '(-2*ee*complex(0,1)*G*I15a22)/3.',
                 order = {'QCD':1,'QED':1})

GC_55 = Coupling(name = 'GC_55',
                 value = 'complex(0,1)*G**2*I15a22',
                 order = {'QCD':2})

GC_56 = Coupling(name = 'GC_56',
                 value = '(2*ee**2*complex(0,1)*I15a33)/9. + (2*ee**2*complex(0,1)*I16a33)/9.',
                 order = {'QED':2})

GC_57 = Coupling(name = 'GC_57',
                 value = '(-2*ee*complex(0,1)*G*I15a33)/3. - (2*ee*complex(0,1)*G*I16a33)/3.',
                 order = {'QCD':1,'QED':1})

GC_58 = Coupling(name = 'GC_58',
                 value = 'complex(0,1)*G**2*I15a33 + complex(0,1)*G**2*I16a33',
                 order = {'QCD':2})

GC_59 = Coupling(name = 'GC_59',
                 value = '(2*ee**2*complex(0,1)*I15a36)/9. + (2*ee**2*complex(0,1)*I16a36)/9.',
                 order = {'QED':2})

GC_60 = Coupling(name = 'GC_60',
                 value = '(-2*ee*complex(0,1)*G*I15a36)/3. - (2*ee*complex(0,1)*G*I16a36)/3.',
                 order = {'QCD':1,'QED':1})

GC_61 = Coupling(name = 'GC_61',
                 value = 'complex(0,1)*G**2*I15a36 + complex(0,1)*G**2*I16a36',
                 order = {'QCD':2})

GC_62 = Coupling(name = 'GC_62',
                 value = '(2*ee**2*complex(0,1)*I16a44)/9.',
                 order = {'QED':2})

GC_63 = Coupling(name = 'GC_63',
                 value = '(-2*ee*complex(0,1)*G*I16a44)/3.',
                 order = {'QCD':1,'QED':1})

GC_64 = Coupling(name = 'GC_64',
                 value = 'complex(0,1)*G**2*I16a44',
                 order = {'QCD':2})

GC_65 = Coupling(name = 'GC_65',
                 value = '(2*ee**2*complex(0,1)*I16a55)/9.',
                 order = {'QED':2})

GC_66 = Coupling(name = 'GC_66',
                 value = '(-2*ee*complex(0,1)*G*I16a55)/3.',
                 order = {'QCD':1,'QED':1})

GC_67 = Coupling(name = 'GC_67',
                 value = 'complex(0,1)*G**2*I16a55',
                 order = {'QCD':2})

GC_68 = Coupling(name = 'GC_68',
                 value = '(2*ee**2*complex(0,1)*I15a63)/9. + (2*ee**2*complex(0,1)*I16a63)/9.',
                 order = {'QED':2})

GC_69 = Coupling(name = 'GC_69',
                 value = '(-2*ee*complex(0,1)*G*I15a63)/3. - (2*ee*complex(0,1)*G*I16a63)/3.',
                 order = {'QCD':1,'QED':1})

GC_70 = Coupling(name = 'GC_70',
                 value = 'complex(0,1)*G**2*I15a63 + complex(0,1)*G**2*I16a63',
                 order = {'QCD':2})

GC_71 = Coupling(name = 'GC_71',
                 value = '(2*ee**2*complex(0,1)*I15a66)/9. + (2*ee**2*complex(0,1)*I16a66)/9.',
                 order = {'QED':2})

GC_72 = Coupling(name = 'GC_72',
                 value = '(-2*ee*complex(0,1)*G*I15a66)/3. - (2*ee*complex(0,1)*G*I16a66)/3.',
                 order = {'QCD':1,'QED':1})

GC_73 = Coupling(name = 'GC_73',
                 value = 'complex(0,1)*G**2*I15a66 + complex(0,1)*G**2*I16a66',
                 order = {'QCD':2})

GC_74 = Coupling(name = 'GC_74',
                 value = 'ee*complex(0,1)*I47a11',
                 order = {'QED':1})

GC_75 = Coupling(name = 'GC_75',
                 value = 'ee*complex(0,1)*I47a22',
                 order = {'QED':1})

GC_76 = Coupling(name = 'GC_76',
                 value = 'ee*complex(0,1)*I47a33 + ee*complex(0,1)*I48a33',
                 order = {'QED':1})

GC_77 = Coupling(name = 'GC_77',
                 value = 'ee*complex(0,1)*I47a36 + ee*complex(0,1)*I48a36',
                 order = {'QED':1})

GC_78 = Coupling(name = 'GC_78',
                 value = 'ee*complex(0,1)*I48a44',
                 order = {'QED':1})

GC_79 = Coupling(name = 'GC_79',
                 value = 'ee*complex(0,1)*I48a55',
                 order = {'QED':1})

GC_80 = Coupling(name = 'GC_80',
                 value = '-(ee*complex(0,1)*I47a63) - ee*complex(0,1)*I48a63',
                 order = {'QED':1})

GC_81 = Coupling(name = 'GC_81',
                 value = 'ee*complex(0,1)*I47a66 + ee*complex(0,1)*I48a66',
                 order = {'QED':1})

GC_82 = Coupling(name = 'GC_82',
                 value = '2*ee**2*complex(0,1)*I53a11',
                 order = {'QED':2})

GC_83 = Coupling(name = 'GC_83',
                 value = '2*ee**2*complex(0,1)*I53a22',
                 order = {'QED':2})

GC_84 = Coupling(name = 'GC_84',
                 value = '2*ee**2*complex(0,1)*I53a33 + 2*ee**2*complex(0,1)*I54a33',
                 order = {'QED':2})

GC_85 = Coupling(name = 'GC_85',
                 value = '2*ee**2*complex(0,1)*I53a36 + 2*ee**2*complex(0,1)*I54a36',
                 order = {'QED':2})

GC_86 = Coupling(name = 'GC_86',
                 value = '2*ee**2*complex(0,1)*I54a44',
                 order = {'QED':2})

GC_87 = Coupling(name = 'GC_87',
                 value = '2*ee**2*complex(0,1)*I54a55',
                 order = {'QED':2})

GC_88 = Coupling(name = 'GC_88',
                 value = '2*ee**2*complex(0,1)*I53a63 + 2*ee**2*complex(0,1)*I54a63',
                 order = {'QED':2})

GC_89 = Coupling(name = 'GC_89',
                 value = '2*ee**2*complex(0,1)*I53a66 + 2*ee**2*complex(0,1)*I54a66',
                 order = {'QED':2})

GC_90 = Coupling(name = 'GC_90',
                 value = '(ee*complex(0,1)*I8a11)/3.',
                 order = {'QED':1})

GC_91 = Coupling(name = 'GC_91',
                 value = '-(complex(0,1)*G*I8a11)',
                 order = {'QCD':1})

GC_92 = Coupling(name = 'GC_92',
                 value = '(ee*complex(0,1)*I8a22)/3.',
                 order = {'QED':1})

GC_93 = Coupling(name = 'GC_93',
                 value = '-(complex(0,1)*G*I8a22)',
                 order = {'QCD':1})

GC_94 = Coupling(name = 'GC_94',
                 value = '(ee*complex(0,1)*I8a33)/3. + (ee*complex(0,1)*I9a33)/3.',
                 order = {'QED':1})

GC_95 = Coupling(name = 'GC_95',
                 value = '-(complex(0,1)*G*I8a33) - complex(0,1)*G*I9a33',
                 order = {'QCD':1})

GC_96 = Coupling(name = 'GC_96',
                 value = '(ee*complex(0,1)*I8a36)/3. + (ee*complex(0,1)*I9a36)/3.',
                 order = {'QED':1})

GC_97 = Coupling(name = 'GC_97',
                 value = '-(complex(0,1)*G*I8a36) - complex(0,1)*G*I9a36',
                 order = {'QCD':1})

GC_98 = Coupling(name = 'GC_98',
                 value = '(ee*complex(0,1)*I9a44)/3.',
                 order = {'QED':1})

GC_99 = Coupling(name = 'GC_99',
                 value = '-(complex(0,1)*G*I9a44)',
                 order = {'QCD':1})

GC_100 = Coupling(name = 'GC_100',
                  value = '(ee*complex(0,1)*I9a55)/3.',
                  order = {'QED':1})

GC_101 = Coupling(name = 'GC_101',
                  value = '-(complex(0,1)*G*I9a55)',
                  order = {'QCD':1})

GC_102 = Coupling(name = 'GC_102',
                  value = '-(ee*complex(0,1)*I8a63)/3. - (ee*complex(0,1)*I9a63)/3.',
                  order = {'QED':1})

GC_103 = Coupling(name = 'GC_103',
                  value = 'complex(0,1)*G*I8a63 + complex(0,1)*G*I9a63',
                  order = {'QCD':1})

GC_104 = Coupling(name = 'GC_104',
                  value = '(ee*complex(0,1)*I8a66)/3. + (ee*complex(0,1)*I9a66)/3.',
                  order = {'QED':1})

GC_105 = Coupling(name = 'GC_105',
                  value = '-(complex(0,1)*G*I8a66) - complex(0,1)*G*I9a66',
                  order = {'QCD':1})

GC_106 = Coupling(name = 'GC_106',
                  value = '-(complex(0,1)*G*Rd1x1*cmath.sqrt(2))',
                  order = {'QCD':1})

GC_107 = Coupling(name = 'GC_107',
                  value = '-(complex(0,1)*G*Rd2x2*cmath.sqrt(2))',
                  order = {'QCD':1})

GC_108 = Coupling(name = 'GC_108',
                  value = '-(complex(0,1)*G*Rd3x3*cmath.sqrt(2))',
                  order = {'QCD':1})

GC_109 = Coupling(name = 'GC_109',
                  value = 'complex(0,1)*G*Rd3x6*cmath.sqrt(2)',
                  order = {'QCD':1})

GC_110 = Coupling(name = 'GC_110',
                  value = 'complex(0,1)*G*Rd4x4*cmath.sqrt(2)',
                  order = {'QCD':1})

GC_111 = Coupling(name = 'GC_111',
                  value = 'complex(0,1)*G*Rd5x5*cmath.sqrt(2)',
                  order = {'QCD':1})

GC_112 = Coupling(name = 'GC_112',
                  value = '-(complex(0,1)*G*Rd6x3*cmath.sqrt(2))',
                  order = {'QCD':1})

GC_113 = Coupling(name = 'GC_113',
                  value = 'complex(0,1)*G*Rd6x6*cmath.sqrt(2)',
                  order = {'QCD':1})

GC_114 = Coupling(name = 'GC_114',
                  value = '-(complex(0,1)*G*Ru1x1*cmath.sqrt(2))',
                  order = {'QCD':1})

GC_115 = Coupling(name = 'GC_115',
                  value = '-(complex(0,1)*G*Ru2x2*cmath.sqrt(2))',
                  order = {'QCD':1})

GC_116 = Coupling(name = 'GC_116',
                  value = '-(complex(0,1)*G*Ru3x3*cmath.sqrt(2))',
                  order = {'QCD':1})

GC_117 = Coupling(name = 'GC_117',
                  value = 'complex(0,1)*G*Ru3x6*cmath.sqrt(2)',
                  order = {'QCD':1})

GC_118 = Coupling(name = 'GC_118',
                  value = 'complex(0,1)*G*Ru4x4*cmath.sqrt(2)',
                  order = {'QCD':1})

GC_119 = Coupling(name = 'GC_119',
                  value = 'complex(0,1)*G*Ru5x5*cmath.sqrt(2)',
                  order = {'QCD':1})

GC_120 = Coupling(name = 'GC_120',
                  value = '-(complex(0,1)*G*Ru6x3*cmath.sqrt(2))',
                  order = {'QCD':1})

GC_121 = Coupling(name = 'GC_121',
                  value = 'complex(0,1)*G*Ru6x6*cmath.sqrt(2)',
                  order = {'QCD':1})

GC_122 = Coupling(name = 'GC_122',
                  value = '-(ee**2*complex(0,1)) + (ee**2*complex(0,1))/sw**2',
                  order = {'QED':2})

GC_123 = Coupling(name = 'GC_123',
                  value = '-(complex(0,1)*I129a33*I130a33) - (ee**2*complex(0,1)*I127a33*I128a33)/(2.*sw**2)',
                  order = {'QED':2})

GC_124 = Coupling(name = 'GC_124',
                  value = '-(complex(0,1)*I129a33*I130a36) - (ee**2*complex(0,1)*I127a36*I128a33)/(2.*sw**2)',
                  order = {'QED':2})

GC_125 = Coupling(name = 'GC_125',
                  value = '-(complex(0,1)*I129a33*I130a63) - (ee**2*complex(0,1)*I127a63*I128a33)/(2.*sw**2)',
                  order = {'QED':2})

GC_126 = Coupling(name = 'GC_126',
                  value = '-(complex(0,1)*I129a33*I130a66) - (ee**2*complex(0,1)*I127a66*I128a33)/(2.*sw**2)',
                  order = {'QED':2})

GC_127 = Coupling(name = 'GC_127',
                  value = '-(complex(0,1)*I129a36*I130a33) - (ee**2*complex(0,1)*I127a33*I128a36)/(2.*sw**2)',
                  order = {'QED':2})

GC_128 = Coupling(name = 'GC_128',
                  value = '-(complex(0,1)*I129a36*I130a36) - (ee**2*complex(0,1)*I127a36*I128a36)/(2.*sw**2)',
                  order = {'QED':2})

GC_129 = Coupling(name = 'GC_129',
                  value = '-(complex(0,1)*I129a36*I130a63) - (ee**2*complex(0,1)*I127a63*I128a36)/(2.*sw**2)',
                  order = {'QED':2})

GC_130 = Coupling(name = 'GC_130',
                  value = '-(complex(0,1)*I129a36*I130a66) - (ee**2*complex(0,1)*I127a66*I128a36)/(2.*sw**2)',
                  order = {'QED':2})

GC_131 = Coupling(name = 'GC_131',
                  value = '-(complex(0,1)*I102a33*I160a33) - (ee**2*complex(0,1)*I159a33*I99a33)/(2.*sw**2)',
                  order = {'QED':2})

GC_132 = Coupling(name = 'GC_132',
                  value = '-(complex(0,1)*I102a33*I160a36) - (ee**2*complex(0,1)*I159a36*I99a33)/(2.*sw**2)',
                  order = {'QED':2})

GC_133 = Coupling(name = 'GC_133',
                  value = '-(complex(0,1)*I102a33*I160a63) - (ee**2*complex(0,1)*I159a63*I99a33)/(2.*sw**2)',
                  order = {'QED':2})

GC_134 = Coupling(name = 'GC_134',
                  value = '-(complex(0,1)*I102a33*I160a66) - (ee**2*complex(0,1)*I159a66*I99a33)/(2.*sw**2)',
                  order = {'QED':2})

GC_135 = Coupling(name = 'GC_135',
                  value = '-(complex(0,1)*I102a36*I160a33) - (ee**2*complex(0,1)*I159a33*I99a36)/(2.*sw**2)',
                  order = {'QED':2})

GC_136 = Coupling(name = 'GC_136',
                  value = '-(complex(0,1)*I102a36*I160a36) - (ee**2*complex(0,1)*I159a36*I99a36)/(2.*sw**2)',
                  order = {'QED':2})

GC_137 = Coupling(name = 'GC_137',
                  value = '-(complex(0,1)*I102a36*I160a63) - (ee**2*complex(0,1)*I159a63*I99a36)/(2.*sw**2)',
                  order = {'QED':2})

GC_138 = Coupling(name = 'GC_138',
                  value = '-(complex(0,1)*I102a36*I160a66) - (ee**2*complex(0,1)*I159a66*I99a36)/(2.*sw**2)',
                  order = {'QED':2})

GC_139 = Coupling(name = 'GC_139',
                  value = '(ee**2*complex(0,1))/(2.*sw**2)',
                  order = {'QED':2})

GC_140 = Coupling(name = 'GC_140',
                  value = '-((ee**2*complex(0,1))/sw**2)',
                  order = {'QED':2})

GC_141 = Coupling(name = 'GC_141',
                  value = '-(ee**2*complex(0,1)*I127a11*I128a11)/(2.*sw**2)',
                  order = {'QED':2})

GC_142 = Coupling(name = 'GC_142',
                  value = '-(ee**2*complex(0,1)*I127a22*I128a11)/(2.*sw**2)',
                  order = {'QED':2})

GC_143 = Coupling(name = 'GC_143',
                  value = '-(ee**2*complex(0,1)*I127a33*I128a11)/(2.*sw**2)',
                  order = {'QED':2})

GC_144 = Coupling(name = 'GC_144',
                  value = '-(ee**2*complex(0,1)*I127a36*I128a11)/(2.*sw**2)',
                  order = {'QED':2})

GC_145 = Coupling(name = 'GC_145',
                  value = '-(ee**2*complex(0,1)*I127a63*I128a11)/(2.*sw**2)',
                  order = {'QED':2})

GC_146 = Coupling(name = 'GC_146',
                  value = '-(ee**2*complex(0,1)*I127a66*I128a11)/(2.*sw**2)',
                  order = {'QED':2})

GC_147 = Coupling(name = 'GC_147',
                  value = '-(ee**2*complex(0,1)*I127a11*I128a22)/(2.*sw**2)',
                  order = {'QED':2})

GC_148 = Coupling(name = 'GC_148',
                  value = '-(ee**2*complex(0,1)*I127a22*I128a22)/(2.*sw**2)',
                  order = {'QED':2})

GC_149 = Coupling(name = 'GC_149',
                  value = '-(ee**2*complex(0,1)*I127a33*I128a22)/(2.*sw**2)',
                  order = {'QED':2})

GC_150 = Coupling(name = 'GC_150',
                  value = '-(ee**2*complex(0,1)*I127a36*I128a22)/(2.*sw**2)',
                  order = {'QED':2})

GC_151 = Coupling(name = 'GC_151',
                  value = '-(ee**2*complex(0,1)*I127a63*I128a22)/(2.*sw**2)',
                  order = {'QED':2})

GC_152 = Coupling(name = 'GC_152',
                  value = '-(ee**2*complex(0,1)*I127a66*I128a22)/(2.*sw**2)',
                  order = {'QED':2})

GC_153 = Coupling(name = 'GC_153',
                  value = '-(ee**2*complex(0,1)*I127a11*I128a33)/(2.*sw**2)',
                  order = {'QED':2})

GC_154 = Coupling(name = 'GC_154',
                  value = '-(ee**2*complex(0,1)*I127a22*I128a33)/(2.*sw**2)',
                  order = {'QED':2})

GC_155 = Coupling(name = 'GC_155',
                  value = '-(ee**2*complex(0,1)*I127a11*I128a36)/(2.*sw**2)',
                  order = {'QED':2})

GC_156 = Coupling(name = 'GC_156',
                  value = '-(ee**2*complex(0,1)*I127a22*I128a36)/(2.*sw**2)',
                  order = {'QED':2})

GC_157 = Coupling(name = 'GC_157',
                  value = '(ee**2*complex(0,1)*I201a11)/(2.*sw**2)',
                  order = {'QED':2})

GC_158 = Coupling(name = 'GC_158',
                  value = '(ee**2*complex(0,1)*I201a22)/(2.*sw**2)',
                  order = {'QED':2})

GC_159 = Coupling(name = 'GC_159',
                  value = '(ee**2*complex(0,1)*I201a33)/(2.*sw**2)',
                  order = {'QED':2})

GC_160 = Coupling(name = 'GC_160',
                  value = '(ee**2*complex(0,1)*I201a36)/(2.*sw**2)',
                  order = {'QED':2})

GC_161 = Coupling(name = 'GC_161',
                  value = '(ee**2*complex(0,1)*I201a63)/(2.*sw**2)',
                  order = {'QED':2})

GC_162 = Coupling(name = 'GC_162',
                  value = '(ee**2*complex(0,1)*I201a66)/(2.*sw**2)',
                  order = {'QED':2})

GC_163 = Coupling(name = 'GC_163',
                  value = '(ee**2*complex(0,1)*I63a11)/(2.*sw**2)',
                  order = {'QED':2})

GC_164 = Coupling(name = 'GC_164',
                  value = '(ee**2*complex(0,1)*I63a22)/(2.*sw**2)',
                  order = {'QED':2})

GC_165 = Coupling(name = 'GC_165',
                  value = '(ee**2*complex(0,1)*I63a33)/(2.*sw**2)',
                  order = {'QED':2})

GC_166 = Coupling(name = 'GC_166',
                  value = '(ee**2*complex(0,1)*I63a36)/(2.*sw**2)',
                  order = {'QED':2})

GC_167 = Coupling(name = 'GC_167',
                  value = '(ee**2*complex(0,1)*I63a63)/(2.*sw**2)',
                  order = {'QED':2})

GC_168 = Coupling(name = 'GC_168',
                  value = '(ee**2*complex(0,1)*I63a66)/(2.*sw**2)',
                  order = {'QED':2})

GC_169 = Coupling(name = 'GC_169',
                  value = '(ee**2*complex(0,1)*I97a11)/(2.*sw**2)',
                  order = {'QED':2})

GC_170 = Coupling(name = 'GC_170',
                  value = '(ee**2*complex(0,1)*I97a22)/(2.*sw**2)',
                  order = {'QED':2})

GC_171 = Coupling(name = 'GC_171',
                  value = '(ee**2*complex(0,1)*I97a33)/(2.*sw**2)',
                  order = {'QED':2})

GC_172 = Coupling(name = 'GC_172',
                  value = '(ee**2*complex(0,1)*I97a36)/(2.*sw**2)',
                  order = {'QED':2})

GC_173 = Coupling(name = 'GC_173',
                  value = '(ee**2*complex(0,1)*I97a63)/(2.*sw**2)',
                  order = {'QED':2})

GC_174 = Coupling(name = 'GC_174',
                  value = '(ee**2*complex(0,1)*I97a66)/(2.*sw**2)',
                  order = {'QED':2})

GC_175 = Coupling(name = 'GC_175',
                  value = '-(ee**2*complex(0,1)*I159a11*I99a11)/(2.*sw**2)',
                  order = {'QED':2})

GC_176 = Coupling(name = 'GC_176',
                  value = '-(ee**2*complex(0,1)*I159a22*I99a11)/(2.*sw**2)',
                  order = {'QED':2})

GC_177 = Coupling(name = 'GC_177',
                  value = '-(ee**2*complex(0,1)*I159a33*I99a11)/(2.*sw**2)',
                  order = {'QED':2})

GC_178 = Coupling(name = 'GC_178',
                  value = '-(ee**2*complex(0,1)*I159a36*I99a11)/(2.*sw**2)',
                  order = {'QED':2})

GC_179 = Coupling(name = 'GC_179',
                  value = '-(ee**2*complex(0,1)*I159a63*I99a11)/(2.*sw**2)',
                  order = {'QED':2})

GC_180 = Coupling(name = 'GC_180',
                  value = '-(ee**2*complex(0,1)*I159a66*I99a11)/(2.*sw**2)',
                  order = {'QED':2})

GC_181 = Coupling(name = 'GC_181',
                  value = '-(ee**2*complex(0,1)*I159a11*I99a22)/(2.*sw**2)',
                  order = {'QED':2})

GC_182 = Coupling(name = 'GC_182',
                  value = '-(ee**2*complex(0,1)*I159a22*I99a22)/(2.*sw**2)',
                  order = {'QED':2})

GC_183 = Coupling(name = 'GC_183',
                  value = '-(ee**2*complex(0,1)*I159a33*I99a22)/(2.*sw**2)',
                  order = {'QED':2})

GC_184 = Coupling(name = 'GC_184',
                  value = '-(ee**2*complex(0,1)*I159a36*I99a22)/(2.*sw**2)',
                  order = {'QED':2})

GC_185 = Coupling(name = 'GC_185',
                  value = '-(ee**2*complex(0,1)*I159a63*I99a22)/(2.*sw**2)',
                  order = {'QED':2})

GC_186 = Coupling(name = 'GC_186',
                  value = '-(ee**2*complex(0,1)*I159a66*I99a22)/(2.*sw**2)',
                  order = {'QED':2})

GC_187 = Coupling(name = 'GC_187',
                  value = '-(ee**2*complex(0,1)*I159a11*I99a33)/(2.*sw**2)',
                  order = {'QED':2})

GC_188 = Coupling(name = 'GC_188',
                  value = '-(ee**2*complex(0,1)*I159a22*I99a33)/(2.*sw**2)',
                  order = {'QED':2})

GC_189 = Coupling(name = 'GC_189',
                  value = '-(ee**2*complex(0,1)*I159a11*I99a36)/(2.*sw**2)',
                  order = {'QED':2})

GC_190 = Coupling(name = 'GC_190',
                  value = '-(ee**2*complex(0,1)*I159a22*I99a36)/(2.*sw**2)',
                  order = {'QED':2})

GC_191 = Coupling(name = 'GC_191',
                  value = '(ee*complex(0,1))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_192 = Coupling(name = 'GC_192',
                  value = '(CKM1x1*ee*complex(0,1))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_193 = Coupling(name = 'GC_193',
                  value = '(CKM2x2*ee*complex(0,1))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_194 = Coupling(name = 'GC_194',
                  value = '(CKM3x3*ee*complex(0,1))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_195 = Coupling(name = 'GC_195',
                  value = '-((cw*ee*complex(0,1))/sw)',
                  order = {'QED':1})

GC_196 = Coupling(name = 'GC_196',
                  value = '(cw*ee*complex(0,1))/sw',
                  order = {'QED':1})

GC_197 = Coupling(name = 'GC_197',
                  value = '(-2*cw*ee**2*complex(0,1))/sw',
                  order = {'QED':2})

GC_198 = Coupling(name = 'GC_198',
                  value = '-((ee**2*complex(0,1)*I105a11)/(sw*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_199 = Coupling(name = 'GC_199',
                  value = '-((ee**2*complex(0,1)*I105a22)/(sw*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_200 = Coupling(name = 'GC_200',
                  value = '-((ee**2*complex(0,1)*I105a33)/(sw*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_201 = Coupling(name = 'GC_201',
                  value = '-((ee**2*complex(0,1)*I105a36)/(sw*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_202 = Coupling(name = 'GC_202',
                  value = '(ee**2*complex(0,1)*I115a11)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_203 = Coupling(name = 'GC_203',
                  value = '(ee*complex(0,1)*G*I115a11*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_204 = Coupling(name = 'GC_204',
                  value = '(ee**2*complex(0,1)*I115a22)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_205 = Coupling(name = 'GC_205',
                  value = '(ee*complex(0,1)*G*I115a22*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_206 = Coupling(name = 'GC_206',
                  value = '(ee**2*complex(0,1)*I115a33)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_207 = Coupling(name = 'GC_207',
                  value = '(ee*complex(0,1)*G*I115a33*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_208 = Coupling(name = 'GC_208',
                  value = '(ee**2*complex(0,1)*I115a36)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_209 = Coupling(name = 'GC_209',
                  value = '(ee*complex(0,1)*G*I115a36*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_210 = Coupling(name = 'GC_210',
                  value = '(ee**2*complex(0,1)*I115a63)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_211 = Coupling(name = 'GC_211',
                  value = '(ee*complex(0,1)*G*I115a63*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_212 = Coupling(name = 'GC_212',
                  value = '(ee**2*complex(0,1)*I115a66)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_213 = Coupling(name = 'GC_213',
                  value = '(ee*complex(0,1)*G*I115a66*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_214 = Coupling(name = 'GC_214',
                  value = '(ee**2*complex(0,1)*I136a11)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_215 = Coupling(name = 'GC_215',
                  value = '(ee*complex(0,1)*G*I136a11*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_216 = Coupling(name = 'GC_216',
                  value = '(ee**2*complex(0,1)*I136a22)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_217 = Coupling(name = 'GC_217',
                  value = '(ee*complex(0,1)*G*I136a22*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_218 = Coupling(name = 'GC_218',
                  value = '(ee**2*complex(0,1)*I136a33)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_219 = Coupling(name = 'GC_219',
                  value = '(ee*complex(0,1)*G*I136a33*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_220 = Coupling(name = 'GC_220',
                  value = '(ee**2*complex(0,1)*I136a36)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_221 = Coupling(name = 'GC_221',
                  value = '(ee*complex(0,1)*G*I136a36*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_222 = Coupling(name = 'GC_222',
                  value = '(ee**2*complex(0,1)*I136a63)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_223 = Coupling(name = 'GC_223',
                  value = '(ee*complex(0,1)*G*I136a63*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_224 = Coupling(name = 'GC_224',
                  value = '(ee**2*complex(0,1)*I136a66)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_225 = Coupling(name = 'GC_225',
                  value = '(ee*complex(0,1)*G*I136a66*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_226 = Coupling(name = 'GC_226',
                  value = '-((ee*complex(0,1)*I197a11)/(sw*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_227 = Coupling(name = 'GC_227',
                  value = '-((ee*complex(0,1)*I197a22)/(sw*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_228 = Coupling(name = 'GC_228',
                  value = '-((ee*complex(0,1)*I197a33)/(sw*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_229 = Coupling(name = 'GC_229',
                  value = '-((ee*complex(0,1)*I197a36)/(sw*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_230 = Coupling(name = 'GC_230',
                  value = '-((ee*complex(0,1)*I197a63)/(sw*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_231 = Coupling(name = 'GC_231',
                  value = '-((ee*complex(0,1)*I197a66)/(sw*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_232 = Coupling(name = 'GC_232',
                  value = '-((ee*complex(0,1)*I198a11)/(sw*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_233 = Coupling(name = 'GC_233',
                  value = '-((ee*complex(0,1)*I198a22)/(sw*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_234 = Coupling(name = 'GC_234',
                  value = '-((ee*complex(0,1)*I198a33)/(sw*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_235 = Coupling(name = 'GC_235',
                  value = '-((ee*complex(0,1)*I198a36)/(sw*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_236 = Coupling(name = 'GC_236',
                  value = '(ee*complex(0,1)*I199a11)/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_237 = Coupling(name = 'GC_237',
                  value = '(ee*complex(0,1)*I199a22)/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_238 = Coupling(name = 'GC_238',
                  value = '(ee*complex(0,1)*I199a33)/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_239 = Coupling(name = 'GC_239',
                  value = '(ee*complex(0,1)*I199a36)/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_240 = Coupling(name = 'GC_240',
                  value = '(ee*complex(0,1)*I199a63)/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_241 = Coupling(name = 'GC_241',
                  value = '(ee*complex(0,1)*I199a66)/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_242 = Coupling(name = 'GC_242',
                  value = '(ee*complex(0,1)*I200a11)/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_243 = Coupling(name = 'GC_243',
                  value = '(ee*complex(0,1)*I200a22)/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_244 = Coupling(name = 'GC_244',
                  value = '(ee*complex(0,1)*I200a33)/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_245 = Coupling(name = 'GC_245',
                  value = '(ee*complex(0,1)*I200a36)/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_246 = Coupling(name = 'GC_246',
                  value = '-((ee**2*complex(0,1)*I90a11)/(sw*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_247 = Coupling(name = 'GC_247',
                  value = '-((ee**2*complex(0,1)*I90a22)/(sw*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_248 = Coupling(name = 'GC_248',
                  value = '-((ee**2*complex(0,1)*I90a33)/(sw*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_249 = Coupling(name = 'GC_249',
                  value = '-((ee**2*complex(0,1)*I90a36)/(sw*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_250 = Coupling(name = 'GC_250',
                  value = '(ee**2*complex(0,1)*I149a44)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_251 = Coupling(name = 'GC_251',
                  value = '(ee**2*complex(0,1)*I149a55)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_252 = Coupling(name = 'GC_252',
                  value = '(4*ee**2*complex(0,1)*I172a44*I175a44)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_253 = Coupling(name = 'GC_253',
                  value = '(4*ee**2*complex(0,1)*I172a55*I175a55)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_254 = Coupling(name = 'GC_254',
                  value = '(4*ee**2*complex(0,1)*I177a44*I178a44)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_255 = Coupling(name = 'GC_255',
                  value = '(4*ee**2*complex(0,1)*I177a55*I178a55)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_256 = Coupling(name = 'GC_256',
                  value = '(cw*ee**2*complex(0,1)*I197a11)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_257 = Coupling(name = 'GC_257',
                  value = '(cw*ee**2*complex(0,1)*I197a22)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_258 = Coupling(name = 'GC_258',
                  value = '(cw*ee**2*complex(0,1)*I197a33)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_259 = Coupling(name = 'GC_259',
                  value = '(cw*ee**2*complex(0,1)*I197a36)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_260 = Coupling(name = 'GC_260',
                  value = '(cw*ee**2*complex(0,1)*I197a63)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_261 = Coupling(name = 'GC_261',
                  value = '(cw*ee**2*complex(0,1)*I197a66)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_262 = Coupling(name = 'GC_262',
                  value = '-((cw*ee**2*complex(0,1)*I198a11)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_263 = Coupling(name = 'GC_263',
                  value = '-((cw*ee**2*complex(0,1)*I198a22)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_264 = Coupling(name = 'GC_264',
                  value = '-((cw*ee**2*complex(0,1)*I198a33)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_265 = Coupling(name = 'GC_265',
                  value = '-((cw*ee**2*complex(0,1)*I198a36)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_266 = Coupling(name = 'GC_266',
                  value = '(cw*ee**2*complex(0,1)*I199a11)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_267 = Coupling(name = 'GC_267',
                  value = '(cw*ee**2*complex(0,1)*I199a22)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_268 = Coupling(name = 'GC_268',
                  value = '(cw*ee**2*complex(0,1)*I199a33)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_269 = Coupling(name = 'GC_269',
                  value = '(cw*ee**2*complex(0,1)*I199a36)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_270 = Coupling(name = 'GC_270',
                  value = '(cw*ee**2*complex(0,1)*I199a63)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_271 = Coupling(name = 'GC_271',
                  value = '(cw*ee**2*complex(0,1)*I199a66)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_272 = Coupling(name = 'GC_272',
                  value = '-((cw*ee**2*complex(0,1)*I200a11)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_273 = Coupling(name = 'GC_273',
                  value = '-((cw*ee**2*complex(0,1)*I200a22)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_274 = Coupling(name = 'GC_274',
                  value = '-((cw*ee**2*complex(0,1)*I200a33)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_275 = Coupling(name = 'GC_275',
                  value = '-((cw*ee**2*complex(0,1)*I200a36)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_276 = Coupling(name = 'GC_276',
                  value = '(ee**2*complex(0,1)*I30a44*I33a44)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_277 = Coupling(name = 'GC_277',
                  value = '(ee**2*complex(0,1)*I30a55*I33a55)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_278 = Coupling(name = 'GC_278',
                  value = '(ee**2*complex(0,1)*I35a44*I36a44)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_279 = Coupling(name = 'GC_279',
                  value = '(ee**2*complex(0,1)*I35a55*I36a55)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_280 = Coupling(name = 'GC_280',
                  value = '-(ee**2*complex(0,1)*I162a44*I63a11)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_281 = Coupling(name = 'GC_281',
                  value = '-(ee**2*complex(0,1)*I162a55*I63a11)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_282 = Coupling(name = 'GC_282',
                  value = '-(ee**2*complex(0,1)*I162a44*I63a22)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_283 = Coupling(name = 'GC_283',
                  value = '-(ee**2*complex(0,1)*I162a55*I63a22)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_284 = Coupling(name = 'GC_284',
                  value = '-(ee**2*complex(0,1)*I65a44)/(6.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_285 = Coupling(name = 'GC_285',
                  value = '(ee**2*complex(0,1)*I161a11*I65a44)/(18.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_286 = Coupling(name = 'GC_286',
                  value = '(ee**2*complex(0,1)*I161a22*I65a44)/(18.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_287 = Coupling(name = 'GC_287',
                  value = '(-2*ee**2*complex(0,1)*I162a44*I65a44)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_288 = Coupling(name = 'GC_288',
                  value = '(-2*ee**2*complex(0,1)*I162a55*I65a44)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_289 = Coupling(name = 'GC_289',
                  value = '-(ee**2*complex(0,1)*I64a11*I65a44)/(6.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_290 = Coupling(name = 'GC_290',
                  value = '-(ee**2*complex(0,1)*I64a22*I65a44)/(6.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_291 = Coupling(name = 'GC_291',
                  value = '-(ee**2*complex(0,1)*I65a55)/(6.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_292 = Coupling(name = 'GC_292',
                  value = '(ee**2*complex(0,1)*I161a11*I65a55)/(18.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_293 = Coupling(name = 'GC_293',
                  value = '(ee**2*complex(0,1)*I161a22*I65a55)/(18.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_294 = Coupling(name = 'GC_294',
                  value = '(-2*ee**2*complex(0,1)*I162a44*I65a55)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_295 = Coupling(name = 'GC_295',
                  value = '(-2*ee**2*complex(0,1)*I162a55*I65a55)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_296 = Coupling(name = 'GC_296',
                  value = '-(ee**2*complex(0,1)*I64a11*I65a55)/(6.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_297 = Coupling(name = 'GC_297',
                  value = '-(ee**2*complex(0,1)*I64a22*I65a55)/(6.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_298 = Coupling(name = 'GC_298',
                  value = '(ee**2*complex(0,1)*I63a11*I66a44)/(6.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_299 = Coupling(name = 'GC_299',
                  value = '(ee**2*complex(0,1)*I63a22*I66a44)/(6.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_300 = Coupling(name = 'GC_300',
                  value = '(ee**2*complex(0,1)*I65a44*I66a44)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_301 = Coupling(name = 'GC_301',
                  value = '(ee**2*complex(0,1)*I65a55*I66a44)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_302 = Coupling(name = 'GC_302',
                  value = '(ee**2*complex(0,1)*I63a11*I66a55)/(6.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_303 = Coupling(name = 'GC_303',
                  value = '(ee**2*complex(0,1)*I63a22*I66a55)/(6.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_304 = Coupling(name = 'GC_304',
                  value = '(ee**2*complex(0,1)*I65a44*I66a55)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_305 = Coupling(name = 'GC_305',
                  value = '(ee**2*complex(0,1)*I65a55*I66a55)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_306 = Coupling(name = 'GC_306',
                  value = '(ee**2*complex(0,1)*I162a44*I97a11)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_307 = Coupling(name = 'GC_307',
                  value = '(ee**2*complex(0,1)*I162a55*I97a11)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_308 = Coupling(name = 'GC_308',
                  value = '(ee**2*complex(0,1)*I162a44*I97a22)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_309 = Coupling(name = 'GC_309',
                  value = '(ee**2*complex(0,1)*I162a55*I97a22)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_310 = Coupling(name = 'GC_310',
                  value = '-(ee**2*complex(0,1)*I98a44)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_311 = Coupling(name = 'GC_311',
                  value = '(ee**2*complex(0,1)*I161a11*I98a44)/(6.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_312 = Coupling(name = 'GC_312',
                  value = '(ee**2*complex(0,1)*I161a22*I98a44)/(6.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_313 = Coupling(name = 'GC_313',
                  value = '(-2*ee**2*complex(0,1)*I162a44*I98a44)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_314 = Coupling(name = 'GC_314',
                  value = '(-2*ee**2*complex(0,1)*I162a55*I98a44)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_315 = Coupling(name = 'GC_315',
                  value = '-(ee**2*complex(0,1)*I98a55)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_316 = Coupling(name = 'GC_316',
                  value = '(ee**2*complex(0,1)*I161a11*I98a55)/(6.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_317 = Coupling(name = 'GC_317',
                  value = '(ee**2*complex(0,1)*I161a22*I98a55)/(6.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_318 = Coupling(name = 'GC_318',
                  value = '(-2*ee**2*complex(0,1)*I162a44*I98a55)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_319 = Coupling(name = 'GC_319',
                  value = '(-2*ee**2*complex(0,1)*I162a55*I98a55)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_320 = Coupling(name = 'GC_320',
                  value = '(cw*ee*complex(0,1)*NN1x1*Rd4x4*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_321 = Coupling(name = 'GC_321',
                  value = '(cw*ee*complex(0,1)*NN2x1*Rd4x4*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_322 = Coupling(name = 'GC_322',
                  value = '(cw*ee*complex(0,1)*NN3x1*Rd4x4*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_323 = Coupling(name = 'GC_323',
                  value = '(cw*ee*complex(0,1)*NN4x1*Rd4x4*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_324 = Coupling(name = 'GC_324',
                  value = '(cw*ee*complex(0,1)*NN1x1*Rd5x5*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_325 = Coupling(name = 'GC_325',
                  value = '(cw*ee*complex(0,1)*NN2x1*Rd5x5*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_326 = Coupling(name = 'GC_326',
                  value = '(cw*ee*complex(0,1)*NN3x1*Rd5x5*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_327 = Coupling(name = 'GC_327',
                  value = '(cw*ee*complex(0,1)*NN4x1*Rd5x5*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_328 = Coupling(name = 'GC_328',
                  value = '(cw*ee*complex(0,1)*NN1x1*Rl4x4*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_329 = Coupling(name = 'GC_329',
                  value = '(cw*ee*complex(0,1)*NN2x1*Rl4x4*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_330 = Coupling(name = 'GC_330',
                  value = '(cw*ee*complex(0,1)*NN3x1*Rl4x4*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_331 = Coupling(name = 'GC_331',
                  value = '(cw*ee*complex(0,1)*NN4x1*Rl4x4*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_332 = Coupling(name = 'GC_332',
                  value = '(cw*ee*complex(0,1)*NN1x1*Rl5x5*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_333 = Coupling(name = 'GC_333',
                  value = '(cw*ee*complex(0,1)*NN2x1*Rl5x5*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_334 = Coupling(name = 'GC_334',
                  value = '(cw*ee*complex(0,1)*NN3x1*Rl5x5*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_335 = Coupling(name = 'GC_335',
                  value = '(cw*ee*complex(0,1)*NN4x1*Rl5x5*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_336 = Coupling(name = 'GC_336',
                  value = '(-2*cw*ee*complex(0,1)*NN1x1*Ru4x4*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_337 = Coupling(name = 'GC_337',
                  value = '(-2*cw*ee*complex(0,1)*NN2x1*Ru4x4*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_338 = Coupling(name = 'GC_338',
                  value = '(-2*cw*ee*complex(0,1)*NN3x1*Ru4x4*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_339 = Coupling(name = 'GC_339',
                  value = '(-2*cw*ee*complex(0,1)*NN4x1*Ru4x4*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_340 = Coupling(name = 'GC_340',
                  value = '(-2*cw*ee*complex(0,1)*NN1x1*Ru5x5*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_341 = Coupling(name = 'GC_341',
                  value = '(-2*cw*ee*complex(0,1)*NN2x1*Ru5x5*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_342 = Coupling(name = 'GC_342',
                  value = '(-2*cw*ee*complex(0,1)*NN3x1*Ru5x5*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_343 = Coupling(name = 'GC_343',
                  value = '(-2*cw*ee*complex(0,1)*NN4x1*Ru5x5*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_344 = Coupling(name = 'GC_344',
                  value = '(ee**2*complex(0,1))/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_345 = Coupling(name = 'GC_345',
                  value = '-(ee**2*complex(0,1))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_346 = Coupling(name = 'GC_346',
                  value = '(ee**2*complex(0,1))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_347 = Coupling(name = 'GC_347',
                  value = '-(cw*ee*complex(0,1))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_348 = Coupling(name = 'GC_348',
                  value = '(cw*ee*complex(0,1))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_349 = Coupling(name = 'GC_349',
                  value = '-(cw*ee*complex(0,1)*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_350 = Coupling(name = 'GC_350',
                  value = '(2*cw*ee*complex(0,1)*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_351 = Coupling(name = 'GC_351',
                  value = '-((cw*ee*complex(0,1)*sw)/((-1 + sw)*(1 + sw)))',
                  order = {'QED':1})

GC_352 = Coupling(name = 'GC_352',
                  value = '(8*cw*ee**2*complex(0,1)*I114a44*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_353 = Coupling(name = 'GC_353',
                  value = '(4*cw*ee*complex(0,1)*G*I114a44*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_354 = Coupling(name = 'GC_354',
                  value = '(8*cw*ee**2*complex(0,1)*I114a55*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_355 = Coupling(name = 'GC_355',
                  value = '(4*cw*ee*complex(0,1)*G*I114a55*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_356 = Coupling(name = 'GC_356',
                  value = '(-2*cw*ee*complex(0,1)*I203a44*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_357 = Coupling(name = 'GC_357',
                  value = '(-2*cw*ee*complex(0,1)*I203a55*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_358 = Coupling(name = 'GC_358',
                  value = '(2*cw*ee**2*complex(0,1)*I48a44*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_359 = Coupling(name = 'GC_359',
                  value = '(2*cw*ee**2*complex(0,1)*I48a55*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_360 = Coupling(name = 'GC_360',
                  value = '(cw*ee*complex(0,1)*I65a44*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_361 = Coupling(name = 'GC_361',
                  value = '(cw*ee*complex(0,1)*I65a55*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_362 = Coupling(name = 'GC_362',
                  value = '(cw*ee*complex(0,1)*I98a44*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_363 = Coupling(name = 'GC_363',
                  value = '(cw*ee*complex(0,1)*I98a55*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_364 = Coupling(name = 'GC_364',
                  value = '(2*cw*ee**2*complex(0,1)*I9a44*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_365 = Coupling(name = 'GC_365',
                  value = '(-2*cw*ee*complex(0,1)*G*I9a44*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_366 = Coupling(name = 'GC_366',
                  value = '(2*cw*ee**2*complex(0,1)*I9a55*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_367 = Coupling(name = 'GC_367',
                  value = '(-2*cw*ee*complex(0,1)*G*I9a55*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_368 = Coupling(name = 'GC_368',
                  value = '(-8*ee**2*complex(0,1)*I203a44*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_369 = Coupling(name = 'GC_369',
                  value = '(-8*ee**2*complex(0,1)*I203a55*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_370 = Coupling(name = 'GC_370',
                  value = '(-2*ee**2*complex(0,1)*I65a44*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_371 = Coupling(name = 'GC_371',
                  value = '(-2*ee**2*complex(0,1)*I65a55*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_372 = Coupling(name = 'GC_372',
                  value = '(-2*ee**2*complex(0,1)*I98a44*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_373 = Coupling(name = 'GC_373',
                  value = '(-2*ee**2*complex(0,1)*I98a55*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_374 = Coupling(name = 'GC_374',
                  value = '-(ee**2*complex(0,1)*I169a11*I172a44)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a11*I174a44)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a11*I176a44)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a11*I177a44)/(36.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_375 = Coupling(name = 'GC_375',
                  value = '-(ee**2*complex(0,1)*I169a22*I172a44)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a22*I174a44)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a22*I176a44)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a22*I177a44)/(36.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_376 = Coupling(name = 'GC_376',
                  value = '-(ee**2*complex(0,1)*I169a11*I172a55)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a11*I174a55)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a11*I176a55)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a11*I177a55)/(36.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_377 = Coupling(name = 'GC_377',
                  value = '-(ee**2*complex(0,1)*I169a22*I172a55)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a22*I174a55)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a22*I176a55)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a22*I177a55)/(36.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_378 = Coupling(name = 'GC_378',
                  value = '-(ee**2*complex(0,1)*I169a33*I172a44)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a33*I174a44)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a44*I175a33)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a33*I175a44)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a33*I176a44)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a33*I177a44)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a44*I178a33)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a33*I178a44)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_379 = Coupling(name = 'GC_379',
                  value = '-(ee**2*complex(0,1)*I169a36*I172a44)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a36*I174a44)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a44*I175a36)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a36*I175a44)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a36*I176a44)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a36*I177a44)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a44*I178a36)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a36*I178a44)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_380 = Coupling(name = 'GC_380',
                  value = '-(ee**2*complex(0,1)*I169a33*I172a55)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a33*I174a55)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a55*I175a33)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a33*I175a55)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a33*I176a55)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a33*I177a55)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a55*I178a33)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a33*I178a55)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_381 = Coupling(name = 'GC_381',
                  value = '-(ee**2*complex(0,1)*I169a36*I172a55)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a36*I174a55)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a55*I175a36)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a36*I175a55)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a36*I176a55)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a36*I177a55)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a55*I178a36)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a36*I178a55)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_382 = Coupling(name = 'GC_382',
                  value = '(ee**2*complex(0,1)*I172a55*I175a44)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a44*I175a55)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a55*I178a44)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a44*I178a55)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_383 = Coupling(name = 'GC_383',
                  value = '-(ee**2*complex(0,1)*I169a63*I172a44)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a63*I174a44)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a63*I175a44)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a44*I175a63)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a63*I176a44)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a63*I177a44)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a63*I178a44)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a44*I178a63)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_384 = Coupling(name = 'GC_384',
                  value = '-(ee**2*complex(0,1)*I169a63*I172a55)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a63*I174a55)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a63*I175a55)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a55*I175a63)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a63*I176a55)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a63*I177a55)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a63*I178a55)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a55*I178a63)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_385 = Coupling(name = 'GC_385',
                  value = '-(ee**2*complex(0,1)*I169a66*I172a44)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a66*I174a44)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a66*I175a44)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a44*I175a66)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a66*I176a44)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a66*I177a44)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a66*I178a44)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a44*I178a66)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_386 = Coupling(name = 'GC_386',
                  value = '-(ee**2*complex(0,1)*I169a66*I172a55)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a66*I174a55)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a66*I175a55)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a55*I175a66)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a66*I176a55)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a66*I177a55)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a66*I178a55)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a55*I178a66)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_387 = Coupling(name = 'GC_387',
                  value = '(ee**2*complex(0,1)*I30a55*I33a44)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a44*I33a55)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a55*I36a44)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a44*I36a55)/(36.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_388 = Coupling(name = 'GC_388',
                  value = '-(ee**2*complex(0,1)*I162a44*I63a33)/(9.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a44*I65a33)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_389 = Coupling(name = 'GC_389',
                  value = '-(ee**2*complex(0,1)*I162a55*I63a33)/(9.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a55*I65a33)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_390 = Coupling(name = 'GC_390',
                  value = '-(ee**2*complex(0,1)*I162a44*I63a36)/(9.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a44*I65a36)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_391 = Coupling(name = 'GC_391',
                  value = '-(ee**2*complex(0,1)*I162a55*I63a36)/(9.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a55*I65a36)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_392 = Coupling(name = 'GC_392',
                  value = '(ee**2*complex(0,1)*I161a33*I65a44)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a33*I65a44)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_393 = Coupling(name = 'GC_393',
                  value = '(ee**2*complex(0,1)*I161a36*I65a44)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a36*I65a44)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_394 = Coupling(name = 'GC_394',
                  value = '(ee**2*complex(0,1)*I161a63*I65a44)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a63*I65a44)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_395 = Coupling(name = 'GC_395',
                  value = '(ee**2*complex(0,1)*I161a66*I65a44)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a66*I65a44)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_396 = Coupling(name = 'GC_396',
                  value = '(ee**2*complex(0,1)*I161a33*I65a55)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a33*I65a55)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_397 = Coupling(name = 'GC_397',
                  value = '(ee**2*complex(0,1)*I161a36*I65a55)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a36*I65a55)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_398 = Coupling(name = 'GC_398',
                  value = '(ee**2*complex(0,1)*I161a63*I65a55)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a63*I65a55)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_399 = Coupling(name = 'GC_399',
                  value = '(ee**2*complex(0,1)*I161a66*I65a55)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a66*I65a55)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_400 = Coupling(name = 'GC_400',
                  value = '-(ee**2*complex(0,1)*I162a44*I63a63)/(9.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a44*I65a63)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_401 = Coupling(name = 'GC_401',
                  value = '-(ee**2*complex(0,1)*I162a55*I63a63)/(9.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a55*I65a63)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_402 = Coupling(name = 'GC_402',
                  value = '-(ee**2*complex(0,1)*I162a44*I63a66)/(9.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a44*I65a66)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_403 = Coupling(name = 'GC_403',
                  value = '-(ee**2*complex(0,1)*I162a55*I63a66)/(9.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a55*I65a66)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_404 = Coupling(name = 'GC_404',
                  value = '-(ee**2*complex(0,1)*I64a33*I65a44)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a44*I66a33)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_405 = Coupling(name = 'GC_405',
                  value = '-(ee**2*complex(0,1)*I64a33*I65a55)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a55*I66a33)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_406 = Coupling(name = 'GC_406',
                  value = '-(ee**2*complex(0,1)*I64a36*I65a44)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a44*I66a36)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_407 = Coupling(name = 'GC_407',
                  value = '-(ee**2*complex(0,1)*I64a36*I65a55)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a55*I66a36)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_408 = Coupling(name = 'GC_408',
                  value = '(ee**2*complex(0,1)*I63a33*I66a44)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a33*I66a44)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_409 = Coupling(name = 'GC_409',
                  value = '(ee**2*complex(0,1)*I63a36*I66a44)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a36*I66a44)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_410 = Coupling(name = 'GC_410',
                  value = '(ee**2*complex(0,1)*I63a63*I66a44)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a63*I66a44)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_411 = Coupling(name = 'GC_411',
                  value = '(ee**2*complex(0,1)*I63a66*I66a44)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a66*I66a44)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_412 = Coupling(name = 'GC_412',
                  value = '(ee**2*complex(0,1)*I63a33*I66a55)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a33*I66a55)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_413 = Coupling(name = 'GC_413',
                  value = '(ee**2*complex(0,1)*I63a36*I66a55)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a36*I66a55)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_414 = Coupling(name = 'GC_414',
                  value = '(ee**2*complex(0,1)*I63a63*I66a55)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a63*I66a55)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_415 = Coupling(name = 'GC_415',
                  value = '(ee**2*complex(0,1)*I63a66*I66a55)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a66*I66a55)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_416 = Coupling(name = 'GC_416',
                  value = '-(ee**2*complex(0,1)*I64a63*I65a44)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a44*I66a63)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_417 = Coupling(name = 'GC_417',
                  value = '-(ee**2*complex(0,1)*I64a63*I65a55)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a55*I66a63)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_418 = Coupling(name = 'GC_418',
                  value = '-(ee**2*complex(0,1)*I64a66*I65a44)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a44*I66a66)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_419 = Coupling(name = 'GC_419',
                  value = '-(ee**2*complex(0,1)*I64a66*I65a55)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a55*I66a66)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_420 = Coupling(name = 'GC_420',
                  value = '-(ee**2*complex(0,1)*I75a11*I76a44)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a11*I77a44)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a11*I79a44)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a11*I80a44)/(8.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_421 = Coupling(name = 'GC_421',
                  value = '-(ee**2*complex(0,1)*I75a22*I76a44)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a22*I77a44)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a22*I79a44)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a22*I80a44)/(8.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_422 = Coupling(name = 'GC_422',
                  value = '-(ee**2*complex(0,1)*I75a11*I76a55)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a11*I77a55)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a11*I79a55)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a11*I80a55)/(8.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_423 = Coupling(name = 'GC_423',
                  value = '-(ee**2*complex(0,1)*I75a22*I76a55)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a22*I77a55)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a22*I79a55)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a22*I80a55)/(8.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_424 = Coupling(name = 'GC_424',
                  value = '-(ee**2*complex(0,1)*I75a33*I76a44)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a33*I77a44)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a44*I78a33)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a33*I78a44)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a33*I79a44)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a33*I80a44)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a44*I81a33)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a33*I81a44)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_425 = Coupling(name = 'GC_425',
                  value = '-(ee**2*complex(0,1)*I75a36*I76a44)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a36*I77a44)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a44*I78a36)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a36*I78a44)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a36*I79a44)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a36*I80a44)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a44*I81a36)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a36*I81a44)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_426 = Coupling(name = 'GC_426',
                  value = '(ee**2*complex(0,1)*I77a44*I78a44)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a44*I81a44)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_427 = Coupling(name = 'GC_427',
                  value = '-(ee**2*complex(0,1)*I75a33*I76a55)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a33*I77a55)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a55*I78a33)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a33*I78a55)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a33*I79a55)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a33*I80a55)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a55*I81a33)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a33*I81a55)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_428 = Coupling(name = 'GC_428',
                  value = '-(ee**2*complex(0,1)*I75a36*I76a55)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a36*I77a55)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a55*I78a36)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a36*I78a55)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a36*I79a55)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a36*I80a55)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a55*I81a36)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a36*I81a55)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_429 = Coupling(name = 'GC_429',
                  value = '(ee**2*complex(0,1)*I77a55*I78a44)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a44*I78a55)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a55*I81a44)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a44*I81a55)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_430 = Coupling(name = 'GC_430',
                  value = '(ee**2*complex(0,1)*I77a55*I78a55)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a55*I81a55)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_431 = Coupling(name = 'GC_431',
                  value = '-(ee**2*complex(0,1)*I75a63*I76a44)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a63*I77a44)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a63*I78a44)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a44*I78a63)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a63*I79a44)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a63*I80a44)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a63*I81a44)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a44*I81a63)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_432 = Coupling(name = 'GC_432',
                  value = '-(ee**2*complex(0,1)*I75a63*I76a55)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a63*I77a55)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a63*I78a55)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a55*I78a63)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a63*I79a55)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a63*I80a55)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a63*I81a55)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a55*I81a63)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_433 = Coupling(name = 'GC_433',
                  value = '-(ee**2*complex(0,1)*I75a66*I76a44)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a66*I77a44)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a66*I78a44)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a44*I78a66)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a66*I79a44)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a66*I80a44)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a66*I81a44)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a44*I81a66)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_434 = Coupling(name = 'GC_434',
                  value = '-(ee**2*complex(0,1)*I75a66*I76a55)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a66*I77a55)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a66*I78a55)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a55*I78a66)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a66*I79a55)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a66*I80a55)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a66*I81a55)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a55*I81a66)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_435 = Coupling(name = 'GC_435',
                  value = '(ee**2*complex(0,1)*I27a11*I30a44)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a11*I32a44)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a11*I35a44)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a44*I8a11)/(72.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_436 = Coupling(name = 'GC_436',
                  value = '(ee**2*complex(0,1)*I27a11*I30a55)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a11*I32a55)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a11*I35a55)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a55*I8a11)/(72.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_437 = Coupling(name = 'GC_437',
                  value = '(ee**2*complex(0,1)*I27a22*I30a44)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a22*I32a44)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a22*I35a44)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a44*I8a22)/(72.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_438 = Coupling(name = 'GC_438',
                  value = '(ee**2*complex(0,1)*I27a22*I30a55)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a22*I32a55)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a22*I35a55)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a55*I8a22)/(72.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_439 = Coupling(name = 'GC_439',
                  value = '(ee**2*complex(0,1)*I27a33*I30a44)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a33*I32a44)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a44*I33a33)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a33*I33a44)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a33*I35a44)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a44*I36a33)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a33*I36a44)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a44*I8a33)/(72.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_440 = Coupling(name = 'GC_440',
                  value = '(ee**2*complex(0,1)*I27a33*I30a55)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a33*I32a55)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a55*I33a33)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a33*I33a55)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a33*I35a55)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a55*I36a33)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a33*I36a55)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a55*I8a33)/(72.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_441 = Coupling(name = 'GC_441',
                  value = '(ee**2*complex(0,1)*I27a36*I30a44)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a36*I32a44)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a44*I33a36)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a36*I33a44)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a36*I35a44)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a44*I36a36)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a36*I36a44)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a44*I8a36)/(72.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_442 = Coupling(name = 'GC_442',
                  value = '(ee**2*complex(0,1)*I27a36*I30a55)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a36*I32a55)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a55*I33a36)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a36*I33a55)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a36*I35a55)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a55*I36a36)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a36*I36a55)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a55*I8a36)/(72.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_443 = Coupling(name = 'GC_443',
                  value = '(ee**2*complex(0,1)*I27a63*I30a44)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a63*I32a44)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a63*I33a44)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a44*I33a63)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a63*I35a44)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a63*I36a44)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a44*I36a63)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a44*I8a63)/(72.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_444 = Coupling(name = 'GC_444',
                  value = '(ee**2*complex(0,1)*I27a63*I30a55)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a63*I32a55)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a63*I33a55)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a55*I33a63)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a63*I35a55)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a63*I36a55)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a55*I36a63)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a55*I8a63)/(72.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_445 = Coupling(name = 'GC_445',
                  value = '(ee**2*complex(0,1)*I27a66*I30a44)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a66*I32a44)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a66*I33a44)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a44*I33a66)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a66*I35a44)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a66*I36a44)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a44*I36a66)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a44*I8a66)/(72.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_446 = Coupling(name = 'GC_446',
                  value = '(ee**2*complex(0,1)*I27a66*I30a55)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a66*I32a55)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a66*I33a55)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a55*I33a66)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a66*I35a55)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a66*I36a55)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a55*I36a66)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a55*I8a66)/(72.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_447 = Coupling(name = 'GC_447',
                  value = '(ee**2*complex(0,1)*I162a44*I97a33)/(3.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a44*I98a33)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_448 = Coupling(name = 'GC_448',
                  value = '(ee**2*complex(0,1)*I162a55*I97a33)/(3.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a55*I98a33)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_449 = Coupling(name = 'GC_449',
                  value = '(ee**2*complex(0,1)*I162a44*I97a36)/(3.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a44*I98a36)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_450 = Coupling(name = 'GC_450',
                  value = '(ee**2*complex(0,1)*I162a55*I97a36)/(3.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a55*I98a36)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_451 = Coupling(name = 'GC_451',
                  value = '(ee**2*complex(0,1)*I161a33*I98a44)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a33*I98a44)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_452 = Coupling(name = 'GC_452',
                  value = '(ee**2*complex(0,1)*I161a36*I98a44)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a36*I98a44)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_453 = Coupling(name = 'GC_453',
                  value = '(ee**2*complex(0,1)*I161a63*I98a44)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a63*I98a44)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_454 = Coupling(name = 'GC_454',
                  value = '(ee**2*complex(0,1)*I161a66*I98a44)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a66*I98a44)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_455 = Coupling(name = 'GC_455',
                  value = '(ee**2*complex(0,1)*I161a33*I98a55)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a33*I98a55)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_456 = Coupling(name = 'GC_456',
                  value = '(ee**2*complex(0,1)*I161a36*I98a55)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a36*I98a55)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_457 = Coupling(name = 'GC_457',
                  value = '(ee**2*complex(0,1)*I161a63*I98a55)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a63*I98a55)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_458 = Coupling(name = 'GC_458',
                  value = '(ee**2*complex(0,1)*I161a66*I98a55)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a66*I98a55)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_459 = Coupling(name = 'GC_459',
                  value = '(ee**2*complex(0,1)*I162a44*I97a63)/(3.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a44*I98a63)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_460 = Coupling(name = 'GC_460',
                  value = '(ee**2*complex(0,1)*I162a55*I97a63)/(3.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a55*I98a63)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_461 = Coupling(name = 'GC_461',
                  value = '(ee**2*complex(0,1)*I162a44*I97a66)/(3.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a44*I98a66)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_462 = Coupling(name = 'GC_462',
                  value = '(ee**2*complex(0,1)*I162a55*I97a66)/(3.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a55*I98a66)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_463 = Coupling(name = 'GC_463',
                  value = '-(ee**2*complex(0,1)*I148a11)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_464 = Coupling(name = 'GC_464',
                  value = '-(ee**2*complex(0,1)*I148a22)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_465 = Coupling(name = 'GC_465',
                  value = '-(ee**2*complex(0,1)*I148a33)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a33)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a33)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_466 = Coupling(name = 'GC_466',
                  value = '-(ee**2*complex(0,1)*I148a36)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a36)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a36)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_467 = Coupling(name = 'GC_467',
                  value = '-(ee**2*complex(0,1)*I148a63)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a63)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a63)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_468 = Coupling(name = 'GC_468',
                  value = '-(ee**2*complex(0,1)*I148a66)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a66)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a66)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_469 = Coupling(name = 'GC_469',
                  value = '-(ee**2*complex(0,1)*I159a11*I163a11)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a11*I163a11)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_470 = Coupling(name = 'GC_470',
                  value = '-(ee**2*complex(0,1)*I159a22*I163a11)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a22*I163a11)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_471 = Coupling(name = 'GC_471',
                  value = '-(ee**2*complex(0,1)*I159a33*I163a11)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a33*I163a11)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_472 = Coupling(name = 'GC_472',
                  value = '-(ee**2*complex(0,1)*I159a36*I163a11)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a36*I163a11)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_473 = Coupling(name = 'GC_473',
                  value = '-(ee**2*complex(0,1)*I159a63*I163a11)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a63*I163a11)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_474 = Coupling(name = 'GC_474',
                  value = '-(ee**2*complex(0,1)*I159a66*I163a11)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a66*I163a11)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_475 = Coupling(name = 'GC_475',
                  value = '-(ee**2*complex(0,1)*I159a11*I163a22)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a11*I163a22)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_476 = Coupling(name = 'GC_476',
                  value = '-(ee**2*complex(0,1)*I159a22*I163a22)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a22*I163a22)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_477 = Coupling(name = 'GC_477',
                  value = '-(ee**2*complex(0,1)*I159a33*I163a22)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a33*I163a22)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_478 = Coupling(name = 'GC_478',
                  value = '-(ee**2*complex(0,1)*I159a36*I163a22)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a36*I163a22)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_479 = Coupling(name = 'GC_479',
                  value = '-(ee**2*complex(0,1)*I159a63*I163a22)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a63*I163a22)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_480 = Coupling(name = 'GC_480',
                  value = '-(ee**2*complex(0,1)*I159a66*I163a22)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a66*I163a22)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_481 = Coupling(name = 'GC_481',
                  value = '-(ee**2*complex(0,1)*I159a11*I163a33)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a11*I163a33)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_482 = Coupling(name = 'GC_482',
                  value = '-(ee**2*complex(0,1)*I159a22*I163a33)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a22*I163a33)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_483 = Coupling(name = 'GC_483',
                  value = '-(ee**2*complex(0,1)*I159a11*I163a36)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a11*I163a36)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_484 = Coupling(name = 'GC_484',
                  value = '-(ee**2*complex(0,1)*I159a22*I163a36)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a22*I163a36)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_485 = Coupling(name = 'GC_485',
                  value = '-(ee**2*complex(0,1)*I159a11*I163a63)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a11*I163a63)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_486 = Coupling(name = 'GC_486',
                  value = '-(ee**2*complex(0,1)*I159a22*I163a63)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a22*I163a63)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_487 = Coupling(name = 'GC_487',
                  value = '-(ee**2*complex(0,1)*I159a11*I163a66)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a11*I163a66)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_488 = Coupling(name = 'GC_488',
                  value = '-(ee**2*complex(0,1)*I159a22*I163a66)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a22*I163a66)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_489 = Coupling(name = 'GC_489',
                  value = '(-2*ee**2*complex(0,1)*I168a11*I169a11)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a11*I169a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_490 = Coupling(name = 'GC_490',
                  value = '(-2*ee**2*complex(0,1)*I168a22*I169a22)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a22*I169a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_491 = Coupling(name = 'GC_491',
                  value = '(-2*ee**2*complex(0,1)*I170a11*I171a11)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I170a11*I171a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_492 = Coupling(name = 'GC_492',
                  value = '-(ee**2*complex(0,1)*I168a22*I169a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I168a11*I169a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a22*I171a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a11*I171a22)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a22*I169a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I168a11*I169a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a22*I171a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a11*I171a22)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_493 = Coupling(name = 'GC_493',
                  value = '(-2*ee**2*complex(0,1)*I170a22*I171a22)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I170a22*I171a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_494 = Coupling(name = 'GC_494',
                  value = '-(ee**2*complex(0,1)*I168a33*I169a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I168a11*I169a33)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a33*I171a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a11*I171a33)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a11*I172a33)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a11*I174a33)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a11*I176a33)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a11*I177a33)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a33*I169a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I168a11*I169a33)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a33*I171a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a11*I171a33)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_495 = Coupling(name = 'GC_495',
                  value = '-(ee**2*complex(0,1)*I168a33*I169a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I168a22*I169a33)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a33*I171a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a22*I171a33)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a22*I172a33)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a22*I174a33)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a22*I176a33)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a22*I177a33)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a33*I169a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I168a22*I169a33)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a33*I171a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a22*I171a33)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_496 = Coupling(name = 'GC_496',
                  value = '-(ee**2*complex(0,1)*I168a36*I169a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I168a11*I169a36)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a36*I171a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a11*I171a36)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a11*I172a36)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a11*I174a36)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a11*I176a36)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a11*I177a36)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a36*I169a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I168a11*I169a36)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a36*I171a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a11*I171a36)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_497 = Coupling(name = 'GC_497',
                  value = '-(ee**2*complex(0,1)*I168a36*I169a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I168a22*I169a36)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a36*I171a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a22*I171a36)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a22*I172a36)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a22*I174a36)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a22*I176a36)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a22*I177a36)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a36*I169a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I168a22*I169a36)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a36*I171a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a22*I171a36)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_498 = Coupling(name = 'GC_498',
                  value = '-(ee**2*complex(0,1)*I168a63*I169a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I168a11*I169a63)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a63*I171a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a11*I171a63)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a11*I172a63)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a11*I174a63)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a11*I176a63)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a11*I177a63)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a63*I169a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I168a11*I169a63)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a63*I171a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a11*I171a63)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_499 = Coupling(name = 'GC_499',
                  value = '-(ee**2*complex(0,1)*I168a63*I169a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I168a22*I169a63)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a63*I171a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a22*I171a63)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a22*I172a63)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a22*I174a63)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a22*I176a63)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a22*I177a63)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a63*I169a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I168a22*I169a63)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a63*I171a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a22*I171a63)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_500 = Coupling(name = 'GC_500',
                  value = '-(ee**2*complex(0,1)*I168a66*I169a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I168a11*I169a66)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a66*I171a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a11*I171a66)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a11*I172a66)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a11*I174a66)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a11*I176a66)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a11*I177a66)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a66*I169a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I168a11*I169a66)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a66*I171a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a11*I171a66)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_501 = Coupling(name = 'GC_501',
                  value = '-(ee**2*complex(0,1)*I168a66*I169a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I168a22*I169a66)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a66*I171a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a22*I171a66)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a22*I172a66)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a22*I174a66)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a22*I176a66)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a22*I177a66)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a66*I169a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I168a22*I169a66)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a66*I171a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a22*I171a66)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_502 = Coupling(name = 'GC_502',
                  value = '(-2*ee**2*complex(0,1)*I26a11*I27a11)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a11*I27a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_503 = Coupling(name = 'GC_503',
                  value = '(-2*ee**2*complex(0,1)*I26a22*I27a22)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a22*I27a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_504 = Coupling(name = 'GC_504',
                  value = '(-2*ee**2*complex(0,1)*I28a11*I29a11)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a11*I29a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_505 = Coupling(name = 'GC_505',
                  value = '-(ee**2*complex(0,1)*I26a22*I27a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I26a11*I27a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a22*I29a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a11*I29a22)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a22*I27a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I26a11*I27a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a22*I29a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a11*I29a22)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_506 = Coupling(name = 'GC_506',
                  value = '(-2*ee**2*complex(0,1)*I28a22*I29a22)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a22*I29a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_507 = Coupling(name = 'GC_507',
                  value = '-(ee**2*complex(0,1)*I26a33*I27a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I26a11*I27a33)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a33*I29a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a11*I29a33)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a11*I30a33)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a11*I32a33)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a11*I35a33)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a33*I8a11)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a33*I27a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I26a11*I27a33)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a33*I29a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a11*I29a33)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_508 = Coupling(name = 'GC_508',
                  value = '-(ee**2*complex(0,1)*I26a33*I27a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I26a22*I27a33)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a33*I29a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a22*I29a33)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a22*I30a33)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a22*I32a33)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a22*I35a33)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a33*I8a22)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a33*I27a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I26a22*I27a33)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a33*I29a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a22*I29a33)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_509 = Coupling(name = 'GC_509',
                  value = '-(ee**2*complex(0,1)*I26a36*I27a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I26a11*I27a36)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a36*I29a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a11*I29a36)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a11*I30a36)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a11*I32a36)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a11*I35a36)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a36*I8a11)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a36*I27a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I26a11*I27a36)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a36*I29a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a11*I29a36)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_510 = Coupling(name = 'GC_510',
                  value = '-(ee**2*complex(0,1)*I26a36*I27a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I26a22*I27a36)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a36*I29a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a22*I29a36)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a22*I30a36)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a22*I32a36)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a22*I35a36)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a36*I8a22)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a36*I27a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I26a22*I27a36)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a36*I29a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a22*I29a36)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_511 = Coupling(name = 'GC_511',
                  value = '-(ee**2*complex(0,1)*I26a63*I27a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I26a11*I27a63)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a63*I29a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a11*I29a63)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a11*I30a63)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a11*I32a63)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a11*I35a63)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a63*I8a11)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a63*I27a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I26a11*I27a63)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a63*I29a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a11*I29a63)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_512 = Coupling(name = 'GC_512',
                  value = '-(ee**2*complex(0,1)*I26a63*I27a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I26a22*I27a63)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a63*I29a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a22*I29a63)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a22*I30a63)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a22*I32a63)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a22*I35a63)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a63*I8a22)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a63*I27a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I26a22*I27a63)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a63*I29a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a22*I29a63)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_513 = Coupling(name = 'GC_513',
                  value = '-(ee**2*complex(0,1)*I26a66*I27a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I26a11*I27a66)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a66*I29a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a11*I29a66)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a11*I30a66)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a11*I32a66)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a11*I35a66)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a66*I8a11)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a66*I27a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I26a11*I27a66)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a66*I29a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a11*I29a66)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_514 = Coupling(name = 'GC_514',
                  value = '-(ee**2*complex(0,1)*I26a66*I27a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I26a22*I27a66)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a66*I29a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a22*I29a66)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a22*I30a66)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a22*I32a66)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a22*I35a66)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a66*I8a22)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a66*I27a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I26a22*I27a66)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a66*I29a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a22*I29a66)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_515 = Coupling(name = 'GC_515',
                  value = '(ee**2*complex(0,1)*I63a11)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I63a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_516 = Coupling(name = 'GC_516',
                  value = '(5*ee**2*complex(0,1)*I161a11*I63a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a11*I63a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_517 = Coupling(name = 'GC_517',
                  value = '(5*ee**2*complex(0,1)*I161a22*I63a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a22*I63a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_518 = Coupling(name = 'GC_518',
                  value = '(5*ee**2*complex(0,1)*I161a33*I63a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a33*I63a11)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a33*I63a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_519 = Coupling(name = 'GC_519',
                  value = '(5*ee**2*complex(0,1)*I161a36*I63a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a36*I63a11)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a36*I63a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_520 = Coupling(name = 'GC_520',
                  value = '(5*ee**2*complex(0,1)*I161a63*I63a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a63*I63a11)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a63*I63a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_521 = Coupling(name = 'GC_521',
                  value = '(5*ee**2*complex(0,1)*I161a66*I63a11)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a66*I63a11)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a66*I63a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_522 = Coupling(name = 'GC_522',
                  value = '(ee**2*complex(0,1)*I63a22)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I63a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_523 = Coupling(name = 'GC_523',
                  value = '(5*ee**2*complex(0,1)*I161a11*I63a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a11*I63a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_524 = Coupling(name = 'GC_524',
                  value = '(5*ee**2*complex(0,1)*I161a22*I63a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a22*I63a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_525 = Coupling(name = 'GC_525',
                  value = '(5*ee**2*complex(0,1)*I161a33*I63a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a33*I63a22)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a33*I63a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_526 = Coupling(name = 'GC_526',
                  value = '(5*ee**2*complex(0,1)*I161a36*I63a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a36*I63a22)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a36*I63a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_527 = Coupling(name = 'GC_527',
                  value = '(5*ee**2*complex(0,1)*I161a63*I63a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a63*I63a22)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a63*I63a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_528 = Coupling(name = 'GC_528',
                  value = '(5*ee**2*complex(0,1)*I161a66*I63a22)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a66*I63a22)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a66*I63a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_529 = Coupling(name = 'GC_529',
                  value = '(ee**2*complex(0,1)*I63a33)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I65a33)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I63a33)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_530 = Coupling(name = 'GC_530',
                  value = '(5*ee**2*complex(0,1)*I161a11*I63a33)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a11*I65a33)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a11*I63a33)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_531 = Coupling(name = 'GC_531',
                  value = '(5*ee**2*complex(0,1)*I161a22*I63a33)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a22*I65a33)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a22*I63a33)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_532 = Coupling(name = 'GC_532',
                  value = '(5*ee**2*complex(0,1)*I161a33*I63a33)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a33*I63a33)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a33*I65a33)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a33*I65a33)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a33*I63a33)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_533 = Coupling(name = 'GC_533',
                  value = '(5*ee**2*complex(0,1)*I161a36*I63a33)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a36*I63a33)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a36*I65a33)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a36*I65a33)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a36*I63a33)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_534 = Coupling(name = 'GC_534',
                  value = '(5*ee**2*complex(0,1)*I161a63*I63a33)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a63*I63a33)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a63*I65a33)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a63*I65a33)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a63*I63a33)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_535 = Coupling(name = 'GC_535',
                  value = '(5*ee**2*complex(0,1)*I161a66*I63a33)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a66*I63a33)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a66*I65a33)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a66*I65a33)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a66*I63a33)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_536 = Coupling(name = 'GC_536',
                  value = '(ee**2*complex(0,1)*I63a36)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I65a36)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I63a36)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_537 = Coupling(name = 'GC_537',
                  value = '(5*ee**2*complex(0,1)*I161a11*I63a36)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a11*I65a36)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a11*I63a36)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_538 = Coupling(name = 'GC_538',
                  value = '(5*ee**2*complex(0,1)*I161a22*I63a36)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a22*I65a36)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a22*I63a36)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_539 = Coupling(name = 'GC_539',
                  value = '(5*ee**2*complex(0,1)*I161a33*I63a36)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a33*I63a36)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a33*I65a36)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a33*I65a36)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a33*I63a36)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_540 = Coupling(name = 'GC_540',
                  value = '(5*ee**2*complex(0,1)*I161a36*I63a36)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a36*I63a36)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a36*I65a36)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a36*I65a36)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a36*I63a36)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_541 = Coupling(name = 'GC_541',
                  value = '(5*ee**2*complex(0,1)*I161a63*I63a36)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a63*I63a36)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a63*I65a36)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a63*I65a36)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a63*I63a36)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_542 = Coupling(name = 'GC_542',
                  value = '(5*ee**2*complex(0,1)*I161a66*I63a36)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a66*I63a36)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a66*I65a36)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a66*I65a36)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a66*I63a36)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_543 = Coupling(name = 'GC_543',
                  value = '(ee**2*complex(0,1)*I63a63)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I65a63)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I63a63)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_544 = Coupling(name = 'GC_544',
                  value = '(5*ee**2*complex(0,1)*I161a11*I63a63)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a11*I65a63)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a11*I63a63)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_545 = Coupling(name = 'GC_545',
                  value = '(5*ee**2*complex(0,1)*I161a22*I63a63)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a22*I65a63)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a22*I63a63)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_546 = Coupling(name = 'GC_546',
                  value = '(5*ee**2*complex(0,1)*I161a33*I63a63)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a33*I63a63)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a33*I65a63)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a33*I65a63)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a33*I63a63)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_547 = Coupling(name = 'GC_547',
                  value = '(5*ee**2*complex(0,1)*I161a36*I63a63)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a36*I63a63)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a36*I65a63)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a36*I65a63)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a36*I63a63)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_548 = Coupling(name = 'GC_548',
                  value = '(5*ee**2*complex(0,1)*I161a63*I63a63)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a63*I63a63)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a63*I65a63)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a63*I65a63)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a63*I63a63)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_549 = Coupling(name = 'GC_549',
                  value = '(5*ee**2*complex(0,1)*I161a66*I63a63)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a66*I63a63)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a66*I65a63)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a66*I65a63)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a66*I63a63)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_550 = Coupling(name = 'GC_550',
                  value = '(ee**2*complex(0,1)*I63a66)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I65a66)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I63a66)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_551 = Coupling(name = 'GC_551',
                  value = '(5*ee**2*complex(0,1)*I161a11*I63a66)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a11*I65a66)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a11*I63a66)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_552 = Coupling(name = 'GC_552',
                  value = '(5*ee**2*complex(0,1)*I161a22*I63a66)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a22*I65a66)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a22*I63a66)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_553 = Coupling(name = 'GC_553',
                  value = '(5*ee**2*complex(0,1)*I161a33*I63a66)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a33*I63a66)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a33*I65a66)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a33*I65a66)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a33*I63a66)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_554 = Coupling(name = 'GC_554',
                  value = '(5*ee**2*complex(0,1)*I161a36*I63a66)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a36*I63a66)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a36*I65a66)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a36*I65a66)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a36*I63a66)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_555 = Coupling(name = 'GC_555',
                  value = '(5*ee**2*complex(0,1)*I161a63*I63a66)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a63*I63a66)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a63*I65a66)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a63*I65a66)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a63*I63a66)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_556 = Coupling(name = 'GC_556',
                  value = '(5*ee**2*complex(0,1)*I161a66*I63a66)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I162a66*I63a66)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a66*I65a66)/(18.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a66*I65a66)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a66*I63a66)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_557 = Coupling(name = 'GC_557',
                  value = '-(ee**2*complex(0,1)*I63a11*I64a11)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a11*I64a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_558 = Coupling(name = 'GC_558',
                  value = '-(ee**2*complex(0,1)*I63a22*I64a11)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a22*I64a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_559 = Coupling(name = 'GC_559',
                  value = '-(ee**2*complex(0,1)*I63a33*I64a11)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a11*I65a33)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a33*I64a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_560 = Coupling(name = 'GC_560',
                  value = '-(ee**2*complex(0,1)*I63a36*I64a11)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a11*I65a36)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a36*I64a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_561 = Coupling(name = 'GC_561',
                  value = '-(ee**2*complex(0,1)*I63a63*I64a11)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a11*I65a63)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a63*I64a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_562 = Coupling(name = 'GC_562',
                  value = '-(ee**2*complex(0,1)*I63a66*I64a11)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a11*I65a66)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a66*I64a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_563 = Coupling(name = 'GC_563',
                  value = '-(ee**2*complex(0,1)*I63a11*I64a22)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a11*I64a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_564 = Coupling(name = 'GC_564',
                  value = '-(ee**2*complex(0,1)*I63a22*I64a22)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a22*I64a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_565 = Coupling(name = 'GC_565',
                  value = '-(ee**2*complex(0,1)*I63a33*I64a22)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a22*I65a33)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a33*I64a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_566 = Coupling(name = 'GC_566',
                  value = '-(ee**2*complex(0,1)*I63a36*I64a22)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a22*I65a36)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a36*I64a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_567 = Coupling(name = 'GC_567',
                  value = '-(ee**2*complex(0,1)*I63a63*I64a22)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a22*I65a63)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a63*I64a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_568 = Coupling(name = 'GC_568',
                  value = '-(ee**2*complex(0,1)*I63a66*I64a22)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a22*I65a66)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a66*I64a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_569 = Coupling(name = 'GC_569',
                  value = '-(ee**2*complex(0,1)*I63a11*I64a33)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a11*I66a33)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a11*I64a33)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_570 = Coupling(name = 'GC_570',
                  value = '-(ee**2*complex(0,1)*I63a22*I64a33)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a22*I66a33)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a22*I64a33)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_571 = Coupling(name = 'GC_571',
                  value = '-(ee**2*complex(0,1)*I63a11*I64a36)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a11*I66a36)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a11*I64a36)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_572 = Coupling(name = 'GC_572',
                  value = '-(ee**2*complex(0,1)*I63a22*I64a36)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a22*I66a36)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a22*I64a36)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_573 = Coupling(name = 'GC_573',
                  value = '-(ee**2*complex(0,1)*I63a11*I64a63)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a11*I66a63)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a11*I64a63)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_574 = Coupling(name = 'GC_574',
                  value = '-(ee**2*complex(0,1)*I63a22*I64a63)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a22*I66a63)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a22*I64a63)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_575 = Coupling(name = 'GC_575',
                  value = '-(ee**2*complex(0,1)*I63a11*I64a66)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a11*I66a66)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a11*I64a66)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_576 = Coupling(name = 'GC_576',
                  value = '-(ee**2*complex(0,1)*I63a22*I64a66)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a22*I66a66)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a22*I64a66)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_577 = Coupling(name = 'GC_577',
                  value = '(ee**2*complex(0,1)*I71a11*I72a11)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a11*I74a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_578 = Coupling(name = 'GC_578',
                  value = '(ee**2*complex(0,1)*I71a22*I72a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I71a11*I72a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a22*I74a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a11*I74a22)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_579 = Coupling(name = 'GC_579',
                  value = '(ee**2*complex(0,1)*I71a22*I72a22)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a22*I74a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_580 = Coupling(name = 'GC_580',
                  value = '-(ee**2*complex(0,1)*I75a11*I76a33)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a11*I77a33)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a11*I79a33)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a11*I80a33)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I71a33*I72a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I71a11*I72a33)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a33*I74a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a11*I74a33)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_581 = Coupling(name = 'GC_581',
                  value = '-(ee**2*complex(0,1)*I75a22*I76a33)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a22*I77a33)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a22*I79a33)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a22*I80a33)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I71a33*I72a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I71a22*I72a33)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a33*I74a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a22*I74a33)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_582 = Coupling(name = 'GC_582',
                  value = '-(ee**2*complex(0,1)*I75a11*I76a36)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a11*I77a36)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a11*I79a36)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a11*I80a36)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I71a36*I72a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I71a11*I72a36)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a36*I74a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a11*I74a36)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_583 = Coupling(name = 'GC_583',
                  value = '-(ee**2*complex(0,1)*I75a22*I76a36)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a22*I77a36)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a22*I79a36)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a22*I80a36)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I71a36*I72a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I71a22*I72a36)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a36*I74a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a22*I74a36)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_584 = Coupling(name = 'GC_584',
                  value = '-(ee**2*complex(0,1)*I75a11*I76a63)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a11*I77a63)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a11*I79a63)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a11*I80a63)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I71a63*I72a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I71a11*I72a63)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a63*I74a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a11*I74a63)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_585 = Coupling(name = 'GC_585',
                  value = '-(ee**2*complex(0,1)*I75a22*I76a63)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a22*I77a63)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a22*I79a63)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a22*I80a63)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I71a63*I72a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I71a22*I72a63)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a63*I74a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a22*I74a63)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_586 = Coupling(name = 'GC_586',
                  value = '-(ee**2*complex(0,1)*I75a11*I76a66)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a11*I77a66)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a11*I79a66)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a11*I80a66)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I71a66*I72a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I71a11*I72a66)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a66*I74a11)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a11*I74a66)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_587 = Coupling(name = 'GC_587',
                  value = '-(ee**2*complex(0,1)*I75a22*I76a66)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a22*I77a66)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a22*I79a66)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a22*I80a66)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I71a66*I72a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I71a22*I72a66)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a66*I74a22)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a22*I74a66)/(16.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_588 = Coupling(name = 'GC_588',
                  value = '(ee**2*complex(0,1)*I97a11)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_589 = Coupling(name = 'GC_589',
                  value = '(ee**2*complex(0,1)*I161a11*I97a11)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a11*I97a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_590 = Coupling(name = 'GC_590',
                  value = '(ee**2*complex(0,1)*I161a22*I97a11)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a22*I97a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_591 = Coupling(name = 'GC_591',
                  value = '(ee**2*complex(0,1)*I161a33*I97a11)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a33*I97a11)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a33*I97a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_592 = Coupling(name = 'GC_592',
                  value = '(ee**2*complex(0,1)*I161a36*I97a11)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a36*I97a11)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a36*I97a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_593 = Coupling(name = 'GC_593',
                  value = '(ee**2*complex(0,1)*I161a63*I97a11)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a63*I97a11)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a63*I97a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_594 = Coupling(name = 'GC_594',
                  value = '(ee**2*complex(0,1)*I161a66*I97a11)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a66*I97a11)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a66*I97a11)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_595 = Coupling(name = 'GC_595',
                  value = '(ee**2*complex(0,1)*I97a22)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_596 = Coupling(name = 'GC_596',
                  value = '(ee**2*complex(0,1)*I161a11*I97a22)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a11*I97a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_597 = Coupling(name = 'GC_597',
                  value = '(ee**2*complex(0,1)*I161a22*I97a22)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a22*I97a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_598 = Coupling(name = 'GC_598',
                  value = '(ee**2*complex(0,1)*I161a33*I97a22)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a33*I97a22)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a33*I97a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_599 = Coupling(name = 'GC_599',
                  value = '(ee**2*complex(0,1)*I161a36*I97a22)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a36*I97a22)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a36*I97a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_600 = Coupling(name = 'GC_600',
                  value = '(ee**2*complex(0,1)*I161a63*I97a22)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a63*I97a22)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a63*I97a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_601 = Coupling(name = 'GC_601',
                  value = '(ee**2*complex(0,1)*I161a66*I97a22)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a66*I97a22)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a66*I97a22)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_602 = Coupling(name = 'GC_602',
                  value = '(ee**2*complex(0,1)*I97a33)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I98a33)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97a33)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_603 = Coupling(name = 'GC_603',
                  value = '(ee**2*complex(0,1)*I161a11*I97a33)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a11*I98a33)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a11*I97a33)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_604 = Coupling(name = 'GC_604',
                  value = '(ee**2*complex(0,1)*I161a22*I97a33)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a22*I98a33)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a22*I97a33)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_605 = Coupling(name = 'GC_605',
                  value = '(ee**2*complex(0,1)*I161a33*I97a33)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a33*I97a33)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a33*I98a33)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a33*I98a33)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a33*I97a33)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_606 = Coupling(name = 'GC_606',
                  value = '(ee**2*complex(0,1)*I161a36*I97a33)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a36*I97a33)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a36*I98a33)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a36*I98a33)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a36*I97a33)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_607 = Coupling(name = 'GC_607',
                  value = '(ee**2*complex(0,1)*I161a63*I97a33)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a63*I97a33)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a63*I98a33)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a63*I98a33)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a63*I97a33)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_608 = Coupling(name = 'GC_608',
                  value = '(ee**2*complex(0,1)*I161a66*I97a33)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a66*I97a33)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a66*I98a33)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a66*I98a33)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a66*I97a33)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_609 = Coupling(name = 'GC_609',
                  value = '(ee**2*complex(0,1)*I97a36)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I98a36)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97a36)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_610 = Coupling(name = 'GC_610',
                  value = '(ee**2*complex(0,1)*I161a11*I97a36)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a11*I98a36)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a11*I97a36)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_611 = Coupling(name = 'GC_611',
                  value = '(ee**2*complex(0,1)*I161a22*I97a36)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a22*I98a36)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a22*I97a36)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_612 = Coupling(name = 'GC_612',
                  value = '(ee**2*complex(0,1)*I161a33*I97a36)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a33*I97a36)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a33*I98a36)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a33*I98a36)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a33*I97a36)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_613 = Coupling(name = 'GC_613',
                  value = '(ee**2*complex(0,1)*I161a36*I97a36)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a36*I97a36)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a36*I98a36)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a36*I98a36)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a36*I97a36)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_614 = Coupling(name = 'GC_614',
                  value = '(ee**2*complex(0,1)*I161a63*I97a36)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a63*I97a36)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a63*I98a36)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a63*I98a36)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a63*I97a36)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_615 = Coupling(name = 'GC_615',
                  value = '(ee**2*complex(0,1)*I161a66*I97a36)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a66*I97a36)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a66*I98a36)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a66*I98a36)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a66*I97a36)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_616 = Coupling(name = 'GC_616',
                  value = '(ee**2*complex(0,1)*I97a63)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I98a63)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97a63)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_617 = Coupling(name = 'GC_617',
                  value = '(ee**2*complex(0,1)*I161a11*I97a63)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a11*I98a63)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a11*I97a63)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_618 = Coupling(name = 'GC_618',
                  value = '(ee**2*complex(0,1)*I161a22*I97a63)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a22*I98a63)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a22*I97a63)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_619 = Coupling(name = 'GC_619',
                  value = '(ee**2*complex(0,1)*I161a33*I97a63)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a33*I97a63)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a33*I98a63)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a33*I98a63)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a33*I97a63)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_620 = Coupling(name = 'GC_620',
                  value = '(ee**2*complex(0,1)*I161a36*I97a63)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a36*I97a63)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a36*I98a63)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a36*I98a63)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a36*I97a63)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_621 = Coupling(name = 'GC_621',
                  value = '(ee**2*complex(0,1)*I161a63*I97a63)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a63*I97a63)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a63*I98a63)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a63*I98a63)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a63*I97a63)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_622 = Coupling(name = 'GC_622',
                  value = '(ee**2*complex(0,1)*I161a66*I97a63)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a66*I97a63)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a66*I98a63)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a66*I98a63)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a66*I97a63)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_623 = Coupling(name = 'GC_623',
                  value = '(ee**2*complex(0,1)*I97a66)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I98a66)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97a66)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_624 = Coupling(name = 'GC_624',
                  value = '(ee**2*complex(0,1)*I161a11*I97a66)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a11*I98a66)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a11*I97a66)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_625 = Coupling(name = 'GC_625',
                  value = '(ee**2*complex(0,1)*I161a22*I97a66)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a22*I98a66)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a22*I97a66)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_626 = Coupling(name = 'GC_626',
                  value = '(ee**2*complex(0,1)*I161a33*I97a66)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a33*I97a66)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a33*I98a66)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a33*I98a66)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a33*I97a66)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_627 = Coupling(name = 'GC_627',
                  value = '(ee**2*complex(0,1)*I161a36*I97a66)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a36*I97a66)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a36*I98a66)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a36*I98a66)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a36*I97a66)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_628 = Coupling(name = 'GC_628',
                  value = '(ee**2*complex(0,1)*I161a63*I97a66)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a63*I97a66)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a63*I98a66)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a63*I98a66)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a63*I97a66)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_629 = Coupling(name = 'GC_629',
                  value = '(ee**2*complex(0,1)*I161a66*I97a66)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I162a66*I97a66)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I161a66*I98a66)/(6.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I162a66*I98a66)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I161a66*I97a66)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_630 = Coupling(name = 'GC_630',
                  value = '(ee**2*complex(0,1)*I97a11)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I100a11*I99a11)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97a11)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I100a11*I99a11)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_631 = Coupling(name = 'GC_631',
                  value = '-(ee**2*complex(0,1)*I100a22*I99a11)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I100a22*I99a11)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_632 = Coupling(name = 'GC_632',
                  value = '-(ee**2*complex(0,1)*I100a33*I99a11)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I100a33*I99a11)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_633 = Coupling(name = 'GC_633',
                  value = '-(ee**2*complex(0,1)*I100a36*I99a11)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I100a36*I99a11)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_634 = Coupling(name = 'GC_634',
                  value = '-(ee**2*complex(0,1)*I100a11*I99a22)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I100a11*I99a22)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_635 = Coupling(name = 'GC_635',
                  value = '(ee**2*complex(0,1)*I97a22)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I100a22*I99a22)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97a22)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I100a22*I99a22)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_636 = Coupling(name = 'GC_636',
                  value = '-(ee**2*complex(0,1)*I100a33*I99a22)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I100a33*I99a22)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_637 = Coupling(name = 'GC_637',
                  value = '-(ee**2*complex(0,1)*I100a36*I99a22)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I100a36*I99a22)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_638 = Coupling(name = 'GC_638',
                  value = '-(ee**2*complex(0,1)*I100a11*I99a33)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I100a11*I99a33)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_639 = Coupling(name = 'GC_639',
                  value = '-(ee**2*complex(0,1)*I100a22*I99a33)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I100a22*I99a33)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_640 = Coupling(name = 'GC_640',
                  value = '-(ee**2*complex(0,1)*I100a11*I99a36)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I100a11*I99a36)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_641 = Coupling(name = 'GC_641',
                  value = '-(ee**2*complex(0,1)*I100a22*I99a36)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I100a22*I99a36)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_642 = Coupling(name = 'GC_642',
                  value = '(-2*cw*ee**2*complex(0,1)*I113a11)/(3.*(-1 + sw)*sw*(1 + sw)) + (8*cw*ee**2*complex(0,1)*I113a11*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_643 = Coupling(name = 'GC_643',
                  value = '-((cw*ee*complex(0,1)*G*I113a11)/((-1 + sw)*sw*(1 + sw))) + (4*cw*ee*complex(0,1)*G*I113a11*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_644 = Coupling(name = 'GC_644',
                  value = '(-2*cw*ee**2*complex(0,1)*I113a22)/(3.*(-1 + sw)*sw*(1 + sw)) + (8*cw*ee**2*complex(0,1)*I113a22*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_645 = Coupling(name = 'GC_645',
                  value = '-((cw*ee*complex(0,1)*G*I113a22)/((-1 + sw)*sw*(1 + sw))) + (4*cw*ee*complex(0,1)*G*I113a22*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_646 = Coupling(name = 'GC_646',
                  value = '(-2*cw*ee**2*complex(0,1)*I113a33)/(3.*(-1 + sw)*sw*(1 + sw)) + (8*cw*ee**2*complex(0,1)*I113a33*sw)/(9.*(-1 + sw)*(1 + sw)) + (8*cw*ee**2*complex(0,1)*I114a33*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_647 = Coupling(name = 'GC_647',
                  value = '-((cw*ee*complex(0,1)*G*I113a33)/((-1 + sw)*sw*(1 + sw))) + (4*cw*ee*complex(0,1)*G*I113a33*sw)/(3.*(-1 + sw)*(1 + sw)) + (4*cw*ee*complex(0,1)*G*I114a33*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_648 = Coupling(name = 'GC_648',
                  value = '(-2*cw*ee**2*complex(0,1)*I113a36)/(3.*(-1 + sw)*sw*(1 + sw)) + (8*cw*ee**2*complex(0,1)*I113a36*sw)/(9.*(-1 + sw)*(1 + sw)) + (8*cw*ee**2*complex(0,1)*I114a36*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_649 = Coupling(name = 'GC_649',
                  value = '-((cw*ee*complex(0,1)*G*I113a36)/((-1 + sw)*sw*(1 + sw))) + (4*cw*ee*complex(0,1)*G*I113a36*sw)/(3.*(-1 + sw)*(1 + sw)) + (4*cw*ee*complex(0,1)*G*I114a36*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_650 = Coupling(name = 'GC_650',
                  value = '(-2*cw*ee**2*complex(0,1)*I113a63)/(3.*(-1 + sw)*sw*(1 + sw)) + (8*cw*ee**2*complex(0,1)*I113a63*sw)/(9.*(-1 + sw)*(1 + sw)) + (8*cw*ee**2*complex(0,1)*I114a63*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_651 = Coupling(name = 'GC_651',
                  value = '-((cw*ee*complex(0,1)*G*I113a63)/((-1 + sw)*sw*(1 + sw))) + (4*cw*ee*complex(0,1)*G*I113a63*sw)/(3.*(-1 + sw)*(1 + sw)) + (4*cw*ee*complex(0,1)*G*I114a63*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_652 = Coupling(name = 'GC_652',
                  value = '(-2*cw*ee**2*complex(0,1)*I113a66)/(3.*(-1 + sw)*sw*(1 + sw)) + (8*cw*ee**2*complex(0,1)*I113a66*sw)/(9.*(-1 + sw)*(1 + sw)) + (8*cw*ee**2*complex(0,1)*I114a66*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_653 = Coupling(name = 'GC_653',
                  value = '-((cw*ee*complex(0,1)*G*I113a66)/((-1 + sw)*sw*(1 + sw))) + (4*cw*ee*complex(0,1)*G*I113a66*sw)/(3.*(-1 + sw)*(1 + sw)) + (4*cw*ee*complex(0,1)*G*I114a66*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_654 = Coupling(name = 'GC_654',
                  value = '(cw*ee*complex(0,1)*I201a11)/(2.*(-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*I201a11*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_655 = Coupling(name = 'GC_655',
                  value = '(cw*ee*complex(0,1)*I201a22)/(2.*(-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*I201a22*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_656 = Coupling(name = 'GC_656',
                  value = '(cw*ee*complex(0,1)*I201a33)/(2.*(-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*I201a33*sw)/(3.*(-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*I203a33*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_657 = Coupling(name = 'GC_657',
                  value = '(cw*ee*complex(0,1)*I201a36)/(2.*(-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*I201a36*sw)/(3.*(-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*I203a36*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_658 = Coupling(name = 'GC_658',
                  value = '-(cw*ee*complex(0,1)*I201a63)/(2.*(-1 + sw)*sw*(1 + sw)) + (2*cw*ee*complex(0,1)*I201a63*sw)/(3.*(-1 + sw)*(1 + sw)) + (2*cw*ee*complex(0,1)*I203a63*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_659 = Coupling(name = 'GC_659',
                  value = '(cw*ee*complex(0,1)*I201a66)/(2.*(-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*I201a66*sw)/(3.*(-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*I203a66*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_660 = Coupling(name = 'GC_660',
                  value = '-((cw*ee**2*complex(0,1)*I47a11)/((-1 + sw)*sw*(1 + sw))) + (2*cw*ee**2*complex(0,1)*I47a11*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_661 = Coupling(name = 'GC_661',
                  value = '-((cw*ee**2*complex(0,1)*I47a22)/((-1 + sw)*sw*(1 + sw))) + (2*cw*ee**2*complex(0,1)*I47a22*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_662 = Coupling(name = 'GC_662',
                  value = '-((cw*ee**2*complex(0,1)*I47a33)/((-1 + sw)*sw*(1 + sw))) + (2*cw*ee**2*complex(0,1)*I47a33*sw)/((-1 + sw)*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I48a33*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_663 = Coupling(name = 'GC_663',
                  value = '-((cw*ee**2*complex(0,1)*I47a36)/((-1 + sw)*sw*(1 + sw))) + (2*cw*ee**2*complex(0,1)*I47a36*sw)/((-1 + sw)*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I48a36*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_664 = Coupling(name = 'GC_664',
                  value = '-((cw*ee**2*complex(0,1)*I47a63)/((-1 + sw)*sw*(1 + sw))) + (2*cw*ee**2*complex(0,1)*I47a63*sw)/((-1 + sw)*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I48a63*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_665 = Coupling(name = 'GC_665',
                  value = '-((cw*ee**2*complex(0,1)*I47a66)/((-1 + sw)*sw*(1 + sw))) + (2*cw*ee**2*complex(0,1)*I47a66*sw)/((-1 + sw)*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I48a66*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_666 = Coupling(name = 'GC_666',
                  value = '-(cw*ee*complex(0,1)*I63a11)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*I63a11*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_667 = Coupling(name = 'GC_667',
                  value = '-(cw*ee*complex(0,1)*I63a22)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*I63a22*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_668 = Coupling(name = 'GC_668',
                  value = '-(cw*ee*complex(0,1)*I63a33)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*I63a33*sw)/(3.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*I65a33*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_669 = Coupling(name = 'GC_669',
                  value = '-(cw*ee*complex(0,1)*I63a36)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*I63a36*sw)/(3.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*I65a36*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_670 = Coupling(name = 'GC_670',
                  value = '(cw*ee*complex(0,1)*I63a63)/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*I63a63*sw)/(3.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*I65a63*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_671 = Coupling(name = 'GC_671',
                  value = '-(cw*ee*complex(0,1)*I63a66)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*I63a66*sw)/(3.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*I65a66*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_672 = Coupling(name = 'GC_672',
                  value = '-(cw*ee**2*complex(0,1)*I8a11)/(3.*(-1 + sw)*sw*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I8a11*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_673 = Coupling(name = 'GC_673',
                  value = '(cw*ee*complex(0,1)*G*I8a11)/((-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*G*I8a11*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_674 = Coupling(name = 'GC_674',
                  value = '-(cw*ee**2*complex(0,1)*I8a22)/(3.*(-1 + sw)*sw*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I8a22*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_675 = Coupling(name = 'GC_675',
                  value = '(cw*ee*complex(0,1)*G*I8a22)/((-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*G*I8a22*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_676 = Coupling(name = 'GC_676',
                  value = '-(cw*ee*complex(0,1)*I97a11)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*I97a11*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_677 = Coupling(name = 'GC_677',
                  value = '-(cw*ee*complex(0,1)*I97a22)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*I97a22*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_678 = Coupling(name = 'GC_678',
                  value = '-(cw*ee*complex(0,1)*I97a33)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*I97a33*sw)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*I98a33*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_679 = Coupling(name = 'GC_679',
                  value = '-(cw*ee*complex(0,1)*I97a36)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*I97a36*sw)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*I98a36*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_680 = Coupling(name = 'GC_680',
                  value = '(cw*ee*complex(0,1)*I97a63)/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*I97a63*sw)/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*I98a63*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_681 = Coupling(name = 'GC_681',
                  value = '-(cw*ee*complex(0,1)*I97a66)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*I97a66*sw)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*I98a66*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_682 = Coupling(name = 'GC_682',
                  value = '-(cw*ee**2*complex(0,1)*I8a33)/(3.*(-1 + sw)*sw*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I8a33*sw)/(9.*(-1 + sw)*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I9a33*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_683 = Coupling(name = 'GC_683',
                  value = '(cw*ee*complex(0,1)*G*I8a33)/((-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*G*I8a33*sw)/(3.*(-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*G*I9a33*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_684 = Coupling(name = 'GC_684',
                  value = '-(cw*ee**2*complex(0,1)*I8a36)/(3.*(-1 + sw)*sw*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I8a36*sw)/(9.*(-1 + sw)*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I9a36*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_685 = Coupling(name = 'GC_685',
                  value = '(cw*ee*complex(0,1)*G*I8a36)/((-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*G*I8a36*sw)/(3.*(-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*G*I9a36*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_686 = Coupling(name = 'GC_686',
                  value = '-(cw*ee**2*complex(0,1)*I8a63)/(3.*(-1 + sw)*sw*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I8a63*sw)/(9.*(-1 + sw)*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I9a63*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_687 = Coupling(name = 'GC_687',
                  value = '(cw*ee*complex(0,1)*G*I8a63)/((-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*G*I8a63*sw)/(3.*(-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*G*I9a63*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_688 = Coupling(name = 'GC_688',
                  value = '-(cw*ee**2*complex(0,1)*I8a66)/(3.*(-1 + sw)*sw*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I8a66*sw)/(9.*(-1 + sw)*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I9a66*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_689 = Coupling(name = 'GC_689',
                  value = '(cw*ee*complex(0,1)*G*I8a66)/((-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*G*I8a66*sw)/(3.*(-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*G*I9a66*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_690 = Coupling(name = 'GC_690',
                  value = '(complex(0,1)*I101a33*I102a33)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I97a33)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I98a33)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I100a33*I99a33)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97a33)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I100a33*I99a33)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I101a33*I102a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_691 = Coupling(name = 'GC_691',
                  value = '(complex(0,1)*I101a36*I102a33)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I97a63)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I98a63)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I100a36*I99a33)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97a63)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I100a36*I99a33)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I101a36*I102a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_692 = Coupling(name = 'GC_692',
                  value = '(complex(0,1)*I101a33*I102a36)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I97a36)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I98a36)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I100a33*I99a36)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97a36)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I100a33*I99a36)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I101a33*I102a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_693 = Coupling(name = 'GC_693',
                  value = '(complex(0,1)*I101a36*I102a36)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I97a66)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I98a66)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I100a36*I99a36)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97a66)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I100a36*I99a36)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I101a36*I102a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_694 = Coupling(name = 'GC_694',
                  value = '-(ee**2*complex(0,1)*I159a33*I163a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I164a33*I165a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I166a33*I167a33)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a33*I163a33)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I164a33*I165a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I166a33*I167a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_695 = Coupling(name = 'GC_695',
                  value = '-(ee**2*complex(0,1)*I159a33*I163a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I164a33*I165a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I166a36*I167a33)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a33*I163a36)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I164a33*I165a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I166a36*I167a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_696 = Coupling(name = 'GC_696',
                  value = '-(ee**2*complex(0,1)*I159a33*I163a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I164a33*I165a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I166a63*I167a33)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a33*I163a63)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I164a33*I165a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I166a63*I167a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_697 = Coupling(name = 'GC_697',
                  value = '-(ee**2*complex(0,1)*I159a33*I163a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I164a33*I165a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I166a66*I167a33)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a33*I163a66)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I164a33*I165a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I166a66*I167a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_698 = Coupling(name = 'GC_698',
                  value = '-(ee**2*complex(0,1)*I159a36*I163a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I164a36*I165a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I166a33*I167a36)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a36*I163a33)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I164a36*I165a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I166a33*I167a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_699 = Coupling(name = 'GC_699',
                  value = '-(ee**2*complex(0,1)*I159a36*I163a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I164a36*I165a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I166a36*I167a36)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a36*I163a36)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I164a36*I165a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I166a36*I167a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_700 = Coupling(name = 'GC_700',
                  value = '-(ee**2*complex(0,1)*I159a36*I163a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I164a36*I165a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I166a63*I167a36)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a36*I163a63)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I164a36*I165a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I166a63*I167a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_701 = Coupling(name = 'GC_701',
                  value = '-(ee**2*complex(0,1)*I159a36*I163a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I164a36*I165a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I166a66*I167a36)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a36*I163a66)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I164a36*I165a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I166a66*I167a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_702 = Coupling(name = 'GC_702',
                  value = '-(ee**2*complex(0,1)*I159a63*I163a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I164a63*I165a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I166a33*I167a63)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a63*I163a33)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I164a63*I165a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I166a33*I167a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_703 = Coupling(name = 'GC_703',
                  value = '-(ee**2*complex(0,1)*I159a63*I163a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I164a63*I165a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I166a36*I167a63)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a63*I163a36)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I164a63*I165a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I166a36*I167a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_704 = Coupling(name = 'GC_704',
                  value = '-(ee**2*complex(0,1)*I159a63*I163a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I164a63*I165a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I166a63*I167a63)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a63*I163a63)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I164a63*I165a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I166a63*I167a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_705 = Coupling(name = 'GC_705',
                  value = '-(ee**2*complex(0,1)*I159a63*I163a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I164a63*I165a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I166a66*I167a63)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a63*I163a66)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I164a63*I165a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I166a66*I167a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_706 = Coupling(name = 'GC_706',
                  value = '-(ee**2*complex(0,1)*I159a66*I163a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I164a66*I165a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I166a33*I167a66)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a66*I163a33)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I164a66*I165a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I166a33*I167a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_707 = Coupling(name = 'GC_707',
                  value = '-(ee**2*complex(0,1)*I159a66*I163a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I164a66*I165a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I166a36*I167a66)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a66*I163a36)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I164a66*I165a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I166a36*I167a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_708 = Coupling(name = 'GC_708',
                  value = '-(ee**2*complex(0,1)*I159a66*I163a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I164a66*I165a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I166a63*I167a66)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a66*I163a63)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I164a66*I165a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I166a63*I167a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_709 = Coupling(name = 'GC_709',
                  value = '-(ee**2*complex(0,1)*I159a66*I163a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I164a66*I165a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I166a66*I167a66)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I159a66*I163a66)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I164a66*I165a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I166a66*I167a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_710 = Coupling(name = 'GC_710',
                  value = '(complex(0,1)*G**2*I168a11*I169a11)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a11*I169a11*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_711 = Coupling(name = 'GC_711',
                  value = '(complex(0,1)*G**2*I168a22*I169a22)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a22*I169a22*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_712 = Coupling(name = 'GC_712',
                  value = '(complex(0,1)*G**2*I170a11*I171a11)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a11*I171a11*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_713 = Coupling(name = 'GC_713',
                  value = '(complex(0,1)*G**2*I168a22*I169a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I168a11*I169a22)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a22*I171a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a11*I171a22)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a22*I169a11*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a11*I169a22*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a22*I171a11*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a11*I171a22*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_714 = Coupling(name = 'GC_714',
                  value = '(complex(0,1)*G**2*I170a22*I171a22)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a22*I171a22*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_715 = Coupling(name = 'GC_715',
                  value = '(complex(0,1)*G**2*I172a44*I175a44)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a44*I175a44*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_716 = Coupling(name = 'GC_716',
                  value = '(complex(0,1)*G**2*I172a55*I175a55)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a55*I175a55*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_717 = Coupling(name = 'GC_717',
                  value = '(complex(0,1)*G**2*I168a33*I169a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a33*I172a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a33*I175a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a33*I176a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a33*I169a33*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a33*I172a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a33*I175a33*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a33*I176a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_718 = Coupling(name = 'GC_718',
                  value = '(complex(0,1)*G**2*I168a36*I169a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a36*I172a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a36*I175a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a36*I176a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a36*I169a36*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a36*I172a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a36*I175a36*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a36*I176a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_719 = Coupling(name = 'GC_719',
                  value = '(complex(0,1)*G**2*I168a63*I169a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a63*I172a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a63*I175a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a63*I176a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a63*I169a63*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a63*I172a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a63*I175a63*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a63*I176a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_720 = Coupling(name = 'GC_720',
                  value = '(complex(0,1)*G**2*I168a66*I169a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a66*I172a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a66*I175a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a66*I176a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a66*I169a66*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a66*I172a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a66*I175a66*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a66*I176a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_721 = Coupling(name = 'GC_721',
                  value = '(complex(0,1)*G**2*I168a33*I169a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I168a11*I169a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a33*I171a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a11*I171a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a11*I172a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a11*I174a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a11*I176a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a11*I177a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a33*I169a11*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a11*I169a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a33*I171a11*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a11*I171a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a11*I172a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a11*I174a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a11*I176a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a11*I177a33*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_722 = Coupling(name = 'GC_722',
                  value = '(complex(0,1)*G**2*I168a33*I169a22)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I168a22*I169a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a33*I171a22)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a22*I171a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a22*I172a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a22*I174a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a22*I176a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a22*I177a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a33*I169a22*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a22*I169a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a33*I171a22*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a22*I171a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a22*I172a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a22*I174a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a22*I176a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a22*I177a33*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_723 = Coupling(name = 'GC_723',
                  value = '(complex(0,1)*G**2*I168a36*I169a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I168a11*I169a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a36*I171a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a11*I171a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a11*I172a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a11*I174a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a11*I176a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a11*I177a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a36*I169a11*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a11*I169a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a36*I171a11*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a11*I171a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a11*I172a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a11*I174a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a11*I176a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a11*I177a36*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_724 = Coupling(name = 'GC_724',
                  value = '(complex(0,1)*G**2*I168a36*I169a22)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I168a22*I169a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a36*I171a22)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a22*I171a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a22*I172a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a22*I174a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a22*I176a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a22*I177a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a36*I169a22*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a22*I169a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a36*I171a22*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a22*I171a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a22*I172a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a22*I174a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a22*I176a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a22*I177a36*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_725 = Coupling(name = 'GC_725',
                  value = '-(complex(0,1)*G**2*I169a11*I172a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a11*I174a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a11*I176a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a11*I177a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a11*I172a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a11*I174a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a11*I176a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a11*I177a44*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_726 = Coupling(name = 'GC_726',
                  value = '-(complex(0,1)*G**2*I169a22*I172a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a22*I174a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a22*I176a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a22*I177a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a22*I172a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a22*I174a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a22*I176a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a22*I177a44*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_727 = Coupling(name = 'GC_727',
                  value = '-(complex(0,1)*G**2*I169a11*I172a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a11*I174a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a11*I176a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a11*I177a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a11*I172a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a11*I174a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a11*I176a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a11*I177a55*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_728 = Coupling(name = 'GC_728',
                  value = '-(complex(0,1)*G**2*I169a22*I172a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a22*I174a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a22*I176a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a22*I177a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a22*I172a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a22*I174a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a22*I176a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a22*I177a55*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_729 = Coupling(name = 'GC_729',
                  value = '(complex(0,1)*G**2*I168a63*I169a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I168a11*I169a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a63*I171a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a11*I171a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a11*I172a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a11*I174a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a11*I176a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a11*I177a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a63*I169a11*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a11*I169a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a63*I171a11*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a11*I171a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a11*I172a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a11*I174a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a11*I176a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a11*I177a63*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_730 = Coupling(name = 'GC_730',
                  value = '(complex(0,1)*G**2*I168a63*I169a22)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I168a22*I169a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a63*I171a22)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a22*I171a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a22*I172a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a22*I174a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a22*I176a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a22*I177a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a63*I169a22*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a22*I169a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a63*I171a22*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a22*I171a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a22*I172a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a22*I174a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a22*I176a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a22*I177a63*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_731 = Coupling(name = 'GC_731',
                  value = '(complex(0,1)*G**2*I168a66*I169a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I168a11*I169a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a66*I171a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a11*I171a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a11*I172a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a11*I174a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a11*I176a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a11*I177a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a66*I169a11*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a11*I169a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a66*I171a11*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a11*I171a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a11*I172a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a11*I174a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a11*I176a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a11*I177a66*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_732 = Coupling(name = 'GC_732',
                  value = '(complex(0,1)*G**2*I168a66*I169a22)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I168a22*I169a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a66*I171a22)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a22*I171a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a22*I172a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a22*I174a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a22*I176a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a22*I177a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a66*I169a22*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a22*I169a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a66*I171a22*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a22*I171a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a22*I172a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a22*I174a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a22*I176a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a22*I177a66*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_733 = Coupling(name = 'GC_733',
                  value = '(complex(0,1)*G**2*I170a33*I171a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a33*I174a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a33*I177a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a33*I178a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a33*I171a33*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a33*I174a33*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a33*I177a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a33*I178a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_734 = Coupling(name = 'GC_734',
                  value = '(complex(0,1)*G**2*I168a36*I169a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a33*I171a36)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a33*I172a36)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a36*I174a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a36*I175a33)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a36*I176a33)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a33*I177a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a36*I178a33)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a36*I169a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a33*I171a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a33*I172a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a36*I174a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a36*I175a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a36*I176a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a33*I177a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a36*I178a33*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_735 = Coupling(name = 'GC_735',
                  value = '(complex(0,1)*G**2*I168a33*I169a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a33*I171a63)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a63*I172a33)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a63*I174a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a33*I175a63)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a33*I176a63)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a33*I177a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a63*I178a33)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a33*I169a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a33*I171a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a63*I172a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a63*I174a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a33*I175a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a33*I176a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a33*I177a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a63*I178a33*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_736 = Coupling(name = 'GC_736',
                  value = '(complex(0,1)*G**2*I168a33*I169a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a36*I171a33)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a36*I172a33)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a33*I174a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a33*I175a36)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a33*I176a36)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a36*I177a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a33*I178a36)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a33*I169a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a36*I171a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a36*I172a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a33*I174a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a33*I175a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a33*I176a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a36*I177a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a33*I178a36*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_737 = Coupling(name = 'GC_737',
                  value = '(complex(0,1)*G**2*I170a36*I171a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a36*I174a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a36*I177a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a36*I178a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a36*I171a36*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a36*I174a36*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a36*I177a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a36*I178a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_738 = Coupling(name = 'GC_738',
                  value = '(complex(0,1)*G**2*I168a36*I169a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a36*I171a66)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a66*I172a36)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a66*I174a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a36*I175a66)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a36*I176a66)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a36*I177a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a66*I178a36)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a36*I169a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a36*I171a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a66*I172a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a66*I174a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a36*I175a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a36*I176a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a36*I177a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a66*I178a36*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_739 = Coupling(name = 'GC_739',
                  value = '-(complex(0,1)*G**2*I169a33*I172a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a33*I174a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a44*I175a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a33*I175a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a33*I176a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a33*I177a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a44*I178a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a33*I178a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a33*I172a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a33*I174a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a44*I175a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a33*I175a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a33*I176a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a33*I177a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a44*I178a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a33*I178a44*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_740 = Coupling(name = 'GC_740',
                  value = '-(complex(0,1)*G**2*I169a36*I172a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a36*I174a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a44*I175a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a36*I175a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a36*I176a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a36*I177a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a44*I178a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a36*I178a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a36*I172a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a36*I174a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a44*I175a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a36*I175a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a36*I176a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a36*I177a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a44*I178a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a36*I178a44*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_741 = Coupling(name = 'GC_741',
                  value = '(complex(0,1)*G**2*I177a44*I178a44)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a44*I178a44*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_742 = Coupling(name = 'GC_742',
                  value = '-(complex(0,1)*G**2*I169a33*I172a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a33*I174a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a55*I175a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a33*I175a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a33*I176a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a33*I177a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a55*I178a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a33*I178a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a33*I172a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a33*I174a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a55*I175a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a33*I175a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a33*I176a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a33*I177a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a55*I178a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a33*I178a55*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_743 = Coupling(name = 'GC_743',
                  value = '-(complex(0,1)*G**2*I169a36*I172a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a36*I174a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a55*I175a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a36*I175a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a36*I176a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a36*I177a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a55*I178a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a36*I178a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a36*I172a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a36*I174a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a55*I175a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a36*I175a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a36*I176a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a36*I177a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a55*I178a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a36*I178a55*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_744 = Coupling(name = 'GC_744',
                  value = '(complex(0,1)*G**2*I172a55*I175a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a44*I175a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a55*I178a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a44*I178a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a55*I175a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a44*I175a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a55*I178a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a44*I178a55*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_745 = Coupling(name = 'GC_745',
                  value = '(complex(0,1)*G**2*I177a55*I178a55)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a55*I178a55*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_746 = Coupling(name = 'GC_746',
                  value = '(complex(0,1)*G**2*I168a63*I169a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a63*I171a33)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a33*I172a63)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a33*I174a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a63*I175a33)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a63*I176a33)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a63*I177a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a33*I178a63)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a63*I169a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a63*I171a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a33*I172a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a33*I174a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a63*I175a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a63*I176a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a63*I177a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a33*I178a63*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_747 = Coupling(name = 'GC_747',
                  value = '(complex(0,1)*G**2*I168a63*I169a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I168a36*I169a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a63*I171a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a36*I171a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a63*I172a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a36*I172a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a63*I174a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a36*I174a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a63*I175a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a36*I175a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a63*I176a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a36*I176a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a63*I177a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a36*I177a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a63*I178a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a36*I178a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a63*I169a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a36*I169a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a63*I171a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a36*I171a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a63*I172a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a36*I172a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a63*I174a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a36*I174a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a63*I175a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a36*I175a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a63*I176a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a36*I176a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a63*I177a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a36*I177a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a63*I178a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a36*I178a63*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_748 = Coupling(name = 'GC_748',
                  value = '-(complex(0,1)*G**2*I169a63*I172a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a63*I174a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a63*I175a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a44*I175a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a63*I176a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a63*I177a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a63*I178a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a44*I178a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a63*I172a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a63*I174a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a63*I175a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a44*I175a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a63*I176a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a63*I177a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a63*I178a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a44*I178a63*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_749 = Coupling(name = 'GC_749',
                  value = '-(complex(0,1)*G**2*I169a63*I172a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a63*I174a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a63*I175a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a55*I175a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a63*I176a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a63*I177a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a63*I178a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a55*I178a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a63*I172a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a63*I174a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a63*I175a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a55*I175a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a63*I176a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a63*I177a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a63*I178a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a55*I178a63*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_750 = Coupling(name = 'GC_750',
                  value = '(complex(0,1)*G**2*I170a63*I171a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a63*I174a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a63*I177a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a63*I178a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a63*I171a63*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a63*I174a63*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a63*I177a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a63*I178a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_751 = Coupling(name = 'GC_751',
                  value = '(complex(0,1)*G**2*I168a66*I169a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a63*I171a66)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a63*I172a66)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a66*I174a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a66*I175a63)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a66*I176a63)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a63*I177a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a66*I178a63)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a66*I169a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a63*I171a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a63*I172a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a66*I174a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a66*I175a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a66*I176a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a63*I177a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a66*I178a63*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_752 = Coupling(name = 'GC_752',
                  value = '(complex(0,1)*G**2*I168a66*I169a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I168a33*I169a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a66*I171a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a33*I171a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a66*I172a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a33*I172a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a66*I174a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a33*I174a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a66*I175a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a33*I175a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a66*I176a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a33*I176a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a66*I177a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a33*I177a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a66*I178a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a33*I178a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a66*I169a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a33*I169a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a66*I171a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a33*I171a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a66*I172a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a33*I172a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a66*I174a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a33*I174a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a66*I175a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a33*I175a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a66*I176a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a33*I176a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a66*I177a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a33*I177a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a66*I178a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a33*I178a66*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_753 = Coupling(name = 'GC_753',
                  value = '(complex(0,1)*G**2*I168a66*I169a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a66*I171a36)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a36*I172a66)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a36*I174a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a66*I175a36)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a66*I176a36)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a66*I177a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a36*I178a66)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a66*I169a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a66*I171a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a36*I172a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a36*I174a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a66*I175a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a66*I176a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a66*I177a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a36*I178a66*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_754 = Coupling(name = 'GC_754',
                  value = '-(complex(0,1)*G**2*I169a66*I172a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a66*I174a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a66*I175a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a44*I175a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a66*I176a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a66*I177a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a66*I178a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a44*I178a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a66*I172a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a66*I174a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a66*I175a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a44*I175a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a66*I176a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a66*I177a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a66*I178a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a44*I178a66*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_755 = Coupling(name = 'GC_755',
                  value = '-(complex(0,1)*G**2*I169a66*I172a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a66*I174a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a66*I175a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a55*I175a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a66*I176a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a66*I177a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a66*I178a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a55*I178a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a66*I172a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a66*I174a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a66*I175a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a55*I175a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a66*I176a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a66*I177a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a66*I178a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a55*I178a66*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_756 = Coupling(name = 'GC_756',
                  value = '(complex(0,1)*G**2*I168a63*I169a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a66*I171a63)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I169a66*I172a63)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a63*I174a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I172a63*I175a66)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I113a63*I176a66)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a66*I177a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a63*I178a66)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I168a63*I169a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a66*I171a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I169a66*I172a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a63*I174a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I172a63*I175a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I113a63*I176a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a66*I177a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a63*I178a66*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_757 = Coupling(name = 'GC_757',
                  value = '(complex(0,1)*G**2*I170a66*I171a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I173a66*I174a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a66*I177a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I177a66*I178a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I170a66*I171a66*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I173a66*I174a66*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I170a66*I177a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I177a66*I178a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_758 = Coupling(name = 'GC_758',
                  value = '(-2*ee**2*complex(0,1)*I168a33*I169a33)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a33*I172a33)/(9.*(-1 + sw)*(1 + sw)) + (4*ee**2*complex(0,1)*I172a33*I175a33)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a33*I176a33)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I179a33*I180a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I183a33*I184a33)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a33*I169a33)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I179a33*I180a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I183a33*I184a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_759 = Coupling(name = 'GC_759',
                  value = '(-2*ee**2*complex(0,1)*I168a36*I169a36)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a36*I172a36)/(9.*(-1 + sw)*(1 + sw)) + (4*ee**2*complex(0,1)*I172a36*I175a36)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a36*I176a36)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I179a36*I180a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I183a36*I184a36)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a36*I169a36)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I179a36*I180a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I183a36*I184a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_760 = Coupling(name = 'GC_760',
                  value = '(-2*ee**2*complex(0,1)*I168a63*I169a63)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a63*I172a63)/(9.*(-1 + sw)*(1 + sw)) + (4*ee**2*complex(0,1)*I172a63*I175a63)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a63*I176a63)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I179a63*I180a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I183a63*I184a63)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a63*I169a63)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I179a63*I180a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I183a63*I184a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_761 = Coupling(name = 'GC_761',
                  value = '(-2*ee**2*complex(0,1)*I168a66*I169a66)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a66*I172a66)/(9.*(-1 + sw)*(1 + sw)) + (4*ee**2*complex(0,1)*I172a66*I175a66)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a66*I176a66)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I179a66*I180a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I183a66*I184a66)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a66*I169a66)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I179a66*I180a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I183a66*I184a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_762 = Coupling(name = 'GC_762',
                  value = '(-2*ee**2*complex(0,1)*I170a33*I171a33)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a33*I174a33)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a33*I177a33)/(9.*(-1 + sw)*(1 + sw)) + (4*ee**2*complex(0,1)*I177a33*I178a33)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I181a33*I182a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I185a33*I186a33)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I170a33*I171a33)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I181a33*I182a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I185a33*I186a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_763 = Coupling(name = 'GC_763',
                  value = '-(ee**2*complex(0,1)*I168a33*I169a36)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a36*I171a33)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a36*I172a33)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a33*I174a36)/(18.*(-1 + sw)*(1 + sw)) + (2*ee**2*complex(0,1)*I172a33*I175a36)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a33*I176a36)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a36*I177a33)/(18.*(-1 + sw)*(1 + sw)) + (2*ee**2*complex(0,1)*I177a33*I178a36)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I179a36*I180a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I181a33*I182a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I183a33*I184a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I185a36*I186a33)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a33*I169a36)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a36*I171a33)/(8.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I179a36*I180a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I181a33*I182a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I183a33*I184a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I185a36*I186a33*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_764 = Coupling(name = 'GC_764',
                  value = '-(ee**2*complex(0,1)*I168a63*I169a33)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a63*I171a33)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a33*I172a63)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a33*I174a63)/(18.*(-1 + sw)*(1 + sw)) + (2*ee**2*complex(0,1)*I172a63*I175a33)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a63*I176a33)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a63*I177a33)/(18.*(-1 + sw)*(1 + sw)) + (2*ee**2*complex(0,1)*I177a33*I178a63)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I179a33*I180a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I181a33*I182a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I183a63*I184a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I185a63*I186a33)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a63*I169a33)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a63*I171a33)/(8.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I179a33*I180a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I181a33*I182a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I183a63*I184a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I185a63*I186a33*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_765 = Coupling(name = 'GC_765',
                  value = '-(ee**2*complex(0,1)*I168a36*I169a33)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a33*I171a36)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a33*I172a36)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a36*I174a33)/(18.*(-1 + sw)*(1 + sw)) + (2*ee**2*complex(0,1)*I172a36*I175a33)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a36*I176a33)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a33*I177a36)/(18.*(-1 + sw)*(1 + sw)) + (2*ee**2*complex(0,1)*I177a36*I178a33)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I179a33*I180a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I181a36*I182a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I183a36*I184a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I185a33*I186a36)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a36*I169a33)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a33*I171a36)/(8.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I179a33*I180a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I181a36*I182a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I183a36*I184a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I185a33*I186a36*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_766 = Coupling(name = 'GC_766',
                  value = '(-2*ee**2*complex(0,1)*I170a36*I171a36)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a36*I174a36)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a36*I177a36)/(9.*(-1 + sw)*(1 + sw)) + (4*ee**2*complex(0,1)*I177a36*I178a36)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I181a36*I182a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I185a36*I186a36)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I170a36*I171a36)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I181a36*I182a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I185a36*I186a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_767 = Coupling(name = 'GC_767',
                  value = '-(ee**2*complex(0,1)*I168a66*I169a36)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a66*I171a36)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a36*I172a66)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a36*I174a66)/(18.*(-1 + sw)*(1 + sw)) + (2*ee**2*complex(0,1)*I172a66*I175a36)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a66*I176a36)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a66*I177a36)/(18.*(-1 + sw)*(1 + sw)) + (2*ee**2*complex(0,1)*I177a36*I178a66)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I179a36*I180a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I181a36*I182a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I183a66*I184a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I185a66*I186a36)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a66*I169a36)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a66*I171a36)/(8.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I179a36*I180a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I181a36*I182a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I183a66*I184a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I185a66*I186a36*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_768 = Coupling(name = 'GC_768',
                  value = '-(ee**2*complex(0,1)*I168a33*I169a63)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a33*I171a63)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a63*I172a33)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a63*I174a33)/(18.*(-1 + sw)*(1 + sw)) + (2*ee**2*complex(0,1)*I172a33*I175a63)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a33*I176a63)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a33*I177a63)/(18.*(-1 + sw)*(1 + sw)) + (2*ee**2*complex(0,1)*I177a63*I178a33)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I179a63*I180a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I181a63*I182a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I183a33*I184a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I185a33*I186a63)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a33*I169a63)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a33*I171a63)/(8.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I179a63*I180a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I181a63*I182a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I183a33*I184a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I185a33*I186a63*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_769 = Coupling(name = 'GC_769',
                  value = '-(ee**2*complex(0,1)*I168a63*I169a36)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I168a36*I169a63)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a63*I171a36)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a36*I171a63)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a63*I172a36)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a36*I172a63)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a63*I174a36)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a36*I174a63)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a63*I175a36)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a36*I175a63)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a63*I176a36)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a36*I176a63)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a63*I177a36)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a36*I177a63)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a63*I178a36)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a36*I178a63)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I179a63*I180a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I179a36*I180a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I181a63*I182a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I181a36*I182a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I183a63*I184a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I183a36*I184a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I185a63*I186a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I185a36*I186a63)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a63*I169a36)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I168a36*I169a63)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a63*I171a36)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a36*I171a63)/(16.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I179a63*I180a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I179a36*I180a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I181a63*I182a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I181a36*I182a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I183a63*I184a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I183a36*I184a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I185a63*I186a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I185a36*I186a63*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_770 = Coupling(name = 'GC_770',
                  value = '(-2*ee**2*complex(0,1)*I170a63*I171a63)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a63*I174a63)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a63*I177a63)/(9.*(-1 + sw)*(1 + sw)) + (4*ee**2*complex(0,1)*I177a63*I178a63)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I181a63*I182a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I185a63*I186a63)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I170a63*I171a63)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I181a63*I182a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I185a63*I186a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_771 = Coupling(name = 'GC_771',
                  value = '-(ee**2*complex(0,1)*I168a63*I169a66)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a66*I171a63)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a66*I172a63)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a63*I174a66)/(18.*(-1 + sw)*(1 + sw)) + (2*ee**2*complex(0,1)*I172a63*I175a66)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a63*I176a66)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a66*I177a63)/(18.*(-1 + sw)*(1 + sw)) + (2*ee**2*complex(0,1)*I177a63*I178a66)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I179a66*I180a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I181a63*I182a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I183a63*I184a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I185a66*I186a63)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a63*I169a66)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a66*I171a63)/(8.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I179a66*I180a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I181a63*I182a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I183a63*I184a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I185a66*I186a63*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_772 = Coupling(name = 'GC_772',
                  value = '-(ee**2*complex(0,1)*I168a66*I169a33)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I168a33*I169a66)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a66*I171a33)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a33*I171a66)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a66*I172a33)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a33*I172a66)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a66*I174a33)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a33*I174a66)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a66*I175a33)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I172a33*I175a66)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a66*I176a33)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a33*I176a66)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a66*I177a33)/(36.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a33*I177a66)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a66*I178a33)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I177a33*I178a66)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I179a66*I180a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I179a33*I180a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I181a66*I182a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I181a33*I182a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I183a66*I184a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I183a33*I184a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I185a66*I186a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I185a33*I186a66)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a66*I169a33)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I168a33*I169a66)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a66*I171a33)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a33*I171a66)/(16.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I179a66*I180a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I179a33*I180a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I181a66*I182a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I181a33*I182a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I183a66*I184a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I183a33*I184a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I185a66*I186a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I185a33*I186a66*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_773 = Coupling(name = 'GC_773',
                  value = '-(ee**2*complex(0,1)*I168a36*I169a66)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a36*I171a66)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a66*I172a36)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a66*I174a36)/(18.*(-1 + sw)*(1 + sw)) + (2*ee**2*complex(0,1)*I172a36*I175a66)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a36*I176a66)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a36*I177a66)/(18.*(-1 + sw)*(1 + sw)) + (2*ee**2*complex(0,1)*I177a66*I178a36)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I179a66*I180a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I181a66*I182a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I183a36*I184a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I185a36*I186a66)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a36*I169a66)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a36*I171a66)/(8.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I179a66*I180a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I181a66*I182a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I183a36*I184a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I185a36*I186a66*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_774 = Coupling(name = 'GC_774',
                  value = '-(ee**2*complex(0,1)*I168a66*I169a63)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a63*I171a66)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I169a63*I172a66)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a66*I174a63)/(18.*(-1 + sw)*(1 + sw)) + (2*ee**2*complex(0,1)*I172a66*I175a63)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a66*I176a63)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a63*I177a66)/(18.*(-1 + sw)*(1 + sw)) + (2*ee**2*complex(0,1)*I177a66*I178a63)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I179a63*I180a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I181a66*I182a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I183a66*I184a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I185a63*I186a66)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I168a66*I169a63)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I170a63*I171a66)/(8.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I179a63*I180a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I181a66*I182a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I183a66*I184a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I185a63*I186a66*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_775 = Coupling(name = 'GC_775',
                  value = '(-2*ee**2*complex(0,1)*I170a66*I171a66)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I173a66*I174a66)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I170a66*I177a66)/(9.*(-1 + sw)*(1 + sw)) + (4*ee**2*complex(0,1)*I177a66*I178a66)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I181a66*I182a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I185a66*I186a66)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I170a66*I171a66)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I181a66*I182a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I185a66*I186a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_776 = Coupling(name = 'GC_776',
                  value = '(4*ee**2*complex(0,1)*I201a11)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I201a11)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (8*ee**2*complex(0,1)*I201a11*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_777 = Coupling(name = 'GC_777',
                  value = '(4*ee**2*complex(0,1)*I201a22)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I201a22)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (8*ee**2*complex(0,1)*I201a22*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_778 = Coupling(name = 'GC_778',
                  value = '(4*ee**2*complex(0,1)*I201a33)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I201a33)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (8*ee**2*complex(0,1)*I201a33*sw**2)/(9.*(-1 + sw)*(1 + sw)) - (8*ee**2*complex(0,1)*I203a33*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_779 = Coupling(name = 'GC_779',
                  value = '(4*ee**2*complex(0,1)*I201a36)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I201a36)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (8*ee**2*complex(0,1)*I201a36*sw**2)/(9.*(-1 + sw)*(1 + sw)) - (8*ee**2*complex(0,1)*I203a36*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_780 = Coupling(name = 'GC_780',
                  value = '(4*ee**2*complex(0,1)*I201a63)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I201a63)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (8*ee**2*complex(0,1)*I201a63*sw**2)/(9.*(-1 + sw)*(1 + sw)) - (8*ee**2*complex(0,1)*I203a63*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_781 = Coupling(name = 'GC_781',
                  value = '(4*ee**2*complex(0,1)*I201a66)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I201a66)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (8*ee**2*complex(0,1)*I201a66*sw**2)/(9.*(-1 + sw)*(1 + sw)) - (8*ee**2*complex(0,1)*I203a66*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_782 = Coupling(name = 'GC_782',
                  value = '(complex(0,1)*G**2*I26a11*I27a11)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a11*I27a11*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_783 = Coupling(name = 'GC_783',
                  value = '(complex(0,1)*G**2*I26a22*I27a22)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a22*I27a22*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_784 = Coupling(name = 'GC_784',
                  value = '(complex(0,1)*G**2*I28a11*I29a11)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a11*I29a11*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_785 = Coupling(name = 'GC_785',
                  value = '(complex(0,1)*G**2*I26a22*I27a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I26a11*I27a22)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a22*I29a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a11*I29a22)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a22*I27a11*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a11*I27a22*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a22*I29a11*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a11*I29a22*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_786 = Coupling(name = 'GC_786',
                  value = '(complex(0,1)*G**2*I28a22*I29a22)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a22*I29a22*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_787 = Coupling(name = 'GC_787',
                  value = '(complex(0,1)*G**2*I30a44*I33a44)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a44*I33a44*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_788 = Coupling(name = 'GC_788',
                  value = '(complex(0,1)*G**2*I30a55*I33a55)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a55*I33a55*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_789 = Coupling(name = 'GC_789',
                  value = '(complex(0,1)*G**2*I28a33*I29a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a33*I32a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a33*I35a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a33*I36a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a33*I29a33*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a33*I32a33*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a33*I35a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a33*I36a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_790 = Coupling(name = 'GC_790',
                  value = '(complex(0,1)*G**2*I28a36*I29a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a36*I32a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a36*I35a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a36*I36a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a36*I29a36*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a36*I32a36*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a36*I35a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a36*I36a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_791 = Coupling(name = 'GC_791',
                  value = '(complex(0,1)*G**2*I35a44*I36a44)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a44*I36a44*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_792 = Coupling(name = 'GC_792',
                  value = '(complex(0,1)*G**2*I30a55*I33a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a44*I33a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a55*I36a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a44*I36a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a55*I33a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a44*I33a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a55*I36a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a44*I36a55*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_793 = Coupling(name = 'GC_793',
                  value = '(complex(0,1)*G**2*I35a55*I36a55)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a55*I36a55*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_794 = Coupling(name = 'GC_794',
                  value = '(complex(0,1)*G**2*I28a63*I29a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a63*I32a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a63*I35a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a63*I36a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a63*I29a63*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a63*I32a63*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a63*I35a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a63*I36a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_795 = Coupling(name = 'GC_795',
                  value = '(complex(0,1)*G**2*I28a66*I29a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a66*I32a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a66*I35a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a66*I36a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a66*I29a66*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a66*I32a66*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a66*I35a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a66*I36a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_796 = Coupling(name = 'GC_796',
                  value = '(-2*ee**2*complex(0,1)*I26a33*I27a33)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a33*I30a33)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a33*I33a33)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I37a33*I38a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I41a33*I42a33)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a33*I8a33)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a33*I27a33)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I37a33*I38a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I41a33*I42a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_797 = Coupling(name = 'GC_797',
                  value = '(-2*ee**2*complex(0,1)*I26a36*I27a36)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a36*I30a36)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a36*I33a36)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I37a36*I38a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I41a36*I42a36)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a36*I8a36)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a36*I27a36)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I37a36*I38a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I41a36*I42a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_798 = Coupling(name = 'GC_798',
                  value = '(-2*ee**2*complex(0,1)*I26a63*I27a63)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a63*I30a63)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a63*I33a63)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I37a63*I38a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I41a63*I42a63)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a63*I8a63)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a63*I27a63)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I37a63*I38a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I41a63*I42a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_799 = Coupling(name = 'GC_799',
                  value = '(-2*ee**2*complex(0,1)*I26a66*I27a66)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a66*I30a66)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a66*I33a66)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I37a66*I38a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I41a66*I42a66)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a66*I8a66)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a66*I27a66)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I37a66*I38a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I41a66*I42a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_800 = Coupling(name = 'GC_800',
                  value = '(-2*ee**2*complex(0,1)*I28a33*I29a33)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a33*I32a33)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a33*I35a33)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a33*I36a33)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I39a33*I40a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I43a33*I44a33)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a33*I29a33)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I39a33*I40a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I43a33*I44a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_801 = Coupling(name = 'GC_801',
                  value = '-(ee**2*complex(0,1)*I26a33*I27a36)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a36*I29a33)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a36*I30a33)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a33*I32a36)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a33*I33a36)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a36*I35a33)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a33*I36a36)/(18.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I37a36*I38a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I39a33*I40a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I41a33*I42a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I43a36*I44a33)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a36*I8a33)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a33*I27a36)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a36*I29a33)/(8.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I37a36*I38a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I39a33*I40a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I41a33*I42a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I43a36*I44a33*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_802 = Coupling(name = 'GC_802',
                  value = '-(ee**2*complex(0,1)*I26a63*I27a33)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a63*I29a33)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a33*I30a63)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a33*I32a63)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a63*I33a33)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a63*I35a33)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a33*I36a63)/(18.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I37a33*I38a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I39a33*I40a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I41a63*I42a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I43a63*I44a33)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a33*I8a63)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a63*I27a33)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a63*I29a33)/(8.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I37a33*I38a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I39a33*I40a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I41a63*I42a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I43a63*I44a33*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_803 = Coupling(name = 'GC_803',
                  value = '-(ee**2*complex(0,1)*I26a36*I27a33)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a33*I29a36)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a33*I30a36)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a36*I32a33)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a36*I33a33)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a33*I35a36)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a36*I36a33)/(18.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I37a33*I38a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I39a36*I40a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I41a36*I42a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I43a33*I44a36)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a33*I8a36)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a36*I27a33)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a33*I29a36)/(8.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I37a33*I38a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I39a36*I40a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I41a36*I42a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I43a33*I44a36*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_804 = Coupling(name = 'GC_804',
                  value = '(-2*ee**2*complex(0,1)*I28a36*I29a36)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a36*I32a36)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a36*I35a36)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a36*I36a36)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I39a36*I40a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I43a36*I44a36)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a36*I29a36)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I39a36*I40a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I43a36*I44a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_805 = Coupling(name = 'GC_805',
                  value = '-(ee**2*complex(0,1)*I26a66*I27a36)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a66*I29a36)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a36*I30a66)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a36*I32a66)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a66*I33a36)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a66*I35a36)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a36*I36a66)/(18.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I37a36*I38a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I39a36*I40a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I41a66*I42a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I43a66*I44a36)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a36*I8a66)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a66*I27a36)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a66*I29a36)/(8.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I37a36*I38a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I39a36*I40a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I41a66*I42a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I43a66*I44a36*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_806 = Coupling(name = 'GC_806',
                  value = '-(ee**2*complex(0,1)*I26a33*I27a63)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a33*I29a63)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a63*I30a33)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a63*I32a33)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a33*I33a63)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a33*I35a63)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a63*I36a33)/(18.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I37a63*I38a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I39a63*I40a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I41a33*I42a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I43a33*I44a63)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a63*I8a33)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a33*I27a63)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a33*I29a63)/(8.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I37a63*I38a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I39a63*I40a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I41a33*I42a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I43a33*I44a63*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_807 = Coupling(name = 'GC_807',
                  value = '-(ee**2*complex(0,1)*I26a63*I27a36)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I26a36*I27a63)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a63*I29a36)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a36*I29a63)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a63*I30a36)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a36*I30a63)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a63*I32a36)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a36*I32a63)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a63*I33a36)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a36*I33a63)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a63*I35a36)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a36*I35a63)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a63*I36a36)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a36*I36a63)/(36.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I37a63*I38a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I37a36*I38a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I39a63*I40a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I39a36*I40a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I41a63*I42a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I41a36*I42a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I43a63*I44a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I43a36*I44a63)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a63*I8a36)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a36*I8a63)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a63*I27a36)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I26a36*I27a63)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a63*I29a36)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a36*I29a63)/(16.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I37a63*I38a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I37a36*I38a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I39a63*I40a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I39a36*I40a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I41a63*I42a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I41a36*I42a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I43a63*I44a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I43a36*I44a63*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_808 = Coupling(name = 'GC_808',
                  value = '(-2*ee**2*complex(0,1)*I28a63*I29a63)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a63*I32a63)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a63*I35a63)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a63*I36a63)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I39a63*I40a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I43a63*I44a63)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a63*I29a63)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I39a63*I40a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I43a63*I44a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_809 = Coupling(name = 'GC_809',
                  value = '-(ee**2*complex(0,1)*I26a63*I27a66)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a66*I29a63)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a66*I30a63)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a63*I32a66)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a63*I33a66)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a66*I35a63)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a63*I36a66)/(18.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I37a66*I38a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I39a63*I40a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I41a63*I42a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I43a66*I44a63)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a66*I8a63)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a63*I27a66)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a66*I29a63)/(8.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I37a66*I38a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I39a63*I40a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I41a63*I42a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I43a66*I44a63*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_810 = Coupling(name = 'GC_810',
                  value = '-(ee**2*complex(0,1)*I26a66*I27a33)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I26a33*I27a66)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a66*I29a33)/(18.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a33*I29a66)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a66*I30a33)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a33*I30a66)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a66*I32a33)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a33*I32a66)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a66*I33a33)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a33*I33a66)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a66*I35a33)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a33*I35a66)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a66*I36a33)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a33*I36a66)/(36.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I37a66*I38a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I37a33*I38a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I39a66*I40a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I39a33*I40a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I41a66*I42a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I41a33*I42a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I43a66*I44a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I43a33*I44a66)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a66*I8a33)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a33*I8a66)/(72.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a66*I27a33)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I26a33*I27a66)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a66*I29a33)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a33*I29a66)/(16.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I37a66*I38a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I37a33*I38a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I39a66*I40a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I39a33*I40a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I41a66*I42a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I41a33*I42a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I43a66*I44a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I43a33*I44a66*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_811 = Coupling(name = 'GC_811',
                  value = '-(ee**2*complex(0,1)*I26a36*I27a66)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a36*I29a66)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a66*I30a36)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a66*I32a36)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a36*I33a66)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a36*I35a66)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a66*I36a36)/(18.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I37a66*I38a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I39a66*I40a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I41a36*I42a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I43a36*I44a66)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a66*I8a36)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a36*I27a66)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a36*I29a66)/(8.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I37a66*I38a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I39a66*I40a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I41a36*I42a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I43a36*I44a66*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_812 = Coupling(name = 'GC_812',
                  value = '-(ee**2*complex(0,1)*I26a66*I27a63)/(9.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I28a63*I29a66)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I27a63*I30a66)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a66*I32a63)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I30a66*I33a63)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a63*I35a66)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a66*I36a63)/(18.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I37a63*I38a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I39a66*I40a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I41a66*I42a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I43a63*I44a66)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I34a63*I8a66)/(36.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I26a66*I27a63)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I28a63*I29a66)/(8.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I37a63*I38a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I39a66*I40a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I41a66*I42a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I43a63*I44a66*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_813 = Coupling(name = 'GC_813',
                  value = '(-2*ee**2*complex(0,1)*I28a66*I29a66)/(9.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I31a66*I32a66)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a66*I35a66)/(18.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I35a66*I36a66)/(9.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I39a66*I40a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I43a66*I44a66)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I28a66*I29a66)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I39a66*I40a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I43a66*I44a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_814 = Coupling(name = 'GC_814',
                  value = '(2*ee**2*complex(0,1)*I63a11)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I63a11)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I63a11*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_815 = Coupling(name = 'GC_815',
                  value = '(complex(0,1)*G**2*I161a11*I63a11)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a11*I63a11*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_816 = Coupling(name = 'GC_816',
                  value = '(complex(0,1)*G**2*I161a22*I63a11)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a22*I63a11*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_817 = Coupling(name = 'GC_817',
                  value = '(complex(0,1)*G**2*I161a33*I63a11)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a33*I63a11)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a33*I63a11*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a33*I63a11*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_818 = Coupling(name = 'GC_818',
                  value = '(complex(0,1)*G**2*I161a36*I63a11)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a36*I63a11)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a36*I63a11*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a36*I63a11*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_819 = Coupling(name = 'GC_819',
                  value = '-((complex(0,1)*G**2*I162a44*I63a11)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I162a44*I63a11*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_820 = Coupling(name = 'GC_820',
                  value = '-((complex(0,1)*G**2*I162a55*I63a11)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I162a55*I63a11*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_821 = Coupling(name = 'GC_821',
                  value = '(complex(0,1)*G**2*I161a63*I63a11)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a63*I63a11)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a63*I63a11*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a63*I63a11*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_822 = Coupling(name = 'GC_822',
                  value = '(complex(0,1)*G**2*I161a66*I63a11)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a66*I63a11)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a66*I63a11*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a66*I63a11*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_823 = Coupling(name = 'GC_823',
                  value = '(2*ee**2*complex(0,1)*I63a22)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I63a22)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I63a22*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_824 = Coupling(name = 'GC_824',
                  value = '(complex(0,1)*G**2*I161a11*I63a22)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a11*I63a22*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_825 = Coupling(name = 'GC_825',
                  value = '(complex(0,1)*G**2*I161a22*I63a22)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a22*I63a22*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_826 = Coupling(name = 'GC_826',
                  value = '(complex(0,1)*G**2*I161a33*I63a22)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a33*I63a22)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a33*I63a22*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a33*I63a22*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_827 = Coupling(name = 'GC_827',
                  value = '(complex(0,1)*G**2*I161a36*I63a22)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a36*I63a22)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a36*I63a22*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a36*I63a22*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_828 = Coupling(name = 'GC_828',
                  value = '-((complex(0,1)*G**2*I162a44*I63a22)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I162a44*I63a22*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_829 = Coupling(name = 'GC_829',
                  value = '-((complex(0,1)*G**2*I162a55*I63a22)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I162a55*I63a22*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_830 = Coupling(name = 'GC_830',
                  value = '(complex(0,1)*G**2*I161a63*I63a22)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a63*I63a22)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a63*I63a22*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a63*I63a22*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_831 = Coupling(name = 'GC_831',
                  value = '(complex(0,1)*G**2*I161a66*I63a22)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a66*I63a22)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a66*I63a22*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a66*I63a22*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_832 = Coupling(name = 'GC_832',
                  value = '(2*ee**2*complex(0,1)*I63a33)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I63a33)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I63a33*sw**2)/(9.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I65a33*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_833 = Coupling(name = 'GC_833',
                  value = '(complex(0,1)*G**2*I161a11*I63a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a11*I65a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a11*I63a33*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a11*I65a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_834 = Coupling(name = 'GC_834',
                  value = '(complex(0,1)*G**2*I161a22*I63a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a22*I65a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a22*I63a33*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a22*I65a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_835 = Coupling(name = 'GC_835',
                  value = '(complex(0,1)*G**2*I161a33*I63a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a33*I63a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a33*I65a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a33*I65a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a33*I63a33*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a33*I63a33*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a33*I65a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a33*I65a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_836 = Coupling(name = 'GC_836',
                  value = '(complex(0,1)*G**2*I161a36*I63a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a36*I63a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a36*I65a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a36*I65a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a36*I63a33*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a36*I63a33*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a36*I65a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a36*I65a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_837 = Coupling(name = 'GC_837',
                  value = '-((complex(0,1)*G**2*I162a44*I63a33)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I162a44*I65a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a44*I63a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a44*I65a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_838 = Coupling(name = 'GC_838',
                  value = '-((complex(0,1)*G**2*I162a55*I63a33)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I162a55*I65a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a55*I63a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a55*I65a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_839 = Coupling(name = 'GC_839',
                  value = '(complex(0,1)*G**2*I161a63*I63a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a63*I63a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a63*I65a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a63*I65a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a63*I63a33*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a63*I63a33*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a63*I65a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a63*I65a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_840 = Coupling(name = 'GC_840',
                  value = '(complex(0,1)*G**2*I161a66*I63a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a66*I63a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a66*I65a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a66*I65a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a66*I63a33*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a66*I63a33*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a66*I65a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a66*I65a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_841 = Coupling(name = 'GC_841',
                  value = '(2*ee**2*complex(0,1)*I63a36)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I63a36)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I63a36*sw**2)/(9.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I65a36*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_842 = Coupling(name = 'GC_842',
                  value = '(complex(0,1)*G**2*I161a11*I63a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a11*I65a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a11*I63a36*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a11*I65a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_843 = Coupling(name = 'GC_843',
                  value = '(complex(0,1)*G**2*I161a22*I63a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a22*I65a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a22*I63a36*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a22*I65a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_844 = Coupling(name = 'GC_844',
                  value = '(complex(0,1)*G**2*I161a33*I63a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a33*I63a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a33*I65a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a33*I65a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a33*I63a36*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a33*I63a36*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a33*I65a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a33*I65a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_845 = Coupling(name = 'GC_845',
                  value = '(complex(0,1)*G**2*I161a36*I63a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a36*I63a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a36*I65a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a36*I65a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a36*I63a36*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a36*I63a36*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a36*I65a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a36*I65a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_846 = Coupling(name = 'GC_846',
                  value = '-((complex(0,1)*G**2*I162a44*I63a36)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I162a44*I65a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a44*I63a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a44*I65a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_847 = Coupling(name = 'GC_847',
                  value = '-((complex(0,1)*G**2*I162a55*I63a36)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I162a55*I65a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a55*I63a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a55*I65a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_848 = Coupling(name = 'GC_848',
                  value = '(complex(0,1)*G**2*I161a63*I63a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a63*I63a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a63*I65a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a63*I65a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a63*I63a36*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a63*I63a36*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a63*I65a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a63*I65a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_849 = Coupling(name = 'GC_849',
                  value = '(complex(0,1)*G**2*I161a66*I63a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a66*I63a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a66*I65a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a66*I65a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a66*I63a36*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a66*I63a36*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a66*I65a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a66*I65a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_850 = Coupling(name = 'GC_850',
                  value = '-((complex(0,1)*G**2*I161a11*I65a44)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I161a11*I65a44*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_851 = Coupling(name = 'GC_851',
                  value = '-((complex(0,1)*G**2*I161a22*I65a44)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I161a22*I65a44*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_852 = Coupling(name = 'GC_852',
                  value = '-((complex(0,1)*G**2*I161a33*I65a44)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I162a33*I65a44)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a33*I65a44*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a33*I65a44*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_853 = Coupling(name = 'GC_853',
                  value = '-((complex(0,1)*G**2*I161a36*I65a44)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I162a36*I65a44)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a36*I65a44*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a36*I65a44*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_854 = Coupling(name = 'GC_854',
                  value = '(complex(0,1)*G**2*I162a44*I65a44)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a44*I65a44*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_855 = Coupling(name = 'GC_855',
                  value = '(complex(0,1)*G**2*I162a55*I65a44)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a55*I65a44*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_856 = Coupling(name = 'GC_856',
                  value = '-((complex(0,1)*G**2*I161a63*I65a44)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I162a63*I65a44)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a63*I65a44*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a63*I65a44*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_857 = Coupling(name = 'GC_857',
                  value = '-((complex(0,1)*G**2*I161a66*I65a44)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I162a66*I65a44)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a66*I65a44*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a66*I65a44*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_858 = Coupling(name = 'GC_858',
                  value = '-((complex(0,1)*G**2*I161a11*I65a55)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I161a11*I65a55*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_859 = Coupling(name = 'GC_859',
                  value = '-((complex(0,1)*G**2*I161a22*I65a55)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I161a22*I65a55*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_860 = Coupling(name = 'GC_860',
                  value = '-((complex(0,1)*G**2*I161a33*I65a55)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I162a33*I65a55)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a33*I65a55*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a33*I65a55*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_861 = Coupling(name = 'GC_861',
                  value = '-((complex(0,1)*G**2*I161a36*I65a55)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I162a36*I65a55)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a36*I65a55*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a36*I65a55*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_862 = Coupling(name = 'GC_862',
                  value = '(complex(0,1)*G**2*I162a44*I65a55)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a44*I65a55*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_863 = Coupling(name = 'GC_863',
                  value = '(complex(0,1)*G**2*I162a55*I65a55)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a55*I65a55*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_864 = Coupling(name = 'GC_864',
                  value = '-((complex(0,1)*G**2*I161a63*I65a55)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I162a63*I65a55)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a63*I65a55*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a63*I65a55*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_865 = Coupling(name = 'GC_865',
                  value = '-((complex(0,1)*G**2*I161a66*I65a55)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I162a66*I65a55)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a66*I65a55*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a66*I65a55*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_866 = Coupling(name = 'GC_866',
                  value = '(2*ee**2*complex(0,1)*I63a63)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I63a63)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I63a63*sw**2)/(9.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I65a63*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_867 = Coupling(name = 'GC_867',
                  value = '(complex(0,1)*G**2*I161a11*I63a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a11*I65a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a11*I63a63*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a11*I65a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_868 = Coupling(name = 'GC_868',
                  value = '(complex(0,1)*G**2*I161a22*I63a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a22*I65a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a22*I63a63*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a22*I65a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_869 = Coupling(name = 'GC_869',
                  value = '(complex(0,1)*G**2*I161a33*I63a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a33*I63a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a33*I65a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a33*I65a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a33*I63a63*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a33*I63a63*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a33*I65a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a33*I65a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_870 = Coupling(name = 'GC_870',
                  value = '(complex(0,1)*G**2*I161a36*I63a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a36*I63a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a36*I65a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a36*I65a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a36*I63a63*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a36*I63a63*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a36*I65a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a36*I65a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_871 = Coupling(name = 'GC_871',
                  value = '-((complex(0,1)*G**2*I162a44*I63a63)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I162a44*I65a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a44*I63a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a44*I65a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_872 = Coupling(name = 'GC_872',
                  value = '-((complex(0,1)*G**2*I162a55*I63a63)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I162a55*I65a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a55*I63a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a55*I65a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_873 = Coupling(name = 'GC_873',
                  value = '(complex(0,1)*G**2*I161a63*I63a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a63*I63a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a63*I65a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a63*I65a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a63*I63a63*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a63*I63a63*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a63*I65a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a63*I65a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_874 = Coupling(name = 'GC_874',
                  value = '(complex(0,1)*G**2*I161a66*I63a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a66*I63a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a66*I65a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a66*I65a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a66*I63a63*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a66*I63a63*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a66*I65a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a66*I65a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_875 = Coupling(name = 'GC_875',
                  value = '(2*ee**2*complex(0,1)*I63a66)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I63a66)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I63a66*sw**2)/(9.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I65a66*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_876 = Coupling(name = 'GC_876',
                  value = '(complex(0,1)*G**2*I161a11*I63a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a11*I65a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a11*I63a66*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a11*I65a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_877 = Coupling(name = 'GC_877',
                  value = '(complex(0,1)*G**2*I161a22*I63a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a22*I65a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a22*I63a66*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a22*I65a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_878 = Coupling(name = 'GC_878',
                  value = '(complex(0,1)*G**2*I161a33*I63a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a33*I63a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a33*I65a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a33*I65a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a33*I63a66*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a33*I63a66*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a33*I65a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a33*I65a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_879 = Coupling(name = 'GC_879',
                  value = '(complex(0,1)*G**2*I161a36*I63a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a36*I63a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a36*I65a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a36*I65a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a36*I63a66*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a36*I63a66*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a36*I65a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a36*I65a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_880 = Coupling(name = 'GC_880',
                  value = '-((complex(0,1)*G**2*I162a44*I63a66)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I162a44*I65a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a44*I63a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a44*I65a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_881 = Coupling(name = 'GC_881',
                  value = '-((complex(0,1)*G**2*I162a55*I63a66)/((-1 + sw)*(1 + sw))) + (complex(0,1)*G**2*I162a55*I65a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a55*I63a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a55*I65a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_882 = Coupling(name = 'GC_882',
                  value = '(complex(0,1)*G**2*I161a63*I63a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a63*I63a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a63*I65a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a63*I65a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a63*I63a66*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a63*I63a66*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a63*I65a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a63*I65a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_883 = Coupling(name = 'GC_883',
                  value = '(complex(0,1)*G**2*I161a66*I63a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a66*I63a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a66*I65a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a66*I65a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I161a66*I63a66*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I162a66*I63a66*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I161a66*I65a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I162a66*I65a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_884 = Coupling(name = 'GC_884',
                  value = '-(ee**2*complex(0,1)*I63a33*I64a33)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a33*I65a33)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a33*I66a33)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a33*I66a33)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67a33*I68a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I69a33*I70a33)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a33*I64a33)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I67a33*I68a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I69a33*I70a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_885 = Coupling(name = 'GC_885',
                  value = '-(ee**2*complex(0,1)*I63a36*I64a33)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a33*I65a36)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a36*I66a33)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a36*I66a33)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67a33*I68a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I69a36*I70a33)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a36*I64a33)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I67a33*I68a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I69a36*I70a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_886 = Coupling(name = 'GC_886',
                  value = '-(ee**2*complex(0,1)*I63a63*I64a33)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a33*I65a63)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a63*I66a33)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a63*I66a33)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67a33*I68a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I69a63*I70a33)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a63*I64a33)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I67a33*I68a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I69a63*I70a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_887 = Coupling(name = 'GC_887',
                  value = '-(ee**2*complex(0,1)*I63a66*I64a33)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a33*I65a66)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a66*I66a33)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a66*I66a33)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67a33*I68a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I69a66*I70a33)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a66*I64a33)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I67a33*I68a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I69a66*I70a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_888 = Coupling(name = 'GC_888',
                  value = '-(ee**2*complex(0,1)*I63a33*I64a36)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a36*I65a33)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a33*I66a36)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a33*I66a36)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67a36*I68a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I69a33*I70a36)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a33*I64a36)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I67a36*I68a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I69a33*I70a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_889 = Coupling(name = 'GC_889',
                  value = '-(ee**2*complex(0,1)*I63a36*I64a36)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a36*I65a36)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a36*I66a36)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a36*I66a36)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67a36*I68a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I69a36*I70a36)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a36*I64a36)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I67a36*I68a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I69a36*I70a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_890 = Coupling(name = 'GC_890',
                  value = '-(ee**2*complex(0,1)*I63a63*I64a36)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a36*I65a63)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a63*I66a36)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a63*I66a36)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67a36*I68a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I69a63*I70a36)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a63*I64a36)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I67a36*I68a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I69a63*I70a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_891 = Coupling(name = 'GC_891',
                  value = '-(ee**2*complex(0,1)*I63a66*I64a36)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a36*I65a66)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a66*I66a36)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a66*I66a36)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67a36*I68a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I69a66*I70a36)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a66*I64a36)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I67a36*I68a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I69a66*I70a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_892 = Coupling(name = 'GC_892',
                  value = '-(ee**2*complex(0,1)*I63a33*I64a63)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a63*I65a33)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a33*I66a63)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a33*I66a63)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67a63*I68a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I69a33*I70a63)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a33*I64a63)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I67a63*I68a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I69a33*I70a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_893 = Coupling(name = 'GC_893',
                  value = '-(ee**2*complex(0,1)*I63a36*I64a63)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a63*I65a36)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a36*I66a63)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a36*I66a63)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67a63*I68a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I69a36*I70a63)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a36*I64a63)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I67a63*I68a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I69a36*I70a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_894 = Coupling(name = 'GC_894',
                  value = '-(ee**2*complex(0,1)*I63a63*I64a63)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a63*I65a63)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a63*I66a63)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a63*I66a63)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67a63*I68a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I69a63*I70a63)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a63*I64a63)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I67a63*I68a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I69a63*I70a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_895 = Coupling(name = 'GC_895',
                  value = '-(ee**2*complex(0,1)*I63a66*I64a63)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a63*I65a66)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a66*I66a63)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a66*I66a63)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67a63*I68a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I69a66*I70a63)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a66*I64a63)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I67a63*I68a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I69a66*I70a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_896 = Coupling(name = 'GC_896',
                  value = '-(ee**2*complex(0,1)*I63a33*I64a66)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a66*I65a33)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a33*I66a66)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a33*I66a66)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67a66*I68a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I69a33*I70a66)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a33*I64a66)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I67a66*I68a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I69a33*I70a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_897 = Coupling(name = 'GC_897',
                  value = '-(ee**2*complex(0,1)*I63a36*I64a66)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a66*I65a36)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a36*I66a66)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a36*I66a66)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67a66*I68a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I69a36*I70a66)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a36*I64a66)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I67a66*I68a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I69a36*I70a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_898 = Coupling(name = 'GC_898',
                  value = '-(ee**2*complex(0,1)*I63a63*I64a66)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a66*I65a63)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a63*I66a66)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a63*I66a66)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67a66*I68a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I69a63*I70a66)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a63*I64a66)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I67a66*I68a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I69a63*I70a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_899 = Coupling(name = 'GC_899',
                  value = '-(ee**2*complex(0,1)*I63a66*I64a66)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I64a66*I65a66)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a66*I66a66)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I65a66*I66a66)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67a66*I68a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I69a66*I70a66)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I63a66*I64a66)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I67a66*I68a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I69a66*I70a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_900 = Coupling(name = 'GC_900',
                  value = '-(ee**2*complex(0,1)*I75a33*I76a33)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a33*I77a33)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a33*I78a33)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a33*I79a33)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a33*I80a33)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a33*I81a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I82a33*I83a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I84a33*I85a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I86a33*I87a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I88a33*I89a33)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I71a33*I72a33)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a33*I74a33)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I82a33*I83a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I84a33*I85a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I86a33*I87a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I88a33*I89a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_901 = Coupling(name = 'GC_901',
                  value = '-(ee**2*complex(0,1)*I75a36*I76a33)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I75a33*I76a36)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a36*I77a33)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a33*I77a36)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a36*I78a33)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a33*I78a36)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a36*I79a33)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a33*I79a36)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a36*I80a33)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a33*I80a36)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a36*I81a33)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a33*I81a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I82a36*I83a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I82a33*I83a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I84a36*I85a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I84a33*I85a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I86a36*I87a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I86a33*I87a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I88a36*I89a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I88a33*I89a36)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I71a36*I72a33)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I71a33*I72a36)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a36*I74a33)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a33*I74a36)/(8.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I82a36*I83a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I82a33*I83a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I84a36*I85a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I84a33*I85a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I86a36*I87a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I86a33*I87a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I88a36*I89a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I88a33*I89a36*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_902 = Coupling(name = 'GC_902',
                  value = '-(ee**2*complex(0,1)*I75a36*I76a36)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a36*I77a36)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a36*I78a36)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a36*I79a36)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a36*I80a36)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a36*I81a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I82a36*I83a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I84a36*I85a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I86a36*I87a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I88a36*I89a36)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I71a36*I72a36)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a36*I74a36)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I82a36*I83a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I84a36*I85a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I86a36*I87a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I88a36*I89a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_903 = Coupling(name = 'GC_903',
                  value = '-(ee**2*complex(0,1)*I75a63*I76a33)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I75a33*I76a63)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a63*I77a33)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a33*I77a63)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a63*I78a33)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a33*I78a63)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a63*I79a33)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a33*I79a63)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a63*I80a33)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a33*I80a63)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a63*I81a33)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a33*I81a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I82a63*I83a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I82a33*I83a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I84a63*I85a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I84a33*I85a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I86a63*I87a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I86a33*I87a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I88a63*I89a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I88a33*I89a63)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I71a63*I72a33)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I71a33*I72a63)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a63*I74a33)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a33*I74a63)/(8.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I82a63*I83a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I82a33*I83a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I84a63*I85a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I84a33*I85a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I86a63*I87a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I86a33*I87a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I88a63*I89a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I88a33*I89a63*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_904 = Coupling(name = 'GC_904',
                  value = '-(ee**2*complex(0,1)*I75a63*I76a63)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a63*I77a63)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a63*I78a63)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a63*I79a63)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a63*I80a63)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a63*I81a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I82a63*I83a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I84a63*I85a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I86a63*I87a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I88a63*I89a63)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I71a63*I72a63)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a63*I74a63)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I82a63*I83a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I84a63*I85a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I86a63*I87a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I88a63*I89a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_905 = Coupling(name = 'GC_905',
                  value = '-(ee**2*complex(0,1)*I75a66*I76a33)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I75a63*I76a36)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I75a36*I76a63)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I75a33*I76a66)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a66*I77a33)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a63*I77a36)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a36*I77a63)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a33*I77a66)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a66*I78a33)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a63*I78a36)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a36*I78a63)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a33*I78a66)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a66*I79a33)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a63*I79a36)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a36*I79a63)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a33*I79a66)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a66*I80a33)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a63*I80a36)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a36*I80a63)/(8.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a33*I80a66)/(8.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a66*I81a33)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a63*I81a36)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a36*I81a63)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a33*I81a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I82a66*I83a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I82a63*I83a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I82a36*I83a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I82a33*I83a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I84a66*I85a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I84a63*I85a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I84a36*I85a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I84a33*I85a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I86a66*I87a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I86a63*I87a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I86a36*I87a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I86a33*I87a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I88a66*I89a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I88a63*I89a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I88a36*I89a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I88a33*I89a66)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I71a66*I72a33)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I71a63*I72a36)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I71a36*I72a63)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I71a33*I72a66)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a66*I74a33)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a63*I74a36)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a36*I74a63)/(16.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a33*I74a66)/(16.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I82a66*I83a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I82a63*I83a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I82a36*I83a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I82a33*I83a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I84a66*I85a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I84a63*I85a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I84a36*I85a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I84a33*I85a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I86a66*I87a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I86a63*I87a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I86a36*I87a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I86a33*I87a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I88a66*I89a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I88a63*I89a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I88a36*I89a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I88a33*I89a66*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_906 = Coupling(name = 'GC_906',
                  value = '-(ee**2*complex(0,1)*I75a66*I76a36)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I75a36*I76a66)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a66*I77a36)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a36*I77a66)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a66*I78a36)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a36*I78a66)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a66*I79a36)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a36*I79a66)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a66*I80a36)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a36*I80a66)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a66*I81a36)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a36*I81a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I82a66*I83a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I82a36*I83a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I84a66*I85a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I84a36*I85a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I86a66*I87a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I86a36*I87a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I88a66*I89a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I88a36*I89a66)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I71a66*I72a36)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I71a36*I72a66)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a66*I74a36)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a36*I74a66)/(8.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I82a66*I83a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I82a36*I83a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I84a66*I85a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I84a36*I85a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I86a66*I87a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I86a36*I87a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I88a66*I89a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I88a36*I89a66*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_907 = Coupling(name = 'GC_907',
                  value = '-(ee**2*complex(0,1)*I75a66*I76a63)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I75a63*I76a66)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a66*I77a63)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a63*I77a66)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a66*I78a63)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a63*I78a66)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a66*I79a63)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a63*I79a66)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a66*I80a63)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a63*I80a66)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a66*I81a63)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a63*I81a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I82a66*I83a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I82a63*I83a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I84a66*I85a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I84a63*I85a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I86a66*I87a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I86a63*I87a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I88a66*I89a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I88a63*I89a66)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I71a66*I72a63)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I71a63*I72a66)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a66*I74a63)/(8.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a63*I74a66)/(8.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I82a66*I83a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I82a63*I83a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I84a66*I85a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I84a63*I85a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I86a66*I87a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I86a63*I87a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I88a66*I89a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I88a63*I89a66*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_908 = Coupling(name = 'GC_908',
                  value = '-(ee**2*complex(0,1)*I75a66*I76a66)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I72a66*I77a66)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I77a66*I78a66)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a66*I79a66)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I73a66*I80a66)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I80a66*I81a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I82a66*I83a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I84a66*I85a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I86a66*I87a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I88a66*I89a66)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I71a66*I72a66)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I73a66*I74a66)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I82a66*I83a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I84a66*I85a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I86a66*I87a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I88a66*I89a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_909 = Coupling(name = 'GC_909',
                  value = '(complex(0,1)*G**2*I26a33*I27a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I26a11*I27a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a33*I29a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a11*I29a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a11*I30a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a11*I32a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a11*I35a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a33*I8a11)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a33*I27a11*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a11*I27a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a33*I29a11*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a11*I29a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a11*I30a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a11*I32a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a11*I35a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a33*I8a11*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_910 = Coupling(name = 'GC_910',
                  value = '(complex(0,1)*G**2*I26a36*I27a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I26a11*I27a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a36*I29a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a11*I29a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a11*I30a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a11*I32a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a11*I35a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a36*I8a11)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a36*I27a11*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a11*I27a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a36*I29a11*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a11*I29a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a11*I30a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a11*I32a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a11*I35a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a36*I8a11*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_911 = Coupling(name = 'GC_911',
                  value = '-(complex(0,1)*G**2*I27a11*I30a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a11*I32a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a11*I35a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a44*I8a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a11*I30a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a11*I32a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a11*I35a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a44*I8a11*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_912 = Coupling(name = 'GC_912',
                  value = '-(complex(0,1)*G**2*I27a11*I30a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a11*I32a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a11*I35a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a55*I8a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a11*I30a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a11*I32a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a11*I35a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a55*I8a11*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_913 = Coupling(name = 'GC_913',
                  value = '(complex(0,1)*G**2*I26a63*I27a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I26a11*I27a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a63*I29a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a11*I29a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a11*I30a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a11*I32a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a11*I35a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a63*I8a11)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a63*I27a11*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a11*I27a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a63*I29a11*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a11*I29a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a11*I30a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a11*I32a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a11*I35a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a63*I8a11*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_914 = Coupling(name = 'GC_914',
                  value = '(complex(0,1)*G**2*I26a66*I27a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I26a11*I27a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a66*I29a11)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a11*I29a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a11*I30a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a11*I32a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a11*I35a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a66*I8a11)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a66*I27a11*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a11*I27a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a66*I29a11*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a11*I29a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a11*I30a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a11*I32a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a11*I35a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a66*I8a11*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_915 = Coupling(name = 'GC_915',
                  value = '(complex(0,1)*G**2*I26a33*I27a22)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I26a22*I27a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a33*I29a22)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a22*I29a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a22*I30a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a22*I32a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a22*I35a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a33*I8a22)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a33*I27a22*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a22*I27a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a33*I29a22*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a22*I29a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a22*I30a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a22*I32a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a22*I35a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a33*I8a22*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_916 = Coupling(name = 'GC_916',
                  value = '(complex(0,1)*G**2*I26a36*I27a22)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I26a22*I27a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a36*I29a22)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a22*I29a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a22*I30a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a22*I32a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a22*I35a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a36*I8a22)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a36*I27a22*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a22*I27a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a36*I29a22*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a22*I29a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a22*I30a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a22*I32a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a22*I35a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a36*I8a22*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_917 = Coupling(name = 'GC_917',
                  value = '-(complex(0,1)*G**2*I27a22*I30a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a22*I32a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a22*I35a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a44*I8a22)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a22*I30a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a22*I32a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a22*I35a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a44*I8a22*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_918 = Coupling(name = 'GC_918',
                  value = '-(complex(0,1)*G**2*I27a22*I30a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a22*I32a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a22*I35a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a55*I8a22)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a22*I30a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a22*I32a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a22*I35a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a55*I8a22*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_919 = Coupling(name = 'GC_919',
                  value = '(complex(0,1)*G**2*I26a63*I27a22)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I26a22*I27a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a63*I29a22)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a22*I29a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a22*I30a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a22*I32a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a22*I35a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a63*I8a22)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a63*I27a22*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a22*I27a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a63*I29a22*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a22*I29a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a22*I30a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a22*I32a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a22*I35a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a63*I8a22*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_920 = Coupling(name = 'GC_920',
                  value = '(complex(0,1)*G**2*I26a66*I27a22)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I26a22*I27a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a66*I29a22)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a22*I29a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a22*I30a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a22*I32a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a22*I35a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a66*I8a22)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a66*I27a22*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a22*I27a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a66*I29a22*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a22*I29a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a22*I30a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a22*I32a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a22*I35a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a66*I8a22*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_921 = Coupling(name = 'GC_921',
                  value = '(complex(0,1)*G**2*I26a33*I27a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a33*I30a33)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a33*I33a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a33*I8a33)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a33*I27a33*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a33*I30a33*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a33*I33a33*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a33*I8a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_922 = Coupling(name = 'GC_922',
                  value = '(complex(0,1)*G**2*I26a33*I27a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a36*I29a33)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a36*I30a33)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a33*I32a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a33*I33a36)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a36*I35a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a33*I36a36)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a36*I8a33)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a33*I27a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a36*I29a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a36*I30a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a33*I32a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a33*I33a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a36*I35a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a33*I36a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a36*I8a33*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_923 = Coupling(name = 'GC_923',
                  value = '-(complex(0,1)*G**2*I27a33*I30a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a33*I32a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a44*I33a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a33*I33a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a33*I35a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a44*I36a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a33*I36a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a44*I8a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a33*I30a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a33*I32a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a44*I33a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a33*I33a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a33*I35a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a44*I36a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a33*I36a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a44*I8a33*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_924 = Coupling(name = 'GC_924',
                  value = '-(complex(0,1)*G**2*I27a33*I30a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a33*I32a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a55*I33a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a33*I33a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a33*I35a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a55*I36a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a33*I36a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a55*I8a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a33*I30a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a33*I32a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a55*I33a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a33*I33a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a33*I35a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a55*I36a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a33*I36a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a55*I8a33*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_925 = Coupling(name = 'GC_925',
                  value = '(complex(0,1)*G**2*I26a33*I27a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a33*I29a63)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a63*I30a33)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a63*I32a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a33*I33a63)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a33*I35a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a63*I36a33)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a63*I8a33)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a33*I27a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a33*I29a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a63*I30a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a63*I32a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a33*I33a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a33*I35a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a63*I36a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a63*I8a33*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_926 = Coupling(name = 'GC_926',
                  value = '(complex(0,1)*G**2*I26a36*I27a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a33*I29a36)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a33*I30a36)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a36*I32a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a36*I33a33)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a33*I35a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a36*I36a33)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a33*I8a36)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a36*I27a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a33*I29a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a33*I30a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a36*I32a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a36*I33a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a33*I35a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a36*I36a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a33*I8a36*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_927 = Coupling(name = 'GC_927',
                  value = '(complex(0,1)*G**2*I26a36*I27a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a36*I30a36)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a36*I33a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a36*I8a36)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a36*I27a36*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a36*I30a36*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a36*I33a36*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a36*I8a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_928 = Coupling(name = 'GC_928',
                  value = '-(complex(0,1)*G**2*I27a36*I30a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a36*I32a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a44*I33a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a36*I33a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a36*I35a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a44*I36a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a36*I36a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a44*I8a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a36*I30a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a36*I32a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a44*I33a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a36*I33a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a36*I35a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a44*I36a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a36*I36a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a44*I8a36*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_929 = Coupling(name = 'GC_929',
                  value = '-(complex(0,1)*G**2*I27a36*I30a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a36*I32a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a55*I33a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a36*I33a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a36*I35a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a55*I36a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a36*I36a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a55*I8a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a36*I30a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a36*I32a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a55*I33a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a36*I33a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a36*I35a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a55*I36a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a36*I36a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a55*I8a36*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_930 = Coupling(name = 'GC_930',
                  value = '(complex(0,1)*G**2*I26a36*I27a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a36*I29a66)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a66*I30a36)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a66*I32a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a36*I33a66)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a36*I35a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a66*I36a36)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a66*I8a36)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a36*I27a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a36*I29a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a66*I30a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a66*I32a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a36*I33a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a36*I35a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a66*I36a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a66*I8a36*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_931 = Coupling(name = 'GC_931',
                  value = '(complex(0,1)*G**2*I26a63*I27a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a63*I29a33)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a33*I30a63)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a33*I32a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a63*I33a33)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a63*I35a33)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a33*I36a63)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a33*I8a63)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a63*I27a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a63*I29a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a33*I30a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a33*I32a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a63*I33a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a63*I35a33*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a33*I36a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a33*I8a63*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_932 = Coupling(name = 'GC_932',
                  value = '(complex(0,1)*G**2*I26a63*I27a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I26a36*I27a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a63*I29a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a36*I29a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a63*I30a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a36*I30a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a63*I32a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a36*I32a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a63*I33a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a36*I33a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a63*I35a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a36*I35a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a63*I36a36)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a36*I36a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a63*I8a36)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a36*I8a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a63*I27a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a36*I27a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a63*I29a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a36*I29a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a63*I30a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a36*I30a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a63*I32a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a36*I32a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a63*I33a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a36*I33a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a63*I35a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a36*I35a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a63*I36a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a36*I36a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a63*I8a36*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a36*I8a63*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_933 = Coupling(name = 'GC_933',
                  value = '-(complex(0,1)*G**2*I27a63*I30a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a63*I32a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a63*I33a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a44*I33a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a63*I35a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a63*I36a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a44*I36a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a44*I8a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a63*I30a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a63*I32a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a63*I33a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a44*I33a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a63*I35a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a63*I36a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a44*I36a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a44*I8a63*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_934 = Coupling(name = 'GC_934',
                  value = '-(complex(0,1)*G**2*I27a63*I30a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a63*I32a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a63*I33a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a55*I33a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a63*I35a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a63*I36a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a55*I36a63)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a55*I8a63)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a63*I30a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a63*I32a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a63*I33a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a55*I33a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a63*I35a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a63*I36a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a55*I36a63*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a55*I8a63*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_935 = Coupling(name = 'GC_935',
                  value = '(complex(0,1)*G**2*I26a63*I27a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a63*I30a63)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a63*I33a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a63*I8a63)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a63*I27a63*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a63*I30a63*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a63*I33a63*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a63*I8a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_936 = Coupling(name = 'GC_936',
                  value = '(complex(0,1)*G**2*I26a63*I27a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a66*I29a63)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a66*I30a63)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a63*I32a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a63*I33a66)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a66*I35a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a63*I36a66)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a66*I8a63)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a63*I27a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a66*I29a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a66*I30a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a63*I32a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a63*I33a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a66*I35a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a63*I36a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a66*I8a63*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_937 = Coupling(name = 'GC_937',
                  value = '(complex(0,1)*G**2*I26a66*I27a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I26a33*I27a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a66*I29a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a33*I29a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a66*I30a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a33*I30a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a66*I32a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a33*I32a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a66*I33a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a33*I33a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a66*I35a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a33*I35a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a66*I36a33)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a33*I36a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a66*I8a33)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a33*I8a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a66*I27a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a33*I27a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a66*I29a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a33*I29a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a66*I30a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a33*I30a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a66*I32a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a33*I32a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a66*I33a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a33*I33a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a66*I35a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a33*I35a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a66*I36a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a33*I36a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a66*I8a33*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a33*I8a66*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_938 = Coupling(name = 'GC_938',
                  value = '(complex(0,1)*G**2*I26a66*I27a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a66*I29a36)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a36*I30a66)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a36*I32a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a66*I33a36)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a66*I35a36)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a36*I36a66)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a36*I8a66)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a66*I27a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a66*I29a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a36*I30a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a36*I32a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a66*I33a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a66*I35a36*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a36*I36a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a36*I8a66*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_939 = Coupling(name = 'GC_939',
                  value = '-(complex(0,1)*G**2*I27a66*I30a44)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a66*I32a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a66*I33a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a44*I33a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a66*I35a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a66*I36a44)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a44*I36a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a44*I8a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a66*I30a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a66*I32a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a66*I33a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a44*I33a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a66*I35a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a66*I36a44*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a44*I36a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a44*I8a66*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_940 = Coupling(name = 'GC_940',
                  value = '-(complex(0,1)*G**2*I27a66*I30a55)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a66*I32a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a66*I33a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a55*I33a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a66*I35a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a66*I36a55)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a55*I36a66)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a55*I8a66)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a66*I30a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a66*I32a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a66*I33a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a55*I33a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a66*I35a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a66*I36a55*sw**2)/(4.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a55*I36a66*sw**2)/(4.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a55*I8a66*sw**2)/(4.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_941 = Coupling(name = 'GC_941',
                  value = '(complex(0,1)*G**2*I26a66*I27a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a63*I29a66)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a63*I30a66)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I31a66*I32a63)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a66*I33a63)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a63*I35a66)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I35a66*I36a63)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a63*I8a66)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a66*I27a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I28a63*I29a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a63*I30a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I31a66*I32a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a66*I33a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I28a63*I35a66*sw**2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I35a66*I36a63*sw**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a63*I8a66*sw**2)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_942 = Coupling(name = 'GC_942',
                  value = '(complex(0,1)*G**2*I26a66*I27a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I27a66*I30a66)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I30a66*I33a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I34a66*I8a66)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I26a66*I27a66*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I27a66*I30a66*sw**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*G**2*I30a66*I33a66*sw**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*G**2*I34a66*I8a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QCD':2})

GC_943 = Coupling(name = 'GC_943',
                  value = '(2*ee**2*complex(0,1)*I97a11)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97a11)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I97a11*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_944 = Coupling(name = 'GC_944',
                  value = '(2*ee**2*complex(0,1)*I97a22)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97a22)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I97a22*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_945 = Coupling(name = 'GC_945',
                  value = '(2*ee**2*complex(0,1)*I97a33)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97a33)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I97a33*sw**2)/((-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I98a33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_946 = Coupling(name = 'GC_946',
                  value = '(2*ee**2*complex(0,1)*I97a36)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97a36)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I97a36*sw**2)/((-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I98a36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_947 = Coupling(name = 'GC_947',
                  value = '(2*ee**2*complex(0,1)*I97a63)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97a63)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I97a63*sw**2)/((-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I98a63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_948 = Coupling(name = 'GC_948',
                  value = '(2*ee**2*complex(0,1)*I97a66)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97a66)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I97a66*sw**2)/((-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I98a66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_949 = Coupling(name = 'GC_949',
                  value = '(complex(0,1)*I10a33*NN1x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I10a33*NN1x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*Rd3x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_950 = Coupling(name = 'GC_950',
                  value = '(complex(0,1)*I10a36*NN1x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I10a36*NN1x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*Rd6x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_951 = Coupling(name = 'GC_951',
                  value = '(complex(0,1)*I49a33*NN1x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I49a33*NN1x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*Rl3x6*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_952 = Coupling(name = 'GC_952',
                  value = '(complex(0,1)*I49a36*NN1x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I49a36*NN1x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*Rl6x6*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_953 = Coupling(name = 'GC_953',
                  value = '(complex(0,1)*I131a33*NN1x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I131a33*NN1x4*sw**2)/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*NN1x1*Ru3x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_954 = Coupling(name = 'GC_954',
                  value = '(complex(0,1)*I131a36*NN1x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I131a36*NN1x4*sw**2)/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*NN1x1*Ru6x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_955 = Coupling(name = 'GC_955',
                  value = '(complex(0,1)*I10a33*NN2x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I10a33*NN2x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*Rd3x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_956 = Coupling(name = 'GC_956',
                  value = '(complex(0,1)*I10a36*NN2x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I10a36*NN2x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*Rd6x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_957 = Coupling(name = 'GC_957',
                  value = '(complex(0,1)*I49a33*NN2x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I49a33*NN2x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*Rl3x6*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_958 = Coupling(name = 'GC_958',
                  value = '(complex(0,1)*I49a36*NN2x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I49a36*NN2x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*Rl6x6*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_959 = Coupling(name = 'GC_959',
                  value = '(complex(0,1)*I131a33*NN2x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I131a33*NN2x4*sw**2)/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*NN2x1*Ru3x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_960 = Coupling(name = 'GC_960',
                  value = '(complex(0,1)*I131a36*NN2x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I131a36*NN2x4*sw**2)/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*NN2x1*Ru6x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_961 = Coupling(name = 'GC_961',
                  value = '(complex(0,1)*I10a33*NN3x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I10a33*NN3x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*Rd3x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_962 = Coupling(name = 'GC_962',
                  value = '(complex(0,1)*I10a36*NN3x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I10a36*NN3x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*Rd6x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_963 = Coupling(name = 'GC_963',
                  value = '(complex(0,1)*I49a33*NN3x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I49a33*NN3x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*Rl3x6*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_964 = Coupling(name = 'GC_964',
                  value = '(complex(0,1)*I49a36*NN3x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I49a36*NN3x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*Rl6x6*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_965 = Coupling(name = 'GC_965',
                  value = '(complex(0,1)*I131a33*NN3x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I131a33*NN3x4*sw**2)/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*NN3x1*Ru3x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_966 = Coupling(name = 'GC_966',
                  value = '(complex(0,1)*I131a36*NN3x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I131a36*NN3x4*sw**2)/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*NN3x1*Ru6x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_967 = Coupling(name = 'GC_967',
                  value = '(complex(0,1)*I10a33*NN4x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I10a33*NN4x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*Rd3x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_968 = Coupling(name = 'GC_968',
                  value = '(complex(0,1)*I10a36*NN4x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I10a36*NN4x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*Rd6x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_969 = Coupling(name = 'GC_969',
                  value = '(complex(0,1)*I49a33*NN4x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I49a33*NN4x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*Rl3x6*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_970 = Coupling(name = 'GC_970',
                  value = '(complex(0,1)*I49a36*NN4x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I49a36*NN4x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*Rl6x6*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_971 = Coupling(name = 'GC_971',
                  value = '(complex(0,1)*I131a33*NN4x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I131a33*NN4x4*sw**2)/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*NN4x1*Ru3x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_972 = Coupling(name = 'GC_972',
                  value = '(complex(0,1)*I131a36*NN4x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I131a36*NN4x4*sw**2)/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*NN4x1*Ru6x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_973 = Coupling(name = 'GC_973',
                  value = '-((ee*complex(0,1)*I187a11*UU1x1)/sw)',
                  order = {'QED':1})

GC_974 = Coupling(name = 'GC_974',
                  value = '-((ee*complex(0,1)*I187a22*UU1x1)/sw)',
                  order = {'QED':1})

GC_975 = Coupling(name = 'GC_975',
                  value = '-((ee*complex(0,1)*I190a11*UU1x1)/sw)',
                  order = {'QED':1})

GC_976 = Coupling(name = 'GC_976',
                  value = '-((ee*complex(0,1)*I190a22*UU1x1)/sw)',
                  order = {'QED':1})

GC_977 = Coupling(name = 'GC_977',
                  value = 'complex(0,1)*I104a33*UU1x2',
                  order = {'QED':1})

GC_978 = Coupling(name = 'GC_978',
                  value = 'complex(0,1)*I134a33*UU1x2',
                  order = {'QED':1})

GC_979 = Coupling(name = 'GC_979',
                  value = 'complex(0,1)*I134a36*UU1x2',
                  order = {'QED':1})

GC_980 = Coupling(name = 'GC_980',
                  value = '-((ee*complex(0,1)*I187a33*UU1x1)/sw) + complex(0,1)*I188a33*UU1x2',
                  order = {'QED':1})

GC_981 = Coupling(name = 'GC_981',
                  value = '-((ee*complex(0,1)*I187a36*UU1x1)/sw) + complex(0,1)*I188a36*UU1x2',
                  order = {'QED':1})

GC_982 = Coupling(name = 'GC_982',
                  value = '-((ee*complex(0,1)*I190a33*UU1x1)/sw) + complex(0,1)*I191a33*UU1x2',
                  order = {'QED':1})

GC_983 = Coupling(name = 'GC_983',
                  value = '-((ee*complex(0,1)*I190a36*UU1x1)/sw) + complex(0,1)*I191a36*UU1x2',
                  order = {'QED':1})

GC_984 = Coupling(name = 'GC_984',
                  value = '-((ee*complex(0,1)*I187a11*UU2x1)/sw)',
                  order = {'QED':1})

GC_985 = Coupling(name = 'GC_985',
                  value = '-((ee*complex(0,1)*I187a22*UU2x1)/sw)',
                  order = {'QED':1})

GC_986 = Coupling(name = 'GC_986',
                  value = '-((ee*complex(0,1)*I190a11*UU2x1)/sw)',
                  order = {'QED':1})

GC_987 = Coupling(name = 'GC_987',
                  value = '-((ee*complex(0,1)*I190a22*UU2x1)/sw)',
                  order = {'QED':1})

GC_988 = Coupling(name = 'GC_988',
                  value = 'complex(0,1)*I104a33*UU2x2',
                  order = {'QED':1})

GC_989 = Coupling(name = 'GC_989',
                  value = 'complex(0,1)*I134a33*UU2x2',
                  order = {'QED':1})

GC_990 = Coupling(name = 'GC_990',
                  value = 'complex(0,1)*I134a36*UU2x2',
                  order = {'QED':1})

GC_991 = Coupling(name = 'GC_991',
                  value = '-((ee*complex(0,1)*I187a33*UU2x1)/sw) + complex(0,1)*I188a33*UU2x2',
                  order = {'QED':1})

GC_992 = Coupling(name = 'GC_992',
                  value = '-((ee*complex(0,1)*I187a36*UU2x1)/sw) + complex(0,1)*I188a36*UU2x2',
                  order = {'QED':1})

GC_993 = Coupling(name = 'GC_993',
                  value = '-((ee*complex(0,1)*I190a33*UU2x1)/sw) + complex(0,1)*I191a33*UU2x2',
                  order = {'QED':1})

GC_994 = Coupling(name = 'GC_994',
                  value = '-((ee*complex(0,1)*I190a36*UU2x1)/sw) + complex(0,1)*I191a36*UU2x2',
                  order = {'QED':1})

GC_995 = Coupling(name = 'GC_995',
                  value = '-((ee*complex(0,1)*I192a11*VV1x1)/sw)',
                  order = {'QED':1})

GC_996 = Coupling(name = 'GC_996',
                  value = '-((ee*complex(0,1)*I192a22*VV1x1)/sw)',
                  order = {'QED':1})

GC_997 = Coupling(name = 'GC_997',
                  value = '-((ee*complex(0,1)*I192a33*VV1x1)/sw)',
                  order = {'QED':1})

GC_998 = Coupling(name = 'GC_998',
                  value = '-((ee*complex(0,1)*I194a11*VV1x1)/sw)',
                  order = {'QED':1})

GC_999 = Coupling(name = 'GC_999',
                  value = '-((ee*complex(0,1)*I194a22*VV1x1)/sw)',
                  order = {'QED':1})

GC_1000 = Coupling(name = 'GC_1000',
                   value = 'complex(0,1)*I13a33*VV1x2',
                   order = {'QED':1})

GC_1001 = Coupling(name = 'GC_1001',
                   value = 'complex(0,1)*I13a36*VV1x2',
                   order = {'QED':1})

GC_1002 = Coupling(name = 'GC_1002',
                   value = '-((ee*complex(0,1)*I194a33*VV1x1)/sw) + complex(0,1)*I195a33*VV1x2',
                   order = {'QED':1})

GC_1003 = Coupling(name = 'GC_1003',
                   value = '-((ee*complex(0,1)*I194a36*VV1x1)/sw) + complex(0,1)*I195a36*VV1x2',
                   order = {'QED':1})

GC_1004 = Coupling(name = 'GC_1004',
                   value = '-((ee*complex(0,1)*I192a11*VV2x1)/sw)',
                   order = {'QED':1})

GC_1005 = Coupling(name = 'GC_1005',
                   value = '-((ee*complex(0,1)*I192a22*VV2x1)/sw)',
                   order = {'QED':1})

GC_1006 = Coupling(name = 'GC_1006',
                   value = '-((ee*complex(0,1)*I192a33*VV2x1)/sw)',
                   order = {'QED':1})

GC_1007 = Coupling(name = 'GC_1007',
                   value = '-((ee*complex(0,1)*I194a11*VV2x1)/sw)',
                   order = {'QED':1})

GC_1008 = Coupling(name = 'GC_1008',
                   value = '-((ee*complex(0,1)*I194a22*VV2x1)/sw)',
                   order = {'QED':1})

GC_1009 = Coupling(name = 'GC_1009',
                   value = 'complex(0,1)*I13a33*VV2x2',
                   order = {'QED':1})

GC_1010 = Coupling(name = 'GC_1010',
                   value = 'complex(0,1)*I13a36*VV2x2',
                   order = {'QED':1})

GC_1011 = Coupling(name = 'GC_1011',
                   value = '-((ee*complex(0,1)*I194a33*VV2x1)/sw) + complex(0,1)*I195a33*VV2x2',
                   order = {'QED':1})

GC_1012 = Coupling(name = 'GC_1012',
                   value = '-((ee*complex(0,1)*I194a36*VV2x1)/sw) + complex(0,1)*I195a36*VV2x2',
                   order = {'QED':1})

GC_1013 = Coupling(name = 'GC_1013',
                   value = '(ee*complex(0,1)*complexconjugate(CKM1x1))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1014 = Coupling(name = 'GC_1014',
                   value = '(ee*complex(0,1)*complexconjugate(CKM2x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1015 = Coupling(name = 'GC_1015',
                   value = '(ee*complex(0,1)*complexconjugate(CKM3x3))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1016 = Coupling(name = 'GC_1016',
                   value = '(cw*ee*complex(0,1)*Rd1x1*complexconjugate(NN1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd1x1*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd1x1*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1017 = Coupling(name = 'GC_1017',
                   value = '(cw*ee*complex(0,1)*Rd2x2*complexconjugate(NN1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd2x2*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd2x2*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1018 = Coupling(name = 'GC_1018',
                   value = '-((cw*ee*complex(0,1)*Rl1x1*complexconjugate(NN1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*Rl1x1*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl1x1*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1019 = Coupling(name = 'GC_1019',
                   value = '-((cw*ee*complex(0,1)*Rl2x2*complexconjugate(NN1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*Rl2x2*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl2x2*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1020 = Coupling(name = 'GC_1020',
                   value = '-((cw*ee*complex(0,1)*Rn1x1*complexconjugate(NN1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn1x1*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn1x1*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1021 = Coupling(name = 'GC_1021',
                   value = '-((cw*ee*complex(0,1)*Rn2x2*complexconjugate(NN1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn2x2*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn2x2*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1022 = Coupling(name = 'GC_1022',
                   value = '-((cw*ee*complex(0,1)*Rn3x3*complexconjugate(NN1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn3x3*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn3x3*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1023 = Coupling(name = 'GC_1023',
                   value = '(cw*ee*complex(0,1)*Ru1x1*complexconjugate(NN1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru1x1*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru1x1*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1024 = Coupling(name = 'GC_1024',
                   value = '(cw*ee*complex(0,1)*Ru2x2*complexconjugate(NN1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru2x2*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru2x2*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1025 = Coupling(name = 'GC_1025',
                   value = '(complex(0,1)*I11a33*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I11a33*sw**2*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Rd3x3*complexconjugate(NN1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd3x3*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd3x3*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1026 = Coupling(name = 'GC_1026',
                   value = '(complex(0,1)*I11a36*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I11a36*sw**2*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Rd6x3*complexconjugate(NN1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd6x3*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd6x3*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1027 = Coupling(name = 'GC_1027',
                   value = '(complex(0,1)*I50a33*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I50a33*sw**2*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*Rl3x3*complexconjugate(NN1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rl3x3*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl3x3*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1028 = Coupling(name = 'GC_1028',
                   value = '(complex(0,1)*I50a36*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I50a36*sw**2*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*Rl6x3*complexconjugate(NN1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rl6x3*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl6x3*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1029 = Coupling(name = 'GC_1029',
                   value = '-((ee*complex(0,1)*UU1x1*complexconjugate(NN1x2))/sw) - (ee*complex(0,1)*UU1x2*complexconjugate(NN1x3))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1030 = Coupling(name = 'GC_1030',
                   value = '-((ee*complex(0,1)*UU2x1*complexconjugate(NN1x2))/sw) - (ee*complex(0,1)*UU2x2*complexconjugate(NN1x3))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1031 = Coupling(name = 'GC_1031',
                   value = '-(cw*ee*complex(0,1)*NN1x3*complexconjugate(NN1x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*complexconjugate(NN1x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                   order = {'QED':1})

GC_1032 = Coupling(name = 'GC_1032',
                   value = '(cw*ee*complex(0,1)*NN2x3*complexconjugate(NN1x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*NN2x4*complexconjugate(NN1x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                   order = {'QED':1})

GC_1033 = Coupling(name = 'GC_1033',
                   value = '(cw*ee*complex(0,1)*NN3x3*complexconjugate(NN1x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*NN3x4*complexconjugate(NN1x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                   order = {'QED':1})

GC_1034 = Coupling(name = 'GC_1034',
                   value = '(cw*ee*complex(0,1)*NN4x3*complexconjugate(NN1x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*NN4x4*complexconjugate(NN1x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                   order = {'QED':1})

GC_1035 = Coupling(name = 'GC_1035',
                   value = '(complex(0,1)*I132a33*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I132a33*sw**2*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Ru3x3*complexconjugate(NN1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru3x3*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru3x3*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1036 = Coupling(name = 'GC_1036',
                   value = '(complex(0,1)*I132a36*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I132a36*sw**2*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Ru6x3*complexconjugate(NN1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru6x3*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru6x3*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1037 = Coupling(name = 'GC_1037',
                   value = '-((ee*complex(0,1)*VV1x1*complexconjugate(NN1x2))/sw) + (ee*complex(0,1)*VV1x2*complexconjugate(NN1x4))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1038 = Coupling(name = 'GC_1038',
                   value = '-((ee*complex(0,1)*VV2x1*complexconjugate(NN1x2))/sw) + (ee*complex(0,1)*VV2x2*complexconjugate(NN1x4))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1039 = Coupling(name = 'GC_1039',
                   value = '(cw*ee*complex(0,1)*Rd1x1*complexconjugate(NN2x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd1x1*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd1x1*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1040 = Coupling(name = 'GC_1040',
                   value = '(cw*ee*complex(0,1)*Rd2x2*complexconjugate(NN2x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd2x2*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd2x2*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1041 = Coupling(name = 'GC_1041',
                   value = '-((cw*ee*complex(0,1)*Rl1x1*complexconjugate(NN2x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*Rl1x1*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl1x1*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1042 = Coupling(name = 'GC_1042',
                   value = '-((cw*ee*complex(0,1)*Rl2x2*complexconjugate(NN2x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*Rl2x2*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl2x2*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1043 = Coupling(name = 'GC_1043',
                   value = '-((cw*ee*complex(0,1)*Rn1x1*complexconjugate(NN2x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn1x1*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn1x1*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1044 = Coupling(name = 'GC_1044',
                   value = '-((cw*ee*complex(0,1)*Rn2x2*complexconjugate(NN2x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn2x2*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn2x2*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1045 = Coupling(name = 'GC_1045',
                   value = '-((cw*ee*complex(0,1)*Rn3x3*complexconjugate(NN2x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn3x3*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn3x3*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1046 = Coupling(name = 'GC_1046',
                   value = '(cw*ee*complex(0,1)*Ru1x1*complexconjugate(NN2x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru1x1*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru1x1*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1047 = Coupling(name = 'GC_1047',
                   value = '(cw*ee*complex(0,1)*Ru2x2*complexconjugate(NN2x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru2x2*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru2x2*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1048 = Coupling(name = 'GC_1048',
                   value = '(complex(0,1)*I11a33*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I11a33*sw**2*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Rd3x3*complexconjugate(NN2x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd3x3*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd3x3*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1049 = Coupling(name = 'GC_1049',
                   value = '(complex(0,1)*I11a36*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I11a36*sw**2*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Rd6x3*complexconjugate(NN2x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd6x3*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd6x3*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1050 = Coupling(name = 'GC_1050',
                   value = '(complex(0,1)*I50a33*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I50a33*sw**2*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*Rl3x3*complexconjugate(NN2x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rl3x3*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl3x3*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1051 = Coupling(name = 'GC_1051',
                   value = '(complex(0,1)*I50a36*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I50a36*sw**2*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*Rl6x3*complexconjugate(NN2x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rl6x3*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl6x3*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1052 = Coupling(name = 'GC_1052',
                   value = '-((ee*complex(0,1)*UU1x1*complexconjugate(NN2x2))/sw) - (ee*complex(0,1)*UU1x2*complexconjugate(NN2x3))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1053 = Coupling(name = 'GC_1053',
                   value = '-((ee*complex(0,1)*UU2x1*complexconjugate(NN2x2))/sw) - (ee*complex(0,1)*UU2x2*complexconjugate(NN2x3))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1054 = Coupling(name = 'GC_1054',
                   value = '-(cw*ee*complex(0,1)*NN1x3*complexconjugate(NN2x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*complexconjugate(NN2x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                   order = {'QED':1})

GC_1055 = Coupling(name = 'GC_1055',
                   value = '-(cw*ee*complex(0,1)*NN2x3*complexconjugate(NN2x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN2x4*complexconjugate(NN2x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                   order = {'QED':1})

GC_1056 = Coupling(name = 'GC_1056',
                   value = '(cw*ee*complex(0,1)*NN3x3*complexconjugate(NN2x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*NN3x4*complexconjugate(NN2x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                   order = {'QED':1})

GC_1057 = Coupling(name = 'GC_1057',
                   value = '(cw*ee*complex(0,1)*NN4x3*complexconjugate(NN2x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*NN4x4*complexconjugate(NN2x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                   order = {'QED':1})

GC_1058 = Coupling(name = 'GC_1058',
                   value = '(complex(0,1)*I132a33*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I132a33*sw**2*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Ru3x3*complexconjugate(NN2x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru3x3*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru3x3*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1059 = Coupling(name = 'GC_1059',
                   value = '(complex(0,1)*I132a36*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I132a36*sw**2*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Ru6x3*complexconjugate(NN2x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru6x3*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru6x3*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1060 = Coupling(name = 'GC_1060',
                   value = '-((ee*complex(0,1)*VV1x1*complexconjugate(NN2x2))/sw) + (ee*complex(0,1)*VV1x2*complexconjugate(NN2x4))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1061 = Coupling(name = 'GC_1061',
                   value = '-((ee*complex(0,1)*VV2x1*complexconjugate(NN2x2))/sw) + (ee*complex(0,1)*VV2x2*complexconjugate(NN2x4))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1062 = Coupling(name = 'GC_1062',
                   value = '(cw*ee*complex(0,1)*Rd1x1*complexconjugate(NN3x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd1x1*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd1x1*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1063 = Coupling(name = 'GC_1063',
                   value = '(cw*ee*complex(0,1)*Rd2x2*complexconjugate(NN3x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd2x2*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd2x2*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1064 = Coupling(name = 'GC_1064',
                   value = '-((cw*ee*complex(0,1)*Rl1x1*complexconjugate(NN3x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*Rl1x1*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl1x1*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1065 = Coupling(name = 'GC_1065',
                   value = '-((cw*ee*complex(0,1)*Rl2x2*complexconjugate(NN3x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*Rl2x2*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl2x2*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1066 = Coupling(name = 'GC_1066',
                   value = '-((cw*ee*complex(0,1)*Rn1x1*complexconjugate(NN3x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn1x1*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn1x1*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1067 = Coupling(name = 'GC_1067',
                   value = '-((cw*ee*complex(0,1)*Rn2x2*complexconjugate(NN3x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn2x2*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn2x2*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1068 = Coupling(name = 'GC_1068',
                   value = '-((cw*ee*complex(0,1)*Rn3x3*complexconjugate(NN3x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn3x3*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn3x3*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1069 = Coupling(name = 'GC_1069',
                   value = '(cw*ee*complex(0,1)*Ru1x1*complexconjugate(NN3x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru1x1*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru1x1*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1070 = Coupling(name = 'GC_1070',
                   value = '(cw*ee*complex(0,1)*Ru2x2*complexconjugate(NN3x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru2x2*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru2x2*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1071 = Coupling(name = 'GC_1071',
                   value = '(complex(0,1)*I11a33*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I11a33*sw**2*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Rd3x3*complexconjugate(NN3x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd3x3*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd3x3*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1072 = Coupling(name = 'GC_1072',
                   value = '(complex(0,1)*I11a36*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I11a36*sw**2*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Rd6x3*complexconjugate(NN3x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd6x3*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd6x3*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1073 = Coupling(name = 'GC_1073',
                   value = '(complex(0,1)*I50a33*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I50a33*sw**2*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*Rl3x3*complexconjugate(NN3x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rl3x3*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl3x3*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1074 = Coupling(name = 'GC_1074',
                   value = '(complex(0,1)*I50a36*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I50a36*sw**2*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*Rl6x3*complexconjugate(NN3x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rl6x3*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl6x3*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1075 = Coupling(name = 'GC_1075',
                   value = '-((ee*complex(0,1)*UU1x1*complexconjugate(NN3x2))/sw) - (ee*complex(0,1)*UU1x2*complexconjugate(NN3x3))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1076 = Coupling(name = 'GC_1076',
                   value = '-((ee*complex(0,1)*UU2x1*complexconjugate(NN3x2))/sw) - (ee*complex(0,1)*UU2x2*complexconjugate(NN3x3))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1077 = Coupling(name = 'GC_1077',
                   value = '-(cw*ee*complex(0,1)*NN1x3*complexconjugate(NN3x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*complexconjugate(NN3x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                   order = {'QED':1})

GC_1078 = Coupling(name = 'GC_1078',
                   value = '-(cw*ee*complex(0,1)*NN2x3*complexconjugate(NN3x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN2x4*complexconjugate(NN3x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                   order = {'QED':1})

GC_1079 = Coupling(name = 'GC_1079',
                   value = '-(cw*ee*complex(0,1)*NN3x3*complexconjugate(NN3x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN3x4*complexconjugate(NN3x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                   order = {'QED':1})

GC_1080 = Coupling(name = 'GC_1080',
                   value = '(cw*ee*complex(0,1)*NN4x3*complexconjugate(NN3x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*NN4x4*complexconjugate(NN3x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                   order = {'QED':1})

GC_1081 = Coupling(name = 'GC_1081',
                   value = '(complex(0,1)*I132a33*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I132a33*sw**2*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Ru3x3*complexconjugate(NN3x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru3x3*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru3x3*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1082 = Coupling(name = 'GC_1082',
                   value = '(complex(0,1)*I132a36*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I132a36*sw**2*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Ru6x3*complexconjugate(NN3x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru6x3*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru6x3*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1083 = Coupling(name = 'GC_1083',
                   value = '-((ee*complex(0,1)*VV1x1*complexconjugate(NN3x2))/sw) + (ee*complex(0,1)*VV1x2*complexconjugate(NN3x4))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1084 = Coupling(name = 'GC_1084',
                   value = '-((ee*complex(0,1)*VV2x1*complexconjugate(NN3x2))/sw) + (ee*complex(0,1)*VV2x2*complexconjugate(NN3x4))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1085 = Coupling(name = 'GC_1085',
                   value = '(cw*ee*complex(0,1)*Rd1x1*complexconjugate(NN4x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd1x1*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd1x1*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1086 = Coupling(name = 'GC_1086',
                   value = '(cw*ee*complex(0,1)*Rd2x2*complexconjugate(NN4x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd2x2*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd2x2*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1087 = Coupling(name = 'GC_1087',
                   value = '-((cw*ee*complex(0,1)*Rl1x1*complexconjugate(NN4x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*Rl1x1*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl1x1*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1088 = Coupling(name = 'GC_1088',
                   value = '-((cw*ee*complex(0,1)*Rl2x2*complexconjugate(NN4x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*Rl2x2*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl2x2*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1089 = Coupling(name = 'GC_1089',
                   value = '-((cw*ee*complex(0,1)*Rn1x1*complexconjugate(NN4x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn1x1*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn1x1*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1090 = Coupling(name = 'GC_1090',
                   value = '-((cw*ee*complex(0,1)*Rn2x2*complexconjugate(NN4x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn2x2*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn2x2*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1091 = Coupling(name = 'GC_1091',
                   value = '-((cw*ee*complex(0,1)*Rn3x3*complexconjugate(NN4x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn3x3*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn3x3*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1092 = Coupling(name = 'GC_1092',
                   value = '(cw*ee*complex(0,1)*Ru1x1*complexconjugate(NN4x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru1x1*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru1x1*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1093 = Coupling(name = 'GC_1093',
                   value = '(cw*ee*complex(0,1)*Ru2x2*complexconjugate(NN4x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru2x2*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru2x2*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1094 = Coupling(name = 'GC_1094',
                   value = '(complex(0,1)*I11a33*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I11a33*sw**2*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Rd3x3*complexconjugate(NN4x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd3x3*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd3x3*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1095 = Coupling(name = 'GC_1095',
                   value = '(complex(0,1)*I11a36*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I11a36*sw**2*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Rd6x3*complexconjugate(NN4x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd6x3*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd6x3*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1096 = Coupling(name = 'GC_1096',
                   value = '(complex(0,1)*I50a33*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I50a33*sw**2*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*Rl3x3*complexconjugate(NN4x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rl3x3*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl3x3*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1097 = Coupling(name = 'GC_1097',
                   value = '(complex(0,1)*I50a36*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I50a36*sw**2*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*Rl6x3*complexconjugate(NN4x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rl6x3*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl6x3*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1098 = Coupling(name = 'GC_1098',
                   value = '-((ee*complex(0,1)*UU1x1*complexconjugate(NN4x2))/sw) - (ee*complex(0,1)*UU1x2*complexconjugate(NN4x3))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1099 = Coupling(name = 'GC_1099',
                   value = '-((ee*complex(0,1)*UU2x1*complexconjugate(NN4x2))/sw) - (ee*complex(0,1)*UU2x2*complexconjugate(NN4x3))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1100 = Coupling(name = 'GC_1100',
                   value = '-(cw*ee*complex(0,1)*NN1x3*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                   order = {'QED':1})

GC_1101 = Coupling(name = 'GC_1101',
                   value = '-(cw*ee*complex(0,1)*NN2x3*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN2x4*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                   order = {'QED':1})

GC_1102 = Coupling(name = 'GC_1102',
                   value = '-(cw*ee*complex(0,1)*NN3x3*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN3x4*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                   order = {'QED':1})

GC_1103 = Coupling(name = 'GC_1103',
                   value = '-(cw*ee*complex(0,1)*NN4x3*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN4x4*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                   order = {'QED':1})

GC_1104 = Coupling(name = 'GC_1104',
                   value = '(complex(0,1)*I132a33*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I132a33*sw**2*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Ru3x3*complexconjugate(NN4x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru3x3*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru3x3*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1105 = Coupling(name = 'GC_1105',
                   value = '(complex(0,1)*I132a36*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I132a36*sw**2*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Ru6x3*complexconjugate(NN4x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru6x3*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru6x3*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1106 = Coupling(name = 'GC_1106',
                   value = '-((ee*complex(0,1)*VV1x1*complexconjugate(NN4x2))/sw) + (ee*complex(0,1)*VV1x2*complexconjugate(NN4x4))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1107 = Coupling(name = 'GC_1107',
                   value = '-((ee*complex(0,1)*VV2x1*complexconjugate(NN4x2))/sw) + (ee*complex(0,1)*VV2x2*complexconjugate(NN4x4))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1108 = Coupling(name = 'GC_1108',
                   value = '-(complex(0,1)*G*complexconjugate(Rd1x1)*cmath.sqrt(2))',
                   order = {'QCD':1})

GC_1109 = Coupling(name = 'GC_1109',
                   value = '(cw*ee*complex(0,1)*NN1x1*complexconjugate(Rd1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*complexconjugate(Rd1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rd1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1110 = Coupling(name = 'GC_1110',
                   value = '(cw*ee*complex(0,1)*NN2x1*complexconjugate(Rd1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*complexconjugate(Rd1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rd1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1111 = Coupling(name = 'GC_1111',
                   value = '(cw*ee*complex(0,1)*NN3x1*complexconjugate(Rd1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*complexconjugate(Rd1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rd1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1112 = Coupling(name = 'GC_1112',
                   value = '(cw*ee*complex(0,1)*NN4x1*complexconjugate(Rd1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*complexconjugate(Rd1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rd1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1113 = Coupling(name = 'GC_1113',
                   value = '-(complex(0,1)*G*complexconjugate(Rd2x2)*cmath.sqrt(2))',
                   order = {'QCD':1})

GC_1114 = Coupling(name = 'GC_1114',
                   value = '(cw*ee*complex(0,1)*NN1x1*complexconjugate(Rd2x2))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*complexconjugate(Rd2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rd2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1115 = Coupling(name = 'GC_1115',
                   value = '(cw*ee*complex(0,1)*NN2x1*complexconjugate(Rd2x2))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*complexconjugate(Rd2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rd2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1116 = Coupling(name = 'GC_1116',
                   value = '(cw*ee*complex(0,1)*NN3x1*complexconjugate(Rd2x2))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*complexconjugate(Rd2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rd2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1117 = Coupling(name = 'GC_1117',
                   value = '(cw*ee*complex(0,1)*NN4x1*complexconjugate(Rd2x2))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*complexconjugate(Rd2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rd2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1118 = Coupling(name = 'GC_1118',
                   value = '-(complex(0,1)*G*complexconjugate(Rd3x3)*cmath.sqrt(2))',
                   order = {'QCD':1})

GC_1119 = Coupling(name = 'GC_1119',
                   value = '(complex(0,1)*I6a33*NN1x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I6a33*NN1x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*complexconjugate(Rd3x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*complexconjugate(Rd3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rd3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1120 = Coupling(name = 'GC_1120',
                   value = '(complex(0,1)*I6a33*NN2x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I6a33*NN2x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*complexconjugate(Rd3x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*complexconjugate(Rd3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rd3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1121 = Coupling(name = 'GC_1121',
                   value = '(complex(0,1)*I6a33*NN3x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I6a33*NN3x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*complexconjugate(Rd3x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*complexconjugate(Rd3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rd3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1122 = Coupling(name = 'GC_1122',
                   value = '(complex(0,1)*I6a33*NN4x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I6a33*NN4x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*complexconjugate(Rd3x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*complexconjugate(Rd3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rd3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1123 = Coupling(name = 'GC_1123',
                   value = 'complex(0,1)*G*complexconjugate(Rd3x6)*cmath.sqrt(2)',
                   order = {'QCD':1})

GC_1124 = Coupling(name = 'GC_1124',
                   value = '(complex(0,1)*I7a33*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I7a33*sw**2*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Rd3x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1125 = Coupling(name = 'GC_1125',
                   value = '(complex(0,1)*I7a33*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I7a33*sw**2*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Rd3x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1126 = Coupling(name = 'GC_1126',
                   value = '(complex(0,1)*I7a33*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I7a33*sw**2*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Rd3x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1127 = Coupling(name = 'GC_1127',
                   value = '(complex(0,1)*I7a33*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I7a33*sw**2*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Rd3x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1128 = Coupling(name = 'GC_1128',
                   value = 'complex(0,1)*G*complexconjugate(Rd4x4)*cmath.sqrt(2)',
                   order = {'QCD':1})

GC_1129 = Coupling(name = 'GC_1129',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Rd4x4)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1130 = Coupling(name = 'GC_1130',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Rd4x4)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1131 = Coupling(name = 'GC_1131',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Rd4x4)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1132 = Coupling(name = 'GC_1132',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Rd4x4)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1133 = Coupling(name = 'GC_1133',
                   value = 'complex(0,1)*G*complexconjugate(Rd5x5)*cmath.sqrt(2)',
                   order = {'QCD':1})

GC_1134 = Coupling(name = 'GC_1134',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Rd5x5)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1135 = Coupling(name = 'GC_1135',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Rd5x5)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1136 = Coupling(name = 'GC_1136',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Rd5x5)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1137 = Coupling(name = 'GC_1137',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Rd5x5)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1138 = Coupling(name = 'GC_1138',
                   value = '-(complex(0,1)*G*complexconjugate(Rd6x3)*cmath.sqrt(2))',
                   order = {'QCD':1})

GC_1139 = Coupling(name = 'GC_1139',
                   value = '(complex(0,1)*I6a36*NN1x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I6a36*NN1x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*complexconjugate(Rd6x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*complexconjugate(Rd6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rd6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1140 = Coupling(name = 'GC_1140',
                   value = '(complex(0,1)*I6a36*NN2x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I6a36*NN2x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*complexconjugate(Rd6x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*complexconjugate(Rd6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rd6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1141 = Coupling(name = 'GC_1141',
                   value = '(complex(0,1)*I6a36*NN3x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I6a36*NN3x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*complexconjugate(Rd6x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*complexconjugate(Rd6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rd6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1142 = Coupling(name = 'GC_1142',
                   value = '(complex(0,1)*I6a36*NN4x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I6a36*NN4x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*complexconjugate(Rd6x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*complexconjugate(Rd6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rd6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1143 = Coupling(name = 'GC_1143',
                   value = 'complex(0,1)*G*complexconjugate(Rd6x6)*cmath.sqrt(2)',
                   order = {'QCD':1})

GC_1144 = Coupling(name = 'GC_1144',
                   value = '(complex(0,1)*I7a36*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I7a36*sw**2*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Rd6x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1145 = Coupling(name = 'GC_1145',
                   value = '(complex(0,1)*I7a36*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I7a36*sw**2*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Rd6x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1146 = Coupling(name = 'GC_1146',
                   value = '(complex(0,1)*I7a36*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I7a36*sw**2*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Rd6x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1147 = Coupling(name = 'GC_1147',
                   value = '(complex(0,1)*I7a36*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I7a36*sw**2*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Rd6x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1148 = Coupling(name = 'GC_1148',
                   value = '-((cw*ee*complex(0,1)*NN1x1*complexconjugate(Rl1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*NN1x2*complexconjugate(Rl1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rl1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1149 = Coupling(name = 'GC_1149',
                   value = '-((cw*ee*complex(0,1)*NN2x1*complexconjugate(Rl1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*NN2x2*complexconjugate(Rl1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rl1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1150 = Coupling(name = 'GC_1150',
                   value = '-((cw*ee*complex(0,1)*NN3x1*complexconjugate(Rl1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*NN3x2*complexconjugate(Rl1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rl1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1151 = Coupling(name = 'GC_1151',
                   value = '-((cw*ee*complex(0,1)*NN4x1*complexconjugate(Rl1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*NN4x2*complexconjugate(Rl1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rl1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1152 = Coupling(name = 'GC_1152',
                   value = '-((cw*ee*complex(0,1)*NN1x1*complexconjugate(Rl2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*NN1x2*complexconjugate(Rl2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rl2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1153 = Coupling(name = 'GC_1153',
                   value = '-((cw*ee*complex(0,1)*NN2x1*complexconjugate(Rl2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*NN2x2*complexconjugate(Rl2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rl2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1154 = Coupling(name = 'GC_1154',
                   value = '-((cw*ee*complex(0,1)*NN3x1*complexconjugate(Rl2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*NN3x2*complexconjugate(Rl2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rl2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1155 = Coupling(name = 'GC_1155',
                   value = '-((cw*ee*complex(0,1)*NN4x1*complexconjugate(Rl2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*NN4x2*complexconjugate(Rl2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rl2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1156 = Coupling(name = 'GC_1156',
                   value = '(complex(0,1)*I45a33*NN1x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I45a33*NN1x3*sw**2)/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*complexconjugate(Rl3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*complexconjugate(Rl3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rl3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1157 = Coupling(name = 'GC_1157',
                   value = '(complex(0,1)*I45a33*NN2x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I45a33*NN2x3*sw**2)/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN2x1*complexconjugate(Rl3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*complexconjugate(Rl3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rl3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1158 = Coupling(name = 'GC_1158',
                   value = '(complex(0,1)*I45a33*NN3x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I45a33*NN3x3*sw**2)/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN3x1*complexconjugate(Rl3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*complexconjugate(Rl3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rl3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1159 = Coupling(name = 'GC_1159',
                   value = '(complex(0,1)*I45a33*NN4x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I45a33*NN4x3*sw**2)/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN4x1*complexconjugate(Rl3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*complexconjugate(Rl3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rl3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1160 = Coupling(name = 'GC_1160',
                   value = '(complex(0,1)*I46a33*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I46a33*sw**2*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Rl3x6)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1161 = Coupling(name = 'GC_1161',
                   value = '(complex(0,1)*I46a33*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I46a33*sw**2*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Rl3x6)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1162 = Coupling(name = 'GC_1162',
                   value = '(complex(0,1)*I46a33*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I46a33*sw**2*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Rl3x6)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1163 = Coupling(name = 'GC_1163',
                   value = '(complex(0,1)*I46a33*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I46a33*sw**2*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Rl3x6)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1164 = Coupling(name = 'GC_1164',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Rl4x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1165 = Coupling(name = 'GC_1165',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Rl4x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1166 = Coupling(name = 'GC_1166',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Rl4x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1167 = Coupling(name = 'GC_1167',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Rl4x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1168 = Coupling(name = 'GC_1168',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Rl5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1169 = Coupling(name = 'GC_1169',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Rl5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1170 = Coupling(name = 'GC_1170',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Rl5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1171 = Coupling(name = 'GC_1171',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Rl5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1172 = Coupling(name = 'GC_1172',
                   value = '(complex(0,1)*I45a36*NN1x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I45a36*NN1x3*sw**2)/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*complexconjugate(Rl6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*complexconjugate(Rl6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rl6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1173 = Coupling(name = 'GC_1173',
                   value = '(complex(0,1)*I45a36*NN2x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I45a36*NN2x3*sw**2)/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN2x1*complexconjugate(Rl6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*complexconjugate(Rl6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rl6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1174 = Coupling(name = 'GC_1174',
                   value = '(complex(0,1)*I45a36*NN3x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I45a36*NN3x3*sw**2)/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN3x1*complexconjugate(Rl6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*complexconjugate(Rl6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rl6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1175 = Coupling(name = 'GC_1175',
                   value = '(complex(0,1)*I45a36*NN4x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I45a36*NN4x3*sw**2)/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN4x1*complexconjugate(Rl6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*complexconjugate(Rl6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rl6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1176 = Coupling(name = 'GC_1176',
                   value = '(complex(0,1)*I46a36*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I46a36*sw**2*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Rl6x6)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1177 = Coupling(name = 'GC_1177',
                   value = '(complex(0,1)*I46a36*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I46a36*sw**2*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Rl6x6)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1178 = Coupling(name = 'GC_1178',
                   value = '(complex(0,1)*I46a36*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I46a36*sw**2*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Rl6x6)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1179 = Coupling(name = 'GC_1179',
                   value = '(complex(0,1)*I46a36*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I46a36*sw**2*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Rl6x6)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1180 = Coupling(name = 'GC_1180',
                   value = '-((cw*ee*complex(0,1)*NN1x1*complexconjugate(Rn1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN1x2*complexconjugate(Rn1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rn1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1181 = Coupling(name = 'GC_1181',
                   value = '-((cw*ee*complex(0,1)*NN2x1*complexconjugate(Rn1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN2x2*complexconjugate(Rn1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rn1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1182 = Coupling(name = 'GC_1182',
                   value = '-((cw*ee*complex(0,1)*NN3x1*complexconjugate(Rn1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN3x2*complexconjugate(Rn1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rn1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1183 = Coupling(name = 'GC_1183',
                   value = '-((cw*ee*complex(0,1)*NN4x1*complexconjugate(Rn1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN4x2*complexconjugate(Rn1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rn1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1184 = Coupling(name = 'GC_1184',
                   value = '-((cw*ee*complex(0,1)*NN1x1*complexconjugate(Rn2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN1x2*complexconjugate(Rn2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rn2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1185 = Coupling(name = 'GC_1185',
                   value = '-((cw*ee*complex(0,1)*NN2x1*complexconjugate(Rn2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN2x2*complexconjugate(Rn2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rn2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1186 = Coupling(name = 'GC_1186',
                   value = '-((cw*ee*complex(0,1)*NN3x1*complexconjugate(Rn2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN3x2*complexconjugate(Rn2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rn2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1187 = Coupling(name = 'GC_1187',
                   value = '-((cw*ee*complex(0,1)*NN4x1*complexconjugate(Rn2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN4x2*complexconjugate(Rn2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rn2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1188 = Coupling(name = 'GC_1188',
                   value = '-((cw*ee*complex(0,1)*NN1x1*complexconjugate(Rn3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN1x2*complexconjugate(Rn3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rn3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1189 = Coupling(name = 'GC_1189',
                   value = '-((cw*ee*complex(0,1)*NN2x1*complexconjugate(Rn3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN2x2*complexconjugate(Rn3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rn3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1190 = Coupling(name = 'GC_1190',
                   value = '-((cw*ee*complex(0,1)*NN3x1*complexconjugate(Rn3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN3x2*complexconjugate(Rn3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rn3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1191 = Coupling(name = 'GC_1191',
                   value = '-((cw*ee*complex(0,1)*NN4x1*complexconjugate(Rn3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN4x2*complexconjugate(Rn3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rn3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1192 = Coupling(name = 'GC_1192',
                   value = '-(complex(0,1)*G*complexconjugate(Ru1x1)*cmath.sqrt(2))',
                   order = {'QCD':1})

GC_1193 = Coupling(name = 'GC_1193',
                   value = '(cw*ee*complex(0,1)*NN1x1*complexconjugate(Ru1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*complexconjugate(Ru1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*complexconjugate(Ru1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1194 = Coupling(name = 'GC_1194',
                   value = '(cw*ee*complex(0,1)*NN2x1*complexconjugate(Ru1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*complexconjugate(Ru1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*complexconjugate(Ru1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1195 = Coupling(name = 'GC_1195',
                   value = '(cw*ee*complex(0,1)*NN3x1*complexconjugate(Ru1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*complexconjugate(Ru1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*complexconjugate(Ru1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1196 = Coupling(name = 'GC_1196',
                   value = '(cw*ee*complex(0,1)*NN4x1*complexconjugate(Ru1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*complexconjugate(Ru1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*complexconjugate(Ru1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1197 = Coupling(name = 'GC_1197',
                   value = '-(complex(0,1)*G*complexconjugate(Ru2x2)*cmath.sqrt(2))',
                   order = {'QCD':1})

GC_1198 = Coupling(name = 'GC_1198',
                   value = '(cw*ee*complex(0,1)*NN1x1*complexconjugate(Ru2x2))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*complexconjugate(Ru2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*complexconjugate(Ru2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1199 = Coupling(name = 'GC_1199',
                   value = '(cw*ee*complex(0,1)*NN2x1*complexconjugate(Ru2x2))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*complexconjugate(Ru2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*complexconjugate(Ru2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1200 = Coupling(name = 'GC_1200',
                   value = '(cw*ee*complex(0,1)*NN3x1*complexconjugate(Ru2x2))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*complexconjugate(Ru2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*complexconjugate(Ru2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1201 = Coupling(name = 'GC_1201',
                   value = '(cw*ee*complex(0,1)*NN4x1*complexconjugate(Ru2x2))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*complexconjugate(Ru2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*complexconjugate(Ru2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1202 = Coupling(name = 'GC_1202',
                   value = '-(complex(0,1)*G*complexconjugate(Ru3x3)*cmath.sqrt(2))',
                   order = {'QCD':1})

GC_1203 = Coupling(name = 'GC_1203',
                   value = '(complex(0,1)*I111a33*NN1x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I111a33*NN1x4*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*complexconjugate(Ru3x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*complexconjugate(Ru3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*complexconjugate(Ru3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1204 = Coupling(name = 'GC_1204',
                   value = '(complex(0,1)*I111a33*NN2x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I111a33*NN2x4*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*complexconjugate(Ru3x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*complexconjugate(Ru3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*complexconjugate(Ru3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1205 = Coupling(name = 'GC_1205',
                   value = '(complex(0,1)*I111a33*NN3x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I111a33*NN3x4*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*complexconjugate(Ru3x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*complexconjugate(Ru3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*complexconjugate(Ru3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1206 = Coupling(name = 'GC_1206',
                   value = '(complex(0,1)*I111a33*NN4x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I111a33*NN4x4*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*complexconjugate(Ru3x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*complexconjugate(Ru3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*complexconjugate(Ru3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1207 = Coupling(name = 'GC_1207',
                   value = 'complex(0,1)*G*complexconjugate(Ru3x6)*cmath.sqrt(2)',
                   order = {'QCD':1})

GC_1208 = Coupling(name = 'GC_1208',
                   value = '(complex(0,1)*I112a33*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I112a33*sw**2*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Ru3x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1209 = Coupling(name = 'GC_1209',
                   value = '(complex(0,1)*I112a33*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I112a33*sw**2*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Ru3x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1210 = Coupling(name = 'GC_1210',
                   value = '(complex(0,1)*I112a33*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I112a33*sw**2*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Ru3x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1211 = Coupling(name = 'GC_1211',
                   value = '(complex(0,1)*I112a33*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I112a33*sw**2*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Ru3x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1212 = Coupling(name = 'GC_1212',
                   value = 'complex(0,1)*G*complexconjugate(Ru4x4)*cmath.sqrt(2)',
                   order = {'QCD':1})

GC_1213 = Coupling(name = 'GC_1213',
                   value = '(-2*cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Ru4x4)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1214 = Coupling(name = 'GC_1214',
                   value = '(-2*cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Ru4x4)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1215 = Coupling(name = 'GC_1215',
                   value = '(-2*cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Ru4x4)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1216 = Coupling(name = 'GC_1216',
                   value = '(-2*cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Ru4x4)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1217 = Coupling(name = 'GC_1217',
                   value = 'complex(0,1)*G*complexconjugate(Ru5x5)*cmath.sqrt(2)',
                   order = {'QCD':1})

GC_1218 = Coupling(name = 'GC_1218',
                   value = '(-2*cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Ru5x5)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1219 = Coupling(name = 'GC_1219',
                   value = '(-2*cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Ru5x5)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1220 = Coupling(name = 'GC_1220',
                   value = '(-2*cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Ru5x5)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1221 = Coupling(name = 'GC_1221',
                   value = '(-2*cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Ru5x5)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1222 = Coupling(name = 'GC_1222',
                   value = '-(complex(0,1)*G*complexconjugate(Ru6x3)*cmath.sqrt(2))',
                   order = {'QCD':1})

GC_1223 = Coupling(name = 'GC_1223',
                   value = '(complex(0,1)*I111a36*NN1x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I111a36*NN1x4*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*complexconjugate(Ru6x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*complexconjugate(Ru6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*complexconjugate(Ru6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1224 = Coupling(name = 'GC_1224',
                   value = '(complex(0,1)*I111a36*NN2x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I111a36*NN2x4*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*complexconjugate(Ru6x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*complexconjugate(Ru6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*complexconjugate(Ru6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1225 = Coupling(name = 'GC_1225',
                   value = '(complex(0,1)*I111a36*NN3x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I111a36*NN3x4*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*complexconjugate(Ru6x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*complexconjugate(Ru6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*complexconjugate(Ru6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1226 = Coupling(name = 'GC_1226',
                   value = '(complex(0,1)*I111a36*NN4x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I111a36*NN4x4*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*complexconjugate(Ru6x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*complexconjugate(Ru6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*complexconjugate(Ru6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1227 = Coupling(name = 'GC_1227',
                   value = 'complex(0,1)*G*complexconjugate(Ru6x6)*cmath.sqrt(2)',
                   order = {'QCD':1})

GC_1228 = Coupling(name = 'GC_1228',
                   value = '(complex(0,1)*I112a36*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I112a36*sw**2*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Ru6x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1229 = Coupling(name = 'GC_1229',
                   value = '(complex(0,1)*I112a36*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I112a36*sw**2*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Ru6x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1230 = Coupling(name = 'GC_1230',
                   value = '(complex(0,1)*I112a36*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I112a36*sw**2*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Ru6x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1231 = Coupling(name = 'GC_1231',
                   value = '(complex(0,1)*I112a36*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I112a36*sw**2*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Ru6x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1232 = Coupling(name = 'GC_1232',
                   value = '-((ee*complex(0,1)*I12a11*complexconjugate(UU1x1))/sw)',
                   order = {'QED':1})

GC_1233 = Coupling(name = 'GC_1233',
                   value = '-((ee*complex(0,1)*I12a22*complexconjugate(UU1x1))/sw)',
                   order = {'QED':1})

GC_1234 = Coupling(name = 'GC_1234',
                   value = '-((ee*complex(0,1)*I51a11*complexconjugate(UU1x1))/sw)',
                   order = {'QED':1})

GC_1235 = Coupling(name = 'GC_1235',
                   value = '-((ee*complex(0,1)*I51a22*complexconjugate(UU1x1))/sw)',
                   order = {'QED':1})

GC_1236 = Coupling(name = 'GC_1236',
                   value = 'complex(0,1)*I193a33*complexconjugate(UU1x2)',
                   order = {'QED':1})

GC_1237 = Coupling(name = 'GC_1237',
                   value = 'complex(0,1)*I196a33*complexconjugate(UU1x2)',
                   order = {'QED':1})

GC_1238 = Coupling(name = 'GC_1238',
                   value = 'complex(0,1)*I196a36*complexconjugate(UU1x2)',
                   order = {'QED':1})

GC_1239 = Coupling(name = 'GC_1239',
                   value = '-((ee*complex(0,1)*I12a33*complexconjugate(UU1x1))/sw) + complex(0,1)*I14a33*complexconjugate(UU1x2)',
                   order = {'QED':1})

GC_1240 = Coupling(name = 'GC_1240',
                   value = '-((ee*complex(0,1)*I12a36*complexconjugate(UU1x1))/sw) + complex(0,1)*I14a36*complexconjugate(UU1x2)',
                   order = {'QED':1})

GC_1241 = Coupling(name = 'GC_1241',
                   value = '-((ee*complex(0,1)*I51a33*complexconjugate(UU1x1))/sw) + complex(0,1)*I52a33*complexconjugate(UU1x2)',
                   order = {'QED':1})

GC_1242 = Coupling(name = 'GC_1242',
                   value = '-((ee*complex(0,1)*I51a36*complexconjugate(UU1x1))/sw) + complex(0,1)*I52a36*complexconjugate(UU1x2)',
                   order = {'QED':1})

GC_1243 = Coupling(name = 'GC_1243',
                   value = '-((ee*complex(0,1)*NN1x2*complexconjugate(UU1x1))/sw) - (ee*complex(0,1)*NN1x3*complexconjugate(UU1x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1244 = Coupling(name = 'GC_1244',
                   value = '-((ee*complex(0,1)*NN2x2*complexconjugate(UU1x1))/sw) - (ee*complex(0,1)*NN2x3*complexconjugate(UU1x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1245 = Coupling(name = 'GC_1245',
                   value = '-((ee*complex(0,1)*NN3x2*complexconjugate(UU1x1))/sw) - (ee*complex(0,1)*NN3x3*complexconjugate(UU1x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1246 = Coupling(name = 'GC_1246',
                   value = '-((ee*complex(0,1)*NN4x2*complexconjugate(UU1x1))/sw) - (ee*complex(0,1)*NN4x3*complexconjugate(UU1x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1247 = Coupling(name = 'GC_1247',
                   value = 'ee*complex(0,1)*UU1x1*complexconjugate(UU1x1) + ee*complex(0,1)*UU1x2*complexconjugate(UU1x2)',
                   order = {'QED':1})

GC_1248 = Coupling(name = 'GC_1248',
                   value = '-((cw*ee*complex(0,1)*UU1x1*complexconjugate(UU1x1))/((-1 + sw)*sw*(1 + sw))) + (cw*ee*complex(0,1)*sw*UU1x1*complexconjugate(UU1x1))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*UU1x2*complexconjugate(UU1x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*sw*UU1x2*complexconjugate(UU1x2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1249 = Coupling(name = 'GC_1249',
                   value = 'ee*complex(0,1)*UU2x1*complexconjugate(UU1x1) + ee*complex(0,1)*UU2x2*complexconjugate(UU1x2)',
                   order = {'QED':1})

GC_1250 = Coupling(name = 'GC_1250',
                   value = '-((cw*ee*complex(0,1)*UU2x1*complexconjugate(UU1x1))/((-1 + sw)*sw*(1 + sw))) + (cw*ee*complex(0,1)*sw*UU2x1*complexconjugate(UU1x1))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*UU2x2*complexconjugate(UU1x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*sw*UU2x2*complexconjugate(UU1x2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1251 = Coupling(name = 'GC_1251',
                   value = '-((ee*complex(0,1)*I12a11*complexconjugate(UU2x1))/sw)',
                   order = {'QED':1})

GC_1252 = Coupling(name = 'GC_1252',
                   value = '-((ee*complex(0,1)*I12a22*complexconjugate(UU2x1))/sw)',
                   order = {'QED':1})

GC_1253 = Coupling(name = 'GC_1253',
                   value = '-((ee*complex(0,1)*I51a11*complexconjugate(UU2x1))/sw)',
                   order = {'QED':1})

GC_1254 = Coupling(name = 'GC_1254',
                   value = '-((ee*complex(0,1)*I51a22*complexconjugate(UU2x1))/sw)',
                   order = {'QED':1})

GC_1255 = Coupling(name = 'GC_1255',
                   value = 'complex(0,1)*I193a33*complexconjugate(UU2x2)',
                   order = {'QED':1})

GC_1256 = Coupling(name = 'GC_1256',
                   value = 'complex(0,1)*I196a33*complexconjugate(UU2x2)',
                   order = {'QED':1})

GC_1257 = Coupling(name = 'GC_1257',
                   value = 'complex(0,1)*I196a36*complexconjugate(UU2x2)',
                   order = {'QED':1})

GC_1258 = Coupling(name = 'GC_1258',
                   value = '-((ee*complex(0,1)*I12a33*complexconjugate(UU2x1))/sw) + complex(0,1)*I14a33*complexconjugate(UU2x2)',
                   order = {'QED':1})

GC_1259 = Coupling(name = 'GC_1259',
                   value = '-((ee*complex(0,1)*I12a36*complexconjugate(UU2x1))/sw) + complex(0,1)*I14a36*complexconjugate(UU2x2)',
                   order = {'QED':1})

GC_1260 = Coupling(name = 'GC_1260',
                   value = '-((ee*complex(0,1)*I51a33*complexconjugate(UU2x1))/sw) + complex(0,1)*I52a33*complexconjugate(UU2x2)',
                   order = {'QED':1})

GC_1261 = Coupling(name = 'GC_1261',
                   value = '-((ee*complex(0,1)*I51a36*complexconjugate(UU2x1))/sw) + complex(0,1)*I52a36*complexconjugate(UU2x2)',
                   order = {'QED':1})

GC_1262 = Coupling(name = 'GC_1262',
                   value = '-((ee*complex(0,1)*NN1x2*complexconjugate(UU2x1))/sw) - (ee*complex(0,1)*NN1x3*complexconjugate(UU2x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1263 = Coupling(name = 'GC_1263',
                   value = '-((ee*complex(0,1)*NN2x2*complexconjugate(UU2x1))/sw) - (ee*complex(0,1)*NN2x3*complexconjugate(UU2x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1264 = Coupling(name = 'GC_1264',
                   value = '-((ee*complex(0,1)*NN3x2*complexconjugate(UU2x1))/sw) - (ee*complex(0,1)*NN3x3*complexconjugate(UU2x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1265 = Coupling(name = 'GC_1265',
                   value = '-((ee*complex(0,1)*NN4x2*complexconjugate(UU2x1))/sw) - (ee*complex(0,1)*NN4x3*complexconjugate(UU2x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1266 = Coupling(name = 'GC_1266',
                   value = 'ee*complex(0,1)*UU1x1*complexconjugate(UU2x1) + ee*complex(0,1)*UU1x2*complexconjugate(UU2x2)',
                   order = {'QED':1})

GC_1267 = Coupling(name = 'GC_1267',
                   value = '-((cw*ee*complex(0,1)*UU1x1*complexconjugate(UU2x1))/((-1 + sw)*sw*(1 + sw))) + (cw*ee*complex(0,1)*sw*UU1x1*complexconjugate(UU2x1))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*UU1x2*complexconjugate(UU2x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*sw*UU1x2*complexconjugate(UU2x2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1268 = Coupling(name = 'GC_1268',
                   value = 'ee*complex(0,1)*UU2x1*complexconjugate(UU2x1) + ee*complex(0,1)*UU2x2*complexconjugate(UU2x2)',
                   order = {'QED':1})

GC_1269 = Coupling(name = 'GC_1269',
                   value = '-((cw*ee*complex(0,1)*UU2x1*complexconjugate(UU2x1))/((-1 + sw)*sw*(1 + sw))) + (cw*ee*complex(0,1)*sw*UU2x1*complexconjugate(UU2x1))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*UU2x2*complexconjugate(UU2x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*sw*UU2x2*complexconjugate(UU2x2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1270 = Coupling(name = 'GC_1270',
                   value = '-((ee*complex(0,1)*I103a11*complexconjugate(VV1x1))/sw)',
                   order = {'QED':1})

GC_1271 = Coupling(name = 'GC_1271',
                   value = '-((ee*complex(0,1)*I103a22*complexconjugate(VV1x1))/sw)',
                   order = {'QED':1})

GC_1272 = Coupling(name = 'GC_1272',
                   value = '-((ee*complex(0,1)*I103a33*complexconjugate(VV1x1))/sw)',
                   order = {'QED':1})

GC_1273 = Coupling(name = 'GC_1273',
                   value = '-((ee*complex(0,1)*I133a11*complexconjugate(VV1x1))/sw)',
                   order = {'QED':1})

GC_1274 = Coupling(name = 'GC_1274',
                   value = '-((ee*complex(0,1)*I133a22*complexconjugate(VV1x1))/sw)',
                   order = {'QED':1})

GC_1275 = Coupling(name = 'GC_1275',
                   value = 'complex(0,1)*I189a33*complexconjugate(VV1x2)',
                   order = {'QED':1})

GC_1276 = Coupling(name = 'GC_1276',
                   value = 'complex(0,1)*I189a36*complexconjugate(VV1x2)',
                   order = {'QED':1})

GC_1277 = Coupling(name = 'GC_1277',
                   value = '-((ee*complex(0,1)*I133a33*complexconjugate(VV1x1))/sw) + complex(0,1)*I135a33*complexconjugate(VV1x2)',
                   order = {'QED':1})

GC_1278 = Coupling(name = 'GC_1278',
                   value = '-((ee*complex(0,1)*I133a36*complexconjugate(VV1x1))/sw) + complex(0,1)*I135a36*complexconjugate(VV1x2)',
                   order = {'QED':1})

GC_1279 = Coupling(name = 'GC_1279',
                   value = '-((ee*complex(0,1)*NN1x2*complexconjugate(VV1x1))/sw) + (ee*complex(0,1)*NN1x4*complexconjugate(VV1x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1280 = Coupling(name = 'GC_1280',
                   value = '-((ee*complex(0,1)*NN2x2*complexconjugate(VV1x1))/sw) + (ee*complex(0,1)*NN2x4*complexconjugate(VV1x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1281 = Coupling(name = 'GC_1281',
                   value = '-((ee*complex(0,1)*NN3x2*complexconjugate(VV1x1))/sw) + (ee*complex(0,1)*NN3x4*complexconjugate(VV1x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1282 = Coupling(name = 'GC_1282',
                   value = '-((ee*complex(0,1)*NN4x2*complexconjugate(VV1x1))/sw) + (ee*complex(0,1)*NN4x4*complexconjugate(VV1x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1283 = Coupling(name = 'GC_1283',
                   value = 'ee*complex(0,1)*VV1x1*complexconjugate(VV1x1) + ee*complex(0,1)*VV1x2*complexconjugate(VV1x2)',
                   order = {'QED':1})

GC_1284 = Coupling(name = 'GC_1284',
                   value = '-((cw*ee*complex(0,1)*VV1x1*complexconjugate(VV1x1))/((-1 + sw)*sw*(1 + sw))) + (cw*ee*complex(0,1)*sw*VV1x1*complexconjugate(VV1x1))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*VV1x2*complexconjugate(VV1x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*sw*VV1x2*complexconjugate(VV1x2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1285 = Coupling(name = 'GC_1285',
                   value = 'ee*complex(0,1)*VV2x1*complexconjugate(VV1x1) + ee*complex(0,1)*VV2x2*complexconjugate(VV1x2)',
                   order = {'QED':1})

GC_1286 = Coupling(name = 'GC_1286',
                   value = '-((cw*ee*complex(0,1)*VV2x1*complexconjugate(VV1x1))/((-1 + sw)*sw*(1 + sw))) + (cw*ee*complex(0,1)*sw*VV2x1*complexconjugate(VV1x1))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*VV2x2*complexconjugate(VV1x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*sw*VV2x2*complexconjugate(VV1x2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1287 = Coupling(name = 'GC_1287',
                   value = '-((ee*complex(0,1)*I103a11*complexconjugate(VV2x1))/sw)',
                   order = {'QED':1})

GC_1288 = Coupling(name = 'GC_1288',
                   value = '-((ee*complex(0,1)*I103a22*complexconjugate(VV2x1))/sw)',
                   order = {'QED':1})

GC_1289 = Coupling(name = 'GC_1289',
                   value = '-((ee*complex(0,1)*I103a33*complexconjugate(VV2x1))/sw)',
                   order = {'QED':1})

GC_1290 = Coupling(name = 'GC_1290',
                   value = '-((ee*complex(0,1)*I133a11*complexconjugate(VV2x1))/sw)',
                   order = {'QED':1})

GC_1291 = Coupling(name = 'GC_1291',
                   value = '-((ee*complex(0,1)*I133a22*complexconjugate(VV2x1))/sw)',
                   order = {'QED':1})

GC_1292 = Coupling(name = 'GC_1292',
                   value = 'complex(0,1)*I189a33*complexconjugate(VV2x2)',
                   order = {'QED':1})

GC_1293 = Coupling(name = 'GC_1293',
                   value = 'complex(0,1)*I189a36*complexconjugate(VV2x2)',
                   order = {'QED':1})

GC_1294 = Coupling(name = 'GC_1294',
                   value = '-((ee*complex(0,1)*I133a33*complexconjugate(VV2x1))/sw) + complex(0,1)*I135a33*complexconjugate(VV2x2)',
                   order = {'QED':1})

GC_1295 = Coupling(name = 'GC_1295',
                   value = '-((ee*complex(0,1)*I133a36*complexconjugate(VV2x1))/sw) + complex(0,1)*I135a36*complexconjugate(VV2x2)',
                   order = {'QED':1})

GC_1296 = Coupling(name = 'GC_1296',
                   value = '-((ee*complex(0,1)*NN1x2*complexconjugate(VV2x1))/sw) + (ee*complex(0,1)*NN1x4*complexconjugate(VV2x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1297 = Coupling(name = 'GC_1297',
                   value = '-((ee*complex(0,1)*NN2x2*complexconjugate(VV2x1))/sw) + (ee*complex(0,1)*NN2x4*complexconjugate(VV2x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1298 = Coupling(name = 'GC_1298',
                   value = '-((ee*complex(0,1)*NN3x2*complexconjugate(VV2x1))/sw) + (ee*complex(0,1)*NN3x4*complexconjugate(VV2x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1299 = Coupling(name = 'GC_1299',
                   value = '-((ee*complex(0,1)*NN4x2*complexconjugate(VV2x1))/sw) + (ee*complex(0,1)*NN4x4*complexconjugate(VV2x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1300 = Coupling(name = 'GC_1300',
                   value = 'ee*complex(0,1)*VV1x1*complexconjugate(VV2x1) + ee*complex(0,1)*VV1x2*complexconjugate(VV2x2)',
                   order = {'QED':1})

GC_1301 = Coupling(name = 'GC_1301',
                   value = '-((cw*ee*complex(0,1)*VV1x1*complexconjugate(VV2x1))/((-1 + sw)*sw*(1 + sw))) + (cw*ee*complex(0,1)*sw*VV1x1*complexconjugate(VV2x1))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*VV1x2*complexconjugate(VV2x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*sw*VV1x2*complexconjugate(VV2x2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1302 = Coupling(name = 'GC_1302',
                   value = 'ee*complex(0,1)*VV2x1*complexconjugate(VV2x1) + ee*complex(0,1)*VV2x2*complexconjugate(VV2x2)',
                   order = {'QED':1})

GC_1303 = Coupling(name = 'GC_1303',
                   value = '-((cw*ee*complex(0,1)*VV2x1*complexconjugate(VV2x1))/((-1 + sw)*sw*(1 + sw))) + (cw*ee*complex(0,1)*sw*VV2x1*complexconjugate(VV2x1))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*VV2x2*complexconjugate(VV2x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*sw*VV2x2*complexconjugate(VV2x2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1304 = Coupling(name = 'GC_1304',
                   value = '-((complex(0,1)*yd3x3*cmath.cos(alp))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1305 = Coupling(name = 'GC_1305',
                   value = '-((complex(0,1)*ye3x3*cmath.cos(alp))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1306 = Coupling(name = 'GC_1306',
                   value = '-((complex(0,1)*yu3x3*cmath.cos(alp))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1307 = Coupling(name = 'GC_1307',
                   value = '-((complex(0,1)*complexconjugate(yd3x3)*cmath.cos(alp))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1308 = Coupling(name = 'GC_1308',
                   value = '-((complex(0,1)*complexconjugate(ye3x3)*cmath.cos(alp))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1309 = Coupling(name = 'GC_1309',
                   value = '-((complex(0,1)*complexconjugate(yu3x3)*cmath.cos(alp))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1310 = Coupling(name = 'GC_1310',
                   value = 'complex(0,1)*I1a33*cmath.cos(beta)',
                   order = {'QED':1})

GC_1311 = Coupling(name = 'GC_1311',
                   value = '-(complex(0,1)*I202a33*cmath.cos(beta))',
                   order = {'QED':1})

GC_1312 = Coupling(name = 'GC_1312',
                   value = '-(complex(0,1)*I2a33*cmath.cos(beta))',
                   order = {'QED':1})

GC_1313 = Coupling(name = 'GC_1313',
                   value = '-(complex(0,1)*I3a33*cmath.cos(beta))',
                   order = {'QED':1})

GC_1314 = Coupling(name = 'GC_1314',
                   value = 'complex(0,1)*I4a33*cmath.cos(beta)',
                   order = {'QED':1})

GC_1315 = Coupling(name = 'GC_1315',
                   value = '-(complex(0,1)*I5a33*cmath.cos(beta))',
                   order = {'QED':1})

GC_1316 = Coupling(name = 'GC_1316',
                   value = '-((yd3x3*cmath.cos(beta))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1317 = Coupling(name = 'GC_1317',
                   value = '-((ye3x3*cmath.cos(beta))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1318 = Coupling(name = 'GC_1318',
                   value = '(yu3x3*cmath.cos(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1319 = Coupling(name = 'GC_1319',
                   value = '(complexconjugate(yd3x3)*cmath.cos(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1320 = Coupling(name = 'GC_1320',
                   value = '(complexconjugate(ye3x3)*cmath.cos(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1321 = Coupling(name = 'GC_1321',
                   value = '-((complexconjugate(yu3x3)*cmath.cos(beta))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1322 = Coupling(name = 'GC_1322',
                   value = '-((ee*complex(0,1)*NN1x3*UU1x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw))) + (ee*complex(0,1)*NN1x3*sw*UU1x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*UU1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*UU1x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*UU1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1323 = Coupling(name = 'GC_1323',
                   value = '-((ee*complex(0,1)*NN2x3*UU1x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw))) + (ee*complex(0,1)*NN2x3*sw*UU1x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*UU1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*UU1x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*UU1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1324 = Coupling(name = 'GC_1324',
                   value = '-((ee*complex(0,1)*NN3x3*UU1x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw))) + (ee*complex(0,1)*NN3x3*sw*UU1x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*UU1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*UU1x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*UU1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1325 = Coupling(name = 'GC_1325',
                   value = '-((ee*complex(0,1)*NN4x3*UU1x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw))) + (ee*complex(0,1)*NN4x3*sw*UU1x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*UU1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*UU1x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*UU1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1326 = Coupling(name = 'GC_1326',
                   value = '-((ee*complex(0,1)*NN1x3*UU2x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw))) + (ee*complex(0,1)*NN1x3*sw*UU2x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*UU2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*UU2x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*UU2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1327 = Coupling(name = 'GC_1327',
                   value = '-((ee*complex(0,1)*NN2x3*UU2x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw))) + (ee*complex(0,1)*NN2x3*sw*UU2x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*UU2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*UU2x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*UU2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1328 = Coupling(name = 'GC_1328',
                   value = '-((ee*complex(0,1)*NN3x3*UU2x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw))) + (ee*complex(0,1)*NN3x3*sw*UU2x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*UU2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*UU2x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*UU2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1329 = Coupling(name = 'GC_1329',
                   value = '-((ee*complex(0,1)*NN4x3*UU2x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw))) + (ee*complex(0,1)*NN4x3*sw*UU2x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*UU2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*UU2x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*UU2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1330 = Coupling(name = 'GC_1330',
                   value = '(ee*complex(0,1)*NN1x4*VV1x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x4*sw*VV1x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*VV1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*VV1x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*VV1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1331 = Coupling(name = 'GC_1331',
                   value = '(ee*complex(0,1)*NN2x4*VV1x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x4*sw*VV1x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*VV1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*VV1x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*VV1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1332 = Coupling(name = 'GC_1332',
                   value = '(ee*complex(0,1)*NN3x4*VV1x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x4*sw*VV1x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*VV1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*VV1x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*VV1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1333 = Coupling(name = 'GC_1333',
                   value = '(ee*complex(0,1)*NN4x4*VV1x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN4x4*sw*VV1x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*VV1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*VV1x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*VV1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1334 = Coupling(name = 'GC_1334',
                   value = '(ee*complex(0,1)*NN1x4*VV2x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x4*sw*VV2x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*VV2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*VV2x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*VV2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1335 = Coupling(name = 'GC_1335',
                   value = '(ee*complex(0,1)*NN2x4*VV2x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x4*sw*VV2x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*VV2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*VV2x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*VV2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1336 = Coupling(name = 'GC_1336',
                   value = '(ee*complex(0,1)*NN3x4*VV2x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x4*sw*VV2x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*VV2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*VV2x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*VV2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1337 = Coupling(name = 'GC_1337',
                   value = '(ee*complex(0,1)*NN4x4*VV2x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN4x4*sw*VV2x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*VV2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*VV2x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*VV2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1338 = Coupling(name = 'GC_1338',
                   value = '-((ee*complex(0,1)*complexconjugate(NN1x3)*complexconjugate(UU1x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw))) + (ee*complex(0,1)*sw*complexconjugate(NN1x3)*complexconjugate(UU1x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1339 = Coupling(name = 'GC_1339',
                   value = '-((ee*complex(0,1)*complexconjugate(NN2x3)*complexconjugate(UU1x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw))) + (ee*complex(0,1)*sw*complexconjugate(NN2x3)*complexconjugate(UU1x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1340 = Coupling(name = 'GC_1340',
                   value = '-((ee*complex(0,1)*complexconjugate(NN3x3)*complexconjugate(UU1x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw))) + (ee*complex(0,1)*sw*complexconjugate(NN3x3)*complexconjugate(UU1x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN3x2)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN3x2)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1341 = Coupling(name = 'GC_1341',
                   value = '-((ee*complex(0,1)*complexconjugate(NN4x3)*complexconjugate(UU1x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw))) + (ee*complex(0,1)*sw*complexconjugate(NN4x3)*complexconjugate(UU1x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN4x2)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN4x2)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1342 = Coupling(name = 'GC_1342',
                   value = '-((ee*complex(0,1)*complexconjugate(NN1x3)*complexconjugate(UU2x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw))) + (ee*complex(0,1)*sw*complexconjugate(NN1x3)*complexconjugate(UU2x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1343 = Coupling(name = 'GC_1343',
                   value = '-((ee*complex(0,1)*complexconjugate(NN2x3)*complexconjugate(UU2x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw))) + (ee*complex(0,1)*sw*complexconjugate(NN2x3)*complexconjugate(UU2x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1344 = Coupling(name = 'GC_1344',
                   value = '-((ee*complex(0,1)*complexconjugate(NN3x3)*complexconjugate(UU2x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw))) + (ee*complex(0,1)*sw*complexconjugate(NN3x3)*complexconjugate(UU2x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN3x2)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN3x2)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1345 = Coupling(name = 'GC_1345',
                   value = '-((ee*complex(0,1)*complexconjugate(NN4x3)*complexconjugate(UU2x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw))) + (ee*complex(0,1)*sw*complexconjugate(NN4x3)*complexconjugate(UU2x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN4x2)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN4x2)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1346 = Coupling(name = 'GC_1346',
                   value = '(ee*complex(0,1)*complexconjugate(NN1x4)*complexconjugate(VV1x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN1x4)*complexconjugate(VV1x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1347 = Coupling(name = 'GC_1347',
                   value = '(ee*complex(0,1)*complexconjugate(NN2x4)*complexconjugate(VV1x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN2x4)*complexconjugate(VV1x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1348 = Coupling(name = 'GC_1348',
                   value = '(ee*complex(0,1)*complexconjugate(NN3x4)*complexconjugate(VV1x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN3x4)*complexconjugate(VV1x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN3x2)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN3x2)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1349 = Coupling(name = 'GC_1349',
                   value = '(ee*complex(0,1)*complexconjugate(NN4x4)*complexconjugate(VV1x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN4x4)*complexconjugate(VV1x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN4x2)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN4x2)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1350 = Coupling(name = 'GC_1350',
                   value = '(ee*complex(0,1)*complexconjugate(NN1x4)*complexconjugate(VV2x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN1x4)*complexconjugate(VV2x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1351 = Coupling(name = 'GC_1351',
                   value = '(ee*complex(0,1)*complexconjugate(NN2x4)*complexconjugate(VV2x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN2x4)*complexconjugate(VV2x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1352 = Coupling(name = 'GC_1352',
                   value = '(ee*complex(0,1)*complexconjugate(NN3x4)*complexconjugate(VV2x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN3x4)*complexconjugate(VV2x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN3x2)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN3x2)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1353 = Coupling(name = 'GC_1353',
                   value = '(ee*complex(0,1)*complexconjugate(NN4x4)*complexconjugate(VV2x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN4x4)*complexconjugate(VV2x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN4x2)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN4x2)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1354 = Coupling(name = 'GC_1354',
                   value = '(complex(0,1)*yd3x3*cmath.sin(alp))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1355 = Coupling(name = 'GC_1355',
                   value = '(complex(0,1)*ye3x3*cmath.sin(alp))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1356 = Coupling(name = 'GC_1356',
                   value = '-((complex(0,1)*yu3x3*cmath.sin(alp))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1357 = Coupling(name = 'GC_1357',
                   value = '(complex(0,1)*complexconjugate(yd3x3)*cmath.sin(alp))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1358 = Coupling(name = 'GC_1358',
                   value = '(complex(0,1)*complexconjugate(ye3x3)*cmath.sin(alp))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1359 = Coupling(name = 'GC_1359',
                   value = '-((complex(0,1)*complexconjugate(yu3x3)*cmath.sin(alp))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1360 = Coupling(name = 'GC_1360',
                   value = '(-2*ee**2*complex(0,1)*I149a44*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1361 = Coupling(name = 'GC_1361',
                   value = '(-2*ee**2*complex(0,1)*I149a55*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1362 = Coupling(name = 'GC_1362',
                   value = '(ee**2*complex(0,1)*I16a44*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1363 = Coupling(name = 'GC_1363',
                   value = '(ee**2*complex(0,1)*I16a55*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1364 = Coupling(name = 'GC_1364',
                   value = '(ee**2*complex(0,1)*I54a44*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1365 = Coupling(name = 'GC_1365',
                   value = '(ee**2*complex(0,1)*I54a55*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1366 = Coupling(name = 'GC_1366',
                   value = '-(ee**2*complex(0,1)*cmath.cos(alp)*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1367 = Coupling(name = 'GC_1367',
                   value = '(cw*ee*complex(0,1)*NN1x1*NN1x4*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN1x4*cmath.cos(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN1x4*sw*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN1x3*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN1x3*cmath.sin(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN1x3*sw*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1368 = Coupling(name = 'GC_1368',
                   value = '-((cw*ee*complex(0,1)*NN1x1*NN1x3*cmath.cos(alp))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*NN1x2*NN1x3*cmath.cos(alp))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN1x3*sw*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN1x4*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN1x4*cmath.sin(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN1x4*sw*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1369 = Coupling(name = 'GC_1369',
                   value = '(cw*ee*complex(0,1)*NN1x4*NN2x1*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN2x4*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x4*NN2x2*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN2x4*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x4*NN2x2*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN2x4*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x3*NN2x1*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN2x3*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x3*NN2x2*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN2x3*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x3*NN2x2*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN2x3*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1370 = Coupling(name = 'GC_1370',
                   value = '(cw*ee*complex(0,1)*NN2x1*NN2x4*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN2x4*cmath.cos(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN2x4*sw*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*NN2x3*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN2x3*cmath.sin(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN2x3*sw*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1371 = Coupling(name = 'GC_1371',
                   value = '-(cw*ee*complex(0,1)*NN1x3*NN2x1*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*NN2x3*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x3*NN2x2*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN2x3*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x3*NN2x2*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN2x3*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*NN2x1*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN2x4*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x4*NN2x2*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN2x4*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x4*NN2x2*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN2x4*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1372 = Coupling(name = 'GC_1372',
                   value = '-((cw*ee*complex(0,1)*NN2x1*NN2x3*cmath.cos(alp))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*NN2x2*NN2x3*cmath.cos(alp))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN2x3*sw*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*NN2x4*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN2x4*cmath.sin(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN2x4*sw*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1373 = Coupling(name = 'GC_1373',
                   value = '(cw*ee*complex(0,1)*NN1x4*NN3x1*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN3x4*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x4*NN3x2*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN3x4*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x4*NN3x2*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN3x4*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x3*NN3x1*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN3x3*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x3*NN3x2*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN3x3*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x3*NN3x2*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN3x3*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1374 = Coupling(name = 'GC_1374',
                   value = '(cw*ee*complex(0,1)*NN2x4*NN3x1*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*NN3x4*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x4*NN3x2*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN3x4*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x4*NN3x2*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN3x4*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x3*NN3x1*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*NN3x3*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x3*NN3x2*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN3x3*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x3*NN3x2*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN3x3*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1375 = Coupling(name = 'GC_1375',
                   value = '(cw*ee*complex(0,1)*NN3x1*NN3x4*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN3x4*cmath.cos(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN3x4*sw*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*NN3x3*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN3x3*cmath.sin(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN3x3*sw*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1376 = Coupling(name = 'GC_1376',
                   value = '-(cw*ee*complex(0,1)*NN1x3*NN3x1*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*NN3x3*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x3*NN3x2*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN3x3*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x3*NN3x2*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN3x3*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*NN3x1*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN3x4*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x4*NN3x2*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN3x4*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x4*NN3x2*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN3x4*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1377 = Coupling(name = 'GC_1377',
                   value = '-(cw*ee*complex(0,1)*NN2x3*NN3x1*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN2x1*NN3x3*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x3*NN3x2*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN3x3*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x3*NN3x2*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN3x3*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x4*NN3x1*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*NN3x4*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x4*NN3x2*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN3x4*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x4*NN3x2*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN3x4*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1378 = Coupling(name = 'GC_1378',
                   value = '-((cw*ee*complex(0,1)*NN3x1*NN3x3*cmath.cos(alp))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*NN3x2*NN3x3*cmath.cos(alp))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN3x3*sw*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*NN3x4*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN3x4*cmath.sin(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN3x4*sw*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1379 = Coupling(name = 'GC_1379',
                   value = '(cw*ee*complex(0,1)*NN1x4*NN4x1*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN4x4*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x4*NN4x2*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN4x4*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x4*NN4x2*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN4x4*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x3*NN4x1*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN4x3*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x3*NN4x2*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN4x3*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x3*NN4x2*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN4x3*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1380 = Coupling(name = 'GC_1380',
                   value = '(cw*ee*complex(0,1)*NN2x4*NN4x1*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*NN4x4*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x4*NN4x2*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN4x4*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x4*NN4x2*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN4x4*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x3*NN4x1*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*NN4x3*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x3*NN4x2*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN4x3*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x3*NN4x2*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN4x3*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1381 = Coupling(name = 'GC_1381',
                   value = '(cw*ee*complex(0,1)*NN3x4*NN4x1*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*NN4x4*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x4*NN4x2*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN4x4*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x4*NN4x2*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN4x4*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x3*NN4x1*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*NN4x3*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x3*NN4x2*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN4x3*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x3*NN4x2*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN4x3*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1382 = Coupling(name = 'GC_1382',
                   value = '(cw*ee*complex(0,1)*NN4x1*NN4x4*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN4x2*NN4x4*cmath.cos(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN4x2*NN4x4*sw*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*NN4x3*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN4x2*NN4x3*cmath.sin(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN4x2*NN4x3*sw*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1383 = Coupling(name = 'GC_1383',
                   value = '-(cw*ee*complex(0,1)*NN1x3*NN4x1*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*NN4x3*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x3*NN4x2*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN4x3*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x3*NN4x2*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN4x3*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*NN4x1*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN4x4*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x4*NN4x2*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN4x4*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x4*NN4x2*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN4x4*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1384 = Coupling(name = 'GC_1384',
                   value = '-(cw*ee*complex(0,1)*NN2x3*NN4x1*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN2x1*NN4x3*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x3*NN4x2*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN4x3*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x3*NN4x2*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN4x3*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x4*NN4x1*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*NN4x4*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x4*NN4x2*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN4x4*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x4*NN4x2*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN4x4*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1385 = Coupling(name = 'GC_1385',
                   value = '-(cw*ee*complex(0,1)*NN3x3*NN4x1*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN3x1*NN4x3*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN3x3*NN4x2*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN4x3*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x3*NN4x2*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN4x3*sw*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x4*NN4x1*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*NN4x4*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x4*NN4x2*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN4x4*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x4*NN4x2*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN4x4*sw*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1386 = Coupling(name = 'GC_1386',
                   value = '-((cw*ee*complex(0,1)*NN4x1*NN4x3*cmath.cos(alp))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*NN4x2*NN4x3*cmath.cos(alp))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN4x2*NN4x3*sw*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*NN4x4*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN4x2*NN4x4*cmath.sin(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN4x2*NN4x4*sw*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1387 = Coupling(name = 'GC_1387',
                   value = '-(ee**2*complex(0,1)*vu*cmath.cos(alp))/(4.*sw**2) + (ee**2*complex(0,1)*vd*cmath.sin(alp))/(4.*sw**2)',
                   order = {'QED':1})

GC_1388 = Coupling(name = 'GC_1388',
                   value = '(ee**2*complex(0,1)*vu*cmath.cos(alp))/(2.*sw**2) - (ee**2*complex(0,1)*vd*cmath.sin(alp))/(2.*sw**2)',
                   order = {'QED':1})

GC_1389 = Coupling(name = 'GC_1389',
                   value = '-(ee**2*complex(0,1)*vu*cmath.cos(alp))/2. - (cw**2*ee**2*complex(0,1)*vu*cmath.cos(alp))/(4.*sw**2) - (ee**2*complex(0,1)*sw**2*vu*cmath.cos(alp))/(4.*cw**2) + (ee**2*complex(0,1)*vd*cmath.sin(alp))/2. + (cw**2*ee**2*complex(0,1)*vd*cmath.sin(alp))/(4.*sw**2) + (ee**2*complex(0,1)*sw**2*vd*cmath.sin(alp))/(4.*cw**2)',
                   order = {'QED':1})

GC_1390 = Coupling(name = 'GC_1390',
                   value = '-(ee**2*complex(0,1)*I114a44*vu*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I114a44*vd*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1391 = Coupling(name = 'GC_1391',
                   value = '-(ee**2*complex(0,1)*I114a55*vu*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I114a55*vd*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1392 = Coupling(name = 'GC_1392',
                   value = '(ee**2*complex(0,1)*I48a44*vu*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I48a44*vd*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1393 = Coupling(name = 'GC_1393',
                   value = '(ee**2*complex(0,1)*I48a55*vu*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I48a55*vd*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1394 = Coupling(name = 'GC_1394',
                   value = '(ee**2*complex(0,1)*I9a44*vu*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I9a44*vd*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1395 = Coupling(name = 'GC_1395',
                   value = '(ee**2*complex(0,1)*I9a55*vu*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I9a55*vd*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1396 = Coupling(name = 'GC_1396',
                   value = '-(ee**2*complex(0,1)*vu*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vd*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1397 = Coupling(name = 'GC_1397',
                   value = '-(ee**2*complex(0,1)*vu*cmath.cos(alp))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vd*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1398 = Coupling(name = 'GC_1398',
                   value = '(ee**2*complex(0,1)*I113a11*vu*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a11*vu*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I113a11*vd*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a11*vd*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1399 = Coupling(name = 'GC_1399',
                   value = '(ee**2*complex(0,1)*I113a22*vu*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a22*vu*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I113a22*vd*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a22*vd*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1400 = Coupling(name = 'GC_1400',
                   value = '-(ee**2*complex(0,1)*I47a11*vu*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I47a11*vu*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I47a11*vd*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I47a11*vd*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1401 = Coupling(name = 'GC_1401',
                   value = '-(ee**2*complex(0,1)*I47a22*vu*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I47a22*vu*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I47a22*vd*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I47a22*vd*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1402 = Coupling(name = 'GC_1402',
                   value = '-(ee**2*complex(0,1)*I8a11*vu*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a11*vu*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I8a11*vd*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a11*vd*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1403 = Coupling(name = 'GC_1403',
                   value = '-(ee**2*complex(0,1)*I8a22*vu*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a22*vu*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I8a22*vd*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a22*vd*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1404 = Coupling(name = 'GC_1404',
                   value = '-(ee**2*complex(0,1)*I8a33*vu*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I9a33*vu*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a33*vu*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I18a33*MUH*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I18a33*MUH*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I21a33*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I21a33*sw**2*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I20a33*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I22a33*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I8a33*vd*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I9a33*vd*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a33*vd*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I20a33*sw**2*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I22a33*sw**2*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I17a33*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I19a33*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I17a33*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I19a33*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1405 = Coupling(name = 'GC_1405',
                   value = '-(ee**2*complex(0,1)*I8a36*vu*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I9a36*vu*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a36*vu*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I18a36*MUH*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I18a36*MUH*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I21a36*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I21a36*sw**2*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I20a36*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I22a36*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I8a36*vd*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I9a36*vd*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a36*vd*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I20a36*sw**2*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I22a36*sw**2*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I17a36*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I19a36*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I17a36*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I19a36*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1406 = Coupling(name = 'GC_1406',
                   value = '-(ee**2*complex(0,1)*I8a63*vu*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I9a63*vu*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a63*vu*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I18a63*MUH*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I18a63*MUH*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I21a63*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I21a63*sw**2*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I20a63*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I22a63*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I8a63*vd*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I9a63*vd*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a63*vd*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I20a63*sw**2*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I22a63*sw**2*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I17a63*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I19a63*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I17a63*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I19a63*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1407 = Coupling(name = 'GC_1407',
                   value = '-(ee**2*complex(0,1)*I8a66*vu*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I9a66*vu*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a66*vu*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I18a66*MUH*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I18a66*MUH*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I21a66*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I21a66*sw**2*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I20a66*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I22a66*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I8a66*vd*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I9a66*vd*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a66*vd*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I20a66*sw**2*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I22a66*sw**2*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I17a66*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I19a66*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I17a66*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I19a66*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1408 = Coupling(name = 'GC_1408',
                   value = '-(ee**2*complex(0,1)*I47a33*vu*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I48a33*vu*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I47a33*vu*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I56a33*MUH*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I56a33*MUH*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I59a33*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I59a33*sw**2*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee**2*complex(0,1)*I47a33*vd*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I48a33*vd*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I58a33*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I60a33*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I47a33*vd*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I58a33*sw**2*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I60a33*sw**2*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I55a33*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I57a33*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I55a33*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I57a33*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1409 = Coupling(name = 'GC_1409',
                   value = '-(ee**2*complex(0,1)*I47a36*vu*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I48a36*vu*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I47a36*vu*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I56a36*MUH*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I56a36*MUH*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I59a36*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I59a36*sw**2*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee**2*complex(0,1)*I47a36*vd*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I48a36*vd*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I58a36*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I60a36*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I47a36*vd*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I58a36*sw**2*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I60a36*sw**2*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I55a36*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I57a36*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I55a36*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I57a36*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1410 = Coupling(name = 'GC_1410',
                   value = '-(ee**2*complex(0,1)*I47a63*vu*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I48a63*vu*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I47a63*vu*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I56a63*MUH*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I56a63*MUH*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I59a63*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I59a63*sw**2*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee**2*complex(0,1)*I47a63*vd*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I48a63*vd*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I58a63*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I60a63*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I47a63*vd*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I58a63*sw**2*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I60a63*sw**2*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I55a63*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I57a63*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I55a63*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I57a63*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1411 = Coupling(name = 'GC_1411',
                   value = '-(ee**2*complex(0,1)*I47a66*vu*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I48a66*vu*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I47a66*vu*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I56a66*MUH*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I56a66*MUH*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I59a66*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I59a66*sw**2*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee**2*complex(0,1)*I47a66*vd*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I48a66*vd*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I58a66*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I60a66*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I47a66*vd*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I58a66*sw**2*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I60a66*sw**2*vd*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I55a66*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I57a66*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I55a66*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I57a66*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1412 = Coupling(name = 'GC_1412',
                   value = '-(ee**2*complex(0,1)*vd*cmath.cos(alp))/(4.*sw**2) - (ee**2*complex(0,1)*vu*cmath.sin(alp))/(4.*sw**2)',
                   order = {'QED':1})

GC_1413 = Coupling(name = 'GC_1413',
                   value = '(ee**2*complex(0,1)*vd*cmath.cos(alp))/(2.*sw**2) + (ee**2*complex(0,1)*vu*cmath.sin(alp))/(2.*sw**2)',
                   order = {'QED':1})

GC_1414 = Coupling(name = 'GC_1414',
                   value = '-(ee**2*complex(0,1)*vd*cmath.cos(alp))/2. - (cw**2*ee**2*complex(0,1)*vd*cmath.cos(alp))/(4.*sw**2) - (ee**2*complex(0,1)*sw**2*vd*cmath.cos(alp))/(4.*cw**2) - (ee**2*complex(0,1)*vu*cmath.sin(alp))/2. - (cw**2*ee**2*complex(0,1)*vu*cmath.sin(alp))/(4.*sw**2) - (ee**2*complex(0,1)*sw**2*vu*cmath.sin(alp))/(4.*cw**2)',
                   order = {'QED':1})

GC_1415 = Coupling(name = 'GC_1415',
                   value = '(ee**2*complex(0,1)*I114a44*vd*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I114a44*vu*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1416 = Coupling(name = 'GC_1416',
                   value = '(ee**2*complex(0,1)*I114a55*vd*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I114a55*vu*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1417 = Coupling(name = 'GC_1417',
                   value = '-(ee**2*complex(0,1)*I48a44*vd*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I48a44*vu*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1418 = Coupling(name = 'GC_1418',
                   value = '-(ee**2*complex(0,1)*I48a55*vd*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I48a55*vu*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1419 = Coupling(name = 'GC_1419',
                   value = '-(ee**2*complex(0,1)*I9a44*vd*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I9a44*vu*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1420 = Coupling(name = 'GC_1420',
                   value = '-(ee**2*complex(0,1)*I9a55*vd*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I9a55*vu*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1421 = Coupling(name = 'GC_1421',
                   value = '(ee**2*complex(0,1)*vd*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vu*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1422 = Coupling(name = 'GC_1422',
                   value = '-(ee**2*complex(0,1)*vd*cmath.cos(alp))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vu*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1423 = Coupling(name = 'GC_1423',
                   value = '-(ee**2*complex(0,1)*I113a11*vd*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I113a11*vd*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I113a11*vu*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a11*vu*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1424 = Coupling(name = 'GC_1424',
                   value = '-(ee**2*complex(0,1)*I113a22*vd*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I113a22*vd*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I113a22*vu*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a22*vu*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1425 = Coupling(name = 'GC_1425',
                   value = '(ee**2*complex(0,1)*I47a11*vd*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a11*vd*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I47a11*vu*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I47a11*vu*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1426 = Coupling(name = 'GC_1426',
                   value = '(ee**2*complex(0,1)*I47a22*vd*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a22*vd*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I47a22*vu*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I47a22*vu*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1427 = Coupling(name = 'GC_1427',
                   value = '(ee**2*complex(0,1)*I8a11*vd*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I8a11*vd*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I8a11*vu*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a11*vu*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1428 = Coupling(name = 'GC_1428',
                   value = '(ee**2*complex(0,1)*I8a22*vd*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I8a22*vd*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I8a22*vu*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a22*vu*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1429 = Coupling(name = 'GC_1429',
                   value = '-(ee**2*complex(0,1)*I113a33*vd*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I114a33*vd*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I113a33*vd*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I150a33*MUH*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I150a33*MUH*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I153a33*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I153a33*sw**2*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee**2*complex(0,1)*I113a33*vu*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I114a33*vu*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I154a33*vu*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I155a33*vu*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a33*vu*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I154a33*sw**2*vu*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I155a33*sw**2*vu*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I151a33*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I152a33*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I151a33*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I152a33*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1430 = Coupling(name = 'GC_1430',
                   value = '-(ee**2*complex(0,1)*I113a36*vd*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I114a36*vd*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I113a36*vd*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I150a36*MUH*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I150a36*MUH*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I153a36*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I153a36*sw**2*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee**2*complex(0,1)*I113a36*vu*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I114a36*vu*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I154a36*vu*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I155a36*vu*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a36*vu*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I154a36*sw**2*vu*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I155a36*sw**2*vu*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I151a36*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I152a36*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I151a36*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I152a36*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1431 = Coupling(name = 'GC_1431',
                   value = '-(ee**2*complex(0,1)*I113a63*vd*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I114a63*vd*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I113a63*vd*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I150a63*MUH*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I150a63*MUH*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I153a63*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I153a63*sw**2*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee**2*complex(0,1)*I113a63*vu*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I114a63*vu*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I154a63*vu*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I155a63*vu*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a63*vu*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I154a63*sw**2*vu*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I155a63*sw**2*vu*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I151a63*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I152a63*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I151a63*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I152a63*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1432 = Coupling(name = 'GC_1432',
                   value = '-(ee**2*complex(0,1)*I113a66*vd*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I114a66*vd*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I113a66*vd*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I150a66*MUH*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I150a66*MUH*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I153a66*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I153a66*sw**2*complexconjugate(MUH)*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee**2*complex(0,1)*I113a66*vu*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I114a66*vu*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I154a66*vu*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I155a66*vu*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a66*vu*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I154a66*sw**2*vu*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I155a66*sw**2*vu*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I151a66*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I152a66*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I151a66*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I152a66*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1433 = Coupling(name = 'GC_1433',
                   value = '-((ee*complex(0,1)*UU1x1*VV1x2*cmath.cos(alp))/(sw*cmath.sqrt(2))) + (ee*complex(0,1)*UU1x2*VV1x1*cmath.sin(alp))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1434 = Coupling(name = 'GC_1434',
                   value = '-((ee*complex(0,1)*UU2x1*VV1x2*cmath.cos(alp))/(sw*cmath.sqrt(2))) + (ee*complex(0,1)*UU2x2*VV1x1*cmath.sin(alp))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1435 = Coupling(name = 'GC_1435',
                   value = '-((ee*complex(0,1)*UU1x2*VV1x1*cmath.cos(alp))/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*UU1x1*VV1x2*cmath.sin(alp))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1436 = Coupling(name = 'GC_1436',
                   value = '-((ee*complex(0,1)*UU2x2*VV1x1*cmath.cos(alp))/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*UU2x1*VV1x2*cmath.sin(alp))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1437 = Coupling(name = 'GC_1437',
                   value = '-((ee*complex(0,1)*UU1x1*VV2x2*cmath.cos(alp))/(sw*cmath.sqrt(2))) + (ee*complex(0,1)*UU1x2*VV2x1*cmath.sin(alp))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1438 = Coupling(name = 'GC_1438',
                   value = '-((ee*complex(0,1)*UU2x1*VV2x2*cmath.cos(alp))/(sw*cmath.sqrt(2))) + (ee*complex(0,1)*UU2x2*VV2x1*cmath.sin(alp))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1439 = Coupling(name = 'GC_1439',
                   value = '-((ee*complex(0,1)*UU1x2*VV2x1*cmath.cos(alp))/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*UU1x1*VV2x2*cmath.sin(alp))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1440 = Coupling(name = 'GC_1440',
                   value = '-((ee*complex(0,1)*UU2x2*VV2x1*cmath.cos(alp))/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*UU2x1*VV2x2*cmath.sin(alp))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1441 = Coupling(name = 'GC_1441',
                   value = '(ee**2*complex(0,1)*I113a33*vu*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I114a33*vu*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I154a33*vu*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I155a33*vu*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a33*vu*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I154a33*sw**2*vu*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I155a33*sw**2*vu*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I151a33*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I152a33*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I151a33*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I152a33*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee**2*complex(0,1)*I113a33*vd*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I114a33*vd*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a33*vd*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I150a33*MUH*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I150a33*MUH*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I153a33*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I153a33*sw**2*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1442 = Coupling(name = 'GC_1442',
                   value = '(ee**2*complex(0,1)*I113a36*vu*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I114a36*vu*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I154a36*vu*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I155a36*vu*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a36*vu*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I154a36*sw**2*vu*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I155a36*sw**2*vu*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I151a36*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I152a36*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I151a36*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I152a36*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee**2*complex(0,1)*I113a36*vd*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I114a36*vd*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a36*vd*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I150a36*MUH*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I150a36*MUH*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I153a36*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I153a36*sw**2*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1443 = Coupling(name = 'GC_1443',
                   value = '(ee**2*complex(0,1)*I113a63*vu*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I114a63*vu*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I154a63*vu*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I155a63*vu*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a63*vu*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I154a63*sw**2*vu*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I155a63*sw**2*vu*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I151a63*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I152a63*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I151a63*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I152a63*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee**2*complex(0,1)*I113a63*vd*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I114a63*vd*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a63*vd*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I150a63*MUH*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I150a63*MUH*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I153a63*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I153a63*sw**2*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1444 = Coupling(name = 'GC_1444',
                   value = '(ee**2*complex(0,1)*I113a66*vu*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I114a66*vu*cmath.cos(alp))/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I154a66*vu*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I155a66*vu*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a66*vu*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I154a66*sw**2*vu*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I155a66*sw**2*vu*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I151a66*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I152a66*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I151a66*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I152a66*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee**2*complex(0,1)*I113a66*vd*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I114a66*vd*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I113a66*vd*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I150a66*MUH*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I150a66*MUH*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I153a66*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I153a66*sw**2*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1445 = Coupling(name = 'GC_1445',
                   value = '(complex(0,1)*I20a33*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I22a33*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a33*vd*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I9a33*vd*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I8a33*vd*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I20a33*sw**2*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I22a33*sw**2*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I17a33*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I19a33*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I17a33*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I19a33*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee**2*complex(0,1)*I8a33*vu*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I9a33*vu*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a33*vu*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I18a33*MUH*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I18a33*MUH*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I21a33*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I21a33*sw**2*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1446 = Coupling(name = 'GC_1446',
                   value = '(complex(0,1)*I20a36*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I22a36*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a36*vd*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I9a36*vd*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I8a36*vd*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I20a36*sw**2*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I22a36*sw**2*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I17a36*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I19a36*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I17a36*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I19a36*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee**2*complex(0,1)*I8a36*vu*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I9a36*vu*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a36*vu*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I18a36*MUH*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I18a36*MUH*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I21a36*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I21a36*sw**2*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1447 = Coupling(name = 'GC_1447',
                   value = '(complex(0,1)*I20a63*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I22a63*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a63*vd*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I9a63*vd*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I8a63*vd*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I20a63*sw**2*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I22a63*sw**2*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I17a63*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I19a63*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I17a63*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I19a63*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee**2*complex(0,1)*I8a63*vu*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I9a63*vu*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a63*vu*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I18a63*MUH*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I18a63*MUH*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I21a63*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I21a63*sw**2*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1448 = Coupling(name = 'GC_1448',
                   value = '(complex(0,1)*I20a66*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I22a66*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a66*vd*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I9a66*vd*cmath.cos(alp))/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I8a66*vd*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I20a66*sw**2*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I22a66*sw**2*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I17a66*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I19a66*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I17a66*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I19a66*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee**2*complex(0,1)*I8a66*vu*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I9a66*vu*cmath.sin(alp))/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I8a66*vu*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I18a66*MUH*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I18a66*MUH*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I21a66*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I21a66*sw**2*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1449 = Coupling(name = 'GC_1449',
                   value = '(ee**2*complex(0,1)*I47a33*vd*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I48a33*vd*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I58a33*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I60a33*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a33*vd*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I58a33*sw**2*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I60a33*sw**2*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I55a33*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I57a33*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I55a33*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I57a33*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee**2*complex(0,1)*I47a33*vu*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I48a33*vu*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I47a33*vu*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I56a33*MUH*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I56a33*MUH*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I59a33*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I59a33*sw**2*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1450 = Coupling(name = 'GC_1450',
                   value = '(ee**2*complex(0,1)*I47a36*vd*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I48a36*vd*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I58a36*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I60a36*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a36*vd*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I58a36*sw**2*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I60a36*sw**2*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I55a36*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I57a36*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I55a36*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I57a36*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee**2*complex(0,1)*I47a36*vu*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I48a36*vu*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I47a36*vu*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I56a36*MUH*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I56a36*MUH*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I59a36*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I59a36*sw**2*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1451 = Coupling(name = 'GC_1451',
                   value = '(ee**2*complex(0,1)*I47a63*vd*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I48a63*vd*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I58a63*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I60a63*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a63*vd*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I58a63*sw**2*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I60a63*sw**2*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I55a63*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I57a63*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I55a63*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I57a63*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee**2*complex(0,1)*I47a63*vu*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I48a63*vu*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I47a63*vu*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I56a63*MUH*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I56a63*MUH*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I59a63*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I59a63*sw**2*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1452 = Coupling(name = 'GC_1452',
                   value = '(ee**2*complex(0,1)*I47a66*vd*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I48a66*vd*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I58a66*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I60a66*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I47a66*vd*cmath.cos(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I58a66*sw**2*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I60a66*sw**2*vd*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I55a66*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I57a66*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I55a66*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I57a66*sw**2*cmath.cos(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee**2*complex(0,1)*I47a66*vu*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I48a66*vu*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I47a66*vu*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I56a66*MUH*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I56a66*MUH*sw**2*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I59a66*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I59a66*sw**2*complexconjugate(MUH)*cmath.sin(alp))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1453 = Coupling(name = 'GC_1453',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(NN1x4)*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(NN1x4)*cmath.cos(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(NN1x4)*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(NN1x3)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(NN1x3)*cmath.sin(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(NN1x3)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1454 = Coupling(name = 'GC_1454',
                   value = '-((cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(NN1x3)*cmath.cos(alp))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(NN1x3)*cmath.cos(alp))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(NN1x3)*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(NN1x4)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(NN1x4)*cmath.sin(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(NN1x4)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1455 = Coupling(name = 'GC_1455',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN1x4)*complexconjugate(NN2x1)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x4)*complexconjugate(NN2x2)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x4)*complexconjugate(NN2x2)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(NN2x4)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(NN2x4)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(NN2x4)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x3)*complexconjugate(NN2x1)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x3)*complexconjugate(NN2x2)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x3)*complexconjugate(NN2x2)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(NN2x3)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(NN2x3)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(NN2x3)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1456 = Coupling(name = 'GC_1456',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(NN2x4)*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(NN2x4)*cmath.cos(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(NN2x4)*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(NN2x3)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(NN2x3)*cmath.sin(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(NN2x3)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1457 = Coupling(name = 'GC_1457',
                   value = '-(cw*ee*complex(0,1)*complexconjugate(NN1x3)*complexconjugate(NN2x1)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN1x3)*complexconjugate(NN2x2)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN1x3)*complexconjugate(NN2x2)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(NN2x3)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(NN2x3)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(NN2x3)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x4)*complexconjugate(NN2x1)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x4)*complexconjugate(NN2x2)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x4)*complexconjugate(NN2x2)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(NN2x4)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(NN2x4)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(NN2x4)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1458 = Coupling(name = 'GC_1458',
                   value = '-((cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(NN2x3)*cmath.cos(alp))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(NN2x3)*cmath.cos(alp))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(NN2x3)*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(NN2x4)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(NN2x4)*cmath.sin(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(NN2x4)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1459 = Coupling(name = 'GC_1459',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN1x4)*complexconjugate(NN3x1)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x4)*complexconjugate(NN3x2)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x4)*complexconjugate(NN3x2)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(NN3x4)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(NN3x4)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(NN3x4)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x3)*complexconjugate(NN3x1)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x3)*complexconjugate(NN3x2)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x3)*complexconjugate(NN3x2)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(NN3x3)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(NN3x3)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(NN3x3)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1460 = Coupling(name = 'GC_1460',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN2x4)*complexconjugate(NN3x1)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN2x4)*complexconjugate(NN3x2)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN2x4)*complexconjugate(NN3x2)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(NN3x4)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(NN3x4)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(NN3x4)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x3)*complexconjugate(NN3x1)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN2x3)*complexconjugate(NN3x2)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN2x3)*complexconjugate(NN3x2)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(NN3x3)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(NN3x3)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(NN3x3)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1461 = Coupling(name = 'GC_1461',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(NN3x4)*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN3x2)*complexconjugate(NN3x4)*cmath.cos(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN3x2)*complexconjugate(NN3x4)*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(NN3x3)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN3x2)*complexconjugate(NN3x3)*cmath.sin(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN3x2)*complexconjugate(NN3x3)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1462 = Coupling(name = 'GC_1462',
                   value = '-(cw*ee*complex(0,1)*complexconjugate(NN1x3)*complexconjugate(NN3x1)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN1x3)*complexconjugate(NN3x2)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN1x3)*complexconjugate(NN3x2)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(NN3x3)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(NN3x3)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(NN3x3)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x4)*complexconjugate(NN3x1)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x4)*complexconjugate(NN3x2)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x4)*complexconjugate(NN3x2)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(NN3x4)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(NN3x4)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(NN3x4)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1463 = Coupling(name = 'GC_1463',
                   value = '-(cw*ee*complex(0,1)*complexconjugate(NN2x3)*complexconjugate(NN3x1)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN2x3)*complexconjugate(NN3x2)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN2x3)*complexconjugate(NN3x2)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(NN3x3)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(NN3x3)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(NN3x3)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x4)*complexconjugate(NN3x1)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN2x4)*complexconjugate(NN3x2)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN2x4)*complexconjugate(NN3x2)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(NN3x4)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(NN3x4)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(NN3x4)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1464 = Coupling(name = 'GC_1464',
                   value = '-((cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(NN3x3)*cmath.cos(alp))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*complexconjugate(NN3x2)*complexconjugate(NN3x3)*cmath.cos(alp))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN3x2)*complexconjugate(NN3x3)*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(NN3x4)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN3x2)*complexconjugate(NN3x4)*cmath.sin(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN3x2)*complexconjugate(NN3x4)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1465 = Coupling(name = 'GC_1465',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN1x4)*complexconjugate(NN4x1)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x4)*complexconjugate(NN4x2)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x4)*complexconjugate(NN4x2)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(NN4x4)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(NN4x4)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(NN4x4)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x3)*complexconjugate(NN4x1)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x3)*complexconjugate(NN4x2)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x3)*complexconjugate(NN4x2)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(NN4x3)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(NN4x3)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(NN4x3)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1466 = Coupling(name = 'GC_1466',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN2x4)*complexconjugate(NN4x1)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN2x4)*complexconjugate(NN4x2)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN2x4)*complexconjugate(NN4x2)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(NN4x4)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(NN4x4)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(NN4x4)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x3)*complexconjugate(NN4x1)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN2x3)*complexconjugate(NN4x2)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN2x3)*complexconjugate(NN4x2)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(NN4x3)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(NN4x3)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(NN4x3)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1467 = Coupling(name = 'GC_1467',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN3x4)*complexconjugate(NN4x1)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN3x4)*complexconjugate(NN4x2)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN3x4)*complexconjugate(NN4x2)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(NN4x4)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN3x2)*complexconjugate(NN4x4)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN3x2)*complexconjugate(NN4x4)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x3)*complexconjugate(NN4x1)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN3x3)*complexconjugate(NN4x2)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN3x3)*complexconjugate(NN4x2)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(NN4x3)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN3x2)*complexconjugate(NN4x3)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN3x2)*complexconjugate(NN4x3)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1468 = Coupling(name = 'GC_1468',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(NN4x4)*cmath.cos(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN4x2)*complexconjugate(NN4x4)*cmath.cos(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN4x2)*complexconjugate(NN4x4)*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(NN4x3)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN4x2)*complexconjugate(NN4x3)*cmath.sin(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN4x2)*complexconjugate(NN4x3)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1469 = Coupling(name = 'GC_1469',
                   value = '-(cw*ee*complex(0,1)*complexconjugate(NN1x3)*complexconjugate(NN4x1)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN1x3)*complexconjugate(NN4x2)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN1x3)*complexconjugate(NN4x2)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(NN4x3)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(NN4x3)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(NN4x3)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x4)*complexconjugate(NN4x1)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x4)*complexconjugate(NN4x2)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x4)*complexconjugate(NN4x2)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(NN4x4)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(NN4x4)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(NN4x4)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1470 = Coupling(name = 'GC_1470',
                   value = '-(cw*ee*complex(0,1)*complexconjugate(NN2x3)*complexconjugate(NN4x1)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN2x3)*complexconjugate(NN4x2)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN2x3)*complexconjugate(NN4x2)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(NN4x3)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(NN4x3)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(NN4x3)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x4)*complexconjugate(NN4x1)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN2x4)*complexconjugate(NN4x2)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN2x4)*complexconjugate(NN4x2)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(NN4x4)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(NN4x4)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(NN4x4)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1471 = Coupling(name = 'GC_1471',
                   value = '-(cw*ee*complex(0,1)*complexconjugate(NN3x3)*complexconjugate(NN4x1)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN3x3)*complexconjugate(NN4x2)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN3x3)*complexconjugate(NN4x2)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(NN4x3)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN3x2)*complexconjugate(NN4x3)*cmath.cos(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN3x2)*complexconjugate(NN4x3)*cmath.cos(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x4)*complexconjugate(NN4x1)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN3x4)*complexconjugate(NN4x2)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN3x4)*complexconjugate(NN4x2)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(NN4x4)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN3x2)*complexconjugate(NN4x4)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN3x2)*complexconjugate(NN4x4)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1472 = Coupling(name = 'GC_1472',
                   value = '-((cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(NN4x3)*cmath.cos(alp))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*complexconjugate(NN4x2)*complexconjugate(NN4x3)*cmath.cos(alp))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN4x2)*complexconjugate(NN4x3)*cmath.cos(alp))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(NN4x4)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*complexconjugate(NN4x2)*complexconjugate(NN4x4)*cmath.sin(alp))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*complexconjugate(NN4x2)*complexconjugate(NN4x4)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1473 = Coupling(name = 'GC_1473',
                   value = '-((ee*complex(0,1)*complexconjugate(UU1x1)*complexconjugate(VV1x2)*cmath.cos(alp))/(sw*cmath.sqrt(2))) + (ee*complex(0,1)*complexconjugate(UU1x2)*complexconjugate(VV1x1)*cmath.sin(alp))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1474 = Coupling(name = 'GC_1474',
                   value = '-((ee*complex(0,1)*complexconjugate(UU2x1)*complexconjugate(VV1x2)*cmath.cos(alp))/(sw*cmath.sqrt(2))) + (ee*complex(0,1)*complexconjugate(UU2x2)*complexconjugate(VV1x1)*cmath.sin(alp))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1475 = Coupling(name = 'GC_1475',
                   value = '-((ee*complex(0,1)*complexconjugate(UU1x2)*complexconjugate(VV1x1)*cmath.cos(alp))/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*complexconjugate(UU1x1)*complexconjugate(VV1x2)*cmath.sin(alp))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1476 = Coupling(name = 'GC_1476',
                   value = '-((ee*complex(0,1)*complexconjugate(UU2x2)*complexconjugate(VV1x1)*cmath.cos(alp))/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*complexconjugate(UU2x1)*complexconjugate(VV1x2)*cmath.sin(alp))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1477 = Coupling(name = 'GC_1477',
                   value = '-((ee*complex(0,1)*complexconjugate(UU1x1)*complexconjugate(VV2x2)*cmath.cos(alp))/(sw*cmath.sqrt(2))) + (ee*complex(0,1)*complexconjugate(UU1x2)*complexconjugate(VV2x1)*cmath.sin(alp))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1478 = Coupling(name = 'GC_1478',
                   value = '-((ee*complex(0,1)*complexconjugate(UU2x1)*complexconjugate(VV2x2)*cmath.cos(alp))/(sw*cmath.sqrt(2))) + (ee*complex(0,1)*complexconjugate(UU2x2)*complexconjugate(VV2x1)*cmath.sin(alp))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1479 = Coupling(name = 'GC_1479',
                   value = '-((ee*complex(0,1)*complexconjugate(UU1x2)*complexconjugate(VV2x1)*cmath.cos(alp))/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*complexconjugate(UU1x1)*complexconjugate(VV2x2)*cmath.sin(alp))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1480 = Coupling(name = 'GC_1480',
                   value = '-((ee*complex(0,1)*complexconjugate(UU2x2)*complexconjugate(VV2x1)*cmath.cos(alp))/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*complexconjugate(UU2x1)*complexconjugate(VV2x2)*cmath.sin(alp))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1481 = Coupling(name = 'GC_1481',
                   value = '(2*ee**2*complex(0,1)*I148a11*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a11*cmath.cos(alp)*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1482 = Coupling(name = 'GC_1482',
                   value = '(2*ee**2*complex(0,1)*I148a22*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a22*cmath.cos(alp)*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1483 = Coupling(name = 'GC_1483',
                   value = '-(ee**2*complex(0,1)*I15a11*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a11*cmath.cos(alp)*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1484 = Coupling(name = 'GC_1484',
                   value = '-(ee**2*complex(0,1)*I15a22*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a22*cmath.cos(alp)*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1485 = Coupling(name = 'GC_1485',
                   value = '-((ee**2*complex(0,1)*I53a11*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw))) + (ee**2*complex(0,1)*I53a11*cmath.cos(alp)*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1486 = Coupling(name = 'GC_1486',
                   value = '-((ee**2*complex(0,1)*I53a22*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw))) + (ee**2*complex(0,1)*I53a22*cmath.cos(alp)*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1487 = Coupling(name = 'GC_1487',
                   value = '(2*ee**2*complex(0,1)*I148a33*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I149a33*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a33*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a33*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a33*cmath.cos(alp)*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a33*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a33*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1488 = Coupling(name = 'GC_1488',
                   value = '(2*ee**2*complex(0,1)*I148a36*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I149a36*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a36*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a36*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a36*cmath.cos(alp)*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a36*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a36*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1489 = Coupling(name = 'GC_1489',
                   value = '(2*ee**2*complex(0,1)*I148a63*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I149a63*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a63*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a63*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a63*cmath.cos(alp)*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a63*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a63*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1490 = Coupling(name = 'GC_1490',
                   value = '(2*ee**2*complex(0,1)*I148a66*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I149a66*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a66*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a66*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a66*cmath.cos(alp)*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a66*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a66*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1491 = Coupling(name = 'GC_1491',
                   value = '-(ee**2*complex(0,1)*I15a33*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a33*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I23a33*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a33*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a33*cmath.cos(alp)*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I23a33*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a33*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1492 = Coupling(name = 'GC_1492',
                   value = '-(ee**2*complex(0,1)*I15a36*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a36*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I23a36*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a36*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a36*cmath.cos(alp)*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I23a36*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a36*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1493 = Coupling(name = 'GC_1493',
                   value = '-(ee**2*complex(0,1)*I15a63*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a63*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I23a63*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a63*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a63*cmath.cos(alp)*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I23a63*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a63*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1494 = Coupling(name = 'GC_1494',
                   value = '-(ee**2*complex(0,1)*I15a66*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a66*cmath.cos(alp)*cmath.sin(alp))/(3.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I23a66*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a66*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a66*cmath.cos(alp)*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I23a66*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a66*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1495 = Coupling(name = 'GC_1495',
                   value = '-((ee**2*complex(0,1)*I53a33*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw))) + (ee**2*complex(0,1)*I54a33*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I61a33*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a33*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a33*cmath.cos(alp)*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I61a33*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a33*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1496 = Coupling(name = 'GC_1496',
                   value = '-((ee**2*complex(0,1)*I53a36*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw))) + (ee**2*complex(0,1)*I54a36*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I61a36*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a36*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a36*cmath.cos(alp)*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I61a36*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a36*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1497 = Coupling(name = 'GC_1497',
                   value = '-((ee**2*complex(0,1)*I53a63*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw))) + (ee**2*complex(0,1)*I54a63*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I61a63*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a63*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a63*cmath.cos(alp)*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I61a63*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a63*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1498 = Coupling(name = 'GC_1498',
                   value = '-((ee**2*complex(0,1)*I53a66*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw))) + (ee**2*complex(0,1)*I54a66*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I61a66*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a66*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a66*cmath.cos(alp)*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I61a66*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a66*sw**2*cmath.cos(alp)*cmath.sin(alp))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1499 = Coupling(name = 'GC_1499',
                   value = '(ee**2*complex(0,1)*cmath.cos(alp)**2)/(2.*sw**2) + (ee**2*complex(0,1)*cmath.sin(alp)**2)/(2.*sw**2)',
                   order = {'QED':2})

GC_1500 = Coupling(name = 'GC_1500',
                   value = '(ee**2*complex(0,1)*I149a44*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a44*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1501 = Coupling(name = 'GC_1501',
                   value = '-(ee**2*complex(0,1)*I149a44*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a44*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1502 = Coupling(name = 'GC_1502',
                   value = '(ee**2*complex(0,1)*I149a55*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a55*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1503 = Coupling(name = 'GC_1503',
                   value = '-(ee**2*complex(0,1)*I149a55*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a55*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1504 = Coupling(name = 'GC_1504',
                   value = '(ee**2*complex(0,1)*I16a44*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a44*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1505 = Coupling(name = 'GC_1505',
                   value = '-(ee**2*complex(0,1)*I16a44*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a44*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1506 = Coupling(name = 'GC_1506',
                   value = '(ee**2*complex(0,1)*I16a55*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a55*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1507 = Coupling(name = 'GC_1507',
                   value = '-(ee**2*complex(0,1)*I16a55*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a55*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1508 = Coupling(name = 'GC_1508',
                   value = '(ee**2*complex(0,1)*I54a44*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I54a44*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1509 = Coupling(name = 'GC_1509',
                   value = '-(ee**2*complex(0,1)*I54a44*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a44*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1510 = Coupling(name = 'GC_1510',
                   value = '(ee**2*complex(0,1)*I54a55*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I54a55*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1511 = Coupling(name = 'GC_1511',
                   value = '-(ee**2*complex(0,1)*I54a55*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a55*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1512 = Coupling(name = 'GC_1512',
                   value = '(ee**2*complex(0,1)*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1513 = Coupling(name = 'GC_1513',
                   value = '-(ee**2*complex(0,1)*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1514 = Coupling(name = 'GC_1514',
                   value = '-(ee**2*complex(0,1)*cmath.cos(alp)**2)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.sin(alp)**2)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1515 = Coupling(name = 'GC_1515',
                   value = '-(ee**2*complex(0,1)*I148a11*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a11*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I148a11*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a11*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1516 = Coupling(name = 'GC_1516',
                   value = '(ee**2*complex(0,1)*I148a11*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a11*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I148a11*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a11*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1517 = Coupling(name = 'GC_1517',
                   value = '-(ee**2*complex(0,1)*I148a22*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a22*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I148a22*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a22*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1518 = Coupling(name = 'GC_1518',
                   value = '(ee**2*complex(0,1)*I148a22*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a22*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I148a22*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a22*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1519 = Coupling(name = 'GC_1519',
                   value = '(ee**2*complex(0,1)*I148a33*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a33*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a33*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a33*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a33*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a33*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a33*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a33*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a33*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a33*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1520 = Coupling(name = 'GC_1520',
                   value = '(ee**2*complex(0,1)*I148a36*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a36*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a36*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a36*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a36*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a36*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a36*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a36*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a36*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a36*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1521 = Coupling(name = 'GC_1521',
                   value = '(ee**2*complex(0,1)*I148a63*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a63*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a63*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a63*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a63*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a63*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a63*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a63*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a63*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a63*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1522 = Coupling(name = 'GC_1522',
                   value = '(ee**2*complex(0,1)*I148a66*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a66*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a66*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a66*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a66*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a66*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a66*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a66*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a66*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a66*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1523 = Coupling(name = 'GC_1523',
                   value = '-(ee**2*complex(0,1)*I15a11*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a11*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I15a11*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a11*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1524 = Coupling(name = 'GC_1524',
                   value = '(ee**2*complex(0,1)*I15a11*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a11*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I15a11*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a11*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1525 = Coupling(name = 'GC_1525',
                   value = '-(ee**2*complex(0,1)*I15a22*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a22*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I15a22*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a22*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1526 = Coupling(name = 'GC_1526',
                   value = '(ee**2*complex(0,1)*I15a22*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a22*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I15a22*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a22*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1527 = Coupling(name = 'GC_1527',
                   value = '(ee**2*complex(0,1)*I15a33*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a33*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I23a33*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a33*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a33*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I23a33*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a33*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a33*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a33*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a33*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1528 = Coupling(name = 'GC_1528',
                   value = '(ee**2*complex(0,1)*I15a36*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a36*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I23a36*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a36*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a36*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I23a36*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a36*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a36*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a36*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a36*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1529 = Coupling(name = 'GC_1529',
                   value = '(ee**2*complex(0,1)*I15a63*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a63*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I23a63*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a63*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a63*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I23a63*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a63*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a63*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a63*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a63*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1530 = Coupling(name = 'GC_1530',
                   value = '(ee**2*complex(0,1)*I15a66*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a66*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I23a66*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a66*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a66*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I23a66*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a66*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a66*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a66*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a66*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1531 = Coupling(name = 'GC_1531',
                   value = '-(ee**2*complex(0,1)*I53a11*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a11*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I53a11*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a11*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1532 = Coupling(name = 'GC_1532',
                   value = '(ee**2*complex(0,1)*I53a11*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a11*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I53a11*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a11*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1533 = Coupling(name = 'GC_1533',
                   value = '-(ee**2*complex(0,1)*I53a22*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a22*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I53a22*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a22*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1534 = Coupling(name = 'GC_1534',
                   value = '(ee**2*complex(0,1)*I53a22*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a22*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I53a22*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a22*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1535 = Coupling(name = 'GC_1535',
                   value = '(ee**2*complex(0,1)*I53a33*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I54a33*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I61a33*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a33*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a33*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I61a33*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a33*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a33*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a33*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a33*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1536 = Coupling(name = 'GC_1536',
                   value = '(ee**2*complex(0,1)*I53a36*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I54a36*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I61a36*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a36*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a36*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I61a36*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a36*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a36*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a36*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a36*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1537 = Coupling(name = 'GC_1537',
                   value = '(ee**2*complex(0,1)*I53a63*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I54a63*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I61a63*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a63*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a63*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I61a63*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a63*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a63*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a63*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a63*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1538 = Coupling(name = 'GC_1538',
                   value = '(ee**2*complex(0,1)*I53a66*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I54a66*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I61a66*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a66*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a66*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I61a66*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a66*sw**2*cmath.cos(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a66*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a66*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a66*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1539 = Coupling(name = 'GC_1539',
                   value = '-(ee**2*complex(0,1)*I148a33*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a33*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a33*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I148a33*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a33*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a33*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a33*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a33*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a33*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a33*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1540 = Coupling(name = 'GC_1540',
                   value = '-(ee**2*complex(0,1)*I148a36*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a36*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a36*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I148a36*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a36*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a36*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a36*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a36*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a36*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a36*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1541 = Coupling(name = 'GC_1541',
                   value = '-(ee**2*complex(0,1)*I148a63*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a63*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a63*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I148a63*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a63*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a63*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a63*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a63*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a63*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a63*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1542 = Coupling(name = 'GC_1542',
                   value = '-(ee**2*complex(0,1)*I148a66*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a66*cmath.cos(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a66*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I148a66*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a66*cmath.sin(alp)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a66*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a66*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a66*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a66*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a66*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1543 = Coupling(name = 'GC_1543',
                   value = '-(ee**2*complex(0,1)*I15a33*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a33*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a33*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I15a33*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a33*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I23a33*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a33*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a33*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I23a33*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a33*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1544 = Coupling(name = 'GC_1544',
                   value = '-(ee**2*complex(0,1)*I15a36*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a36*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a36*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I15a36*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a36*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I23a36*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a36*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a36*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I23a36*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a36*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1545 = Coupling(name = 'GC_1545',
                   value = '-(ee**2*complex(0,1)*I15a63*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a63*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a63*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I15a63*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a63*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I23a63*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a63*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a63*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I23a63*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a63*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1546 = Coupling(name = 'GC_1546',
                   value = '-(ee**2*complex(0,1)*I15a66*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a66*cmath.cos(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a66*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I15a66*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a66*cmath.sin(alp)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I23a66*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a66*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a66*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I23a66*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a66*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1547 = Coupling(name = 'GC_1547',
                   value = '-(ee**2*complex(0,1)*I53a33*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a33*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a33*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I53a33*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I54a33*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I61a33*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a33*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a33*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I61a33*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a33*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1548 = Coupling(name = 'GC_1548',
                   value = '-(ee**2*complex(0,1)*I53a36*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a36*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a36*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I53a36*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I54a36*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I61a36*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a36*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a36*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I61a36*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a36*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1549 = Coupling(name = 'GC_1549',
                   value = '-(ee**2*complex(0,1)*I53a63*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a63*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a63*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I53a63*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I54a63*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I61a63*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a63*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a63*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I61a63*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a63*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1550 = Coupling(name = 'GC_1550',
                   value = '-(ee**2*complex(0,1)*I53a66*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a66*cmath.cos(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a66*cmath.cos(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I53a66*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I54a66*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I61a66*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a66*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a66*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I61a66*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a66*sw**2*cmath.sin(alp)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1551 = Coupling(name = 'GC_1551',
                   value = '-(ee**2*complex(0,1)*vu*cmath.cos(alp)**3)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (5*ee**2*complex(0,1)*vd*cmath.cos(alp)**2*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (5*ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vd*cmath.sin(alp)**3)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1552 = Coupling(name = 'GC_1552',
                   value = '(3*ee**2*complex(0,1)*vu*cmath.cos(alp)**3)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (3*ee**2*complex(0,1)*vd*cmath.cos(alp)**2*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (3*ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (3*ee**2*complex(0,1)*vd*cmath.sin(alp)**3)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1553 = Coupling(name = 'GC_1553',
                   value = '-(ee**2*complex(0,1)*vd*cmath.cos(alp)**3)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (5*ee**2*complex(0,1)*vu*cmath.cos(alp)**2*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (5*ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vu*cmath.sin(alp)**3)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1554 = Coupling(name = 'GC_1554',
                   value = '(3*ee**2*complex(0,1)*vd*cmath.cos(alp)**3)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (3*ee**2*complex(0,1)*vu*cmath.cos(alp)**2*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (3*ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (3*ee**2*complex(0,1)*vu*cmath.sin(alp)**3)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1555 = Coupling(name = 'GC_1555',
                   value = '(3*ee**2*complex(0,1)*cmath.cos(alp)**3*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (3*ee**2*complex(0,1)*cmath.cos(alp)*cmath.sin(alp)**3)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1556 = Coupling(name = 'GC_1556',
                   value = '(-3*ee**2*complex(0,1)*cmath.cos(alp)**3*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (3*ee**2*complex(0,1)*cmath.cos(alp)*cmath.sin(alp)**3)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1557 = Coupling(name = 'GC_1557',
                   value = '-(ee**2*complex(0,1)*cmath.cos(alp)**4)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (5*ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.sin(alp)**2)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.sin(alp)**4)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1558 = Coupling(name = 'GC_1558',
                   value = '(3*ee**2*complex(0,1)*cmath.cos(alp)**4)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (3*ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.sin(alp)**2)/(2.*(-1 + sw)*sw**2*(1 + sw)) + (3*ee**2*complex(0,1)*cmath.sin(alp)**4)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1559 = Coupling(name = 'GC_1559',
                   value = 'complex(0,1)*I1a33*cmath.sin(beta)',
                   order = {'QED':1})

GC_1560 = Coupling(name = 'GC_1560',
                   value = 'complex(0,1)*I202a33*cmath.sin(beta)',
                   order = {'QED':1})

GC_1561 = Coupling(name = 'GC_1561',
                   value = 'complex(0,1)*I2a33*cmath.sin(beta)',
                   order = {'QED':1})

GC_1562 = Coupling(name = 'GC_1562',
                   value = 'complex(0,1)*I3a33*cmath.sin(beta)',
                   order = {'QED':1})

GC_1563 = Coupling(name = 'GC_1563',
                   value = 'complex(0,1)*I4a33*cmath.sin(beta)',
                   order = {'QED':1})

GC_1564 = Coupling(name = 'GC_1564',
                   value = 'complex(0,1)*I5a33*cmath.sin(beta)',
                   order = {'QED':1})

GC_1565 = Coupling(name = 'GC_1565',
                   value = '(yd3x3*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1566 = Coupling(name = 'GC_1566',
                   value = '(ye3x3*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1567 = Coupling(name = 'GC_1567',
                   value = '(yu3x3*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1568 = Coupling(name = 'GC_1568',
                   value = '-((complexconjugate(yd3x3)*cmath.sin(beta))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1569 = Coupling(name = 'GC_1569',
                   value = '-((complexconjugate(ye3x3)*cmath.sin(beta))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1570 = Coupling(name = 'GC_1570',
                   value = '-((complexconjugate(yu3x3)*cmath.sin(beta))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1571 = Coupling(name = 'GC_1571',
                   value = '-((ee**2*I109a11*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2)))',
                   order = {'QED':2})

GC_1572 = Coupling(name = 'GC_1572',
                   value = '-((ee**2*I109a22*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2)))',
                   order = {'QED':2})

GC_1573 = Coupling(name = 'GC_1573',
                   value = '(ee**2*I123a11*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1574 = Coupling(name = 'GC_1574',
                   value = '(ee**2*I123a22*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1575 = Coupling(name = 'GC_1575',
                   value = '-((ee**2*I144a11*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2)))',
                   order = {'QED':2})

GC_1576 = Coupling(name = 'GC_1576',
                   value = '-((ee**2*I144a22*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2)))',
                   order = {'QED':2})

GC_1577 = Coupling(name = 'GC_1577',
                   value = '(ee**2*I94a11*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1578 = Coupling(name = 'GC_1578',
                   value = '(ee**2*I94a22*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1579 = Coupling(name = 'GC_1579',
                   value = '(-2*ee**2*complex(0,1)*I149a44*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1580 = Coupling(name = 'GC_1580',
                   value = '(-2*ee**2*complex(0,1)*I149a55*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1581 = Coupling(name = 'GC_1581',
                   value = '(ee**2*complex(0,1)*I16a44*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1582 = Coupling(name = 'GC_1582',
                   value = '(ee**2*complex(0,1)*I16a55*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1583 = Coupling(name = 'GC_1583',
                   value = '(ee**2*complex(0,1)*I54a44*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1584 = Coupling(name = 'GC_1584',
                   value = '(ee**2*complex(0,1)*I54a55*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1585 = Coupling(name = 'GC_1585',
                   value = '-(ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1586 = Coupling(name = 'GC_1586',
                   value = '-(ee**2*complex(0,1)*I53a11*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1587 = Coupling(name = 'GC_1587',
                   value = '-(ee**2*complex(0,1)*I53a22*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1588 = Coupling(name = 'GC_1588',
                   value = '(ee**2*complex(0,1)*cmath.cos(alp)*cmath.cos(beta)*cmath.sin(alp)*cmath.sin(beta))/((-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1589 = Coupling(name = 'GC_1589',
                   value = '(2*ee**2*complex(0,1)*cmath.cos(beta)**2*cmath.sin(beta)**2)/((-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1590 = Coupling(name = 'GC_1590',
                   value = '(I150a33*MUH*cmath.cos(beta))/cmath.sqrt(2) - (I153a33*complexconjugate(MUH)*cmath.cos(beta))/cmath.sqrt(2) - (I151a33*cmath.sin(beta))/cmath.sqrt(2) + (I152a33*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1591 = Coupling(name = 'GC_1591',
                   value = '(I150a36*MUH*cmath.cos(beta))/cmath.sqrt(2) - (I153a36*complexconjugate(MUH)*cmath.cos(beta))/cmath.sqrt(2) - (I151a36*cmath.sin(beta))/cmath.sqrt(2) + (I152a36*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1592 = Coupling(name = 'GC_1592',
                   value = '(I150a63*MUH*cmath.cos(beta))/cmath.sqrt(2) - (I153a63*complexconjugate(MUH)*cmath.cos(beta))/cmath.sqrt(2) - (I151a63*cmath.sin(beta))/cmath.sqrt(2) + (I152a63*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1593 = Coupling(name = 'GC_1593',
                   value = '(I150a66*MUH*cmath.cos(beta))/cmath.sqrt(2) - (I153a66*complexconjugate(MUH)*cmath.cos(beta))/cmath.sqrt(2) - (I151a66*cmath.sin(beta))/cmath.sqrt(2) + (I152a66*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1594 = Coupling(name = 'GC_1594',
                   value = '-((I18a33*MUH*cmath.cos(beta))/cmath.sqrt(2)) + (I21a33*complexconjugate(MUH)*cmath.cos(beta))/cmath.sqrt(2) - (I17a33*cmath.sin(beta))/cmath.sqrt(2) + (I19a33*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1595 = Coupling(name = 'GC_1595',
                   value = '-((I18a36*MUH*cmath.cos(beta))/cmath.sqrt(2)) + (I21a36*complexconjugate(MUH)*cmath.cos(beta))/cmath.sqrt(2) - (I17a36*cmath.sin(beta))/cmath.sqrt(2) + (I19a36*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1596 = Coupling(name = 'GC_1596',
                   value = '-((I18a63*MUH*cmath.cos(beta))/cmath.sqrt(2)) + (I21a63*complexconjugate(MUH)*cmath.cos(beta))/cmath.sqrt(2) - (I17a63*cmath.sin(beta))/cmath.sqrt(2) + (I19a63*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1597 = Coupling(name = 'GC_1597',
                   value = '-((I18a66*MUH*cmath.cos(beta))/cmath.sqrt(2)) + (I21a66*complexconjugate(MUH)*cmath.cos(beta))/cmath.sqrt(2) - (I17a66*cmath.sin(beta))/cmath.sqrt(2) + (I19a66*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1598 = Coupling(name = 'GC_1598',
                   value = '-((I56a33*MUH*cmath.cos(beta))/cmath.sqrt(2)) + (I59a33*complexconjugate(MUH)*cmath.cos(beta))/cmath.sqrt(2) - (I55a33*cmath.sin(beta))/cmath.sqrt(2) + (I57a33*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1599 = Coupling(name = 'GC_1599',
                   value = '-((I56a36*MUH*cmath.cos(beta))/cmath.sqrt(2)) + (I59a36*complexconjugate(MUH)*cmath.cos(beta))/cmath.sqrt(2) - (I55a36*cmath.sin(beta))/cmath.sqrt(2) + (I57a36*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1600 = Coupling(name = 'GC_1600',
                   value = '-((I56a63*MUH*cmath.cos(beta))/cmath.sqrt(2)) + (I59a63*complexconjugate(MUH)*cmath.cos(beta))/cmath.sqrt(2) - (I55a63*cmath.sin(beta))/cmath.sqrt(2) + (I57a63*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1601 = Coupling(name = 'GC_1601',
                   value = '-((I56a66*MUH*cmath.cos(beta))/cmath.sqrt(2)) + (I59a66*complexconjugate(MUH)*cmath.cos(beta))/cmath.sqrt(2) - (I55a66*cmath.sin(beta))/cmath.sqrt(2) + (I57a66*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1602 = Coupling(name = 'GC_1602',
                   value = '-((cw*ee*NN1x1*NN1x4*cmath.cos(beta))/((-1 + sw)*(1 + sw))) + (ee*NN1x2*NN1x4*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*NN1x2*NN1x4*sw*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*NN1x1*NN1x3*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee*NN1x2*NN1x3*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) + (ee*NN1x2*NN1x3*sw*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1603 = Coupling(name = 'GC_1603',
                   value = '-((cw*ee*NN1x1*NN1x3*cmath.cos(beta))/((-1 + sw)*(1 + sw))) + (ee*NN1x2*NN1x3*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*NN1x2*NN1x3*sw*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*NN1x1*NN1x4*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee*NN1x2*NN1x4*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*NN1x2*NN1x4*sw*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1604 = Coupling(name = 'GC_1604',
                   value = '-(cw*ee*NN1x4*NN2x1*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x1*NN2x4*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x4*NN2x2*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x2*NN2x4*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x4*NN2x2*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x2*NN2x4*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN1x3*NN2x1*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN1x1*NN2x3*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x3*NN2x2*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x2*NN2x3*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x3*NN2x2*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x2*NN2x3*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1605 = Coupling(name = 'GC_1605',
                   value = '-((cw*ee*NN2x1*NN2x4*cmath.cos(beta))/((-1 + sw)*(1 + sw))) + (ee*NN2x2*NN2x4*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*NN2x2*NN2x4*sw*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*NN2x1*NN2x3*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee*NN2x2*NN2x3*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) + (ee*NN2x2*NN2x3*sw*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1606 = Coupling(name = 'GC_1606',
                   value = '-(cw*ee*NN1x3*NN2x1*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x1*NN2x3*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x3*NN2x2*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x2*NN2x3*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x3*NN2x2*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x2*NN2x3*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x4*NN2x1*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x1*NN2x4*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x4*NN2x2*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x2*NN2x4*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x4*NN2x2*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x2*NN2x4*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1607 = Coupling(name = 'GC_1607',
                   value = '-((cw*ee*NN2x1*NN2x3*cmath.cos(beta))/((-1 + sw)*(1 + sw))) + (ee*NN2x2*NN2x3*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*NN2x2*NN2x3*sw*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*NN2x1*NN2x4*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee*NN2x2*NN2x4*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*NN2x2*NN2x4*sw*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1608 = Coupling(name = 'GC_1608',
                   value = '-(cw*ee*NN1x4*NN3x1*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x1*NN3x4*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x4*NN3x2*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x2*NN3x4*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x4*NN3x2*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x2*NN3x4*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN1x3*NN3x1*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN1x1*NN3x3*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x3*NN3x2*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x2*NN3x3*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x3*NN3x2*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x2*NN3x3*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1609 = Coupling(name = 'GC_1609',
                   value = '-(cw*ee*NN2x4*NN3x1*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN2x1*NN3x4*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN2x4*NN3x2*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN2x2*NN3x4*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN2x4*NN3x2*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN2x2*NN3x4*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN2x3*NN3x1*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN2x1*NN3x3*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN2x3*NN3x2*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN2x2*NN3x3*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN2x3*NN3x2*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN2x2*NN3x3*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1610 = Coupling(name = 'GC_1610',
                   value = '-((cw*ee*NN3x1*NN3x4*cmath.cos(beta))/((-1 + sw)*(1 + sw))) + (ee*NN3x2*NN3x4*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*NN3x2*NN3x4*sw*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*NN3x1*NN3x3*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee*NN3x2*NN3x3*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) + (ee*NN3x2*NN3x3*sw*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1611 = Coupling(name = 'GC_1611',
                   value = '-(cw*ee*NN1x3*NN3x1*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x1*NN3x3*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x3*NN3x2*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x2*NN3x3*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x3*NN3x2*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x2*NN3x3*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x4*NN3x1*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x1*NN3x4*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x4*NN3x2*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x2*NN3x4*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x4*NN3x2*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x2*NN3x4*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1612 = Coupling(name = 'GC_1612',
                   value = '-(cw*ee*NN2x3*NN3x1*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN2x1*NN3x3*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN2x3*NN3x2*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN2x2*NN3x3*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN2x3*NN3x2*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN2x2*NN3x3*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN2x4*NN3x1*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN2x1*NN3x4*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN2x4*NN3x2*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN2x2*NN3x4*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN2x4*NN3x2*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN2x2*NN3x4*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1613 = Coupling(name = 'GC_1613',
                   value = '-((cw*ee*NN3x1*NN3x3*cmath.cos(beta))/((-1 + sw)*(1 + sw))) + (ee*NN3x2*NN3x3*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*NN3x2*NN3x3*sw*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*NN3x1*NN3x4*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee*NN3x2*NN3x4*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*NN3x2*NN3x4*sw*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1614 = Coupling(name = 'GC_1614',
                   value = '-(cw*ee*NN1x4*NN4x1*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x1*NN4x4*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x4*NN4x2*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x2*NN4x4*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x4*NN4x2*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x2*NN4x4*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN1x3*NN4x1*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN1x1*NN4x3*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x3*NN4x2*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x2*NN4x3*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x3*NN4x2*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x2*NN4x3*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1615 = Coupling(name = 'GC_1615',
                   value = '-(cw*ee*NN2x4*NN4x1*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN2x1*NN4x4*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN2x4*NN4x2*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN2x2*NN4x4*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN2x4*NN4x2*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN2x2*NN4x4*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN2x3*NN4x1*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN2x1*NN4x3*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN2x3*NN4x2*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN2x2*NN4x3*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN2x3*NN4x2*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN2x2*NN4x3*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1616 = Coupling(name = 'GC_1616',
                   value = '-(cw*ee*NN3x4*NN4x1*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN3x1*NN4x4*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN3x4*NN4x2*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN3x2*NN4x4*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN3x4*NN4x2*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN3x2*NN4x4*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN3x3*NN4x1*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN3x1*NN4x3*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN3x3*NN4x2*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN3x2*NN4x3*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN3x3*NN4x2*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN3x2*NN4x3*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1617 = Coupling(name = 'GC_1617',
                   value = '-((cw*ee*NN4x1*NN4x4*cmath.cos(beta))/((-1 + sw)*(1 + sw))) + (ee*NN4x2*NN4x4*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*NN4x2*NN4x4*sw*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*NN4x1*NN4x3*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee*NN4x2*NN4x3*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) + (ee*NN4x2*NN4x3*sw*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1618 = Coupling(name = 'GC_1618',
                   value = '-(cw*ee*NN1x3*NN4x1*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x1*NN4x3*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x3*NN4x2*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x2*NN4x3*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x3*NN4x2*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x2*NN4x3*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x4*NN4x1*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x1*NN4x4*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x4*NN4x2*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x2*NN4x4*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x4*NN4x2*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x2*NN4x4*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1619 = Coupling(name = 'GC_1619',
                   value = '-(cw*ee*NN2x3*NN4x1*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN2x1*NN4x3*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN2x3*NN4x2*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN2x2*NN4x3*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN2x3*NN4x2*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN2x2*NN4x3*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN2x4*NN4x1*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN2x1*NN4x4*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN2x4*NN4x2*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN2x2*NN4x4*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN2x4*NN4x2*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN2x2*NN4x4*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1620 = Coupling(name = 'GC_1620',
                   value = '-(cw*ee*NN3x3*NN4x1*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN3x1*NN4x3*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN3x3*NN4x2*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN3x2*NN4x3*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN3x3*NN4x2*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN3x2*NN4x3*sw*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN3x4*NN4x1*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN3x1*NN4x4*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*NN3x4*NN4x2*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN3x2*NN4x4*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN3x4*NN4x2*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*NN3x2*NN4x4*sw*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1621 = Coupling(name = 'GC_1621',
                   value = '-((cw*ee*NN4x1*NN4x3*cmath.cos(beta))/((-1 + sw)*(1 + sw))) + (ee*NN4x2*NN4x3*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*NN4x2*NN4x3*sw*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*NN4x1*NN4x4*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee*NN4x2*NN4x4*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*NN4x2*NN4x4*sw*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1622 = Coupling(name = 'GC_1622',
                   value = '(ee*complex(0,1)*NN1x3*UU1x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x3*sw*UU1x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*UU1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*UU1x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*sw*UU1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1623 = Coupling(name = 'GC_1623',
                   value = '(ee*complex(0,1)*NN2x3*UU1x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x3*sw*UU1x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN2x1*UU1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*UU1x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*sw*UU1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1624 = Coupling(name = 'GC_1624',
                   value = '(ee*complex(0,1)*NN3x3*UU1x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x3*sw*UU1x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN3x1*UU1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*UU1x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*sw*UU1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1625 = Coupling(name = 'GC_1625',
                   value = '(ee*complex(0,1)*NN4x3*UU1x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN4x3*sw*UU1x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN4x1*UU1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*UU1x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*sw*UU1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1626 = Coupling(name = 'GC_1626',
                   value = '(ee*complex(0,1)*NN1x3*UU2x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x3*sw*UU2x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*UU2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*UU2x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*sw*UU2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1627 = Coupling(name = 'GC_1627',
                   value = '(ee*complex(0,1)*NN2x3*UU2x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x3*sw*UU2x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN2x1*UU2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*UU2x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*sw*UU2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1628 = Coupling(name = 'GC_1628',
                   value = '(ee*complex(0,1)*NN3x3*UU2x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x3*sw*UU2x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN3x1*UU2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*UU2x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*sw*UU2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1629 = Coupling(name = 'GC_1629',
                   value = '(ee*complex(0,1)*NN4x3*UU2x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN4x3*sw*UU2x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN4x1*UU2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*UU2x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*sw*UU2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1630 = Coupling(name = 'GC_1630',
                   value = '(ee**2*vu*cmath.cos(beta))/(4.*sw**2) - (ee**2*vd*cmath.sin(beta))/(4.*sw**2)',
                   order = {'QED':1})

GC_1631 = Coupling(name = 'GC_1631',
                   value = '-(ee**2*vu*cmath.cos(beta))/(4.*sw**2) + (ee**2*vd*cmath.sin(beta))/(4.*sw**2)',
                   order = {'QED':1})

GC_1632 = Coupling(name = 'GC_1632',
                   value = '(ee**2*complex(0,1)*vu*cmath.cos(beta))/(4.*cw) + (cw*ee**2*complex(0,1)*vu*cmath.cos(beta))/(4.*sw**2) - (ee**2*complex(0,1)*vd*cmath.sin(beta))/(4.*cw) - (cw*ee**2*complex(0,1)*vd*cmath.sin(beta))/(4.*sw**2)',
                   order = {'QED':1})

GC_1633 = Coupling(name = 'GC_1633',
                   value = '(ee**2*complex(0,1)*vu*cmath.cos(beta))/(4.*cw) - (cw*ee**2*complex(0,1)*vu*cmath.cos(beta))/(4.*sw**2) - (ee**2*complex(0,1)*vd*cmath.sin(beta))/(4.*cw) + (cw*ee**2*complex(0,1)*vd*cmath.sin(beta))/(4.*sw**2)',
                   order = {'QED':1})

GC_1634 = Coupling(name = 'GC_1634',
                   value = '-(ee**2*complex(0,1)*I105a11*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I105a11*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1635 = Coupling(name = 'GC_1635',
                   value = '-(ee**2*complex(0,1)*I105a22*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I105a22*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1636 = Coupling(name = 'GC_1636',
                   value = 'complex(0,1)*I107a33*MUH*cmath.cos(beta) - (ee**2*complex(0,1)*I105a33*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I106a33*cmath.sin(beta) + (complex(0,1)*I108a33*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I105a33*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1637 = Coupling(name = 'GC_1637',
                   value = 'complex(0,1)*I107a36*MUH*cmath.cos(beta) - (ee**2*complex(0,1)*I105a36*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I106a36*cmath.sin(beta) + (complex(0,1)*I108a36*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I105a36*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1638 = Coupling(name = 'GC_1638',
                   value = '-(ee**2*complex(0,1)*I115a11*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I115a11*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1639 = Coupling(name = 'GC_1639',
                   value = '-(ee**2*complex(0,1)*I115a22*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I115a22*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1640 = Coupling(name = 'GC_1640',
                   value = '-(ee**2*complex(0,1)*I136a11*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I136a11*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1641 = Coupling(name = 'GC_1641',
                   value = '-(ee**2*complex(0,1)*I136a22*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I136a22*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1642 = Coupling(name = 'GC_1642',
                   value = '-(ee**2*complex(0,1)*I90a11*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I90a11*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1643 = Coupling(name = 'GC_1643',
                   value = '-(ee**2*complex(0,1)*I90a22*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I90a22*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1644 = Coupling(name = 'GC_1644',
                   value = 'complex(0,1)*I93a33*complexconjugate(MUH)*cmath.cos(beta) - (ee**2*complex(0,1)*I90a33*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I91a33*cmath.sin(beta) + (complex(0,1)*I92a33*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I90a33*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1645 = Coupling(name = 'GC_1645',
                   value = 'complex(0,1)*I93a36*complexconjugate(MUH)*cmath.cos(beta) - (ee**2*complex(0,1)*I90a36*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I91a36*cmath.sin(beta) + (complex(0,1)*I92a36*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I90a36*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1646 = Coupling(name = 'GC_1646',
                   value = '(ee**2*complex(0,1)*vu*cmath.cos(beta))/(2.*sw) - (ee**2*complex(0,1)*vd*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':1})

GC_1647 = Coupling(name = 'GC_1647',
                   value = '-(ee**2*complex(0,1)*vu*cmath.cos(beta))/(2.*sw) + (ee**2*complex(0,1)*vd*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':1})

GC_1648 = Coupling(name = 'GC_1648',
                   value = '(cw*ee**2*complex(0,1)*vu*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee**2*complex(0,1)*vd*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1649 = Coupling(name = 'GC_1649',
                   value = 'complex(0,1)*I117a33*cmath.cos(beta) + complex(0,1)*I120a33*complexconjugate(MUH)*cmath.cos(beta) + (complex(0,1)*I121a33*vd*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I122a33*vu*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I115a33*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I118a33*cmath.sin(beta) + complex(0,1)*I116a33*MUH*cmath.sin(beta) + (complex(0,1)*I119a33*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I115a33*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I121a33*vu*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1650 = Coupling(name = 'GC_1650',
                   value = 'complex(0,1)*I117a36*cmath.cos(beta) + complex(0,1)*I120a36*complexconjugate(MUH)*cmath.cos(beta) + (complex(0,1)*I121a36*vd*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I122a36*vu*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I115a36*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I118a36*cmath.sin(beta) + complex(0,1)*I116a36*MUH*cmath.sin(beta) + (complex(0,1)*I119a36*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I115a36*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I121a36*vu*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1651 = Coupling(name = 'GC_1651',
                   value = 'complex(0,1)*I117a63*cmath.cos(beta) + complex(0,1)*I120a63*complexconjugate(MUH)*cmath.cos(beta) + (complex(0,1)*I121a63*vd*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I122a63*vu*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I115a63*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I118a63*cmath.sin(beta) + complex(0,1)*I116a63*MUH*cmath.sin(beta) + (complex(0,1)*I119a63*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I115a63*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I121a63*vu*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1652 = Coupling(name = 'GC_1652',
                   value = 'complex(0,1)*I117a66*cmath.cos(beta) + complex(0,1)*I120a66*complexconjugate(MUH)*cmath.cos(beta) + (complex(0,1)*I121a66*vd*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I122a66*vu*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I115a66*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I118a66*cmath.sin(beta) + complex(0,1)*I116a66*MUH*cmath.sin(beta) + (complex(0,1)*I119a66*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I115a66*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I121a66*vu*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1653 = Coupling(name = 'GC_1653',
                   value = '-(ee**2*vd*cmath.cos(beta))/(4.*sw**2) - (ee**2*vu*cmath.sin(beta))/(4.*sw**2)',
                   order = {'QED':1})

GC_1654 = Coupling(name = 'GC_1654',
                   value = '(ee**2*vd*cmath.cos(beta))/(4.*sw**2) + (ee**2*vu*cmath.sin(beta))/(4.*sw**2)',
                   order = {'QED':1})

GC_1655 = Coupling(name = 'GC_1655',
                   value = '(ee**2*complex(0,1)*vd*cmath.cos(beta))/(4.*cw) - (cw*ee**2*complex(0,1)*vd*cmath.cos(beta))/(4.*sw**2) + (ee**2*complex(0,1)*vu*cmath.sin(beta))/(4.*cw) - (cw*ee**2*complex(0,1)*vu*cmath.sin(beta))/(4.*sw**2)',
                   order = {'QED':1})

GC_1656 = Coupling(name = 'GC_1656',
                   value = '(ee**2*complex(0,1)*vd*cmath.cos(beta))/(4.*cw) + (cw*ee**2*complex(0,1)*vd*cmath.cos(beta))/(4.*sw**2) + (ee**2*complex(0,1)*vu*cmath.sin(beta))/(4.*cw) + (cw*ee**2*complex(0,1)*vu*cmath.sin(beta))/(4.*sw**2)',
                   order = {'QED':1})

GC_1657 = Coupling(name = 'GC_1657',
                   value = '(ee**2*complex(0,1)*I105a11*vd*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I105a11*vu*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1658 = Coupling(name = 'GC_1658',
                   value = '(ee**2*complex(0,1)*I105a22*vd*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I105a22*vu*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1659 = Coupling(name = 'GC_1659',
                   value = '-(complex(0,1)*I106a33*cmath.cos(beta)) - (complex(0,1)*I108a33*vd*cmath.cos(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I105a33*vd*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I107a33*MUH*cmath.sin(beta) - (ee**2*complex(0,1)*I105a33*vu*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1660 = Coupling(name = 'GC_1660',
                   value = '-(complex(0,1)*I106a36*cmath.cos(beta)) - (complex(0,1)*I108a36*vd*cmath.cos(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I105a36*vd*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I107a36*MUH*cmath.sin(beta) - (ee**2*complex(0,1)*I105a36*vu*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1661 = Coupling(name = 'GC_1661',
                   value = '(ee**2*complex(0,1)*I115a11*vd*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I115a11*vu*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1662 = Coupling(name = 'GC_1662',
                   value = '(ee**2*complex(0,1)*I115a22*vd*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I115a22*vu*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1663 = Coupling(name = 'GC_1663',
                   value = '(ee**2*complex(0,1)*I136a11*vd*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I136a11*vu*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1664 = Coupling(name = 'GC_1664',
                   value = '(ee**2*complex(0,1)*I136a22*vd*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I136a22*vu*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1665 = Coupling(name = 'GC_1665',
                   value = '-(complex(0,1)*I137a33*cmath.cos(beta)) - complex(0,1)*I141a33*complexconjugate(MUH)*cmath.cos(beta) - (complex(0,1)*I140a33*vd*cmath.cos(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I136a33*vd*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I143a33*vu*cmath.cos(beta))/cmath.sqrt(2) + complex(0,1)*I139a33*cmath.sin(beta) + complex(0,1)*I138a33*MUH*cmath.sin(beta) + (complex(0,1)*I143a33*vd*cmath.sin(beta))/cmath.sqrt(2) + (complex(0,1)*I142a33*vu*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I136a33*vu*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1666 = Coupling(name = 'GC_1666',
                   value = '-(complex(0,1)*I137a36*cmath.cos(beta)) - complex(0,1)*I141a36*complexconjugate(MUH)*cmath.cos(beta) - (complex(0,1)*I140a36*vd*cmath.cos(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I136a36*vd*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I143a36*vu*cmath.cos(beta))/cmath.sqrt(2) + complex(0,1)*I139a36*cmath.sin(beta) + complex(0,1)*I138a36*MUH*cmath.sin(beta) + (complex(0,1)*I143a36*vd*cmath.sin(beta))/cmath.sqrt(2) + (complex(0,1)*I142a36*vu*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I136a36*vu*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1667 = Coupling(name = 'GC_1667',
                   value = '-(complex(0,1)*I137a63*cmath.cos(beta)) - complex(0,1)*I141a63*complexconjugate(MUH)*cmath.cos(beta) - (complex(0,1)*I140a63*vd*cmath.cos(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I136a63*vd*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I143a63*vu*cmath.cos(beta))/cmath.sqrt(2) + complex(0,1)*I139a63*cmath.sin(beta) + complex(0,1)*I138a63*MUH*cmath.sin(beta) + (complex(0,1)*I143a63*vd*cmath.sin(beta))/cmath.sqrt(2) + (complex(0,1)*I142a63*vu*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I136a63*vu*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1668 = Coupling(name = 'GC_1668',
                   value = '-(complex(0,1)*I137a66*cmath.cos(beta)) - complex(0,1)*I141a66*complexconjugate(MUH)*cmath.cos(beta) - (complex(0,1)*I140a66*vd*cmath.cos(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I136a66*vd*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I143a66*vu*cmath.cos(beta))/cmath.sqrt(2) + complex(0,1)*I139a66*cmath.sin(beta) + complex(0,1)*I138a66*MUH*cmath.sin(beta) + (complex(0,1)*I143a66*vd*cmath.sin(beta))/cmath.sqrt(2) + (complex(0,1)*I142a66*vu*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I136a66*vu*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1669 = Coupling(name = 'GC_1669',
                   value = '(ee**2*complex(0,1)*I90a11*vd*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I90a11*vu*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1670 = Coupling(name = 'GC_1670',
                   value = '(ee**2*complex(0,1)*I90a22*vd*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I90a22*vu*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1671 = Coupling(name = 'GC_1671',
                   value = '-(ee**2*complex(0,1)*vd*cmath.cos(beta))/(2.*sw) - (ee**2*complex(0,1)*vu*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':1})

GC_1672 = Coupling(name = 'GC_1672',
                   value = '(ee**2*complex(0,1)*vd*cmath.cos(beta))/(2.*sw) + (ee**2*complex(0,1)*vu*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':1})

GC_1673 = Coupling(name = 'GC_1673',
                   value = '(cw*ee**2*complex(0,1)*vd*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee**2*complex(0,1)*vu*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1674 = Coupling(name = 'GC_1674',
                   value = '(ee*UU1x1*VV1x2*cmath.cos(beta))/(sw*cmath.sqrt(2)) + (ee*UU1x2*VV1x1*cmath.sin(beta))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1675 = Coupling(name = 'GC_1675',
                   value = '(ee*UU2x1*VV1x2*cmath.cos(beta))/(sw*cmath.sqrt(2)) + (ee*UU2x2*VV1x1*cmath.sin(beta))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1676 = Coupling(name = 'GC_1676',
                   value = '(ee*complex(0,1)*NN1x4*VV1x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x4*sw*VV1x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*VV1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*VV1x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*VV1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1677 = Coupling(name = 'GC_1677',
                   value = '(ee*complex(0,1)*NN2x4*VV1x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x4*sw*VV1x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*VV1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*VV1x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*VV1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1678 = Coupling(name = 'GC_1678',
                   value = '(ee*complex(0,1)*NN3x4*VV1x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x4*sw*VV1x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*VV1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*VV1x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*VV1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1679 = Coupling(name = 'GC_1679',
                   value = '(ee*complex(0,1)*NN4x4*VV1x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN4x4*sw*VV1x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*VV1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*VV1x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*VV1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1680 = Coupling(name = 'GC_1680',
                   value = '-((ee*UU1x2*VV1x1*cmath.cos(beta))/(sw*cmath.sqrt(2))) + (ee*UU1x1*VV1x2*cmath.sin(beta))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1681 = Coupling(name = 'GC_1681',
                   value = '-((ee*UU2x2*VV1x1*cmath.cos(beta))/(sw*cmath.sqrt(2))) + (ee*UU2x1*VV1x2*cmath.sin(beta))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1682 = Coupling(name = 'GC_1682',
                   value = '(ee*UU1x1*VV2x2*cmath.cos(beta))/(sw*cmath.sqrt(2)) + (ee*UU1x2*VV2x1*cmath.sin(beta))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1683 = Coupling(name = 'GC_1683',
                   value = '(ee*UU2x1*VV2x2*cmath.cos(beta))/(sw*cmath.sqrt(2)) + (ee*UU2x2*VV2x1*cmath.sin(beta))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1684 = Coupling(name = 'GC_1684',
                   value = '(ee*complex(0,1)*NN1x4*VV2x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x4*sw*VV2x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*VV2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*VV2x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*VV2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1685 = Coupling(name = 'GC_1685',
                   value = '(ee*complex(0,1)*NN2x4*VV2x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x4*sw*VV2x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*VV2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*VV2x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*VV2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1686 = Coupling(name = 'GC_1686',
                   value = '(ee*complex(0,1)*NN3x4*VV2x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x4*sw*VV2x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*VV2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*VV2x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*VV2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1687 = Coupling(name = 'GC_1687',
                   value = '(ee*complex(0,1)*NN4x4*VV2x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN4x4*sw*VV2x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*VV2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*VV2x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*VV2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1688 = Coupling(name = 'GC_1688',
                   value = '-((ee*UU1x2*VV2x1*cmath.cos(beta))/(sw*cmath.sqrt(2))) + (ee*UU1x1*VV2x2*cmath.sin(beta))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1689 = Coupling(name = 'GC_1689',
                   value = '-((ee*UU2x2*VV2x1*cmath.cos(beta))/(sw*cmath.sqrt(2))) + (ee*UU2x1*VV2x2*cmath.sin(beta))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1690 = Coupling(name = 'GC_1690',
                   value = '-(complex(0,1)*I118a33*cmath.cos(beta)) - complex(0,1)*I116a33*MUH*cmath.cos(beta) - (complex(0,1)*I119a33*vd*cmath.cos(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I115a33*vd*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I121a33*vu*cmath.cos(beta))/cmath.sqrt(2) + complex(0,1)*I117a33*cmath.sin(beta) + complex(0,1)*I120a33*complexconjugate(MUH)*cmath.sin(beta) + (complex(0,1)*I121a33*vd*cmath.sin(beta))/cmath.sqrt(2) + (complex(0,1)*I122a33*vu*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I115a33*vu*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1691 = Coupling(name = 'GC_1691',
                   value = '-(complex(0,1)*I118a36*cmath.cos(beta)) - complex(0,1)*I116a36*MUH*cmath.cos(beta) - (complex(0,1)*I119a36*vd*cmath.cos(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I115a36*vd*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I121a36*vu*cmath.cos(beta))/cmath.sqrt(2) + complex(0,1)*I117a36*cmath.sin(beta) + complex(0,1)*I120a36*complexconjugate(MUH)*cmath.sin(beta) + (complex(0,1)*I121a36*vd*cmath.sin(beta))/cmath.sqrt(2) + (complex(0,1)*I122a36*vu*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I115a36*vu*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1692 = Coupling(name = 'GC_1692',
                   value = '-(complex(0,1)*I118a63*cmath.cos(beta)) - complex(0,1)*I116a63*MUH*cmath.cos(beta) - (complex(0,1)*I119a63*vd*cmath.cos(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I115a63*vd*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I121a63*vu*cmath.cos(beta))/cmath.sqrt(2) + complex(0,1)*I117a63*cmath.sin(beta) + complex(0,1)*I120a63*complexconjugate(MUH)*cmath.sin(beta) + (complex(0,1)*I121a63*vd*cmath.sin(beta))/cmath.sqrt(2) + (complex(0,1)*I122a63*vu*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I115a63*vu*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1693 = Coupling(name = 'GC_1693',
                   value = '-(complex(0,1)*I118a66*cmath.cos(beta)) - complex(0,1)*I116a66*MUH*cmath.cos(beta) - (complex(0,1)*I119a66*vd*cmath.cos(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I115a66*vd*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I121a66*vu*cmath.cos(beta))/cmath.sqrt(2) + complex(0,1)*I117a66*cmath.sin(beta) + complex(0,1)*I120a66*complexconjugate(MUH)*cmath.sin(beta) + (complex(0,1)*I121a66*vd*cmath.sin(beta))/cmath.sqrt(2) + (complex(0,1)*I122a66*vu*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I115a66*vu*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1694 = Coupling(name = 'GC_1694',
                   value = 'complex(0,1)*I139a33*cmath.cos(beta) + complex(0,1)*I138a33*MUH*cmath.cos(beta) + (complex(0,1)*I143a33*vd*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I142a33*vu*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I136a33*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I137a33*cmath.sin(beta) + complex(0,1)*I141a33*complexconjugate(MUH)*cmath.sin(beta) + (complex(0,1)*I140a33*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I136a33*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I143a33*vu*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1695 = Coupling(name = 'GC_1695',
                   value = 'complex(0,1)*I139a36*cmath.cos(beta) + complex(0,1)*I138a36*MUH*cmath.cos(beta) + (complex(0,1)*I143a36*vd*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I142a36*vu*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I136a36*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I137a36*cmath.sin(beta) + complex(0,1)*I141a36*complexconjugate(MUH)*cmath.sin(beta) + (complex(0,1)*I140a36*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I136a36*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I143a36*vu*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1696 = Coupling(name = 'GC_1696',
                   value = 'complex(0,1)*I139a63*cmath.cos(beta) + complex(0,1)*I138a63*MUH*cmath.cos(beta) + (complex(0,1)*I143a63*vd*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I142a63*vu*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I136a63*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I137a63*cmath.sin(beta) + complex(0,1)*I141a63*complexconjugate(MUH)*cmath.sin(beta) + (complex(0,1)*I140a63*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I136a63*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I143a63*vu*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1697 = Coupling(name = 'GC_1697',
                   value = 'complex(0,1)*I139a66*cmath.cos(beta) + complex(0,1)*I138a66*MUH*cmath.cos(beta) + (complex(0,1)*I143a66*vd*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I142a66*vu*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I136a66*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I137a66*cmath.sin(beta) + complex(0,1)*I141a66*complexconjugate(MUH)*cmath.sin(beta) + (complex(0,1)*I140a66*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I136a66*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I143a66*vu*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1698 = Coupling(name = 'GC_1698',
                   value = '-((I151a33*cmath.cos(beta))/cmath.sqrt(2)) + (I152a33*cmath.cos(beta))/cmath.sqrt(2) - (I150a33*MUH*cmath.sin(beta))/cmath.sqrt(2) + (I153a33*complexconjugate(MUH)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1699 = Coupling(name = 'GC_1699',
                   value = '-((I151a36*cmath.cos(beta))/cmath.sqrt(2)) + (I152a36*cmath.cos(beta))/cmath.sqrt(2) - (I150a36*MUH*cmath.sin(beta))/cmath.sqrt(2) + (I153a36*complexconjugate(MUH)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1700 = Coupling(name = 'GC_1700',
                   value = '-((I151a63*cmath.cos(beta))/cmath.sqrt(2)) + (I152a63*cmath.cos(beta))/cmath.sqrt(2) - (I150a63*MUH*cmath.sin(beta))/cmath.sqrt(2) + (I153a63*complexconjugate(MUH)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1701 = Coupling(name = 'GC_1701',
                   value = '-((I151a66*cmath.cos(beta))/cmath.sqrt(2)) + (I152a66*cmath.cos(beta))/cmath.sqrt(2) - (I150a66*MUH*cmath.sin(beta))/cmath.sqrt(2) + (I153a66*complexconjugate(MUH)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1702 = Coupling(name = 'GC_1702',
                   value = '(I17a33*cmath.cos(beta))/cmath.sqrt(2) - (I19a33*cmath.cos(beta))/cmath.sqrt(2) - (I18a33*MUH*cmath.sin(beta))/cmath.sqrt(2) + (I21a33*complexconjugate(MUH)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1703 = Coupling(name = 'GC_1703',
                   value = '(I17a36*cmath.cos(beta))/cmath.sqrt(2) - (I19a36*cmath.cos(beta))/cmath.sqrt(2) - (I18a36*MUH*cmath.sin(beta))/cmath.sqrt(2) + (I21a36*complexconjugate(MUH)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1704 = Coupling(name = 'GC_1704',
                   value = '(I17a63*cmath.cos(beta))/cmath.sqrt(2) - (I19a63*cmath.cos(beta))/cmath.sqrt(2) - (I18a63*MUH*cmath.sin(beta))/cmath.sqrt(2) + (I21a63*complexconjugate(MUH)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1705 = Coupling(name = 'GC_1705',
                   value = '(I17a66*cmath.cos(beta))/cmath.sqrt(2) - (I19a66*cmath.cos(beta))/cmath.sqrt(2) - (I18a66*MUH*cmath.sin(beta))/cmath.sqrt(2) + (I21a66*complexconjugate(MUH)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1706 = Coupling(name = 'GC_1706',
                   value = '(I55a33*cmath.cos(beta))/cmath.sqrt(2) - (I57a33*cmath.cos(beta))/cmath.sqrt(2) - (I56a33*MUH*cmath.sin(beta))/cmath.sqrt(2) + (I59a33*complexconjugate(MUH)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1707 = Coupling(name = 'GC_1707',
                   value = '(I55a36*cmath.cos(beta))/cmath.sqrt(2) - (I57a36*cmath.cos(beta))/cmath.sqrt(2) - (I56a36*MUH*cmath.sin(beta))/cmath.sqrt(2) + (I59a36*complexconjugate(MUH)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1708 = Coupling(name = 'GC_1708',
                   value = '(I55a63*cmath.cos(beta))/cmath.sqrt(2) - (I57a63*cmath.cos(beta))/cmath.sqrt(2) - (I56a63*MUH*cmath.sin(beta))/cmath.sqrt(2) + (I59a63*complexconjugate(MUH)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1709 = Coupling(name = 'GC_1709',
                   value = '(I55a66*cmath.cos(beta))/cmath.sqrt(2) - (I57a66*cmath.cos(beta))/cmath.sqrt(2) - (I56a66*MUH*cmath.sin(beta))/cmath.sqrt(2) + (I59a66*complexconjugate(MUH)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1710 = Coupling(name = 'GC_1710',
                   value = '-(complex(0,1)*I91a33*cmath.cos(beta)) - (complex(0,1)*I92a33*vd*cmath.cos(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I90a33*vd*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I93a33*complexconjugate(MUH)*cmath.sin(beta) - (ee**2*complex(0,1)*I90a33*vu*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1711 = Coupling(name = 'GC_1711',
                   value = '-(complex(0,1)*I91a36*cmath.cos(beta)) - (complex(0,1)*I92a36*vd*cmath.cos(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I90a36*vd*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I93a36*complexconjugate(MUH)*cmath.sin(beta) - (ee**2*complex(0,1)*I90a36*vu*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1712 = Coupling(name = 'GC_1712',
                   value = '(cw*ee*complexconjugate(NN1x1)*complexconjugate(NN1x4)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x2)*complexconjugate(NN1x4)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x2)*complexconjugate(NN1x4)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complexconjugate(NN1x1)*complexconjugate(NN1x3)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee*complexconjugate(NN1x2)*complexconjugate(NN1x3)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*sw*complexconjugate(NN1x2)*complexconjugate(NN1x3)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1713 = Coupling(name = 'GC_1713',
                   value = '(cw*ee*complexconjugate(NN1x1)*complexconjugate(NN1x3)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x2)*complexconjugate(NN1x3)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x2)*complexconjugate(NN1x3)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN1x1)*complexconjugate(NN1x4)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x2)*complexconjugate(NN1x4)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x2)*complexconjugate(NN1x4)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1714 = Coupling(name = 'GC_1714',
                   value = '(cw*ee*complexconjugate(NN1x4)*complexconjugate(NN2x1)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x4)*complexconjugate(NN2x2)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x4)*complexconjugate(NN2x2)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN1x1)*complexconjugate(NN2x4)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x2)*complexconjugate(NN2x4)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x2)*complexconjugate(NN2x4)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complexconjugate(NN1x3)*complexconjugate(NN2x1)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*complexconjugate(NN1x3)*complexconjugate(NN2x2)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*complexconjugate(NN1x3)*complexconjugate(NN2x2)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complexconjugate(NN1x1)*complexconjugate(NN2x3)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*complexconjugate(NN1x2)*complexconjugate(NN2x3)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*complexconjugate(NN1x2)*complexconjugate(NN2x3)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1715 = Coupling(name = 'GC_1715',
                   value = '(cw*ee*complexconjugate(NN2x1)*complexconjugate(NN2x4)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN2x2)*complexconjugate(NN2x4)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN2x2)*complexconjugate(NN2x4)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complexconjugate(NN2x1)*complexconjugate(NN2x3)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee*complexconjugate(NN2x2)*complexconjugate(NN2x3)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*sw*complexconjugate(NN2x2)*complexconjugate(NN2x3)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1716 = Coupling(name = 'GC_1716',
                   value = '(cw*ee*complexconjugate(NN1x3)*complexconjugate(NN2x1)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x3)*complexconjugate(NN2x2)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x3)*complexconjugate(NN2x2)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN1x1)*complexconjugate(NN2x3)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x2)*complexconjugate(NN2x3)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x2)*complexconjugate(NN2x3)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN1x4)*complexconjugate(NN2x1)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x4)*complexconjugate(NN2x2)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x4)*complexconjugate(NN2x2)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN1x1)*complexconjugate(NN2x4)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x2)*complexconjugate(NN2x4)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x2)*complexconjugate(NN2x4)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1717 = Coupling(name = 'GC_1717',
                   value = '(cw*ee*complexconjugate(NN2x1)*complexconjugate(NN2x3)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN2x2)*complexconjugate(NN2x3)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN2x2)*complexconjugate(NN2x3)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN2x1)*complexconjugate(NN2x4)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN2x2)*complexconjugate(NN2x4)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN2x2)*complexconjugate(NN2x4)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1718 = Coupling(name = 'GC_1718',
                   value = '(cw*ee*complexconjugate(NN1x4)*complexconjugate(NN3x1)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x4)*complexconjugate(NN3x2)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x4)*complexconjugate(NN3x2)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN1x1)*complexconjugate(NN3x4)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x2)*complexconjugate(NN3x4)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x2)*complexconjugate(NN3x4)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complexconjugate(NN1x3)*complexconjugate(NN3x1)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*complexconjugate(NN1x3)*complexconjugate(NN3x2)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*complexconjugate(NN1x3)*complexconjugate(NN3x2)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complexconjugate(NN1x1)*complexconjugate(NN3x3)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*complexconjugate(NN1x2)*complexconjugate(NN3x3)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*complexconjugate(NN1x2)*complexconjugate(NN3x3)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1719 = Coupling(name = 'GC_1719',
                   value = '(cw*ee*complexconjugate(NN2x4)*complexconjugate(NN3x1)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN2x4)*complexconjugate(NN3x2)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN2x4)*complexconjugate(NN3x2)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN2x1)*complexconjugate(NN3x4)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN2x2)*complexconjugate(NN3x4)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN2x2)*complexconjugate(NN3x4)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complexconjugate(NN2x3)*complexconjugate(NN3x1)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*complexconjugate(NN2x3)*complexconjugate(NN3x2)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*complexconjugate(NN2x3)*complexconjugate(NN3x2)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complexconjugate(NN2x1)*complexconjugate(NN3x3)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*complexconjugate(NN2x2)*complexconjugate(NN3x3)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*complexconjugate(NN2x2)*complexconjugate(NN3x3)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1720 = Coupling(name = 'GC_1720',
                   value = '(cw*ee*complexconjugate(NN3x1)*complexconjugate(NN3x4)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN3x2)*complexconjugate(NN3x4)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN3x2)*complexconjugate(NN3x4)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complexconjugate(NN3x1)*complexconjugate(NN3x3)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee*complexconjugate(NN3x2)*complexconjugate(NN3x3)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*sw*complexconjugate(NN3x2)*complexconjugate(NN3x3)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1721 = Coupling(name = 'GC_1721',
                   value = '(cw*ee*complexconjugate(NN1x3)*complexconjugate(NN3x1)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x3)*complexconjugate(NN3x2)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x3)*complexconjugate(NN3x2)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN1x1)*complexconjugate(NN3x3)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x2)*complexconjugate(NN3x3)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x2)*complexconjugate(NN3x3)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN1x4)*complexconjugate(NN3x1)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x4)*complexconjugate(NN3x2)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x4)*complexconjugate(NN3x2)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN1x1)*complexconjugate(NN3x4)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x2)*complexconjugate(NN3x4)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x2)*complexconjugate(NN3x4)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1722 = Coupling(name = 'GC_1722',
                   value = '(cw*ee*complexconjugate(NN2x3)*complexconjugate(NN3x1)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN2x3)*complexconjugate(NN3x2)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN2x3)*complexconjugate(NN3x2)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN2x1)*complexconjugate(NN3x3)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN2x2)*complexconjugate(NN3x3)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN2x2)*complexconjugate(NN3x3)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN2x4)*complexconjugate(NN3x1)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN2x4)*complexconjugate(NN3x2)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN2x4)*complexconjugate(NN3x2)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN2x1)*complexconjugate(NN3x4)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN2x2)*complexconjugate(NN3x4)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN2x2)*complexconjugate(NN3x4)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1723 = Coupling(name = 'GC_1723',
                   value = '(cw*ee*complexconjugate(NN3x1)*complexconjugate(NN3x3)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN3x2)*complexconjugate(NN3x3)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN3x2)*complexconjugate(NN3x3)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN3x1)*complexconjugate(NN3x4)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN3x2)*complexconjugate(NN3x4)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN3x2)*complexconjugate(NN3x4)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1724 = Coupling(name = 'GC_1724',
                   value = '(cw*ee*complexconjugate(NN1x4)*complexconjugate(NN4x1)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x4)*complexconjugate(NN4x2)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x4)*complexconjugate(NN4x2)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN1x1)*complexconjugate(NN4x4)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x2)*complexconjugate(NN4x4)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x2)*complexconjugate(NN4x4)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complexconjugate(NN1x3)*complexconjugate(NN4x1)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*complexconjugate(NN1x3)*complexconjugate(NN4x2)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*complexconjugate(NN1x3)*complexconjugate(NN4x2)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complexconjugate(NN1x1)*complexconjugate(NN4x3)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*complexconjugate(NN1x2)*complexconjugate(NN4x3)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*complexconjugate(NN1x2)*complexconjugate(NN4x3)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1725 = Coupling(name = 'GC_1725',
                   value = '(cw*ee*complexconjugate(NN2x4)*complexconjugate(NN4x1)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN2x4)*complexconjugate(NN4x2)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN2x4)*complexconjugate(NN4x2)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN2x1)*complexconjugate(NN4x4)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN2x2)*complexconjugate(NN4x4)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN2x2)*complexconjugate(NN4x4)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complexconjugate(NN2x3)*complexconjugate(NN4x1)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*complexconjugate(NN2x3)*complexconjugate(NN4x2)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*complexconjugate(NN2x3)*complexconjugate(NN4x2)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complexconjugate(NN2x1)*complexconjugate(NN4x3)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*complexconjugate(NN2x2)*complexconjugate(NN4x3)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*complexconjugate(NN2x2)*complexconjugate(NN4x3)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1726 = Coupling(name = 'GC_1726',
                   value = '(cw*ee*complexconjugate(NN3x4)*complexconjugate(NN4x1)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN3x4)*complexconjugate(NN4x2)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN3x4)*complexconjugate(NN4x2)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN3x1)*complexconjugate(NN4x4)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN3x2)*complexconjugate(NN4x4)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN3x2)*complexconjugate(NN4x4)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complexconjugate(NN3x3)*complexconjugate(NN4x1)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*complexconjugate(NN3x3)*complexconjugate(NN4x2)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*complexconjugate(NN3x3)*complexconjugate(NN4x2)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complexconjugate(NN3x1)*complexconjugate(NN4x3)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee*complexconjugate(NN3x2)*complexconjugate(NN4x3)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*complexconjugate(NN3x2)*complexconjugate(NN4x3)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1727 = Coupling(name = 'GC_1727',
                   value = '(cw*ee*complexconjugate(NN4x1)*complexconjugate(NN4x4)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN4x2)*complexconjugate(NN4x4)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN4x2)*complexconjugate(NN4x4)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complexconjugate(NN4x1)*complexconjugate(NN4x3)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee*complexconjugate(NN4x2)*complexconjugate(NN4x3)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*sw*complexconjugate(NN4x2)*complexconjugate(NN4x3)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1728 = Coupling(name = 'GC_1728',
                   value = '(cw*ee*complexconjugate(NN1x3)*complexconjugate(NN4x1)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x3)*complexconjugate(NN4x2)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x3)*complexconjugate(NN4x2)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN1x1)*complexconjugate(NN4x3)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x2)*complexconjugate(NN4x3)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x2)*complexconjugate(NN4x3)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN1x4)*complexconjugate(NN4x1)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x4)*complexconjugate(NN4x2)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x4)*complexconjugate(NN4x2)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN1x1)*complexconjugate(NN4x4)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN1x2)*complexconjugate(NN4x4)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN1x2)*complexconjugate(NN4x4)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1729 = Coupling(name = 'GC_1729',
                   value = '(cw*ee*complexconjugate(NN2x3)*complexconjugate(NN4x1)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN2x3)*complexconjugate(NN4x2)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN2x3)*complexconjugate(NN4x2)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN2x1)*complexconjugate(NN4x3)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN2x2)*complexconjugate(NN4x3)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN2x2)*complexconjugate(NN4x3)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN2x4)*complexconjugate(NN4x1)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN2x4)*complexconjugate(NN4x2)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN2x4)*complexconjugate(NN4x2)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN2x1)*complexconjugate(NN4x4)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN2x2)*complexconjugate(NN4x4)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN2x2)*complexconjugate(NN4x4)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1730 = Coupling(name = 'GC_1730',
                   value = '(cw*ee*complexconjugate(NN3x3)*complexconjugate(NN4x1)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN3x3)*complexconjugate(NN4x2)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN3x3)*complexconjugate(NN4x2)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN3x1)*complexconjugate(NN4x3)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN3x2)*complexconjugate(NN4x3)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN3x2)*complexconjugate(NN4x3)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN3x4)*complexconjugate(NN4x1)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN3x4)*complexconjugate(NN4x2)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN3x4)*complexconjugate(NN4x2)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN3x1)*complexconjugate(NN4x4)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN3x2)*complexconjugate(NN4x4)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN3x2)*complexconjugate(NN4x4)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1731 = Coupling(name = 'GC_1731',
                   value = '(cw*ee*complexconjugate(NN4x1)*complexconjugate(NN4x3)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN4x2)*complexconjugate(NN4x3)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN4x2)*complexconjugate(NN4x3)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complexconjugate(NN4x1)*complexconjugate(NN4x4)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee*complexconjugate(NN4x2)*complexconjugate(NN4x4)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*complexconjugate(NN4x2)*complexconjugate(NN4x4)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1732 = Coupling(name = 'GC_1732',
                   value = '(ee*complex(0,1)*complexconjugate(NN1x3)*complexconjugate(UU1x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN1x3)*complexconjugate(UU1x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1733 = Coupling(name = 'GC_1733',
                   value = '(ee*complex(0,1)*complexconjugate(NN2x3)*complexconjugate(UU1x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN2x3)*complexconjugate(UU1x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1734 = Coupling(name = 'GC_1734',
                   value = '(ee*complex(0,1)*complexconjugate(NN3x3)*complexconjugate(UU1x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN3x3)*complexconjugate(UU1x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*complexconjugate(NN3x2)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*sw*complexconjugate(NN3x2)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1735 = Coupling(name = 'GC_1735',
                   value = '(ee*complex(0,1)*complexconjugate(NN4x3)*complexconjugate(UU1x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN4x3)*complexconjugate(UU1x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*complexconjugate(NN4x2)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*sw*complexconjugate(NN4x2)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1736 = Coupling(name = 'GC_1736',
                   value = '(ee*complex(0,1)*complexconjugate(NN1x3)*complexconjugate(UU2x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN1x3)*complexconjugate(UU2x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1737 = Coupling(name = 'GC_1737',
                   value = '(ee*complex(0,1)*complexconjugate(NN2x3)*complexconjugate(UU2x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN2x3)*complexconjugate(UU2x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1738 = Coupling(name = 'GC_1738',
                   value = '(ee*complex(0,1)*complexconjugate(NN3x3)*complexconjugate(UU2x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN3x3)*complexconjugate(UU2x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*complexconjugate(NN3x2)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*sw*complexconjugate(NN3x2)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1739 = Coupling(name = 'GC_1739',
                   value = '(ee*complex(0,1)*complexconjugate(NN4x3)*complexconjugate(UU2x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN4x3)*complexconjugate(UU2x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*complexconjugate(NN4x2)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*sw*complexconjugate(NN4x2)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1740 = Coupling(name = 'GC_1740',
                   value = '-((ee*complexconjugate(UU1x1)*complexconjugate(VV1x2)*cmath.cos(beta))/(sw*cmath.sqrt(2))) - (ee*complexconjugate(UU1x2)*complexconjugate(VV1x1)*cmath.sin(beta))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1741 = Coupling(name = 'GC_1741',
                   value = '-((ee*complexconjugate(UU2x1)*complexconjugate(VV1x2)*cmath.cos(beta))/(sw*cmath.sqrt(2))) - (ee*complexconjugate(UU2x2)*complexconjugate(VV1x1)*cmath.sin(beta))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1742 = Coupling(name = 'GC_1742',
                   value = '(ee*complex(0,1)*complexconjugate(NN1x4)*complexconjugate(VV1x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN1x4)*complexconjugate(VV1x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1743 = Coupling(name = 'GC_1743',
                   value = '(ee*complex(0,1)*complexconjugate(NN2x4)*complexconjugate(VV1x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN2x4)*complexconjugate(VV1x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1744 = Coupling(name = 'GC_1744',
                   value = '(ee*complex(0,1)*complexconjugate(NN3x4)*complexconjugate(VV1x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN3x4)*complexconjugate(VV1x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN3x2)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN3x2)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1745 = Coupling(name = 'GC_1745',
                   value = '(ee*complex(0,1)*complexconjugate(NN4x4)*complexconjugate(VV1x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN4x4)*complexconjugate(VV1x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN4x2)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN4x2)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1746 = Coupling(name = 'GC_1746',
                   value = '(ee*complexconjugate(UU1x2)*complexconjugate(VV1x1)*cmath.cos(beta))/(sw*cmath.sqrt(2)) - (ee*complexconjugate(UU1x1)*complexconjugate(VV1x2)*cmath.sin(beta))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1747 = Coupling(name = 'GC_1747',
                   value = '(ee*complexconjugate(UU2x2)*complexconjugate(VV1x1)*cmath.cos(beta))/(sw*cmath.sqrt(2)) - (ee*complexconjugate(UU2x1)*complexconjugate(VV1x2)*cmath.sin(beta))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1748 = Coupling(name = 'GC_1748',
                   value = '-((ee*complexconjugate(UU1x1)*complexconjugate(VV2x2)*cmath.cos(beta))/(sw*cmath.sqrt(2))) - (ee*complexconjugate(UU1x2)*complexconjugate(VV2x1)*cmath.sin(beta))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1749 = Coupling(name = 'GC_1749',
                   value = '-((ee*complexconjugate(UU2x1)*complexconjugate(VV2x2)*cmath.cos(beta))/(sw*cmath.sqrt(2))) - (ee*complexconjugate(UU2x2)*complexconjugate(VV2x1)*cmath.sin(beta))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1750 = Coupling(name = 'GC_1750',
                   value = '(ee*complex(0,1)*complexconjugate(NN1x4)*complexconjugate(VV2x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN1x4)*complexconjugate(VV2x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1751 = Coupling(name = 'GC_1751',
                   value = '(ee*complex(0,1)*complexconjugate(NN2x4)*complexconjugate(VV2x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN2x4)*complexconjugate(VV2x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1752 = Coupling(name = 'GC_1752',
                   value = '(ee*complex(0,1)*complexconjugate(NN3x4)*complexconjugate(VV2x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN3x4)*complexconjugate(VV2x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN3x2)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN3x2)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1753 = Coupling(name = 'GC_1753',
                   value = '(ee*complex(0,1)*complexconjugate(NN4x4)*complexconjugate(VV2x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN4x4)*complexconjugate(VV2x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN4x2)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN4x2)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1754 = Coupling(name = 'GC_1754',
                   value = '(ee*complexconjugate(UU1x2)*complexconjugate(VV2x1)*cmath.cos(beta))/(sw*cmath.sqrt(2)) - (ee*complexconjugate(UU1x1)*complexconjugate(VV2x2)*cmath.sin(beta))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1755 = Coupling(name = 'GC_1755',
                   value = '(ee*complexconjugate(UU2x2)*complexconjugate(VV2x1)*cmath.cos(beta))/(sw*cmath.sqrt(2)) - (ee*complexconjugate(UU2x1)*complexconjugate(VV2x2)*cmath.sin(beta))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1756 = Coupling(name = 'GC_1756',
                   value = '-(ee**2*complex(0,1)*I109a11*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I109a11*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1757 = Coupling(name = 'GC_1757',
                   value = '-(ee**2*complex(0,1)*I109a22*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I109a22*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1758 = Coupling(name = 'GC_1758',
                   value = '(complex(0,1)*I110a33*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) - (ee**2*complex(0,1)*I109a33*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I109a33*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1759 = Coupling(name = 'GC_1759',
                   value = '-(ee**2*complex(0,1)*I109a33*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I110a33*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I109a33*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1760 = Coupling(name = 'GC_1760',
                   value = '(complex(0,1)*I110a36*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) - (ee**2*complex(0,1)*I109a36*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I109a36*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1761 = Coupling(name = 'GC_1761',
                   value = '-(ee**2*complex(0,1)*I109a36*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I110a36*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I109a36*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1762 = Coupling(name = 'GC_1762',
                   value = '-(ee**2*complex(0,1)*I123a11*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I123a11*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1763 = Coupling(name = 'GC_1763',
                   value = '-(ee**2*complex(0,1)*I123a22*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I123a22*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1764 = Coupling(name = 'GC_1764',
                   value = '-(ee**2*complex(0,1)*I144a11*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I144a11*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1765 = Coupling(name = 'GC_1765',
                   value = '-(ee**2*complex(0,1)*I144a22*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I144a22*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1766 = Coupling(name = 'GC_1766',
                   value = '-(ee**2*complex(0,1)*I94a11*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I94a11*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1767 = Coupling(name = 'GC_1767',
                   value = '-(ee**2*complex(0,1)*I94a22*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I94a22*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1768 = Coupling(name = 'GC_1768',
                   value = '(complex(0,1)*I95a33*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) - (ee**2*complex(0,1)*I94a33*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I94a33*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1769 = Coupling(name = 'GC_1769',
                   value = '-(ee**2*complex(0,1)*I94a33*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I95a33*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I94a33*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1770 = Coupling(name = 'GC_1770',
                   value = '(complex(0,1)*I95a36*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) - (ee**2*complex(0,1)*I94a36*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I94a36*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1771 = Coupling(name = 'GC_1771',
                   value = '-(ee**2*complex(0,1)*I94a36*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I95a36*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I94a36*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1772 = Coupling(name = 'GC_1772',
                   value = '(ee*complex(0,1)*cmath.cos(beta)*cmath.sin(alp))/(2.*sw) - (ee*complex(0,1)*cmath.cos(alp)*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':1})

GC_1773 = Coupling(name = 'GC_1773',
                   value = '-(ee*complex(0,1)*cmath.cos(beta)*cmath.sin(alp))/(2.*sw) + (ee*complex(0,1)*cmath.cos(alp)*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':1})

GC_1774 = Coupling(name = 'GC_1774',
                   value = '(ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(alp))/(2.*sw) - (ee**2*complex(0,1)*cmath.cos(alp)*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':2})

GC_1775 = Coupling(name = 'GC_1775',
                   value = '-(ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(alp))/(2.*sw) + (ee**2*complex(0,1)*cmath.cos(alp)*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':2})

GC_1776 = Coupling(name = 'GC_1776',
                   value = '(cw*ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee**2*complex(0,1)*cmath.cos(alp)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1777 = Coupling(name = 'GC_1777',
                   value = '-(cw*ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee**2*complex(0,1)*cmath.cos(alp)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1778 = Coupling(name = 'GC_1778',
                   value = '(cw*ee*cmath.cos(beta)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*cmath.cos(alp)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw))',
                   order = {'QED':1})

GC_1779 = Coupling(name = 'GC_1779',
                   value = '-(cw*ee*cmath.cos(beta)*cmath.sin(alp))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*cmath.cos(alp)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw))',
                   order = {'QED':1})

GC_1780 = Coupling(name = 'GC_1780',
                   value = '(I110a33*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*I109a33*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1781 = Coupling(name = 'GC_1781',
                   value = '(I110a36*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*I109a36*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1782 = Coupling(name = 'GC_1782',
                   value = '-((I95a33*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2)) + (ee**2*I94a33*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1783 = Coupling(name = 'GC_1783',
                   value = '-((I95a36*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2)) + (ee**2*I94a36*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1784 = Coupling(name = 'GC_1784',
                   value = '-((ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))) + (ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1785 = Coupling(name = 'GC_1785',
                   value = '(2*ee**2*complex(0,1)*I148a11*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a11*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1786 = Coupling(name = 'GC_1786',
                   value = '-(ee**2*complex(0,1)*I148a11*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a11*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1787 = Coupling(name = 'GC_1787',
                   value = '(2*ee**2*complex(0,1)*I148a22*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a22*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1788 = Coupling(name = 'GC_1788',
                   value = '-(ee**2*complex(0,1)*I148a22*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a22*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1789 = Coupling(name = 'GC_1789',
                   value = '(2*ee**2*complex(0,1)*I15a11*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a11*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1790 = Coupling(name = 'GC_1790',
                   value = '-(ee**2*complex(0,1)*I15a11*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a11*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1791 = Coupling(name = 'GC_1791',
                   value = '(2*ee**2*complex(0,1)*I15a22*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a22*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1792 = Coupling(name = 'GC_1792',
                   value = '-(ee**2*complex(0,1)*I15a22*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a22*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1793 = Coupling(name = 'GC_1793',
                   value = '-((ee**2*complex(0,1)*I53a11*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))) + (ee**2*complex(0,1)*I53a11*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1794 = Coupling(name = 'GC_1794',
                   value = '-((ee**2*complex(0,1)*I53a22*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))) + (ee**2*complex(0,1)*I53a22*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1795 = Coupling(name = 'GC_1795',
                   value = '(2*ee**2*complex(0,1)*I148a33*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I149a33*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a33*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a33*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a33*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a33*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a33*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1796 = Coupling(name = 'GC_1796',
                   value = '(2*ee**2*complex(0,1)*I148a36*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I149a36*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a36*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a36*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a36*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a36*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a36*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1797 = Coupling(name = 'GC_1797',
                   value = '(2*ee**2*complex(0,1)*I148a63*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I149a63*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a63*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a63*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a63*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a63*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a63*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1798 = Coupling(name = 'GC_1798',
                   value = '(2*ee**2*complex(0,1)*I148a66*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I149a66*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a66*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a66*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a66*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a66*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a66*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1799 = Coupling(name = 'GC_1799',
                   value = '-(ee**2*complex(0,1)*I148a33*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I149a33*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I157a33*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I158a33*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a33*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I157a33*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I158a33*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1800 = Coupling(name = 'GC_1800',
                   value = '-(ee**2*complex(0,1)*I148a36*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I149a36*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I157a36*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I158a36*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a36*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I157a36*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I158a36*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1801 = Coupling(name = 'GC_1801',
                   value = '-(ee**2*complex(0,1)*I148a63*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I149a63*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I157a63*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I158a63*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a63*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I157a63*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I158a63*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1802 = Coupling(name = 'GC_1802',
                   value = '-(ee**2*complex(0,1)*I148a66*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I149a66*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I157a66*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I158a66*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a66*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I157a66*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I158a66*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1803 = Coupling(name = 'GC_1803',
                   value = '-(ee**2*complex(0,1)*I15a33*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a33*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I23a33*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a33*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a33*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I23a33*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a33*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1804 = Coupling(name = 'GC_1804',
                   value = '-(ee**2*complex(0,1)*I15a36*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a36*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I23a36*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a36*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a36*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I23a36*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a36*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1805 = Coupling(name = 'GC_1805',
                   value = '-(ee**2*complex(0,1)*I15a63*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a63*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I23a63*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a63*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a63*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I23a63*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a63*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1806 = Coupling(name = 'GC_1806',
                   value = '-(ee**2*complex(0,1)*I15a66*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a66*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I23a66*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a66*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a66*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I23a66*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a66*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1807 = Coupling(name = 'GC_1807',
                   value = '(2*ee**2*complex(0,1)*I15a33*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a33*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I24a33*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I25a33*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a33*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I24a33*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I25a33*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1808 = Coupling(name = 'GC_1808',
                   value = '(2*ee**2*complex(0,1)*I15a36*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a36*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I24a36*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I25a36*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a36*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I24a36*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I25a36*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1809 = Coupling(name = 'GC_1809',
                   value = '(2*ee**2*complex(0,1)*I15a63*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a63*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I24a63*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I25a63*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a63*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I24a63*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I25a63*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1810 = Coupling(name = 'GC_1810',
                   value = '(2*ee**2*complex(0,1)*I15a66*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a66*cmath.cos(beta)*cmath.sin(beta))/(3.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I24a66*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I25a66*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a66*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I24a66*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I25a66*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1811 = Coupling(name = 'GC_1811',
                   value = '(ee**2*complex(0,1)*I54a33*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a33*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a33*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I62a33*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1812 = Coupling(name = 'GC_1812',
                   value = '-((ee**2*complex(0,1)*I53a33*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))) + (ee**2*complex(0,1)*I54a33*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I61a33*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a33*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a33*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I61a33*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a33*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1813 = Coupling(name = 'GC_1813',
                   value = '(ee**2*complex(0,1)*I54a36*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a36*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a36*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I62a36*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1814 = Coupling(name = 'GC_1814',
                   value = '-((ee**2*complex(0,1)*I53a36*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))) + (ee**2*complex(0,1)*I54a36*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I61a36*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a36*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a36*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I61a36*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a36*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1815 = Coupling(name = 'GC_1815',
                   value = '(ee**2*complex(0,1)*I54a63*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a63*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a63*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I62a63*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1816 = Coupling(name = 'GC_1816',
                   value = '-((ee**2*complex(0,1)*I53a63*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))) + (ee**2*complex(0,1)*I54a63*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I61a63*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a63*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a63*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I61a63*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a63*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1817 = Coupling(name = 'GC_1817',
                   value = '(ee**2*complex(0,1)*I54a66*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a66*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a66*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I62a66*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1818 = Coupling(name = 'GC_1818',
                   value = '-((ee**2*complex(0,1)*I53a66*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))) + (ee**2*complex(0,1)*I54a66*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I61a66*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a66*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a66*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I61a66*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a66*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1819 = Coupling(name = 'GC_1819',
                   value = '-((ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))) - (complex(0,1)*I96a33*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*I96a33*sw**2*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1820 = Coupling(name = 'GC_1820',
                   value = '-((complex(0,1)*I125a33*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2)) + (complex(0,1)*I124a33*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a33*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I126a33*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a33*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I125a33*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1821 = Coupling(name = 'GC_1821',
                   value = '(complex(0,1)*I125a33*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I126a33*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a33*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I124a33*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a33*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I125a33*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1822 = Coupling(name = 'GC_1822',
                   value = '-((complex(0,1)*I125a36*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2)) + (complex(0,1)*I124a36*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a36*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I126a36*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a36*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I125a36*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1823 = Coupling(name = 'GC_1823',
                   value = '(complex(0,1)*I125a36*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I126a36*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a36*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I124a36*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a36*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I125a36*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1824 = Coupling(name = 'GC_1824',
                   value = '-((complex(0,1)*I125a63*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2)) + (complex(0,1)*I124a63*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a63*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I126a63*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a63*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I125a63*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1825 = Coupling(name = 'GC_1825',
                   value = '(complex(0,1)*I125a63*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I126a63*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a63*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I124a63*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a63*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I125a63*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1826 = Coupling(name = 'GC_1826',
                   value = '-((complex(0,1)*I125a66*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2)) + (complex(0,1)*I124a66*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a66*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I126a66*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a66*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I125a66*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1827 = Coupling(name = 'GC_1827',
                   value = '(complex(0,1)*I125a66*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I126a66*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a66*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I124a66*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a66*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I125a66*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1828 = Coupling(name = 'GC_1828',
                   value = '-((complex(0,1)*I147a33*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2)) + (complex(0,1)*I145a33*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a33*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I146a33*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a33*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I147a33*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1829 = Coupling(name = 'GC_1829',
                   value = '(complex(0,1)*I147a33*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I146a33*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a33*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I145a33*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a33*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I147a33*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1830 = Coupling(name = 'GC_1830',
                   value = '-((complex(0,1)*I147a36*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2)) + (complex(0,1)*I145a36*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a36*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I146a36*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a36*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I147a36*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1831 = Coupling(name = 'GC_1831',
                   value = '(complex(0,1)*I147a36*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I146a36*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a36*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I145a36*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a36*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I147a36*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1832 = Coupling(name = 'GC_1832',
                   value = '-((complex(0,1)*I147a63*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2)) + (complex(0,1)*I145a63*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a63*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I146a63*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a63*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I147a63*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1833 = Coupling(name = 'GC_1833',
                   value = '(complex(0,1)*I147a63*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I146a63*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a63*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I145a63*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a63*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I147a63*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1834 = Coupling(name = 'GC_1834',
                   value = '-((complex(0,1)*I147a66*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2)) + (complex(0,1)*I145a66*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a66*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I146a66*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a66*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I147a66*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1835 = Coupling(name = 'GC_1835',
                   value = '(complex(0,1)*I147a66*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I146a66*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a66*cmath.cos(beta)*cmath.sin(alp))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I145a66*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a66*cmath.cos(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I147a66*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1836 = Coupling(name = 'GC_1836',
                   value = '(ee**2*complex(0,1)*I109a11*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I109a11*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1837 = Coupling(name = 'GC_1837',
                   value = '-(ee**2*complex(0,1)*I109a11*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + (ee**2*complex(0,1)*I109a11*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1838 = Coupling(name = 'GC_1838',
                   value = '(ee**2*complex(0,1)*I109a22*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I109a22*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1839 = Coupling(name = 'GC_1839',
                   value = '-(ee**2*complex(0,1)*I109a22*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + (ee**2*complex(0,1)*I109a22*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1840 = Coupling(name = 'GC_1840',
                   value = '-((complex(0,1)*I110a33*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2)) + (ee**2*complex(0,1)*I109a33*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I109a33*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1841 = Coupling(name = 'GC_1841',
                   value = '-(ee**2*complex(0,1)*I109a33*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I110a33*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I109a33*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1842 = Coupling(name = 'GC_1842',
                   value = '-((complex(0,1)*I110a36*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2)) + (ee**2*complex(0,1)*I109a36*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I109a36*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1843 = Coupling(name = 'GC_1843',
                   value = '-(ee**2*complex(0,1)*I109a36*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I110a36*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I109a36*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1844 = Coupling(name = 'GC_1844',
                   value = '(ee**2*complex(0,1)*I123a11*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I123a11*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1845 = Coupling(name = 'GC_1845',
                   value = '-(ee**2*complex(0,1)*I123a11*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + (ee**2*complex(0,1)*I123a11*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1846 = Coupling(name = 'GC_1846',
                   value = '(ee**2*complex(0,1)*I123a22*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I123a22*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1847 = Coupling(name = 'GC_1847',
                   value = '-(ee**2*complex(0,1)*I123a22*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + (ee**2*complex(0,1)*I123a22*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1848 = Coupling(name = 'GC_1848',
                   value = '-((complex(0,1)*I124a33*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2)) + (ee**2*complex(0,1)*I123a33*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I125a33*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) + (complex(0,1)*I125a33*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) + (complex(0,1)*I126a33*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a33*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1849 = Coupling(name = 'GC_1849',
                   value = '(complex(0,1)*I126a33*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a33*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I125a33*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) + (complex(0,1)*I125a33*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (complex(0,1)*I124a33*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I123a33*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1850 = Coupling(name = 'GC_1850',
                   value = '-((complex(0,1)*I124a36*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2)) + (ee**2*complex(0,1)*I123a36*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I125a36*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) + (complex(0,1)*I125a36*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) + (complex(0,1)*I126a36*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a36*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1851 = Coupling(name = 'GC_1851',
                   value = '(complex(0,1)*I126a36*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a36*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I125a36*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) + (complex(0,1)*I125a36*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (complex(0,1)*I124a36*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I123a36*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1852 = Coupling(name = 'GC_1852',
                   value = '-((complex(0,1)*I124a63*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2)) + (ee**2*complex(0,1)*I123a63*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I125a63*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) + (complex(0,1)*I125a63*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) + (complex(0,1)*I126a63*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a63*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1853 = Coupling(name = 'GC_1853',
                   value = '(complex(0,1)*I126a63*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a63*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I125a63*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) + (complex(0,1)*I125a63*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (complex(0,1)*I124a63*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I123a63*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1854 = Coupling(name = 'GC_1854',
                   value = '-((complex(0,1)*I124a66*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2)) + (ee**2*complex(0,1)*I123a66*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I125a66*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) + (complex(0,1)*I125a66*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) + (complex(0,1)*I126a66*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a66*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1855 = Coupling(name = 'GC_1855',
                   value = '(complex(0,1)*I126a66*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I123a66*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I125a66*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) + (complex(0,1)*I125a66*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (complex(0,1)*I124a66*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I123a66*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1856 = Coupling(name = 'GC_1856',
                   value = '(ee**2*complex(0,1)*I144a11*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I144a11*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1857 = Coupling(name = 'GC_1857',
                   value = '-(ee**2*complex(0,1)*I144a11*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + (ee**2*complex(0,1)*I144a11*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1858 = Coupling(name = 'GC_1858',
                   value = '(ee**2*complex(0,1)*I144a22*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I144a22*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1859 = Coupling(name = 'GC_1859',
                   value = '-(ee**2*complex(0,1)*I144a22*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + (ee**2*complex(0,1)*I144a22*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1860 = Coupling(name = 'GC_1860',
                   value = '-((complex(0,1)*I145a33*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2)) + (ee**2*complex(0,1)*I144a33*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I147a33*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) + (complex(0,1)*I147a33*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) + (complex(0,1)*I146a33*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a33*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1861 = Coupling(name = 'GC_1861',
                   value = '(complex(0,1)*I146a33*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a33*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I147a33*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) + (complex(0,1)*I147a33*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (complex(0,1)*I145a33*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I144a33*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1862 = Coupling(name = 'GC_1862',
                   value = '-((complex(0,1)*I145a36*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2)) + (ee**2*complex(0,1)*I144a36*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I147a36*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) + (complex(0,1)*I147a36*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) + (complex(0,1)*I146a36*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a36*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1863 = Coupling(name = 'GC_1863',
                   value = '(complex(0,1)*I146a36*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a36*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I147a36*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) + (complex(0,1)*I147a36*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (complex(0,1)*I145a36*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I144a36*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1864 = Coupling(name = 'GC_1864',
                   value = '-((complex(0,1)*I145a63*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2)) + (ee**2*complex(0,1)*I144a63*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I147a63*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) + (complex(0,1)*I147a63*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) + (complex(0,1)*I146a63*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a63*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1865 = Coupling(name = 'GC_1865',
                   value = '(complex(0,1)*I146a63*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a63*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I147a63*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) + (complex(0,1)*I147a63*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (complex(0,1)*I145a63*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I144a63*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1866 = Coupling(name = 'GC_1866',
                   value = '-((complex(0,1)*I145a66*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2)) + (ee**2*complex(0,1)*I144a66*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I147a66*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) + (complex(0,1)*I147a66*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) + (complex(0,1)*I146a66*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a66*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1867 = Coupling(name = 'GC_1867',
                   value = '(complex(0,1)*I146a66*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I144a66*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I147a66*cmath.cos(beta)*cmath.sin(alp))/cmath.sqrt(2) + (complex(0,1)*I147a66*cmath.cos(alp)*cmath.sin(beta))/cmath.sqrt(2) - (complex(0,1)*I145a66*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I144a66*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1868 = Coupling(name = 'GC_1868',
                   value = '(ee**2*complex(0,1)*I94a11*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I94a11*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1869 = Coupling(name = 'GC_1869',
                   value = '-(ee**2*complex(0,1)*I94a11*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + (ee**2*complex(0,1)*I94a11*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1870 = Coupling(name = 'GC_1870',
                   value = '(ee**2*complex(0,1)*I94a22*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I94a22*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1871 = Coupling(name = 'GC_1871',
                   value = '-(ee**2*complex(0,1)*I94a22*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + (ee**2*complex(0,1)*I94a22*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1872 = Coupling(name = 'GC_1872',
                   value = '-((complex(0,1)*I95a33*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2)) + (ee**2*complex(0,1)*I94a33*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I94a33*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1873 = Coupling(name = 'GC_1873',
                   value = '-(ee**2*complex(0,1)*I94a33*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I95a33*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I94a33*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1874 = Coupling(name = 'GC_1874',
                   value = '-((complex(0,1)*I95a36*cmath.cos(alp)*cmath.cos(beta))/cmath.sqrt(2)) + (ee**2*complex(0,1)*I94a36*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I94a36*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1875 = Coupling(name = 'GC_1875',
                   value = '-(ee**2*complex(0,1)*I94a36*cmath.cos(alp)*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (complex(0,1)*I95a36*cmath.sin(alp)*cmath.sin(beta))/cmath.sqrt(2) + (ee**2*complex(0,1)*I94a36*cmath.sin(alp)*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1876 = Coupling(name = 'GC_1876',
                   value = '-(ee*complex(0,1)*cmath.cos(alp)*cmath.cos(beta))/(2.*sw) - (ee*complex(0,1)*cmath.sin(alp)*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':1})

GC_1877 = Coupling(name = 'GC_1877',
                   value = '(ee*complex(0,1)*cmath.cos(alp)*cmath.cos(beta))/(2.*sw) + (ee*complex(0,1)*cmath.sin(alp)*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':1})

GC_1878 = Coupling(name = 'GC_1878',
                   value = '(ee**2*complex(0,1)*cmath.cos(alp)*cmath.cos(beta))/(2.*sw) + (ee**2*complex(0,1)*cmath.sin(alp)*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':2})

GC_1879 = Coupling(name = 'GC_1879',
                   value = '(cw*ee**2*complex(0,1)*cmath.cos(alp)*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee**2*complex(0,1)*cmath.sin(alp)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1880 = Coupling(name = 'GC_1880',
                   value = '(cw*ee*cmath.cos(alp)*cmath.cos(beta))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*cmath.sin(alp)*cmath.sin(beta))/(2.*(-1 + sw)*sw*(1 + sw))',
                   order = {'QED':1})

GC_1881 = Coupling(name = 'GC_1881',
                   value = '(ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vd*cmath.cos(beta)*cmath.sin(alp)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1882 = Coupling(name = 'GC_1882',
                   value = '-(ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vu*cmath.cos(beta)*cmath.sin(alp)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1883 = Coupling(name = 'GC_1883',
                   value = '(ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(alp)**2*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1884 = Coupling(name = 'GC_1884',
                   value = '-(ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(alp)**2*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1885 = Coupling(name = 'GC_1885',
                   value = '-(ee*complex(0,1)*cmath.cos(beta)**2) - ee*complex(0,1)*cmath.sin(beta)**2',
                   order = {'QED':1})

GC_1886 = Coupling(name = 'GC_1886',
                   value = '2*ee**2*complex(0,1)*cmath.cos(beta)**2 + 2*ee**2*complex(0,1)*cmath.sin(beta)**2',
                   order = {'QED':2})

GC_1887 = Coupling(name = 'GC_1887',
                   value = '-((I125a33*cmath.cos(beta)**2)/cmath.sqrt(2)) - (I124a33*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) - (I126a33*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) + (ee**2*I123a33*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2)) - (I125a33*cmath.sin(beta)**2)/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1888 = Coupling(name = 'GC_1888',
                   value = '(I125a33*cmath.cos(beta)**2)/cmath.sqrt(2) - (I124a33*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) - (I126a33*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) + (ee**2*I123a33*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2)) + (I125a33*cmath.sin(beta)**2)/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1889 = Coupling(name = 'GC_1889',
                   value = '-((I125a36*cmath.cos(beta)**2)/cmath.sqrt(2)) - (I124a36*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) - (I126a36*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) + (ee**2*I123a36*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2)) - (I125a36*cmath.sin(beta)**2)/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1890 = Coupling(name = 'GC_1890',
                   value = '(I125a36*cmath.cos(beta)**2)/cmath.sqrt(2) - (I124a36*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) - (I126a36*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) + (ee**2*I123a36*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2)) + (I125a36*cmath.sin(beta)**2)/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1891 = Coupling(name = 'GC_1891',
                   value = '-((I125a63*cmath.cos(beta)**2)/cmath.sqrt(2)) - (I124a63*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) - (I126a63*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) + (ee**2*I123a63*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2)) - (I125a63*cmath.sin(beta)**2)/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1892 = Coupling(name = 'GC_1892',
                   value = '(I125a63*cmath.cos(beta)**2)/cmath.sqrt(2) - (I124a63*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) - (I126a63*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) + (ee**2*I123a63*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2)) + (I125a63*cmath.sin(beta)**2)/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1893 = Coupling(name = 'GC_1893',
                   value = '-((I125a66*cmath.cos(beta)**2)/cmath.sqrt(2)) - (I124a66*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) - (I126a66*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) + (ee**2*I123a66*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2)) - (I125a66*cmath.sin(beta)**2)/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1894 = Coupling(name = 'GC_1894',
                   value = '(I125a66*cmath.cos(beta)**2)/cmath.sqrt(2) - (I124a66*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) - (I126a66*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) + (ee**2*I123a66*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2)) + (I125a66*cmath.sin(beta)**2)/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1895 = Coupling(name = 'GC_1895',
                   value = '-((I147a33*cmath.cos(beta)**2)/cmath.sqrt(2)) + (I145a33*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) + (I146a33*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*I144a33*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2)) - (I147a33*cmath.sin(beta)**2)/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1896 = Coupling(name = 'GC_1896',
                   value = '(I147a33*cmath.cos(beta)**2)/cmath.sqrt(2) + (I145a33*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) + (I146a33*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*I144a33*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2)) + (I147a33*cmath.sin(beta)**2)/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1897 = Coupling(name = 'GC_1897',
                   value = '-((I147a36*cmath.cos(beta)**2)/cmath.sqrt(2)) + (I145a36*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) + (I146a36*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*I144a36*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2)) - (I147a36*cmath.sin(beta)**2)/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1898 = Coupling(name = 'GC_1898',
                   value = '(I147a36*cmath.cos(beta)**2)/cmath.sqrt(2) + (I145a36*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) + (I146a36*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*I144a36*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2)) + (I147a36*cmath.sin(beta)**2)/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1899 = Coupling(name = 'GC_1899',
                   value = '-((I147a63*cmath.cos(beta)**2)/cmath.sqrt(2)) + (I145a63*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) + (I146a63*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*I144a63*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2)) - (I147a63*cmath.sin(beta)**2)/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1900 = Coupling(name = 'GC_1900',
                   value = '(I147a63*cmath.cos(beta)**2)/cmath.sqrt(2) + (I145a63*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) + (I146a63*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*I144a63*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2)) + (I147a63*cmath.sin(beta)**2)/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1901 = Coupling(name = 'GC_1901',
                   value = '-((I147a66*cmath.cos(beta)**2)/cmath.sqrt(2)) + (I145a66*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) + (I146a66*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*I144a66*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2)) - (I147a66*cmath.sin(beta)**2)/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1902 = Coupling(name = 'GC_1902',
                   value = '(I147a66*cmath.cos(beta)**2)/cmath.sqrt(2) + (I145a66*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) + (I146a66*cmath.cos(beta)*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*I144a66*cmath.cos(beta)*cmath.sin(beta))/(sw**2*cmath.sqrt(2)) + (I147a66*cmath.sin(beta)**2)/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1903 = Coupling(name = 'GC_1903',
                   value = '(ee**2*complex(0,1)*cmath.cos(beta)**2)/(2.*sw**2) + (ee**2*complex(0,1)*cmath.sin(beta)**2)/(2.*sw**2)',
                   order = {'QED':2})

GC_1904 = Coupling(name = 'GC_1904',
                   value = '(ee**2*I109a11*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) - (ee**2*I109a11*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1905 = Coupling(name = 'GC_1905',
                   value = '-(ee**2*I109a11*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) + (ee**2*I109a11*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1906 = Coupling(name = 'GC_1906',
                   value = '(ee**2*I109a22*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) - (ee**2*I109a22*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1907 = Coupling(name = 'GC_1907',
                   value = '-(ee**2*I109a22*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) + (ee**2*I109a22*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1908 = Coupling(name = 'GC_1908',
                   value = '-((I110a33*cmath.cos(beta)**2)/cmath.sqrt(2)) + (ee**2*I109a33*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) - (ee**2*I109a33*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1909 = Coupling(name = 'GC_1909',
                   value = '-(ee**2*I109a33*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) - (I110a33*cmath.sin(beta)**2)/cmath.sqrt(2) + (ee**2*I109a33*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1910 = Coupling(name = 'GC_1910',
                   value = '-((I110a36*cmath.cos(beta)**2)/cmath.sqrt(2)) + (ee**2*I109a36*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) - (ee**2*I109a36*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1911 = Coupling(name = 'GC_1911',
                   value = '-(ee**2*I109a36*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) - (I110a36*cmath.sin(beta)**2)/cmath.sqrt(2) + (ee**2*I109a36*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1912 = Coupling(name = 'GC_1912',
                   value = '(ee**2*I123a11*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) - (ee**2*I123a11*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1913 = Coupling(name = 'GC_1913',
                   value = '-(ee**2*I123a11*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) + (ee**2*I123a11*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1914 = Coupling(name = 'GC_1914',
                   value = '(ee**2*I123a22*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) - (ee**2*I123a22*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1915 = Coupling(name = 'GC_1915',
                   value = '-(ee**2*I123a22*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) + (ee**2*I123a22*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1916 = Coupling(name = 'GC_1916',
                   value = '-((I126a33*cmath.cos(beta)**2)/cmath.sqrt(2)) + (ee**2*I123a33*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) + (I124a33*cmath.sin(beta)**2)/cmath.sqrt(2) - (ee**2*I123a33*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1917 = Coupling(name = 'GC_1917',
                   value = '(I124a33*cmath.cos(beta)**2)/cmath.sqrt(2) - (ee**2*I123a33*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) - (I126a33*cmath.sin(beta)**2)/cmath.sqrt(2) + (ee**2*I123a33*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1918 = Coupling(name = 'GC_1918',
                   value = '-((I126a36*cmath.cos(beta)**2)/cmath.sqrt(2)) + (ee**2*I123a36*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) + (I124a36*cmath.sin(beta)**2)/cmath.sqrt(2) - (ee**2*I123a36*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1919 = Coupling(name = 'GC_1919',
                   value = '(I124a36*cmath.cos(beta)**2)/cmath.sqrt(2) - (ee**2*I123a36*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) - (I126a36*cmath.sin(beta)**2)/cmath.sqrt(2) + (ee**2*I123a36*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1920 = Coupling(name = 'GC_1920',
                   value = '-((I126a63*cmath.cos(beta)**2)/cmath.sqrt(2)) + (ee**2*I123a63*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) + (I124a63*cmath.sin(beta)**2)/cmath.sqrt(2) - (ee**2*I123a63*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1921 = Coupling(name = 'GC_1921',
                   value = '(I124a63*cmath.cos(beta)**2)/cmath.sqrt(2) - (ee**2*I123a63*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) - (I126a63*cmath.sin(beta)**2)/cmath.sqrt(2) + (ee**2*I123a63*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1922 = Coupling(name = 'GC_1922',
                   value = '-((I126a66*cmath.cos(beta)**2)/cmath.sqrt(2)) + (ee**2*I123a66*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) + (I124a66*cmath.sin(beta)**2)/cmath.sqrt(2) - (ee**2*I123a66*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1923 = Coupling(name = 'GC_1923',
                   value = '(I124a66*cmath.cos(beta)**2)/cmath.sqrt(2) - (ee**2*I123a66*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) - (I126a66*cmath.sin(beta)**2)/cmath.sqrt(2) + (ee**2*I123a66*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1924 = Coupling(name = 'GC_1924',
                   value = '(ee**2*I144a11*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) - (ee**2*I144a11*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1925 = Coupling(name = 'GC_1925',
                   value = '-(ee**2*I144a11*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) + (ee**2*I144a11*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1926 = Coupling(name = 'GC_1926',
                   value = '(ee**2*I144a22*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) - (ee**2*I144a22*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1927 = Coupling(name = 'GC_1927',
                   value = '-(ee**2*I144a22*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) + (ee**2*I144a22*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1928 = Coupling(name = 'GC_1928',
                   value = '-((I145a33*cmath.cos(beta)**2)/cmath.sqrt(2)) + (ee**2*I144a33*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) + (I146a33*cmath.sin(beta)**2)/cmath.sqrt(2) - (ee**2*I144a33*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1929 = Coupling(name = 'GC_1929',
                   value = '(I146a33*cmath.cos(beta)**2)/cmath.sqrt(2) - (ee**2*I144a33*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) - (I145a33*cmath.sin(beta)**2)/cmath.sqrt(2) + (ee**2*I144a33*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1930 = Coupling(name = 'GC_1930',
                   value = '-((I145a36*cmath.cos(beta)**2)/cmath.sqrt(2)) + (ee**2*I144a36*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) + (I146a36*cmath.sin(beta)**2)/cmath.sqrt(2) - (ee**2*I144a36*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1931 = Coupling(name = 'GC_1931',
                   value = '(I146a36*cmath.cos(beta)**2)/cmath.sqrt(2) - (ee**2*I144a36*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) - (I145a36*cmath.sin(beta)**2)/cmath.sqrt(2) + (ee**2*I144a36*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1932 = Coupling(name = 'GC_1932',
                   value = '-((I145a63*cmath.cos(beta)**2)/cmath.sqrt(2)) + (ee**2*I144a63*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) + (I146a63*cmath.sin(beta)**2)/cmath.sqrt(2) - (ee**2*I144a63*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1933 = Coupling(name = 'GC_1933',
                   value = '(I146a63*cmath.cos(beta)**2)/cmath.sqrt(2) - (ee**2*I144a63*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) - (I145a63*cmath.sin(beta)**2)/cmath.sqrt(2) + (ee**2*I144a63*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1934 = Coupling(name = 'GC_1934',
                   value = '-((I145a66*cmath.cos(beta)**2)/cmath.sqrt(2)) + (ee**2*I144a66*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) + (I146a66*cmath.sin(beta)**2)/cmath.sqrt(2) - (ee**2*I144a66*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1935 = Coupling(name = 'GC_1935',
                   value = '(I146a66*cmath.cos(beta)**2)/cmath.sqrt(2) - (ee**2*I144a66*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) - (I145a66*cmath.sin(beta)**2)/cmath.sqrt(2) + (ee**2*I144a66*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1936 = Coupling(name = 'GC_1936',
                   value = '(ee**2*I94a11*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) - (ee**2*I94a11*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1937 = Coupling(name = 'GC_1937',
                   value = '-(ee**2*I94a11*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) + (ee**2*I94a11*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1938 = Coupling(name = 'GC_1938',
                   value = '(ee**2*I94a22*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) - (ee**2*I94a22*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1939 = Coupling(name = 'GC_1939',
                   value = '-(ee**2*I94a22*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) + (ee**2*I94a22*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1940 = Coupling(name = 'GC_1940',
                   value = '(ee**2*I94a33*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) + (I95a33*cmath.sin(beta)**2)/cmath.sqrt(2) - (ee**2*I94a33*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1941 = Coupling(name = 'GC_1941',
                   value = '(I95a33*cmath.cos(beta)**2)/cmath.sqrt(2) - (ee**2*I94a33*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) + (ee**2*I94a33*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1942 = Coupling(name = 'GC_1942',
                   value = '(ee**2*I94a36*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) + (I95a36*cmath.sin(beta)**2)/cmath.sqrt(2) - (ee**2*I94a36*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1943 = Coupling(name = 'GC_1943',
                   value = '(I95a36*cmath.cos(beta)**2)/cmath.sqrt(2) - (ee**2*I94a36*cmath.cos(beta)**2)/(2.*sw**2*cmath.sqrt(2)) + (ee**2*I94a36*cmath.sin(beta)**2)/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':2})

GC_1944 = Coupling(name = 'GC_1944',
                   value = '(ee*cmath.cos(beta)**2)/(2.*sw) + (ee*cmath.sin(beta)**2)/(2.*sw)',
                   order = {'QED':1})

GC_1945 = Coupling(name = 'GC_1945',
                   value = '-(ee**2*cmath.cos(beta)**2)/(2.*sw) - (ee**2*cmath.sin(beta)**2)/(2.*sw)',
                   order = {'QED':2})

GC_1946 = Coupling(name = 'GC_1946',
                   value = '(ee**2*cmath.cos(beta)**2)/(2.*sw) + (ee**2*cmath.sin(beta)**2)/(2.*sw)',
                   order = {'QED':2})

GC_1947 = Coupling(name = 'GC_1947',
                   value = '-(cw*ee**2*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee**2*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1948 = Coupling(name = 'GC_1948',
                   value = '(cw*ee**2*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee**2*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1949 = Coupling(name = 'GC_1949',
                   value = '(ee**2*complex(0,1)*I149a44*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a44*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1950 = Coupling(name = 'GC_1950',
                   value = '-(ee**2*complex(0,1)*I149a44*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a44*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1951 = Coupling(name = 'GC_1951',
                   value = '(ee**2*complex(0,1)*I149a55*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a55*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1952 = Coupling(name = 'GC_1952',
                   value = '-(ee**2*complex(0,1)*I149a55*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a55*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1953 = Coupling(name = 'GC_1953',
                   value = '(ee**2*complex(0,1)*I16a44*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a44*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1954 = Coupling(name = 'GC_1954',
                   value = '-(ee**2*complex(0,1)*I16a44*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a44*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1955 = Coupling(name = 'GC_1955',
                   value = '(ee**2*complex(0,1)*I16a55*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a55*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1956 = Coupling(name = 'GC_1956',
                   value = '-(ee**2*complex(0,1)*I16a55*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a55*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1957 = Coupling(name = 'GC_1957',
                   value = '(ee**2*complex(0,1)*I54a44*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I54a44*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1958 = Coupling(name = 'GC_1958',
                   value = '-(ee**2*complex(0,1)*I54a44*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a44*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1959 = Coupling(name = 'GC_1959',
                   value = '(ee**2*complex(0,1)*I54a55*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I54a55*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1960 = Coupling(name = 'GC_1960',
                   value = '-(ee**2*complex(0,1)*I54a55*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a55*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1961 = Coupling(name = 'GC_1961',
                   value = '(ee**2*complex(0,1)*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1962 = Coupling(name = 'GC_1962',
                   value = '-(ee**2*complex(0,1)*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1963 = Coupling(name = 'GC_1963',
                   value = '-(ee**2*complex(0,1)*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1964 = Coupling(name = 'GC_1964',
                   value = '(ee**2*complex(0,1)*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1965 = Coupling(name = 'GC_1965',
                   value = '(ee**2*complex(0,1)*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I96a33*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I96a33*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1966 = Coupling(name = 'GC_1966',
                   value = '-(ee**2*complex(0,1)*cmath.cos(beta)**2)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.sin(beta)**2)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1967 = Coupling(name = 'GC_1967',
                   value = '-(ee**2*complex(0,1)*I148a11*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a11*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I148a11*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a11*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1968 = Coupling(name = 'GC_1968',
                   value = '-(ee**2*complex(0,1)*I148a11*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a11*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I148a11*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a11*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1969 = Coupling(name = 'GC_1969',
                   value = '(ee**2*complex(0,1)*I148a11*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a11*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I148a11*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a11*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1970 = Coupling(name = 'GC_1970',
                   value = '(ee**2*complex(0,1)*I148a11*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a11*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I148a11*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a11*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1971 = Coupling(name = 'GC_1971',
                   value = '-(ee**2*complex(0,1)*I148a22*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a22*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I148a22*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a22*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1972 = Coupling(name = 'GC_1972',
                   value = '-(ee**2*complex(0,1)*I148a22*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a22*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I148a22*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a22*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1973 = Coupling(name = 'GC_1973',
                   value = '(ee**2*complex(0,1)*I148a22*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a22*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I148a22*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a22*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1974 = Coupling(name = 'GC_1974',
                   value = '(ee**2*complex(0,1)*I148a22*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a22*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I148a22*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a22*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1975 = Coupling(name = 'GC_1975',
                   value = '(ee**2*complex(0,1)*I148a33*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a33*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a33*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a33*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a33*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a33*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a33*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a33*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a33*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a33*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1976 = Coupling(name = 'GC_1976',
                   value = '(ee**2*complex(0,1)*I148a36*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a36*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a36*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a36*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a36*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a36*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a36*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a36*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a36*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a36*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1977 = Coupling(name = 'GC_1977',
                   value = '(ee**2*complex(0,1)*I148a63*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a63*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a63*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a63*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a63*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a63*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a63*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a63*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a63*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a63*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1978 = Coupling(name = 'GC_1978',
                   value = '(ee**2*complex(0,1)*I148a66*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a66*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a66*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a66*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a66*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a66*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a66*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a66*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a66*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a66*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1979 = Coupling(name = 'GC_1979',
                   value = '-(ee**2*complex(0,1)*I15a11*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a11*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I15a11*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a11*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1980 = Coupling(name = 'GC_1980',
                   value = '-(ee**2*complex(0,1)*I15a11*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a11*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I15a11*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a11*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1981 = Coupling(name = 'GC_1981',
                   value = '(ee**2*complex(0,1)*I15a11*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a11*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I15a11*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a11*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1982 = Coupling(name = 'GC_1982',
                   value = '(ee**2*complex(0,1)*I15a11*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a11*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I15a11*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a11*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1983 = Coupling(name = 'GC_1983',
                   value = '-(ee**2*complex(0,1)*I15a22*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a22*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I15a22*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a22*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1984 = Coupling(name = 'GC_1984',
                   value = '-(ee**2*complex(0,1)*I15a22*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a22*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I15a22*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a22*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1985 = Coupling(name = 'GC_1985',
                   value = '(ee**2*complex(0,1)*I15a22*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a22*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I15a22*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a22*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1986 = Coupling(name = 'GC_1986',
                   value = '(ee**2*complex(0,1)*I15a22*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a22*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I15a22*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a22*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1987 = Coupling(name = 'GC_1987',
                   value = '(ee**2*complex(0,1)*I15a33*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a33*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I23a33*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a33*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a33*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I23a33*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a33*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a33*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a33*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a33*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1988 = Coupling(name = 'GC_1988',
                   value = '(ee**2*complex(0,1)*I15a36*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a36*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I23a36*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a36*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a36*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I23a36*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a36*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a36*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a36*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a36*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1989 = Coupling(name = 'GC_1989',
                   value = '(ee**2*complex(0,1)*I15a63*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a63*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I23a63*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a63*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a63*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I23a63*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a63*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a63*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a63*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a63*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1990 = Coupling(name = 'GC_1990',
                   value = '(ee**2*complex(0,1)*I15a66*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a66*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I23a66*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a66*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a66*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I23a66*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a66*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a66*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a66*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a66*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1991 = Coupling(name = 'GC_1991',
                   value = '(ee**2*complex(0,1)*I53a11*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I53a11*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1992 = Coupling(name = 'GC_1992',
                   value = '-(ee**2*complex(0,1)*I53a11*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a11*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I53a11*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a11*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1993 = Coupling(name = 'GC_1993',
                   value = '-(ee**2*complex(0,1)*I53a11*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I53a11*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1994 = Coupling(name = 'GC_1994',
                   value = '(ee**2*complex(0,1)*I53a11*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a11*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I53a11*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a11*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1995 = Coupling(name = 'GC_1995',
                   value = '(ee**2*complex(0,1)*I53a22*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I53a22*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1996 = Coupling(name = 'GC_1996',
                   value = '-(ee**2*complex(0,1)*I53a22*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a22*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I53a22*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a22*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1997 = Coupling(name = 'GC_1997',
                   value = '-(ee**2*complex(0,1)*I53a22*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I53a22*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1998 = Coupling(name = 'GC_1998',
                   value = '(ee**2*complex(0,1)*I53a22*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a22*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I53a22*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a22*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_1999 = Coupling(name = 'GC_1999',
                   value = '-(ee**2*complex(0,1)*I54a33*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I62a33*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a33*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I62a33*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a33*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a33*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2000 = Coupling(name = 'GC_2000',
                   value = '(ee**2*complex(0,1)*I53a33*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I54a33*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I61a33*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a33*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a33*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I61a33*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a33*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a33*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a33*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a33*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2001 = Coupling(name = 'GC_2001',
                   value = '-(ee**2*complex(0,1)*I54a36*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I62a36*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a36*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I62a36*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a36*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a36*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2002 = Coupling(name = 'GC_2002',
                   value = '(ee**2*complex(0,1)*I53a36*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I54a36*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I61a36*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a36*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a36*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I61a36*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a36*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a36*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a36*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a36*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2003 = Coupling(name = 'GC_2003',
                   value = '-(ee**2*complex(0,1)*I54a63*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I62a63*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a63*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I62a63*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a63*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a63*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2004 = Coupling(name = 'GC_2004',
                   value = '(ee**2*complex(0,1)*I53a63*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I54a63*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I61a63*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a63*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a63*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I61a63*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a63*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a63*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a63*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a63*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2005 = Coupling(name = 'GC_2005',
                   value = '-(ee**2*complex(0,1)*I54a66*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I62a66*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a66*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I62a66*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a66*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a66*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2006 = Coupling(name = 'GC_2006',
                   value = '(ee**2*complex(0,1)*I53a66*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I54a66*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I61a66*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a66*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a66*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I61a66*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a66*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a66*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a66*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a66*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2007 = Coupling(name = 'GC_2007',
                   value = '(cw*ee*complex(0,1)*cmath.cos(beta)**2)/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*sw*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*cmath.sin(beta)**2)/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*sw*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_2008 = Coupling(name = 'GC_2008',
                   value = '-((cw*ee**2*complex(0,1)*cmath.cos(beta)**2)/((-1 + sw)*sw*(1 + sw))) + (2*cw*ee**2*complex(0,1)*sw*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (cw*ee**2*complex(0,1)*cmath.sin(beta)**2)/((-1 + sw)*sw*(1 + sw)) + (2*cw*ee**2*complex(0,1)*sw*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2009 = Coupling(name = 'GC_2009',
                   value = '(2*ee**2*complex(0,1)*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(beta)**2)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (2*ee**2*complex(0,1)*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.sin(beta)**2)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2010 = Coupling(name = 'GC_2010',
                   value = '(ee**2*complex(0,1)*I148a33*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a33*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I158a33*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a33*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I158a33*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a33*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a33*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I157a33*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a33*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I157a33*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2011 = Coupling(name = 'GC_2011',
                   value = '-(ee**2*complex(0,1)*I148a33*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a33*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a33*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I148a33*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a33*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a33*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a33*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a33*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a33*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a33*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2012 = Coupling(name = 'GC_2012',
                   value = '(ee**2*complex(0,1)*I148a36*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a36*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I158a36*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a36*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I158a36*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a36*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a36*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I157a36*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a36*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I157a36*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2013 = Coupling(name = 'GC_2013',
                   value = '-(ee**2*complex(0,1)*I148a36*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a36*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a36*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I148a36*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a36*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a36*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a36*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a36*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a36*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a36*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2014 = Coupling(name = 'GC_2014',
                   value = '(ee**2*complex(0,1)*I148a63*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a63*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I158a63*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a63*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I158a63*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a63*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a63*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I157a63*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a63*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I157a63*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2015 = Coupling(name = 'GC_2015',
                   value = '-(ee**2*complex(0,1)*I148a63*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a63*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a63*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I148a63*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a63*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a63*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a63*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a63*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a63*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a63*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2016 = Coupling(name = 'GC_2016',
                   value = '(ee**2*complex(0,1)*I148a66*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a66*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I158a66*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a66*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I158a66*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a66*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a66*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I157a66*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a66*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I157a66*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2017 = Coupling(name = 'GC_2017',
                   value = '-(ee**2*complex(0,1)*I148a66*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a66*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a66*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I148a66*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a66*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I156a66*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I157a66*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a66*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I156a66*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I157a66*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2018 = Coupling(name = 'GC_2018',
                   value = '-(ee**2*complex(0,1)*I148a33*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a33*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I157a33*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a33*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I157a33*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a33*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a33*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I158a33*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a33*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I158a33*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2019 = Coupling(name = 'GC_2019',
                   value = '-(ee**2*complex(0,1)*I148a36*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a36*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I157a36*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a36*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I157a36*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a36*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a36*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I158a36*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a36*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I158a36*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2020 = Coupling(name = 'GC_2020',
                   value = '-(ee**2*complex(0,1)*I148a63*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a63*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I157a63*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a63*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I157a63*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a63*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a63*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I158a63*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a63*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I158a63*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2021 = Coupling(name = 'GC_2021',
                   value = '-(ee**2*complex(0,1)*I148a66*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I149a66*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I157a66*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a66*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I157a66*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I148a66*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I149a66*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I158a66*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I148a66*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I158a66*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2022 = Coupling(name = 'GC_2022',
                   value = '(ee**2*complex(0,1)*I15a33*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a33*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I25a33*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a33*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I25a33*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a33*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a33*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I24a33*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a33*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I24a33*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2023 = Coupling(name = 'GC_2023',
                   value = '-(ee**2*complex(0,1)*I15a33*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a33*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a33*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I15a33*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a33*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I23a33*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a33*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a33*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I23a33*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a33*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2024 = Coupling(name = 'GC_2024',
                   value = '(ee**2*complex(0,1)*I15a36*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a36*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I25a36*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a36*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I25a36*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a36*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a36*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I24a36*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a36*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I24a36*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2025 = Coupling(name = 'GC_2025',
                   value = '-(ee**2*complex(0,1)*I15a36*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a36*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a36*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I15a36*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a36*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I23a36*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a36*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a36*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I23a36*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a36*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2026 = Coupling(name = 'GC_2026',
                   value = '(ee**2*complex(0,1)*I15a63*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a63*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I25a63*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a63*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I25a63*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a63*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a63*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I24a63*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a63*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I24a63*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2027 = Coupling(name = 'GC_2027',
                   value = '-(ee**2*complex(0,1)*I15a63*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a63*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a63*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I15a63*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a63*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I23a63*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a63*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a63*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I23a63*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a63*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2028 = Coupling(name = 'GC_2028',
                   value = '(ee**2*complex(0,1)*I15a66*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a66*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I25a66*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a66*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I25a66*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a66*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a66*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I24a66*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a66*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I24a66*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2029 = Coupling(name = 'GC_2029',
                   value = '-(ee**2*complex(0,1)*I15a66*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a66*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a66*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I15a66*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a66*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I23a66*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I24a66*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a66*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I23a66*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I24a66*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2030 = Coupling(name = 'GC_2030',
                   value = '-(ee**2*complex(0,1)*I15a33*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a33*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I24a33*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a33*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I24a33*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a33*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a33*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I25a33*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a33*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I25a33*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2031 = Coupling(name = 'GC_2031',
                   value = '-(ee**2*complex(0,1)*I15a36*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a36*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I24a36*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a36*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I24a36*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a36*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a36*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I25a36*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a36*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I25a36*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2032 = Coupling(name = 'GC_2032',
                   value = '-(ee**2*complex(0,1)*I15a63*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a63*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I24a63*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a63*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I24a63*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a63*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a63*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I25a63*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a63*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I25a63*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2033 = Coupling(name = 'GC_2033',
                   value = '-(ee**2*complex(0,1)*I15a66*cmath.cos(beta)**2)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I16a66*cmath.cos(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I24a66*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a66*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I24a66*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I15a66*cmath.sin(beta)**2)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I16a66*cmath.sin(beta)**2)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I25a66*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I15a66*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I25a66*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2034 = Coupling(name = 'GC_2034',
                   value = '(ee**2*complex(0,1)*I54a33*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a33*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I54a33*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I62a33*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a33*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I62a33*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2035 = Coupling(name = 'GC_2035',
                   value = '-(ee**2*complex(0,1)*I53a33*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a33*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a33*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I53a33*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I54a33*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I61a33*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a33*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a33*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I61a33*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a33*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2036 = Coupling(name = 'GC_2036',
                   value = '(ee**2*complex(0,1)*I54a36*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a36*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I54a36*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I62a36*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a36*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I62a36*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2037 = Coupling(name = 'GC_2037',
                   value = '-(ee**2*complex(0,1)*I53a36*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a36*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a36*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I53a36*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I54a36*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I61a36*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a36*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a36*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I61a36*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a36*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2038 = Coupling(name = 'GC_2038',
                   value = '(ee**2*complex(0,1)*I54a63*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a63*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I54a63*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I62a63*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a63*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I62a63*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2039 = Coupling(name = 'GC_2039',
                   value = '-(ee**2*complex(0,1)*I53a63*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a63*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a63*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I53a63*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I54a63*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I61a63*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a63*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a63*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I61a63*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a63*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2040 = Coupling(name = 'GC_2040',
                   value = '(ee**2*complex(0,1)*I54a66*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a66*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I54a66*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I62a66*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a66*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I62a66*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2041 = Coupling(name = 'GC_2041',
                   value = '-(ee**2*complex(0,1)*I53a66*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I54a66*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I53a66*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I53a66*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I54a66*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I61a66*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I62a66*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I53a66*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I61a66*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I62a66*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2042 = Coupling(name = 'GC_2042',
                   value = '-(ee**2*complex(0,1)*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I96a33*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I96a33*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2043 = Coupling(name = 'GC_2043',
                   value = '(ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vd*cmath.cos(beta)**2*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vd*cmath.sin(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_2044 = Coupling(name = 'GC_2044',
                   value = '(ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vd*cmath.cos(beta)**2*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*vd*cmath.cos(beta)**2*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vu*cmath.cos(beta)*cmath.sin(alp)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*vu*cmath.cos(beta)*cmath.sin(alp)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vd*cmath.sin(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_2045 = Coupling(name = 'GC_2045',
                   value = '-(ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vd*cmath.cos(beta)**2*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vu*cmath.cos(beta)*cmath.sin(alp)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*vu*cmath.cos(beta)*cmath.sin(alp)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vd*cmath.sin(alp)*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*vd*cmath.sin(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_2046 = Coupling(name = 'GC_2046',
                   value = '-(ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vd*cmath.cos(beta)**2*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vd*cmath.sin(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_2047 = Coupling(name = 'GC_2047',
                   value = '(ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.cos(beta)**2)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vd*cmath.cos(beta)**2*cmath.sin(alp))/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*vd*cmath.cos(beta)**2*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*vu*cmath.cos(beta)*cmath.sin(alp)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vd*cmath.sin(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*vd*cmath.sin(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_2048 = Coupling(name = 'GC_2048',
                   value = '-(ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vu*cmath.cos(beta)**2*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vu*cmath.sin(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_2049 = Coupling(name = 'GC_2049',
                   value = '(ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.cos(beta)**2)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vu*cmath.cos(beta)**2*cmath.sin(alp))/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*vu*cmath.cos(beta)**2*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*vd*cmath.cos(beta)*cmath.sin(alp)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vu*cmath.sin(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*vu*cmath.sin(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_2050 = Coupling(name = 'GC_2050',
                   value = '(ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vu*cmath.cos(beta)**2*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vu*cmath.sin(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_2051 = Coupling(name = 'GC_2051',
                   value = '(ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vu*cmath.cos(beta)**2*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*vu*cmath.cos(beta)**2*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vd*cmath.cos(beta)*cmath.sin(alp)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*vd*cmath.cos(beta)*cmath.sin(alp)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vu*cmath.sin(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_2052 = Coupling(name = 'GC_2052',
                   value = '-(ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vu*cmath.cos(beta)**2*cmath.sin(alp))/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*vu*cmath.cos(alp)*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vd*cmath.cos(beta)*cmath.sin(alp)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*vd*cmath.cos(beta)*cmath.sin(alp)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*vd*cmath.cos(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*vu*cmath.sin(alp)*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*vu*cmath.sin(alp)*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_2053 = Coupling(name = 'GC_2053',
                   value = '(ee**2*complex(0,1)*cmath.cos(alp)*cmath.cos(beta)**2*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(alp)**2*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(alp)**2*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(alp)*cmath.sin(alp)*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2054 = Coupling(name = 'GC_2054',
                   value = '-(ee**2*complex(0,1)*cmath.cos(alp)*cmath.cos(beta)**2*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(alp)**2*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(alp)**2*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(alp)*cmath.sin(alp)*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_2055 = Coupling(name = 'GC_2055',
                   value = '(ee**2*complex(0,1)*cmath.cos(alp)*cmath.cos(beta)**2*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(alp)*cmath.sin(alp)*cmath.sin(beta)**2)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2056 = Coupling(name = 'GC_2056',
                   value = '-(ee**2*complex(0,1)*cmath.cos(alp)*cmath.cos(beta)**2*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(alp)*cmath.cos(beta)**2*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(alp)**2*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(alp)*cmath.sin(alp)*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(alp)*cmath.sin(alp)*cmath.sin(beta)**2)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2057 = Coupling(name = 'GC_2057',
                   value = '-(ee**2*complex(0,1)*cmath.cos(alp)*cmath.cos(beta)**2*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(alp)*cmath.sin(alp)*cmath.sin(beta)**2)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2058 = Coupling(name = 'GC_2058',
                   value = '(ee**2*complex(0,1)*cmath.cos(alp)*cmath.cos(beta)**2*cmath.sin(alp))/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(alp)*cmath.cos(beta)**2*cmath.sin(alp))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(alp)**2*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(alp)*cmath.sin(alp)*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(alp)*cmath.sin(alp)*cmath.sin(beta)**2)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2059 = Coupling(name = 'GC_2059',
                   value = '-(ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(beta)**2*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.sin(alp)**2*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2060 = Coupling(name = 'GC_2060',
                   value = '(ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.cos(beta)**2)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(beta)**2*cmath.sin(alp)**2)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(beta)**2*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(alp)*cmath.cos(beta)*cmath.sin(alp)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.sin(beta)**2)/(4.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.sin(alp)**2*cmath.sin(beta)**2)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.sin(alp)**2*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2061 = Coupling(name = 'GC_2061',
                   value = '(ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(beta)**2*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.sin(alp)**2*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2062 = Coupling(name = 'GC_2062',
                   value = '(ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(beta)**2*cmath.sin(alp)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(beta)**2*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(alp)*cmath.cos(beta)*cmath.sin(alp)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(alp)*cmath.cos(beta)*cmath.sin(alp)*cmath.sin(beta))/((-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.sin(alp)**2*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2063 = Coupling(name = 'GC_2063',
                   value = '-(ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(beta)**2*cmath.sin(alp)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(alp)*cmath.cos(beta)*cmath.sin(alp)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(alp)*cmath.cos(beta)*cmath.sin(alp)*cmath.sin(beta))/((-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(alp)**2*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.sin(alp)**2*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.sin(alp)**2*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2064 = Coupling(name = 'GC_2064',
                   value = '(ee**2*vu*cmath.cos(beta)**3)/(4.*sw**2) - (ee**2*vd*cmath.cos(beta)**2*cmath.sin(beta))/(4.*sw**2) + (ee**2*vu*cmath.cos(beta)*cmath.sin(beta)**2)/(4.*sw**2) - (ee**2*vd*cmath.sin(beta)**3)/(4.*sw**2)',
                   order = {'QED':1})

GC_2065 = Coupling(name = 'GC_2065',
                   value = '-(ee**2*vu*cmath.cos(beta)**3)/(4.*sw**2) + (ee**2*vd*cmath.cos(beta)**2*cmath.sin(beta))/(4.*sw**2) - (ee**2*vu*cmath.cos(beta)*cmath.sin(beta)**2)/(4.*sw**2) + (ee**2*vd*cmath.sin(beta)**3)/(4.*sw**2)',
                   order = {'QED':1})

GC_2066 = Coupling(name = 'GC_2066',
                   value = '-(ee**2*vd*cmath.cos(beta)**3)/(4.*sw**2) - (ee**2*vu*cmath.cos(beta)**2*cmath.sin(beta))/(4.*sw**2) - (ee**2*vd*cmath.cos(beta)*cmath.sin(beta)**2)/(4.*sw**2) - (ee**2*vu*cmath.sin(beta)**3)/(4.*sw**2)',
                   order = {'QED':1})

GC_2067 = Coupling(name = 'GC_2067',
                   value = '(ee**2*vd*cmath.cos(beta)**3)/(4.*sw**2) + (ee**2*vu*cmath.cos(beta)**2*cmath.sin(beta))/(4.*sw**2) + (ee**2*vd*cmath.cos(beta)*cmath.sin(beta)**2)/(4.*sw**2) + (ee**2*vu*cmath.sin(beta)**3)/(4.*sw**2)',
                   order = {'QED':1})

GC_2068 = Coupling(name = 'GC_2068',
                   value = '(ee**2*cmath.cos(beta)**3*cmath.sin(alp))/(4.*sw**2) - (ee**2*cmath.cos(alp)*cmath.cos(beta)**2*cmath.sin(beta))/(4.*sw**2) + (ee**2*cmath.cos(beta)*cmath.sin(alp)*cmath.sin(beta)**2)/(4.*sw**2) - (ee**2*cmath.cos(alp)*cmath.sin(beta)**3)/(4.*sw**2)',
                   order = {'QED':2})

GC_2069 = Coupling(name = 'GC_2069',
                   value = '-(ee**2*cmath.cos(beta)**3*cmath.sin(alp))/(4.*sw**2) + (ee**2*cmath.cos(alp)*cmath.cos(beta)**2*cmath.sin(beta))/(4.*sw**2) - (ee**2*cmath.cos(beta)*cmath.sin(alp)*cmath.sin(beta)**2)/(4.*sw**2) + (ee**2*cmath.cos(alp)*cmath.sin(beta)**3)/(4.*sw**2)',
                   order = {'QED':2})

GC_2070 = Coupling(name = 'GC_2070',
                   value = '(ee**2*complex(0,1)*cmath.cos(beta)**3*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(beta)**3)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2071 = Coupling(name = 'GC_2071',
                   value = '-(ee**2*complex(0,1)*cmath.cos(beta)**3*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(beta)**3)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2072 = Coupling(name = 'GC_2072',
                   value = '(ee**2*complex(0,1)*cmath.cos(beta)**3*cmath.sin(beta))/((-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(beta)**3)/((-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2073 = Coupling(name = 'GC_2073',
                   value = '-((ee**2*complex(0,1)*cmath.cos(beta)**3*cmath.sin(beta))/((-1 + sw)*sw**2*(1 + sw))) + (ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(beta)**3)/((-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2074 = Coupling(name = 'GC_2074',
                   value = '(3*ee**2*complex(0,1)*cmath.cos(beta)**3*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) - (3*ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(beta)**3)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2075 = Coupling(name = 'GC_2075',
                   value = '(-3*ee**2*complex(0,1)*cmath.cos(beta)**3*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (3*ee**2*complex(0,1)*cmath.cos(beta)*cmath.sin(beta)**3)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2076 = Coupling(name = 'GC_2076',
                   value = '-(ee**2*cmath.cos(alp)*cmath.cos(beta)**3)/(4.*sw**2) - (ee**2*cmath.cos(beta)**2*cmath.sin(alp)*cmath.sin(beta))/(4.*sw**2) - (ee**2*cmath.cos(alp)*cmath.cos(beta)*cmath.sin(beta)**2)/(4.*sw**2) - (ee**2*cmath.sin(alp)*cmath.sin(beta)**3)/(4.*sw**2)',
                   order = {'QED':2})

GC_2077 = Coupling(name = 'GC_2077',
                   value = '(ee**2*cmath.cos(alp)*cmath.cos(beta)**3)/(4.*sw**2) + (ee**2*cmath.cos(beta)**2*cmath.sin(alp)*cmath.sin(beta))/(4.*sw**2) + (ee**2*cmath.cos(alp)*cmath.cos(beta)*cmath.sin(beta)**2)/(4.*sw**2) + (ee**2*cmath.sin(alp)*cmath.sin(beta)**3)/(4.*sw**2)',
                   order = {'QED':2})

GC_2078 = Coupling(name = 'GC_2078',
                   value = '-(ee**2*complex(0,1)*cmath.cos(beta)**4)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (3*ee**2*complex(0,1)*cmath.cos(beta)**2*cmath.sin(beta)**2)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.sin(beta)**4)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2079 = Coupling(name = 'GC_2079',
                   value = '-(ee**2*complex(0,1)*cmath.cos(beta)**4)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (5*ee**2*complex(0,1)*cmath.cos(beta)**2*cmath.sin(beta)**2)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.sin(beta)**4)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2080 = Coupling(name = 'GC_2080',
                   value = '(ee**2*complex(0,1)*cmath.cos(beta)**4)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(beta)**4)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(beta)**2*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(beta)**2*cmath.sin(beta)**2)/(2.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.sin(beta)**4)/(4.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.sin(beta)**4)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2081 = Coupling(name = 'GC_2081',
                   value = '(ee**2*complex(0,1)*cmath.cos(beta)**4)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(beta)**2*cmath.sin(beta)**2)/(2.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.sin(beta)**4)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2082 = Coupling(name = 'GC_2082',
                   value = '-(ee**2*complex(0,1)*cmath.cos(beta)**4)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.cos(beta)**4)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(beta)**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) + (3*ee**2*complex(0,1)*cmath.cos(beta)**2*cmath.sin(beta)**2)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.sin(beta)**4)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*cmath.sin(beta)**4)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2083 = Coupling(name = 'GC_2083',
                   value = '(ee**2*complex(0,1)*cmath.cos(beta)**4)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(beta)**2*cmath.sin(beta)**2)/((-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*cmath.sin(beta)**4)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

GC_2084 = Coupling(name = 'GC_2084',
                   value = '(3*ee**2*complex(0,1)*cmath.cos(beta)**4)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (3*ee**2*complex(0,1)*cmath.cos(beta)**2*cmath.sin(beta)**2)/(2.*(-1 + sw)*sw**2*(1 + sw)) + (3*ee**2*complex(0,1)*cmath.sin(beta)**4)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':2})

