# pdb swaxs

# [a1, a2, a3, a4, b1, b2, b3, b4, c]
const COEFS = Dict{String, Array{Float64, 1}}(
"H" => [0.489918, 0.262003, 0.196767, 0.049879, 20.6593, 7.74039, 49.5519, 2.20159, 0.001305],
"H1-" => [0.897661, 0.565616, 0.415815, 0.116973, 53.1368, 15.187, 186.576, 3.56709, 0.002389],
"He" => [0.8734, 0.6309, 0.3112, 0.178, 9.1037, 3.3568, 22.9276, 0.9821, 0.0064],
"Li" => [1.1282, 0.7508, 0.6175, 0.4653, 3.9546, 1.0524, 85.3905, 168.261, 0.0377],
"Li1+" => [0.6968, 0.7888, 0.3414, 0.1563, 4.6237, 1.9557, 0.6316, 10.0953, 0.0167],
"Be" => [1.5919, 1.1278, 0.5391, 0.7029, 43.6427, 1.8623, 103.483, 0.542, 0.0385],
"Be2+" => [6.2603, 0.8849, 0.7993, 0.1647, 0.0027, 0.8313, 2.2758, 5.1146, -6.1092],
"B" => [2.0545, 1.3326, 1.0979, 0.7068, 23.2185, 1.021, 60.3498, 0.1403, -0.1932],
"C" => [2.31, 1.02, 1.5886, 0.865, 20.8439, 10.2075, 0.5687, 51.6512, 0.2156],
"Cval" => [2.26069, 1.56165, 1.05075, 0.839259, 22.6907, 0.656665, 9.75618, 55.5949, 0.286977],
"N" => [12.2126, 3.1322, 2.0125, 1.1663, 0.0057, 9.8933, 28.9975, 0.5826, -11.529],
"O" => [3.0485, 2.2868, 1.5463, 0.867, 13.2771, 5.7011, 0.3239, 32.9089, 0.2508],
"O1-" => [4.1916, 1.63969, 1.52673, -20.307, 12.8573, 4.17236, 47.0179, -0.01404, 21.9412],
"F" => [3.5392, 2.6412, 1.517, 1.0243, 10.2825, 4.2944, 0.2615, 26.1476, 0.2776],
"F1-" => [3.6322, 3.51057, 1.26064, 0.940706, 5.27756, 14.7353, 0.442258, 47.3437, 0.653396],
"Ne" => [3.9553, 3.1125, 1.4546, 1.1251, 8.4042, 3.4262, 0.2306, 21.7184, 0.3515],
"Na" => [4.7626, 3.1736, 1.2674, 1.1128, 3.285, 8.8422, 0.3136, 129.424, 0.676],
"Na1+" => [3.2565, 3.9362, 1.3998, 1.0032, 2.6671, 6.1153, 0.2001, 14.039, 0.404],
"NA" => [3.2565, 3.9362, 1.3998, 1.0032, 2.6671, 6.1153, 0.2001, 14.039, 0.404],
"Mg" => [5.4204, 2.1735, 1.2269, 2.3073, 2.8275, 79.2611, 0.3808, 7.1937, 0.8584],
"Mg2+" => [3.4988, 3.8378, 1.3284, 0.8497, 2.1676, 4.7542, 0.185, 10.1411, 0.4853],
"MG" => [3.4988, 3.8378, 1.3284, 0.8497, 2.1676, 4.7542, 0.185, 10.1411, 0.4853],
"Al" => [6.4202, 1.9002, 1.5936, 1.9646, 3.0387, 0.7426, 31.5472, 85.0886, 1.1151],
"Al3+" => [4.17448, 3.3876, 1.20296, 0.528137, 1.93816, 4.14553, 0.228753, 8.28524, 0.706786],
"Siv" => [6.2915, 3.0353, 1.9891, 1.541, 2.4386, 32.3337, 0.6785, 81.6937, 1.1407],
"Sival" => [5.66269, 3.07164, 2.62446, 1.3932, 2.6652, 38.6634, 0.916946, 93.5458, 1.24707],
"Si4+" => [4.43918, 3.20345, 1.19453, 0.41653, 1.64167, 3.43757, 0.2149, 6.65365, 0.746297],
"P" => [6.4345, 4.1791, 1.78, 1.4908, 1.9067, 27.157, 0.526, 68.1645, 1.1149],
"S" => [6.9053, 5.2034, 1.4379, 1.5863, 1.4679, 22.2151, 0.2536, 56.172, 0.8669],
"Cl" => [11.4604, 7.1964, 6.2556, 1.6455, 0.0104, 1.1662, 18.5194, 47.7784, -9.5574],
"Cl-" => [18.2915, 7.2084, 6.5337, 2.3386, 0.0066, 1.1717, 19.5424, 60.4486, -16.378],
"CL" => [18.2915, 7.2084, 6.5337, 2.3386, 0.0066, 1.1717, 19.5424, 60.4486, -16.378],
"Ar" => [7.4845, 6.7723, 0.6539, 1.6442, 0.9072, 14.8407, 43.8983, 33.3929, 1.4445],
"K" => [8.2186, 7.4398, 1.0519, 0.8659, 12.7949, 0.7748, 213.187, 41.6841, 1.4228],
"K+" => [7.9578, 7.4917, 6.359, 1.1915, 12.6331, 0.7674, -0.002, 31.9128, -4.9978],
"Ca" => [8.6266, 7.3873, 1.5899, 1.0211, 10.4421, 0.6599, 85.7484, 178.437, 1.3751],
"Ca2+" => [15.6348, 7.9518, 8.4372, 0.8537, -0.0074, 0.6089, 10.3116, 25.9905, -14.875],
"Sc" => [9.189, 7.3679, 1.6409, 1.468, 9.0213, 0.5729, 136.108, 51.3531, 1.3329],
"Sc3+" => [13.4008, 8.0273, 1.65943, 1.57936, 0.29854, 7.9629, -0.28604, 16.0662, -6.6667],
"Ti" => [9.7595, 7.3558, 1.6991, 1.9021, 7.8508, 0.5, 35.6338, 116.105, 1.2807],
"Ti2+" => [9.11423, 7.62174, 2.2793, 0.087899, 7.5243, 0.457585, 19.5361, 61.6558, 0.897155],
"Ti3+" => [17.7344, 8.73816, 5.25691, 1.92134, 0.22061, 7.04716, -0.15762, 15.9768, -14.652],
"Ti4+" => [19.5114, 8.23473, 2.01341, 1.5208, 0.178847, 6.67018, -0.29263, 12.9464, -13.28],
"V" => [10.2971, 7.3511, 2.0703, 2.0571, 6.8657, 0.4385, 26.8938, 102.478, 1.2199],
"V2+" => [10.106, 7.3541, 2.2884, 0.0223, 6.8818, 0.4409, 20.3004, 115.122, 1.2298],
"V3+" => [9.43141, 7.7419, 2.15343, 0.016865, 6.39535, 0.383349, 15.1908, 63.969, 0.656565],
"V5+" => [15.6887, 8.14208, 2.03081, -9.576, 0.679003, 5.40135, 9.97278, 0.940464, 1.7143],
"Cr" => [10.6406, 7.3537, 3.324, 1.4922, 6.1038, 0.392, 20.2626, 98.7399, 1.1832],
"Cr2+" => [9.54034, 7.7509, 3.58274, 0.509107, 5.66078, 0.344261, 13.3075, 32.4224, 0.616898],
"Cr3+" => [9.6809, 7.81136, 2.87603, 0.113575, 5.59463, 0.334393, 12.8288, 32.8761, 0.518275],
"Mn" => [11.2819, 7.3573, 3.0193, 2.2441, 5.3409, 0.3432, 17.8674, 83.7543, 1.0896],
"Mn2+" => [10.8061, 7.362, 3.5268, 0.2184, 5.2796, 0.3435, 14.343, 41.3235, 1.0874],
"Mn3+" => [9.84521, 7.87194, 3.56531, 0.323613, 4.91797, 0.294393, 10.8171, 24.1281, 0.393974],
"Mn4+" => [9.96253, 7.97057, 2.76067, 0.054447, 4.8485, 0.283303, 10.4852, 27.573, 0.251877],
"Fe" => [11.7695, 7.3573, 3.5222, 2.3045, 4.7611, 0.3072, 15.3535, 76.8805, 1.0369],
"Fe2+" => [11.0424, 7.374, 4.1346, 0.4399, 4.6538, 0.3053, 12.0546, 31.2809, 1.0097],
"Fe3+" => [11.1764, 7.3863, 3.3948, 0.0724, 4.6147, 0.3005, 11.6729, 38.5566, 0.9707],
"Co" => [12.2841, 7.3409, 4.0034, 2.3488, 4.2791, 0.2784, 13.5359, 71.1692, 1.0118],
"Co2+" => [11.2296, 7.3883, 4.7393, 0.7108, 4.1231, 0.2726, 10.2443, 25.6466, 0.9324],
"Co3+" => [10.338, 7.88173, 4.76795, 0.725591, 3.90969, 0.238668, 8.35583, 18.3491, 0.286667],
"Ni" => [12.8376, 7.292, 4.4438, 2.38, 3.8785, 0.2565, 12.1763, 66.3421, 1.0341],
"Ni2+" => [11.4166, 7.4005, 5.3442, 0.9773, 3.6766, 0.2449, 8.873, 22.1626, 0.8614],
"Ni3+" => [10.7806, 7.75868, 5.22746, 0.847114, 3.5477, 0.22314, 7.64468, 16.9673, 0.386044],
"Cu" => [13.338, 7.1676, 5.6158, 1.6735, 3.5828, 0.247, 11.3966, 64.8126, 1.191],
"Cu1+" => [11.9475, 7.3573, 6.2455, 1.5578, 3.3669, 0.2274, 8.6625, 25.8487, 0.89],
"Cu2+" => [11.8168, 7.11181, 5.78135, 1.14523, 3.37484, 0.244078, 7.9876, 19.897, 1.14431],
"Zn" => [14.0743, 7.0318, 5.1652, 2.41, 3.2655, 0.2333, 10.3163, 58.7097, 1.3041],
"Zn2+" => [11.9719, 7.3862, 6.4668, 1.394, 2.9946, 0.2031, 7.0826, 18.0995, 0.7807],
"Ga" => [15.2354, 6.7006, 4.3591, 2.9623, 3.0669, 0.2412, 10.7805, 61.4135, 1.7189],
"Ga3+" => [12.692, 6.69883, 6.06692, 1.0066, 2.81262, 0.22789, 6.36441, 14.4122, 1.53545],
"Ge" => [16.0816, 6.3747, 3.7068, 3.683, 2.8509, 0.2516, 11.4468, 54.7625, 2.1313],
"Ge4+" => [12.9172, 6.70003, 6.06791, 0.859041, 2.53718, 0.205855, 5.47913, 11.603, 1.45572],
"As" => [16.6723, 6.0701, 3.4313, 4.2779, 2.6345, 0.2647, 12.9479, 47.7972, 2.531],
"Se" => [17.0006, 5.8196, 3.9731, 4.3543, 2.4098, 0.2726, 15.2372, 43.8163, 2.8409],
"Br" => [17.1789, 5.2358, 5.6377, 3.9851, 2.1723, 16.5796, 0.2609, 41.4328, 2.9557],
"Br1-" => [17.1718, 6.3338, 5.5754, 3.7272, 2.2059, 19.3345, 0.2871, 58.1535, 3.1776],
"Kr" => [17.3555, 6.7286, 5.5493, 3.5375, 1.9384, 16.5623, 0.2261, 39.3972, 2.825],
"Rb" => [17.1784, 9.6435, 5.1399, 1.5292, 1.7888, 17.3151, 0.2748, 164.934, 3.4873],
"Rb1+" => [17.5816, 7.6598, 5.8981, 2.7817, 1.7139, 14.7957, 0.1603, 31.2087, 2.0782],
"Sr" => [17.5663, 9.8184, 5.422, 2.6694, 1.5564, 14.0988, 0.1664, 132.376, 2.5064],
"Sr2+" => [18.0874, 8.1373, 2.5654, -34.193, 1.4907, 12.6963, 24.5651, -0.0138, 41.4025],
"Y" => [17.776, 10.2946, 5.72629, 3.26588, 1.4029, 12.8006, 0.125599, 104.354, 1.91213],
"Y3+" => [17.9268, 9.1531, 1.76795, -33.108, 1.35417, 11.2145, 22.6599, -0.01319, 40.2602],
"Zr" => [17.8765, 10.948, 5.41732, 3.65721, 1.27618, 11.916, 0.117622, 87.6627, 2.06929],
"Zr4+" => [18.1668, 10.0562, 1.01118, -2.6479, 1.2148, 10.1483, 21.6054, -0.10276, 9.41454],
"Nb" => [17.6142, 12.0144, 4.04183, 3.53346, 1.18865, 11.766, 0.204785, 69.7957, 3.75591],
"Nb3+" => [19.8812, 18.0653, 11.0177, 1.94715, 0.019175, 1.13305, 10.1621, 28.3389, -12.912],
"Nb5+" => [17.9163, 13.3417, 10.799, 0.337905, 1.12446, 0.028781, 9.28206, 25.7228, -6.3934],
"Mo" => [3.7025, 17.2356, 12.8876, 3.7429, 0.2772, 1.0958, 11.004, 61.6584, 4.3875],
"Mo3+" => [21.1664, 18.2017, 11.7423, 2.30951, 0.014734, 1.03031, 9.53659, 26.6307, -14.421],
"Mo5+" => [21.0149, 18.0992, 11.4632, 0.740625, 0.014345, 1.02238, 8.78809, 23.3452, -14.316],
"Mo6+" => [17.8871, 11.175, 6.57891, 0, 1.03649, 8.48061, 0.058881, 0, 0.344941],
"Tc" => [19.1301, 11.0948, 4.64901, 2.71263, 0.864132, 8.14487, 21.5707, 86.8472, 5.40428],
"Ru" => [19.2674, 12.9182, 4.86337, 1.56756, 0.80852, 8.43467, 24.7997, 94.2928, 5.37874],
"Ru3+" => [18.5638, 13.2885, 9.32602, 3.00964, 0.847329, 8.37164, 0.017662, 22.887, -3.1892],
"Ru4+" => [18.5003, 13.1787, 4.71304, 2.18535, 0.844582, 8.12534, 0.36495, 20.8504, 1.42357],
"Rh" => [19.2957, 14.3501, 4.73425, 1.28918, 0.751536, 8.21758, 25.8749, 98.6062, 5.328],
"Rh3+" => [18.8785, 14.1259, 3.32515, -6.1989, 0.764252, 7.84438, 21.2487, -0.01036, 11.8678],
"Rh4+" => [18.8545, 13.9806, 2.53464, -5.6526, 0.760825, 7.62436, 19.3317, -0.0102, 11.2835],
"Pd" => [19.3319, 15.5017, 5.29537, 0.605844, 0.698655, 7.98929, 25.2052, 76.8986, 5.26593],
"Pd2+" => [19.1701, 15.2096, 4.32234, 0, 0.696219, 7.55573, 22.5057, 0, 5.2916],
"Pd4+" => [19.2493, 14.79, 2.89289, -7.9492, 0.683839, 7.14833, 17.9144, 0.005127, 13.0174],
"Ag" => [19.2808, 16.6885, 4.8045, 1.0463, 0.6446, 7.4726, 24.6605, 99.8156, 5.179],
"Ag1+" => [19.1812, 15.9719, 5.27475, 0.357534, 0.646179, 7.19123, 21.7326, 66.1147, 5.21572],
"Ag2+" => [19.1643, 16.2456, 4.3709, 0, 0.645643, 7.18544, 21.4072, 0, 5.21404],
"Cd" => [19.2214, 17.6444, 4.461, 1.6029, 0.5946, 6.9089, 24.7008, 87.4825, 5.0694],
"Cd2+" => [19.1514, 17.2535, 4.47128, 0, 0.597922, 6.80639, 20.2521, 0, 5.11937],
"In" => [19.1624, 18.5596, 4.2948, 2.0396, 0.5476, 6.3776, 25.8499, 92.8029, 4.9391],
"In3+" => [19.1045, 18.1108, 3.78897, 0, 0.551522, 6.3247, 17.3595, 0, 4.99635],
"Sn" => [19.1889, 19.1005, 4.4585, 2.4663, 5.8303, 0.5031, 26.8909, 83.9571, 4.7821],
"Sn2+" => [19.1094, 19.0548, 4.5648, 0.487, 0.5036, 5.8378, 23.3752, 62.2061, 4.7861],
"Sn4+" => [18.9333, 19.7131, 3.4182, 0.0193, 5.764, 0.4655, 14.0049, -0.7583, 3.9182],
"Sb" => [19.6418, 19.0455, 5.0371, 2.6827, 5.3034, 0.4607, 27.9074, 75.2825, 4.5909],
"Sb3+" => [18.9755, 18.933, 5.10789, 0.288753, 0.467196, 5.22126, 19.5902, 55.5113, 4.69626],
"Sb5+" => [19.8685, 19.0302, 2.41253, 0, 5.44853, 0.467973, 14.1259, 0, 4.69263],
"Te" => [19.9644, 19.0138, 6.14487, 2.5239, 4.81742, 0.420885, 28.5284, 70.8403, 4.352],
"I" => [20.1472, 18.9949, 7.5138, 2.2735, 4.347, 0.3814, 27.766, 66.8776, 4.0712],
"I1-" => [20.2332, 18.997, 7.8069, 2.8868, 4.3579, 0.3815, 29.5259, 84.9304, 4.0714],
"Xe" => [20.2933, 19.0298, 8.9767, 1.99, 3.9282, 0.344, 26.4659, 64.2658, 3.7118],
"Cs" => [20.3892, 19.1062, 10.662, 1.4953, 3.569, 0.3107, 24.3879, 213.904, 3.3352],
"Cs1+" => [20.3524, 19.1278, 10.2821, 0.9615, 3.552, 0.3086, 23.7128, 59.4565, 3.2791],
"Ba" => [20.3361, 19.297, 10.888, 2.6959, 3.216, 0.2756, 20.2073, 167.202, 2.7731],
"Ba2+" => [20.1807, 19.1136, 10.9054, 0.77634, 3.21367, 0.28331, 20.0558, 51.746, 3.02902],
"La" => [20.578, 19.599, 11.3727, 3.28719, 2.94817, 0.244475, 18.7726, 133.124, 2.14678],
"La3+" => [20.2489, 19.3763, 11.6323, 0.336048, 2.9207, 0.250698, 17.8211, 54.9453, 2.4086],
"Ce" => [21.1671, 19.7695, 11.8513, 3.33049, 2.81219, 0.226836, 17.6083, 127.113, 1.86264],
"Ce3+" => [20.8036, 19.559, 11.9369, 0.612376, 2.77691, 0.23154, 16.5408, 43.1692, 2.09013],
"Ce4+" => [20.3235, 19.8186, 12.1233, 0.144583, 2.65941, 0.21885, 15.7992, 62.2355, 1.5918],
"Pr" => [22.044, 19.6697, 12.3856, 2.82428, 2.77393, 0.222087, 16.7669, 143.644, 2.0583],
"Pr3+" => [21.3727, 19.7491, 12.1329, 0.97518, 2.6452, 0.214299, 15.323, 36.4065, 1.77132],
"Pr4+" => [20.9413, 20.0539, 12.4668, 0.296689, 2.54467, 0.202481, 14.8137, 45.4643, 1.24285],
"Nd" => [22.6845, 19.6847, 12.774, 2.85137, 2.66248, 0.210628, 15.885, 137.903, 1.98486],
"Nd3+" => [21.961, 19.9339, 12.12, 1.51031, 2.52722, 0.199237, 14.1783, 30.8717, 1.47588],
"Pm" => [23.3405, 19.6095, 13.1235, 2.87516, 2.5627, 0.202088, 15.1009, 132.721, 2.02876],
"Pm3+" => [22.5527, 20.1108, 12.0671, 2.07492, 2.4174, 0.185769, 13.1275, 27.4491, 1.19499],
"Sm" => [24.0042, 19.4258, 13.4396, 2.89604, 2.47274, 0.196451, 14.3996, 128.007, 2.20963],
"Sm3+" => [23.1504, 20.2599, 11.9202, 2.71488, 2.31641, 0.174081, 12.1571, 24.8242, 0.954586],
"Eu" => [24.6274, 19.0886, 13.7603, 2.9227, 2.3879, 0.1942, 13.7546, 123.174, 2.5745],
"Eu2+" => [24.0063, 19.9504, 11.8034, 3.87243, 2.27783, 0.17353, 11.6096, 26.5156, 1.36389],
"Eu3+" => [23.7497, 20.3745, 11.8509, 3.26503, 2.22258, 0.16394, 11.311, 22.9966, 0.759344],
"Gd" => [25.0709, 19.0798, 13.8518, 3.54545, 2.25341, 0.181951, 12.9331, 101.398, 2.4196],
"Gd3+" => [24.3466, 20.4208, 11.8708, 3.7149, 2.13553, 0.155525, 10.5782, 21.7029, 0.645089],
"Tb" => [25.8976, 18.2185, 14.3167, 2.95354, 2.24256, 0.196143, 12.6648, 115.362, 3.58324],
"Tb3+" => [24.9559, 20.3271, 12.2471, 3.773, 2.05601, 0.149525, 10.0499, 21.2773, 0.691967],
"Dy" => [26.507, 17.6383, 14.5596, 2.96577, 2.1802, 0.202172, 12.1899, 111.874, 4.29728],
"Dy3+" => [25.5395, 20.2861, 11.9812, 4.50073, 1.9804, 0.143384, 9.34972, 19.581, 0.68969],
"Ho" => [26.9049, 17.294, 14.5583, 3.63837, 2.07051, 0.19794, 11.4407, 92.6566, 4.56796],
"Ho3+" => [26.1296, 20.0994, 11.9788, 4.93676, 1.91072, 0.139358, 8.80018, 18.5908, 0.852795],
"Er" => [27.6563, 16.4285, 14.9779, 2.98233, 2.07356, 0.223545, 11.3604, 105.703, 5.92046],
"Er3+" => [26.722, 19.7748, 12.1506, 5.17379, 1.84659, 0.13729, 8.36225, 17.8974, 1.17613],
"Tm" => [28.1819, 15.8851, 15.1542, 2.98706, 2.02859, 0.238849, 10.9975, 102.961, 6.75621],
"Tm3+" => [27.3083, 19.332, 12.3339, 5.38348, 1.78711, 0.136974, 7.96778, 17.2922, 1.63929],
"Yb" => [28.6641, 15.4345, 15.3087, 2.98963, 1.9889, 0.257119, 10.6647, 100.417, 7.56672],
"Yb2+" => [28.1209, 17.6817, 13.3335, 5.14657, 1.78503, 0.15997, 8.18304, 20.39, 3.70983],
"Yb3+" => [27.8917, 18.7614, 12.6072, 5.47647, 1.73272, 0.13879, 7.64412, 16.8153, 2.26001],
"Lu" => [28.9476, 15.2208, 15.1, 3.71601, 1.90182, 9.98519, 0.261033, 84.3298, 7.97628],
"Lu3+" => [28.4628, 18.121, 12.8429, 5.59415, 1.68216, 0.142292, 7.33727, 16.3535, 2.97573],
"Hf" => [29.144, 15.1726, 14.7586, 4.30013, 1.83262, 9.5999, 0.275116, 72.029, 8.58154],
"Hf4+" => [28.8131, 18.4601, 12.7285, 5.59927, 1.59136, 0.128903, 6.76232, 14.0366, 2.39699],
"Ta" => [29.2024, 15.2293, 14.5135, 4.76492, 1.77333, 9.37046, 0.295977, 63.3644, 9.24354],
"Ta5+" => [29.1587, 18.8407, 12.8268, 5.38695, 1.50711, 0.116741, 6.31524, 12.4244, 1.78555],
"W" => [29.0818, 15.43, 14.4327, 5.11982, 1.72029, 9.2259, 0.321703, 57.056, 9.8875],
"W6+" => [29.4936, 19.3763, 13.0544, 5.06412, 1.42755, 0.104621, 5.93667, 11.1972, 1.01074],
"Re" => [28.7621, 15.7189, 14.5564, 5.44174, 1.67191, 9.09227, 0.3505, 52.0861, 10.472],
"Os" => [28.1894, 16.155, 14.9305, 5.67589, 1.62903, 8.97948, 0.382661, 48.1647, 11.0005],
"Os4+" => [30.419, 15.2637, 14.7458, 5.06795, 1.37113, 6.84706, 0.165191, 18.003, 6.49804],
"Ir" => [27.3049, 16.7296, 15.6115, 5.83377, 1.59279, 8.86553, 0.417916, 45.0011, 11.4722],
"Ir3+" => [30.4156, 15.862, 13.6145, 5.82008, 1.34323, 7.10909, 0.204633, 20.3254, 8.27903],
"Ir4+" => [30.7058, 15.5512, 14.2326, 5.53672, 1.30923, 6.71983, 0.167252, 17.4911, 6.96824],
"Pt" => [27.0059, 17.7639, 15.7131, 5.7837, 1.51293, 8.81174, 0.424593, 38.6103, 11.6883],
"Pt2+" => [29.8429, 16.7224, 13.2153, 6.35234, 1.32927, 7.38979, 0.263297, 22.9426, 9.85329],
"Pt4+" => [30.9612, 15.9829, 13.7348, 5.92034, 1.24813, 6.60834, 0.16864, 16.9392, 7.39534],
"Au" => [16.8819, 18.5913, 25.5582, 5.86, 0.4611, 8.6216, 1.4826, 36.3956, 12.0658],
"Au1+" => [28.0109, 17.8204, 14.3359, 6.58077, 1.35321, 7.7395, 0.356752, 26.4043, 11.2299],
"Au3+" => [30.6886, 16.9029, 12.7801, 6.52354, 1.2199, 6.82872, 0.212867, 18.659, 9.0968],
"Hg" => [20.6809, 19.0417, 21.6575, 5.9676, 0.545, 8.4484, 1.5729, 38.3246, 12.6089],
"Hg1+" => [25.0853, 18.4973, 16.8883, 6.48216, 1.39507, 7.65105, 0.443378, 28.2262, 12.0205],
"Hg2+" => [29.5641, 18.06, 12.8374, 6.89912, 1.21152, 7.05639, 0.284738, 20.7482, 10.6268],
"Tl" => [27.5446, 19.1584, 15.538, 5.52593, 0.65515, 8.70751, 1.96347, 45.8149, 13.1746],
"Tl1+" => [21.3985, 20.4723, 18.7478, 6.82847, 1.4711, 0.517394, 7.43463, 28.8482, 12.5258],
"Tl3+" => [30.8695, 18.3481, 11.9328, 7.00574, 1.1008, 6.53852, 0.219074, 17.2114, 9.8027],
"Pb" => [31.0617, 13.0637, 18.442, 5.9696, 0.6902, 2.3576, 8.618, 47.2579, 13.4118],
"Pb2+" => [21.7886, 19.5682, 19.1406, 7.01107, 1.3366, 0.488383, 6.7727, 23.8132, 12.4734],
"Pb4+" => [32.1244, 18.8003, 12.0175, 6.96886, 1.00566, 6.10926, 0.147041, 14.714, 8.08428],
"Bi" => [33.3689, 12.951, 16.5877, 6.4692, 0.704, 2.9238, 8.7937, 48.0093, 13.5782],
"Bi3+" => [21.8053, 19.5026, 19.1053, 7.10295, 1.2356, 6.24149, 0.469999, 20.3185, 12.4711],
"Bi5+" => [33.5364, 25.0946, 19.2497, 6.91555, 0.91654, 0.39042, 5.71414, 12.8285, -6.7994],
"Po" => [34.6726, 15.4733, 13.1138, 7.02588, 0.700999, 3.55078, 9.55642, 47.0045, 13.677],
"At" => [35.3163, 19.0211, 9.49887, 7.42518, 0.68587, 3.97458, 11.3824, 45.4715, 13.7108],
"Rn" => [35.5631, 21.2816, 8.0037, 7.4433, 0.6631, 4.0691, 14.0422, 44.2473, 13.6905],
"Fr" => [35.9299, 23.0547, 12.1439, 2.11253, 0.646453, 4.17619, 23.1052, 150.645, 13.7247],
"Ra" => [35.763, 22.9064, 12.4739, 3.21097, 0.616341, 3.87135, 19.9887, 142.325, 13.6211],
"Ra2+" => [35.215, 21.67, 7.91342, 7.65078, 0.604909, 3.5767, 12.601, 29.8436, 13.5431],
"Ac" => [35.6597, 23.1032, 12.5977, 4.08655, 0.589092, 3.65155, 18.599, 117.02, 13.5266],
"Ac3+" => [35.1736, 22.1112, 8.19216, 7.05545, 0.579689, 3.41437, 12.9187, 25.9443, 13.4637],
"Th" => [35.5645, 23.4219, 12.7473, 4.80703, 0.563359, 3.46204, 17.8309, 99.1722, 13.4314],
"Th4+" => [35.1007, 22.4418, 9.78554, 5.29444, 0.555054, 3.24498, 13.4661, 23.9533, 13.376],
"Pa" => [35.8847, 23.2948, 14.1891, 4.17287, 0.547751, 3.41519, 16.9235, 105.251, 13.4287],
"U" => [36.0228, 23.4128, 14.9491, 4.188, 0.5293, 3.3253, 16.0927, 100.613, 13.3966],
"U3+" => [35.5747, 22.5259, 12.2165, 5.37073, 0.52048, 3.12293, 12.7148, 26.3394, 13.3092],
"U4+" => [35.3715, 22.5326, 12.0291, 4.7984, 0.516598, 3.05053, 12.5723, 23.4582, 13.2671],
"U6+" => [34.8509, 22.7584, 14.0099, 1.21457, 0.507079, 2.8903, 13.1767, 25.2017, 13.1665],
"Np" => [36.1874, 23.5964, 15.6402, 4.1855, 0.511929, 3.25396, 15.3622, 97.4908, 13.3573],
"Np3+" => [35.7074, 22.613, 12.9898, 5.43227, 0.502322, 3.03807, 12.1449, 25.4928, 13.2544],
"Np4+" => [35.5103, 22.5787, 12.7766, 4.92159, 0.498626, 2.96627, 11.9484, 22.7502, 13.2116],
"Np6+" => [35.0136, 22.7286, 14.3884, 1.75669, 0.48981, 2.81099, 12.33, 22.6581, 13.113],
"Pu" => [36.5254, 23.8083, 16.7707, 3.47947, 0.499384, 3.26371, 14.9455, 105.98, 13.3812],
"Pu3+" => [35.84, 22.7169, 13.5807, 5.66016, 0.484938, 2.96118, 11.5331, 24.3992, 13.1991],
"Pu4+" => [35.6493, 22.646, 13.3595, 5.18831, 0.481422, 2.8902, 11.316, 21.8301, 13.1555],
"Pu6+" => [35.1736, 22.7181, 14.7635, 2.28678, 0.473204, 2.73848, 11.553, 20.9303, 13.0582],
"Am" => [36.6706, 24.0992, 17.3415, 3.49331, 0.483629, 3.20647, 14.3136, 102.273, 13.3592],
"Cm" => [36.6488, 24.4096, 17.399, 4.21665, 0.465154, 3.08997, 13.4346, 88.4834, 13.2887],
"Bk" => [36.7881, 24.7736, 17.8919, 4.23284, 0.451018, 3.04619, 12.8946, 86.003, 13.2754],
"Cf" => [36.9185, 25.1995, 18.3317, 4.24391, 0.437533, 3.00775, 12.4044, 83.7881, 13.2674]);


const TAA = Tuple{Array{String,1},Array{Float64,2}};
const AFT = typeof(COEFS);


function AtomAFF(atomid::Array{String, 1}, q::AFloat; coefs::AFT=COEFS)

    uniqueid = unique(atomid);
    idx = [id in collect(keys(coefs)) for id in uniqueid];
    uniqueid = uniqueid[idx];
    s = q/(4*pi);

    "O" in uniqueid ? nothing : push!(uniqueid, "O");
    "H" in uniqueid ? nothing : push!(uniqueid, "H");

    asf = Dict{String, Float64}();

    for id in uniqueid
        a, b, c = coefs[id][1:4], coefs[id][5:8], coefs[id][9];
        fq = c + a' * exp.(-b * s^2);
        merge!(asf, Dict{String, Float64}(id => fq));
    end

    # For solvent
    fo = asf["O"] * (1. + 0.12 * exp(-0.5 * (q / 2.2) ^ 2));
    fh = asf["H"] * (1. - 0.48 * exp(-0.5 * (q / 2.2) ^ 2));
    merge!(asf, Dict{String, Float64}("SOL-O" => fo, "SOL-H" => fh));

    # convert all the atoms to the ASFs
    n = length(atomid);
    atomaff = zeros(n);
    [atomaff[i] = asf[atomid[i]] for i = 1:n];
    return atomaff;
end







function SimplyPDB(filename::AS; coefs::AFT=COEFS, waters::Bool=true, ions::Bool=true)::TAA

    atomnames = collect(keys(coefs));
    IONS = atomnames[endswith.(atomnames, "+") .| endswith.(atomnames, "-")];
    append!(IONS, ["NA"; "K"; "CL"; "MG"; "MGH"]);
    SOL = ["SOL"; "WAT"; "HOH"];

    # Read pdb file
    f = open(filename, "r");
    lines = readlines(f);
    close(f);
    atoms = lines[startswith.(lines, "ATOM")];
    n = length(atoms);

    # trying to locate the coordinates
    sp = split(atoms[1]);
    id = 0;
    for s in sp
        id = id + 1;
        try parse(Int64, s)
            continue;
        catch er1
            if isa(er1, ArgumentError)
                try parse(Float64, s)
                    break;
                catch er2
                    isa(er2, ArgumentError) ? continue : throw(er2);
                end
            else
                throw(er1);
            end
        end
    end


    # forming the debye matrix
    mat = Array{Float64, 2}(undef, 0, 3);
    atomid = Array{String, 1}(undef, 0);
    for i = 1:n
        sp = convert.(String, split(atoms[i]));
        if (sp[4] in IONS)
            if ions
                push!(atomid, sp[end]);
                # with chain name or not
                z = parse(Float64, sp[id+2]);
                mat = z == 1.0 || z == 0.0 ? vcat(mat, parse.(Float64, sp[id-1:id+1])') : vcat(mat, parse.(Float64, sp[id:id+2])');
            end
        elseif (sp[4] in SOL)
            if waters
                push!(atomid, "SOL-" * sp[end]);
                # with chain name or not
                z = parse(Float64, sp[id+2]);
                mat = z == 1.0 || z == 0.0 ? vcat(mat, parse.(Float64, sp[id-1:id+1])') : vcat(mat, parse.(Float64, sp[id:id+2])');
            end
        else
            push!(atomid, "" * sp[end]);
            # with chain name or not
            z = parse(Float64, sp[id+2]);
            mat = z == 1.0 || z == 0.0 ? vcat(mat, parse.(Float64, sp[id-1:id+1])') : vcat(mat, parse.(Float64, sp[id:id+2])');
        end
    end
    return atomid, mat;
end




# multi dispatch of _DQ
function _DQ(sysA::TAA, sysB::TAA, q::AFloat; J::Int64=1500)

    atomidA, pdbmatA = sysA;
    atomidB, pdbmatB = sysB;
    atomaffA = AtomAFF(atomidA, q);
    atomaffB = AtomAFF(atomidB, q);

    qmat = _Orie(q, J);

    # Calculate A, B, D
    qrA = pdbmatA * qmat';
    qrB = pdbmatB * qmat';
    A = [atomaffA' * cos.(qrA); -atomaffA' * sin.(qrA)];
    B = [atomaffB' * cos.(qrB); -atomaffB' * sin.(qrB)];

    d = A .- B;
    D = mean(sum(d.^2, dims=1));

    return D;
end



# multi dispatch of _DQ
function _DQ(sysA::TAA, q::AFloat; J::Int64=1500)

    atomidA, pdbmatA = sysA;
    atomaffA = AtomAFF(atomidA, q);

    qmat = _Orie(q, J);

    # Calculate A, B, D
    qrA = pdbmatA * qmat';
    A = [atomaffA' * cos.(qrA); -atomaffA' * sin.(qrA)];

    D = mean(sum(A.^2, dims=1));

    return D;
end




function PDBSWAXS(solutefn::AS, solventfn::AS, q::AVec; J::Int64=1500, waters=true, ions=true)

    solute = SimplyPDB(solutefn; waters=waters, ions=ions);
    solvent = SimplyPDB(solventfn);

    intensity = pmap(x -> _DQ(solute, solvent, x; J=J), q, distributed=true);
    return intensity;
end





function PDBSWAXS(pdbfn::AS, q::AVec; J::Int64=1500, waters=true, ions=true)
    pdb = SimplyPDB(pdbfn; waters=waters, ions=ions);

    intensity = pmap(x -> _DQ(pdb, x; J=J), q, distributed=true);
    return intensity;
end


##new branch
