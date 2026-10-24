# The finite-difference gradient qualification (W10, §7 item 12(5))

Produced by `test/newcore/gates/fd_gradient_qualification.jl` at commit ded1543 (tree b71a0e1). Julia 1.13.0, 1 thread; load average [13.5, 13.4, 12.7] at the start and [13.8, 14.1, 13.5] at the end; 639 s. V-BO-01 at 60 PMMA cells, the temperature form, the case's solve and keywords (abstol scaled with reltol). J: the mean square misfit of the mass-loss rate against Gpyro V2 [(g m⁻² s⁻¹)²]. D(h) = (J(h) − J(−h))/(2h) ≈ dJ/d ln p. The script's header and `fd_convergence.jl` define the rule.

## The verdict

Per parameter, at reltol 1.0e-8: the converged window (three consecutive steps whose two Richardson extrapolants lie within 5 % of each other), the converged dJ/d ln p (the extrapolant of the window's smaller pair), its error (the difference of the two extrapolants) and per unit of the parameter; the reltol 1.0e-7 and 1.0e-8 differences within 5 % at each step of the window; the usable steps at each tolerance (the difference within 5 % of the converged value); the verdict; and the literal reading (the three largest steps within 5 % of one another, and its value).

| parameter | J(0) | window | dJ/d ln p | error | per unit | tolerances agree | usable at reltol 1.0e-6 | usable at reltol 1.0e-7 | usable at reltol 1.0e-8 | qualified | literal reading |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PMMA absorption coefficient α | 19.4426 | 1e-3, 2e-3, 4e-3 | -24.653 | 0.014 | -10.39 per m²/kg | pass | 1e-3, 2e-3, 4e-3 | 1e-3, 2e-3, 4e-3 | 1e-3, 2e-3, 4e-3 | pass | **fail** (the three largest disagree) |
| incident flux q | 19.4426 | 2e-3, 4e-3, 8e-3 | 957.7 | 4.7 | 0.01915 per W/m² | pass | 2e-3, 4e-3, 8e-3 | 2e-3, 4e-3, 8e-3 | 2e-3, 4e-3, 8e-3 | pass | pass, 957.7 |
| pre-exponential factor A | 19.4426 | 1e-3, 2e-3, 4e-3 | 32.567 | 0.13 | 2.412e-10 per s⁻¹ | pass | 1e-3, 2e-3 | 1e-3, 2e-3 | 1e-3, 2e-3 | pass | **fail** (the three largest disagree) |

## The differences

D(h) at each step and tolerance; the Richardson extrapolants of successive pairs; the tolerance's own window and converged value; the spread of the three largest differences ((max − min)/min |D|, the literal reading's); the runs' invocations and their accepted steps at ε = 0.

| parameter | reltol | D(1e-3) | D(2e-3) | D(4e-3) | D(8e-3) | extrapolants | window | converged | spread of the three largest | invocations | steps at ε = 0 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PMMA absorption coefficient α | 1.0e-6 | -25.168 | -24.209 | -25.415 | -28.473 | -25.487, -23.807, -24.396 | 2e-3, 4e-3, 8e-3 | -23.807 | 0.18 | 1 | 1116 |
| PMMA absorption coefficient α | 1.0e-7 | -24.694 | -24.834 | -25.419 | -28.496 | -24.647, -24.639, -24.393 | 1e-3, 2e-3, 4e-3 | -24.647 | 0.15 | 1 | 2019 |
| PMMA absorption coefficient α | 1.0e-8 | -24.698 | -24.832 | -25.411 | -28.487 | -24.653, -24.639, -24.386 | 1e-3, 2e-3, 4e-3 | -24.653 | 0.15 | 1 | 3127 |
| incident flux q | 1.0e-6 | 899.96 | 954.64 | 945.97 | 925.89 | 881.74, 957.52, 952.66 | 2e-3, 4e-3, 8e-3 | 957.52 | 0.031 | 1 | 1116 |
| incident flux q | 1.0e-7 | 899.69 | 954.84 | 946.31 | 926.18 | 881.31, 957.68, 953.02 | 2e-3, 4e-3, 8e-3 | 957.68 | 0.031 | 1 | 2019 |
| incident flux q | 1.0e-8 | 899.65 | 954.86 | 946.33 | 926.19 | 881.25, 957.7, 953.04 | 2e-3, 4e-3, 8e-3 | 957.7 | 0.031 | 1 | 3127 |
| pre-exponential factor A | 1.0e-6 | 32.925 | 32.866 | 35.359 | 41.811 | 32.945, 32.035, 33.208 | 1e-3, 2e-3, 4e-3 | 32.945 | 0.27 | 1 | 1116 |
| pre-exponential factor A | 1.0e-7 | 32.697 | 33.054 | 35.007 | 41.488 | 32.578, 32.403, 32.846 | 1e-3, 2e-3, 4e-3 | 32.578 | 0.26 | 1 | 2019 |
| pre-exponential factor A | 1.0e-8 | 32.695 | 33.079 | 34.999 | 41.487 | 32.567, 32.439, 32.837 | 1e-3, 2e-3, 4e-3 | 32.567 | 0.25 | 1 | 3127 |

## The runs

| parameter | reltol | ε | retcode | J | invocations (first at) | steps (rejected) | wall [s] |
|---|---|---|---|---|---|---|---|
| PMMA absorption coefficient α | 1.0e-6 | 0.0 | Success | 19.43916117 | 1 (357.7987) | 1116 (228) | 3.4 |
| PMMA absorption coefficient α | 1.0e-6 | 0.001 | Success | 19.41476822 | 1 (357.7983) | 1104 (203) | 3.54 |
| PMMA absorption coefficient α | 1.0e-6 | -0.001 | Success | 19.4651034 | 1 (357.7973) | 1137 (248) | 3.71 |
| PMMA absorption coefficient α | 1.0e-6 | 0.002 | Success | 19.39181751 | 1 (357.7965) | 1141 (251) | 3.89 |
| PMMA absorption coefficient α | 1.0e-6 | -0.002 | Success | 19.48865321 | 1 (357.7983) | 1107 (230) | 4.05 |
| PMMA absorption coefficient α | 1.0e-6 | 0.004 | Success | 19.33746244 | 1 (357.7967) | 1112 (233) | 3.67 |
| PMMA absorption coefficient α | 1.0e-6 | -0.004 | Success | 19.54078641 | 1 (357.7965) | 1146 (279) | 4.05 |
| PMMA absorption coefficient α | 1.0e-6 | 0.008 | Success | 19.19132216 | 1 (357.7954) | 1152 (300) | 4.17 |
| PMMA absorption coefficient α | 1.0e-6 | -0.008 | Success | 19.64689607 | 1 (357.7959) | 1119 (235) | 4.22 |
| PMMA absorption coefficient α | 1.0e-7 | 0.0 | Success | 19.44250571 | 1 (357.7958) | 2019 (820) | 6.47 |
| PMMA absorption coefficient α | 1.0e-7 | 0.001 | Success | 19.41774408 | 1 (357.7958) | 1987 (753) | 6.35 |
| PMMA absorption coefficient α | 1.0e-7 | -0.001 | Success | 19.4671324 | 1 (357.7957) | 1991 (765) | 6.77 |
| PMMA absorption coefficient α | 1.0e-7 | 0.002 | Success | 19.39271772 | 1 (357.7957) | 2003 (773) | 5.42 |
| PMMA absorption coefficient α | 1.0e-7 | -0.002 | Success | 19.49205462 | 1 (357.7955) | 2020 (806) | 5.92 |
| PMMA absorption coefficient α | 1.0e-7 | 0.004 | Success | 19.33949899 | 1 (357.7953) | 1976 (718) | 8.79 |
| PMMA absorption coefficient α | 1.0e-7 | -0.004 | Success | 19.54284826 | 1 (357.7948) | 1998 (780) | 6.26 |
| PMMA absorption coefficient α | 1.0e-7 | 0.008 | Success | 19.19610185 | 1 (357.7932) | 2021 (830) | 7.83 |
| PMMA absorption coefficient α | 1.0e-7 | -0.008 | Success | 19.65204148 | 1 (357.7923) | 2009 (794) | 7.47 |
| PMMA absorption coefficient α | 1.0e-8 | 0.0 | Success | 19.442583 | 1 (357.7957) | 3127 (1167) | 11.5 |
| PMMA absorption coefficient α | 1.0e-8 | 0.001 | Success | 19.41786123 | 1 (357.7957) | 3098 (1117) | 12.3 |
| PMMA absorption coefficient α | 1.0e-8 | -0.001 | Success | 19.46725633 | 1 (357.7956) | 3117 (1062) | 7.37 |
| PMMA absorption coefficient α | 1.0e-8 | 0.002 | Success | 19.39278944 | 1 (357.7956) | 3117 (1119) | 9.35 |
| PMMA absorption coefficient α | 1.0e-8 | -0.002 | Success | 19.49211683 | 1 (357.7954) | 3122 (1155) | 10.4 |
| PMMA absorption coefficient α | 1.0e-8 | 0.004 | Success | 19.33966719 | 1 (357.7952) | 3118 (1103) | 10.6 |
| PMMA absorption coefficient α | 1.0e-8 | -0.004 | Success | 19.54295688 | 1 (357.7947) | 3104 (1215) | 10.0 |
| PMMA absorption coefficient α | 1.0e-8 | 0.008 | Success | 19.19629932 | 1 (357.7931) | 3118 (1191) | 10.4 |
| PMMA absorption coefficient α | 1.0e-8 | -0.008 | Success | 19.65209804 | 1 (357.7922) | 3094 (1072) | 13.7 |
| incident flux q | 1.0e-6 | 0.0 | Success | 19.43916117 | 1 (357.7987) | 1116 (228) | 4.9 |
| incident flux q | 1.0e-6 | 0.001 | Success | 20.18640534 | 1 (357.3412) | 1136 (244) | 3.99 |
| incident flux q | 1.0e-6 | -0.001 | Success | 18.38648147 | 1 (358.2571) | 1126 (214) | 3.72 |
| incident flux q | 1.0e-6 | 0.002 | Success | 21.35955479 | 1 (356.8847) | 1137 (245) | 4.34 |
| incident flux q | 1.0e-6 | -0.002 | Success | 17.54100965 | 1 (358.7161) | 1112 (213) | 3.85 |
| incident flux q | 1.0e-6 | 0.004 | Success | 23.24762053 | 1 (355.9751) | 1139 (263) | 4.07 |
| incident flux q | 1.0e-6 | -0.004 | Success | 15.67985109 | 1 (359.6381) | 1152 (289) | 3.65 |
| incident flux q | 1.0e-6 | 0.008 | Success | 26.90019813 | 1 (354.171) | 1107 (211) | 3.66 |
| incident flux q | 1.0e-6 | -0.008 | Success | 12.08590323 | 1 (361.4984) | 1125 (263) | 3.63 |
| incident flux q | 1.0e-7 | 0.0 | Success | 19.44250571 | 1 (357.7958) | 2019 (820) | 7.19 |
| incident flux q | 1.0e-7 | 0.001 | Success | 20.19351535 | 1 (357.3384) | 1983 (766) | 7.7 |
| incident flux q | 1.0e-7 | -0.001 | Success | 18.39413536 | 1 (358.2543) | 2009 (808) | 9.37 |
| incident flux q | 1.0e-7 | 0.002 | Success | 21.36318625 | 1 (356.8822) | 2002 (775) | 6.99 |
| incident flux q | 1.0e-7 | -0.002 | Success | 17.54383138 | 1 (358.7139) | 2006 (791) | 8.73 |
| incident flux q | 1.0e-7 | 0.004 | Success | 23.25243912 | 1 (355.9731) | 1984 (733) | 9.17 |
| incident flux q | 1.0e-7 | -0.004 | Success | 15.68195331 | 1 (359.6365) | 1994 (771) | 5.1 |
| incident flux q | 1.0e-7 | 0.008 | Success | 26.90940451 | 1 (354.1679) | 1992 (715) | 5.04 |
| incident flux q | 1.0e-7 | -0.008 | Success | 12.09056458 | 1 (361.4954) | 2023 (878) | 5.23 |
| incident flux q | 1.0e-8 | 0.0 | Success | 19.442583 | 1 (357.7957) | 3127 (1167) | 7.97 |
| incident flux q | 1.0e-8 | 0.001 | Success | 20.1936887 | 1 (357.3384) | 3101 (1124) | 7.4 |
| incident flux q | 1.0e-8 | -0.001 | Success | 18.39439022 | 1 (358.2542) | 3112 (1146) | 8.72 |
| incident flux q | 1.0e-8 | 0.002 | Success | 21.36336428 | 1 (356.8821) | 3111 (1184) | 11.1 |
| incident flux q | 1.0e-8 | -0.002 | Success | 17.54392036 | 1 (358.7138) | 3075 (1052) | 9.94 |
| incident flux q | 1.0e-8 | 0.004 | Success | 23.25270468 | 1 (355.9729) | 3129 (1120) | 11.2 |
| incident flux q | 1.0e-8 | -0.004 | Success | 15.68206387 | 1 (359.6364) | 3118 (1134) | 11.0 |
| incident flux q | 1.0e-8 | 0.008 | Success | 26.90970722 | 1 (354.1678) | 3130 (1152) | 12.7 |
| incident flux q | 1.0e-8 | -0.008 | Success | 12.0906312 | 1 (361.4954) | 3099 (1098) | 7.55 |
| pre-exponential factor A | 1.0e-6 | 0.0 | Success | 19.43916117 | 1 (357.7987) | 1116 (228) | 2.7 |
| pre-exponential factor A | 1.0e-6 | 0.001 | Success | 19.47212961 | 1 (357.7704) | 1098 (213) | 3.56 |
| pre-exponential factor A | 1.0e-6 | -0.001 | Success | 19.40627898 | 1 (357.8267) | 1121 (229) | 4.91 |
| pre-exponential factor A | 1.0e-6 | 0.002 | Success | 19.50374867 | 1 (357.7434) | 1097 (204) | 2.68 |
| pre-exponential factor A | 1.0e-6 | -0.002 | Success | 19.37228395 | 1 (357.8545) | 1097 (188) | 4.05 |
| pre-exponential factor A | 1.0e-6 | 0.004 | Success | 19.57365718 | 1 (357.6863) | 1115 (211) | 3.22 |
| pre-exponential factor A | 1.0e-6 | -0.004 | Success | 19.29078717 | 1 (357.9111) | 1096 (215) | 2.85 |
| pre-exponential factor A | 1.0e-6 | 0.008 | Success | 19.72531019 | 1 (357.5756) | 1113 (211) | 3.45 |
| pre-exponential factor A | 1.0e-6 | -0.008 | Success | 19.05632647 | 1 (358.0239) | 1111 (232) | 3.52 |
| pre-exponential factor A | 1.0e-7 | 0.0 | Success | 19.44250571 | 1 (357.7958) | 2019 (820) | 7.61 |
| pre-exponential factor A | 1.0e-7 | 0.001 | Success | 19.47507602 | 1 (357.7678) | 2020 (826) | 7.07 |
| pre-exponential factor A | 1.0e-7 | -0.001 | Success | 19.40968199 | 1 (357.8238) | 2012 (828) | 6.62 |
| pre-exponential factor A | 1.0e-7 | 0.002 | Success | 19.50794611 | 1 (357.7399) | 1999 (779) | 7.15 |
| pre-exponential factor A | 1.0e-7 | -0.002 | Success | 19.37573161 | 1 (357.8518) | 2003 (785) | 5.04 |
| pre-exponential factor A | 1.0e-7 | 0.004 | Success | 19.5763196 | 1 (357.6841) | 2045 (833) | 6.18 |
| pre-exponential factor A | 1.0e-7 | -0.004 | Success | 19.29626743 | 1 (357.908) | 1992 (768) | 8.24 |
| pre-exponential factor A | 1.0e-7 | 0.008 | Success | 19.72919791 | 1 (357.573) | 1994 (774) | 6.0 |
| pre-exponential factor A | 1.0e-7 | -0.008 | Success | 19.06538922 | 1 (358.0206) | 2001 (799) | 6.36 |
| pre-exponential factor A | 1.0e-8 | 0.0 | Success | 19.442583 | 1 (357.7957) | 3127 (1167) | 13.3 |
| pre-exponential factor A | 1.0e-8 | 0.001 | Success | 19.4751578 | 1 (357.7678) | 3096 (1098) | 10.4 |
| pre-exponential factor A | 1.0e-8 | -0.001 | Success | 19.40976818 | 1 (357.8237) | 3125 (1128) | 8.95 |
| pre-exponential factor A | 1.0e-8 | 0.002 | Success | 19.50806918 | 1 (357.7398) | 3097 (1082) | 7.79 |
| pre-exponential factor A | 1.0e-8 | -0.002 | Success | 19.37575301 | 1 (357.8518) | 3103 (1065) | 8.63 |
| pre-exponential factor A | 1.0e-8 | 0.004 | Success | 19.57636275 | 1 (357.6841) | 3122 (1163) | 14.1 |
| pre-exponential factor A | 1.0e-8 | -0.004 | Success | 19.29636948 | 1 (357.9079) | 3107 (1098) | 12.3 |
| pre-exponential factor A | 1.0e-8 | 0.008 | Success | 19.72929705 | 1 (357.5729) | 3119 (1189) | 10.5 |
| pre-exponential factor A | 1.0e-8 | -0.008 | Success | 19.06550859 | 1 (358.0206) | 3120 (1142) | 9.87 |
