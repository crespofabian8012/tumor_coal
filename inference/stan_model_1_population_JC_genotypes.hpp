
//#include <stan⁩/lib⁩/stan_math⁩/stan/math/prim/fun/Eigen.hpp>
//#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/src/stan/model/model_header.hpp>

namespace stan_model_1_population_JC_genotypes_model_namespace {

using stan::io::dump;
using stan::model::assign;
using stan::model::index_uni;
using stan::model::index_max;
using stan::model::index_min;
using stan::model::index_min_max;
using stan::model::index_multi;
using stan::model::index_omni;
using stan::model::model_base_crtp;
using stan::model::rvalue;
using namespace stan::math;


stan::math::profile_map profiles__;
static constexpr std::array<const char*, 867> locations_array__ =
{" (found before start of program)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1106, column 2 to column 47)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1108, column 2 to column 25)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1109, column 2 to column 59)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1110, column 2 to column 24)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1112, column 2 to column 38)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1117, column 2 to column 86)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1118, column 2 to column 142)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1121, column 2 to column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1122, column 2 to column 21)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1123, column 2 to column 18)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1124, column 2 to column 19)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1125, column 22 to column 39)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1125, column 40 to column 41)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1125, column 2 to column 43)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1126, column 27 to column 42)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1126, column 2 to column 44)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1127, column 9 to column 24)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1127, column 2 to column 41)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1128, column 9 to column 28)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1128, column 2 to column 37)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1129, column 9 to column 24)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1129, column 2 to column 39)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1131, column 2 to column 52)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1133, column 2 to column 52)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1135, column 2 to column 25)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1137, column 2 to column 95)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1140, column 2 to column 85)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1142, column 2 to column 29)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1145, column 2 to line 1147, column 52)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1149, column 2 to column 140)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1152, column 6 to column 63)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1153, column 6 to column 43)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1155, column 8 to column 37)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1154, column 20 to line 1156, column 7)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1154, column 6 to line 1156, column 7)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1151, column 29 to line 1157, column 5)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1151, column 1 to line 1157, column 5)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1163, column 10 to column 141)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1164, column 10 to column 141)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1167, column 10 to column 108)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1162, column 38 to line 1168, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1162, column 2 to line 1168, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1171, column 9 to column 140)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1169, column 22 to line 1172, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1169, column 8 to line 1172, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1173, column 8 to column 58)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1161, column 17 to line 1174, column 4)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1161, column 2 to line 1174, column 4)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1032, column 2 to column 31)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1033, column 2 to column 17)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1034, column 2 to column 33)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1035, column 2 to column 18)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1036, column 41 to column 58)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1036, column 59 to column 60)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1036, column 2 to column 62)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1037, column 50 to column 69)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1037, column 2 to column 73)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1038, column 55 to column 72)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1038, column 2 to column 74)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1042, column 4 to column 48)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1043, column 29 to column 50)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1043, column 4 to column 52)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1045, column 41 to column 62)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1045, column 4 to column 65)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1046, column 37 to column 56)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1046, column 57 to column 58)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1046, column 4 to column 60)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1047, column 4 to column 18)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1048, column 11 to column 30)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1048, column 32 to column 53)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1048, column 4 to column 134)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1049, column 39 to column 60)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1049, column 4 to column 98)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1050, column 11 to column 26)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1050, column 28 to column 47)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1050, column 4 to column 132)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1051, column 12 to column 33)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1051, column 4 to column 37)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1052, column 10 to column 31)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1052, column 4 to column 33)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1053, column 10 to column 27)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1053, column 4 to column 29)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1056, column 7 to column 33)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1057, column 7 to column 45)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1055, column 35 to line 1058, column 5)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1055, column 4 to line 1058, column 5)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1060, column 7 to column 116)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1061, column 7 to column 154)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1062, column 7 to column 55)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1063, column 8 to column 51)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1064, column 8 to column 49)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1059, column 40 to line 1065, column 4)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1059, column 7 to line 1065, column 4)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1069, column 11 to column 101)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1068, column 44 to line 1070, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1068, column 8 to line 1070, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1073, column 11 to column 99)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1071, column 44 to line 1074, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1071, column 8 to line 1074, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1067, column 3 to line 1075, column 5)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1066, column 1 to line 1075, column 5)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1077, column 4 to column 50)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1078, column 4 to column 50)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1079, column 4 to column 50)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1084, column 17 to column 65)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1085, column 17 to column 48)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1086, column 17 to column 82)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1082, column 26 to line 1088, column 11)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1082, column 10 to line 1088, column 11)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1081, column 36 to line 1089, column 7)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1081, column 4 to line 1089, column 7)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1093, column 12 to column 60)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1095, column 14 to column 107)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1094, column 50 to line 1096, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1094, column 12 to line 1096, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1098, column 17 to column 136)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1097, column 49 to line 1100, column 12)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1097, column 11 to line 1100, column 12)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1101, column 11 to column 50)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1091, column 40 to line 1102, column 7)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1091, column 6 to line 1102, column 7)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1112, column 10 to column 27)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1117, column 19 to column 38)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 3, column 5 to column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 4, column 5 to column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 7, column 10 to column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 6, column 5 to line 8, column 6)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 5, column 5 to line 8, column 6)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 12, column 9 to column 18)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 10, column 8 to column 18)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 9, column 5 to line 12, column 18)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 13, column 6 to column 21)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 2, column 39 to line 14, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 16, column 5 to column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 17, column 5 to column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 20, column 10 to column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 19, column 5 to line 21, column 6)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 18, column 5 to line 21, column 6)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 25, column 9 to column 18)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 23, column 8 to column 18)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 22, column 5 to line 25, column 18)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 26, column 6 to column 21)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 15, column 33 to line 27, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 29, column 5 to column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 30, column 5 to column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 33, column 10 to column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 32, column 5 to line 34, column 6)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 31, column 5 to line 34, column 6)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 38, column 9 to column 18)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 36, column 8 to column 18)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 35, column 5 to line 38, column 18)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 39, column 6 to column 21)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 28, column 47 to line 40, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 42, column 18 to column 19)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 42, column 7 to column 21)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 43, column 7 to column 25)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 45, column 12 to column 35)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 46, column 12 to column 34)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 44, column 21 to line 47, column 8)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 44, column 7 to line 47, column 8)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 48, column 7 to column 22)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 41, column 52 to line 49, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 51, column 7 to column 24)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 52, column 7 to column 22)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 53, column 7 to column 25)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 64, column 14 to column 35)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 65, column 14 to column 41)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 63, column 16 to line 66, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 60, column 14 to column 64)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 61, column 14 to column 41)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 59, column 72 to line 62, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 59, column 17 to line 66, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 56, column 14 to column 64)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 57, column 14 to column 35)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 55, column 67 to line 58, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 55, column 12 to line 66, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 54, column 44 to line 67, column 8)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 54, column 7 to line 67, column 8)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 68, column 7 to column 22)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 50, column 76 to line 69, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 71, column 4 to column 12)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 72, column 4 to column 10)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 74, column 12 to column 23)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 73, column 149 to line 75, column 8)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 73, column 4 to line 75, column 8)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 76, column 4 to column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 70, column 162 to line 77, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 80, column 4 to column 160)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 81, column 11 to column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 81, column 4 to column 88)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 82, column 4 to column 19)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 78, column 148 to line 83, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 86, column 4 to column 21)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 87, column 4 to column 12)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 88, column 4 to column 25)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 89, column 4 to column 10)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 95, column 14 to column 130)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 94, column 16 to line 96, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 92, column 14 to column 80)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 91, column 23 to line 93, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 91, column 12 to line 96, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 97, column 12 to column 89)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 98, column 12 to column 23)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 90, column 149 to line 99, column 8)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 90, column 4 to line 99, column 8)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 100, column 4 to column 19)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 85, column 154 to line 101, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 106, column 4 to column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 107, column 4 to column 33)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 108, column 4 to column 17)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 109, column 4 to column 18)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 110, column 4 to column 23)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 111, column 4 to column 24)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 112, column 4 to column 17)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 113, column 4 to column 22)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 115, column 12 to column 47)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 115, column 4 to column 73)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 116, column 12 to column 47)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 116, column 4 to column 75)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 121, column 27 to column 52)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 122, column 27 to column 53)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 123, column 27 to column 121)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 124, column 27 to column 123)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 132, column 32 to line 141, column 29)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 127, column 30 to column 78)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 128, column 30 to column 76)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 129, column 30 to column 132)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 130, column 30 to column 56)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 126, column 51 to line 131, column 29)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 126, column 27 to line 141, column 29)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 148, column 31 to line 157, column 29)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 143, column 28 to column 76)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 144, column 28 to column 76)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 145, column 28 to column 130)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 146, column 28 to column 54)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 142, column 51 to line 147, column 28)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 142, column 26 to line 157, column 29)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 120, column 97 to line 158, column 24)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 120, column 4 to line 158, column 24)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 165, column 4 to column 38)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 104, column 102 to line 166, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 170, column 20 to column 32)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 171, column 20 to column 49)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 172, column 20 to column 33)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 173, column 20 to column 34)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 174, column 20 to column 39)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 175, column 20 to column 40)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 176, column 20 to column 33)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 177, column 20 to column 38)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 178, column 20 to line 179, column 103)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 182, column 52 to column 87)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 182, column 20 to column 143)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 183, column 28 to column 63)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 183, column 20 to column 89)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 184, column 28 to column 63)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 184, column 20 to column 91)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 187, column 27 to column 52)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 188, column 27 to column 53)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 189, column 27 to column 121)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 190, column 27 to column 123)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 201, column 34 to column 82)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 202, column 34 to column 80)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 203, column 34 to column 80)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 204, column 34 to column 60)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 200, column 69 to line 206, column 32)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 200, column 31 to line 206, column 32)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 198, column 32 to line 207, column 29)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 193, column 30 to column 78)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 194, column 30 to column 76)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 195, column 30 to column 132)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 196, column 30 to column 56)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 192, column 51 to line 197, column 29)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 192, column 27 to line 207, column 29)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 217, column 34 to column 82)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 218, column 34 to column 81)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 219, column 34 to column 80)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 220, column 34 to column 60)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 216, column 70 to line 222, column 32)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 216, column 31 to line 222, column 32)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 214, column 31 to line 223, column 29)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 209, column 28 to column 76)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 210, column 28 to column 76)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 211, column 28 to column 130)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 212, column 28 to column 54)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 208, column 51 to line 213, column 28)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 208, column 26 to line 223, column 29)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 186, column 113 to line 224, column 24)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 186, column 20 to line 224, column 24)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 232, column 14 to column 44)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 230, column 15 to column 53)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 229, column 12 to line 232, column 44)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 235, column 12 to column 27)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 168, column 62 to line 236, column 11)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 242, column 11 to column 60)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 242, column 4 to column 87)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 243, column 4 to column 12)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 244, column 4 to column 109)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 245, column 23 to column 34)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 245, column 4 to column 36)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 246, column 4 to column 17)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 247, column 4 to column 19)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 248, column 4 to column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 249, column 4 to column 25)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 250, column 4 to column 27)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 251, column 4 to column 22)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 254, column 4 to column 59)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 255, column 4 to column 30)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 282, column 15 to column 116)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 284, column 18 to column 73)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 285, column 18 to column 44)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 283, column 66 to line 286, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 283, column 16 to line 286, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 288, column 19 to column 74)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 289, column 19 to column 45)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 287, column 65 to line 290, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 287, column 15 to line 290, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 291, column 15 to column 52)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 281, column 15 to line 292, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 280, column 9 to line 292, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 279, column 6 to line 293, column 7)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 266, column 15 to column 116)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 268, column 18 to column 73)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 269, column 18 to column 44)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 267, column 65 to line 270, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 267, column 15 to line 270, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 272, column 19 to column 74)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 273, column 19 to column 45)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 271, column 65 to line 274, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 271, column 15 to line 274, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 275, column 15 to column 52)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 265, column 15 to line 276, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 262, column 10 to line 276, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 261, column 6 to line 277, column 7)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 260, column 6 to line 293, column 7)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 296, column 6 to column 33)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 259, column 82 to line 297, column 7)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 259, column 4 to line 297, column 7)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 306, column 4 to column 79)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 241, column 79 to line 307, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 313, column 15 to column 32)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 313, column 4 to column 34)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 314, column 4 to column 17)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 315, column 4 to column 18)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 316, column 4 to column 108)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 317, column 23 to column 36)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 317, column 4 to column 38)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 318, column 23 to column 36)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 318, column 4 to column 38)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 319, column 4 to column 17)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 320, column 4 to column 19)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 321, column 4 to column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 322, column 4 to column 25)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 323, column 4 to column 22)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 324, column 36 to column 49)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 324, column 4 to column 83)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 325, column 4 to column 14)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 327, column 4 to column 59)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 328, column 4 to column 30)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 364, column 18 to column 63)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 365, column 18 to column 28)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 363, column 19 to line 366, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 360, column 18 to column 73)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 361, column 18 to column 44)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 359, column 64 to line 362, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 359, column 14 to line 366, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 372, column 18 to column 63)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 373, column 18 to column 28)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 371, column 19 to line 374, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 368, column 19 to column 74)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 369, column 19 to column 45)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 367, column 65 to line 370, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 367, column 15 to line 374, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 358, column 12 to line 375, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 357, column 9 to line 375, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 356, column 6 to line 376, column 7)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 342, column 18 to column 63)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 343, column 18 to column 28)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 341, column 19 to line 344, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 338, column 18 to column 73)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 339, column 18 to column 44)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 337, column 65 to line 340, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 337, column 15 to line 344, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 350, column 18 to column 63)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 351, column 18 to column 28)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 349, column 19 to line 352, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 346, column 19 to column 74)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 347, column 19 to column 45)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 345, column 65 to line 348, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 345, column 15 to line 352, column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 336, column 12 to line 353, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 333, column 10 to line 353, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 332, column 6 to line 354, column 7)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 331, column 6 to line 376, column 7)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 379, column 6 to column 33)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 330, column 58 to line 381, column 5)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 330, column 4 to line 381, column 5)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 382, column 4 to column 92)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 383, column 4 to column 19)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 311, column 60 to line 384, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 389, column 13 to column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 389, column 6 to column 29)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 390, column 6 to column 14)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 401, column 21 to column 61)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 400, column 32 to line 402, column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 400, column 16 to line 402, column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 403, column 15 to column 30)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 399, column 17 to line 404, column 14)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 395, column 21 to column 61)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 394, column 35 to line 396, column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 394, column 19 to line 396, column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 397, column 18 to column 72)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 393, column 53 to line 398, column 14)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 393, column 14 to line 404, column 14)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 405, column 11 to column 17)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 392, column 20 to line 406, column 10)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 392, column 9 to line 406, column 10)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 391, column 19 to line 407, column 7)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 391, column 6 to line 407, column 7)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 408, column 5 to column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 387, column 51 to line 409, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 411, column 5 to column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 412, column 5 to column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 413, column 5 to column 27)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 416, column 12 to column 47)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 415, column 62 to line 417, column 11)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 415, column 10 to line 417, column 11)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 418, column 10 to column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 414, column 87 to line 419, column 6)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 414, column 5 to line 419, column 6)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 423, column 9 to column 31)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 421, column 8 to column 18)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 420, column 5 to line 423, column 31)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 424, column 6 to column 21)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 410, column 114 to line 425, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 431, column 10 to column 30)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 432, column 10 to column 40)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 433, column 10 to column 53)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 434, column 10 to column 59)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 435, column 10 to column 17)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 436, column 10 to column 25)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 437, column 10 to column 28)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 443, column 11 to column 43)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 444, column 11 to column 22)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 445, column 11 to column 37)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 446, column 11 to column 34)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 447, column 11 to column 41)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 448, column 11 to column 51)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 449, column 11 to column 41)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 450, column 11 to column 51)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 452, column 11 to column 62)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 453, column 11 to column 30)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 442, column 13 to line 454, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 439, column 12 to column 45)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 440, column 12 to column 42)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 438, column 19 to line 441, column 11)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 438, column 10 to line 454, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 430, column 9 to line 455, column 10)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 461, column 8 to column 15)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 462, column 8 to column 15)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 463, column 8 to column 15)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 464, column 8 to column 15)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 465, column 8 to column 15)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 466, column 8 to column 15)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 467, column 8 to column 23)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 468, column 8 to column 38)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 469, column 8 to column 23)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 470, column 8 to column 23)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 471, column 8 to column 25)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 472, column 8 to column 27)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 504, column 10 to column 28)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 505, column 10 to column 42)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 506, column 10 to column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 507, column 10 to column 35)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 508, column 10 to column 33)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 509, column 10 to column 27)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 510, column 10 to column 33)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 511, column 10 to column 48)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 512, column 10 to column 49)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 503, column 12 to line 513, column 11)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 474, column 11 to column 44)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 475, column 11 to column 45)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 476, column 11 to column 21)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 477, column 11 to column 35)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 478, column 11 to column 21)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 479, column 11 to column 38)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 480, column 11 to column 100)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 481, column 11 to column 77)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 482, column 12 to column 51)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 486, column 12 to column 32)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 487, column 12 to column 28)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 488, column 12 to column 29)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 489, column 12 to column 28)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 490, column 12 to column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 491, column 12 to column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 492, column 12 to column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 493, column 12 to column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 494, column 12 to column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 495, column 12 to column 26)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 496, column 12 to column 35)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 497, column 12 to column 30)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 498, column 12 to column 72)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 499, column 12 to column 28)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 500, column 12 to column 51)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 485, column 8 to line 501, column 12)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 484, column 8 to line 501, column 12)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 473, column 17 to line 502, column 8)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 473, column 8 to line 513, column 11)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 514, column 8 to column 29)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 460, column 9 to line 515, column 10)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 520, column 16 to column 43)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 520, column 8 to column 52)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 525, column 14 to column 27)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 523, column 15 to column 100)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 522, column 11 to line 525, column 27)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 521, column 47 to line 526, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 521, column 8 to line 526, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 527, column 8 to column 23)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 518, column 16 to line 528, column 7)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 530, column 4 to column 24)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 531, column 4 to column 23)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 532, column 4 to column 25)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 534, column 6 to column 73)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 535, column 6 to column 77)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 536, column 6 to column 111)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 533, column 27 to line 537, column 5)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 533, column 4 to line 537, column 5)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 538, column 4 to column 19)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 529, column 120 to line 539, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 541, column 8 to column 63)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 542, column 8 to column 15)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 543, column 8 to column 21)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 544, column 8 to column 53)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 545, column 8 to column 18)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 546, column 8 to column 74)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 547, column 8 to column 78)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 548, column 8 to column 37)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 553, column 12 to column 46)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 552, column 12 to line 554, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 550, column 10 to column 65)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 549, column 22 to line 551, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 549, column 8 to line 554, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 555, column 8 to column 24)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 540, column 65 to line 556, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 558, column 6 to column 18)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 559, column 6 to column 17)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 560, column 6 to column 17)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 561, column 6 to column 17)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 562, column 6 to column 19)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 563, column 6 to column 32)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 564, column 6 to column 28)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 565, column 6 to column 24)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 566, column 6 to column 47)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 567, column 6 to column 49)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 569, column 6 to column 24)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 570, column 6 to column 21)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 557, column 76 to line 571, column 4)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 573, column 6 to column 18)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 574, column 6 to column 75)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 575, column 6 to column 21)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 572, column 77 to line 576, column 4)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 578, column 5 to column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 579, column 5 to column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 580, column 5 to column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 581, column 5 to column 17)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 582, column 5 to column 35)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 583, column 5 to column 27)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 584, column 5 to column 23)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 585, column 5 to column 39)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 587, column 5 to column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 577, column 75 to line 588, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 590, column 5 to column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 591, column 5 to column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 592, column 5 to column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 593, column 5 to column 17)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 594, column 5 to column 35)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 595, column 5 to column 27)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 596, column 5 to column 23)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 597, column 5 to column 33)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 598, column 4 to column 24)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 589, column 78 to line 599, column 2)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 601, column 4 to column 15)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 602, column 4 to column 15)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 603, column 4 to column 15)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 604, column 4 to column 16)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 605, column 4 to column 67)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 606, column 4 to column 26)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 600, column 77 to line 607, column 2)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 613, column 4 to column 53)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 614, column 4 to column 35)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 615, column 4 to column 40)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 616, column 4 to column 37)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 617, column 4 to column 47)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 618, column 4 to column 47)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 619, column 4 to column 36)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 620, column 4 to column 66)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 621, column 4 to column 14)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 625, column 8 to column 53)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 623, column 8 to column 39)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 622, column 4 to line 625, column 53)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 626, column 4 to column 17)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 611, column 16 to line 627, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 633, column 9 to column 24)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 634, column 9 to column 26)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 635, column 9 to column 19)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 636, column 9 to column 30)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 637, column 9 to column 19)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 638, column 9 to column 21)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 639, column 16 to column 17)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 639, column 9 to column 123)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 640, column 9 to column 159)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 644, column 18 to column 179)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 645, column 18 to column 78)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 646, column 18 to column 64)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 647, column 18 to column 80)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 648, column 18 to column 52)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 643, column 37 to line 649, column 17)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 643, column 16 to line 649, column 17)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 641, column 76 to line 650, column 10)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 641, column 9 to line 650, column 10)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 652, column 9 to column 48)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 653, column 9 to column 38)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 654, column 9 to column 59)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 655, column 9 to column 39)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 657, column 9 to column 41)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 658, column 8 to column 23)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 631, column 128 to line 659, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 663, column 4 to column 23)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 664, column 4 to column 22)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 665, column 4 to column 21)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 666, column 4 to column 42)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 667, column 4 to column 19)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 668, column 4 to column 45)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 669, column 11 to column 72)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 669, column 4 to line 670, column 122)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 671, column 4 to column 46)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 676, column 4 to column 30)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 678, column 6 to column 48)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 681, column 8 to column 23)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 680, column 6 to line 681, column 23)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 683, column 30 to column 39)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 683, column 6 to column 39)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 697, column 12 to column 69)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 698, column 12 to column 79)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 703, column 19 to column 133)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 702, column 16 to line 704, column 15)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 700, column 14 to column 51)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 699, column 43 to line 701, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 699, column 12 to line 704, column 15)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 705, column 13 to line 706, column 52)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 707, column 14 to column 41)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 708, column 14 to column 85)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 709, column 14 to column 61)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 696, column 28 to line 710, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 696, column 8 to line 710, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 694, column 10 to line 711, column 10)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 688, column 12 to column 144)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 689, column 12 to column 39)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 690, column 12 to column 48)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 691, column 12 to column 59)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 687, column 27 to line 692, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 687, column 8 to line 692, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 685, column 125 to line 693, column 7)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 685, column 6 to line 711, column 10)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 677, column 45 to line 712, column 7)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 677, column 4 to line 712, column 7)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 713, column 6 to column 22)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 661, column 138 to line 714, column 5)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 717, column 4 to column 44)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 718, column 4 to column 25)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 719, column 4 to column 22)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 720, column 4 to column 33)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 721, column 4 to column 28)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 722, column 4 to column 32)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 723, column 4 to column 39)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 725, column 4 to column 32)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 726, column 11 to column 30)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 726, column 4 to column 78)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 728, column 4 to column 128)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 732, column 8 to column 41)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 733, column 8 to column 45)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 734, column 8 to column 38)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 735, column 8 to column 22)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 737, column 19 to column 32)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 737, column 12 to column 58)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 738, column 19 to column 32)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 738, column 12 to column 43)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 739, column 19 to column 32)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 739, column 12 to column 40)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 741, column 19 to column 32)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 741, column 12 to column 48)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 744, column 15 to column 46)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 745, column 15 to column 41)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 746, column 15 to column 59)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 743, column 35 to line 747, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 743, column 12 to line 747, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 750, column 14 to column 63)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 751, column 14 to line 754, column 100)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 749, column 39 to line 756, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 749, column 12 to line 756, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 758, column 12 to line 759, column 120)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 760, column 12 to column 116)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 761, column 12 to column 40)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 763, column 12 to column 38)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 736, column 17 to line 764, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 736, column 8 to line 764, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 731, column 42 to line 765, column 5)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 731, column 4 to line 765, column 5)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 766, column 4 to column 58)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 716, column 93 to line 767, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 771, column 4 to column 21)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 772, column 4 to column 26)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 773, column 4 to column 26)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 774, column 4 to column 36)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 775, column 4 to column 24)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 776, column 4 to column 29)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 779, column 4 to column 22)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 782, column 9 to column 21)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 783, column 9 to column 46)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 784, column 23 to column 38)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 784, column 9 to column 72)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 785, column 16 to column 29)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 785, column 9 to column 72)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 788, column 11 to column 138)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 789, column 11 to line 791, column 55)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 787, column 18 to line 793, column 10)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 787, column 9 to line 793, column 10)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 781, column 42 to line 795, column 5)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 781, column 4 to line 795, column 5)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 796, column 4 to column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 769, column 144 to line 797, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 800, column 6 to column 19)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 801, column 6 to column 83)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 806, column 9 to column 56)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 805, column 6 to line 806, column 56)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 807, column 6 to column 21)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 799, column 140 to line 808, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 813, column 6 to column 23)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 814, column 13 to column 35)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 814, column 6 to column 44)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 818, column 8 to line 819, column 57)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 820, column 8 to column 62)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 821, column 17 to column 30)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 821, column 33 to column 46)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 821, column 8 to column 52)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 822, column 8 to column 66)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 826, column 8 to column 114)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 828, column 8 to column 42)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 831, column 8 to column 42)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 833, column 8 to column 38)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 817, column 40 to line 835, column 6)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 817, column 6 to line 835, column 6)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 836, column 4 to column 34)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 837, column 4 to column 19)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 812, column 102 to line 838, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 845, column 11 to column 12)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 845, column 4 to column 114)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 846, column 41 to column 42)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 846, column 4 to column 44)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 847, column 4 to column 22)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 848, column 4 to column 23)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 849, column 4 to column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 850, column 4 to column 42)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 851, column 4 to column 12)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 852, column 4 to column 21)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 853, column 4 to column 17)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 854, column 21 to column 22)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 854, column 4 to column 24)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 856, column 4 to column 109)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 857, column 4 to column 35)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 860, column 4 to column 15)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 862, column 10 to column 121)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 863, column 10 to column 138)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 865, column 11 to column 146)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 866, column 10 to column 58)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 867, column 10 to column 126)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 868, column 24 to column 41)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 868, column 10 to column 43)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 869, column 10 to column 41)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 870, column 17 to column 31)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 870, column 10 to column 87)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 872, column 17 to column 24)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 872, column 10 to column 74)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 873, column 10 to column 160)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 916, column 15 to line 918, column 169)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 924, column 15 to column 184)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 929, column 15 to column 198)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 931, column 15 to line 933, column 51)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 937, column 15 to line 939, column 148)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 911, column 14 to line 942, column 11)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 878, column 14 to line 880, column 80)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 883, column 15 to line 885, column 217)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 889, column 15 to column 156)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 892, column 15 to column 198)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 893, column 15 to column 135)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 895, column 15 to line 896, column 56)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 900, column 15 to line 902, column 51)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 906, column 14 to line 908, column 155)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 876, column 22 to line 910, column 11)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 876, column 10 to line 942, column 11)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 944, column 14 to line 946, column 208)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 943, column 33 to line 947, column 11)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 943, column 10 to line 947, column 11)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 861, column 19 to line 948, column 8)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 861, column 4 to line 948, column 8)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 949, column 4 to column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 844, column 4 to line 950, column 3)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 954, column 11 to column 32)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 954, column 4 to column 86)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 955, column 12 to column 31)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 955, column 4 to column 83)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 957, column 4 to column 14)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 958, column 4 to column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 959, column 4 to column 14)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 960, column 4 to column 17)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 961, column 4 to column 18)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 962, column 4 to column 26)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 963, column 4 to column 26)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 964, column 4 to column 36)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 965, column 26 to column 27)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 965, column 4 to column 44)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 966, column 4 to column 29)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 967, column 4 to column 21)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 968, column 4 to column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 969, column 4 to column 22)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 970, column 4 to column 65)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 972, column 9 to column 17)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 978, column 17 to column 84)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 980, column 17 to column 27)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 976, column 13 to line 981, column 14)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 975, column 13 to line 981, column 14)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 982, column 13 to column 19)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 974, column 9 to line 983, column 10)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 973, column 9 to line 983, column 10)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 971, column 17 to line 984, column 5)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 971, column 4 to line 984, column 5)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 986, column 4 to column 10)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 989, column 8 to column 27)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 990, column 8 to column 28)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 998, column 11 to column 80)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1005, column 16 to column 62)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1004, column 16 to line 1007, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1001, column 15 to column 82)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1002, column 15 to column 22)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1000, column 80 to line 1003, column 14)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1000, column 11 to line 1007, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 997, column 12 to line 1009, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 994, column 12 to column 53)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 995, column 13 to column 20)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 992, column 38 to line 996, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 992, column 8 to line 1009, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1015, column 11 to column 82)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1023, column 16 to column 62)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1022, column 15 to line 1025, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1019, column 14 to column 82)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1020, column 14 to column 21)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1018, column 82 to line 1021, column 12)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1018, column 11 to line 1025, column 13)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1014, column 12 to line 1026, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1011, column 11 to column 52)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1012, column 11 to column 18)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1010, column 39 to line 1013, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1010, column 8 to line 1026, column 9)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 988, column 4 to line 1027, column 5)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 987, column 4 to line 1027, column 5)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 1028, column 4 to column 27)",
 " (in 'stan_model_1_population_JC_genotypes.stan', line 952, column 62 to line 1029, column 3)"};


int
find_integer(const int& x, const std::vector<int>& values,
             std::ostream* pstream__) {
  using local_scalar_t__ = double;
  int current_statement__ = 0;
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    int i;
    i = std::numeric_limits<int>::min();
    
    current_statement__ = 124;
    i = 1;
    int result;
    result = std::numeric_limits<int>::min();
    
    current_statement__ = 125;
    result = -1;
    current_statement__ = 128;
    while ((primitive_value(logical_lte(i, num_elements(values))) &&
           primitive_value(
           logical_gt(
             stan::math::abs((rvalue(values, "values", index_uni(i)) - x)),
             0.00001)))) {
      current_statement__ = 126;
      i = (i + 1);
    }
    current_statement__ = 131;
    if ((primitive_value(logical_eq(num_elements(values), 0)) ||
        primitive_value(logical_gt(i, num_elements(values))))) {
      current_statement__ = 130;
      result = -1;
    } else {
      current_statement__ = 129;
      result = i;
    }
    current_statement__ = 132;
    return result;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct find_integer_functor__ {
int
operator()(const int& x, const std::vector<int>& values,
           std::ostream* pstream__)  const
{
return find_integer(x, values, pstream__);
}
};

template <typename T0__, typename T1__>
int
find(const T0__& x, const T1__& values_arg__, std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__,
          stan::value_type_t<T1__>>;
  int current_statement__ = 0;
  const auto& values = to_ref(values_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    int i;
    i = std::numeric_limits<int>::min();
    
    current_statement__ = 134;
    i = 1;
    int result;
    result = std::numeric_limits<int>::min();
    
    current_statement__ = 135;
    result = -1;
    current_statement__ = 138;
    while ((primitive_value(logical_lte(i, rows(values))) && primitive_value(
           logical_gt(
             stan::math::fabs((rvalue(values, "values", index_uni(i)) - x)),
             0.00001)))) {
      current_statement__ = 136;
      i = (i + 1);
    }
    current_statement__ = 141;
    if ((primitive_value(logical_eq(rows(values), 0)) || primitive_value(
        logical_gt(i, rows(values))))) {
      current_statement__ = 140;
      result = -1;
    } else {
      current_statement__ = 139;
      result = i;
    }
    current_statement__ = 142;
    return result;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct find_functor__ {
template <typename T0__, typename T1__>
int
operator()(const T0__& x, const T1__& values, std::ostream* pstream__)  const
{
return find(x, values, pstream__);
}
};

template <typename T0__, typename T1__>
int
find_sorted(const T0__& x, const T1__& sorted_values_arg__,
            std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__,
          stan::value_type_t<T1__>>;
  int current_statement__ = 0;
  const auto& sorted_values = to_ref(sorted_values_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    int i;
    i = std::numeric_limits<int>::min();
    
    current_statement__ = 144;
    i = 1;
    int result;
    result = std::numeric_limits<int>::min();
    
    current_statement__ = 145;
    result = -1;
    current_statement__ = 148;
    while ((primitive_value((primitive_value(
           logical_lte(i, rows(sorted_values))) && primitive_value(
           logical_gt(
             stan::math::fabs(
               (rvalue(sorted_values, "sorted_values", index_uni(i)) - x)),
             0.00001)))) && primitive_value(
           logical_lt(rvalue(sorted_values, "sorted_values", index_uni(i)),
             x)))) {
      current_statement__ = 146;
      i = (i + 1);
    }
    current_statement__ = 151;
    if ((primitive_value((primitive_value(logical_eq(rows(sorted_values), 0))
        || primitive_value(logical_gt(i, rows(sorted_values))))) ||
        primitive_value((primitive_value((primitive_value(logical_gte(i, 1))
        && primitive_value(logical_lte(i, rows(sorted_values))))) &&
        primitive_value(
        logical_gt(rvalue(sorted_values, "sorted_values", index_uni(i)), x)))))) {
      current_statement__ = 150;
      result = -1;
    } else {
      current_statement__ = 149;
      result = i;
    }
    current_statement__ = 152;
    return result;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct find_sorted_functor__ {
template <typename T0__, typename T1__>
int
operator()(const T0__& x, const T1__& sorted_values, std::ostream* pstream__)  const
{
return find_sorted(x, sorted_values, pstream__);
}
};

std::vector<int>
cumulative_sum_integer(const std::vector<int>& values, const int& N,
                       std::ostream* pstream__) {
  using local_scalar_t__ = double;
  int current_statement__ = 0;
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    current_statement__ = 154;
    validate_non_negative_index("result", "N", N);
    std::vector<int> result;
    result = std::vector<int>(N, std::numeric_limits<int>::min());
    
    int current_sum;
    current_sum = std::numeric_limits<int>::min();
    
    current_statement__ = 156;
    current_sum = 0;
    current_statement__ = 160;
    for (int i = 1; i <= N; ++i) {
      current_statement__ = 157;
      current_sum = (current_sum + rvalue(values, "values", index_uni(i)));
      current_statement__ = 158;
      assign(result, current_sum, "assigning variable result", index_uni(i));
    }
    current_statement__ = 161;
    return result;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct cumulative_sum_integer_functor__ {
std::vector<int>
operator()(const std::vector<int>& values, const int& N,
           std::ostream* pstream__)  const
{
return cumulative_sum_integer(values, N, pstream__);
}
};

template <typename T0__, typename T1__>
Eigen::Matrix<stan::promote_args_t<stan::value_type_t<T0__>,
stan::value_type_t<T1__>>, -1, 1>
difference_sorted_lists(const T0__& sorted_list_arg__,
                        const T1__& sorted_sublist_arg__,
                        std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<stan::value_type_t<T0__>,
          stan::value_type_t<T1__>>;
  int current_statement__ = 0;
  const auto& sorted_list = to_ref(sorted_list_arg__);
  const auto& sorted_sublist = to_ref(sorted_sublist_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    Eigen::Matrix<local_scalar_t__, -1, 1> result;
    result = Eigen::Matrix<local_scalar_t__, -1, 1>(0);
    stan::math::fill(result, DUMMY_VAR__);
    
    int idx_list;
    idx_list = std::numeric_limits<int>::min();
    
    current_statement__ = 164;
    idx_list = 1;
    int idx_sublist;
    idx_sublist = std::numeric_limits<int>::min();
    
    current_statement__ = 165;
    idx_sublist = 1;
    current_statement__ = 178;
    while (logical_lte(idx_list, rows(sorted_list))) {
      current_statement__ = 176;
      if (logical_lt(rvalue(sorted_list, "sorted_list", index_uni(idx_list)),
            rvalue(sorted_sublist, "sorted_sublist", index_uni(idx_sublist)))) {
        current_statement__ = 173;
        assign(result,
          append_row(stan::model::deep_copy(result),
            rvalue(sorted_list, "sorted_list", index_uni(idx_list))),
          "assigning variable result");
        current_statement__ = 174;
        idx_list = (idx_list + 1);
      } else {
        current_statement__ = 172;
        if (logical_gt(
              rvalue(sorted_list, "sorted_list", index_uni(idx_list)),
              rvalue(sorted_sublist, "sorted_sublist",
                index_uni(idx_sublist)))) {
          current_statement__ = 169;
          assign(result,
            append_row(stan::model::deep_copy(result),
              rvalue(sorted_list, "sorted_list", index_uni(idx_list))),
            "assigning variable result");
          current_statement__ = 170;
          idx_sublist = (idx_sublist + 1);
        } else {
          current_statement__ = 166;
          idx_list = (idx_list + 1);
          current_statement__ = 167;
          idx_sublist = (idx_sublist + 1);
        }
      }
    }
    current_statement__ = 179;
    return result;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct difference_sorted_lists_functor__ {
template <typename T0__, typename T1__>
Eigen::Matrix<stan::promote_args_t<stan::value_type_t<T0__>,
stan::value_type_t<T1__>>, -1, 1>
operator()(const T0__& sorted_list, const T1__& sorted_sublist,
           std::ostream* pstream__)  const
{
return difference_sorted_lists(sorted_list, sorted_sublist, pstream__);
}
};

template <typename T0__, typename T1__>
int
get_index_coal_time_above_MRCA_in_oldest_population_time_less_than(const T0__& time_in_oldest_population_time,
                                                                   const T1__& sorted_coal_times_in_oldest_population_time_arg__,
                                                                   std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__,
          stan::value_type_t<T1__>>;
  int current_statement__ = 0;
  const auto& sorted_coal_times_in_oldest_population_time = to_ref(sorted_coal_times_in_oldest_population_time_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    int idx;
    idx = std::numeric_limits<int>::min();
    
    current_statement__ = 182;
    idx = 1;
    current_statement__ = 185;
    while ((primitive_value(
           logical_lte(idx,
             rows(sorted_coal_times_in_oldest_population_time))) &&
           primitive_value(
           logical_lt(
             rvalue(sorted_coal_times_in_oldest_population_time,
               "sorted_coal_times_in_oldest_population_time", index_uni(idx)),
             time_in_oldest_population_time)))) {
      current_statement__ = 183;
      idx = (idx + 1);
    }
    current_statement__ = 186;
    return idx;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct get_index_coal_time_above_MRCA_in_oldest_population_time_less_than_functor__ {
template <typename T0__, typename T1__>
int
operator()(const T0__& time_in_oldest_population_time,
           const T1__& sorted_coal_times_in_oldest_population_time,
           std::ostream* pstream__)  const
{
return get_index_coal_time_above_MRCA_in_oldest_population_time_less_than(
         time_in_oldest_population_time,
         sorted_coal_times_in_oldest_population_time, pstream__);
}
};

template <typename T0__, typename T1__>
Eigen::Matrix<stan::promote_args_t<T0__,
stan::value_type_t<T1__>>, -1, 1>
get_coal_times_in_oldest_population_time_less_than(const T0__& time_in_oldest_population_time,
                                                   const T1__& sorted_coal_times_in_oldest_population_time_arg__,
                                                   std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__,
          stan::value_type_t<T1__>>;
  int current_statement__ = 0;
  const auto& sorted_coal_times_in_oldest_population_time = to_ref(sorted_coal_times_in_oldest_population_time_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    int idx;
    idx = std::numeric_limits<int>::min();
    
    current_statement__ = 188;
    idx = get_index_coal_time_above_MRCA_in_oldest_population_time_less_than(
            time_in_oldest_population_time,
            sorted_coal_times_in_oldest_population_time, pstream__);
    current_statement__ = 189;
    validate_non_negative_index("result", "idx - 1", (idx - 1));
    Eigen::Matrix<local_scalar_t__, -1, 1> result;
    result = Eigen::Matrix<local_scalar_t__, -1, 1>((idx - 1));
    stan::math::fill(result, DUMMY_VAR__);
    
    current_statement__ = 190;
    assign(result,
      segment(sorted_coal_times_in_oldest_population_time, 1, (idx - 1)),
      "assigning variable result");
    current_statement__ = 191;
    return result;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct get_coal_times_in_oldest_population_time_less_than_functor__ {
template <typename T0__, typename T1__>
Eigen::Matrix<stan::promote_args_t<T0__,
stan::value_type_t<T1__>>, -1, 1>
operator()(const T0__& time_in_oldest_population_time,
           const T1__& sorted_coal_times_in_oldest_population_time,
           std::ostream* pstream__)  const
{
return get_coal_times_in_oldest_population_time_less_than(
         time_in_oldest_population_time,
         sorted_coal_times_in_oldest_population_time, pstream__);
}
};

template <typename T0__, typename T1__>
Eigen::Matrix<stan::promote_args_t<T0__,
stan::value_type_t<T1__>>, -1, 1>
get_inter_coal_times_in_oldest_population_time_less_than(const T0__& time_in_oldest_population_time,
                                                         const T1__& sorted_coal_times_in_oldest_population_time_arg__,
                                                         std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__,
          stan::value_type_t<T1__>>;
  int current_statement__ = 0;
  const auto& sorted_coal_times_in_oldest_population_time = to_ref(sorted_coal_times_in_oldest_population_time_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    Eigen::Matrix<local_scalar_t__, -1, 1> result;
    result = Eigen::Matrix<local_scalar_t__, -1, 1>(0);
    stan::math::fill(result, DUMMY_VAR__);
    
    int idx;
    idx = std::numeric_limits<int>::min();
    
    local_scalar_t__ inter_coal_time;
    inter_coal_time = DUMMY_VAR__;
    
    current_statement__ = 196;
    idx = 1;
    current_statement__ = 205;
    while ((primitive_value(
           logical_lte(idx,
             rows(sorted_coal_times_in_oldest_population_time))) &&
           primitive_value(
           logical_lt(
             rvalue(sorted_coal_times_in_oldest_population_time,
               "sorted_coal_times_in_oldest_population_time", index_uni(idx)),
             time_in_oldest_population_time)))) {
      current_statement__ = 201;
      if (logical_eq(idx, 1)) {
        current_statement__ = 199;
        inter_coal_time = rvalue(sorted_coal_times_in_oldest_population_time,
                            "sorted_coal_times_in_oldest_population_time",
                            index_uni(idx));
      } else {
        current_statement__ = 197;
        inter_coal_time = (rvalue(
                             sorted_coal_times_in_oldest_population_time,
                             "sorted_coal_times_in_oldest_population_time",
                             index_uni(idx)) -
                            rvalue(
                              sorted_coal_times_in_oldest_population_time,
                              "sorted_coal_times_in_oldest_population_time",
                              index_uni((idx - 1))));
      }
      current_statement__ = 202;
      assign(result,
        append_row(stan::model::deep_copy(result),
          rvalue(sorted_coal_times_in_oldest_population_time,
            "sorted_coal_times_in_oldest_population_time", index_uni(idx))),
        "assigning variable result");
      current_statement__ = 203;
      idx = (idx + 1);
    }
    current_statement__ = 206;
    return result;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct get_inter_coal_times_in_oldest_population_time_less_than_functor__ {
template <typename T0__, typename T1__>
Eigen::Matrix<stan::promote_args_t<T0__,
stan::value_type_t<T1__>>, -1, 1>
operator()(const T0__& time_in_oldest_population_time,
           const T1__& sorted_coal_times_in_oldest_population_time,
           std::ostream* pstream__)  const
{
return get_inter_coal_times_in_oldest_population_time_less_than(
         time_in_oldest_population_time,
         sorted_coal_times_in_oldest_population_time, pstream__);
}
};

template <typename T0__, typename T1__>
Eigen::Matrix<stan::promote_args_t<T0__,
stan::value_type_t<T1__>>, -1, 1>
get_candidate_time_MRCAs(const T0__& time_origin_in_oldest_population_time,
                         const T1__& sorted_coal_times_in_oldest_population_time_arg__,
                         const std::vector<std::vector<int>>& topology,
                         const std::vector<int>& number_of_tips_below,
                         const int& idx_first_coal_time_above_torigin,
                         std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__,
          stan::value_type_t<T1__>>;
  int current_statement__ = 0;
  const auto& sorted_coal_times_in_oldest_population_time = to_ref(sorted_coal_times_in_oldest_population_time_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    local_scalar_t__ result;
    result = DUMMY_VAR__;
    
    int number_candidate_MRCAs;
    number_candidate_MRCAs = std::numeric_limits<int>::min();
    
    current_statement__ = 209;
    number_candidate_MRCAs = 0;
    int idx_left;
    idx_left = std::numeric_limits<int>::min();
    
    int idx_right;
    idx_right = std::numeric_limits<int>::min();
    
    int idx_left_below;
    idx_left_below = std::numeric_limits<int>::min();
    
    int idx_right_below;
    idx_right_below = std::numeric_limits<int>::min();
    
    int idx_MRCA;
    idx_MRCA = std::numeric_limits<int>::min();
    
    int current_pos;
    current_pos = std::numeric_limits<int>::min();
    
    current_statement__ = 215;
    current_pos = 1;
    current_statement__ = 216;
    validate_non_negative_index("indexes_candidate_MRCAs",
                                "idx_first_coal_time_above_torigin - 1",
                                (idx_first_coal_time_above_torigin - 1));
    Eigen::Matrix<local_scalar_t__, -1, 1> indexes_candidate_MRCAs;
    indexes_candidate_MRCAs = Eigen::Matrix<local_scalar_t__, -1, 1>(
      (idx_first_coal_time_above_torigin - 1));
    stan::math::fill(indexes_candidate_MRCAs, DUMMY_VAR__);
    
    current_statement__ = 218;
    validate_non_negative_index("coal_time_candidate_MRCAs",
                                "idx_first_coal_time_above_torigin - 1",
                                (idx_first_coal_time_above_torigin - 1));
    Eigen::Matrix<local_scalar_t__, -1, 1> coal_time_candidate_MRCAs;
    coal_time_candidate_MRCAs = Eigen::Matrix<local_scalar_t__, -1, 1>(
      (idx_first_coal_time_above_torigin - 1));
    stan::math::fill(coal_time_candidate_MRCAs, DUMMY_VAR__);
    
    current_statement__ = 239;
    for (int i = idx_first_coal_time_above_torigin;
         i <= rows(sorted_coal_times_in_oldest_population_time); ++i) {
      current_statement__ = 220;
      idx_left = rvalue(topology, "topology", index_uni(i), index_uni(2));
      current_statement__ = 221;
      idx_right = rvalue(topology, "topology", index_uni(i), index_uni(3));
      current_statement__ = 222;
      idx_left_below = find_integer(idx_left,
                         rvalue(topology, "topology",
                           index_min_max(1, (idx_first_coal_time_above_torigin
                                              - 1)), index_uni(1)), pstream__);
      current_statement__ = 223;
      idx_right_below = find_integer(idx_right,
                          rvalue(topology, "topology",
                            index_min_max(1, (idx_first_coal_time_above_torigin
                                               - 1)), index_uni(1)), pstream__);
      current_statement__ = 230;
      if (logical_neq(idx_left_below, -1)) {
        current_statement__ = 225;
        number_candidate_MRCAs = (number_candidate_MRCAs + 1);
        current_statement__ = 226;
        assign(indexes_candidate_MRCAs, idx_left,
          "assigning variable indexes_candidate_MRCAs", index_uni(current_pos));
        current_statement__ = 227;
        assign(coal_time_candidate_MRCAs,
          rvalue(sorted_coal_times_in_oldest_population_time,
            "sorted_coal_times_in_oldest_population_time",
            index_uni(idx_left_below)),
          "assigning variable coal_time_candidate_MRCAs", index_uni(current_pos));
        current_statement__ = 228;
        current_pos = (current_pos + 1);
      } else {
        
      }
      current_statement__ = 237;
      if (logical_neq(idx_right_below, -1)) {
        current_statement__ = 232;
        number_candidate_MRCAs = (number_candidate_MRCAs + 1);
        current_statement__ = 233;
        assign(indexes_candidate_MRCAs, idx_right,
          "assigning variable indexes_candidate_MRCAs", index_uni(current_pos));
        current_statement__ = 234;
        assign(coal_time_candidate_MRCAs,
          rvalue(sorted_coal_times_in_oldest_population_time,
            "sorted_coal_times_in_oldest_population_time",
            index_uni(idx_right_below)),
          "assigning variable coal_time_candidate_MRCAs", index_uni(current_pos));
        current_statement__ = 235;
        current_pos = (current_pos + 1);
      } else {
        
      }
    }
    current_statement__ = 240;
    return coal_time_candidate_MRCAs;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct get_candidate_time_MRCAs_functor__ {
template <typename T0__, typename T1__>
Eigen::Matrix<stan::promote_args_t<T0__,
stan::value_type_t<T1__>>, -1, 1>
operator()(const T0__& time_origin_in_oldest_population_time,
           const T1__& sorted_coal_times_in_oldest_population_time,
           const std::vector<std::vector<int>>& topology,
           const std::vector<int>& number_of_tips_below,
           const int& idx_first_coal_time_above_torigin,
           std::ostream* pstream__)  const
{
return get_candidate_time_MRCAs(time_origin_in_oldest_population_time,
         sorted_coal_times_in_oldest_population_time, topology,
         number_of_tips_below, idx_first_coal_time_above_torigin, pstream__);
}
};

template <typename T0__, typename T1__, typename T2__>
stan::promote_args_t<T0__, T1__,
stan::value_type_t<T2__>>
random_time_MRCA_below_time_origin(const T0__& value,
                                   const T1__& time_origin_in_oldest_population_time,
                                   const T2__& sorted_coal_times_in_oldest_population_time_arg__,
                                   const std::vector<std::vector<int>>& topology,
                                   const std::vector<int>& number_of_tips_below,
                                   std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__,
          T1__,
          stan::value_type_t<T2__>>;
  int current_statement__ = 0;
  const auto& sorted_coal_times_in_oldest_population_time = to_ref(sorted_coal_times_in_oldest_population_time_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    local_scalar_t__ result;
    result = DUMMY_VAR__;
    
    int number_candidate_MRCAs;
    number_candidate_MRCAs = std::numeric_limits<int>::min();
    
    current_statement__ = 243;
    number_candidate_MRCAs = 0;
    int idx_left;
    idx_left = std::numeric_limits<int>::min();
    
    int idx_right;
    idx_right = std::numeric_limits<int>::min();
    
    int idx_left_below;
    idx_left_below = std::numeric_limits<int>::min();
    
    int idx_right_below;
    idx_right_below = std::numeric_limits<int>::min();
    
    int idx_MRCA;
    idx_MRCA = std::numeric_limits<int>::min();
    
    int current_pos;
    current_pos = std::numeric_limits<int>::min();
    
    current_statement__ = 249;
    current_pos = 1;
    int idx_first_coal_time_above_torigin;
    idx_first_coal_time_above_torigin = std::numeric_limits<int>::min();
    
    current_statement__ = 250;
    idx_first_coal_time_above_torigin = get_index_coal_time_above_MRCA_in_oldest_population_time_less_than(
                                          time_origin_in_oldest_population_time,
                                          sorted_coal_times_in_oldest_population_time, pstream__);
    current_statement__ = 251;
    validate_non_negative_index("indexes_nodes_below_torigin",
                                "idx_first_coal_time_above_torigin - 1",
                                (idx_first_coal_time_above_torigin - 1));
    std::vector<int> indexes_nodes_below_torigin;
    indexes_nodes_below_torigin = std::vector<int>((idx_first_coal_time_above_torigin
                                                     - 1), std::numeric_limits<int>::min());
    
    
    current_statement__ = 252;
    assign(indexes_nodes_below_torigin,
      rvalue(topology, "topology",
        index_min_max(1, (idx_first_coal_time_above_torigin - 1)),
          index_uni(1)), "assigning variable indexes_nodes_below_torigin");
    current_statement__ = 253;
    validate_non_negative_index("indexes_candidate_MRCAs",
                                "idx_first_coal_time_above_torigin - 1",
                                (idx_first_coal_time_above_torigin - 1));
    Eigen::Matrix<local_scalar_t__, -1, 1> indexes_candidate_MRCAs;
    indexes_candidate_MRCAs = Eigen::Matrix<local_scalar_t__, -1, 1>(
      (idx_first_coal_time_above_torigin - 1));
    stan::math::fill(indexes_candidate_MRCAs, DUMMY_VAR__);
    
    current_statement__ = 255;
    validate_non_negative_index("coal_time_candidate_MRCAs",
                                "idx_first_coal_time_above_torigin - 1",
                                (idx_first_coal_time_above_torigin - 1));
    Eigen::Matrix<local_scalar_t__, -1, 1> coal_time_candidate_MRCAs;
    coal_time_candidate_MRCAs = Eigen::Matrix<local_scalar_t__, -1, 1>(
      (idx_first_coal_time_above_torigin - 1));
    stan::math::fill(coal_time_candidate_MRCAs, DUMMY_VAR__);
    
    current_statement__ = 288;
    for (int i = idx_first_coal_time_above_torigin;
         i <= rows(sorted_coal_times_in_oldest_population_time); ++i) {
      current_statement__ = 257;
      idx_left = rvalue(topology, "topology", index_uni(i), index_uni(2));
      current_statement__ = 258;
      idx_right = rvalue(topology, "topology", index_uni(i), index_uni(3));
      current_statement__ = 259;
      idx_left_below = find_integer(idx_left,
                         rvalue(topology, "topology",
                           index_min_max(1, (idx_first_coal_time_above_torigin
                                              - 1)), index_uni(1)), pstream__);
      current_statement__ = 260;
      idx_right_below = find_integer(idx_right,
                          rvalue(topology, "topology",
                            index_min_max(1, (idx_first_coal_time_above_torigin
                                               - 1)), index_uni(1)), pstream__);
      current_statement__ = 273;
      if (logical_neq(idx_left_below, -1)) {
        current_statement__ = 268;
        number_candidate_MRCAs = (number_candidate_MRCAs + 1);
        current_statement__ = 269;
        assign(indexes_candidate_MRCAs, idx_left,
          "assigning variable indexes_candidate_MRCAs", index_uni(current_pos));
        current_statement__ = 270;
        assign(coal_time_candidate_MRCAs,
          rvalue(sorted_coal_times_in_oldest_population_time,
            "sorted_coal_times_in_oldest_population_time",
            index_uni(idx_left_below)),
          "assigning variable coal_time_candidate_MRCAs", index_uni(current_pos));
        current_statement__ = 271;
        current_pos = (current_pos + 1);
      } else {
        current_statement__ = 266;
        if (logical_eq(
              rvalue(number_of_tips_below, "number_of_tips_below",
                index_uni(idx_left)), 1)) {
          current_statement__ = 261;
          number_candidate_MRCAs = (number_candidate_MRCAs + 1);
          current_statement__ = 262;
          assign(indexes_candidate_MRCAs, idx_left,
            "assigning variable indexes_candidate_MRCAs", index_uni(current_pos));
          current_statement__ = 263;
          assign(coal_time_candidate_MRCAs, 0.0,
            "assigning variable coal_time_candidate_MRCAs", index_uni(current_pos));
          current_statement__ = 264;
          current_pos = (current_pos + 1);
        }
      }
      current_statement__ = 286;
      if (logical_neq(idx_right_below, -1)) {
        current_statement__ = 281;
        number_candidate_MRCAs = (number_candidate_MRCAs + 1);
        current_statement__ = 282;
        assign(indexes_candidate_MRCAs, idx_right,
          "assigning variable indexes_candidate_MRCAs", index_uni(current_pos));
        current_statement__ = 283;
        assign(coal_time_candidate_MRCAs,
          rvalue(sorted_coal_times_in_oldest_population_time,
            "sorted_coal_times_in_oldest_population_time",
            index_uni(idx_right_below)),
          "assigning variable coal_time_candidate_MRCAs", index_uni(current_pos));
        current_statement__ = 284;
        current_pos = (current_pos + 1);
      } else {
        current_statement__ = 279;
        if (logical_eq(
              rvalue(number_of_tips_below, "number_of_tips_below",
                index_uni(idx_right)), 1)) {
          current_statement__ = 274;
          number_candidate_MRCAs = (number_candidate_MRCAs + 1);
          current_statement__ = 275;
          assign(indexes_candidate_MRCAs, idx_right,
            "assigning variable indexes_candidate_MRCAs", index_uni(current_pos));
          current_statement__ = 276;
          assign(coal_time_candidate_MRCAs, 0.0,
            "assigning variable coal_time_candidate_MRCAs", index_uni(current_pos));
          current_statement__ = 277;
          current_pos = (current_pos + 1);
        }
      }
    }
    current_statement__ = 291;
    if (logical_gt(number_candidate_MRCAs, 0)) {
      current_statement__ = 290;
      result = -stan::math::log(number_candidate_MRCAs);
    } else {
      current_statement__ = 289;
      std::stringstream errmsg_stream__;
      errmsg_stream__ << "no candidate MRCAs!";
      throw std::domain_error(errmsg_stream__.str());
    }
    current_statement__ = 292;
    return result;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct random_time_MRCA_below_time_origin_functor__ {
template <typename T0__, typename T1__, typename T2__>
stan::promote_args_t<T0__, T1__,
stan::value_type_t<T2__>>
operator()(const T0__& value,
           const T1__& time_origin_in_oldest_population_time,
           const T2__& sorted_coal_times_in_oldest_population_time,
           const std::vector<std::vector<int>>& topology,
           const std::vector<int>& number_of_tips_below,
           std::ostream* pstream__)  const
{
return random_time_MRCA_below_time_origin(value,
         time_origin_in_oldest_population_time,
         sorted_coal_times_in_oldest_population_time, topology,
         number_of_tips_below, pstream__);
}
};

template <typename T4__, typename T5__, typename T6__, typename T7__>
Eigen::Matrix<stan::promote_args_t<T4__, stan::value_type_t<T5__>, T6__,
stan::value_type_t<T7__>>, -1, 1>
get_list_coal_times_in_oldest_population_time_below_time_origin(const int& idx_pop,
                                                                const int& order,
                                                                const int& total_sample_size,
                                                                const int& number_internal_nodes_under_MRCA,
                                                                const T4__& coal_time_MRCA_in_oldest_population,
                                                                const T5__& inmigrants_MRCA_times_in_model_time_oldest_population_arg__,
                                                                const T6__& time_origin_in_oldest_population_time,
                                                                const T7__& sorted_coal_times_in_oldest_population_time_arg__,
                                                                const std::vector<std::vector<int>>& topology,
                                                                const std::vector<int>& number_of_coalescent_times_below,
                                                                std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T4__,
          stan::value_type_t<T5__>,
          T6__,
          stan::value_type_t<T7__>>;
  int current_statement__ = 0;
  const auto& inmigrants_MRCA_times_in_model_time_oldest_population = to_ref(inmigrants_MRCA_times_in_model_time_oldest_population_arg__);
  const auto& sorted_coal_times_in_oldest_population_time = to_ref(sorted_coal_times_in_oldest_population_time_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    current_statement__ = 294;
    validate_non_negative_index("subtree_coal_times_below",
                                "rows(sorted_coal_times_in_oldest_population_time)",
                                rows(
                                  sorted_coal_times_in_oldest_population_time));
    Eigen::Matrix<local_scalar_t__, -1, 1> subtree_coal_times_below;
    subtree_coal_times_below = Eigen::Matrix<local_scalar_t__, -1, 1>(
      rows(sorted_coal_times_in_oldest_population_time));
    stan::math::fill(subtree_coal_times_below, DUMMY_VAR__);
    
    int idx;
    idx = std::numeric_limits<int>::min();
    
    int idx_subtree;
    idx_subtree = std::numeric_limits<int>::min();
    
    current_statement__ = 297;
    idx_subtree = find(coal_time_MRCA_in_oldest_population,
                    sorted_coal_times_in_oldest_population_time, pstream__);
    current_statement__ = 298;
    validate_non_negative_index("nodes_to_visit", "idx_subtree", idx_subtree);
    std::vector<int> nodes_to_visit;
    nodes_to_visit = std::vector<int>(idx_subtree, std::numeric_limits<int>::min());
    
    
    int idx_MRCA;
    idx_MRCA = std::numeric_limits<int>::min();
    
    int left_child;
    left_child = std::numeric_limits<int>::min();
    
    int right_child;
    right_child = std::numeric_limits<int>::min();
    
    local_scalar_t__ coal_time_above;
    coal_time_above = DUMMY_VAR__;
    
    int pos_subtree_coal;
    pos_subtree_coal = std::numeric_limits<int>::min();
    
    current_statement__ = 304;
    pos_subtree_coal = 1;
    int pos_subtree;
    pos_subtree = std::numeric_limits<int>::min();
    
    current_statement__ = 305;
    pos_subtree = 1;
    current_statement__ = 306;
    assign(nodes_to_visit,
      rvalue(topology, "topology", index_uni(idx_subtree), index_uni(1)),
      "assigning variable nodes_to_visit", index_uni(pos_subtree));
    current_statement__ = 307;
    pos_subtree = (pos_subtree + 1);
    current_statement__ = 337;
    while ((primitive_value(logical_gte(idx_subtree, 1)) && primitive_value(
           logical_lte(pos_subtree_coal, number_internal_nodes_under_MRCA)))) {
      current_statement__ = 334;
      if (logical_gt(order, 1)) {
        current_statement__ = 332;
        if ((primitive_value((primitive_value(
            logical_gt(
              rows(inmigrants_MRCA_times_in_model_time_oldest_population), 0))
            && primitive_value(
            logical_eq(
              find(
                rvalue(sorted_coal_times_in_oldest_population_time,
                  "sorted_coal_times_in_oldest_population_time",
                  index_uni(idx_subtree)),
                inmigrants_MRCA_times_in_model_time_oldest_population, pstream__),
              -1)))) && primitive_value(
            logical_neq(
              find_integer(
                rvalue(topology, "topology",
                  index_uni(idx_subtree), index_uni(1)),
                nodes_to_visit, pstream__), -1)))) {
          current_statement__ = 321;
          assign(subtree_coal_times_below,
            rvalue(sorted_coal_times_in_oldest_population_time,
              "sorted_coal_times_in_oldest_population_time",
              index_uni(idx_subtree)),
            "assigning variable subtree_coal_times_below", index_uni(pos_subtree_coal));
          current_statement__ = 325;
          if (logical_gt(
                rvalue(topology, "topology",
                  index_uni(idx_subtree), index_uni(2)), total_sample_size)) {
            current_statement__ = 322;
            assign(nodes_to_visit,
              rvalue(topology, "topology",
                index_uni(idx_subtree), index_uni(2)),
              "assigning variable nodes_to_visit", index_uni(pos_subtree));
            current_statement__ = 323;
            pos_subtree = (pos_subtree + 1);
          }
          current_statement__ = 329;
          if (logical_gt(
                rvalue(topology, "topology",
                  index_uni(idx_subtree), index_uni(3)), total_sample_size)) {
            current_statement__ = 326;
            assign(nodes_to_visit,
              rvalue(topology, "topology",
                index_uni(idx_subtree), index_uni(3)),
              "assigning variable nodes_to_visit", index_uni(pos_subtree));
            current_statement__ = 327;
            pos_subtree = (pos_subtree + 1);
          }
          current_statement__ = 330;
          pos_subtree_coal = (pos_subtree_coal + 1);
        }
      } else {
        current_statement__ = 319;
        if (logical_neq(
              find_integer(
                rvalue(topology, "topology",
                  index_uni(idx_subtree), index_uni(1)),
                nodes_to_visit, pstream__), -1)) {
          current_statement__ = 308;
          assign(subtree_coal_times_below,
            rvalue(sorted_coal_times_in_oldest_population_time,
              "sorted_coal_times_in_oldest_population_time",
              index_uni(idx_subtree)),
            "assigning variable subtree_coal_times_below", index_uni(pos_subtree_coal));
          current_statement__ = 312;
          if (logical_gt(
                rvalue(topology, "topology",
                  index_uni(idx_subtree), index_uni(2)), total_sample_size)) {
            current_statement__ = 309;
            assign(nodes_to_visit,
              rvalue(topology, "topology",
                index_uni(idx_subtree), index_uni(2)),
              "assigning variable nodes_to_visit", index_uni(pos_subtree));
            current_statement__ = 310;
            pos_subtree = (pos_subtree + 1);
          }
          current_statement__ = 316;
          if (logical_gt(
                rvalue(topology, "topology",
                  index_uni(idx_subtree), index_uni(3)), total_sample_size)) {
            current_statement__ = 313;
            assign(nodes_to_visit,
              rvalue(topology, "topology",
                index_uni(idx_subtree), index_uni(3)),
              "assigning variable nodes_to_visit", index_uni(pos_subtree));
            current_statement__ = 314;
            pos_subtree = (pos_subtree + 1);
          }
          current_statement__ = 317;
          pos_subtree_coal = (pos_subtree_coal + 1);
        }
      }
      current_statement__ = 335;
      idx_subtree = (idx_subtree - 1);
    }
    current_statement__ = 338;
    return sort_asc(
             segment(subtree_coal_times_below, 1, (pos_subtree_coal - 1)));
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct get_list_coal_times_in_oldest_population_time_below_time_origin_functor__ {
template <typename T4__, typename T5__, typename T6__, typename T7__>
Eigen::Matrix<stan::promote_args_t<T4__, stan::value_type_t<T5__>, T6__,
stan::value_type_t<T7__>>, -1, 1>
operator()(const int& idx_pop, const int& order,
           const int& total_sample_size,
           const int& number_internal_nodes_under_MRCA,
           const T4__& coal_time_MRCA_in_oldest_population,
           const T5__& inmigrants_MRCA_times_in_model_time_oldest_population,
           const T6__& time_origin_in_oldest_population_time,
           const T7__& sorted_coal_times_in_oldest_population_time,
           const std::vector<std::vector<int>>& topology,
           const std::vector<int>& number_of_coalescent_times_below,
           std::ostream* pstream__)  const
{
return get_list_coal_times_in_oldest_population_time_below_time_origin(
         idx_pop, order, total_sample_size, number_internal_nodes_under_MRCA,
         coal_time_MRCA_in_oldest_population,
         inmigrants_MRCA_times_in_model_time_oldest_population,
         time_origin_in_oldest_population_time,
         sorted_coal_times_in_oldest_population_time, topology,
         number_of_coalescent_times_below, pstream__);
}
};

template <typename T2__, typename T3__, typename T4__, typename T5__>
std::vector<int>
get_tips_positions_below_time(const int& order,
                              const int& number_tips_under_MRCA,
                              const T2__& coal_time_MRCA_in_oldest_population,
                              const T3__& inmigrants_MRCA_times_in_model_time_oldest_population_arg__,
                              const T4__& time_origin_in_oldest_population_time,
                              const T5__& sorted_coal_times_in_oldest_population_time_arg__,
                              const std::vector<std::vector<int>>& topology,
                              const int& total_sample_size,
                              std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T2__,
          stan::value_type_t<T3__>,
          T4__,
          stan::value_type_t<T5__>>;
  int current_statement__ = 0;
  const auto& inmigrants_MRCA_times_in_model_time_oldest_population = to_ref(inmigrants_MRCA_times_in_model_time_oldest_population_arg__);
  const auto& sorted_coal_times_in_oldest_population_time = to_ref(sorted_coal_times_in_oldest_population_time_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    current_statement__ = 340;
    validate_non_negative_index("result", "total_sample_size",
                                total_sample_size);
    std::vector<int> result;
    result = std::vector<int>(total_sample_size, std::numeric_limits<int>::min());
    
    
    int idx_left;
    idx_left = std::numeric_limits<int>::min();
    
    int idx_right;
    idx_right = std::numeric_limits<int>::min();
    
    int idx_subtree;
    idx_subtree = std::numeric_limits<int>::min();
    
    current_statement__ = 344;
    idx_subtree = find(coal_time_MRCA_in_oldest_population,
                    sorted_coal_times_in_oldest_population_time, pstream__);
    current_statement__ = 345;
    validate_non_negative_index("tips_positions", "idx_subtree + 1",
                                (idx_subtree + 1));
    std::vector<int> tips_positions;
    tips_positions = std::vector<int>((idx_subtree + 1), std::numeric_limits<int>::min());
    
    
    current_statement__ = 347;
    validate_non_negative_index("nodes_to_visit", "idx_subtree + 1",
                                (idx_subtree + 1));
    std::vector<int> nodes_to_visit;
    nodes_to_visit = std::vector<int>((idx_subtree + 1), std::numeric_limits<int>::min());
    
    
    int idx_MRCA;
    idx_MRCA = std::numeric_limits<int>::min();
    
    int left_child;
    left_child = std::numeric_limits<int>::min();
    
    int right_child;
    right_child = std::numeric_limits<int>::min();
    
    local_scalar_t__ coal_time_above;
    coal_time_above = DUMMY_VAR__;
    
    int pos_subtree;
    pos_subtree = std::numeric_limits<int>::min();
    
    current_statement__ = 353;
    pos_subtree = 1;
    current_statement__ = 354;
    validate_non_negative_index("indexes_nodes_below_torigin",
                                "idx_subtree - 1", (idx_subtree - 1));
    std::vector<int> indexes_nodes_below_torigin;
    indexes_nodes_below_torigin = std::vector<int>((idx_subtree - 1), std::numeric_limits<int>::min());
    
    
    current_statement__ = 355;
    assign(indexes_nodes_below_torigin,
      rvalue(topology, "topology",
        index_min_max(1, (idx_subtree - 1)), index_uni(1)),
      "assigning variable indexes_nodes_below_torigin");
    int pos;
    pos = std::numeric_limits<int>::min();
    
    current_statement__ = 356;
    pos = 1;
    current_statement__ = 357;
    assign(nodes_to_visit,
      rvalue(topology, "topology", index_uni(idx_subtree), index_uni(1)),
      "assigning variable nodes_to_visit", index_uni(pos_subtree));
    current_statement__ = 358;
    pos_subtree = (pos_subtree + 1);
    current_statement__ = 396;
    while ((primitive_value(logical_gte(idx_subtree, 1)) && primitive_value(
           logical_lte(pos, number_tips_under_MRCA)))) {
      current_statement__ = 393;
      if (logical_gt(order, 1)) {
        current_statement__ = 391;
        if ((primitive_value((primitive_value(
            logical_gt(
              rows(inmigrants_MRCA_times_in_model_time_oldest_population), 0))
            && primitive_value(
            logical_eq(
              find(
                rvalue(sorted_coal_times_in_oldest_population_time,
                  "sorted_coal_times_in_oldest_population_time",
                  index_uni(idx_subtree)),
                inmigrants_MRCA_times_in_model_time_oldest_population, pstream__),
              -1)))) && primitive_value(
            logical_neq(
              find_integer(
                rvalue(topology, "topology",
                  index_uni(idx_subtree), index_uni(1)),
                nodes_to_visit, pstream__), -1)))) {
          current_statement__ = 382;
          if (logical_gt(
                rvalue(topology, "topology",
                  index_uni(idx_subtree), index_uni(2)), total_sample_size)) {
            current_statement__ = 379;
            assign(nodes_to_visit,
              rvalue(topology, "topology",
                index_uni(idx_subtree), index_uni(2)),
              "assigning variable nodes_to_visit", index_uni(pos_subtree));
            current_statement__ = 380;
            pos_subtree = (pos_subtree + 1);
          } else {
            current_statement__ = 376;
            assign(tips_positions,
              rvalue(topology, "topology",
                index_uni(idx_subtree), index_uni(2)),
              "assigning variable tips_positions", index_uni(pos));
            current_statement__ = 377;
            pos = (pos + 1);
          }
          current_statement__ = 389;
          if (logical_gt(
                rvalue(topology, "topology",
                  index_uni(idx_subtree), index_uni(3)), total_sample_size)) {
            current_statement__ = 386;
            assign(nodes_to_visit,
              rvalue(topology, "topology",
                index_uni(idx_subtree), index_uni(3)),
              "assigning variable nodes_to_visit", index_uni(pos_subtree));
            current_statement__ = 387;
            pos_subtree = (pos_subtree + 1);
          } else {
            current_statement__ = 383;
            assign(tips_positions,
              rvalue(topology, "topology",
                index_uni(idx_subtree), index_uni(3)),
              "assigning variable tips_positions", index_uni(pos));
            current_statement__ = 384;
            pos = (pos + 1);
          }
        }
      } else {
        current_statement__ = 374;
        if (logical_neq(
              find_integer(
                rvalue(topology, "topology",
                  index_uni(idx_subtree), index_uni(1)),
                nodes_to_visit, pstream__), -1)) {
          current_statement__ = 365;
          if (logical_gt(
                rvalue(topology, "topology",
                  index_uni(idx_subtree), index_uni(2)), total_sample_size)) {
            current_statement__ = 362;
            assign(nodes_to_visit,
              rvalue(topology, "topology",
                index_uni(idx_subtree), index_uni(2)),
              "assigning variable nodes_to_visit", index_uni(pos_subtree));
            current_statement__ = 363;
            pos_subtree = (pos_subtree + 1);
          } else {
            current_statement__ = 359;
            assign(tips_positions,
              rvalue(topology, "topology",
                index_uni(idx_subtree), index_uni(2)),
              "assigning variable tips_positions", index_uni(pos));
            current_statement__ = 360;
            pos = (pos + 1);
          }
          current_statement__ = 372;
          if (logical_gt(
                rvalue(topology, "topology",
                  index_uni(idx_subtree), index_uni(3)), total_sample_size)) {
            current_statement__ = 369;
            assign(nodes_to_visit,
              rvalue(topology, "topology",
                index_uni(idx_subtree), index_uni(3)),
              "assigning variable nodes_to_visit", index_uni(pos_subtree));
            current_statement__ = 370;
            pos_subtree = (pos_subtree + 1);
          } else {
            current_statement__ = 366;
            assign(tips_positions,
              rvalue(topology, "topology",
                index_uni(idx_subtree), index_uni(3)),
              "assigning variable tips_positions", index_uni(pos));
            current_statement__ = 367;
            pos = (pos + 1);
          }
        }
      }
      current_statement__ = 394;
      idx_subtree = (idx_subtree - 1);
    }
    current_statement__ = 397;
    assign(result,
      append_array(
        rvalue(tips_positions, "tips_positions", index_min_max(1, (pos - 1))),
        rep_array(-1, ((total_sample_size - pos) + 1))),
      "assigning variable result");
    current_statement__ = 398;
    return result;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct get_tips_positions_below_time_functor__ {
template <typename T2__, typename T3__, typename T4__, typename T5__>
std::vector<int>
operator()(const int& order, const int& number_tips_under_MRCA,
           const T2__& coal_time_MRCA_in_oldest_population,
           const T3__& inmigrants_MRCA_times_in_model_time_oldest_population,
           const T4__& time_origin_in_oldest_population_time,
           const T5__& sorted_coal_times_in_oldest_population_time,
           const std::vector<std::vector<int>>& topology,
           const int& total_sample_size, std::ostream* pstream__)  const
{
return get_tips_positions_below_time(order, number_tips_under_MRCA,
         coal_time_MRCA_in_oldest_population,
         inmigrants_MRCA_times_in_model_time_oldest_population,
         time_origin_in_oldest_population_time,
         sorted_coal_times_in_oldest_population_time, topology,
         total_sample_size, pstream__);
}
};

template <typename T3__, typename T4__>
Eigen::Matrix<stan::promote_args_t<stan::value_type_t<T3__>,
stan::value_type_t<T4__>>, -1, 1>
get_torigins_of_inmigrant_populations(const int& pos, const int& order,
                                      const int& N,
                                      const T3__& sorted_previous_torigins_in_model_time_oldest_population_arg__,
                                      const T4__& torigins_in_model_time_oldest_population_arg__,
                                      const std::vector<int>& indexes_father_populations,
                                      std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<stan::value_type_t<T3__>,
          stan::value_type_t<T4__>>;
  int current_statement__ = 0;
  const auto& sorted_previous_torigins_in_model_time_oldest_population = to_ref(sorted_previous_torigins_in_model_time_oldest_population_arg__);
  const auto& torigins_in_model_time_oldest_population = to_ref(torigins_in_model_time_oldest_population_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    current_statement__ = 400;
    validate_non_negative_index("result", "order - 1", (order - 1));
    Eigen::Matrix<local_scalar_t__, -1, 1> result;
    result = Eigen::Matrix<local_scalar_t__, -1, 1>((order - 1));
    stan::math::fill(result, DUMMY_VAR__);
    
    int k;
    k = std::numeric_limits<int>::min();
    
    current_statement__ = 402;
    k = 1;
    current_statement__ = 418;
    for (int i = 1; i <= N; ++i) {
      current_statement__ = 416;
      if (logical_neq(i, pos)) {
        current_statement__ = 413;
        if (logical_eq(
              rvalue(indexes_father_populations,
                "indexes_father_populations", index_uni(i)), pos)) {
          current_statement__ = 410;
          if (logical_gt(k, (order - 1))) {
            current_statement__ = 408;
            std::stringstream errmsg_stream__;
            errmsg_stream__ << "index out bounds; found k=";
            errmsg_stream__ << k;
            throw std::domain_error(errmsg_stream__.str());
          }
          current_statement__ = 411;
          assign(result,
            rvalue(torigins_in_model_time_oldest_population,
              "torigins_in_model_time_oldest_population", index_uni(i)),
            "assigning variable result", index_uni(k));
        } else {
          current_statement__ = 405;
          if (logical_gt(k, (order - 1))) {
            current_statement__ = 403;
            std::stringstream errmsg_stream__;
            errmsg_stream__ << "index out bounds; found k=";
            errmsg_stream__ << k;
            throw std::domain_error(errmsg_stream__.str());
          }
          current_statement__ = 406;
          assign(result, -1.0, "assigning variable result", index_uni(k));
        }
        current_statement__ = 414;
        k = (k + 1);
      }
    }
    current_statement__ = 419;
    return result;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct get_torigins_of_inmigrant_populations_functor__ {
template <typename T3__, typename T4__>
Eigen::Matrix<stan::promote_args_t<stan::value_type_t<T3__>,
stan::value_type_t<T4__>>, -1, 1>
operator()(const int& pos, const int& order, const int& N,
           const T3__& sorted_previous_torigins_in_model_time_oldest_population,
           const T4__& torigins_in_model_time_oldest_population,
           const std::vector<int>& indexes_father_populations,
           std::ostream* pstream__)  const
{
return get_torigins_of_inmigrant_populations(pos, order, N,
         sorted_previous_torigins_in_model_time_oldest_population,
         torigins_in_model_time_oldest_population,
         indexes_father_populations, pstream__);
}
};

template <typename T0__, typename T1__, typename T2__>
int
population_of_coalescent_event(const T0__& x,
                               const T1__& grouped_coalescent_times_arg__,
                               const T2__& cumulative_number_coal_times_arg__,
                               std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__,
          stan::value_type_t<T1__>,
          stan::value_type_t<T2__>>;
  int current_statement__ = 0;
  const auto& grouped_coalescent_times = to_ref(grouped_coalescent_times_arg__);
  const auto& cumulative_number_coal_times = to_ref(cumulative_number_coal_times_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    int i;
    i = std::numeric_limits<int>::min();
    
    current_statement__ = 421;
    i = 1;
    int result;
    result = std::numeric_limits<int>::min();
    
    current_statement__ = 422;
    result = -1;
    int idx_population;
    idx_population = std::numeric_limits<int>::min();
    
    current_statement__ = 423;
    idx_population = 1;
    current_statement__ = 429;
    while ((primitive_value(logical_lte(i, rows(grouped_coalescent_times)))
           && primitive_value(
           logical_neq(
             rvalue(grouped_coalescent_times, "grouped_coalescent_times",
               index_uni(i)), x)))) {
      current_statement__ = 426;
      if (logical_eq(i,
            rvalue(cumulative_number_coal_times,
              "cumulative_number_coal_times", index_uni(idx_population)))) {
        current_statement__ = 424;
        idx_population = (idx_population + 1);
      }
      current_statement__ = 427;
      i = (i + 1);
    }
    current_statement__ = 432;
    if ((primitive_value(logical_eq(rows(grouped_coalescent_times), 0)) ||
        primitive_value(logical_gt(i, rows(grouped_coalescent_times))))) {
      current_statement__ = 431;
      result = -1;
    } else {
      current_statement__ = 430;
      result = idx_population;
    }
    current_statement__ = 433;
    return result;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct population_of_coalescent_event_functor__ {
template <typename T0__, typename T1__, typename T2__>
int
operator()(const T0__& x, const T1__& grouped_coalescent_times,
           const T2__& cumulative_number_coal_times, std::ostream* pstream__)  const
{
return population_of_coalescent_event(x, grouped_coalescent_times,
         cumulative_number_coal_times, pstream__);
}
};

template <typename T0__, typename T1__, typename T2__, typename T3__>
stan::promote_args_t<T0__, T1__, T2__,
T3__>
model_time_to_standard_time(const T0__& t, const T1__& torigin,
                            const T2__& delta, const T3__& K,
                            std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__, T1__, T2__, T3__>;
  int current_statement__ = 0;
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    local_scalar_t__ model_time;
    model_time = DUMMY_VAR__;
    
    current_statement__ = 435;
    model_time = 0.0;
    local_scalar_t__ a;
    a = DUMMY_VAR__;
    
    current_statement__ = 436;
    a = (stan::math::exp((delta * t)) - 1.0);
    local_scalar_t__ b;
    b = DUMMY_VAR__;
    
    current_statement__ = 437;
    b = (1.0 - stan::math::exp(((-1.0 * delta) * torigin)));
    local_scalar_t__ c;
    c = DUMMY_VAR__;
    
    current_statement__ = 438;
    c = (1.0 - stan::math::exp(((-1.0 * delta) * (torigin - t))));
    local_scalar_t__ d;
    d = DUMMY_VAR__;
    
    local_scalar_t__ numerator;
    numerator = DUMMY_VAR__;
    
    local_scalar_t__ denominator;
    denominator = DUMMY_VAR__;
    
    current_statement__ = 456;
    if (logical_eq(K, 0)) {
      current_statement__ = 453;
      model_time = ((a * b) / (delta * c));
      current_statement__ = 454;
      return ((2.0 * model_time) / delta);
    } else {
      current_statement__ = 442;
      c = stan::math::exp(((-1.0 * delta) * torigin));
      current_statement__ = 443;
      d = (1.0 - c);
      current_statement__ = 444;
      a = (((K / delta) * d) - c);
      current_statement__ = 445;
      b = (1.0 - ((K / delta) * d));
      current_statement__ = 446;
      numerator = ((a * stan::math::exp((delta * t))) + b);
      current_statement__ = 447;
      denominator = (1 - stan::math::exp((-delta * (torigin - t))));
      current_statement__ = 448;
      numerator = ((a * stan::math::exp((delta * t))) + b);
      current_statement__ = 449;
      denominator = (1 - stan::math::exp((-delta * (torigin - t))));
      current_statement__ = 450;
      model_time = ((2.0 / K) * stan::math::log((numerator / denominator)));
      current_statement__ = 451;
      return model_time;
    }
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct model_time_to_standard_time_functor__ {
template <typename T0__, typename T1__, typename T2__, typename T3__>
stan::promote_args_t<T0__, T1__, T2__,
T3__>
operator()(const T0__& t, const T1__& torigin, const T2__& delta,
           const T3__& K, std::ostream* pstream__)  const
{
return model_time_to_standard_time(t, torigin, delta, K, pstream__);
}
};

template <typename T0__, typename T1__, typename T2__, typename T3__>
stan::promote_args_t<T0__, T1__, T2__,
T3__>
standard_time_to_model_time(const T0__& V, const T1__& torigin,
                            const T2__& delta, const T3__& K,
                            std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__, T1__, T2__, T3__>;
  int current_statement__ = 0;
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    local_scalar_t__ a;
    a = DUMMY_VAR__;
    
    local_scalar_t__ b;
    b = DUMMY_VAR__;
    
    local_scalar_t__ c;
    c = DUMMY_VAR__;
    
    local_scalar_t__ d;
    d = DUMMY_VAR__;
    
    local_scalar_t__ e;
    e = DUMMY_VAR__;
    
    local_scalar_t__ x;
    x = DUMMY_VAR__;
    
    local_scalar_t__ firstTerm;
    firstTerm = DUMMY_VAR__;
    
    local_scalar_t__ secondTerm;
    secondTerm = DUMMY_VAR__;
    
    current_statement__ = 465;
    secondTerm = (1.0 / delta);
    local_scalar_t__ thirdTerm;
    thirdTerm = DUMMY_VAR__;
    
    local_scalar_t__ numerator;
    numerator = DUMMY_VAR__;
    
    local_scalar_t__ denominator;
    denominator = DUMMY_VAR__;
    
    local_scalar_t__ StandardTimeG;
    StandardTimeG = DUMMY_VAR__;
    
    current_statement__ = 507;
    if (logical_eq(K, 0)) {
      current_statement__ = 480;
      a = stan::math::exp(((-1.0 * delta) * torigin));
      current_statement__ = 481;
      b = (((1 - a) * (1 - a)) * (1.0 / a));
      current_statement__ = 482;
      c = (1 - a);
      current_statement__ = 483;
      d = ((1 - a) * (1.0 / a));
      current_statement__ = 484;
      e = (V + d);
      current_statement__ = 485;
      thirdTerm = stan::math::log((1 - (b / e)));
      current_statement__ = 486;
      thirdTerm = stan::math::log(
                    (1 -
                      ((((1 - a) * (1 - a)) * (1.0 / a)) /
                        ((V * delta) + ((1 - a) * (1.0 / a))))));
      current_statement__ = 487;
      thirdTerm = (stan::math::log(((1 + (delta * V)) - a)) -
                    stan::math::log((1 + (((delta * V) - 1) * a))));
      current_statement__ = 488;
      StandardTimeG = (secondTerm * thirdTerm);
      current_statement__ = 505;
      if ((primitive_value(logical_lte(((1 + (delta * V)) - a), 0)) ||
          primitive_value(logical_lte((1 + (((delta * V) - 1) * a)), 0)))) {
        current_statement__ = 489;
        StandardTimeG = 0.0;
        current_statement__ = 490;
        firstTerm = 0.0;
        current_statement__ = 491;
        secondTerm = 0.0;
        current_statement__ = 492;
        thirdTerm = 0.0;
        current_statement__ = 493;
        a = 0.0;
        current_statement__ = 494;
        b = 0.0;
        current_statement__ = 495;
        c = 0.0;
        current_statement__ = 496;
        d = 0.0;
        current_statement__ = 497;
        e = 0.0;
        current_statement__ = 498;
        a = (1 / delta);
        current_statement__ = 499;
        b = stan::math::log((1 + (delta * V)));
        current_statement__ = 500;
        firstTerm = (a * b);
        current_statement__ = 501;
        d = ((((V * V) * delta) *
               stan::math::exp(((-1.0 * delta) * torigin))) / (1 + V));
        current_statement__ = 502;
        secondTerm = d;
        current_statement__ = 503;
        StandardTimeG = (firstTerm - secondTerm);
      }
    } else {
      current_statement__ = 470;
      x = stan::math::exp(((K * V) / 2.0));
      current_statement__ = 471;
      c = stan::math::exp(((-1.0 * delta) * torigin));
      current_statement__ = 472;
      d = (1.0 - c);
      current_statement__ = 473;
      a = (((K / delta) * d) - c);
      current_statement__ = 474;
      b = (1.0 - ((K / delta) * d));
      current_statement__ = 475;
      numerator = (x - b);
      current_statement__ = 476;
      denominator = ((x * c) + a);
      current_statement__ = 477;
      thirdTerm = stan::math::log((numerator / denominator));
      current_statement__ = 478;
      StandardTimeG = (secondTerm * thirdTerm);
    }
    current_statement__ = 508;
    return StandardTimeG;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct standard_time_to_model_time_functor__ {
template <typename T0__, typename T1__, typename T2__, typename T3__>
stan::promote_args_t<T0__, T1__, T2__,
T3__>
operator()(const T0__& V, const T1__& torigin, const T2__& delta,
           const T3__& K, std::ostream* pstream__)  const
{
return standard_time_to_model_time(V, torigin, delta, K, pstream__);
}
};

template <typename T0__, typename T1__, typename T2__, typename T3__>
Eigen::Matrix<stan::promote_args_t<stan::value_type_t<T0__>, T1__, T2__,
T3__>, -1, 1>
from_model_time_to_Kingman_coalescent_time(const T0__& coal_times_model_time_arg__,
                                           const T1__& torigin,
                                           const T2__& delta, const T3__& K,
                                           std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<stan::value_type_t<T0__>,
          T1__,
          T2__,
          T3__>;
  int current_statement__ = 0;
  const auto& coal_times_model_time = to_ref(coal_times_model_time_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    current_statement__ = 510;
    validate_non_negative_index("result", "rows(coal_times_model_time)",
                                rows(coal_times_model_time));
    Eigen::Matrix<local_scalar_t__, -1, 1> result;
    result = Eigen::Matrix<local_scalar_t__, -1, 1>(rows(
                                                      coal_times_model_time));
    stan::math::fill(result, DUMMY_VAR__);
    
    current_statement__ = 516;
    for (int i = 1; i <= rows(coal_times_model_time); ++i) {
      current_statement__ = 514;
      if (logical_negation(
            is_nan(
              rvalue(coal_times_model_time, "coal_times_model_time",
                index_uni(i))))) {
        current_statement__ = 513;
        assign(result,
          model_time_to_standard_time(
            rvalue(coal_times_model_time, "coal_times_model_time",
              index_uni(i)), torigin, delta, K, pstream__),
          "assigning variable result", index_uni(i));
      } else {
        current_statement__ = 512;
        assign(result, -1, "assigning variable result", index_uni(i));
      }
    }
    current_statement__ = 517;
    return result;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct from_model_time_to_Kingman_coalescent_time_functor__ {
template <typename T0__, typename T1__, typename T2__, typename T3__>
Eigen::Matrix<stan::promote_args_t<stan::value_type_t<T0__>, T1__, T2__,
T3__>, -1, 1>
operator()(const T0__& coal_times_model_time, const T1__& torigin,
           const T2__& delta, const T3__& K, std::ostream* pstream__)  const
{
return from_model_time_to_Kingman_coalescent_time(coal_times_model_time,
         torigin, delta, K, pstream__);
}
};

template <typename T0__, typename T1__, typename T3__, typename T4__,
typename T5__>
stan::promote_args_t<T0__, T1__, T3__, T4__,
T5__>
log_lik_no_coalescent_between_times(const T0__& from, const T1__& to,
                                    const int& number_alive_ind,
                                    const T3__& torigin, const T4__& delta,
                                    const T5__& K, std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__,
          T1__,
          T3__,
          T4__,
          T5__>;
  int current_statement__ = 0;
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    local_scalar_t__ result;
    result = DUMMY_VAR__;
    
    current_statement__ = 519;
    result = 0.0;
    local_scalar_t__ cumulative_to;
    cumulative_to = DUMMY_VAR__;
    
    local_scalar_t__ cumulative_from;
    cumulative_from = DUMMY_VAR__;
    
    current_statement__ = 526;
    if (logical_gt(number_alive_ind, 1)) {
      current_statement__ = 522;
      cumulative_to = model_time_to_standard_time(to, torigin, delta,
                        K, pstream__);
      current_statement__ = 523;
      cumulative_from = model_time_to_standard_time(from, torigin, delta,
                          K, pstream__);
      current_statement__ = 524;
      result = (stan::math::log(choose(number_alive_ind, 2)) -
                 ((1.0 * choose(number_alive_ind, 2)) *
                   (cumulative_to - cumulative_from)));
    }
    current_statement__ = 527;
    return result;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct log_lik_no_coalescent_between_times_functor__ {
template <typename T0__, typename T1__, typename T3__, typename T4__,
typename T5__>
stan::promote_args_t<T0__, T1__, T3__, T4__,
T5__>
operator()(const T0__& from, const T1__& to, const int& number_alive_ind,
           const T3__& torigin, const T4__& delta, const T5__& K,
           std::ostream* pstream__)  const
{
return log_lik_no_coalescent_between_times(from, to, number_alive_ind,
         torigin, delta, K, pstream__);
}
};

template <typename T0__, typename RNG>
stan::promote_args_t<T0__>
conditionalDensityTOrigin_rng(const T0__& delta, const int& sample_size,
                              RNG& base_rng__, std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__>;
  int current_statement__ = 0;
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    local_scalar_t__ number_ancestors_population_when_sample_minus_one;
    number_ancestors_population_when_sample_minus_one = DUMMY_VAR__;
    
    local_scalar_t__ U;
    U = DUMMY_VAR__;
    
    local_scalar_t__ torigin;
    torigin = DUMMY_VAR__;
    
    local_scalar_t__ torigin_in_model_time_oldest_population;
    torigin_in_model_time_oldest_population = DUMMY_VAR__;
    
    local_scalar_t__ term;
    term = DUMMY_VAR__;
    
    current_statement__ = 534;
    number_ancestors_population_when_sample_minus_one = (2 * sample_size);
    current_statement__ = 535;
    U = gamma_rng((number_ancestors_population_when_sample_minus_one + 1), 1,
          base_rng__);
    current_statement__ = 536;
    term = (1.0 - (delta / (U + delta)));
    current_statement__ = 541;
    if (logical_lte(term, 0)) {
      current_statement__ = 539;
      std::stringstream errmsg_stream__;
      errmsg_stream__ << "term must not be negative; found term=";
      errmsg_stream__ << term;
      throw std::domain_error(errmsg_stream__.str());
    } else {
      current_statement__ = 537;
      torigin = (-(1.0 / delta) * stan::math::log(term));
    }
    current_statement__ = 542;
    return torigin;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct conditionalDensityTOrigin_rng_functor__ {
template <typename T0__, typename RNG>
stan::promote_args_t<T0__>
operator()(const T0__& delta, const int& sample_size, RNG& base_rng__,
           std::ostream* pstream__)  const
{
return conditionalDensityTOrigin_rng(delta, sample_size, base_rng__,
         pstream__);
}
};

template <typename T0__, typename T1__>
stan::promote_args_t<T0__,
T1__>
conditionalDensityTOrigin_pdf(const T0__& y, const T1__& delta,
                              const int& sample_size, std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__, T1__>;
  int current_statement__ = 0;
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    local_scalar_t__ result;
    result = DUMMY_VAR__;
    
    local_scalar_t__ term1;
    term1 = DUMMY_VAR__;
    
    local_scalar_t__ term2;
    term2 = DUMMY_VAR__;
    
    local_scalar_t__ term3;
    term3 = DUMMY_VAR__;
    
    local_scalar_t__ partial;
    partial = DUMMY_VAR__;
    
    current_statement__ = 549;
    term1 = stan::math::exp(((-1.0 * delta) * y));
    current_statement__ = 550;
    term2 = (delta * term1);
    current_statement__ = 551;
    term3 = (1.0 - term1);
    current_statement__ = 552;
    partial = ((delta * term2) / (term3 * term3));
    current_statement__ = 553;
    partial = (partial * stan::math::exp(((-1 * term2) / term3)));
    current_statement__ = 554;
    result = partial;
    current_statement__ = 555;
    return result;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct conditionalDensityTOrigin_pdf_functor__ {
template <typename T0__, typename T1__>
stan::promote_args_t<T0__,
T1__>
operator()(const T0__& y, const T1__& delta, const int& sample_size,
           std::ostream* pstream__)  const
{
return conditionalDensityTOrigin_pdf(y, delta, sample_size, pstream__);
}
};

template <bool propto__, typename T0__, typename T1__>
stan::promote_args_t<T0__,
T1__>
conditionalDensityTOrigin_lpdf(const T0__& y, const T1__& delta,
                               const int& sample_size,
                               std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__, T1__>;
  int current_statement__ = 0;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    local_scalar_t__ result;
    result = DUMMY_VAR__;
    
    current_statement__ = 558;
    result = stan::math::log(
               conditionalDensityTOrigin_pdf(y, delta,
                 sample_size, pstream__));
    current_statement__ = 559;
    return result;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct conditionalDensityTOrigin_lpdf_functor__ {
template <bool propto__, typename T0__, typename T1__>
stan::promote_args_t<T0__,
T1__>
operator()(const T0__& y, const T1__& delta, const int& sample_size,
           std::ostream* pstream__)  const
{
return conditionalDensityTOrigin_lpdf<propto__>(y, delta, sample_size,
         pstream__);
}
};

template <typename T0__, typename T1__>
stan::promote_args_t<T0__,
T1__>
conditionalDensityTOrigin_cdf(const T0__& y, const T1__& delta,
                              const int& sample_size, std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__, T1__>;
  int current_statement__ = 0;
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    local_scalar_t__ term1;
    term1 = DUMMY_VAR__;
    
    local_scalar_t__ term2;
    term2 = DUMMY_VAR__;
    
    local_scalar_t__ term3;
    term3 = DUMMY_VAR__;
    
    local_scalar_t__ result;
    result = DUMMY_VAR__;
    
    current_statement__ = 565;
    term1 = stan::math::exp(((-1.0 * delta) * y));
    current_statement__ = 566;
    term2 = (delta * term1);
    current_statement__ = 567;
    term3 = (1.0 - term1);
    current_statement__ = 568;
    result = stan::math::exp(((-1.0 * term2) / term3));
    current_statement__ = 569;
    return result;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct conditionalDensityTOrigin_cdf_functor__ {
template <typename T0__, typename T1__>
stan::promote_args_t<T0__,
T1__>
operator()(const T0__& y, const T1__& delta, const int& sample_size,
           std::ostream* pstream__)  const
{
return conditionalDensityTOrigin_cdf(y, delta, sample_size, pstream__);
}
};

template <typename T0__, typename T1__>
stan::promote_args_t<T0__,
T1__>
conditionalDensityTOrigin_lcdf(const T0__& y, const T1__& delta,
                               const int& sample_size,
                               std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__, T1__>;
  int current_statement__ = 0;
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    local_scalar_t__ term1;
    term1 = DUMMY_VAR__;
    
    local_scalar_t__ term2;
    term2 = DUMMY_VAR__;
    
    local_scalar_t__ term3;
    term3 = DUMMY_VAR__;
    
    local_scalar_t__ result;
    result = DUMMY_VAR__;
    
    current_statement__ = 575;
    term1 = stan::math::exp(((-1.0 * delta) * y));
    current_statement__ = 576;
    term2 = (delta * term1);
    current_statement__ = 577;
    term3 = (1.0 - term1);
    current_statement__ = 578;
    result = ((-1.0 * term2) / term3);
    current_statement__ = 579;
    return stan::math::log(result);
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct conditionalDensityTOrigin_lcdf_functor__ {
template <typename T0__, typename T1__>
stan::promote_args_t<T0__,
T1__>
operator()(const T0__& y, const T1__& delta, const int& sample_size,
           std::ostream* pstream__)  const
{
return conditionalDensityTOrigin_lcdf(y, delta, sample_size, pstream__);
}
};

template <typename T0__, typename T1__>
stan::promote_args_t<T0__,
T1__>
conditionalDensityTOrigin_lccdf(const T0__& y, const T1__& delta,
                                const int& sample_size,
                                std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__, T1__>;
  int current_statement__ = 0;
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    local_scalar_t__ term1;
    term1 = DUMMY_VAR__;
    
    local_scalar_t__ term2;
    term2 = DUMMY_VAR__;
    
    local_scalar_t__ term3;
    term3 = DUMMY_VAR__;
    
    local_scalar_t__ result;
    result = DUMMY_VAR__;
    
    current_statement__ = 585;
    result = conditionalDensityTOrigin_cdf(y, delta, sample_size, pstream__);
    current_statement__ = 586;
    return stan::math::log((1 - result));
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct conditionalDensityTOrigin_lccdf_functor__ {
template <typename T0__, typename T1__>
stan::promote_args_t<T0__,
T1__>
operator()(const T0__& y, const T1__& delta, const int& sample_size,
           std::ostream* pstream__)  const
{
return conditionalDensityTOrigin_lccdf(y, delta, sample_size, pstream__);
}
};

template <typename T0__, typename T1__, typename T2__, typename T3__>
stan::promote_args_t<T0__, T1__, T2__,
T3__>
log_h(const T0__& t, const T1__& torigin, const T2__& delta, const T3__& K,
      std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__, T1__, T2__, T3__>;
  int current_statement__ = 0;
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    local_scalar_t__ a;
    a = DUMMY_VAR__;
    
    current_statement__ = 588;
    a = (1.0 - stan::math::exp(((-1.0 * delta) * (torigin - t))));
    local_scalar_t__ first_term;
    first_term = DUMMY_VAR__;
    
    current_statement__ = 589;
    first_term = (2.0 * stan::math::log(a));
    local_scalar_t__ second_term;
    second_term = DUMMY_VAR__;
    
    current_statement__ = 590;
    second_term = ((-1.0 * delta) * t);
    local_scalar_t__ third_term;
    third_term = DUMMY_VAR__;
    
    current_statement__ = 591;
    third_term = stan::math::exp((delta * t));
    local_scalar_t__ above_term;
    above_term = DUMMY_VAR__;
    
    current_statement__ = 592;
    above_term = (first_term + second_term);
    local_scalar_t__ b;
    b = DUMMY_VAR__;
    
    current_statement__ = 593;
    b = (1.0 - stan::math::exp(((-1.0 * delta) * torigin)));
    local_scalar_t__ below_term;
    below_term = DUMMY_VAR__;
    
    current_statement__ = 594;
    below_term = (2.0 * stan::math::log(b));
    local_scalar_t__ extra_term;
    extra_term = DUMMY_VAR__;
    
    current_statement__ = 595;
    extra_term = stan::math::log(
                   (1 + ((((K / delta) * b) * (third_term - 1.0)) / a)));
    local_scalar_t__ logH;
    logH = DUMMY_VAR__;
    
    current_statement__ = 599;
    if (logical_eq(K, 0)) {
      current_statement__ = 598;
      logH = (above_term - below_term);
    } else {
      current_statement__ = 597;
      logH = ((above_term - below_term) + extra_term);
    }
    current_statement__ = 600;
    return logH;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct log_h_functor__ {
template <typename T0__, typename T1__, typename T2__, typename T3__>
stan::promote_args_t<T0__, T1__, T2__,
T3__>
operator()(const T0__& t, const T1__& torigin, const T2__& delta,
           const T3__& K, std::ostream* pstream__)  const
{
return log_h(t, torigin, delta, K, pstream__);
}
};

template <typename T2__, typename T3__, typename T4__, typename T5__,
typename T6__, typename T7__, typename T8__, typename T9__, typename T10__,
typename T11__, typename T13__>
stan::promote_args_t<T2__, T3__, T4__, T5__, T6__, stan::promote_args_t<T7__,
stan::value_type_t<T8__>, stan::value_type_t<T9__>,
stan::value_type_t<T10__>,
T11__, stan::promote_args_t<stan::value_type_t<T13__>>>>
log_prob_father_population(const int& order,
                           const int& index_oldest_population,
                           const T2__& father_torigin,
                           const T3__& father_delta, const T4__& x_father,
                           const T5__& torigin, const T6__& delta,
                           const T7__& x, const T8__& deltas_arg__,
                           const T9__& torigins_in_model_time_arg__,
                           const T10__& torigins_in_model_time_oldest_population_par_arg__,
                           const T11__& K, const int& N,
                           const T13__& pop_sizes_proportion_arg__,
                           std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T2__,
          T3__,
          T4__,
          T5__,
          T6__, stan::promote_args_t<T7__,
          stan::value_type_t<T8__>,
          stan::value_type_t<T9__>,
          stan::value_type_t<T10__>,
          T11__, stan::promote_args_t<stan::value_type_t<T13__>>>>;
  int current_statement__ = 0;
  const auto& deltas = to_ref(deltas_arg__);
  const auto& torigins_in_model_time = to_ref(torigins_in_model_time_arg__);
  const auto& torigins_in_model_time_oldest_population_par = to_ref(torigins_in_model_time_oldest_population_par_arg__);
  const auto& pop_sizes_proportion = to_ref(pop_sizes_proportion_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    local_scalar_t__ numerator;
    numerator = DUMMY_VAR__;
    
    local_scalar_t__ denominator;
    denominator = DUMMY_VAR__;
    
    local_scalar_t__ logh;
    logh = DUMMY_VAR__;
    
    int order_father_pop;
    order_father_pop = std::numeric_limits<int>::min();
    
    local_scalar_t__ temp;
    temp = DUMMY_VAR__;
    
    local_scalar_t__ result;
    result = DUMMY_VAR__;
    
    current_statement__ = 608;
    validate_non_negative_index("ordered_torigins_in_model_time_oldest_population",
                                "N", N);
    Eigen::Matrix<local_scalar_t__, -1, 1> ordered_torigins_in_model_time_oldest_population;
    ordered_torigins_in_model_time_oldest_population = Eigen::Matrix<local_scalar_t__, -1, 1>(N);
    stan::math::fill(ordered_torigins_in_model_time_oldest_population, DUMMY_VAR__);
    
    
    current_statement__ = 609;
    assign(ordered_torigins_in_model_time_oldest_population,
      sort_asc(torigins_in_model_time_oldest_population_par),
      "assigning variable ordered_torigins_in_model_time_oldest_population");
    current_statement__ = 610;
    order_father_pop = find(
                         ((x_father /
                            rvalue(pop_sizes_proportion,
                              "pop_sizes_proportion",
                              index_uni(index_oldest_population))) *
                           father_torigin),
                         ordered_torigins_in_model_time_oldest_population, pstream__);
    current_statement__ = 619;
    if ((primitive_value((primitive_value(logical_neq(order, -1)) &&
        primitive_value(logical_neq(order_father_pop, -1)))) &&
        primitive_value(logical_lt(order, order_father_pop)))) {
      current_statement__ = 617;
      for (int i = (order + 1); i <= N; ++i) {
        int pos;
        pos = std::numeric_limits<int>::min();
        
        current_statement__ = 611;
        pos = find(
                ((rvalue(pop_sizes_proportion, "pop_sizes_proportion",
                    index_uni(index_oldest_population)) /
                   rvalue(pop_sizes_proportion, "pop_sizes_proportion",
                     index_uni(i))) *
                  rvalue(torigins_in_model_time_oldest_population_par,
                    "torigins_in_model_time_oldest_population_par",
                    index_uni(i))), torigins_in_model_time, pstream__);
        current_statement__ = 612;
        denominator = (denominator +
                        stan::math::log(
                          rvalue(pop_sizes_proportion,
                            "pop_sizes_proportion", index_uni(pos))));
        current_statement__ = 613;
        temp = ((torigin * x) /
                 rvalue(pop_sizes_proportion, "pop_sizes_proportion",
                   index_uni(pos)));
        current_statement__ = 614;
        logh = log_h(temp,
                 rvalue(torigins_in_model_time, "torigins_in_model_time",
                   index_uni(pos)), rvalue(deltas, "deltas", index_uni(pos)),
                 K, pstream__);
        current_statement__ = 615;
        denominator = (denominator + temp);
      }
    }
    current_statement__ = 620;
    numerator = (numerator + stan::math::log(x_father));
    current_statement__ = 621;
    temp = ((torigin * x) / x_father);
    current_statement__ = 622;
    logh = log_h(temp, father_torigin, father_delta, K, pstream__);
    current_statement__ = 623;
    numerator = (numerator + temp);
    current_statement__ = 624;
    result = (numerator / denominator);
    current_statement__ = 625;
    return result;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct log_prob_father_population_functor__ {
template <typename T2__, typename T3__, typename T4__, typename T5__,
typename T6__, typename T7__, typename T8__, typename T9__, typename T10__,
typename T11__, typename T13__>
stan::promote_args_t<T2__, T3__, T4__, T5__, T6__, stan::promote_args_t<T7__,
stan::value_type_t<T8__>, stan::value_type_t<T9__>,
stan::value_type_t<T10__>,
T11__, stan::promote_args_t<stan::value_type_t<T13__>>>>
operator()(const int& order, const int& index_oldest_population,
           const T2__& father_torigin, const T3__& father_delta,
           const T4__& x_father, const T5__& torigin, const T6__& delta,
           const T7__& x, const T8__& deltas,
           const T9__& torigins_in_model_time,
           const T10__& torigins_in_model_time_oldest_population_par,
           const T11__& K, const int& N, const T13__& pop_sizes_proportion,
           std::ostream* pstream__)  const
{
return log_prob_father_population(order, index_oldest_population,
         father_torigin, father_delta, x_father, torigin, delta, x, deltas,
         torigins_in_model_time,
         torigins_in_model_time_oldest_population_par, K, N,
         pop_sizes_proportion, pstream__);
}
};

template <typename T2__, typename T3__, typename T4__, typename T6__,
typename T7__>
stan::promote_args_t<T2__, T3__, T4__, stan::value_type_t<T6__>,
stan::value_type_t<T7__>>
log_likelihood_subpopulation(const int& order, const int& sample_size,
                             const T2__& K, const T3__& delta,
                             const T4__& torigin,
                             const std::vector<int>& positions,
                             const T6__& sorted_coal_time_with_previous_torigins_in_model_time_arg__,
                             const T7__& previous_torigins_in_model_time_arg__,
                             std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T2__,
          T3__,
          T4__,
          stan::value_type_t<T6__>,
          stan::value_type_t<T7__>>;
  int current_statement__ = 0;
  const auto& sorted_coal_time_with_previous_torigins_in_model_time = to_ref(sorted_coal_time_with_previous_torigins_in_model_time_arg__);
  const auto& previous_torigins_in_model_time = to_ref(previous_torigins_in_model_time_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    local_scalar_t__ log_lik;
    log_lik = DUMMY_VAR__;
    
    current_statement__ = 627;
    log_lik = 0.0;
    local_scalar_t__ current_time;
    current_time = DUMMY_VAR__;
    
    int alive_cells;
    alive_cells = std::numeric_limits<int>::min();
    
    local_scalar_t__ term_only_after_first_coal_event;
    term_only_after_first_coal_event = DUMMY_VAR__;
    
    local_scalar_t__ zero;
    zero = DUMMY_VAR__;
    
    current_statement__ = 631;
    zero = 0.0;
    int currentCoalescentEventInThisEpoch;
    currentCoalescentEventInThisEpoch = std::numeric_limits<int>::min();
    
    current_statement__ = 632;
    currentCoalescentEventInThisEpoch = -1;
    current_statement__ = 633;
    validate_non_negative_index("padded_coalescent_times",
                                "rows(sorted_coal_time_with_previous_torigins_in_model_time) + 1",
                                (rows(
                                   sorted_coal_time_with_previous_torigins_in_model_time)
                                  + 1));
    Eigen::Matrix<local_scalar_t__, -1, 1> padded_coalescent_times;
    padded_coalescent_times = Eigen::Matrix<local_scalar_t__, -1, 1>(
      (rows(sorted_coal_time_with_previous_torigins_in_model_time) + 1));
    stan::math::fill(padded_coalescent_times, DUMMY_VAR__);
    
    current_statement__ = 634;
    assign(padded_coalescent_times,
      append_row(zero, sorted_coal_time_with_previous_torigins_in_model_time),
      "assigning variable padded_coalescent_times");
    local_scalar_t__ last_event_time_before_migration;
    last_event_time_before_migration = DUMMY_VAR__;
    
    current_statement__ = 635;
    last_event_time_before_migration = 0.0;
    current_statement__ = 636;
    alive_cells = sample_size;
    current_statement__ = 665;
    for (int j = 1; j <= rows(padded_coalescent_times); ++j) {
      current_statement__ = 637;
      current_time = rvalue(padded_coalescent_times,
                       "padded_coalescent_times", index_uni(j));
      current_statement__ = 639;
      if (logical_lt(stan::math::fabs((current_time - torigin)), 0.00000001)) {
        current_statement__ = 638;
        return log_lik;
      }
      current_statement__ = 641;
      if (logical_lt(current_time, 0.0)) {
        continue;
      }
      current_statement__ = 663;
      if ((primitive_value((primitive_value(logical_gt(order, 1)) &&
          primitive_value(
          logical_gt(rows(previous_torigins_in_model_time), 0)))) &&
          primitive_value(
          logical_neq(
            find(current_time, previous_torigins_in_model_time, pstream__),
            -1)))) {
        current_statement__ = 661;
        if (logical_gt(alive_cells, 1)) {
          current_statement__ = 656;
          log_lik = (log_lik +
                      log_lik_no_coalescent_between_times(
                        last_event_time_before_migration, current_time,
                        alive_cells, torigin, delta, K, pstream__));
          current_statement__ = 657;
          alive_cells = (alive_cells + 1);
          current_statement__ = 658;
          currentCoalescentEventInThisEpoch = 1;
          current_statement__ = 659;
          last_event_time_before_migration = current_time;
        }
      } else {
        current_statement__ = 654;
        if (logical_gt(alive_cells, 1)) {
          current_statement__ = 642;
          log_lik = (log_lik +
                      stan::math::log(
                        ((alive_cells * (alive_cells - 1)) / 2.0)));
          current_statement__ = 643;
          log_lik = ((log_lik + stan::math::log(2.0)) -
                      log_h(current_time, torigin, delta, K, pstream__));
          current_statement__ = 648;
          if (logical_lt(stan::math::fabs(current_time), 0.00001)) {
            current_statement__ = 646;
            term_only_after_first_coal_event = 0.0;
          } else {
            current_statement__ = 644;
            term_only_after_first_coal_event = model_time_to_standard_time(
                                                 last_event_time_before_migration,
                                                 torigin, delta,
                                                 K, pstream__);
          }
          current_statement__ = 649;
          log_lik = (log_lik -
                      (((alive_cells * (alive_cells - 1)) / 2.0) *
                        (model_time_to_standard_time(current_time, torigin,
                           delta, K, pstream__) -
                          term_only_after_first_coal_event)));
          current_statement__ = 650;
          alive_cells = (alive_cells - 1);
          current_statement__ = 651;
          currentCoalescentEventInThisEpoch = (currentCoalescentEventInThisEpoch
                                                + 1);
          current_statement__ = 652;
          last_event_time_before_migration = current_time;
        }
      }
    }
    current_statement__ = 666;
    return log_lik;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct log_likelihood_subpopulation_functor__ {
template <typename T2__, typename T3__, typename T4__, typename T6__,
typename T7__>
stan::promote_args_t<T2__, T3__, T4__, stan::value_type_t<T6__>,
stan::value_type_t<T7__>>
operator()(const int& order, const int& sample_size, const T2__& K,
           const T3__& delta, const T4__& torigin,
           const std::vector<int>& positions,
           const T6__& sorted_coal_time_with_previous_torigins_in_model_time,
           const T7__& previous_torigins_in_model_time,
           std::ostream* pstream__)  const
{
return log_likelihood_subpopulation(order, sample_size, K, delta, torigin,
         positions, sorted_coal_time_with_previous_torigins_in_model_time,
         previous_torigins_in_model_time, pstream__);
}
};

template <typename T0__, typename T1__, typename T2__, typename T3__,
typename T5__, typename RNG>
Eigen::Matrix<stan::promote_args_t<stan::value_type_t<T0__>, stan::value_type_t<T1__>,
stan::value_type_t<T2__>, stan::value_type_t<T3__>,
T5__>, -1, 1>
structured_coalescent_rng(const T0__& deltas_arg__,
                          const T1__& torigins_in_model_time_arg__,
                          const T2__& torigins_in_model_time_oldest_population_arg__,
                          const T3__& pop_sizes_proportion_arg__,
                          const std::vector<int>& sample_sizes,
                          const T5__& K, RNG& base_rng__,
                          std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<stan::value_type_t<T0__>,
          stan::value_type_t<T1__>,
          stan::value_type_t<T2__>,
          stan::value_type_t<T3__>,
          T5__>;
  int current_statement__ = 0;
  const auto& deltas = to_ref(deltas_arg__);
  const auto& torigins_in_model_time = to_ref(torigins_in_model_time_arg__);
  const auto& torigins_in_model_time_oldest_population = to_ref(torigins_in_model_time_oldest_population_arg__);
  const auto& pop_sizes_proportion = to_ref(pop_sizes_proportion_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    int total_sample_size;
    total_sample_size = std::numeric_limits<int>::min();
    
    current_statement__ = 668;
    total_sample_size = sum(sample_sizes);
    int N;
    N = std::numeric_limits<int>::min();
    
    current_statement__ = 669;
    N = rows(deltas);
    int current_pos;
    current_pos = std::numeric_limits<int>::min();
    
    current_statement__ = 670;
    current_pos = 1;
    int current_pos_migrations;
    current_pos_migrations = std::numeric_limits<int>::min();
    
    current_statement__ = 671;
    current_pos_migrations = 1;
    int number_inmigrants;
    number_inmigrants = std::numeric_limits<int>::min();
    
    current_statement__ = 672;
    number_inmigrants = 1;
    int index_next_immigrants;
    index_next_immigrants = std::numeric_limits<int>::min();
    
    current_statement__ = 673;
    index_next_immigrants = 1;
    int global_index_next_immigrants;
    global_index_next_immigrants = std::numeric_limits<int>::min();
    
    current_statement__ = 674;
    global_index_next_immigrants = 1;
    int indexes_inmigrants;
    indexes_inmigrants = std::numeric_limits<int>::min();
    
    current_statement__ = 675;
    indexes_inmigrants = 1;
    current_statement__ = 676;
    validate_non_negative_index("coal_event_times_model_time_oldest_population",
                                "total_sample_size - 1",
                                (total_sample_size - 1));
    Eigen::Matrix<local_scalar_t__, -1, 1> coal_event_times_model_time_oldest_population;
    coal_event_times_model_time_oldest_population = Eigen::Matrix<local_scalar_t__, -1, 1>(
      (total_sample_size - 1));
    stan::math::fill(coal_event_times_model_time_oldest_population, DUMMY_VAR__);
    
    
    current_statement__ = 678;
    assign(coal_event_times_model_time_oldest_population,
      rep_vector(
        (rvalue(torigins_in_model_time_oldest_population,
           "torigins_in_model_time_oldest_population", index_uni(1)) + 1),
        (total_sample_size - 1)),
      "assigning variable coal_event_times_model_time_oldest_population");
    current_statement__ = 707;
    for (int i = 1; i <= num_elements(sample_sizes); ++i) {
      int sample_size;
      sample_size = std::numeric_limits<int>::min();
      
      current_statement__ = 679;
      sample_size = rvalue(sample_sizes, "sample_sizes", index_uni(i));
      int current_sample_size;
      current_sample_size = std::numeric_limits<int>::min();
      
      current_statement__ = 680;
      current_sample_size = sample_size;
      int index_father_population;
      index_father_population = std::numeric_limits<int>::min();
      
      current_statement__ = 681;
      index_father_population = 0;
      current_statement__ = 682;
      current_pos = 1;
      current_statement__ = 705;
      if (logical_eq(i, 1)) {
        current_statement__ = 683;
        validate_non_negative_index("number_ancestors_sample",
                                    "sample_size - 1", (sample_size - 1));
        Eigen::Matrix<local_scalar_t__, -1, 1> number_ancestors_sample;
        number_ancestors_sample = Eigen::Matrix<local_scalar_t__, -1, 1>(
          (sample_size - 1));
        stan::math::fill(number_ancestors_sample, DUMMY_VAR__);
        
        current_statement__ = 685;
        validate_non_negative_index("rate_exp", "sample_size - 1",
                                    (sample_size - 1));
        Eigen::Matrix<local_scalar_t__, -1, 1> rate_exp;
        rate_exp = Eigen::Matrix<local_scalar_t__, -1, 1>((sample_size - 1));
        stan::math::fill(rate_exp, DUMMY_VAR__);
        
        current_statement__ = 687;
        validate_non_negative_index("times", "sample_size - 1",
                                    (sample_size - 1));
        Eigen::Matrix<local_scalar_t__, -1, 1> times;
        times = Eigen::Matrix<local_scalar_t__, -1, 1>((sample_size - 1));
        stan::math::fill(times, DUMMY_VAR__);
        
        current_statement__ = 689;
        validate_non_negative_index("cum_sum_times", "sample_size - 1",
                                    (sample_size - 1));
        Eigen::Matrix<local_scalar_t__, -1, 1> cum_sum_times;
        cum_sum_times = Eigen::Matrix<local_scalar_t__, -1, 1>((sample_size -
                                                                 1));
        stan::math::fill(cum_sum_times, DUMMY_VAR__);
        
        current_statement__ = 695;
        for (int j = 2; j <= sample_size; ++j) {
          current_statement__ = 691;
          assign(number_ancestors_sample, j,
            "assigning variable number_ancestors_sample", index_uni((j - 1)));
          current_statement__ = 692;
          assign(rate_exp, ((0.5 * j) * (j - 1)),
            "assigning variable rate_exp", index_uni((j - 1)));
          current_statement__ = 693;
          assign(times,
            exponential_rng(rvalue(rate_exp, "rate_exp", index_uni((j - 1))),
              base_rng__), "assigning variable times", index_uni((j - 1)));
        }
        current_statement__ = 699;
        for (int j = 1; j <= (sample_size - 1); ++j) {
          current_statement__ = 696;
          assign(cum_sum_times,
            sum(rvalue(times, "times", index_min_max(j, (sample_size - 1)))),
            "assigning variable cum_sum_times", index_uni(j));
          current_statement__ = 697;
          assign(coal_event_times_model_time_oldest_population,
            standard_time_to_model_time(
              rvalue(cum_sum_times, "cum_sum_times", index_uni(j)),
              rvalue(torigins_in_model_time, "torigins_in_model_time",
                index_uni(i)), rvalue(deltas, "deltas", index_uni(i)),
              K, pstream__),
            "assigning variable coal_event_times_model_time_oldest_population",
            index_uni((sample_size - j)));
        }
        current_statement__ = 700;
        assign(coal_event_times_model_time_oldest_population,
          multiply(
            (rvalue(pop_sizes_proportion, "pop_sizes_proportion",
               index_uni(i)) /
              rvalue(pop_sizes_proportion, "pop_sizes_proportion",
                index_uni(N))),
            stan::model::deep_copy(
              rvalue(coal_event_times_model_time_oldest_population,
                "coal_event_times_model_time_oldest_population",
                index_min_max(1, (sample_size - 1))))),
          "assigning variable coal_event_times_model_time_oldest_population",
          index_min_max(1, (sample_size - 1)));
        current_statement__ = 701;
        assign(coal_event_times_model_time_oldest_population,
          rvalue(torigins_in_model_time_oldest_population,
            "torigins_in_model_time_oldest_population", index_uni(i)),
          "assigning variable coal_event_times_model_time_oldest_population",
          index_uni(sample_size));
        current_statement__ = 702;
        current_pos = (sample_size + 1);
        current_statement__ = 703;
        index_father_population = 2;
      }
    }
    current_statement__ = 708;
    return coal_event_times_model_time_oldest_population;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct structured_coalescent_rng_functor__ {
template <typename T0__, typename T1__, typename T2__, typename T3__,
typename T5__, typename RNG>
Eigen::Matrix<stan::promote_args_t<stan::value_type_t<T0__>, stan::value_type_t<T1__>,
stan::value_type_t<T2__>, stan::value_type_t<T3__>,
T5__>, -1, 1>
operator()(const T0__& deltas, const T1__& torigins_in_model_time,
           const T2__& torigins_in_model_time_oldest_population,
           const T3__& pop_sizes_proportion,
           const std::vector<int>& sample_sizes, const T5__& K,
           RNG& base_rng__, std::ostream* pstream__)  const
{
return structured_coalescent_rng(deltas, torigins_in_model_time,
         torigins_in_model_time_oldest_population, pop_sizes_proportion,
         sample_sizes, K, base_rng__, pstream__);
}
};

template <bool propto__, typename T0__, typename T1__, typename T2__,
typename T3__, typename T4__, typename T6__>
stan::promote_args_t<stan::value_type_t<T0__>, stan::value_type_t<T1__>,
stan::value_type_t<T2__>, stan::value_type_t<T3__>,
stan::value_type_t<T4__>, stan::promote_args_t<T6__>>
structured_coalescent_lpdf(const T0__& coal_event_times_model_time_oldest_population_arg__,
                           const T1__& deltas_arg__,
                           const T2__& torigins_in_model_time_arg__,
                           const T3__& sorted_torigins_in_model_time_oldest_population_arg__,
                           const T4__& pop_sizes_proportion_arg__,
                           const std::vector<int>& sample_sizes,
                           const T6__& K, std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<stan::value_type_t<T0__>,
          stan::value_type_t<T1__>,
          stan::value_type_t<T2__>,
          stan::value_type_t<T3__>,
          stan::value_type_t<T4__>, stan::promote_args_t<T6__>>;
  int current_statement__ = 0;
  const auto& coal_event_times_model_time_oldest_population = to_ref(coal_event_times_model_time_oldest_population_arg__);
  const auto& deltas = to_ref(deltas_arg__);
  const auto& torigins_in_model_time = to_ref(torigins_in_model_time_arg__);
  const auto& sorted_torigins_in_model_time_oldest_population = to_ref(sorted_torigins_in_model_time_oldest_population_arg__);
  const auto& pop_sizes_proportion = to_ref(pop_sizes_proportion_arg__);
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    local_scalar_t__ log_lik;
    log_lik = DUMMY_VAR__;
    
    current_statement__ = 710;
    log_lik = 0.0;
    int cum_sample_size;
    cum_sample_size = std::numeric_limits<int>::min();
    
    current_statement__ = 711;
    cum_sample_size = 0;
    int cum_number_coal;
    cum_number_coal = std::numeric_limits<int>::min();
    
    current_statement__ = 712;
    cum_number_coal = 0;
    int first_position_coal_event;
    first_position_coal_event = std::numeric_limits<int>::min();
    
    current_statement__ = 713;
    first_position_coal_event = 1;
    int N;
    N = std::numeric_limits<int>::min();
    
    current_statement__ = 714;
    N = rows(deltas);
    int pos_next_torigin;
    pos_next_torigin = std::numeric_limits<int>::min();
    
    current_statement__ = 715;
    pos_next_torigin = 0;
    current_statement__ = 716;
    cum_sample_size = 0;
    current_statement__ = 728;
    for (int i = 1; i <= num_elements(sample_sizes); ++i) {
      int order;
      order = std::numeric_limits<int>::min();
      
      current_statement__ = 717;
      order = i;
      int number_events;
      number_events = std::numeric_limits<int>::min();
      
      current_statement__ = 718;
      number_events = (rvalue(sample_sizes, "sample_sizes", index_uni(i)) -
                        1);
      current_statement__ = 719;
      validate_non_negative_index("positions", "sample_sizes[i]",
                                  rvalue(sample_sizes, "sample_sizes",
                                    index_uni(i)));
      std::vector<int> positions;
      positions = std::vector<int>(rvalue(sample_sizes, "sample_sizes",
                                     index_uni(i)), std::numeric_limits<int>::min());
      
      
      current_statement__ = 720;
      assign(positions,
        rep_array(1, rvalue(sample_sizes, "sample_sizes", index_uni(i))),
        "assigning variable positions");
      current_statement__ = 721;
      validate_non_negative_index("current_coal_event_times_with_immigrants",
                                  "number_events", number_events);
      Eigen::Matrix<local_scalar_t__, -1, 1> current_coal_event_times_with_immigrants;
      current_coal_event_times_with_immigrants = Eigen::Matrix<local_scalar_t__, -1, 1>(number_events);
      stan::math::fill(current_coal_event_times_with_immigrants, DUMMY_VAR__);
      
      current_statement__ = 726;
      if (logical_eq(i, 1)) {
        current_statement__ = 723;
        assign(current_coal_event_times_with_immigrants,
          rvalue(coal_event_times_model_time_oldest_population,
            "coal_event_times_model_time_oldest_population",
            index_min_max(1, number_events)),
          "assigning variable current_coal_event_times_with_immigrants",
          index_min_max(1, number_events));
        current_statement__ = 724;
        log_lik = (log_lik +
                    log_likelihood_subpopulation(i,
                      rvalue(sample_sizes, "sample_sizes", index_uni(i)), K,
                      rvalue(deltas, "deltas", index_uni(i)),
                      rvalue(torigins_in_model_time,
                        "torigins_in_model_time", index_uni(i)), positions,
                      current_coal_event_times_with_immigrants,
                      rep_vector(1, 1), pstream__));
      }
    }
    current_statement__ = 729;
    return log_lik;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct structured_coalescent_lpdf_functor__ {
template <bool propto__, typename T0__, typename T1__, typename T2__,
typename T3__, typename T4__, typename T6__>
stan::promote_args_t<stan::value_type_t<T0__>, stan::value_type_t<T1__>,
stan::value_type_t<T2__>, stan::value_type_t<T3__>,
stan::value_type_t<T4__>, stan::promote_args_t<T6__>>
operator()(const T0__& coal_event_times_model_time_oldest_population,
           const T1__& deltas, const T2__& torigins_in_model_time,
           const T3__& sorted_torigins_in_model_time_oldest_population,
           const T4__& pop_sizes_proportion,
           const std::vector<int>& sample_sizes, const T6__& K,
           std::ostream* pstream__)  const
{
return structured_coalescent_lpdf<propto__>(
         coal_event_times_model_time_oldest_population, deltas,
         torigins_in_model_time,
         sorted_torigins_in_model_time_oldest_population,
         pop_sizes_proportion, sample_sizes, K, pstream__);
}
};

template <typename T0__, typename T1__>
int
get_sample_size(const T0__& sorted_coal_times_arg__,
                const T1__& MRCA_time_in_model_time_oldest_population,
                const std::vector<int>& number_of_tips_below,
                const std::vector<std::vector<int>>& topology,
                std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<stan::value_type_t<T0__>,
          T1__>;
  int current_statement__ = 0;
  const auto& sorted_coal_times = to_ref(sorted_coal_times_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    int result;
    result = std::numeric_limits<int>::min();
    
    current_statement__ = 731;
    result = 0;
    int pos;
    pos = std::numeric_limits<int>::min();
    
    current_statement__ = 732;
    pos = find(MRCA_time_in_model_time_oldest_population,
            sorted_coal_times, pstream__);
    current_statement__ = 734;
    if ((primitive_value(logical_gt(pos, 0)) && primitive_value(
        logical_lte(pos, num_elements(number_of_tips_below))))) {
      current_statement__ = 733;
      result = rvalue(number_of_tips_below, "number_of_tips_below",
                 index_uni(rvalue(topology, "topology",
                             index_uni(pos), index_uni(1))));
    }
    current_statement__ = 735;
    return result;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct get_sample_size_functor__ {
template <typename T0__, typename T1__>
int
operator()(const T0__& sorted_coal_times,
           const T1__& MRCA_time_in_model_time_oldest_population,
           const std::vector<int>& number_of_tips_below,
           const std::vector<std::vector<int>>& topology,
           std::ostream* pstream__)  const
{
return get_sample_size(sorted_coal_times,
         MRCA_time_in_model_time_oldest_population, number_of_tips_below,
         topology, pstream__);
}
};

template <typename T2__, typename T3__, typename T4__, typename T5__,
typename T6__, typename T8__>
Eigen::Matrix<stan::promote_args_t<stan::value_type_t<T2__>, stan::value_type_t<T3__>,
stan::value_type_t<T4__>, T5__,
T6__, stan::promote_args_t<stan::value_type_t<T8__>>>, -1, 1>
compute_probab_vector(const int& total_sample_size,
                      const int& number_candidate_MRCAs,
                      const T2__& candidate_time_MRCAs_in_model_time_arg__,
                      const T3__& candidate_time_MRCAs_in_model_time_oldest_population_arg__,
                      const T4__& Q_coeff_hypoexponential_arg__,
                      const T5__& delta, const T6__& torigin_in_model_time,
                      const std::vector<int>& number_of_tips_below,
                      const T8__& sorted_coal_times_in_oldest_population_time_arg__,
                      const std::vector<std::vector<int>>& topology,
                      std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<stan::value_type_t<T2__>,
          stan::value_type_t<T3__>,
          stan::value_type_t<T4__>,
          T5__,
          T6__, stan::promote_args_t<stan::value_type_t<T8__>>>;
  int current_statement__ = 0;
  const auto& candidate_time_MRCAs_in_model_time = to_ref(candidate_time_MRCAs_in_model_time_arg__);
  const auto& candidate_time_MRCAs_in_model_time_oldest_population = to_ref(candidate_time_MRCAs_in_model_time_oldest_population_arg__);
  const auto& Q_coeff_hypoexponential = to_ref(Q_coeff_hypoexponential_arg__);
  const auto& sorted_coal_times_in_oldest_population_time = to_ref(sorted_coal_times_in_oldest_population_time_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    local_scalar_t__ cum_sum;
    cum_sum = DUMMY_VAR__;
    
    current_statement__ = 737;
    cum_sum = 0.0;
    current_statement__ = 738;
    validate_non_negative_index("result", "number_candidate_MRCAs",
                                number_candidate_MRCAs);
    Eigen::Matrix<local_scalar_t__, -1, 1> result;
    result = Eigen::Matrix<local_scalar_t__, -1, 1>(number_candidate_MRCAs);
    stan::math::fill(result, DUMMY_VAR__);
    
    current_statement__ = 751;
    for (int j = 1; j <= number_candidate_MRCAs; ++j) {
      int sample_size;
      sample_size = std::numeric_limits<int>::min();
      
      current_statement__ = 740;
      sample_size = get_sample_size(
                      sorted_coal_times_in_oldest_population_time,
                      rvalue(
                        candidate_time_MRCAs_in_model_time_oldest_population,
                        "candidate_time_MRCAs_in_model_time_oldest_population",
                        index_uni(j)), number_of_tips_below,
                      topology, pstream__);
      int size_sub_matrix;
      size_sub_matrix = std::numeric_limits<int>::min();
      
      current_statement__ = 741;
      size_sub_matrix = ((total_sample_size - sample_size) + 1);
      current_statement__ = 742;
      validate_non_negative_index("Q", "(sample_size - 1)", (sample_size - 1));
      current_statement__ = 743;
      validate_non_negative_index("Q", "(sample_size - 1)", (sample_size - 1));
      Eigen::Matrix<local_scalar_t__, -1, -1> Q;
      Q = Eigen::Matrix<local_scalar_t__, -1, -1>((sample_size - 1),
        (sample_size - 1));
      stan::math::fill(Q, DUMMY_VAR__);
      
      local_scalar_t__ current_time;
      current_time = DUMMY_VAR__;
      
      current_statement__ = 745;
      current_time = rvalue(candidate_time_MRCAs_in_model_time,
                       "candidate_time_MRCAs_in_model_time", index_uni(j));
      current_statement__ = 746;
      assign(Q,
        rvalue(Q_coeff_hypoexponential, "Q_coeff_hypoexponential",
          index_min_max(size_sub_matrix, (total_sample_size - 1)),
            index_min_max(size_sub_matrix, (total_sample_size - 1))),
        "assigning variable Q");
      current_statement__ = 747;
      assign(Q,
        matrix_exp(multiply(stan::model::deep_copy(Q), current_time)),
        "assigning variable Q");
      current_statement__ = 748;
      assign(result,
        rvalue(Q, "Q", index_uni(1), index_uni((sample_size - 1))),
        "assigning variable result", index_uni(j));
      current_statement__ = 749;
      cum_sum = (cum_sum + rvalue(result, "result", index_uni(j)));
    }
    current_statement__ = 752;
    assign(result, multiply((1.0 / cum_sum), stan::model::deep_copy(result)),
      "assigning variable result");
    current_statement__ = 753;
    return result;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct compute_probab_vector_functor__ {
template <typename T2__, typename T3__, typename T4__, typename T5__,
typename T6__, typename T8__>
Eigen::Matrix<stan::promote_args_t<stan::value_type_t<T2__>, stan::value_type_t<T3__>,
stan::value_type_t<T4__>, T5__,
T6__, stan::promote_args_t<stan::value_type_t<T8__>>>, -1, 1>
operator()(const int& total_sample_size, const int& number_candidate_MRCAs,
           const T2__& candidate_time_MRCAs_in_model_time,
           const T3__& candidate_time_MRCAs_in_model_time_oldest_population,
           const T4__& Q_coeff_hypoexponential, const T5__& delta,
           const T6__& torigin_in_model_time,
           const std::vector<int>& number_of_tips_below,
           const T8__& sorted_coal_times_in_oldest_population_time,
           const std::vector<std::vector<int>>& topology,
           std::ostream* pstream__)  const
{
return compute_probab_vector(total_sample_size, number_candidate_MRCAs,
         candidate_time_MRCAs_in_model_time,
         candidate_time_MRCAs_in_model_time_oldest_population,
         Q_coeff_hypoexponential, delta, torigin_in_model_time,
         number_of_tips_below, sorted_coal_times_in_oldest_population_time,
         topology, pstream__);
}
};

template <typename T2__, typename T3__, typename T4__, typename T5__,
typename T6__, typename T8__, typename T10__, typename T11__>
stan::promote_args_t<T2__, stan::value_type_t<T3__>,
stan::value_type_t<T4__>, stan::value_type_t<T5__>,
stan::value_type_t<T6__>, stan::promote_args_t<stan::value_type_t<T8__>,
stan::value_type_t<T10__>,
T11__>>
log_likelihood(const int& total_sample_size, const int& N, const T2__& K,
               const T3__& deltas_arg__,
               const T4__& torigins_in_model_time_arg__,
               const T5__& torigins_in_model_time_oldest_population_arg__,
               const T6__& sorted_coalescent_times_in_model_time_oldest_population_arg__,
               const std::vector<int>& number_inmigrants_per_population,
               const T8__& MRCA_times_in_model_time_oldest_population_arg__,
               const std::vector<int>& indexes_father_populations,
               const T10__& pop_sizes_proportion_arg__,
               const T11__& max_torigin_in_model_time_oldest_population,
               const std::vector<std::vector<int>>& topology,
               const std::vector<int>& number_of_tips_below,
               const std::vector<int>& number_of_coalescent_times_below,
               std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T2__,
          stan::value_type_t<T3__>,
          stan::value_type_t<T4__>,
          stan::value_type_t<T5__>,
          stan::value_type_t<T6__>, stan::promote_args_t<stan::value_type_t<T8__>,
          stan::value_type_t<T10__>,
          T11__>>;
  int current_statement__ = 0;
  const auto& deltas = to_ref(deltas_arg__);
  const auto& torigins_in_model_time = to_ref(torigins_in_model_time_arg__);
  const auto& torigins_in_model_time_oldest_population = to_ref(torigins_in_model_time_oldest_population_arg__);
  const auto& sorted_coalescent_times_in_model_time_oldest_population = to_ref(sorted_coalescent_times_in_model_time_oldest_population_arg__);
  const auto& MRCA_times_in_model_time_oldest_population = to_ref(MRCA_times_in_model_time_oldest_population_arg__);
  const auto& pop_sizes_proportion = to_ref(pop_sizes_proportion_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    current_statement__ = 755;
    validate_non_negative_index("sorted_torigins_in_model_time_oldest_population",
                                "N", N);
    Eigen::Matrix<local_scalar_t__, -1, 1> sorted_torigins_in_model_time_oldest_population;
    sorted_torigins_in_model_time_oldest_population = Eigen::Matrix<local_scalar_t__, -1, 1>(N);
    stan::math::fill(sorted_torigins_in_model_time_oldest_population, DUMMY_VAR__);
    
    
    current_statement__ = 756;
    assign(sorted_torigins_in_model_time_oldest_population,
      sort_asc(torigins_in_model_time_oldest_population),
      "assigning variable sorted_torigins_in_model_time_oldest_population");
    current_statement__ = 757;
    validate_non_negative_index("number_internal_nodes_under_MRCA", "N", N);
    std::vector<int> number_internal_nodes_under_MRCA;
    number_internal_nodes_under_MRCA = std::vector<int>(N, std::numeric_limits<int>::min());
    
    
    local_scalar_t__ current_time;
    current_time = DUMMY_VAR__;
    
    local_scalar_t__ log_lik;
    log_lik = DUMMY_VAR__;
    
    current_statement__ = 760;
    log_lik = 0.0;
    int alive_cells;
    alive_cells = std::numeric_limits<int>::min();
    
    local_scalar_t__ term_only_after_first_coal_event;
    term_only_after_first_coal_event = DUMMY_VAR__;
    
    int pos;
    pos = std::numeric_limits<int>::min();
    
    local_scalar_t__ random_unif;
    random_unif = DUMMY_VAR__;
    
    int pos_tips;
    pos_tips = std::numeric_limits<int>::min();
    
    current_statement__ = 766;
    validate_non_negative_index("sample_sizes", "N", N);
    std::vector<int> sample_sizes;
    sample_sizes = std::vector<int>(N, std::numeric_limits<int>::min());
    
    local_scalar_t__ current_max_torigin_in_model_time_oldest_population;
    current_max_torigin_in_model_time_oldest_population = DUMMY_VAR__;
    
    current_statement__ = 768;
    current_max_torigin_in_model_time_oldest_population = max(
                                                            torigins_in_model_time_oldest_population);
    int index_oldest_population;
    index_oldest_population = std::numeric_limits<int>::min();
    
    current_statement__ = 769;
    index_oldest_population = N;
    current_statement__ = 770;
    pos_tips = 1;
    current_statement__ = 804;
    for (int i = 1; i <= N; ++i) {
      int order;
      order = std::numeric_limits<int>::min();
      
      current_statement__ = 771;
      order = find(
                rvalue(torigins_in_model_time_oldest_population,
                  "torigins_in_model_time_oldest_population", index_uni(i)),
                sorted_torigins_in_model_time_oldest_population, pstream__);
      int idx_time_MRCA;
      idx_time_MRCA = std::numeric_limits<int>::min();
      
      current_statement__ = 772;
      idx_time_MRCA = find(
                        rvalue(MRCA_times_in_model_time_oldest_population,
                          "MRCA_times_in_model_time_oldest_population",
                          index_uni(i)),
                        sorted_coalescent_times_in_model_time_oldest_population, pstream__);
      int coal_times_size;
      coal_times_size = std::numeric_limits<int>::min();
      
      current_statement__ = 773;
      coal_times_size = (rvalue(number_of_coalescent_times_below,
                           "number_of_coalescent_times_below",
                           index_uni(rvalue(topology, "topology",
                                       index_uni(idx_time_MRCA), index_uni(1))))
                          -
                          sum(
                            rvalue(number_internal_nodes_under_MRCA,
                              "number_internal_nodes_under_MRCA",
                              index_min_max(1, (order - 1)))));
      int all_times_size;
      all_times_size = std::numeric_limits<int>::min();
      
      current_statement__ = 774;
      all_times_size = (coal_times_size + (order - 1));
      int number_tips_under_MRCA;
      number_tips_under_MRCA = std::numeric_limits<int>::min();
      
      current_statement__ = 775;
      number_tips_under_MRCA = (rvalue(number_of_tips_below,
                                  "number_of_tips_below",
                                  index_uni(rvalue(topology, "topology",
                                              index_uni(idx_time_MRCA),
                                                index_uni(1)))) -
                                 sum(
                                   rvalue(number_of_tips_below,
                                     "number_of_tips_below",
                                     index_min_max(1, (order - 1)))));
      current_statement__ = 776;
      validate_non_negative_index("positions", "total_sample_size",
                                  total_sample_size);
      std::vector<int> positions;
      positions = std::vector<int>(total_sample_size, std::numeric_limits<int>::min());
      
      
      int index_coal_time_above_MRCA;
      index_coal_time_above_MRCA = std::numeric_limits<int>::min();
      
      current_statement__ = 779;
      validate_non_negative_index("sorted_coal_time_with_previous_torigins_in_model_time",
                                  "all_times_size", all_times_size);
      Eigen::Matrix<local_scalar_t__, -1, 1> sorted_coal_time_with_previous_torigins_in_model_time;
      sorted_coal_time_with_previous_torigins_in_model_time = Eigen::Matrix<local_scalar_t__, -1, 1>(all_times_size);
      stan::math::fill(sorted_coal_time_with_previous_torigins_in_model_time, DUMMY_VAR__);
      
      
      current_statement__ = 781;
      validate_non_negative_index("torigin_inmigrants_model_time_oldest_population",
                                  "order - 1", (order - 1));
      Eigen::Matrix<local_scalar_t__, -1, 1> torigin_inmigrants_model_time_oldest_population;
      torigin_inmigrants_model_time_oldest_population = Eigen::Matrix<local_scalar_t__, -1, 1>(
        (order - 1));
      stan::math::fill(torigin_inmigrants_model_time_oldest_population, DUMMY_VAR__);
      
      
      current_statement__ = 783;
      assign(number_internal_nodes_under_MRCA,
        (rvalue(number_of_coalescent_times_below,
           "number_of_coalescent_times_below",
           index_uni(rvalue(topology, "topology",
                       index_uni(idx_time_MRCA), index_uni(1)))) -
          sum(
            rvalue(number_internal_nodes_under_MRCA,
              "number_internal_nodes_under_MRCA",
              index_min_max(1, (order - 1))))),
        "assigning variable number_internal_nodes_under_MRCA", index_uni(i));
      current_statement__ = 799;
      if (logical_gt(order, 1)) {
        current_statement__ = 790;
        assign(torigin_inmigrants_model_time_oldest_population,
          get_torigins_of_inmigrant_populations(i, order, N,
            segment(sorted_torigins_in_model_time_oldest_population, 1,
              (order - 1)), torigins_in_model_time_oldest_population,
            indexes_father_populations, pstream__),
          "assigning variable torigin_inmigrants_model_time_oldest_population");
        current_statement__ = 791;
        assign(sorted_coal_time_with_previous_torigins_in_model_time,
          get_list_coal_times_in_oldest_population_time_below_time_origin(i,
            order, total_sample_size,
            rvalue(number_internal_nodes_under_MRCA,
              "number_internal_nodes_under_MRCA", index_uni(i)),
            rvalue(MRCA_times_in_model_time_oldest_population,
              "MRCA_times_in_model_time_oldest_population", index_uni(i)),
            (Eigen::Matrix<local_scalar_t__,-1,1>(1) <<
            rvalue(MRCA_times_in_model_time_oldest_population,
              "MRCA_times_in_model_time_oldest_population", index_uni(1))).finished(),
            rvalue(torigins_in_model_time_oldest_population,
              "torigins_in_model_time_oldest_population", index_uni(i)),
            sorted_coalescent_times_in_model_time_oldest_population,
            topology, number_of_coalescent_times_below, pstream__),
          "assigning variable sorted_coal_time_with_previous_torigins_in_model_time",
          index_min_max(1, coal_times_size));
        current_statement__ = 792;
        assign(sorted_coal_time_with_previous_torigins_in_model_time,
          torigin_inmigrants_model_time_oldest_population,
          "assigning variable sorted_coal_time_with_previous_torigins_in_model_time",
          index_min_max((coal_times_size + 1), all_times_size));
        current_statement__ = 793;
        assign(sorted_coal_time_with_previous_torigins_in_model_time,
          multiply(
            (rvalue(pop_sizes_proportion, "pop_sizes_proportion",
               index_uni(index_oldest_population)) /
              rvalue(pop_sizes_proportion, "pop_sizes_proportion",
                index_uni(i))),
            stan::model::deep_copy(
              sorted_coal_time_with_previous_torigins_in_model_time)),
          "assigning variable sorted_coal_time_with_previous_torigins_in_model_time");
        current_statement__ = 794;
        assign(sorted_coal_time_with_previous_torigins_in_model_time,
          sort_asc(
            stan::model::deep_copy(
              sorted_coal_time_with_previous_torigins_in_model_time)),
          "assigning variable sorted_coal_time_with_previous_torigins_in_model_time");
        current_statement__ = 795;
        assign(sample_sizes,
          (get_sample_size(
             sorted_coalescent_times_in_model_time_oldest_population,
             rvalue(MRCA_times_in_model_time_oldest_population,
               "MRCA_times_in_model_time_oldest_population", index_uni(i)),
             number_of_tips_below, topology, pstream__) -
            sum(
              rvalue(sample_sizes, "sample_sizes",
                index_min_max(1, (order - 1))))),
          "assigning variable sample_sizes", index_uni(i));
        current_statement__ = 796;
        assign(positions,
          get_tips_positions_below_time(order, number_tips_under_MRCA,
            rvalue(MRCA_times_in_model_time_oldest_population,
              "MRCA_times_in_model_time_oldest_population", index_uni(i)),
            (Eigen::Matrix<local_scalar_t__,-1,1>(1) <<
            rvalue(MRCA_times_in_model_time_oldest_population,
              "MRCA_times_in_model_time_oldest_population", index_uni(1))).finished(),
            rvalue(torigins_in_model_time_oldest_population,
              "torigins_in_model_time_oldest_population", index_uni(i)),
            sorted_coalescent_times_in_model_time_oldest_population,
            topology, total_sample_size, pstream__),
          "assigning variable positions");
        current_statement__ = 797;
        log_lik = (log_lik +
                    log_likelihood_subpopulation(order,
                      rvalue(sample_sizes, "sample_sizes", index_uni(i)), K,
                      rvalue(deltas, "deltas", index_uni(i)),
                      rvalue(torigins_in_model_time,
                        "torigins_in_model_time", index_uni(i)), positions,
                      sorted_coal_time_with_previous_torigins_in_model_time,
                      multiply(
                        (rvalue(pop_sizes_proportion, "pop_sizes_proportion",
                           index_uni(index_oldest_population)) /
                          rvalue(pop_sizes_proportion,
                            "pop_sizes_proportion", index_uni(i))),
                        torigin_inmigrants_model_time_oldest_population), pstream__));
      } else {
        current_statement__ = 784;
        assign(sorted_coal_time_with_previous_torigins_in_model_time,
          get_list_coal_times_in_oldest_population_time_below_time_origin(i,
            order, total_sample_size,
            rvalue(number_internal_nodes_under_MRCA,
              "number_internal_nodes_under_MRCA", index_uni(i)),
            rvalue(MRCA_times_in_model_time_oldest_population,
              "MRCA_times_in_model_time_oldest_population", index_uni(i)),
            rep_vector(-1, 1),
            rvalue(torigins_in_model_time_oldest_population,
              "torigins_in_model_time_oldest_population", index_uni(i)),
            sorted_coalescent_times_in_model_time_oldest_population,
            topology, number_of_coalescent_times_below, pstream__),
          "assigning variable sorted_coal_time_with_previous_torigins_in_model_time");
        current_statement__ = 785;
        assign(sample_sizes,
          get_sample_size(
            sorted_coalescent_times_in_model_time_oldest_population,
            rvalue(MRCA_times_in_model_time_oldest_population,
              "MRCA_times_in_model_time_oldest_population", index_uni(i)),
            number_of_tips_below, topology, pstream__),
          "assigning variable sample_sizes", index_uni(i));
        current_statement__ = 786;
        assign(sorted_coal_time_with_previous_torigins_in_model_time,
          multiply(
            (rvalue(pop_sizes_proportion, "pop_sizes_proportion",
               index_uni(index_oldest_population)) /
              rvalue(pop_sizes_proportion, "pop_sizes_proportion",
                index_uni(i))),
            stan::model::deep_copy(
              sorted_coal_time_with_previous_torigins_in_model_time)),
          "assigning variable sorted_coal_time_with_previous_torigins_in_model_time");
        current_statement__ = 787;
        assign(positions,
          get_tips_positions_below_time(order, number_tips_under_MRCA,
            rvalue(MRCA_times_in_model_time_oldest_population,
              "MRCA_times_in_model_time_oldest_population", index_uni(i)),
            rep_vector(-1, 1),
            rvalue(torigins_in_model_time_oldest_population,
              "torigins_in_model_time_oldest_population", index_uni(i)),
            sorted_coalescent_times_in_model_time_oldest_population,
            topology, total_sample_size, pstream__),
          "assigning variable positions");
        current_statement__ = 788;
        log_lik = (log_lik +
                    log_likelihood_subpopulation(order,
                      rvalue(sample_sizes, "sample_sizes", index_uni(i)), K,
                      rvalue(deltas, "deltas", index_uni(i)),
                      rvalue(torigins_in_model_time,
                        "torigins_in_model_time", index_uni(i)), positions,
                      sorted_coal_time_with_previous_torigins_in_model_time,
                      multiply(
                        (rvalue(pop_sizes_proportion, "pop_sizes_proportion",
                           index_uni(index_oldest_population)) /
                          rvalue(pop_sizes_proportion,
                            "pop_sizes_proportion", index_uni(i))),
                        torigins_in_model_time_oldest_population), pstream__));
      }
      current_statement__ = 802;
      if ((primitive_value(logical_lt(order, (N - 1))) && primitive_value(
          logical_gt(N, 2)))) {
        current_statement__ = 800;
        log_lik = (log_lik +
                    log_prob_father_population(order,
                      index_oldest_population,
                      rvalue(torigins_in_model_time,
                        "torigins_in_model_time",
                        index_uni(rvalue(indexes_father_populations,
                                    "indexes_father_populations",
                                    index_uni(i)))),
                      rvalue(deltas, "deltas",
                        index_uni(rvalue(indexes_father_populations,
                                    "indexes_father_populations",
                                    index_uni(i)))),
                      rvalue(pop_sizes_proportion, "pop_sizes_proportion",
                        index_uni(rvalue(indexes_father_populations,
                                    "indexes_father_populations",
                                    index_uni(i)))),
                      rvalue(torigins_in_model_time,
                        "torigins_in_model_time", index_uni(i)),
                      rvalue(deltas, "deltas", index_uni(i)),
                      rvalue(pop_sizes_proportion, "pop_sizes_proportion",
                        index_uni(i)), deltas, torigins_in_model_time,
                      torigins_in_model_time_oldest_population, K, N,
                      pop_sizes_proportion, pstream__));
      }
    }
    current_statement__ = 805;
    return log_lik;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct log_likelihood_functor__ {
template <typename T2__, typename T3__, typename T4__, typename T5__,
typename T6__, typename T8__, typename T10__, typename T11__>
stan::promote_args_t<T2__, stan::value_type_t<T3__>,
stan::value_type_t<T4__>, stan::value_type_t<T5__>,
stan::value_type_t<T6__>, stan::promote_args_t<stan::value_type_t<T8__>,
stan::value_type_t<T10__>,
T11__>>
operator()(const int& total_sample_size, const int& N, const T2__& K,
           const T3__& deltas, const T4__& torigins_in_model_time,
           const T5__& torigins_in_model_time_oldest_population,
           const T6__& sorted_coalescent_times_in_model_time_oldest_population,
           const std::vector<int>& number_inmigrants_per_population,
           const T8__& MRCA_times_in_model_time_oldest_population,
           const std::vector<int>& indexes_father_populations,
           const T10__& pop_sizes_proportion,
           const T11__& max_torigin_in_model_time_oldest_population,
           const std::vector<std::vector<int>>& topology,
           const std::vector<int>& number_of_tips_below,
           const std::vector<int>& number_of_coalescent_times_below,
           std::ostream* pstream__)  const
{
return log_likelihood(total_sample_size, N, K, deltas,
         torigins_in_model_time, torigins_in_model_time_oldest_population,
         sorted_coalescent_times_in_model_time_oldest_population,
         number_inmigrants_per_population,
         MRCA_times_in_model_time_oldest_population,
         indexes_father_populations, pop_sizes_proportion,
         max_torigin_in_model_time_oldest_population, topology,
         number_of_tips_below, number_of_coalescent_times_below, pstream__);
}
};

template <typename T4__, typename T6__>
Eigen::Matrix<stan::promote_args_t<stan::value_type_t<T4__>,
stan::value_type_t<T6__>>, -1, 1>
compute_branch_lengths_from_coal_times(const int& N,
                                       const std::vector<int>& map_internal_node_topology_row,
                                       const int& total_sample_size,
                                       const std::vector<int>& sample_sizes,
                                       const T4__& coal_times_in_model_time_oldest_population_arg__,
                                       const std::vector<std::vector<int>>& topology,
                                       const T6__& torigins_in_model_time_oldest_population_arg__,
                                       std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<stan::value_type_t<T4__>,
          stan::value_type_t<T6__>>;
  int current_statement__ = 0;
  const auto& coal_times_in_model_time_oldest_population = to_ref(coal_times_in_model_time_oldest_population_arg__);
  const auto& torigins_in_model_time_oldest_population = to_ref(torigins_in_model_time_oldest_population_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    current_statement__ = 807;
    validate_non_negative_index("branch_lengths",
                                "2 * total_sample_size - 2",
                                ((2 * total_sample_size) - 2));
    Eigen::Matrix<local_scalar_t__, -1, 1> branch_lengths;
    branch_lengths = Eigen::Matrix<local_scalar_t__, -1, 1>(((2 *
                                                               total_sample_size)
                                                              - 2));
    stan::math::fill(branch_lengths, DUMMY_VAR__);
    
    current_statement__ = 808;
    assign(branch_lengths, rep_vector(0, ((2 * total_sample_size) - 2)),
      "assigning variable branch_lengths");
    current_statement__ = 809;
    validate_non_negative_index("only_coal_times", "total_sample_size - 1",
                                (total_sample_size - 1));
    Eigen::Matrix<local_scalar_t__, -1, 1> only_coal_times;
    only_coal_times = Eigen::Matrix<local_scalar_t__, -1, 1>((total_sample_size
                                                               - 1));
    stan::math::fill(only_coal_times, DUMMY_VAR__);
    
    current_statement__ = 810;
    assign(only_coal_times, rep_vector(0, (total_sample_size - 1)),
      "assigning variable only_coal_times");
    int pos;
    pos = std::numeric_limits<int>::min();
    
    current_statement__ = 811;
    pos = 1;
    int left;
    left = std::numeric_limits<int>::min();
    
    int right;
    right = std::numeric_limits<int>::min();
    
    int idx_left;
    idx_left = std::numeric_limits<int>::min();
    
    int idx_right;
    idx_right = std::numeric_limits<int>::min();
    
    int cum_sample_size;
    cum_sample_size = std::numeric_limits<int>::min();
    
    current_statement__ = 816;
    cum_sample_size = 0;
    int cum_number_coal;
    cum_number_coal = std::numeric_limits<int>::min();
    
    current_statement__ = 817;
    cum_number_coal = 0;
    int first_position_coal_event;
    first_position_coal_event = std::numeric_limits<int>::min();
    
    current_statement__ = 818;
    first_position_coal_event = 1;
    current_statement__ = 819;
    validate_non_negative_index("number_inmigrants", "N", N);
    std::vector<int> number_inmigrants;
    number_inmigrants = std::vector<int>(N, std::numeric_limits<int>::min());
    
    current_statement__ = 820;
    assign(number_inmigrants, rep_array(0, N),
      "assigning variable number_inmigrants");
    int pos_next_torigin;
    pos_next_torigin = std::numeric_limits<int>::min();
    
    current_statement__ = 821;
    pos_next_torigin = 0;
    int is_torigin;
    is_torigin = std::numeric_limits<int>::min();
    
    current_statement__ = 822;
    is_torigin = 1;
    int pos_torigin;
    pos_torigin = std::numeric_limits<int>::min();
    
    int current_pop;
    current_pop = std::numeric_limits<int>::min();
    
    current_statement__ = 824;
    current_pop = 1;
    current_statement__ = 825;
    assign(only_coal_times, coal_times_in_model_time_oldest_population,
      "assigning variable only_coal_times");
    current_statement__ = 835;
    for (int i = 1; i <= N; ++i) {
      int j;
      j = std::numeric_limits<int>::min();
      
      current_statement__ = 826;
      j = 1;
      current_statement__ = 833;
      while ((primitive_value(logical_lte(pos, (total_sample_size - 1))) &&
             primitive_value(
             logical_lt(
               rvalue(coal_times_in_model_time_oldest_population,
                 "coal_times_in_model_time_oldest_population", index_uni(j)),
               rvalue(torigins_in_model_time_oldest_population,
                 "torigins_in_model_time_oldest_population", index_uni(i)))))) {
        current_statement__ = 830;
        if (logical_eq(
              find(
                rvalue(coal_times_in_model_time_oldest_population,
                  "coal_times_in_model_time_oldest_population", index_uni(j)),
                torigins_in_model_time_oldest_population, pstream__),
              -1)) {
          current_statement__ = 827;
          assign(only_coal_times,
            rvalue(coal_times_in_model_time_oldest_population,
              "coal_times_in_model_time_oldest_population", index_uni(j)),
            "assigning variable only_coal_times", index_uni(pos));
          current_statement__ = 828;
          pos = (pos + 1);
        }
        current_statement__ = 831;
        j = (j + 1);
      }
    }
    current_statement__ = 836;
    pos = 1;
    current_statement__ = 864;
    for (int i = 1; i <= (total_sample_size - 1); ++i) {
      current_statement__ = 837;
      left = rvalue(topology, "topology", index_uni(i), index_uni(2));
      current_statement__ = 838;
      right = rvalue(topology, "topology", index_uni(i), index_uni(3));
      current_statement__ = 850;
      if (logical_lte(left, total_sample_size)) {
        current_statement__ = 847;
        assign(branch_lengths,
          rvalue(only_coal_times, "only_coal_times", index_uni(i)),
          "assigning variable branch_lengths", index_uni(pos));
        current_statement__ = 848;
        pos = (pos + 1);
      } else {
        current_statement__ = 839;
        idx_left = find_integer(left,
                     rvalue(topology, "topology",
                       index_min_max(1, ((total_sample_size - 1) - 1)),
                         index_uni(1)), pstream__);
        current_statement__ = 845;
        if ((primitive_value(logical_neq(idx_left, -1)) && primitive_value(
            logical_gt(
              (rvalue(only_coal_times, "only_coal_times", index_uni(i)) -
                rvalue(only_coal_times, "only_coal_times",
                  index_uni(idx_left))), 0)))) {
          current_statement__ = 842;
          assign(branch_lengths,
            (rvalue(only_coal_times, "only_coal_times", index_uni(i)) -
              rvalue(only_coal_times, "only_coal_times", index_uni(idx_left))),
            "assigning variable branch_lengths", index_uni(pos));
          current_statement__ = 843;
          pos = (pos + 1);
        } else {
          current_statement__ = 840;
          std::stringstream errmsg_stream__;
          errmsg_stream__ << "the branch lengths must be positive";
          throw std::domain_error(errmsg_stream__.str());
        }
      }
      current_statement__ = 862;
      if (logical_lte(right, total_sample_size)) {
        current_statement__ = 859;
        assign(branch_lengths,
          rvalue(only_coal_times, "only_coal_times", index_uni(i)),
          "assigning variable branch_lengths", index_uni(pos));
        current_statement__ = 860;
        pos = (pos + 1);
      } else {
        current_statement__ = 851;
        idx_right = find_integer(right,
                      rvalue(topology, "topology",
                        index_min_max(1, ((total_sample_size - 1) - 1)),
                          index_uni(1)), pstream__);
        current_statement__ = 857;
        if ((primitive_value(logical_neq(idx_right, -1)) && primitive_value(
            logical_gt(
              (rvalue(only_coal_times, "only_coal_times", index_uni(i)) -
                rvalue(only_coal_times, "only_coal_times",
                  index_uni(idx_right))), 0)))) {
          current_statement__ = 854;
          assign(branch_lengths,
            (rvalue(only_coal_times, "only_coal_times", index_uni(i)) -
              rvalue(only_coal_times, "only_coal_times",
                index_uni(idx_right))),
            "assigning variable branch_lengths", index_uni(pos));
          current_statement__ = 855;
          pos = (pos + 1);
        } else {
          current_statement__ = 852;
          std::stringstream errmsg_stream__;
          errmsg_stream__ << "the branch lengths must be positive";
          throw std::domain_error(errmsg_stream__.str());
        }
      }
    }
    current_statement__ = 865;
    return branch_lengths;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct compute_branch_lengths_from_coal_times_functor__ {
template <typename T4__, typename T6__>
Eigen::Matrix<stan::promote_args_t<stan::value_type_t<T4__>,
stan::value_type_t<T6__>>, -1, 1>
operator()(const int& N,
           const std::vector<int>& map_internal_node_topology_row,
           const int& total_sample_size,
           const std::vector<int>& sample_sizes,
           const T4__& coal_times_in_model_time_oldest_population,
           const std::vector<std::vector<int>>& topology,
           const T6__& torigins_in_model_time_oldest_population,
           std::ostream* pstream__)  const
{
return compute_branch_lengths_from_coal_times(N,
         map_internal_node_topology_row, total_sample_size, sample_sizes,
         coal_times_in_model_time_oldest_population, topology,
         torigins_in_model_time_oldest_population, pstream__);
}
};

class stan_model_1_population_JC_genotypes_model final : public model_base_crtp<stan_model_1_population_JC_genotypes_model> {

 private:
  double K;
  int N;
  int total_sample_size;
  int L;
  std::vector<std::vector<int>> genotype_tipdata;
  std::vector<std::vector<int>> topology;
  std::vector<int> tip_association;
  int number_branches;
  std::vector<int> number_of_tips_below;
  std::vector<int> number_of_coalescent_times_below;
  std::vector<std::vector<Eigen::Matrix<double, -1, 1>>> genotype_tip_partials;
  int tip_index;
  Eigen::Matrix<double, -1, -1> indexes_nodes_below__;
  std::vector<int> map_internal_node_topology_row;
  Eigen::Matrix<double, -1, -1> coal_times_to_branch_lengths__;
  Eigen::Matrix<double, -1, 1> w__;
  std::vector<int> v;
  std::vector<int> u;
  int coal_event_times_model_time_oldest_population_1dim__;
  Eigen::Map<Eigen::Matrix<double, -1, -1>> indexes_nodes_below{nullptr, 0, 0};
  Eigen::Map<Eigen::Matrix<double, -1, -1>> coal_times_to_branch_lengths{nullptr, 0, 0};
  Eigen::Map<Eigen::Matrix<double, -1, 1>> w{nullptr, 0};
 
 public:
  ~stan_model_1_population_JC_genotypes_model() { }
  
  inline std::string model_name() const final { return "stan_model_1_population_JC_genotypes_model"; }

  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 d11c003", "stancflags = "};
  }
  
  
  stan_model_1_population_JC_genotypes_model(stan::io::var_context& context__,
                                             unsigned int random_seed__ = 0,
                                             std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ =
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static constexpr const char* function__ = "stan_model_1_population_JC_genotypes_model_namespace::stan_model_1_population_JC_genotypes_model";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      current_statement__ = 49;
      context__.validate_dims("data initialization","K","double",
           std::vector<size_t>{});
      K = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 49;
      K = context__.vals_r("K")[(1 - 1)];
      current_statement__ = 49;
      check_greater_or_equal(function__, "K", K, 0.0);
      current_statement__ = 49;
      check_less_or_equal(function__, "K", K, 2.0);
      current_statement__ = 50;
      context__.validate_dims("data initialization","N","int",
           std::vector<size_t>{});
      N = std::numeric_limits<int>::min();
      
      current_statement__ = 50;
      N = context__.vals_i("N")[(1 - 1)];
      current_statement__ = 50;
      check_greater_or_equal(function__, "N", N, 1);
      current_statement__ = 51;
      context__.validate_dims("data initialization","total_sample_size",
          "int", std::vector<size_t>{});
      total_sample_size = std::numeric_limits<int>::min();
      
      current_statement__ = 51;
      total_sample_size = context__.vals_i("total_sample_size")[(1 - 1)];
      current_statement__ = 51;
      check_greater_or_equal(function__, "total_sample_size",
                             total_sample_size, 2);
      current_statement__ = 52;
      context__.validate_dims("data initialization","L","int",
           std::vector<size_t>{});
      L = std::numeric_limits<int>::min();
      
      current_statement__ = 52;
      L = context__.vals_i("L")[(1 - 1)];
      current_statement__ = 52;
      check_greater_or_equal(function__, "L", L, 0);
      current_statement__ = 53;
      validate_non_negative_index("genotype_tipdata", "total_sample_size",
                                  total_sample_size);
      current_statement__ = 54;
      validate_non_negative_index("genotype_tipdata", "L", L);
      current_statement__ = 55;
      context__.validate_dims("data initialization","genotype_tipdata","int",
           std::vector<size_t>{static_cast<size_t>(total_sample_size),
            static_cast<size_t>(L)});
      genotype_tipdata = std::vector<std::vector<int>>(total_sample_size, std::vector<int>(L, std::numeric_limits<int>::min()));
      
      
      {
        std::vector<int> genotype_tipdata_flat__;
        current_statement__ = 55;
        genotype_tipdata_flat__ = context__.vals_i("genotype_tipdata");
        current_statement__ = 55;
        pos__ = 1;
        current_statement__ = 55;
        for (int sym1__ = 1; sym1__ <= L; ++sym1__) {
          current_statement__ = 55;
          for (int sym2__ = 1; sym2__ <= total_sample_size; ++sym2__) {
            current_statement__ = 55;
            assign(genotype_tipdata, genotype_tipdata_flat__[(pos__ - 1)],
              "assigning variable genotype_tipdata", index_uni(sym2__),
                                                       index_uni(sym1__));
            current_statement__ = 55;
            pos__ = (pos__ + 1);
          }
        }
      }
      current_statement__ = 55;
      for (int sym1__ = 1; sym1__ <= total_sample_size; ++sym1__) {
        current_statement__ = 55;
        for (int sym2__ = 1; sym2__ <= L; ++sym2__) {
          current_statement__ = 55;
          check_greater_or_equal(function__,
                                 "genotype_tipdata[sym1__, sym2__]",
                                 genotype_tipdata[(sym1__ - 1)][(sym2__ - 1)],
                                 0);
        }
      }
      current_statement__ = 55;
      for (int sym1__ = 1; sym1__ <= total_sample_size; ++sym1__) {
        current_statement__ = 55;
        for (int sym2__ = 1; sym2__ <= L; ++sym2__) {
          current_statement__ = 55;
          check_less_or_equal(function__, "genotype_tipdata[sym1__, sym2__]",
                              genotype_tipdata[(sym1__ - 1)][(sym2__ - 1)],
                              16);
        }
      }
      current_statement__ = 56;
      validate_non_negative_index("topology", "total_sample_size - 1",
                                  (total_sample_size - 1));
      current_statement__ = 57;
      context__.validate_dims("data initialization","topology","int",
           std::vector<size_t>{static_cast<size_t>((total_sample_size - 1)),
            static_cast<size_t>(3)});
      topology = std::vector<std::vector<int>>((total_sample_size - 1), std::vector<int>(3, std::numeric_limits<int>::min()));
      
      
      {
        std::vector<int> topology_flat__;
        current_statement__ = 57;
        topology_flat__ = context__.vals_i("topology");
        current_statement__ = 57;
        pos__ = 1;
        current_statement__ = 57;
        for (int sym1__ = 1; sym1__ <= 3; ++sym1__) {
          current_statement__ = 57;
          for (int sym2__ = 1; sym2__ <= (total_sample_size - 1); ++sym2__) {
            current_statement__ = 57;
            assign(topology, topology_flat__[(pos__ - 1)],
              "assigning variable topology", index_uni(sym2__),
                                               index_uni(sym1__));
            current_statement__ = 57;
            pos__ = (pos__ + 1);
          }
        }
      }
      current_statement__ = 57;
      for (int sym1__ = 1; sym1__ <= (total_sample_size - 1); ++sym1__) {
        current_statement__ = 57;
        for (int sym2__ = 1; sym2__ <= 3; ++sym2__) {
          current_statement__ = 57;
          check_greater_or_equal(function__, "topology[sym1__, sym2__]",
                                 topology[(sym1__ - 1)][(sym2__ - 1)], 0);
        }
      }
      current_statement__ = 57;
      for (int sym1__ = 1; sym1__ <= (total_sample_size - 1); ++sym1__) {
        current_statement__ = 57;
        for (int sym2__ = 1; sym2__ <= 3; ++sym2__) {
          current_statement__ = 57;
          check_less_or_equal(function__, "topology[sym1__, sym2__]",
                              topology[(sym1__ - 1)][(sym2__ - 1)],
                              (2 * total_sample_size));
        }
      }
      current_statement__ = 58;
      validate_non_negative_index("tip_association", "total_sample_size",
                                  total_sample_size);
      current_statement__ = 59;
      context__.validate_dims("data initialization","tip_association","int",
           std::vector<size_t>{static_cast<size_t>(total_sample_size)});
      tip_association = std::vector<int>(total_sample_size, std::numeric_limits<int>::min());
      
      
      current_statement__ = 59;
      tip_association = context__.vals_i("tip_association");
      current_statement__ = 59;
      for (int sym1__ = 1; sym1__ <= total_sample_size; ++sym1__) {
        current_statement__ = 59;
        check_greater_or_equal(function__, "tip_association[sym1__]",
                               tip_association[(sym1__ - 1)], 0);
      }
      current_statement__ = 59;
      for (int sym1__ = 1; sym1__ <= total_sample_size; ++sym1__) {
        current_statement__ = 59;
        check_less_or_equal(function__, "tip_association[sym1__]",
                            tip_association[(sym1__ - 1)], total_sample_size);
      }
      current_statement__ = 60;
      number_branches = std::numeric_limits<int>::min();
      
      current_statement__ = 60;
      number_branches = ((2 * total_sample_size) - 2);
      current_statement__ = 61;
      validate_non_negative_index("number_of_tips_below",
                                  "2 * total_sample_size - 1",
                                  ((2 * total_sample_size) - 1));
      current_statement__ = 62;
      number_of_tips_below = std::vector<int>(((2 * total_sample_size) - 1), std::numeric_limits<int>::min());
      
      
      current_statement__ = 63;
      validate_non_negative_index("number_of_coalescent_times_below",
                                  "2 * total_sample_size - 1",
                                  ((2 * total_sample_size) - 1));
      current_statement__ = 64;
      number_of_coalescent_times_below = std::vector<int>(((2 *
                                                             total_sample_size)
                                                            - 1), std::numeric_limits<int>::min());
      
      
      current_statement__ = 65;
      validate_non_negative_index("genotype_tip_partials",
                                  "2 * total_sample_size",
                                  (2 * total_sample_size));
      current_statement__ = 66;
      validate_non_negative_index("genotype_tip_partials", "L", L);
      current_statement__ = 67;
      genotype_tip_partials = std::vector<std::vector<Eigen::Matrix<double, -1, 1>>>(
        (2 * total_sample_size), std::vector<Eigen::Matrix<double, -1, 1>>(L, Eigen::Matrix<double, -1, 1>(16)));
      stan::math::fill(genotype_tip_partials, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 68;
      tip_index = std::numeric_limits<int>::min();
      
      current_statement__ = 69;
      validate_non_negative_index("indexes_nodes_below",
                                  "total_sample_size - 1",
                                  (total_sample_size - 1));
      current_statement__ = 70;
      validate_non_negative_index("indexes_nodes_below",
                                  "2 * total_sample_size - 1",
                                  ((2 * total_sample_size) - 1));
      current_statement__ = 71;
      indexes_nodes_below__ = Eigen::Matrix<double, -1, -1>((total_sample_size
                                                              - 1), ((2 *
                                                                    total_sample_size)
                                                                    - 1));
      new (&indexes_nodes_below) Eigen::Map<Eigen::Matrix<double, -1, -1>>(indexes_nodes_below__.data(),
        (total_sample_size - 1), ((2 * total_sample_size) - 1));
      stan::math::fill(indexes_nodes_below, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 71;
      assign(indexes_nodes_below,
        rep_matrix(0, (total_sample_size - 1), ((2 * total_sample_size) - 1)),
        "assigning variable indexes_nodes_below");
      current_statement__ = 72;
      validate_non_negative_index("map_internal_node_topology_row",
                                  "2 * total_sample_size - 1",
                                  ((2 * total_sample_size) - 1));
      current_statement__ = 73;
      map_internal_node_topology_row = std::vector<int>(((2 *
                                                           total_sample_size)
                                                          - 1), std::numeric_limits<int>::min());
      
      
      current_statement__ = 73;
      assign(map_internal_node_topology_row,
        rep_array(0, ((2 * total_sample_size) - 1)),
        "assigning variable map_internal_node_topology_row");
      current_statement__ = 74;
      validate_non_negative_index("coal_times_to_branch_lengths",
                                  "number_branches", number_branches);
      current_statement__ = 75;
      validate_non_negative_index("coal_times_to_branch_lengths",
                                  "total_sample_size - 1",
                                  (total_sample_size - 1));
      current_statement__ = 76;
      coal_times_to_branch_lengths__ = Eigen::Matrix<double, -1, -1>(number_branches,
        (total_sample_size - 1));
      new (&coal_times_to_branch_lengths) Eigen::Map<Eigen::Matrix<double, -1, -1>>(coal_times_to_branch_lengths__.data(), number_branches,
        (total_sample_size - 1));
      stan::math::fill(coal_times_to_branch_lengths, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 76;
      assign(coal_times_to_branch_lengths,
        rep_matrix(0, number_branches, (total_sample_size - 1)),
        "assigning variable coal_times_to_branch_lengths");
      current_statement__ = 77;
      validate_non_negative_index("w", "3 * total_sample_size - 4",
                                  ((3 * total_sample_size) - 4));
      current_statement__ = 78;
      w__ = Eigen::Matrix<double, -1, 1>(((3 * total_sample_size) - 4));
      new (&w) Eigen::Map<Eigen::Matrix<double, -1, 1>>(w__.data(), ((3 *
                                                                    total_sample_size)
                                                                    - 4));
      stan::math::fill(w, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 79;
      validate_non_negative_index("v", "3 * total_sample_size - 4",
                                  ((3 * total_sample_size) - 4));
      current_statement__ = 80;
      v = std::vector<int>(((3 * total_sample_size) - 4), std::numeric_limits<int>::min());
      
      
      current_statement__ = 81;
      validate_non_negative_index("u", "number_branches + 1",
                                  (number_branches + 1));
      current_statement__ = 82;
      u = std::vector<int>((number_branches + 1), std::numeric_limits<int>::min());
      
      
      current_statement__ = 86;
      for (int i = 1; i <= total_sample_size; ++i) {
        current_statement__ = 83;
        assign(number_of_tips_below, 1,
          "assigning variable number_of_tips_below", index_uni(i));
        current_statement__ = 84;
        assign(number_of_coalescent_times_below, 0,
          "assigning variable number_of_coalescent_times_below", index_uni(i));
      }
      current_statement__ = 93;
      for (int i = 1; i <= (total_sample_size - 1); ++i) {
        current_statement__ = 87;
        assign(number_of_tips_below,
          (rvalue(number_of_tips_below, "number_of_tips_below",
             index_uni(rvalue(topology, "topology",
                         index_uni(i), index_uni(2)))) +
            rvalue(number_of_tips_below, "number_of_tips_below",
              index_uni(rvalue(topology, "topology",
                          index_uni(i), index_uni(3))))),
          "assigning variable number_of_tips_below", index_uni(rvalue(
                                                                 topology,
                                                                 "topology",
                                                                 index_uni(i),
                                                                   index_uni(1))));
        current_statement__ = 88;
        assign(number_of_coalescent_times_below,
          ((rvalue(number_of_coalescent_times_below,
              "number_of_coalescent_times_below",
              index_uni(rvalue(topology, "topology",
                          index_uni(i), index_uni(2)))) +
             rvalue(number_of_coalescent_times_below,
               "number_of_coalescent_times_below",
               index_uni(rvalue(topology, "topology",
                           index_uni(i), index_uni(3))))) + 1),
          "assigning variable number_of_coalescent_times_below", index_uni(
                                                                   rvalue(
                                                                    topology,
                                                                    "topology",
                                                                    index_uni(i),
                                                                    index_uni(1))));
        current_statement__ = 89;
        assign(map_internal_node_topology_row, i,
          "assigning variable map_internal_node_topology_row", index_uni(
                                                                 rvalue(
                                                                   topology,
                                                                   "topology",
                                                                   index_uni(i),
                                                                    index_uni(1))));
        current_statement__ = 90;
        assign(coal_times_to_branch_lengths, 1.0,
          "assigning variable coal_times_to_branch_lengths", index_uni(
                                                               ((2 * i) - 1)),
                                                               index_uni(i));
        current_statement__ = 91;
        assign(coal_times_to_branch_lengths, 1.0,
          "assigning variable coal_times_to_branch_lengths", index_uni(
                                                               (2 * i)),
                                                               index_uni(i));
      }
      current_statement__ = 101;
      for (int i = 1; i <= (total_sample_size - 1); ++i) {
        current_statement__ = 96;
        if (logical_gt(
              rvalue(topology, "topology", index_uni(i), index_uni(2)),
              total_sample_size)) {
          current_statement__ = 94;
          assign(coal_times_to_branch_lengths, -1.0,
            "assigning variable coal_times_to_branch_lengths", index_uni(
                                                                 ((2 * i) -
                                                                   1)),
                                                                 index_uni(
                                                                 rvalue(
                                                                   map_internal_node_topology_row,
                                                                   "map_internal_node_topology_row",
                                                                   index_uni(
                                                                    rvalue(
                                                                    topology,
                                                                    "topology",
                                                                    index_uni(i),
                                                                    index_uni(2))))));
        }
        current_statement__ = 99;
        if (logical_gt(
              rvalue(topology, "topology", index_uni(i), index_uni(3)),
              total_sample_size)) {
          current_statement__ = 97;
          assign(coal_times_to_branch_lengths, -1.0,
            "assigning variable coal_times_to_branch_lengths", index_uni(
                                                                 (2 * i)),
                                                                 index_uni(
                                                                 rvalue(
                                                                   map_internal_node_topology_row,
                                                                   "map_internal_node_topology_row",
                                                                   index_uni(
                                                                    rvalue(
                                                                    topology,
                                                                    "topology",
                                                                    index_uni(i),
                                                                    index_uni(3))))));
        }
      }
      current_statement__ = 102;
      assign(w, csr_extract_w(coal_times_to_branch_lengths),
        "assigning variable w");
      current_statement__ = 103;
      assign(v, csr_extract_v(coal_times_to_branch_lengths),
        "assigning variable v");
      current_statement__ = 104;
      assign(u, csr_extract_u(coal_times_to_branch_lengths),
        "assigning variable u");
      current_statement__ = 111;
      for (int n = 1; n <= total_sample_size; ++n) {
        current_statement__ = 109;
        for (int i = 1; i <= L; ++i) {
          current_statement__ = 105;
          assign(genotype_tip_partials, rep_vector(0.0, 16),
            "assigning variable genotype_tip_partials", index_uni(n),
                                                          index_uni(i));
          current_statement__ = 106;
          tip_index = rvalue(tip_association, "tip_association",
                        index_uni(n));
          current_statement__ = 107;
          assign(genotype_tip_partials, 1.0,
            "assigning variable genotype_tip_partials", index_uni(n),
                                                          index_uni(i),
                                                          index_uni(rvalue(
                                                                    genotype_tipdata,
                                                                    "genotype_tipdata",
                                                                    index_uni(tip_index),
                                                                    index_uni(i))));
        }
      }
      current_statement__ = 121;
      for (int i = 1; i <= (total_sample_size - 1); ++i) {
        current_statement__ = 112;
        assign(map_internal_node_topology_row, i,
          "assigning variable map_internal_node_topology_row", index_uni(
                                                                 rvalue(
                                                                   topology,
                                                                   "topology",
                                                                   index_uni(i),
                                                                    index_uni(1))));
        current_statement__ = 115;
        if (logical_gt(
              rvalue(topology, "topology", index_uni(i), index_uni(2)),
              total_sample_size)) {
          current_statement__ = 113;
          assign(indexes_nodes_below,
            stan::model::deep_copy(
              rvalue(indexes_nodes_below, "indexes_nodes_below",
                index_uni(rvalue(map_internal_node_topology_row,
                            "map_internal_node_topology_row",
                            index_uni(rvalue(topology, "topology",
                                        index_uni(i), index_uni(2))))),
                  index_omni())),
            "assigning variable indexes_nodes_below", index_uni(i),
                                                        index_omni());
        }
        current_statement__ = 118;
        if (logical_gt(
              rvalue(topology, "topology", index_uni(i), index_uni(3)),
              total_sample_size)) {
          current_statement__ = 116;
          assign(indexes_nodes_below,
            add(
              stan::model::deep_copy(
                rvalue(indexes_nodes_below, "indexes_nodes_below",
                  index_uni(i), index_omni())),
              stan::model::deep_copy(
                rvalue(indexes_nodes_below, "indexes_nodes_below",
                  index_uni(rvalue(map_internal_node_topology_row,
                              "map_internal_node_topology_row",
                              index_uni(rvalue(topology, "topology",
                                          index_uni(i), index_uni(3))))),
                    index_omni()))),
            "assigning variable indexes_nodes_below", index_uni(i),
                                                        index_omni());
        }
        current_statement__ = 119;
        assign(indexes_nodes_below, 1,
          "assigning variable indexes_nodes_below", index_uni(i),
                                                      index_uni(rvalue(
                                                                  topology,
                                                                  "topology",
                                                                  index_uni(i),
                                                                    index_uni(1))));
      }
      current_statement__ = 122;
      validate_positive_index("simplex1", "total_sample_size",
                              total_sample_size);
      current_statement__ = 123;
      coal_event_times_model_time_oldest_population_1dim__ = std::numeric_limits<int>::min();
      
      
      current_statement__ = 123;
      coal_event_times_model_time_oldest_population_1dim__ = (total_sample_size
                                                               - 1);
      current_statement__ = 123;
      validate_non_negative_index("coal_event_times_model_time_oldest_population",
                                  "total_sample_size - 1",
                                  coal_event_times_model_time_oldest_population_1dim__);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = 1 + 1 + 1 + 1 + (total_sample_size - 1);
    
  }
  
  template <bool propto__, bool jacobian__ , typename VecR, typename VecI,
  stan::require_vector_like_t<VecR>* = nullptr,
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR> log_prob_impl(VecR& params_r__,
                                                 VecI& params_i__,
                                                 std::ostream* pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    static constexpr const char* function__ = "stan_model_1_population_JC_genotypes_model_namespace::log_prob";
    (void) function__;  // suppress unused var warning
    
    try {
      local_scalar_t__ hyper_parameters_exponential;
      hyper_parameters_exponential = DUMMY_VAR__;
      
      current_statement__ = 1;
      hyper_parameters_exponential = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(
                                       0.0, lp__);
      local_scalar_t__ gamma;
      gamma = DUMMY_VAR__;
      
      current_statement__ = 2;
      gamma = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(
                0.01, lp__);
      local_scalar_t__ torigins_in_model_time_oldest_population;
      torigins_in_model_time_oldest_population = DUMMY_VAR__;
      
      current_statement__ = 3;
      torigins_in_model_time_oldest_population = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(
                                                   0.0, lp__);
      local_scalar_t__ theta;
      theta = DUMMY_VAR__;
      
      current_statement__ = 4;
      theta = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(
                0.0, lp__);
      Eigen::Matrix<local_scalar_t__, -1, 1> simplex1;
      simplex1 = Eigen::Matrix<local_scalar_t__, -1, 1>(total_sample_size);
      stan::math::fill(simplex1, DUMMY_VAR__);
      
      current_statement__ = 5;
      simplex1 = in__.template read_constrain_simplex<Eigen::Matrix<local_scalar_t__, -1, 1>, jacobian__>(
                   lp__, total_sample_size);
      Eigen::Matrix<local_scalar_t__, -1, 1> coal_event_times_model_time_oldest_population;
      coal_event_times_model_time_oldest_population = Eigen::Matrix<local_scalar_t__, -1, 1>(coal_event_times_model_time_oldest_population_1dim__);
      stan::math::fill(coal_event_times_model_time_oldest_population, DUMMY_VAR__);
      
      
      current_statement__ = 7;
      assign(coal_event_times_model_time_oldest_population,
        multiply(head(cumulative_sum(simplex1), (total_sample_size - 1)),
          torigins_in_model_time_oldest_population),
        "assigning variable coal_event_times_model_time_oldest_population");
      current_statement__ = 6;
      check_positive_ordered(function__,
                             "coal_event_times_model_time_oldest_population",
                             coal_event_times_model_time_oldest_population);
      {
        local_scalar_t__ a;
        a = DUMMY_VAR__;
        
        local_scalar_t__ a_divide_by_3;
        a_divide_by_3 = DUMMY_VAR__;
        
        Eigen::Matrix<local_scalar_t__, -1, 1> left;
        left = Eigen::Matrix<local_scalar_t__, -1, 1>(16);
        stan::math::fill(left, DUMMY_VAR__);
        
        Eigen::Matrix<local_scalar_t__, -1, 1> right;
        right = Eigen::Matrix<local_scalar_t__, -1, 1>(16);
        stan::math::fill(right, DUMMY_VAR__);
        
        current_statement__ = 12;
        validate_non_negative_index("partials", "total_sample_size",
                                    total_sample_size);
        current_statement__ = 13;
        validate_non_negative_index("partials", "L", L);
        std::vector<std::vector<Eigen::Matrix<local_scalar_t__, -1, 1>>> partials;
        partials = std::vector<std::vector<Eigen::Matrix<local_scalar_t__, -1, 1>>>(total_sample_size, std::vector<Eigen::Matrix<local_scalar_t__, -1, 1>>(L, Eigen::Matrix<local_scalar_t__, -1, 1>(16)));
        stan::math::fill(partials, DUMMY_VAR__);
        
        current_statement__ = 15;
        validate_non_negative_index("p_matrices", "number_branches",
                                    number_branches);
        std::vector<Eigen::Matrix<local_scalar_t__, -1, -1>> p_matrices;
        p_matrices = std::vector<Eigen::Matrix<local_scalar_t__, -1, -1>>(number_branches, Eigen::Matrix<local_scalar_t__, -1, -1>(16, 16));
        stan::math::fill(p_matrices, DUMMY_VAR__);
        
        current_statement__ = 17;
        validate_non_negative_index("branch_lengths", "number_branches",
                                    number_branches);
        Eigen::Matrix<local_scalar_t__, -1, 1> branch_lengths;
        branch_lengths = Eigen::Matrix<local_scalar_t__, -1, 1>(number_branches);
        stan::math::fill(branch_lengths, DUMMY_VAR__);
        
        current_statement__ = 19;
        validate_non_negative_index("J_diag", "total_sample_size - 1",
                                    (total_sample_size - 1));
        Eigen::Matrix<local_scalar_t__, -1, 1> J_diag;
        J_diag = Eigen::Matrix<local_scalar_t__, -1, 1>((total_sample_size -
                                                          1));
        stan::math::fill(J_diag, DUMMY_VAR__);
        
        current_statement__ = 21;
        validate_non_negative_index("J_diag_theta", "number_branches",
                                    number_branches);
        Eigen::Matrix<local_scalar_t__, -1, 1> J_diag_theta;
        J_diag_theta = Eigen::Matrix<local_scalar_t__, -1, 1>(number_branches);
        stan::math::fill(J_diag_theta, DUMMY_VAR__);
        
        current_statement__ = 23;
        lp_accum__.add(
          gamma_lpdf<propto__>(hyper_parameters_exponential, 0.001, 0.001));
        current_statement__ = 24;
        lp_accum__.add(
          exponential_lpdf<propto__>(gamma, hyper_parameters_exponential));
        current_statement__ = 25;
        lp_accum__.add(exponential_lpdf<propto__>(theta, 1));
        current_statement__ = 26;
        lp_accum__.add(
          conditionalDensityTOrigin_lpdf<propto__>(
            torigins_in_model_time_oldest_population, gamma,
            total_sample_size, pstream__));
        current_statement__ = 27;
        assign(J_diag,
          rep_vector(torigins_in_model_time_oldest_population,
            (total_sample_size - 1)), "assigning variable J_diag");
        current_statement__ = 28;
        lp_accum__.add(sum(stan::math::log(J_diag)));
        current_statement__ = 29;
        lp_accum__.add(
          structured_coalescent_lpdf<propto__>(
            coal_event_times_model_time_oldest_population,
            (Eigen::Matrix<local_scalar_t__,-1,1>(1) << gamma).finished(),
            (Eigen::Matrix<local_scalar_t__,-1,1>(1) <<
            torigins_in_model_time_oldest_population).finished(),
            (Eigen::Matrix<local_scalar_t__,-1,1>(1) <<
            torigins_in_model_time_oldest_population).finished(),
            (Eigen::Matrix<double,-1,1>(1) << 1.0).finished(),
            std::vector<int>{total_sample_size}, K, pstream__));
        current_statement__ = 30;
        assign(branch_lengths,
          csr_matrix_times_vector(number_branches, (total_sample_size - 1),
            w, v, u, coal_event_times_model_time_oldest_population),
          "assigning variable branch_lengths");
        current_statement__ = 37;
        for (int b = 1; b <= number_branches; ++b) {
          current_statement__ = 31;
          a = ((1.0 / 16.0) *
                (1 -
                  stan::math::exp(
                    (((-16 *
                        rvalue(branch_lengths, "branch_lengths",
                          index_uni(b))) * theta) / 15.0))));
          current_statement__ = 32;
          assign(p_matrices, rep_matrix(a, 16, 16),
            "assigning variable p_matrices", index_uni(b));
          current_statement__ = 35;
          for (int i = 1; i <= 16; ++i) {
            current_statement__ = 33;
            assign(p_matrices, (1.0 - (15 * a)),
              "assigning variable p_matrices", index_uni(b), index_uni(i),
                                                 index_uni(i));
          }
        }
        current_statement__ = 48;
        for (int i = 1; i <= L; ++i) {
          current_statement__ = 42;
          for (int n = 1; n <= (total_sample_size - 1); ++n) {
            current_statement__ = 38;
            assign(left,
              (
                 stan::math::eval(logical_gt(
                                    rvalue(topology, "topology",
                                      index_uni(n), index_uni(2)),
                                    total_sample_size)) ?
                 stan::math::promote_scalar<local_scalar_t__>(rvalue(
                                                                partials,
                                                                "partials",
                                                                index_uni(
                                                                  (rvalue(
                                                                    topology,
                                                                    "topology",
                                                                    index_uni(n),
                                                                    index_uni(2))
                                                                    -
                                                                    total_sample_size)),
                                                                  index_uni(i)))
                 :
                 stan::math::promote_scalar<local_scalar_t__>(rvalue(
                                                                genotype_tip_partials,
                                                                "genotype_tip_partials",
                                                                index_uni(
                                                                  rvalue(
                                                                    topology,
                                                                    "topology",
                                                                    index_uni(n),
                                                                    index_uni(2))),
                                                                  index_uni(i)))),
              "assigning variable left");
            current_statement__ = 39;
            assign(right,
              (
                 stan::math::eval(logical_gt(
                                    rvalue(topology, "topology",
                                      index_uni(n), index_uni(3)),
                                    total_sample_size)) ?
                 stan::math::promote_scalar<local_scalar_t__>(rvalue(
                                                                partials,
                                                                "partials",
                                                                index_uni(
                                                                  (rvalue(
                                                                    topology,
                                                                    "topology",
                                                                    index_uni(n),
                                                                    index_uni(3))
                                                                    -
                                                                    total_sample_size)),
                                                                  index_uni(i)))
                 :
                 stan::math::promote_scalar<local_scalar_t__>(rvalue(
                                                                genotype_tip_partials,
                                                                "genotype_tip_partials",
                                                                index_uni(
                                                                  rvalue(
                                                                    topology,
                                                                    "topology",
                                                                    index_uni(n),
                                                                    index_uni(3))),
                                                                  index_uni(i)))),
              "assigning variable right");
            current_statement__ = 40;
            assign(partials,
              elt_multiply(
                multiply(
                  rvalue(p_matrices, "p_matrices", index_uni(((2 * n) - 1))),
                  left),
                multiply(
                  rvalue(p_matrices, "p_matrices", index_uni((2 * n))),
                  right)),
              "assigning variable partials", index_uni((rvalue(topology,
                                                          "topology",
                                                          index_uni(n),
                                                            index_uni(1)) -
                                                         total_sample_size)),
                                               index_uni(i));
          }
          current_statement__ = 45;
          for (int j = 1; j <= 16; ++j) {
            current_statement__ = 43;
            assign(partials,
              (rvalue(
                 rvalue(partials, "partials",
  index_uni((rvalue(topology, "topology",
               index_uni((total_sample_size - 1)), index_uni(1)) -
              total_sample_size)), index_uni(i)),
                 "partials[(topology[(total_sample_size - 1), 1] - total_sample_size), i]",
                 index_uni(j)) * (1.0 / 16)),
              "assigning variable partials", index_uni(((2 *
                                                          total_sample_size)
                                                         - total_sample_size)),
                                               index_uni(i), index_uni(j));
          }
          current_statement__ = 46;
          lp_accum__.add(
            stan::math::log(
              sum(
                rvalue(partials, "partials",
                  index_uni(total_sample_size), index_uni(i)))));
        }
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob_impl()
    
  template <typename RNG, typename VecR, typename VecI, typename VecVar,
  stan::require_vector_like_vt<std::is_floating_point, VecR>* = nullptr,
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr,
  stan::require_std_vector_vt<std::is_floating_point, VecVar>* = nullptr>
  inline void write_array_impl(RNG& base_rng__, VecR& params_r__,
                               VecI& params_i__, VecVar& vars__,
                               const bool emit_transformed_parameters__ = true,
                               const bool emit_generated_quantities__ = true,
                               std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    static constexpr bool propto__ = true;
    (void) propto__;
    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    int current_statement__ = 0;
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    constexpr bool jacobian__ = false;
    (void) DUMMY_VAR__;  // suppress unused var warning
    static constexpr const char* function__ = "stan_model_1_population_JC_genotypes_model_namespace::write_array";
    (void) function__;  // suppress unused var warning
    
    try {
      double hyper_parameters_exponential;
      hyper_parameters_exponential = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      hyper_parameters_exponential = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(
                                       0.0, lp__);
      double gamma;
      gamma = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      gamma = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(
                0.01, lp__);
      double torigins_in_model_time_oldest_population;
      torigins_in_model_time_oldest_population = std::numeric_limits<double>::quiet_NaN();
      
      
      current_statement__ = 3;
      torigins_in_model_time_oldest_population = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(
                                                   0.0, lp__);
      double theta;
      theta = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 4;
      theta = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(
                0.0, lp__);
      Eigen::Matrix<double, -1, 1> simplex1;
      simplex1 = Eigen::Matrix<double, -1, 1>(total_sample_size);
      stan::math::fill(simplex1, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 5;
      simplex1 = in__.template read_constrain_simplex<Eigen::Matrix<local_scalar_t__, -1, 1>, jacobian__>(
                   lp__, total_sample_size);
      Eigen::Matrix<double, -1, 1> coal_event_times_model_time_oldest_population;
      coal_event_times_model_time_oldest_population = Eigen::Matrix<double, -1, 1>(coal_event_times_model_time_oldest_population_1dim__);
      stan::math::fill(coal_event_times_model_time_oldest_population, std::numeric_limits<double>::quiet_NaN());
      
      
      out__.write(hyper_parameters_exponential);
      out__.write(gamma);
      out__.write(torigins_in_model_time_oldest_population);
      out__.write(theta);
      out__.write(simplex1);
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      }
      current_statement__ = 7;
      assign(coal_event_times_model_time_oldest_population,
        multiply(head(cumulative_sum(simplex1), (total_sample_size - 1)),
          torigins_in_model_time_oldest_population),
        "assigning variable coal_event_times_model_time_oldest_population");
      current_statement__ = 6;
      check_positive_ordered(function__,
                             "coal_event_times_model_time_oldest_population",
                             coal_event_times_model_time_oldest_population);
      if (emit_transformed_parameters__) {
        out__.write(coal_event_times_model_time_oldest_population);
      }
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    } // write_array_impl()
    
  template <typename VecVar, typename VecI,
  stan::require_std_vector_t<VecVar>* = nullptr,
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void transform_inits_impl(VecVar& params_r__, VecI& params_i__,
                                   VecVar& vars__,
                                   std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      local_scalar_t__ hyper_parameters_exponential;
      hyper_parameters_exponential = DUMMY_VAR__;
      
      hyper_parameters_exponential = in__.read<local_scalar_t__>();
      out__.write_free_lb(0.0, hyper_parameters_exponential);
      local_scalar_t__ gamma;
      gamma = DUMMY_VAR__;
      
      gamma = in__.read<local_scalar_t__>();
      out__.write_free_lb(0.01, gamma);
      local_scalar_t__ torigins_in_model_time_oldest_population;
      torigins_in_model_time_oldest_population = DUMMY_VAR__;
      
      torigins_in_model_time_oldest_population = in__.read<local_scalar_t__>(
                                                   );
      out__.write_free_lb(0.0, torigins_in_model_time_oldest_population);
      local_scalar_t__ theta;
      theta = DUMMY_VAR__;
      
      theta = in__.read<local_scalar_t__>();
      out__.write_free_lb(0.0, theta);
      Eigen::Matrix<local_scalar_t__, -1, 1> simplex1;
      simplex1 = Eigen::Matrix<local_scalar_t__, -1, 1>(total_sample_size);
      stan::math::fill(simplex1, DUMMY_VAR__);
      
      for (int sym1__ = 1; sym1__ <= total_sample_size; ++sym1__) {
        assign(simplex1, in__.read<local_scalar_t__>(),
          "assigning variable simplex1", index_uni(sym1__));
      }
      out__.write_free_simplex(simplex1);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    } // transform_inits_impl()
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__ = std::vector<std::string>{"hyper_parameters_exponential",
      "gamma", "torigins_in_model_time_oldest_population", "theta",
      "simplex1", "coal_event_times_model_time_oldest_population"};
    
    } // get_param_names()
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{},
      std::vector<size_t>{}, std::vector<size_t>{}, std::vector<size_t>{
      }, std::vector<size_t>{static_cast<size_t>(total_sample_size)},
      std::vector<size_t>{
                          static_cast<size_t>(coal_event_times_model_time_oldest_population_1dim__)
                          }};
    
    } // get_dims()
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "hyper_parameters_exponential");
    param_names__.emplace_back(std::string() + "gamma");
    param_names__.emplace_back(std::string() + "torigins_in_model_time_oldest_population");
    param_names__.emplace_back(std::string() + "theta");
    for (int sym1__ = 1; sym1__ <= total_sample_size; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "simplex1" + '.' + std::to_string(sym1__));
      }
    }
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1;
           sym1__ <= coal_event_times_model_time_oldest_population_1dim__;
           ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "coal_event_times_model_time_oldest_population" + '.' + std::to_string(sym1__));
        }
      }
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // constrained_param_names()
    
  inline void unconstrained_param_names(
                                        std::vector<std::string>& param_names__,
                                        bool emit_transformed_parameters__ = true,
                                        bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "hyper_parameters_exponential");
    param_names__.emplace_back(std::string() + "gamma");
    param_names__.emplace_back(std::string() + "torigins_in_model_time_oldest_population");
    param_names__.emplace_back(std::string() + "theta");
    for (int sym1__ = 1; sym1__ <= (total_sample_size - 1); ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "simplex1" + '.' + std::to_string(sym1__));
      }
    }
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1;
           sym1__ <= coal_event_times_model_time_oldest_population_1dim__;
           ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "coal_event_times_model_time_oldest_population" + '.' + std::to_string(sym1__));
        }
      }
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // unconstrained_param_names()
    
  inline std::string get_constrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"hyper_parameters_exponential\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"gamma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"torigins_in_model_time_oldest_population\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"theta\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"simplex1\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(total_sample_size) + "},\"block\":\"parameters\"},{\"name\":\"coal_event_times_model_time_oldest_population\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(coal_event_times_model_time_oldest_population_1dim__) + "},\"block\":\"transformed_parameters\"}]");
    
    } // get_constrained_sizedtypes()
    
  inline std::string get_unconstrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"hyper_parameters_exponential\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"gamma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"torigins_in_model_time_oldest_population\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"theta\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"simplex1\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string((total_sample_size - 1)) + "},\"block\":\"parameters\"},{\"name\":\"coal_event_times_model_time_oldest_population\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(coal_event_times_model_time_oldest_population_1dim__) + "},\"block\":\"transformed_parameters\"}]");
    
    } // get_unconstrained_sizedtypes()
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    inline void write_array(RNG& base_rng,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                            const bool emit_transformed_parameters = true,
                            const bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      const size_t num_params__ =
  ((((1 + 1) + 1) + 1) + total_sample_size);
      const size_t num_transformed = coal_event_times_model_time_oldest_population_1dim__;
      const size_t num_gen_quantities = 0;
      std::vector<double> vars_vec(num_params__
       + (emit_transformed_parameters * num_transformed)
       + (emit_generated_quantities * num_gen_quantities));
      std::vector<int> params_i;
      write_array_impl(base_rng, params_r, params_i, vars_vec,
          emit_transformed_parameters, emit_generated_quantities, pstream);
      vars = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>>(
        vars_vec.data(), vars_vec.size());
    }

    template <typename RNG>
    inline void write_array(RNG& base_rng, std::vector<double>& params_r,
                            std::vector<int>& params_i,
                            std::vector<double>& vars,
                            bool emit_transformed_parameters = true,
                            bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      const size_t num_params__ =
  ((((1 + 1) + 1) + 1) + total_sample_size);
      const size_t num_transformed = coal_event_times_model_time_oldest_population_1dim__;
      const size_t num_gen_quantities = 0;
      vars.resize(num_params__
        + (emit_transformed_parameters * num_transformed)
        + (emit_generated_quantities * num_gen_quantities));
      write_array_impl(base_rng, params_r, params_i, vars, emit_transformed_parameters, emit_generated_quantities, pstream);
    }

    template <bool propto__, bool jacobian__, typename T_>
    inline T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
                       std::ostream* pstream = nullptr) const {
      Eigen::Matrix<int, -1, 1> params_i;
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }

    template <bool propto__, bool jacobian__, typename T__>
    inline T__ log_prob(std::vector<T__>& params_r,
                        std::vector<int>& params_i,
                        std::ostream* pstream = nullptr) const {
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }


    inline void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream = nullptr) const final {
      std::vector<double> params_r_vec(params_r.size());
      std::vector<int> params_i;
      transform_inits(context, params_i, params_r_vec, pstream);
      params_r = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>>(
        params_r_vec.data(), params_r_vec.size());
    }

  inline void transform_inits(const stan::io::var_context& context,
                              std::vector<int>& params_i,
                              std::vector<double>& vars,
                              std::ostream* pstream__ = nullptr) const {
     constexpr std::array<const char*, 5> names__{"hyper_parameters_exponential",
      "gamma", "torigins_in_model_time_oldest_population", "theta",
      "simplex1"};
      const std::array<Eigen::Index, 5> constrain_param_sizes__{1, 1,
       1, 1, total_sample_size};
      const auto num_constrained_params__ = std::accumulate(
        constrain_param_sizes__.begin(), constrain_param_sizes__.end(), 0);
    
     std::vector<double> params_r_flat__(num_constrained_params__);
     Eigen::Index size_iter__ = 0;
     Eigen::Index flat_iter__ = 0;
     for (auto&& param_name__ : names__) {
       const auto param_vec__ = context.vals_r(param_name__);
       for (Eigen::Index i = 0; i < constrain_param_sizes__[size_iter__]; ++i) {
         params_r_flat__[flat_iter__] = param_vec__[i];
         ++flat_iter__;
       }
       ++size_iter__;
     }
     vars.resize(num_params_r__);
     transform_inits_impl(params_r_flat__, params_i, vars, pstream__);
    } // transform_inits()
    
};
}

using stan_model = stan_model_1_population_JC_genotypes_model_namespace::stan_model_1_population_JC_genotypes_model;

#ifndef USING_R

// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}

stan::math::profile_map& get_stan_profile_data() {
  return stan_model_1_population_JC_genotypes_model_namespace::profiles__;
}

#endif




