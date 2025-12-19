#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);

// ################################################

matrix ff1T(matrix, matrix = NAN, matrix = NAN);
matrix ff1R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix gg1R(matrix, matrix = NAN, matrix = NAN);

// ################################################

matrix ff2T(matrix x, matrix ud1, matrix ud2);

matrix df2(double, matrix, matrix = NAN, matrix = NAN);
matrix ff2R(matrix, matrix = NAN, matrix = NAN);

// ################################################

matrix ff3T(matrix x, matrix ud1, matrix ud2);
matrix ff3T_zew(matrix x, matrix ud1, matrix ud2);
matrix ff3T_wew(matrix x, matrix ud1, matrix ud2);

matrix df3(double t, matrix Y, matrix ud1, matrix ud2);
matrix ff3R(matrix x, matrix ud1, matrix ud2);

// ################################################

matrix ff4T(matrix x, matrix ud1, matrix ud2);
matrix gf4T(matrix x, matrix ud1, matrix ud2);
matrix hf4T(matrix x, matrix ud1, matrix ud2);

matrix h_phi_4R(matrix x, matrix ud1, matrix ud2);
matrix ff_4R(matrix x, matrix ud1, matrix ud2);
matrix gf_4R(matrix x, matrix ud1, matrix ud2);
matrix p_od_phi_4R(matrix x, matrix ud1, matrix ud2);