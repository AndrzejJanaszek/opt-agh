#include"user_funs.h"
#include <math.h>	//todo check on linux

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
// #include <cmath>
// #include <corecrt_math_defines.h>
// #include <corecrt_math.h>			// for windows
// #include <corecrt_math_defines.h>	// for windows

matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera warto�� funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera wsp�rz�dne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera warto�� funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	int n = get_len(Y[0]);									// d�ugo�� rozwi�zania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahad�a
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// warto�� funkcji celu (ud1 to za�o�one maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z po�o�enia to pr�dko��
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z pr�dko�ci to przyspieszenie
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2){
	return -cos(0.1*x(0)) * exp(-pow(0.1*x(0) - 2*M_PI , 2)) + 0.002*pow(0.1*x(0), 2);
}

// funkcja zwraca maksymalną temperaturę w trakcie symulacji dla podanego przekroju D_a
matrix gg1R(matrix x, matrix ud1, matrix ud2){
	// wspolczynniki
	constexpr double a = 0.98;
	constexpr double b = 0.63;
	constexpr double g = 9.81;
	
	// warunki początkowe
	// zbiornik A
	double V_a = 5.0;						// [m3]
	constexpr double P_a = 2.0;					// [m2]
	// przeliczenie przekroju z [cm2] na [m2]
	double D_a =  m2d(x) / 10000.0;			// [m2]
	constexpr double T_a = 95.0;				// [95*C]
	
	// zbiornik B
	double V_b = 1.0;						// [m3]
	constexpr double P_b = 1.0;					// [m2]
	// przeliczenie przekroju z [cm2] na [m2]
	constexpr double D_b = 36.5665 / 10000.0;	// [m2]
	double T = 20.0; 						// [*C]
	
	constexpr double T_kran = 20.0;				// [*C]
	// 1L = 1dm3 = 0.001 m3
	constexpr double F_kran = 10 * 0.001; 		//[m3/s]
	
	// czas symulacji [double celowo]
	const double time_start = 0.0;
	const double time_end = 2000.0;
	const double time_step = 1.0;

	double dV_a, dV_b, dT, V_in, T_in;

	double T_max = -300;
	// printf("%lf,%lf,%lf,%lf\n", time, V_a, V_b, T);
	
	for(double time = time_start; time < time_end; time+=time_step) {
		// przeliczenie objętości dla zbiornika A
		if(V_a > 0){
			dV_a = (-a * b * D_a * sqrt(2 * g * V_a / P_a))*time_step;
			V_a = V_a + dV_a;
			if(V_a < 0){
				dV_a = 0;
				V_a = 0;
			}
		}
		// przeliczenie objętości dla zbiornika B
		dV_b = (-a * b * D_b * sqrt(2 * g * V_b / P_b))*time_step;
		V_b = V_b + dV_b;

		// suma wpływających płynów do pojemnika B
		// [] = (z poj. A[m3]) + (z kranu[m3])
		V_in = -dV_a + (F_kran*time_step);

		// Średnia temperatura wpływających płynów do pojemnika B
		T_in = (-dV_a * T_a + (F_kran*time_step) * T_kran) / V_in;

		// zmiana temperatury płynu w pojemniku B
		dT =(V_in / V_b * (T_in - T)) * time_step;	//*dt pominiete bo 1[s]
		V_b += V_in;
		T = T + dT;
		// printf("time: %lf\tT: %lf\n", time, T);

		if(T > T_max){
			T_max = T;
		}

		if(time >= 1000)
		printf("%lf,%lf,%lf,%lf\n", time+1, V_a, V_b, T);
	}
	return matrix(T_max);
}

matrix ff1R(matrix x, matrix ud1, matrix ud2){
	return fabs(m2d(gg1R(x, ud1, ud2)) - 50);
}

matrix ff2T(matrix x, matrix ud1, matrix ud2){
	return pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;
}


matrix df2(double t, matrix Y, matrix ud1, matrix ud2){
	double M = ud2(0) * (ud1(0) - Y(0)) + ud2(1) * (ud1(1) - Y(1));

	constexpr double mr = 1;
	constexpr double mc = 5;
	constexpr double l = 2;
	constexpr double b = 0.25;
	constexpr double I = 1.0 / 3.0 * mr * l*l + mc * l*l;

	matrix dY(2,1);
	dY(0) = Y(1);
	dY(1) = ( M - b*Y(1) ) / I;

	return dY;
}

matrix ff2R(matrix x, matrix ud1, matrix ud2){
	matrix y = 0;

	matrix Y0(2,1);
	Y0(0)=0;
	Y0(1)=0;

	matrix Y_ref(2,1);
	Y_ref(0) = M_PI;
	Y_ref(1) = 0;
	
	matrix *Y = solve_ode(df2, 0 , 0.1, 100, Y0, Y_ref, x);

	int n = get_len(Y[0]);
	
	for(int i = 0; i < n; i++){
		double Mt = x(0)*( Y_ref(0) - Y[1](i,0)) + 
					x(1)*( Y_ref(1) - Y[1](i,1));
		y = y + 
		10 * pow(Y_ref(0) - Y[1](i,0),2) + 
		pow(Y_ref(1) - Y[1](i,1),2) +
		pow(Mt, 2);

		printf("%lf ", Y[1](i,0));
		printf("%lf", Y[1](i,1));
		printf("\n");
	}

	delete[] Y;

	y = 0.1 * y;
	return y;
	// printf("Jakiś śmieszny wynik całki: %lf\n", m2d(y));
	// printf("k1: %lf\n", m2d(x(0)) );
	// printf("k2: %lf\n", m2d(x(1)) );
}

matrix ff3T(matrix x, matrix ud1, matrix ud2){
	double r = std::sqrt( (x(0)/M_PI)*(x(0)/M_PI) + (x(1)/M_PI)*(x(1)/M_PI) );
    double denom = M_PI * r;

    double val = (denom == 0.0)
        ? 1.0           // granica sin(t)/t -> 1 dla t→0
        : std::sin(M_PI * r) / denom;

	matrix res;
	res = val;
	return res;
}
matrix ff3T_zew(matrix x, matrix ud1, matrix ud2){
	double r = std::sqrt( (x(0)/M_PI)*(x(0)/M_PI) + (x(1)/M_PI)*(x(1)/M_PI) );
    double denom = M_PI * r;

    double val = (denom == 0.0)
        ? 1.0           // granica sin(t)/t -> 1 dla t→0
        : std::sin(M_PI * r) / denom;

	matrix res;
	res = val;

	if(-x(0) + 1 > 0){
		res = res + ud2 * pow(-x(0) + 1, 2);
	}
	if(-x(1) + 1 > 0){
		res = res + ud2 * pow(-x(1) + 1, 2);
	}
	if(norm(x) - ud1 > 0){
		res = res + ud2 * pow(norm(x) - ud1, 2);
	}


	return res;
}
matrix ff3T_wew(matrix x, matrix ud1, matrix ud2){
	double r = std::sqrt( (x(0)/M_PI)*(x(0)/M_PI) + (x(1)/M_PI)*(x(1)/M_PI) );
    double denom = M_PI * r;

    double val = (denom == 0.0)
        ? 1.0           // granica sin(t)/t -> 1 dla t→0
        : std::sin(M_PI * r) / denom;

	matrix res;
	res = val;

	// w warunkach daliśmy >= itd żeby nie było błędu dzielenia przez 0
	#define DBL_MAX 1.7976931348623158e+308
	if(-x(0) + 1 >= 0){
		res = DBL_MAX;
	}
	else{
		res = res - ud2 / (-x(0) + 1);
	}
	if(-x(1) + 1 >= 0){
		res = DBL_MAX;
	}
	else{
		res = res - ud2 / (-x(1) + 1);
	}
	if(norm(x) - ud1 >= 0){
		res = DBL_MAX;
	}
	else{
		res = res - ud2 / (norm(x) - ud1);
	}

	return res;
}

matrix df3(double t, matrix Y, matrix ud1, matrix ud2){
	// dane wstepne
	const double r = 0.12;	// 12 cm
	const double m = 0.6; 	// 600g = masa
	const double y0 = 100;

	// stałe
	const double g = 9.81;
	const double C = 0.47;
	const double ro = 1.2;
	const double S = M_PI * r*r;

	double vx =  Y(1);	// dx/dt
	double vy = Y(3);	// dx/dt
	double Dx = 0.5 * C * ro * S * vx * fabs(vx);
	double Dy = 0.5 * C * ro * S * vy * fabs(vy);
	double omega = ud1(1);	// TODO
	double Fmx = ro * vy * omega * M_PI * r*r*r;
	double Fmy = ro * vx * omega * M_PI * r*r*r;

	matrix dY(4,1);
	dY(0) = Y(1);
	dY(1) = (-Dx - Fmx)/m;
	dY(2) = Y(3);
	dY(3) = (-m*g -Dy - Fmy)/m;

	return dY;
}

matrix ff3R(matrix x, matrix ud1, matrix ud2){
	matrix Y0(4,1);
	Y0(0) = 0;
	Y0(1) = x(0);
	Y0(2) = 100;
	Y0(3) = 0;

	matrix *Y = solve_ode(df3, 0 , 0.01, 7, Y0, x, NULL);	//todo

	int i0 = 0;
	int i50 = 0;

	int n = get_len(Y[0]);

	matrix y = 0;

	for(int i = 0; i < n; i++){
		if( fabs( Y[1](i,2) - 50 ) < fabs( Y[1](i50,2) - 50 ) ){	//todo
			i50 = i;
		}
		if( fabs( Y[1](i,2) ) < fabs( Y[1](i0,2) ) ){
			i0 = i;
			y = -Y[1](i0,0);
		}
		if( fabs( x(0) ) - 10 > 0){
			// printf("x0");
			y = y + ud2 * pow( fabs( x(0) ) - 10, 2 );
		}
		if( fabs( x(1) ) - 10 > 0){
			// printf("x1");
			y = y + ud2 * pow( fabs( x(1) ) - 10, 2 );
		}
		if( fabs( Y[1](i50,0) - 5 ) - 2 > 0){
			y = y + ud2 * pow( fabs( Y[1](i50,0) - 5 ) - 2, 2 );
			// y = 99999999;
		}

		printf("%lf ", i*0.01);
		printf("%lf ", Y[1](i,0));
		printf("%lf ", Y[1](i,2));
		printf("\n");
	}

	// printf("Przez kosz: %lf\n",Y[1](i50,0));
	return y;
}

// jeżeli dynamiczne h - dla liniowej phi()
// ud2 => 
matrix ff4T(matrix x, matrix ud1, matrix ud2){

	if(ud2.m < 2){
		return 1.0/6.0 * pow(x(0), 6) - 1.05*pow(x(0), 4) + 2*pow(x(0), 2) + pow(x(1), 2) + x(0)*x(1);
	}
	return ff4T(ud2[0] + x*ud2[1], NULL, NULL);

}
matrix gf4T(matrix x, matrix ud1, matrix ud2){
	matrix res(2,1);
	res(0) = pow(x(0), 5) - 4.20*pow(x(0), 3) + 4*x(0) + x(1);
	res(1) = 2*x(1) + x(0);
	return res;
}
matrix hf4T(matrix x, matrix ud1, matrix ud2){
	matrix res(2,2);
	// [ 5*pow(x(0), 4) - 12.6*pow(x(0), 2) + 4		1 ]
	// [ 1											2 ]

	// dxdx
	res(0,0) = 5*pow(x(0), 4) - 12.6*pow(x(0), 2) + 4;
	// dxdy
	res(0,1) = 1;
	// dydx
	res(1,0) = 1;
	// dydy
	res(1,1) = 2;

	return res;
}

// ud2 to wwektor phi
matrix h_phi_4R(matrix x, matrix ud1, matrix ud2){
	return 1/(1 + exp( m2d(-trans(ud2)*x) ));
}

// phi jest w x
// w ud1 jest xi
// w ud2 jest yi
matrix ff_4R(matrix x, matrix ud1, matrix ud2){
	int m = 100;
	double sum = 0;
	for(int i = 0; i < m; i++){
		sum += ud2(i) * log(h_phi_4R(ud1(i),NULL,x)(0)) + (1 - ud2(i)) * log(1 - h_phi_4R(ud1(i),NULL,x)(0));
	}
	sum = -sum / (double)m;

	matrix res = sum;
	return res;
}
