/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/
#include<random>

#include"opt_alg.h"
#include"matrix.h"
#include"ode_solver.h"
#include"solution.h"
#include"user_funs.h"


void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		// lab0();
		lab2();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;									// dok�adno��
	int Nmax = 10000;										// maksymalna liczba wywo�a� funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5),						// dolne oraz g�rne ograniczenie
		a(2, 1);											// dok�adne rozwi�zanie optymalne
	solution opt;											// rozwi�zanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo�anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik�w

	//Wahadlo
	Nmax = 1000;											// dok�adno��
	epsilon = 1e-2;											// maksymalna liczba wywo�a� funkcji celu
	lb = 0, ub = 5;											// dolne oraz g�rne ograniczenie
	double teta_opt = 1;									// maksymalne wychylenie wahad�a
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo�anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik�w

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	ofstream Sout("symulacja_lab0.csv");					// definiujemy strumie� do pliku .csv
	Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout.close();											// zamykamy strumie�
	Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
}

void lab1()
{
	//#######################################
	//				   STAŁE
	//#######################################
	double epsilon = 1e-2;									// dok�adno��
	int Nmax = 10000;
	const double gamma = 1e-2;

	//#######################################
	//			 FUN TESTOWA
	//#######################################

 	std::random_device rd;                      // ziarno (sprzętowe, jeśli dostępne)
    std::mt19937 gen(rd());                     // generator Mersenne Twister
    std::uniform_real_distribution<double> dist(-100, 100.0); // równomierny rozkład

    double x0 = 0;
	double* m;
	solution s_fib;
	solution s_lag;

	/* for(int i = 0; i < 100; i++){
		solution::clear_calls();
		x0 = dist(gen);
		m = expansion(ff1T, x0, 1, 1.5, 1000);

		solution::clear_calls();
		s_fib = fib(ff1T, m[0], m[1], epsilon);

		int fib_calls = s_fib.f_calls;
		
		solution::clear_calls();
		s_lag = lag(ff1T, m[0], m[1], epsilon, gamma, Nmax);
		int lag_calls = s_lag.f_calls;


		// lokalne = 0
		// globalne = 1
		
		// expansion
		printf("%lf,%lf,%lf,%.0lf,", 
			x0,       // x0
			m[0],     // a
			m[1],     // b
			m[2]      // f_calls
		);

		// fib
		printf("%.15lf,%.15lf,%d,%d,", 
			m2d(s_fib.x),         // x*
			m2d(s_fib.y),         // y*
			fib_calls,   // f_calls
			(s_fib.x > 10) ? 1 : 0 // lokalne/globalne
		);

		// lag
		printf("%.15lf,%.15lf,%d,%d\n", 
			m2d(s_lag.x),         // x*
			m2d(s_lag.y),         // y*
			lag_calls,   // f_calls
			(s_lag.x > 10) ? 1 : 0 // lokalne/globalne
		);

		delete[] m;
	}
 */
	


	gg1R(matrix(20.03293157));

}

void lab2()
{
	// #####################################################
	// ###################### TESTOWA ######################
	// #####################################################
	/* double epsilon = 1e-2;
	int Nmax = 10000;


	
	double xx[2] = {1.0, 0.0};
	matrix x_zero(2, xx);

	// double y = m2d(ff2T(x_zero, NULL, NULL));
	// printf("%lf \n", y);

	// const double s = 0.2;
	// const double alpha = 0.8;
	// solution::clear_calls();
	// solution rozwiazanie = HJ(ff2T, x_zero, s, alpha, epsilon, Nmax, NULL, NULL);

	// printf("x1: %lf \n", rozwiazanie.x(0));
	// printf("x2: %lf \n", rozwiazanie.x(1));
	// printf("y: %lf \n", rozwiazanie.y(0));

	solution::clear_calls();
	double ss[2] = {0.3,0.3};
	matrix s_zero(2, xx);
	const double alpha = 1.2;
	const double beta = 0.5;
	solution rozwiazanie = Rosen(ff2T, x_zero, s_zero, alpha, beta, epsilon, Nmax, NULL, NULL);
	printf("x1: %lf \n", rozwiazanie.x(0));
	printf("x2: %lf \n", rozwiazanie.x(1));
	printf("y: %lf \n", rozwiazanie.y(0)); */

	// #####################################################
	// ###################### RZECZYW ######################
	// #####################################################
	// double epsilon = 1e-2;
	// int Nmax = 10000;
	// // ------------------
	// // warunki początkowe
	// double xx[2] = {3.0, 9.0};
	// matrix x_zero(2, xx);

	// const double s = 0.5;
	// const double alpha = 0.8;

	// solution::clear_calls();
	// solution rozwiazanie_hj = HJ(ff2R, x_zero, s, alpha, epsilon, Nmax, NULL, NULL);

	// printf("%lf ", rozwiazanie_hj.x(0));
	// printf("%lf ", rozwiazanie_hj.x(1));
	// printf("%lf ", rozwiazanie_hj.y(0));
	// printf("%d ", solution::f_calls);

	// double ss[2] = {0.5,0.5};
	// matrix s_zero(2, ss);
	// const double alpha_rosen = 1.2;
	// const double beta = 0.8;
	// solution::clear_calls();
	// solution rozwiazanie_rosen = Rosen(ff2R, x_zero, s_zero, alpha_rosen, beta, epsilon, Nmax, NULL, NULL);

	// printf("%lf  ", rozwiazanie_rosen.x(0));
	// printf("%lf  ", rozwiazanie_rosen.x(1));
	// printf("%lf  ", rozwiazanie_rosen.y(0));
	// printf("%d ", solution::f_calls);

	matrix k(2,1);

	k(0)=3.0023;
	k(1)=10.834834;
	ff2R(k, NULL, NULL);

	// #####################################################
	// ###################### EXCEL ########################
	// #####################################################
// 	double epsilon = 1e-2;
// 	int Nmax = 10000;

// 	std::random_device rd;                      // ziarno (sprzętowe, jeśli dostępne)
//     std::mt19937 gen(rd());                     // generator Mersenne Twister
//     std::uniform_real_distribution<double> dist(-1.0, 1.0); // równomierny rozkład
	
// 	matrix x0(2,1);

// 	const double s = 0.9;
// 	matrix s0(2,1);
// 	s0(0) = s;
// 	s0(1) = s;
// 	const double alpha = 0.3;
// 	const double alpha_rosen = 1.1;
// 	const double beta_rosen = 0.2;

// 	for(int i = 0; i < 1; i++){
// 		x0(0) = dist(gen);
// 		x0(1) = dist(gen);

// 		x0(0) = 0.21;
// 		x0(1) = 0.37;

// 		solution sol;

// // ############################################################

// 		solution::clear_calls();
// 		sol = HJ(ff2T, x0, s, alpha, epsilon, Nmax, NULL, NULL);
// 		/* 
// 		// 0 x1
// 		// 0 x2
// 		printf("%lf ", x0(0));
// 		printf("%lf ", x0(1));

// 		// x1
// 		// x2
// 		// y
// 		printf("%lf ", sol.x(0));
// 		printf("%lf ", sol.x(1));
// 		printf("%lf ", sol.y(0));
		
// 		// f_calls
// 		printf("%d ", sol.f_calls);
// 		printf("%d ", 0); */

// printf("\n");
// // ############################################################

// 		solution::clear_calls();
// 		sol = Rosen(ff2T, x0, s0, alpha_rosen, beta_rosen, epsilon, Nmax, NULL, NULL);
// 		/* // sol = Rosen();

// 		// x1
// 		// x2
// 		// y
// 		printf("%lf ", sol.x(0));
// 		printf("%lf ", sol.x(1));
// 		printf("%lf ", sol.y(0));
		
// 		// f_calls
// 		printf("%d ", sol.f_calls);
// 		printf("%d ", 0);

// 		printf("\n"); */
// 	}

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
