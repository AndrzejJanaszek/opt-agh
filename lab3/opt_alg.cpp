#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wej�ciowe:
	// ff - wska�nik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i g�rne ograniczenie
	// epslion - zak��dana dok�adno�� rozwi�zania
	// Nmax - maksymalna liczba wywo�a� funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);									// losujemy macierz Nx1 stosuj�c rozk�ad jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);// przeskalowywujemy rozwi�zanie do przedzia�u [lb, ub]
			Xopt.fit_fun(ff, ud1, ud2);							// obliczmy warto�� funkcji celu
			if (Xopt.y < epsilon)								// sprawdzmy 1. kryterium stopu
			{
				Xopt.flag = 1;									// flaga = 1 ozancza znalezienie rozwi�zanie z zadan� dok�adno�ci�
				break;
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalne liczby wywo�a� funkcji celu
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		int i = 0;
		double f_od_x_i[3] = {0};
		// x_i f(x_i)
		solution s_i[3];

		double* p = new double[3]{ 0, 0, 0 };
		//Tu wpisz kod funkcji

		solution X0(x0);
		solution X1(x0+d);

		X0.fit_fun(ff, ud1, ud2);
		X1.fit_fun(ff, ud1, ud2);

		// todo
		if (X1.y == X0.y){
			p[0] = m2d(X0.x);
			p[1] = m2d(X1.x);
			p[2] = static_cast<double>(X1.f_calls);

			return p;
		}

		if(X1.y > X0.y){
			d = -d;
			X1.x = X0.x + d;
			X1.fit_fun(ff, ud1, ud2);

			if(X1.y >= X0.y){
				p[0] = m2d(X1.x);
				p[1] = m2d(X0.x)-d;
				p[2] = static_cast<double>(X1.f_calls);

				return p;
			}
		}

		s_i[0] = X0;
		s_i[1] = X1;
		do{
			if(X1.f_calls > Nmax){
				throw std::runtime_error("NA DUZO WYWOLAN\n");
			}
			i++;

			// nadpisujemy wartość dla X1 (bez znaczenia czy x1 czy x0 musi byc jakies solution)
			X1.x = X0.x + pow(alpha, i) * d;
			X1.fit_fun(ff, ud1, ud2);
			s_i[(i+1)%3] = X1;
		}while(s_i[i%3].y >= s_i[(i+1)%3].y);

		if (d > 0){
			p[0] = m2d(s_i[(i-1)%3].x);	// i-1
			p[1] = m2d(s_i[(i+1)%3].x);	// i+1
			p[2] = static_cast<double>(X1.f_calls);

			return p;
		}

		p[0] = m2d(s_i[(i+1)%3].x);	// i+1
		p[1] = m2d(s_i[(i-1)%3].x);	// i-1
		p[2] = static_cast<double>(X1.f_calls);
		// std::cout << "arrytm: " << static_cast<double>(X1.f_calls) + static_cast<double>(X0.f_calls) << "\n";
		// std::cout << "p[2]: " << p[2] << "\n";
		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	solution::clear_calls();

	try
	{
		solution Xopt;
		solution s_left;
		solution s_right;
		
		double c = 0, d = 0;
		int k = 0;

		// ----------------------------------------------------------------------------
		// szukanie k -> ciag fibonaciego
		std::vector<double> fib_n;
		fib_n.push_back(1);	// k = 0
		fib_n.push_back(1);	// k = 1
		for(k = 1; fib_n[k] <= ((b-a)/epsilon); k++){
			fib_n.push_back(fib_n[k] + fib_n[k-1]);
		}
		// ----------------------------------------------------------------------------

		c = b - fib_n[k-1]/fib_n[k] * (b - a);
		d = a + b - c;
		// printf("%d,%lf\n", 0, (b-a));


		for(int i = 0; i <= (k-3); i++){
			s_left.x = c;
			s_right.x = d;
			s_left.fit_fun(ff, ud1, ud2);
			s_right.fit_fun(ff, ud1, ud2);

			if(s_left.y < s_right.y){
				// a = a;
				b = d;
			}
			else{
				// b = b;
				a = c;
			}

			c = b - fib_n[k-i-2]/fib_n[k-i-1] * (b - a);
			d = a + b - c;

			// printf("%d,%lf\n", i+1, (b-a));
		}
		s_left.fit_fun(ff, ud1, ud2);
		Xopt.x = s_left.x;
		Xopt.y = s_left.y;
		// s_left.f_calls += s_right.f_calls;
		return s_left;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	solution::clear_calls();

	try
	{
		int i = 0;
		double d = 0, l = 0, m = 0;
		solution Xopt;
		solution s_a;
		solution s_b;
		solution  s_c;

		solution s_d_i[2];
		s_d_i[0].x = 2000;
		//Tu wpisz kod funkcji
		double c = a + ((b-a)/2.0);
		
		s_a.x = a;
		s_b.x = b;
		s_c.x = c;

		// printf("%d,%lf\n", i, (b-a));
		do{
			i++;
			s_a.fit_fun(ff, ud1, ud2);
			s_b.fit_fun(ff, ud1, ud2);
			s_c.fit_fun(ff, ud1, ud2);

			l = m2d(s_a.y) *( pow(m2d(s_b.x),2) - pow(m2d(s_c.x),2)) + m2d(s_b.y) *( pow(m2d(s_c.x),2) - pow(m2d(s_a.x),2)) + m2d(s_c.y) *( pow(m2d(s_a.x),2) - pow(m2d(s_b.x),2));
			m = m2d(s_a.y) *( m2d(s_b.x) - m2d(s_c.x)) + m2d(s_b.y) *( m2d(s_c.x) - m2d(s_a.x)) + m2d(s_c.y) *( m2d(s_a.x) - m2d(s_b.x));

			if (m <= 0){
				throw std::runtime_error("[Error]: m <= 0 in lag function");
			}

			d = 0.5 * l / m;
			s_d_i[i%2].x = d;
			s_d_i[i%2].fit_fun(ff, ud1, ud2);

			if(s_a.x < s_d_i[i%2].x && s_d_i[i%2].x < s_c.x){
				if(s_d_i[i%2].y < s_c.y){
					// a = a;
					// c = d;
					matrix tmp = s_c.x;
					s_c.x = s_d_i[i%2].x;
					// b = c;
					s_b.x = tmp;
				}
				else{
					s_a.x = s_d_i[i%2].x;
					// c = c
					// b = b
				}
			}
			else if(s_c.x < s_d_i[i%2].x && s_d_i[i%2].x < s_b.x){
				// s_d_i[i%2].fit_fun(ff, ud1, ud2);
				if(s_d_i[i%2].y < s_c.y){
					s_a.x = s_c.x;
					s_c.x = s_d_i[i%2].x;
					// b = b
				}
				else{
					// a = a
					// c = c
					s_b.x = s_d_i[i%2].x;
				}
			}
			else{
				throw std::runtime_error("[Error]: double if: in lag function");
			}

			if(s_a.f_calls > Nmax){
				throw std::runtime_error("[Error]: to many f_calls!");
			}

			// printf("%d,%lf\n", i, m2d((s_b.x)-s_a.x));
		}while((m2d(s_b.x - s_a.x) > epsilon) && (fabs( m2d(s_d_i[i%2].x) - m2d(s_d_i[(i-1)%2].x)) > gamma));

		Xopt.x = s_d_i[i%2].x;
		Xopt.y = s_d_i[i%2].y;
		// Xopt.f_calls = s_a.f_calls + s_b.f_calls + s_c.f_calls + s_d_i[0].f_calls + s_d_i[1].f_calls;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

// double norm(matrix &x, matrix &y) {
//     int* size = get_size(x);
//     int n = size[0];  // liczba wierszy (zakładam, że wektor kolumnowy)
//     delete[] size;

//     double sum = 0.0;
//     for (int i = 0; i < n; i++) {
//         double diff = x(i) - y(i);
//         sum += diff * diff;
//     }

//     return sqrt(sum);
// }

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution::clear_calls();
		// solution Xopt;
		// ---------------
		solution sol_xb, sol_x, sol_xb_try;
		sol_x.x = x0;
		
		do{
			sol_xb.x = sol_x.x;
			sol_x = HJ_trial(ff, sol_xb, s, ud1, ud2);
	
			sol_x.fit_fun(ff, ud1, ud2);
			sol_xb.fit_fun(ff, ud1, ud2);
			if(sol_x.y < sol_xb.y){
				do{
					sol_xb_try.x = sol_xb.x;
					sol_xb.x = sol_x.x;
					sol_x.x = 2*sol_xb.x - sol_xb_try.x;
					printf("x: %lf\n", m2d(sol_x.x(0)));
					printf("xb: %lf\n", m2d(sol_xb.x(0)));
	
					if (norm(sol_x.x - sol_xb.x) < 1e-8)
						break;

					sol_x = HJ_trial(ff, sol_xb, s, ud1, ud2);

	
					if(solution::f_calls > Nmax){
						throw std::runtime_error("Za duzo wywołań w funkcji HJ; (1)");
					}
				}while(!(sol_x.y >= sol_xb.y));
				sol_x.x = sol_xb.x;
			}
			else{
				s = alpha * s;
			}
	
			if(solution::f_calls > Nmax){
				throw std::runtime_error("Za duzo wywołań w funkcji HJ; (2)");
			}
		}while(!(s < epsilon));

		return sol_xb;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try {
		const int WYMIAR = 2; // dynamiczny wymiar
		solution sol_x, sol_trial;

		sol_x.x = XB.x;
		sol_x.fit_fun(ff, ud1, ud2);

		for (int j = 0; j < WYMIAR; j++) {
			// f(x + s*e_j)
			sol_trial.x = sol_x.x;
			sol_trial.x(j) += s;
			sol_trial.fit_fun(ff, ud1, ud2);

			if (sol_trial.y < sol_x.y) {
				sol_x = sol_trial;
			}
			else {
				// f(x - s*e_j)
				sol_trial.x = sol_x.x;
				sol_trial.x(j) -= s;
				sol_trial.fit_fun(ff, ud1, ud2);

				if (sol_trial.y < sol_x.y) {
					sol_x = sol_trial;
				}
			}
		}

		return sol_x;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution::clear_calls();
		// solution Xopt;
		//Tu wpisz kod funkcji

		const int WYMIAR = 2;
		solution sol_xb;
		solution sol_step;
		solution sol_x;

		// todo wypełnić jedynkami na przekatnej
		matrix d = ident_mat(WYMIAR);
		matrix lambda(WYMIAR, 1, 0.0);
		matrix p(WYMIAR, 1, 0.0);
		matrix s(WYMIAR, 1, 0.0);
		s = s0;

		sol_xb.x = x0;

		double max = -1;
		do{
			for(int j = 0; j < WYMIAR; j++){
				sol_step.x = sol_xb.x + m2d(s(j) * d(j));
				sol_step.fit_fun(ff, ud1, ud2);
				sol_xb.fit_fun(ff, ud1, ud2);
				if(sol_step.y < sol_xb.y){
					sol_xb.x = sol_xb.x + s(j) * d(j);
					lambda(j) = lambda(j) + s(j);
					s(j) = alpha * s(j);
				}
				else{
					s(j) = -beta * s(j);
					p(j) = p(j) + 1;
				}
			}

			sol_x.x = sol_xb.x;

			bool breakFlag = false;
			for(int j = 0; j < WYMIAR; j++){
				if(lambda(j) == 0 || p(j) == 0){
					breakFlag = true;
					break;
				}
			}

			if(breakFlag == false){
				// zmiana bazy kierunków D
				matrix temp(WYMIAR, WYMIAR, 0.0);
				// wypełnienie macierzy pomocniczej
				for(int j = 0; j < WYMIAR; j++){
					for(int col = 0; col < j+1; col++){
						// wiersz j
						// kolumna col
						temp(j,col) = lambda(j);	// popraw
					}
				}

				matrix Q(WYMIAR, WYMIAR, 0.0);
				Q = d*temp;

				for(int j = 0; j < WYMIAR; j++){
					matrix v = get_col(Q, j);        // pobranie j-tej kolumny
					double v_norm = norm(v);         // norma euklidesowa
					if(v_norm > 1e-12){              // unikamy dzielenia przez 0
						for(int i = 0; i < WYMIAR; i++){
							d(i, j) = v(i) / v_norm; // normalizacja
						}
					} else {
						// jeśli norma ~0, pozostawiamy kierunek bez zmian
						for(int i = 0; i < WYMIAR; i++){
							d(i, j) = v(i);
						}
					}
				}

				lambda = 0;
				p = 0;
				s = s0;
			}

			if(solution::f_calls > Nmax){
				throw std::runtime_error("Za duzo wywołań w funkcji Rosen");
			}

			max = fabs(s(0));
			for(int i = 1; i < WYMIAR; i++){
				if(fabs(s(i)) > max)
					max = fabs(s(i));
			}
		}while(max > epsilon);


		sol_x.fit_fun(ff, ud1, ud2);
		return sol_x;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
