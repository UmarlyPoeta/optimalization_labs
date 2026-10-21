#include"opt_alg.h"
#include<vector>

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
		Xopt = rand_mat(N);									// losujemy macierz Nx1 stosuj�c rozk�ad jednostajny na przedziale [0,1]
		for (int i = 0; i < N; ++i)							// przeskalowywujemy rozwi�zanie do przedzia�u [lb, ub]
			Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
		Xopt.fit_fun(ff, ud1, ud2);							// obliczamy warto�� funkcji celu
		Xopt.ud = (Xopt.y);									// warto�� funckji celu zapami�tujemy w macierzy ud
		if (Xopt.y < epsilon)								// sprawdzamy 1. kryterium stopu
		{
			Xopt.flag = 1;
			return Xopt;
		}
		solution X;
		while (true)
		{
			X = rand_mat(N);									// losujemy macierz Nx1 stosuj�c rozk�ad jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)							// przeskalowywujemy rozwi�zanie do przedzia�u [lb, ub]
				X.x(i) = (ub(i) - lb(i)) * X.x(i) + lb(i);
			X.fit_fun(ff, ud1, ud2);							// obliczmy warto�� funkcji celu
			if (X.y < Xopt.y)
			{
				Xopt = X;
				if (Xopt.y < epsilon)							// sprawdzmy 1. kryterium stopu
				{
					Xopt.flag = 1;								// flaga = 1 ozancza znalezienie rozwi�zanie z zadan� dok�adno�ci�
					break;
				}
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalnej liczby wywo�a� funkcji celu
				break;
			}
			Xopt.ud.add_row(Xopt.y);							// warto�� funckji dodajemy jako kolejny wiersz macierzy ud
		}
		Xopt.ud.add_row(Xopt.y);
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
		solution::clear_calls();
		solution X0(x0), X1(x0 + d);
		X0.fit_fun(ff, ud1, ud2);
		X1.fit_fun(ff, ud1, ud2);
		if (X0.y == X1.y)
			return new double[2] { x0, x0 + d };
		if (X1.y > X0.y)
		{
			d = -d;
			X1 = solution(x0 + d);
			X1.fit_fun(ff, ud1, ud2);
			if (X1.y >= X0.y)
				return new double[2] { m2d(X1.x), x0 - d };
		}
		vector<solution> X;
		X.push_back(X0);
		X.push_back(X1);
		int i = 0;
		while (true)
		{
			if (solution::f_calls > Nmax)
				throw string("przekroczono maksymalna liczbe wywolan funkcji celu");
			++i;
			solution Xn(x0 + pow(alpha, i) * d);
			Xn.fit_fun(ff, ud1, ud2);
			X.push_back(Xn);
			if (X[i].y <= X[i + 1].y)
				break;
		}
		if (d > 0)
			return new double[2] { m2d(X[i - 1].x), m2d(X[i + 1].x) };
		return new double[2] { m2d(X[i + 1].x), m2d(X[i - 1].x) };
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution::clear_calls();

		vector<double> phi = { 1, 1 };
		while (phi.back() <= (b - a) / epsilon)
			phi.push_back(phi[phi.size() - 1] + phi[phi.size() - 2]);
		int k = static_cast<int>(phi.size()) - 1;

		double c = b - phi[k - 1] / phi[k] * (b - a);
		double d = a + b - c;
		solution C(c), D(d);
		C.fit_fun(ff, ud1, ud2);
		D.fit_fun(ff, ud1, ud2);
		Xopt.ud = matrix(b - a);

		for (int i = 0; i <= k - 3; ++i)
		{
			if (C.y < D.y)
			{
				b = d;  d = c;  D = C;
				c = b - phi[k - i - 2] / phi[k - i - 1] * (b - a);
				C = solution(c);
				C.fit_fun(ff, ud1, ud2);
			}
			else
			{
				a = c;  c = d;  C = D;
				d = a + b - c;
				D = solution(d);
				D.fit_fun(ff, ud1, ud2);
			}
			Xopt.ud.add_row(b - a);
		}
		Xopt.x = C.x;
		Xopt.y = C.y;
		Xopt.flag = 0;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}
}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution::clear_calls();

		double c = (a + b) / 2.0;
		double d = c, d_prev = c;
		solution A(a), B(b), C(c), D(c);
		A.fit_fun(ff, ud1, ud2);
		B.fit_fun(ff, ud1, ud2);
		C.fit_fun(ff, ud1, ud2);
		double fa = m2d(A.y), fb = m2d(B.y), fc = m2d(C.y), fd = fc;
		Xopt.ud = matrix(b - a);
		int i = 0;

		while (true)
		{
			double l = fa * (b * b - c * c) + fb * (c * c - a * a) + fc * (a * a - b * b);
			double m = fa * (b - c) + fb * (c - a) + fc * (a - b);
			if (m <= 0)
				throw string("m <= 0 (parabola nie ma minimum)");

			d_prev = d;
			d = 0.5 * l / m;

			if (a < d && d < c)
			{
				D = solution(d);
				D.fit_fun(ff, ud1, ud2);
				fd = m2d(D.y);
				if (fd < fc) { b = c; fb = fc; c = d; fc = fd; }
				else { a = d; fa = fd; }
			}
			else if (c < d && d < b)
			{
				D = solution(d);
				D.fit_fun(ff, ud1, ud2);
				fd = m2d(D.y);
				if (fd < fc) { a = c; fa = fc; c = d; fc = fd; }
				else { b = d; fb = fd; }
			}
			else
				throw string("d poza przedzialem (a, b)");

			++i;
			Xopt.ud.add_row(b - a);
			if (solution::f_calls > Nmax)
				throw string("przekroczono maksymalna liczbe wywolan funkcji celu");
			if (b - a < epsilon || (i > 1 && fabs(d - d_prev) < gamma))
				break;
		}
		Xopt.x = D.x;
		Xopt.y = D.y;
		Xopt.flag = 0;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
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
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
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
