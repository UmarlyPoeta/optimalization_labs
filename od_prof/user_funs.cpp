#include"user_funs.h"

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

matrix ff1T(matrix x, matrix ud1, matrix ud2)
{
	int n = get_len(x);
	matrix y(n, 1);
	for (int i = 0; i < n; ++i)
		y(i) = -cos(0.1 * x(i)) * std::exp(-1 * pow(0.1 * x(i) - 2 * M_PI, 2)) + 0.002 * pow(0.1 * x(i), 2);
	return y;
}


matrix df1(double t, matrix Y, matrix ud1, matrix ud2)
{
	const double a = 0.98, b = 0.63, g = 9.81;
	const double PA = 1.5, PB = 1.0;
	const double DB = 36.5665e-4;
	const double FBin = 0.01, TBin = 20, TA = 90;
	double DA = m2d(ud1);
	double FAout = a * b * DA * sqrt(2 * g * fmax(Y(0), 0.0) / PA);
	double FBout = a * b * DB * sqrt(2 * g * fmax(Y(1), 0.0) / PB);
	matrix dY(3, 1);
	dY(0) = -FAout;
	dY(1) = FAout + FBin - FBout;
	dY(2) = FBin / Y(1) * (TBin - Y(2)) + FAout / Y(1) * (TA - Y(2));
	return dY;
}

matrix ff1R(matrix x, matrix ud1, matrix ud2)
{
	int n = get_len(x);
	matrix y(n, 1);
	for (int k = 0; k < n; ++k)
	{
		double DA = fmax(x(k), 0.0) * 1e-4;
		matrix Y0(3, 1);
		Y0(0) = 5;  Y0(1) = 1;  Y0(2) = 20;
		matrix* Y = solve_ode(df1, 0, 1, 2000, Y0, matrix(DA));
		int m = get_len(Y[0]);
		double Tmax = Y[1](0, 2);
		for (int i = 1; i < m; ++i)
			if (Y[1](i, 2) > Tmax)
				Tmax = Y[1](i, 2);
		y(k) = fabs(Tmax - 50);
		delete[] Y;
	}
	return y;
}

matrix ff2T(matrix x, matrix ud1, matrix ud2)
{
	double v = pow(x(0), 2) + pow(x(1), 2)
		- cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;
	return matrix(v);
}

// Y = [alpha; omega], ud1 = [k1; k2]
matrix df2(double t, matrix Y, matrix ud1, matrix ud2)
{
	double mr = 1.0, mc = 5.0, l = 1.5, b = 0.5;
	double a_ref = 3.14, w_ref = 0.0;							// wartosc zadana (3.14 - zgodnie z wartoscia kontrolna Q = 453.713)
	double I = mr * l * l / 3.0 + mc * l * l;					// moment bezwladnosci
	double M = ud1(0) * (a_ref - Y(0)) + ud1(1) * (w_ref - Y(1));
	matrix dY(2, 1);
	dY(0) = Y(1);
	dY(1) = (M - b * Y(1)) / I;
	return dY;
}

matrix ff2R(matrix x, matrix ud1, matrix ud2)
{
	double a_ref = 3.14, w_ref = 0.0, dt = 0.1;
	matrix Y0(2, 1, 0.0);
	matrix* Y = solve_ode(df2, 0, dt, 100, Y0, x);
	int n = get_len(Y[0]);
	double Q = 0.0;
	for (int i = 0; i < n; ++i)
	{
		double a = Y[1](i, 0), w = Y[1](i, 1);
		double M = x(0) * (a_ref - a) + x(1) * (w_ref - w);
		Q += (10.0 * pow(a_ref - a, 2) + pow(w_ref - w, 2) + pow(M, 2)) * dt;
	}
	delete[] Y;
	return matrix(Q);
}

// ----------------------------------------------------------------- lab 3

// testowa funkcja celu: f(x1,x2) = sin(pi*r') / (pi*r'),  r' = sqrt((x1/pi)^2 + (x2/pi)^2)
matrix ff3T(matrix x, matrix ud1, matrix ud2)
{
	double t = M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2));
	return matrix(sin(t) / t);
}

// ograniczenia g_i(x) <= 0; ud1 = a
static void g3T(matrix x, double a, double g[3])
{
	g[0] = -x(0) + 1;
	g[1] = -x(1) + 1;
	g[2] = sqrt(pow(x(0), 2) + pow(x(1), 2)) - a;
}

// zewnetrzna funkcja kary: F(x) = f(x) + c * sum(max(0, g_i)^2); ud1 = a, ud2 = c
matrix ff3T_zew(matrix x, matrix ud1, matrix ud2)
{
	double g[3];
	g3T(x, m2d(ud1), g);
	double S = 0.0;
	for (int i = 0; i < 3; ++i)
		S += pow(fmax(0.0, g[i]), 2);
	return ff3T(x) + m2d(ud2) * S;
}

// wewnetrzna funkcja kary: F(x) = f(x) - c * sum(1 / g_i); ud1 = a, ud2 = c
// poza obszarem dopuszczalnym (g_i >= 0) zwracana jest bardzo duza wartosc (bariera)
matrix ff3T_wew(matrix x, matrix ud1, matrix ud2)
{
	double g[3];
	g3T(x, m2d(ud1), g);
	double S = 0.0;
	for (int i = 0; i < 3; ++i)
	{
		if (g[i] >= 0)
			return matrix(1e10);
		S -= 1.0 / g[i];
	}
	return ff3T(x) + m2d(ud2) * S;
}

// model lotu pilki z efektem Magnusa; Y = [x; vx; y; vy], ud1 = omega
matrix df3(double t, matrix Y, matrix ud1, matrix ud2)
{
	const double m = 0.6, r = 0.12, C = 0.47, rho = 1.2, g = 9.81;
	const double S = M_PI * r * r;
	double omega = m2d(ud1);
	double vx = Y(1), vy = Y(3);
	double Dx = 0.5 * C * rho * S * vx * fabs(vx);
	double Dy = 0.5 * C * rho * S * vy * fabs(vy);
	double FMx = rho * vy * omega * M_PI * pow(r, 3);
	double FMy = rho * vx * omega * M_PI * pow(r, 3);
	matrix dY(4, 1);
	dY(0) = vx;
	dY(1) = -(Dx + FMx) / m;
	dY(2) = vy;
	dY(3) = -(Dy + FMy) / m - g;
	return dY;
}

// symulacja lotu pilki dla x = [v0x; omega]; zwraca x_end (x dla y najblizszego 0) oraz x dla y najblizszego 50 m
matrix sym3R(matrix x)
{
	matrix Y0(4, 1, 0.0);
	Y0(1) = x(0);
	Y0(2) = 100;
	matrix* Y = solve_ode(df3, 0, 0.01, 7, Y0, matrix(x(1)));
	int n = get_len(Y[0]);
	int i0 = 0, i50 = 0;
	for (int i = 1; i < n; ++i)
	{
		if (fabs(Y[1](i, 2)) < fabs(Y[1](i0, 2))) i0 = i;
		if (fabs(Y[1](i, 2) - 50) < fabs(Y[1](i50, 2) - 50)) i50 = i;
	}
	matrix wynik(3, 1);
	wynik(0) = Y[1](i0, 0);		// x_end
	wynik(1) = Y[1](i50, 0);	// x dla y = 50 m
	wynik(2) = Y[0](i0);		// czas upadku
	delete[] Y;
	return wynik;
}

// problem rzeczywisty: maksymalizacja x_end -> minimalizacja -x_end z zewnetrzna funkcja kary; ud2 = c
matrix ff3R(matrix x, matrix ud1, matrix ud2)
{
	matrix w = sym3R(x);
	double x50 = w(1);
	double g[6] = {
		x(0) - 10, -10 - x(0),		// v0x w [-10, 10]
		x(1) - 15, -15 - x(1),		// omega w [-15, 15]
		2 - x50, x50 - 8			// x(y = 50) w [2, 8]
	};
	double S = 0.0;
	for (int i = 0; i < 6; ++i)
		S += pow(fmax(0.0, g[i]), 2);
	return matrix(-w(0) + m2d(ud2) * S);
}
