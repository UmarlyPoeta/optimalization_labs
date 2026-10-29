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