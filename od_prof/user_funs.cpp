#include"user_funs.h"

matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1),
		MT = matrix(2, new double[2] { m2d(x), 0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 20, Y0, ud1, MT);
	int n = get_len(Y[0]);
	double teta_max = Y[1](0, 0);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));
	delete[] Y;
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	double m = 1, l = 0.5, b = 0.5, g = 9.81;
	double I = m * pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2)
{
	double t = 0.1 * m2d(x);
	matrix y;
	y = -cos(t) * std::exp(-1 * pow(t - 2 * M_PI, 2)) + 0.002 * pow(t, 2);
	return y;
}

// Y = [V_A; V_B; T_B], ud1 = D_A w m^2
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

// x = D_A w cm^2
matrix ff1R(matrix x, matrix ud1, matrix ud2)
{
	double DA = fmax(m2d(x), 0.0) * 1e-4;
	matrix Y0(3, 1);
	Y0(0) = 5;  Y0(1) = 1;  Y0(2) = 20;
	matrix* Y = solve_ode(df1, 0, 1, 2000, Y0, matrix(DA));
	int n = get_len(Y[0]);
	double Tmax = Y[1](0, 2);
	for (int i = 1; i < n; ++i)
		if (Y[1](i, 2) > Tmax)
			Tmax = Y[1](i, 2);
	matrix y;
	y = fabs(Tmax - 50);
	delete[] Y;
	return y;
}