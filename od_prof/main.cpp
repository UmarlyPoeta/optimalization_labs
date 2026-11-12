/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
Data ostatniej modyfikacji: 15.09.2026
*********************************************/

#include"opt_alg.h"

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
	double epsilon;												// dok�adno��
	int Nmax;													// maksymalna liczba wywo�a� funkcji celu
	matrix lb, ub;												// dolne oraz g�rne ograniczenie
	solution opt;												// rozwi�zanie optymalne znalezione przez algorytm
	if (1) // Funkcja testowa
	{
		epsilon = 1e-2;
		Nmax = 10000;
		lb = matrix(2, 1, -5);
		ub = matrix(2, 1, 5);
		matrix a(2, 1);											// dok�adne rozwi�zanie optymalne
		a(0) = -3;
		a(1) = 4;
		opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo�anie procedury optymalizacji
		cout << opt << endl << endl;							// wypisanie wyniku
		solution::clear_calls();								// wyzerowanie licznik�w
		int* n = get_size(opt.ud);
		cout << n[0] << '\t' << n[1] << endl << endl;
	}
	if (0) // Wahad�o
	{
		Nmax = 1000;
		epsilon = 1e-2;
		lb = 0, ub = 5;
		double teta_opt = 1;									// maksymalne wychylenie wahad�a
		opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo�anie procedury optymalizacji
		cout << opt << endl << endl;							// wypisanie wyniku
		solution::clear_calls();								// wyzerowanie licznik�w

		// Zapis zmian warto�ci funckji celu do pliku
		ofstream Sout("optymalicja_lab0.csv");					// definiujemy strumie� do pliku .csv
		Sout << opt.ud;											// zapisujemy wyniki w pliku
		Sout.close();											// zamykamy strumie�

		// Zapis symulacji do pliku csv
		matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
			MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
		matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi�zujemy r�wnanie r�niczkowe
		Sout.open("symulacja_lab0.csv");						// otwieramy plik .csv
		Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
		Sout.close();											// zamykamy strumie�
		Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
		Y[1].~matrix();
	}
	if (1)
	{
		Nmax = 1000;
		epsilon = 1e-2;
		lb = 0, ub = 5;
		
	}
}


static string num(double v)
{
	ostringstream s;
	s.precision(12);
	s << v;
	string t = s.str();
	for (char& ch : t)
		if (ch == '.') ch = ',';
	return t;
}

// dok�adno��
// maksymalna liczba wywo�a� funkcji celu
// dolne oraz g�rne ograniczenie

void lab1()
{
	double epsilon = 1e-9, gamma = 1e-9, d = 1.0;
	int Nmax = 1000;
	double lb = -100, ub = 100;
	double alpha[3] = { 1.5, 2.0, 2.5 };

	double xg = lb, yg = 1e9;
	for (double x = lb; x <= ub; x += 1e-3)
	{
		double y = m2d(ff1T(matrix(x)));
		if (y < yg) { yg = y; xg = x; }
	}
	cout << "Minimum globalne (siatka): x = " << xg << ", y = " << yg << endl;

	ofstream T1("lab1_tabela1.csv");
	T1 << "alpha;x0;a;b;f_calls_ekspansja;fib_x;fib_y;fib_f_calls;fib_globalne;lag_x;lag_y;lag_f_calls;lag_globalne;lag_blad\n";
	for (int ai = 0; ai < 3; ++ai)
		for (int n = 0; n < 100; ++n)
		{
			double x0 = lb + (ub - lb) * m2d(rand_mat(1, 1));
			double* p = nullptr;
			try
			{
				p = expansion(ff1T, x0, d, alpha[ai], Nmax);
			}
			catch (string)
			{
				T1 << num(alpha[ai]) << ";" << num(x0) << ";blad ekspansji\n";
				continue;
			}
			int f_exp = solution::f_calls;
			double a = p[0], b = p[1];
			delete[] p;

			solution F = fib(ff1T, a, b, epsilon);
			int f_fib = solution::f_calls;

			T1 << num(alpha[ai]) << ";" << num(x0) << ";" << num(a) << ";" << num(b) << ";" << f_exp << ";"
				<< num(m2d(F.x)) << ";" << num(m2d(F.y)) << ";" << f_fib << ";"
				<< (fabs(m2d(F.x) - xg) < 0.01 ? 1 : 0) << ";";
			try
			{
				solution L = lag(ff1T, a, b, epsilon, gamma, Nmax);
				T1 << num(m2d(L.x)) << ";" << num(m2d(L.y)) << ";" << solution::f_calls << ";"
					<< (fabs(m2d(L.x) - xg) < 0.01 ? 1 : 0) << ";0\n";
			}
			catch (string)
			{
				T1 << ";;" << solution::f_calls << ";0;1\n";
			}
		}
	T1.close();

	solution F = fib(ff1T, lb, ub, epsilon);
	cout << "Fib bez ekspansji:\n" << F << endl;
	matrix lenF = F.ud, lenL;
	bool lagOK = true;
	try
	{
		solution L = lag(ff1T, lb, ub, epsilon, gamma, Nmax);
		cout << "Lag bez ekspansji:\n" << L << endl;
		lenL = L.ud;
	}
	catch (string ex)
	{
		lagOK = false;
		cout << "Lag bez ekspansji: blad\n" << ex << endl;
	}
	{
		int nF = get_len(lenF), nL = lagOK ? get_len(lenL) : 0;
		ofstream f("lab1_przedzial.csv");
		f << "iteracja;fib;lag\n";
		for (int i = 0; i < max(nF, nL); ++i)
			f << i << ";" << (i < nF ? num(lenF(i)) : "") << ";" << (i < nL ? num(lenL(i)) : "") << "\n";
	}

	{
		matrix Y0(3, 1);
		Y0(0) = 5; Y0(1) = 1; Y0(2) = 20;

		matrix* Yc = solve_ode(df1, 0, 1, 2000, Y0, matrix(50e-4));
		double Tm = 0;
		for (int i = 0; i < get_len(Yc[0]); ++i) Tm = max(Tm, Yc[1](i, 2));
		cout << "Kontrola modelu (D_A = 50 cm2): max T_B = " << Tm << endl;
		delete[] Yc;

		double epsR = 1e-5, gammaR = 1e-5, alphaR = 1.5, dR = 1.0;
		double x0 = 1 + 99 * m2d(rand_mat(1, 1));
		double* p = expansion(ff1R, x0, dR, alphaR, Nmax);
		double a = p[0], b = p[1];
		int f_exp = solution::f_calls;
		delete[] p;

		solution Fr = fib(ff1R, a, b, epsR);
		int f_fib = solution::f_calls;
		solution Lr;
		int f_lag = 0;
		string lagErr;
		try
		{
			Lr = lag(ff1R, a, b, epsR, gammaR, Nmax);
		}
		catch (string ex) { lagErr = ex; }
		f_lag = solution::f_calls;

		ofstream T3("lab1_tabela3.csv");
		T3 << "metoda;x0;a;b;f_calls_ekspansja;D_A_cm2;y;f_calls;D_A_m2\n";
		T3 << "Fibonacci;" << num(x0) << ";" << num(a) << ";" << num(b) << ";" << f_exp << ";"
			<< num(m2d(Fr.x)) << ";" << num(m2d(Fr.y)) << ";" << f_fib << ";" << num(m2d(Fr.x) * 1e-4) << "\n";
		if (lagErr.empty())
			T3 << "Lagrange;" << num(x0) << ";" << num(a) << ";" << num(b) << ";" << f_exp << ";"
			<< num(m2d(Lr.x)) << ";" << num(m2d(Lr.y)) << ";" << f_lag << ";" << num(m2d(Lr.x) * 1e-4) << "\n";
		else
			cout << "Lagrange (zbiorniki): blad: " << lagErr << endl;
		T3.close();

		cout << "Zbiorniki: x0 = " << x0 << ", [a,b] = [" << a << ", " << b << "]\n";
		cout << "Fibonacci:\n" << Fr << endl;
		if (lagErr.empty()) cout << "Lagrange:\n" << Lr << endl;

		matrix* Y = solve_ode(df1, 0, 1, 2000, Y0, matrix(m2d(Fr.x) * 1e-4));
		ofstream Sout("lab1_symulacja.csv");
		Sout << hcat(Y[0], Y[1]);
		Sout.close();
		delete[] Y;
	}
}

struct Wyn2
{
	bool ok = false;
	double x1 = 0, x2 = 0, y = 0;
	int calls = 0;
	bool glob = false;
	vector<vector<double>> traj;
};

// metoda: 0 = Hooke-Jeeves, 1 = Rosenbrock
static Wyn2 run2(int metoda, matrix x0, double s, double alphaHJ, double alphaR, double betaR, double eps, int Nmax)
{
	Wyn2 w;
	try
	{
		solution S;
		if (metoda == 0)
			S = HJ(ff2T, x0, s, alphaHJ, eps, Nmax);
		else
			S = Rosen(ff2T, x0, matrix(2, 1, s), alphaR, betaR, eps, Nmax);
		w.ok = true;
		w.calls = solution::f_calls;
		w.x1 = S.x(0);
		w.x2 = S.x(1);
		w.y = m2d(S.y);
		w.glob = sqrt(w.x1 * w.x1 + w.x2 * w.x2) < 0.01;		// kryterium minimum globalnego (0,0)
		int* sz = get_size(S.ud);
		int rows = sz[0];
		for (int r = 0; r < rows; ++r)
			w.traj.push_back({ S.ud(r, 0), S.ud(r, 1) });
	}
	catch (string)
	{
		w.ok = false;
		w.calls = solution::f_calls;
	}
	return w;
}


void lab2()
{
	double epsilon = 1e-4;											// dokladnosc
	int Nmax = 10000;												// max liczba wywolan funkcji celu
	double alphaHJ = 0.5;											// wsp. zmniejszania kroku (HJ)
	double alphaR = 2.0, betaR = 0.5;								// wsp. ekspansji i kontrakcji (Rosenbrock)
	double s[3] = { 0.05, 0.25, 0.5 };								// trzy dlugosci kroku

	//funkcja testowa
	double sum[2][3][4] = {};										// [metoda][krok][x1,x2,y,calls]
	int cnt[2][3] = {};												// liczba minimow globalnych
	vector<vector<double>> trajHJ, trajR;							// trajektorie do wykresu (krok nr 2, 1. losowanie)

	ofstream T1("lab2_tabela1.csv");
	T1 << "s;lp;x1_0;x2_0;HJ_x1;HJ_x2;HJ_y;HJ_calls;HJ_global;R_x1;R_x2;R_y;R_calls;R_global\n";
	for (int si = 0; si < 3; ++si)
		for (int n = 0; n < 100; ++n)
		{
			matrix x0(2, 1);
			x0(0) = -1 + 2 * m2d(rand_mat(1, 1));
			x0(1) = -1 + 2 * m2d(rand_mat(1, 1));

			Wyn2 r[2];
			for (int m = 0; m < 2; ++m)
				r[m] = run2(m, x0, s[si], alphaHJ, alphaR, betaR, epsilon, Nmax);

			T1 << num(s[si]) << ";" << n + 1 << ";" << num(x0(0)) << ";" << num(x0(1));
			for (int m = 0; m < 2; ++m)
			{
				if (r[m].ok)
				{
					T1 << ";" << num(r[m].x1) << ";" << num(r[m].x2) << ";" << num(r[m].y) << ";" << r[m].calls
						<< ";" << (r[m].glob ? "TAK" : "NIE");
					if (r[m].glob)
					{
						sum[m][si][0] += r[m].x1;
						sum[m][si][1] += r[m].x2;
						sum[m][si][2] += r[m].y;
						sum[m][si][3] += r[m].calls;
						cnt[m][si]++;
					}
				}
				else
					T1 << ";;;;" << r[m].calls << ";NIE";
			}
			T1 << "\n";

			if (si == 1 && n == 0)
			{
				trajHJ = r[0].traj;
				trajR = r[1].traj;
			}
		}
	T1.close();

	//tabela 2 - wartosci srednie (tylko minima globalne)
	ofstream T2("lab2_tabela2.csv");
	T2 << "s;HJ_x1;HJ_x2;HJ_y;HJ_calls;HJ_liczba_global;R_x1;R_x2;R_y;R_calls;R_liczba_global\n";
	for (int si = 0; si < 3; ++si)
	{
		T2 << num(s[si]);
		for (int m = 0; m < 2; ++m)
		{
			for (int k = 0; k < 4; ++k)
				T2 << ";" << (cnt[m][si] ? num(sum[m][si][k] / cnt[m][si]) : string(""));
			T2 << ";" << cnt[m][si];
		}
		T2 << "\n";
	}
	T2.close();

	//trajektorie do wykresu
	{
		ofstream W("lab2_wykres.csv");
		W << "iteracja;HJ_x1;HJ_x2;R_x1;R_x2\n";
		size_t L = max(trajHJ.size(), trajR.size());
		for (size_t i = 0; i < L; ++i)
		{
			W << i << ";";
			if (i < trajHJ.size()) W << num(trajHJ[i][0]) << ";" << num(trajHJ[i][1]); else W << ";";
			W << ";";
			if (i < trajR.size()) W << num(trajR[i][0]) << ";" << num(trajR[i][1]); else W << ";";
			W << "\n";
		}
		W.close();

		//siatka wartosci funkcji celu do narysowania poziomic (wiersz 1: x1, kolumna 1: x2)
		ofstream G("lab2_poziomice.csv");
		int N = 101;
		G << "x2\\x1";
		for (int i = 0; i < N; ++i) G << ";" << num(-1 + 2.0 * i / (N - 1));
		G << "\n";
		for (int j = 0; j < N; ++j)
		{
			double x2 = -1 + 2.0 * j / (N - 1);
			G << num(x2);
			for (int i = 0; i < N; ++i)
			{
				matrix x(2, 1);
				x(0) = -1 + 2.0 * i / (N - 1);
				x(1) = x2;
				G << ";" << num(m2d(ff2T(x, NAN, NAN)));
			}
			G << "\n";
		}
		G.close();
	}

	//problem rzeczywisty
	{
		matrix kt(2, 1);
		kt(0) = 5; kt(1) = 5;
		cout << "Kontrola modelu: Q(5,5) = " << m2d(ff2R(kt, NAN, NAN)) << " (oczekiwane ok. 453.713)" << endl;

		double sR = 1.0, epsR = 1e-3;
		matrix k0(2, 1);
		k0(0) = 20 * m2d(rand_mat(1, 1));
		k0(1) = 20 * m2d(rand_mat(1, 1));

		solution H, R;
		int cH = 0, cR = 0;
		try { H = HJ(ff2R, k0, sR, alphaHJ, epsR, Nmax); cH = solution::f_calls; }
		catch (string ex) { cout << "HJ (ramie): blad\n" << ex << endl; }
		try { R = Rosen(ff2R, k0, matrix(2, 1, sR), alphaR, betaR, epsR, Nmax); cR = solution::f_calls; }
		catch (string ex) { cout << "Rosen (ramie): blad\n" << ex << endl; }

		cout << "k0 = " << k0(0) << ", " << k0(1) << endl;
		cout << "HJ:\n" << H << endl;
		cout << "Rosenbrock:\n" << R << endl;

		ofstream T3("lab2_tabela3.csv");
		T3 << "s;HJ_k1;HJ_k2;HJ_Q;HJ_calls;R_k1;R_k2;R_Q;R_calls;k1_0;k2_0\n";
		T3 << num(sR) << ";" << num(H.x(0)) << ";" << num(H.x(1)) << ";" << num(m2d(H.y)) << ";" << cH << ";"
			<< num(R.x(0)) << ";" << num(R.x(1)) << ";" << num(m2d(R.y)) << ";" << cR << ";"
			<< num(k0(0)) << ";" << num(k0(1)) << "\n";
		T3.close();

		//symulacja dla optymalnych k1, k2: t; alpha_HJ; omega_HJ; alpha_R; omega_R
		matrix Y0(2, 1, 0.0);
		matrix* YH = solve_ode(df2, 0, 0.1, 100, Y0, H.x);
		matrix* YR = solve_ode(df2, 0, 0.1, 100, Y0, R.x);
		ofstream Sout("lab2_symulacja.csv");
		Sout << hcat(hcat(YH[0], YH[1]), YR[1]);
		Sout.close();
		delete[] YH;
		delete[] YR;
	}


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
