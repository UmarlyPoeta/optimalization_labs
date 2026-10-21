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
		lab1();
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

void lab1()
{
	double epsilon = 1e-9, gamma = 1e-9, d = 1.0;
	int Nmax = 1000;
	double lb = -100, ub = 100;
	double alpha[3] = { 1.5, 2.0, 2.5 };

	// minimum globalne z przeszukania siatki, do oceny czy metoda je znalazla
	double xg = lb, yg = 1e9;
	for (double x = lb; x <= ub; x += 1e-3)
	{
		double y = m2d(ff1T(matrix(x)));
		if (y < yg) { yg = y; xg = x; }
	}
	cout << "Minimum globalne (siatka): x = " << xg << ", y = " << yg << endl;

	// tabela 1
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

	// bez zawezania przedzialu
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

	// problem rzeczywisty
	{
		matrix Y0(3, 1);
		Y0(0) = 5; Y0(1) = 1; Y0(2) = 20;

		// kontrola modelu: dla D_A = 50 cm^2 max T_B powinno wyjsc ok. 61.2859
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

void lab2()
{

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
