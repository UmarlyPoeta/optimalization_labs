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
		lab0();
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

void lab1()
{
	
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
