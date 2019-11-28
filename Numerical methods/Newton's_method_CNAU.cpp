#include <iostream>
#include <iomanip>
#include <math.h>
#include <conio.h>
using namespace std;

const int n = 2;
const double eps1 = 1e-9;
const double eps2 = 1e-9;

void Newton_met(double ** Jac, double*F, double*x,int n, double* dx, int iter);
void M_Gaussa(double **Jac, double *F, int n, double *dx, double*x);
void Equation(double* F, double* x);
void Jacobi(double ** Jac, double *F, double* x,int n);

int main()
{
	int iter;
	double* F; //Вектор невязки
	F = new double[n];

	double* x; //вектор-столбец переменных
	x = new double[n];

	double** Jac; // матрица Якоби
	Jac = new double* [n];
	for (int i = 0; i < n; i++)
	{
		Jac[i] = new double[n+1];
	}

	double* dx; // вектор поправки
	dx = new double[n];

	cout << "Vvedite chislo iteraciy:" << endl;
	cin >> iter;
	cout << endl << "Vvedite " << n << " znacheniy:" << endl;
	for (int i = 0; i < n; i++)
	{
		cin >> x[i];
	}
	cout << "==================================================" << endl;
	Newton_met(Jac,F,x,n,dx,iter);
	cout << "==================================================" << endl;
	system("pause");
	return 0;
}

void Equation(double* F, double* x)
{
	F[0] = sin(x[0] + 1) - x[1] - 1;
	F[1] = 2 * x[0] + cos(x[1]) - 2;
}

void Jacobi(double** Jac, double* F, double* x, int n)
{
	double f1, f2;
	for (int i = 0; i < n;i++)
	{
		for (int j = 0; j < n; j++)
		{
			x[j] += eps1;
			Equation(F, x);
			f1 = F[i];
			x[j] -= eps1;
			Equation(F, x);
			f2 = F[i];
			Jac[i][j] = (f1 - f2) / eps1;
			Jac[i][n] = - F[i];
		}
	}
}

void M_Gaussa(double** Jac, double* F, int n, double* dx, double* x)
{
	//Прямой ход 
	for (int i = 0; i < n; i++)
	{
		double max = abs(Jac[i][i]);
		int my = i;
		for (int t = i; t < n; t++)
			if (abs(Jac[t][i]) > max)
			{
				max = abs(Jac[t][i]);
				my = t;
			}

		//Перемещение строк
		if (my != i)
		{
			double* per = Jac[i];
			Jac[i] = Jac[my];
			Jac[my] = per;
		}

		//деление строки 
		double amain = Jac[i][i];
		for (int z = 0; z < n + 1; z++)
		{
			Jac[i][z] = Jac[i][z] / amain;
		}

		//Вычитание из i=1 и деление строк i-y * на i-y коэффициент соотвествующей строки
		for (int j = i + 1; j < n; j++)
		{
			double b = Jac[j][i];
			for (int z = i; z < n + 1; z++)
				Jac[j][z] = Jac[j][z] - Jac[i][z] * b;
		}
	}

	//Обратный ход
	for (int i = n - 1; i > 0; i--)
	{
		for (int j = i - 1; j >= 0; j--)
			Jac[j][n] = Jac[j][n] - Jac[j][i] * Jac[i][n];
	}

	for (int i = 0; i < n; ++i)
		dx[i] = Jac[i][n];

	for (int i = 0; i < n; i++)
		x[i] += Jac[i][n];
}

void Newton_met(double** Jac, double* F, double* x, int n, double* dx, int iter)
{
	double D1, D2;
	double max;
	int k = 0;
	cout << "k_iter" << setw(12) << "del1" << setw(16) << "del2" << endl;
	cout << "==================================================" << endl;
	while (true)
	{
		Equation(F,x);
		Jacobi(Jac, F,x,n);

		//Вектор невязки
		for (int i = 0; i < n; i++)
		{
			F[i] *= -1;
		}
		M_Gaussa(Jac, F,n, dx, x);
		max = 0;

		Equation(F, x);
		for (int i = 0; i < n; i++)
		{
			if (abs(F[i]) > max)
			{
				max = abs(F[i]);
			}
		}
		D1 = max;
		max = 0;

		for (int i = 0; i < n; i++)
		{
			if (abs(x[i]) < 1 && abs(dx[i]) > max)
				max = abs(dx[i]);
			if (abs(x[i]) >= 1 && abs(dx[i] / x[i]) > max)
				max = abs(dx[i] / x[i]);
		}
		D2 = max;
		cout << setw(4)<<k + 1 << "\t" <<setw(10)<<D1 <<"\t"<< setw(10)<< D2 << endl;
		k++;
		
		if (D1 <= eps2 && D2 <= eps2 || k >= iter)
			break;
	}
	cout << "==================================================" << endl;
	cout << "X1 = " << x[0] << endl;
	cout << "X2 = " << x[1] << endl;
}
