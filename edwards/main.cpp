#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include "main.h"

using namespace std;

#define PI 3.1415926535

double function(double x)
{
	//return (1.0 / (pow(x, 2.0) * pow(-0.126 * x + 0.109, 2.0) * 2.46038 * 2.46038));
	return (1.0 / pow(1.35,3.0)/(pow(x, 2.0) * (0.31829 - 0.90469 * x + 0.64821 * x * x)));

}

double integral(double(*fun)(double x), double a, double b, int n)
{
	double sum, step, result;
	int i;
	sum = (fun(a) + fun(b)) / 2.0;
	step = (b - a) / n;
	for (i = 1; i < n; i++)
	{
		sum += fun(a + i * step);
	}
	result = sum * step;
	return result;
}

void main()
{
	double result;

	ofstream newfile("ellipsoid.txt");//write files
	newfile << fixed << setprecision(17);
	newfile << "PD" << "\t" << "1/X" << "\t" << "X" << "\t" << "PD" << endl;

	for(int i = 0; i < 50; i++)
	{
		result = integral(function, 0.584, 0.58648 + 0.00248 * i, 200*(50-i));
		newfile << 0.58648 + 0.00248 * i << "\t" << result << "\t" << 1.0 / result << "\t" << 0.58648 + 0.00248 * i << endl;
	}
	newfile.close();

	double a = tgamma(0.5 / 1.0);
	double b = tgamma(1.5 / 1.0);
	double vol = 2.0 * a * a * a / b / 3.0;
	vol *= (2.20 * 3.75 * 3.75);

	//double resvol = 4.0 / 3.0 * PI * pow(3.0 * 1.35, 3.0);
	double resvol = pow(1.35, 3.0);

	double ratio = resvol;
	cout << ratio << endl;

	result = integral(function, 0.584, 0.708, 10000);
	//result = integral(function, 0.587, 0.64, 10000);
	cout << result << endl;

	cout << "press any key to continue!\n";
	getchar();
}