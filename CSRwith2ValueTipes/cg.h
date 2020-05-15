#pragma once
#include <string>
#include"chrono"
#include <numeric>
#include <vector>
#include <algorithm>
#include <iomanip>

using namespace std;

namespace cg
{
	size_t n;
	char repNumber{ 1 };
	double alpha;
	double* x, * r, * d, * q;

	//functions
	void setParameters(void);
	void destructParameters(void);
	void printSolution(void) noexcept;

	template<typename T>
	void cgSymmPosDef(unique_ptr<T>&&, double) noexcept;

	template<typename T>
	void cgSymmPosDefPrecond(unique_ptr<T>&&, double) noexcept;

	template<typename T>
	void SteepestDescent(unique_ptr<T>&&, double) noexcept;

	template<typename T>
	void cgDanielSymmPosDef(unique_ptr<T>&&, double) noexcept;

	template<typename Solver, typename T>
	void multiSolve(Solver, unique_ptr<T>&& m, double tolerance) noexcept;
};
void cg::printSolution(void)  noexcept
{
	int quarter = n / 4;
	for (int i = 0; i < quarter; i++)
	{
		printf("x[%3d] = %-15.11f     x[%3d] = %15.11f     x[%3d] = %15.11f     x[%3d] = %2.11f  \n",
			i, x[i], quarter + i, x[quarter + i], 2 * quarter + i,
			x[2 * quarter + i], 3 * quarter + i, x[3 * quarter + i]);
	}
	for (int i = 4 * quarter; i < n; i++)
		printf("%87sx[%3d] = %-15.11f     \n", "", i, x[i]);
}

void cg::setParameters(void)
{
	x = (double*)malloc(n * sizeof(double));
	d = (double*)malloc(n * sizeof(double));
	r = (double*)malloc(n * sizeof(double));
	q = (double*)malloc(n * sizeof(double));
}

void cg::destructParameters(void)
{
	free(x);
	free(r);
	free(d);
	free(q);
}

template<typename T>
void cg::cgSymmPosDef(unique_ptr<T>&& m, double tolerance) noexcept
{
	for (int i = 0; i < n; i++)
	{
		d[i] = r[i] = 1;
		x[i] = 0.;
	}
	double deltaNew(inner_product(r, r + n, r, 0.));

	while (deltaNew >= tolerance * tolerance)
	{
		m->multiplyAndSave(d, q);

		double alpha = deltaNew / inner_product(d, d + n, q, 0.);
		for (int i = 0; i < n; i++)
			x[i] += alpha * d[i];

		for (int i = 0; i < n; i++)
			r[i] -= alpha * q[i];

		double deltaOld(deltaNew);
		deltaNew = inner_product(r, r + n, r, 0.);
		double beta(deltaNew / deltaOld);
		for (int i = 0; i < n; i++)
			d[i] = r[i] + (beta * d[i]);
	}
}

//The simplest preconditioner is used
template<typename T>
void cg::cgSymmPosDefPrecond(unique_ptr<T>&& m, double tolerance) noexcept
{
	double* s = (double*)malloc(n * sizeof(double));
	double* M = (double*)malloc(n * sizeof(double));

	m->Diagonal(M);

	for (int i = 0; i < n; i++)
	{
		r[i] = 1;
		d[i] = r[i] / M[i];
		x[i] = 0.;
	}
	double deltaNew(inner_product(r, r + n, d, 0.));

	while (inner_product(r, r + n, r, 0.) >= tolerance * tolerance)
	{
		m->multiplyAndSave(d, q);

		double alpha = deltaNew / inner_product(d, d + n, q, 0.);
		for (size_t i = 0; i < n; i++)
			x[i] += alpha * d[i];

		for (size_t i = 0; i < n; i++)
			r[i] -= alpha * q[i];

		for (size_t i = 0; i < n; i++)
			s[i] = r[i] / M[i];

		double deltaOld(deltaNew);
		deltaNew = inner_product(r, r + n, s, 0.);
		double beta(deltaNew / deltaOld);
		for (int i = 0; i < n; i++)
			d[i] = s[i] + beta * d[i];
	}
	free(s);
	free(M);
}

template<typename Solver, typename T>
void cg::multiSolve(Solver solve, unique_ptr<T>&& m, double tolerance) noexcept
{
	n = m->n;
	setParameters();
	vector<_int64> repetitions(repNumber);
	//Create expensive objects	
	auto st = chrono::high_resolution_clock::now();
	auto diff = chrono::high_resolution_clock::now() - st;
	auto time = chrono::duration_cast<chrono::microseconds>(diff);

	for (int j = 0; j < repNumber; j++)
	{
		st = chrono::high_resolution_clock::now();

		solve(move(m), tolerance);

		diff = chrono::high_resolution_clock::now() - st;
		time = chrono::duration_cast<chrono::microseconds>(diff);
		repetitions[j] = time.count();
	}
	sort(repetitions.begin(), repetitions.end());
	freopen("results.txt", "w", stdout);
	printSolution();
	cout << "Compute time: " << repetitions[repNumber / 2] << endl;

	destructParameters();
}

template<typename T>
void cg::SteepestDescent(unique_ptr<T>&& m, double tolerance) noexcept
{
	for (int i = 0; i < n; i++)
	{
		r[i] = 1;
		x[i] = 0.;
	}
	double deltaNew(inner_product(r, r + n, r, 0.));

	while (deltaNew >= tolerance * tolerance)
	{
		m->multiplyAndSave(r, q);

		double alpha = deltaNew / inner_product(r, r + n, q, 0.);


		for (int i = 0; i < n; i++)
			x[i] += alpha * r[i];

		for (int i = 0; i < n; i++)
			r[i] -= alpha * q[i];

		deltaNew = inner_product(r, r + n, r, 0.);
	}
}

template<typename T>
void cg::cgDanielSymmPosDef(unique_ptr<T>&& m, double tolerance) noexcept
{
	for (int i = 0; i < n; i++)
	{
		d[i] = r[i] = 1;
		x[i] = 0.;
	}
	double deltaNew(inner_product(r, r + n, r, 0.));
	for (int i = 0; i < n; i++)
		m->multiplyAndSave(d, q);
	double alpha = deltaNew / inner_product(d, d + n, q, 0.);
	for (int i = 0; i < n; i++)
		x[i] += alpha * d[i];
	for (int i = 0; i < n; i++)
		r[i] -= alpha * q[i];
	deltaNew = inner_product(r, r + n, r, 0.);

	while (deltaNew >= tolerance * tolerance)
	{
		double beta(inner_product(q, q + n, r, 0.) / inner_product(d, d + n, q, 0.));
		for (int i = 0; i < n; i++)
			d[i] = r[i] - beta * d[i];

		m->multiplyAndSave(d, q);
		alpha = deltaNew / inner_product(d, d + n, q, 0.);
		for (int i = 0; i < n; i++)
			x[i] += alpha * d[i];

		for (int i = 0; i < n; i++)
			r[i] -= alpha * q[i];
		deltaNew = inner_product(r, r + n, r, 0.);
	}
}