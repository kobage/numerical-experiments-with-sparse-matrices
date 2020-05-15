#pragma once
#include <iostream>
#include <string>
#include<random>
#include"AuxiliaryFunctions.h"
using namespace std;

class simpleAbstractSM
{
public:
	uint32_t n{};
	uint32_t nnzz{};
	virtual ~simpleAbstractSM() {};

	virtual void multiplyAndSave(double*, double*) = 0;
	virtual void Diagonal(double*) = 0;
};

// Straightforward implementation, aimed on the conjugate gradient method
class CSR : public simpleAbstractSM
{
public:
	double* AA;
	uint32_t* JA;
	uint32_t* IA;
	CSR(string);
	CSR(char*& data);
	~CSR();

	virtual void multiplyAndSave(double*, double*);
	virtual void Diagonal(double*);
};

CSR::CSR(string s)
{
	/* open an existing file for reading */
	FILE* infile = fopen(s.c_str(), "r");

	/* declare a file pointer */
	char* buffer;
	long numbytes;

	/* if the file does not exist */
	if (infile == NULL)
		cout << "the file does not exist!" << endl;

	/* Get the number of bytes */
	fseek(infile, 0L, SEEK_END);
	numbytes = ftell(infile);

	/* reset the file position indicator to
	the beginning of the file */
	fseek(infile, 0L, SEEK_SET);

	/* grab sufficient memory for the
	buffer to hold the text */
	buffer = (char*)malloc(numbytes * sizeof(char));

	/* memory error */
	if (buffer == NULL)
		cout << "memory error!" << endl;

	/* copy all the text into the buffer */
	fread(buffer, sizeof(char), numbytes, infile);
	fclose(infile);

	//read quantities
	char* data = buffer;
	n = make_index(data);
	++data;
	nnzz = make_index(data);

	//allocate memory
	AA = (double*)malloc(sizeof(double) * nnzz);
	JA = (uint32_t*)malloc(sizeof(uint32_t) * nnzz);
	IA = (uint32_t*)malloc(sizeof(uint32_t) * (n + 1));
	IA[0] = 0;

	size_t row = 1;
	size_t iaInd = 1;
	size_t c = 0;		//counter for rows

	//temporary storage for the strings from the file
	double v; 	uint32_t i, j;

	bool dummy{ false };
	for (uint32_t ii = 0; ii < nnzz; ++ii)
	{
		++data;
		//read with the predefined order: j,i,v
		//j:
		j = make_index(data);
		++data;
		JA[ii] = j - 1;
		//i:
		i = make_index(data);
		++data;
		//v
		v = my_atof(data, dummy);

		AA[ii] = v;
		if (row != i)
		{
			IA[iaInd] = IA[iaInd - 1] + c;
			++iaInd;
			row = i;
			c = 0;
		}
		++c;
	}
	IA[n] = IA[n - 1] + c;
	free(buffer);
}
CSR::CSR(char*& data)
{
	//quantities
	n = make_index(data);
	++data;
	nnzz = make_index(data);

	//allocate memory
	AA = (double*)malloc(sizeof(double) * nnzz);
	JA = (uint32_t*)malloc(sizeof(uint32_t) * nnzz);
	IA = (uint32_t*)malloc(sizeof(uint32_t) * (n + 1));
	IA[0] = 0;

	size_t row = 1;
	size_t iaInd = 1;
	size_t c = 0;		//counter for rows

	//temporary storage for the strings from the file
	double v; 	uint32_t i, j;

	bool dummy{ false };
	for (uint32_t ii = 0; ii < nnzz; ++ii)
	{
		++data;
		//read with the predefined order: j,i,v
		//j:
		j = make_index(data);
		++data;
		JA[ii] = j - 1;
		//i:
		i = make_index(data);
		++data;
		//v
		v = my_atof(data, dummy);

		AA[ii] = v;
		if (row != i)
		{
			IA[iaInd] = IA[iaInd - 1] + c;
			++iaInd;
			row = i;
			c = 0;
		}
		++c;
	}
	IA[n] = IA[n - 1] + c;
}
CSR::~CSR(void)
{
	free(AA);
	free(IA);
	free(JA);
	cout << "CSR matrix is freed out" << endl;
}

//multiplication: sparse symmetric matrix  - vector 
void CSR::multiplyAndSave(double* a, double* res)
{
	size_t i, j;
	size_t left, right;

	for (i = 0; i < n; ++i)
		res[i] = 0.;

	for (i = 0; i < n; ++i)
	{
		left = IA[i];
		right = IA[i + 1];

		res[i] += AA[left] * a[JA[left]];

		for (j = left + 1; j < right; ++j)
		{
			res[i] += AA[j] * a[JA[j]];
			res[JA[j]] += AA[j] * a[i];
		}
	}
}

void CSR::Diagonal(double* dg)
{
	for (size_t i = 0; i < n; ++i)
	{
		dg[i] = AA[IA[i]];
	}
}

//---------------------------------- low-memory CSR: lm_CSR  -------------------------------
template<typename T>
class lm_CSR : public simpleAbstractSM
{
public:
	double* AAd;
	float* AAs;
	T* IA;
	uint32_t* JAd;
	uint32_t* JAs; 
	lm_CSR(char*&, int32_t, int32_t);
	~lm_CSR();
	virtual void multiplyAndSave(double*, double*);
	virtual void Diagonal(double*);
};

template<typename T>
lm_CSR<T>::lm_CSR(char*& data, int32_t cd, int32_t cs)
{
	//read quantities
	n = make_index(data);
	++data;
	nnzz = make_index(data);

	//allocate memory
	AAd = (double*)malloc(sizeof(double) * cd);
	AAs = (float*)malloc(sizeof(float) * cs);
	JAd = (uint32_t*)malloc(sizeof(uint32_t) * cd);
	JAs = (uint32_t*)malloc(sizeof(uint32_t) * cs);
	IA = (T*)malloc(sizeof(T) * n);

	//fill = reenter
	T rd = 0;		//quantities of s -(single or float) and d (doubles)
	T rs = 0;
	cd = 0;   //indexes of the variables of the corresponding types
	cs = 0;   //
	size_t row = 1;
	double v; 	uint32_t i, j;

	//fill
	bool isFloat;
	uint8_t bits = (8 * sizeof(T)) / 2;			// T bits = 

	for (int ii = 0; ii < nnzz; ++ii)
	{
		++data;
		//read with the predefined order: j,i,v
		//j:
		j = make_index(data);
		++data;
		//i:
		i = make_index(data);
		++data;
		//v
		isFloat = false;
		v = my_atof(data, isFloat);
		if (row != i)
		{
			IA[row - 1] = (rd << bits) + rs;
			row = i;
			rs = 0;
			rd = 0;
		}
		if (isFloat)
		{
			AAs[cs] = v;
			JAs[cs] = j - 1;
			++cs;
			++rs;
		}
		else
		{
			AAd[cd] = v;
			JAd[cd] = j - 1;
			++cd;
			++rd;
		}
	}
	IA[n - 1] = (rd << bits) + rs;

}
template<typename T>
lm_CSR<T>::~lm_CSR(void)
{
	free(AAd);
	free(AAs);
	free(IA);
	free(JAd);
	free(JAs);
	cout << "lm_CSR<> matrix is freed out" << endl;
}

//multiplication: sparse symmetric matrix  - vector 
template<typename T>
void lm_CSR<T>::multiplyAndSave(double* x, double* res)
{
	size_t i, j;

	for (i = 0; i < n; ++i)
		res[i] = 0.;

	uint32_t leftS, rightS{ 0 };
	uint32_t leftD, rightD{ 0 };

	uint8_t bits = (8 * sizeof(T)) / 2;
	T tmp = (2 << (bits - 1));

	for (i = 0; i < n; ++i)
	{
		leftS = rightS;
		leftD = rightD;
		rightS += IA[i] % tmp;
		rightD += IA[i] >> bits;

		if (rightS > 0 && JAs[leftS] == i)
		{
			res[i] += AAs[leftS] * x[JAs[leftS]];
			++leftS;
		}
		else
		{
			res[i] += AAd[leftD] * x[JAd[leftD]];
			++leftD;
		}
		for (j = leftS; j < rightS; ++j)
		{
			res[i] += AAs[j] * x[JAs[j]];
			res[JAs[j]] += AAs[j] * x[i];
		}
		for (j = leftD; j < rightD; ++j)
		{
			res[i] += AAd[j] * x[JAd[j]];
			res[JAd[j]] += AAd[j] * x[i];
		}
	}
}

template<typename T>
void lm_CSR<T>::Diagonal(double* dg)
{
	size_t i, j;
	uint32_t leftS, rightS{ 0 };
	uint32_t leftD, rightD{ 0 };

	uint8_t bits = (8 * sizeof(T)) / 2;
	T tmp = (2 << (bits - 1));

	for (i = 0; i < n; ++i)
	{
		leftS = rightS;
		leftD = rightD;
		rightS += IA[i] % tmp;
		rightD += IA[i] >> bits;

		if (rightS > leftS && JAs[leftS] == i)
		{
			dg[i] = AAs[leftS];
			++leftS;
		}
		else
		{
			dg[i] = AAd[leftD];
			++leftD;
		}
	}
}

//----low-memory CSR: lm_CSR_co - compromise between speed and low memory -----
template<typename T>
class lm_CSR_co : public simpleAbstractSM
{
public:
	double* ADiagonal;
	double* AAd;
	float* AAs;
	T* IA;
	uint32_t* JAd;
	uint32_t* JAs;
	lm_CSR_co(char*&, int32_t, int32_t);
	~lm_CSR_co();
	virtual void multiplyAndSave(double*, double*);
	virtual void Diagonal(double*);
};

template<typename T>
lm_CSR_co<T>::lm_CSR_co(char*& data, int32_t cd, int32_t cs)
{
	//read quantities
	n = make_index(data);
	++data;
	nnzz = make_index(data);

	//allocate memory
	ADiagonal = (double*)malloc(sizeof(double) * n);
	AAd = (double*)malloc(sizeof(double) * cd);
	AAs = (float*)malloc(sizeof(float) * cs);
	JAd = (uint32_t*)malloc(sizeof(uint32_t) * cd);
	JAs = (uint32_t*)malloc(sizeof(uint32_t) * cs);
	IA = (T*)malloc(sizeof(T) * n);

	//fill = re enter
	T rd = 0;		//quantities of s -(single or float) and d (doubles)
	T rs = 0;
	cd = 0;   //შესაბამისი ტიპის ცვლადების ინდექსებია
	cs = 0;   //
	size_t row = 1;
	double v; 	uint32_t i, j;

	//fill matrix
	bool isFloat;
	uint8_t bits = (8 * sizeof(T)) / 2;

	++data;
	//read with the predefined order: j,i,v
	//j:
	j = make_index(data);
	++data;
	//i:
	i = make_index(data);
	++data;
	//v
	isFloat = false;
	v = my_atof(data, isFloat);
	ADiagonal[0] = v;

	for (int ii = 1; ii < nnzz; ++ii)
	{
		++data;
		//read with the predefined order: j,i,v
		//j:
		j = make_index(data);
		++data;
		//i:
		i = make_index(data);
		++data;
		//v
		isFloat = false;
		v = my_atof(data, isFloat);
		if (row != i)
		{
			IA[row - 1] = (rd << bits) + rs;
			ADiagonal[row] = v;
			row = i;
			rs = 0;
			rd = 0;
			continue;
		}
		if (isFloat)
		{
			AAs[cs] = v;
			JAs[cs] = j - 1;
			++cs;
			++rs;
		}
		else
		{
			AAd[cd] = v;
			JAd[cd] = j - 1;
			++cd;
			++rd;
		}
	}
	IA[n - 1] = (rd << bits) + rs;
}
template<typename T>
lm_CSR_co<T>::~lm_CSR_co(void)
{
	free(ADiagonal);
	free(AAd);
	free(AAs);
	free(IA);
	free(JAd);
	free(JAs);
	cout << "lm_CSR<> matrix is freed out" << endl;
}

//multiplication: sparse symmetric matrix  - vector 
template<typename T>
void lm_CSR_co<T>::multiplyAndSave(double* x, double* res)
{
	size_t i, j;

	for (i = 0; i < n; ++i)
		res[i] = ADiagonal[i] * x[i];

	uint32_t leftS, rightS{ 0 };
	uint32_t leftD, rightD{ 0 };

	uint8_t bits = (8 * sizeof(T)) / 2;
	T tmp = (2 << (bits - 1));
	for (i = 0; i < n; ++i)
	{
		leftS = rightS;
		leftD = rightD;
		rightS += IA[i] % tmp;
		rightD += IA[i] >> bits;
		for (j = leftS; j < rightS; ++j)
		{
			res[i] += AAs[j] * x[JAs[j]];
			res[JAs[j]] += AAs[j] * x[i];
		}
		for (j = leftD; j < rightD; ++j)
		{
			res[i] += AAd[j] * x[JAd[j]];
			res[JAd[j]] += AAd[j] * x[i];
		}
	}
}

template<typename T>
void lm_CSR_co<T>::Diagonal(double* dg)
{
	size_t i;
	for (i = 0; i < n; ++i)
		dg[i] = ADiagonal[i];
}

//----------------------------  GenericCSR  ---------------------------
class GenericCSR : public simpleAbstractSM
{
public:
	unique_ptr<simpleAbstractSM> mp;
	GenericCSR(string s)
	{
		/* open an existing file for reading */
		FILE* infile = fopen(s.c_str(), "r");

		/* declare a file pointer */
		char* buffer;
		long numbytes;

		/* if the file does not exist */
		if (infile == NULL)
			cout << "the file does not exist!" << endl;

		/* Get the number of bytes */
		fseek(infile, 0L, SEEK_END);
		numbytes = ftell(infile);

		/* reset the file position indicator to
		the beginning of the file */
		fseek(infile, 0L, SEEK_SET);

		/* grab sufficient memory for the
		buffer to hold the text */
		buffer = (char*)malloc(numbytes * sizeof(char));

		/* memory error */
		if (buffer == NULL)
			cout << "memory error!" << endl;

		/* copy all the text into the buffer */
		fread(buffer, sizeof(char), numbytes, infile);
		fclose(infile);

		//read quantities
		char* data = buffer;
		char* p = data;

		n = make_index(data);
		++data;
		nnzz = make_index(data);

		//the re enter point
		uint32_t cd = 0;		//counter for doubles
		uint32_t cs = 0;		//counter for singles

		uint16_t rd = 0;		//quantities of s -(single or float)
		uint16_t rs = 0;		// and d (doubles) in the row

		size_t row = 1;
		//temporary storage for the strings from the file
		double v; 	uint32_t i, j;

		//for the lengthiest row
		uint16_t maxLength{ 0 };

		//starting count of singles and doubles 
		//find the lengthiest row
		bool isFloat;

		++data;
		//read: j,i,v
		//j: - jump ----------------- 
		while (*data != ' ')  ++data;
		++data;
		//i:- jump -----------------  
		while (*data != ' ')  ++data;
		++data;
		//v -  
		isFloat = false;
		v = my_atof(data, isFloat);
		
		for (int ii = 1; ii < nnzz; ++ii)
		{
			++data;
			//read with the predefined order: j,i,v
			//j: - jump ----------------- 
			while (*data != ' ')  ++data;
			++data;
			//i:
			i = make_index(data);
			++data;
			//v -  
			isFloat = false;
			v = my_atof(data, isFloat);

			if (row != i)
			{
				row = i;
				rs = 0;
				rd = 0;
				continue;
			}

			if (isFloat)
			{
				++cs;
				++rs;
			}
			else
			{
				++cd;
				++rd;
			}
			if (maxLength < rs) maxLength = rs;
			if (maxLength < rd) maxLength = rd;
		}

		if (maxLength < 16)
		{
			mp = make_unique<lm_CSR_co<uint8_t>>(p, cd, cs);
			return;
		}
		if (maxLength < 255)
		{
			mp = make_unique<lm_CSR_co<uint16_t>>(p, cd, cs);
			return;
		}
		mp = make_unique<lm_CSR_co<uint32_t>>(p, cd, cs);
		free(buffer);
	}


	~GenericCSR() {}
	virtual void multiplyAndSave(double* x, double* res)
	{
		mp->multiplyAndSave(x, res);
	}

	virtual void Diagonal(double* res)
	{
		mp->Diagonal(res);
	}
};