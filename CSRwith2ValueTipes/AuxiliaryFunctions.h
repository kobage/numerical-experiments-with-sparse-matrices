#pragma once
//--------------------- string to numbers functions -------------
/*
* atof.c --
*
*	Source code for the "atof" library procedure.
*
* Copyright 1988 Regents of the University of California
* Permission to use, copy, modify, and distribute this
* software and its documentation for any purpose and without
* fee is hereby granted, provided that the above copyright
* notice appear in all copies.  The University of California
* makes no representations about the suitability of this
* software for any purpose.  It is provided "as is" without
* express or implied warranty.
*
* Modified at 2020 to detect the type (float or double) of
* the real number, which is being read from the file
*/

#include <stdlib.h>
#include <ctype.h>
#include"cfloat"
static int maxExponent = 511;	/* Largest possible base 10 exponent.  Any
								* exponent larger than this will already
								* produce underflow or overflow, so there's
								* no need to worry about additional digits.
								*/
static double powersOf10[] = {	/* Table giving binary powers of 10.  Entry */
	10.,			/* is 10^2^i.  Used to convert decimal */
	100.,			/* exponents into floating-point numbers. */
	1.0e4,
	1.0e8,
	1.0e16,
	1.0e32,
	1.0e64,
	1.0e128,
	1.0e256
};

/*
*----------------------------------------------------------------------
*
* atof --
*
*	This procedure converts a floating-point number from an ASCII
*	decimal representation to internal double-precision format.
*
* Results:
*	The return value is the floating-point equivalent of string.
*	If a terminating character is found before any floating-point
*	digits, then zero is returned.
*
* Side effects:
*	None.
*
* The function is modified at 09.2019
* After returning, the second argument indicates whether the returning number
* is float (isFloat = true) or double (isFloat = false).
* isFloat = true is sufficient condition (not necessary)
*----------------------------------------------------------------------
*/

//parameter bool & isFloat has been added at the modification

double my_atof(char*& string, bool& isFloat)
/* A decimal ASCII floating-point number,
* optionally preceded by white space.
* Must have form "-I.FE-X", where I is the
* integer part of the mantissa, F is the
* fractional part of the mantissa, and X
* is the exponent.  Either of the signs
* may be "+", "-", or omitted.  Either I
* or F may be omitted, or both.  The decimal
* point isn't necessary unless F is present.
* The "E" may actually be an "e".  E and X
* may both be omitted (but not just one).
*/
{
	bool sign, expSign = false;
	double fraction, dblExp, * d;
	char c;
	int exp = 0;		/* Exponent read from "EX" field. */
	int fracExp = 0;		/* Exponent that derives from the fractional
							* part.  Under normal circumstatnces, it is
							* the negative of the number of digits in F.
							* However, if I is very long, the last digits
							* of I get dropped (otherwise a long I with a
							* large negative exponent could cause an
							* unnecessary overflow on I alone).  In this
							* case, fracExp is incremented one for each
							* dropped digit.
							*/
	int mantSize;		/* Number of digits in mantissa. */
	int decPt;			/* Number of mantissa digits BEFORE decimal point. */
	char* pExp;			/* Temporarily holds location of exponent in string. */

						/*
						* Strip off leading blanks and check for a sign.
						*/

	char*& p = string;
	while (isspace(*p)) {
		p += 1;
	}
	if (*p == '-') {
		sign = true;
		p += 1;
	}
	else {
		if (*p == '+') {
			p += 1;
		}
		sign = false;
	}

	/*
	* Count the number of digits in the mantissa (including the decimal
	* point), and also locate the decimal point.
	*/

	decPt = -1;
	for (mantSize = 0; ; mantSize += 1)
	{
		c = *p;
		if (!isdigit(c))
		{
			if ((c != '.') || (decPt >= 0))
				break;
			decPt = mantSize;
		}
		p += 1;
	}

	//First step to distinguish float from double
	if (mantSize < 8) isFloat = true;
	else isFloat = false;

	/*
	* Now suck up the digits in the mantissa.  Use two integers to
	* collect 9 digits each (this is faster than using floating-point).
	* If the mantissa has more than 18 digits, ignore the extras, since
	* they can't affect the value anyway.
	*/

	pExp = p;
	p -= mantSize;
	if (decPt < 0) {
		decPt = mantSize;
	}
	else {
		mantSize -= 1;			/* One of the digits was the point. */
	}
	if (mantSize > 18) {
		fracExp = decPt - 18;
		mantSize = 18;
	}
	else {
		fracExp = decPt - mantSize;
	}
	if (mantSize == 0) {
		isFloat = true;					// it looks like as it is odd
		return 0.0;
	}
	else {
		int frac1, frac2;
		frac1 = 0;
		for (; mantSize > 9; mantSize -= 1)
		{
			c = *p;
			p += 1;
			if (c == '.') {
				c = *p;
				p += 1;
			}
			frac1 = 10 * frac1 + (c - '0');
		}
		frac2 = 0;
		for (; mantSize > 0; mantSize -= 1)
		{
			c = *p;
			p += 1;
			if (c == '.') {
				c = *p;
				p += 1;
			}
			frac2 = 10 * frac2 + (c - '0');
		}
		fraction = (1.0e9 * frac1) + frac2;
	}

	/*
	* Skim off the exponent.
	*/

	p = pExp;
	if ((*p == 'E') || (*p == 'e')) {
		p += 1;
		if (*p == '-') {
			expSign = true;
			p += 1;
		}
		else {
			if (*p == '+') {
				p += 1;
			}
			expSign = false;
		}
		while (isdigit(*p)) {
			exp = exp * 10 + (*p - '0');
			p += 1;
		}
	}
	if (expSign) {
		exp = fracExp - exp;
	}
	else {
		exp = fracExp + exp;
	}

	/*
	* Generate a floating-point number that represents the exponent.
	* Do this by processing the exponent one bit at a time to combine
	* many powers of 2 of 10. Then combine the exponent with the
	* fraction.
	*/

	if (exp < 0) {
		expSign = true;
		exp = -exp;
	}
	else {
		expSign = false;
	}
	if (exp > maxExponent) {
		exp = maxExponent;
	}
	dblExp = 1.0;
	for (d = powersOf10; exp != 0; exp >>= 1, d += 1) {
		if (exp & 01) {
			dblExp *= *d;
		}
	}
	if (expSign) {
		fraction /= dblExp;
	}
	else {
		fraction *= dblExp;
	}

	float mn = FLT_MIN;
	float mx = FLT_MAX;

	if (isFloat && fraction)
	{
		if (mn > fraction || fraction >= mx) isFloat = false;
	}

	if (sign) {
		return -fraction;
	}
	return fraction;
}

//restores index from the ASCII c-styoe string, pointer moves to the end of the string
inline size_t make_index(char*& str)
{
	size_t res = 0;
	while (isdigit(*str))
	{
		res = res * 10 + *str - '0';
		++str;
	}
	return res;
}

inline short myLog2(size_t num)
{
	short res = 0;
	for (short i = 32; i > 0; i >>= 1)
	{
		if ((1LL << (res + i)) <= num)
			res += i;
	}
	return res;
}