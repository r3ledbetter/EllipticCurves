#include "ecm.h"

int DIGITS = 7;

int main() {
	srand(time(NULL));
	mpfr::random(time(NULL));

	// Needs to be at least 3x as big for generating curve around point
	int maxDigits = 3 * DIGITS + 1;
	mpreal::set_default_prec(mpfr::digits2bits(maxDigits));
	cout.precision(3 * DIGITS);

	int leastDigits = 1;
	for (int i = 1; i < DIGITS; ++i) {
		leastDigits *= 10;
	}

	mpreal factor1 = getPrime(DIGITS, 10);
	mpreal factor2 = getPrime(DIGITS, 10);
	while (factor1 < leastDigits) {
		factor1 = getPrime(DIGITS, 10);
	}
	while (factor2 < leastDigits) {
		factor2 = getPrime(DIGITS, 10);
	}

	mpreal number = factor1 * factor2;







	/*
	// Loop until finding a curve and point that "fails"
	while (true) {
		// Generate random point p = (a, b)
		Point p(DIGITS);
		cout << number << endl;
		EllipticCurve E(p, number);

		for (int i = 2; i < 1000; ++i) {
			//p1.print();
			p = E.doubleAndAdd(i, p);
			if (p.infinity) {
				cout << "darn" << endl;
				getchar();
			}
		}
	}
	*/
}
