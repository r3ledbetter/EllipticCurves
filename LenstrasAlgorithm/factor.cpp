#include <vector>
#include <time.h>

#include "ecm.h"

int DIGITS = 4;
vector<int> primePowers;

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

	int bound = 100000;
	vector<int> sieveVector;
	for (int i = 0; i < bound; ++i) {
		sieveVector.push_back(i);
	}
	sieveVector[1] = 0;
	// Sieve for prime numbers less than bound
	for (int p = 2; p < bound; ++p) {
		if (sieveVector[p] != 0) {
			int maxPower = p;
			for (maxPower = p; maxPower < bound; maxPower *= p) {}
			primePowers.push_back(maxPower);
			for (int i = 2 * p; i < bound; i += p) {
				sieveVector[i] = 0;
			}
		}
	}

	cout << "Start time: " << double(clock()) / CLOCKS_PER_SEC << endl;

	// Loop until finding a curve and point that "fails"
	while (true) {
		// Generate random point p = (a, b)
		Point p(DIGITS);
		cout << number << endl;
		EllipticCurve E(p, number);

		for (int i = 0; i < primePowers.size(); ++i) {
			//p1.print();
			p = E.doubleAndAdd(primePowers[i], p);
		}
	}
}
