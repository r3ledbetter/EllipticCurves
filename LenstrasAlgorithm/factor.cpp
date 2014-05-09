#include <vector>
#include <time.h>
#include <utility>	// std::pair
#include <math.h>	// pow
#include <fstream>

#include "ecm.h"

using namespace std;

void factor(const mpreal& number);
void saveResults(int bound, int boundIndex);

int DIGITS = 5;
vector<int> primePowers;
vector<int> boundValues;
vector<mpreal> numbers;

// Store number of curves needed and length of time needed to factor
vector<pair<int, double> > numCurvesAndTime;

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

	for (int i = 0; i < 1000; ++i) {
		// Make sure to generate primes with desired number of digits
		mpreal factor1 = getPrime(DIGITS, 10);
		mpreal factor2 = getPrime(DIGITS, 10);
		while (factor1 < leastDigits) {
			factor1 = getPrime(DIGITS, 10);
		}
		while (factor2 < leastDigits) {
			factor2 = getPrime(DIGITS, 10);
		}
		numbers.push_back(factor1 * factor2);
	}

	// y = 10^(.1x)
	for (double x = 10; x <= 50; ++x) {
		boundValues.push_back(pow(10, 0.1 * x));
	}

	// Test
	for (int i = 0; i < boundValues.size(); ++i) {
		cout << "Beginning test with bound: " << boundValues[i] << endl;
		// Clear for fresh start with new bound
		primePowers.clear();
		numCurvesAndTime.clear();
		
		vector<int> sieveVector;
		for (int j = 0; j < boundValues[i]; ++j) {
			sieveVector.push_back(j);
		}
		sieveVector[1] = 0;
		// Sieve for prime numbers less than boundValues[i]
		for (int p = 2; p < boundValues[i]; ++p) {
			if (sieveVector[p] != 0) {
				int maxPower = p;
				for (maxPower = p; maxPower < boundValues[i]; maxPower *= p) {}
				primePowers.push_back(maxPower);
				for (int j = 2 * p; j < boundValues[i]; j += p) {
					sieveVector[j] = 0;
				}
			}
		}
		cout << "Time till completion:" << endl;
		cout << "          |" << endl;
		for (int j = 0; j < numbers.size(); ++j) {
			factor(numbers[j]);
			if (j % 100 == 0) {
				cout << "*" << flush;
			}
		}
		cout << endl;
		saveResults(boundValues[i], i);
	}
}

void saveResults(int bound, int boundIndex) {
	string path;
	stringstream tostring;
	tostring << boundIndex << "-" << bound;
	path = tostring.str();
	ofstream output;
	path = "../results/" + path + ".txt";
	output.open(path.c_str());
	for (int i = 0; i < numCurvesAndTime.size(); ++i) {
		output << numCurvesAndTime[i].first << " " << numCurvesAndTime[i].second << "\n";
	}
	output.close();
}

void factor(const mpreal& number) {
	clock_t t = clock();

	// Loop until finding a curve and point that "fails"
	for (int numCurves = 1; true; ++numCurves) {
		// Generate random point p = (a, b)
		Point p(DIGITS);
		EllipticCurve E(p, number);

		for (int i = 0; i < primePowers.size(); ++i) {
			p = E.doubleAndAdd(primePowers[i], p);
			// If true, factor was found
			if (p.done) {
				t = clock() - t;
				pair<int, double> record;
				record.first = numCurves;
				record.second = ((double)t) / CLOCKS_PER_SEC;
				numCurvesAndTime.push_back(record);
				return;
			}
		}
	}
}
