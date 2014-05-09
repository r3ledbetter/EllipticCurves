#include "../GeneratePrime.h"

#include <time.h>

using mpfr::mpreal;
using mpfr::floor;
using mpfr::mod;
using mpfr::abs;
using std::vector;
using std::cout;
using std::endl;

// Solve the congruence a = bx (mod n)
mpreal solveCongruence(const mpreal& a, const mpreal& b, const mpreal& n, bool& factorFound);

// On the elliptic curve
struct Point {
	mpreal X;
	mpreal Y;
	bool infinity;
	bool done;	// used to pass information when done factoring

	Point() : infinity(false), done(false) {}
	Point(bool i) : infinity(i), done(false) {}
	Point(int digits) {
		X = getLargeInt(digits);
		Y = getLargeInt(digits);
		infinity = false;
		done = false;
	}
	Point(mpreal x, mpreal y) : X(x), Y(y), infinity(false), done(false) {}
	Point(mpreal x, mpreal y, bool i) : X(x), Y(y), infinity(i), done(false) {}

	void print() {
		if (infinity) {
			std::cout << "Infinity" << std::endl;
		}
		else {
			std::cout << "(" << X << ", " << Y << ")" << std::endl;
		}
	}
};

// Y^2 = X^3 + AX + B
class EllipticCurve {
	mpreal A;
	mpreal B;
	mpreal N;

public:
	EllipticCurve(const mpreal& a, const mpreal& b, const mpreal& n)
		: A(a), B(b), N(n) {}
	EllipticCurve(const Point& p, const mpreal& n);

	Point add(const Point& p1, const Point& p2);

	Point doubleAndAdd(int multiple, const Point& p);

	void print();

};

Point EllipticCurve::doubleAndAdd(int multiple, const Point& p) {
	Point product(true);	// Infinity point
	Point Q(p.X, p.Y);
	while (multiple > 0) {
		if (mod(multiple, 2) == 1) {
			product = add(product, Q);
			if (product.done) {
				return product;
			}
			--multiple;
		}
		else {
			Q = add(Q, Q);
			if (Q.done) {
				return Q;
			}
			multiple = multiple / 2;	// floor -> int
		}
	}
	return product;
}

EllipticCurve::EllipticCurve(const Point& p, const mpreal& n) {
	// Generate elliptic curve containing the point p
	N = n;
	A = getLargeInt(2);	// How many digits for "random" A??
	B = (p.Y * p.Y) - (p.X * p.X * p.X) - (A * p.X);
	B = mod(B, N);
}

Point EllipticCurve::add(const Point& p1, const Point& p2) {
	Point sum;
	mpreal lambda, numerator, denominator;
	bool factorFound = false;

	// Consider reordering to skip several checks..
	if (p1.X == p2.X && p1.Y == -p2.Y) {
		sum.infinity = true;
		return sum;
	}
	else if (p1.infinity) {
		return p2;
	}
	else if (p2.infinity) {
		return p1;
	}
	else if (p1.X == p2.X && p1.Y == p2.Y) {
		numerator = ((3 * p1.X * p1.X) + A);
		denominator = (2* p1.Y);
	}
	else {
		numerator = (p2.Y - p1.Y);
		denominator = (p2.X - p1.X);
	}
	// Working mod N
	numerator = mod(numerator, N);
	denominator = mod(denominator, N);

	// Need to use extended Euclidean algorithm to compute lambda
	lambda = solveCongruence(numerator, denominator, N, factorFound);
	if (factorFound) {
		sum.done = true;
		return sum;
	}

	sum.X = (lambda * lambda) - p1.X - p2.X;
	sum.X = mod(sum.X, N);
	sum.Y = lambda * (p1.X - sum.X) - p1.Y;
	sum.Y = mod(sum.Y, N);
	return sum;
}

void EllipticCurve::print() {
	std::cout << "Y^2 = X^3 + " << A << "X + " << B << "  (mod " << N << ")" << std::endl;
}

// Solve the congruence a = bx (mod n)
mpreal solveCongruence(const mpreal& a, const mpreal& b, const mpreal& n, bool& factorFound) {
	// Values to iteratively solve Euclidean algorithm
	// i.e. find A = gcd(n, b) and x2, y2 such that n*x2 + b*y2 = A
	mpreal A = n, B = b, q, r, xTemp, yTemp;
	mpreal x1 = 0;
	mpreal x2 = 1;
	mpreal y1 = 1;
	mpreal y2 = 0;
	while (B != 0) {
		q = floor(A / B);
		r = A - (q * B);
		A = B;
		B = r;
		xTemp = x2 - (q * x1);
		yTemp = y2 - (q * y1);
		x2 = x1;
		x1 = xTemp;
		y2 = y1;
		y1 = yTemp;
	}
	if (A != 1) {
		//if (A == n) {
		//	cout << "Found n as a factor of n :(" << endl;
		//}
		factorFound = true;
		return A;
	}
	// x2 and y2 now contain values such that n*x2 + b*y2 = A


	mpreal inverse;
	if (abs(x2) > abs(y2)) {
		// x2 is the inverse of b (mod n)
		inverse = x2;
	}
	else if (abs(y2) > abs(x2)) {
		// y2 is the inverse of b (mod n)
		inverse = y2;
	}

	// Working mod n
	if (inverse < 0) {
		inverse += n;
	}
	return mod(inverse * a, n);
}
