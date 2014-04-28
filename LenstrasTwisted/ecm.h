#include <vector>

#include "../GeneratePrime/GeneratePrime.cpp"
#include "../mpreal.h"

using mpfr::mpreal;
using mpfr::floor;
using mpfr::mod;
using mpfr::abs;
using std::vector;
using std::cout;
using std::endl;

// Solve the congruence a = bx (mod n)
mpreal solveCongruence(const mpreal& a, const mpreal& b, const mpreal& n);

// On the elliptic curve
struct Point {
	mpreal X;
	mpreal Y;
	mpreal Z;
	mpreal T;

	Point();
	// Random point
	Point(int digits) {
		X = getLargeInt(digits);
		Y = getLargeInt(digits);
		Z = 1;
		T = X * Y;
	}
	Point(mpreal x, mpreal y) : X(x), Y(y), Z(1), T(X * Y) {}

	void print() {
		cout << "(" << X << " : " << Y << " : "
			<< Z << " : " << T << ")" << endl;
	}

	bool operator==(const Point& p) {
		if (p.X == X && p.Y == Y &&
			p.Z == Z && p.T == T) {
			return true;
		}
		return false;
	}
};

/*** EXTENDED TWISTED EDWARDS CURVE ***/
// {(X : Y : Z : T) Element of P^3
// (x, y) -> (x : y : 1 : xy)
// aX^2 + Y^2 = Z^2 + dT^2 and XY = ZT (choose a = 1)
class EllipticCurve {
	vector<Point> knownPoints;
	mpreal a;
	mpreal d;
	mpreal N;

public:
	EllipticCurve(const mpreal& d, const mpreal& n) : d(d), N(n) {}
	EllipticCurve(const Point& p, const mpreal& n);

	Point add(const Point& p1, const Point& p2);

	Point doubleAndAdd(int multiple, const Point& p);

	void print();

};

Point EllipticCurve::doubleAndAdd(int multiple, const Point& p) {
	Point product(true);	// Infinity point
	Point Q(p.X, p.Y);
	while (multiple > 0) {
		if (multiple % 2 == 1) {
			product = add(product, Q);
		}
		Q = add(Q, Q);
		multiple = multiple / 2;	// floor -> int
	}
	return product;
}

Point EllipticCurve::add(const Point& p1, const Point& p2) {
	Point p3;	// sum

	if (p1.X == p2.X && p1.Y == p2.Y && p1.Z == p2.Z && p1.T == p2.T) {
		mpreal B = mod((p1.X + p1.Y) * (p1.X + p1.Y), N);
		mpreal C = mod(p1.X * p1.X, N);
		mpreal D = mod(p1.Y * p1.Y, N);
		mpreal E = mod(a * C, N);
		mpreal F = mod(E + D, N);
		mpreal H = mod(p1.Z * p1.Z, N);
		mpreal J = mod(F - (2 * H), N);
		p3.X = mod((B - C - D) * J, N);
		p3.Y = mod(F * (E - D), N);
		p3.Z = mod(F * J, N);
	}
	else {
		mpreal A = mod(p1.X * p2.X, N);
		mpreal B = mod(p1.Y * p2.Y, N);
		mpreal C = mod(p1.Z * p2.T, N);
		mpreal D = mod(p1.T * p2.Z, N);
		mpreal E = mod(D + C, N);
		mpreal F = mod((p1.X - p1.Y) * (p2.X - p2.Y) + B - A, N);
		mpreal G = mod(B + (a * A), N);
		mpreal H = mod(D - C, N);
		p3.X = mod(E * F, N);
		p3.Y = mod(G * H, N);
		p3.Z = mod(F * G, N);
		p3.T = mod(E * H, N);
	}



	/*
	mpreal lambda, numerator, denominator;
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
	lambda = solveCongruence(numerator, denominator, N);

	//cout << "numerator: " << numerator << ", denominator: " << denominator << endl;
	//cout << "lambda: " << lambda << endl;

	sum.X = (lambda * lambda) - p1.X - p2.X;
	sum.X = mod(sum.X, N);
	sum.Y = lambda * (p1.X - sum.X) - p1.Y;
	sum.Y = mod(sum.Y, N);
	return sum;
	*/
}

/*
EllipticCurve::EllipticCurve(const Point& p, const mpreal& n) {
	// Generate elliptic curve containing the point p
	N = n;
	A = getLargeInt(2);	// How many digits for "random" A??
	B = (p.Y * p.Y) - (p.X * p.X * p.X) - (A * p.X);
	B = mod(B, N);
}
*/

void EllipticCurve::print() {
	cout << "X^2 + Y^2 = Z^2 + " << d << "T^2  (mod " << N << ")" << endl;
}









// Solve the congruence a = bx (mod n)
mpreal solveCongruence(const mpreal& a, const mpreal& b, const mpreal& n) {
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
		if (A == n) {
			cout << "Found n as a factor of n :(" << endl;
			getchar();
		}
		cout << "Just found a factor!! "<< A << " divides " << n << endl;
		getchar();
	}
	// x2 and y2 now contain values such that n*x2 + b*y2 = A

	//cout << x2 << " " << y2 << endl;

	mpreal inverse;
	if (abs(x2) > abs(y2)) {
		// x2 is the inverse of b (mod n)
		inverse = x2;
	}
	else if (abs(y2) > abs(x2)){
		// y2 is the inverse of b (mod n)
		inverse = y2;
	}

	// Working mod n
	if (inverse < 0) {
		inverse += n;
	}
	return mod(inverse * a, n);
}
