#include <stdio.h>

int main(int argc, char** argv) {
	double e = 1;
	double fact = 1;
	long n = 100;
	for (int i = 1; i <= n; i++) {
		e += 1 / fact;
		fact *= (i + 1);
	}
	printf("e = %.20lf\n", e);
	return 0;
}
