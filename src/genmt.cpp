/**
* @file genmt.cpp
 * @brief IDM - F2 ZZ3
 * @copyright Romain KLODZINSKI - ISIMA F2 ZZ3 - (c) 2024
 * */

#include <sys/types.h>
#include "settings.h"
#include "CLHEP/Random/MTwistEngine.h"


/**
 * @brief generate 'REPLICATES' independent status of the random number generator CLHEP::MTwistEngine in the 'status' directory
* */
int main() {
	// unhashed seed with equal number of 0 and 1
	CLHEP::MTwistEngine mtRng{0b10101010101010101010101010101010};

	for (int i = 0; i < REPLICATES; i++) {
		char s[30];
		sprintf(s, "../status/status-%02d", i);

		printf("Computing status: '%s'...\n", s);

		for (unsigned u = SIZE; u--;) (void)mtRng.flat();

		printf("saved %s\n", s);
		mtRng.saveStatus(s);
	}
}
