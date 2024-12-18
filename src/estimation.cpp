/**
* @file estimation.cpp
 * @brief IDM - F2 ZZ3
 * @copyright Romain KLODZINSKI - ISIMA F2 ZZ3 - (c) 2024
 * */

#include <chrono>
#include <cmath>
#include <numbers>
#include <thread>
#include "settings.h"
#include "CLHEP/Random/MTwistEngine.h"

using std::chrono::high_resolution_clock;
using f64 = double;

/**
 * @brief represent the independent state of a simulation
 */
struct State {
	CLHEP::MTwistEngine rng; ///< the random number generator
	f64 value; ///< the estimated value of π
	std::chrono::duration<f64> duration; ///< the duration of the simulation
	std::thread thread; ///< the thread of the simulation

	/**
	 * @brief run the simulation to estimate the volume of a sphere of radius 1
	 * @param points the number of points to use (Monte Carlo method)
	 */
	void computeSphereVolume(const unsigned points = POINTS) {
		const auto start = high_resolution_clock::now();
		unsigned in = 0;

		for (unsigned i = points; i--;) {
			f64 sum = rng.flat(); // x
			sum *= sum;
			f64 c = rng.flat(); // y
			sum += c * c;
			c = rng.flat(); // z
			sum += c * c;

			if (sum < 1.) in++; // sqrt unnecessary here
		}

		value = 8. * in / points;
		duration = high_resolution_clock::now() - start;
		printf("estimation: %.08f (0x%016llx) in (%4.02f sec)\n",
			   value,
			   *reinterpret_cast<const unsigned long long*>(&value), // real f64 value as hex
			   duration.count());
	}

	/**
	 * @brief run the simulation in its own thread
	 * @param points the number of points to use (Monte Carlo method)
	 */
	void start(const unsigned points = POINTS) {
		thread = std::thread(&State::computeSphereVolume, this, points);
	}

	/**
	 * @brief initialize the random number generator with a status file
	 * @param seq the sequence number of the status file to load
	 */
	void initRng(const int seq) {
		char s[30];
		sprintf(s, "../status/status-%02d", seq);
		printf("loading '%s'...\n", s);

		rng.restoreStatus(s); // load the rng status from file
	}
};


/**
 * @brief print the result of a set of replicated experiments including mean, variance,
 * standard deviation, error against π and confidence interval at 99%
 * @param[in] mean the mean of the set
 * @param[in] variance the variance of the set
 * @param[in] replicates the number of replicates
 * @note reused and reworked from ZZ2 Simulation TP3 code
 * */
static void printResult(const f64 mean, const f64 variance, const unsigned replicates) {
	static const f64 STUDENT[] = {
		// Student law coeff (99% confidence), for k in : [0-30] U {40, 50, 60, 80, 100, 120, inf}
		INFINITY, 63.66, 9.925, 5.841, 4.604, 4.032, 3.707, 3.499, 3.355, 3.25, 3.169, 3.106, 3.055, 3.012, 2.977,
		2.947, 2.921, 2.898, 2.878, 2.861, 2.845, 2.831, 2.819, 2.807, 2.797, 2.787, 2.779, 2.771, 2.763, 2.756, 2.75,
		2.704, 2.678, 2.66, 2.639, 2.626, 2.617, 2.576
	};

	constexpr f64 real_value = 4. * std::numbers::pi / 3.;
	const f64 ub_variance = replicates * variance / (replicates - 1); // unbiased variance
	const f64 error = real_value - mean;
	f64 R = std::sqrt(ub_variance / replicates);

	if (replicates <= 30) R *= STUDENT[replicates]; // replicates=0 will cause R=inf -> interval = read number -> ok
	else if (replicates <= 60) R *= STUDENT[27 + replicates / 10]; // dispatch 10-sized intervals
	else if (replicates < 140) R *= STUDENT[30 + replicates / 20]; // dispatch 20-sized intervals
	else R *= STUDENT[37]; // infinite

	f64 location = (error + R) * 100. / R;
	if (location > 100.) location = 200. - location;

	printf("\nResults for %d replicates:\n", replicates);
	printf("\t- Mean :                         \t%.10f\n", mean);
	printf("\t- Variance :                     \t%.10f\n", variance);
	printf("\t- Unbiased variance :            \t%.10f\n", ub_variance);
	printf("\t- Standard deviation :           \t%.10f\n", std::sqrt(variance));
	printf("\t- Absolute error : 4π/3 - mean : \t%.10f\n", error);
	printf("\t- Relative error : Err / 4π/3 :  \t%.10f %%\n", 100. * error / real_value);
	printf("\t- Standard error :               \t%.10f\n", std::sqrt(variance / replicates)); // stddev / sqrt(n)
	printf("\t- Confidence interval :          \t[ %.10f ; %.10f ]\n", mean - R, mean + R);
	printf("\t- 4π/3 location in interval :    \t%.10f %%\n", location);
	// distance to the center: 100% = center, 0% = limit, -X% = out of interval
	printf("\t- Confidence radius :            \t%.10f\n\n", R);
}


/**
 * @brief main simulation function, start the replicates and show the results
 * */
int main() {
	State states[REPLICATES];
	f64 results[REPLICATES];

	// load the rng status from file
	for (int i = 0; i < REPLICATES; i++) states[i].initRng(i);

	// start the simulation (separated from the previous loop to avoid loading time in the timing)
	for (auto& state : states) state.start();

	printf("\nrunning (thread)...\n");

	f64 mean = 0., squared_mean = 0.;

	for (int i = 0; i < REPLICATES; i++) {
		auto& state = states[i];
		state.thread.join();
		results[i] = state.value;
		mean += state.value;
		squared_mean += state.value * state.value;
	}

	mean /= REPLICATES;
	squared_mean /= REPLICATES;

	printResult(mean, squared_mean - mean * mean, REPLICATES);

	printf("\nrunning (sequential)...\n");

	// load the rng status from file
	for (int i = 0; i < REPLICATES; i++) states[i].initRng(i);

	std::chrono::duration<f64> duration{};

	for (int i = 0; i < REPLICATES; i++) {
		states[i].computeSphereVolume();

		// check if value are exactly the same (need a cast to u64 to avoid test bad optimization)
		if (*reinterpret_cast<unsigned long long*>(&states[i].value) == *reinterpret_cast<unsigned long long*>(&results[
			i])) {
			printf("reproductibility confirmed\n");
		}
		else {
			printf("reproductibility issue %.08f (0x%016llx) vs %.08f (0x%016llx)\n",
				   states[i].value,
				   *reinterpret_cast<const unsigned long long*>(&states[i].value),
				   results[i],
				   *reinterpret_cast<const unsigned long long*>(&results[i])
			);
		}

		duration += states[i].duration;
	}

	printf("Sequential time: %4.02f sec\n", duration.count());

	return 0;
}
