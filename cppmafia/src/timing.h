#ifndef MAFIA_TIMING_H_
#define MAFIA_TIMING_H_

/** @file timing.h a simple timing implementation */
enum TimePhase {
	PhaseNone = -1,
	PhaseReadData = 0, PhaseBuildHisto, PhaseBuildWindows, PhaseBuildBitmaps,
	PhaseFindCdus, // includes deduplication
	PhaseFindDense,
	PhaseFindUnjoined, 
	PhaseBuildGraph,
	PhaseBuildClusters,
	PhaseWriteData,
	PhaseEnd  // ending phase of the program
};

/** ends the old phase and starts the new one */
void start_phase(int phase);

/** gets the time elapsed so far for each phase */
double phase_time(int phase);

/** gets the name of the phase */
const char* phase_name(int phase);

/** prints all timing data */
void print_timing_info();

#endif
