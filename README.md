# FLASHlikeSteadyState 
Computes the steady state of "FLASH-like" sequence (e.g., MPRAGE) in the ideal spoiling condition.

This code allows the efficient solution of the steady state transverse (and longitudinal) signal from sequences that combine FLASH readouts with preparation pulses and delays (e.g., MPRAGE or MP2RAGE). This is done by breaking the sequence into "events": FLASH TR, preparation pulse, or relaxation delay. The simulation assumes perfect spoiling of transverse magnetization after each event. Under this assumption, the steady state condition can be expressed as the solution of a linear equation with a nearly banded structure (one outlier entry). We use the LU decomposition of this matrix structure to find the solution.

The main simulation code is in `SequenceEventVectors.h`, while the LU solver for our specific matrix is in `StructuredLUSolver.h`. For an example of how to specify and simulate an MPRAGE sequence, see the `examples` directory.
