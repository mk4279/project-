#ifndef LNV2011_H_INCLUDED__
#define LNV2011_H_INCLUDED__


// uncomment to test Rq+r, so that the coefficient is not from 0 to q, but TestR to q+TestR;
// *** this option can be use in addition to the MixQandQhalf ONLY!!!
//#define ShiftRqByTestR
#define TestR 100   // for TestR within (0, q]


struct LNV2011Params{
	int mean=0;     			 	// mean for Gaussian distribution; default to 0
	double standard_deviation=1;	// sd for Gaussian distribution; default to 1
	int q_length_in_bits=89;		// determine the length of prime q; 2^q_length_in_bits; at least 15
	NTL::ZZ t=to_ZZ(2);
	unsigned n=16;
};




#endif
