/*
 * BaseCrypto.cpp
 *
 *  Created on: Jul 7, 2015
 *      Author: ph
 */

#include <common/distribution.h>
#include <common/JSONParser.h>
#include <common/libntl_addon.h>
#include <common/util.h>

#include "LNV2011.h"

class BaseCrypto {
public:
	/*
	 * print basic parameters
	 */
	void printParameters(JSONParser *jsonParser){
		printZZ("q", _q);
		jsonParser->addZZToJSON("q", _q);
		//cout<<"after print q"<< endl;

		printZZ_pX("phi", _phi);
		jsonParser->addZZ_pXToJSON("phi", _phi);
		printZZ_pX("e", _e);
		jsonParser->addZZ_pXToJSON("e", _e);
	}


	void generateParameters(LNV2011Params *params, JSONParser *jsonParser){
		if(jsonParser->isLoadedFromFile()){
			_q = jsonParser->readJSONZZ("q");
			ZZ_p::init(_q);
			_phi = jsonParser->readJSONZZ_pX("phi");
			_e = jsonParser->readJSONZZ_pX("e");
		}else{
			// TODO temporary set prime to qMargin bits large than n; define in config.h
			_q = RandomPrime_ZZ(params->q_length_in_bits);
			assert(_q>0);
			//cout<< "_q " << _q << endl;

			ZZ_p::init(_q);  // set modulo p for ZZ_p

			// initialize _f to x^n+1
			SetCoeff(_phi, 0, 1);
			SetCoeff(_phi, params->n, 1);

			while(_e==0)
				_e = generatePolyWithCoeffFromGaussianNoise(params->n, params->mean, params->standard_deviation);
		}
	}


	/*
	 * the prime number q
	 */
	const ZZ& get_q() const {
		return _q;
	}


	/*
	 * get the noise term e
	 */
	const ZZ_pX& get_e() const {
		return _e;
	}

	/*
	 * get the polynomial function f
	 */
	const ZZ_pX& get_phi() const {
		return _phi;
	}


protected:
	//example prime numbers  http://gabaix.us/blog/wp-content/uploads/2012/11/Screen-Shot-2012-11-23-at-4.36.53-PM.png
	ZZ _q;

	ZZ_pX _phi, _e;

};



