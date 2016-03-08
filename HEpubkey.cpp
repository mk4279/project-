/*
 * HEpubkey.cpp
 *
 *  Created on: Jul 1, 2015
 *      Author: ph
 */

#include <common/distribution.h>
#include <common/libntl_addon.h>
#include <common/util.h>
#include "LNV2011.h"
#include "HEpubkey.h"


HEpubkey::HEpubkey(LNV2011Params *params){
	_params=params;
}

HEpubkey::HEpubkey(LNV2011Params *params, const ZZ_pX f, const ZZ q, const ZZ_pX e, const ZZ_pX s){
	_params=params;
	instantiateHEpubkey(f, q, e, s);
}

HEpubkey::HEpubkey(LNV2011Params *params, const ZZ_pX f, const ZZ q, const ZZ_pX s, const int i){
	_params=params;
	instantiateHomomorphismKey(f, q, s, i);
}

HEpubkey::HEpubkey(LNV2011Params *params,const ZZ_pX a1, const ZZ_pX f, const ZZ q, const ZZ_pX e, const ZZ_pX s){
	_params=params;
	instantiateHEpubkey(a1, f, q, e, s);
}

HEpubkey::HEpubkey(LNV2011Params *params, const ZZ_pX a, const ZZ_pX b, const ZZ_pX e){
	_params=params;
	_a0=a;
	_a1=b;
	_e=e;
}

HEpubkey::~HEpubkey() {
	// TODO Auto-generated destructor stub
}


/*
 * a0 = −(a1s + te), a1
 */
void HEpubkey::instantiateHEpubkey(const ZZ_pX f, const ZZ q, const ZZ_pX e, const ZZ_pX s){
	while(_a1==0)
		_a1 = generatePolyWithCoeffFromRq(_params->n);
#ifdef ShiftRqByTestR
	modQ(_a1, q);
#endif

	ZZ_pX a1_s = MulMod_f(_a1, s, f);
	ZZ_pX t_e = MulMod_f(_params->t, e, f);
	//printZZ("t", _params->t);
	_a0 = -(a1_s + t_e);
	_a0 = Mod_f(_a0, f);
}

void HEpubkey::instantiateHEpubkey(const ZZ_pX a1, const ZZ_pX f, const ZZ q, const ZZ_pX e, const ZZ_pX s){
	_a1 = a1;
	//modQ(_a1, q);

	ZZ_pX a1_s = MulMod_f(_a1, s, f);
	ZZ_pX t_e = MulMod_f(_params->t, e, f);
	_a0 = -(a1_s + t_e);
	_a0 = Mod_f(_a0, f);
}

/*
 * ai, bi=-(a*s + t*ei) + t^i * s^2
 */
void HEpubkey::instantiateHomomorphismKey(const ZZ_pX f, const ZZ q, const ZZ_pX s, const int i){
	while(_a1==0)
		_a1 = generatePolyWithCoeffFromRq(_params->n);
#ifdef ShiftRqByTestR
	modQ(_a1, q);
#endif

	while(_e==0)
		_e = generatePolyWithCoeffFromGaussianNoise(_params->n, _params->mean, _params->standard_deviation);

	ZZ_pX a_s = MulMod_f(_a1, s, f);
	ZZ_pX t_e = MulMod_f(_params->t, _e, f);
	ZZ t_i = power(_params->t, i);
	ZZ_pX s_2 = MulMod_f(s, s, f);

	ZZ_pX ti_s2 = MulMod_f(t_i, s_2, f);
	_a0 = -(a_s + t_e) + ti_s2;
	_a0 = Mod_f(_a0, f);
}



const ZZ_pX& HEpubkey::get_a0() const {
	return _a0;
}

const ZZ_pX& HEpubkey::get_a1() const {
	return _a1;
}

const ZZ_pX& HEpubkey::get_HomomorphismKey_e() const {
	return _e;
}
