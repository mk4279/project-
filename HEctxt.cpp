/*
 * HEciphertext.cpp
 *
 *  Created on: Jul 21, 2015
 *      Author: ph
 */

#include "HEctxt.h"

HEctxt::HEctxt() {
}

HEctxt::HEctxt(const ZZ_pX c0, const ZZ_pX c1) {
	_c0=c0;
	_c1=c1;
}

HEctxt::HEctxt(const ZZ_pX c0, const ZZ_pX c1, const ZZ_pX c2) {
	_c0=c0;
	_c1=c1;
	_c2=c2;
}


HEctxt::~HEctxt() {
	// TODO Auto-generated destructor stub
}

