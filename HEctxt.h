/*
 * HEciphertext.h
 *
 *  Created on: Jul 21, 2015
 *      Author: ph
 */

#ifndef LNV2011_HECTXT_H_
#define LNV2011_HECTXT_H_

#include <NTL/ZZ_pX.h>

using namespace std;
using namespace NTL;


class HEctxt {
public:
	HEctxt();
	HEctxt(const ZZ_pX c0, const ZZ_pX c1);
	HEctxt(const ZZ_pX c0, const ZZ_pX c1, const ZZ_pX c2);
	virtual ~HEctxt();

	const ZZ_pX& getC0() const {
		return _c0;
	}

	const ZZ_pX& getC1() const {
		return _c1;
	}

	const ZZ_pX& getC2() const {
		return _c2;
	}


private:
	ZZ_pX _c0, _c1, _c2;
};




#endif /* LNV2011_HECTXT_H_ */
