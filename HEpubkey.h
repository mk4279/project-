/*
 * HEpubkey.h
 *
 *  Created on: Jul 1, 2015
 *      Author: ph
 */

#ifndef LNV2011_HEPUBKEY_H_
#define LNV2011_HEPUBKEY_H_


using namespace std;
using namespace NTL;

class HEpubkey {
protected:
	// define the public key pair
	ZZ_pX _a0, _a1, _e;
	LNV2011Params *_params;

	void instantiateHEpubkey(const ZZ_pX f, const ZZ q, const ZZ_pX e, const ZZ_pX s);
	void instantiateHomomorphismKey(const ZZ_pX f, const ZZ q, const ZZ_pX s, const int i);
	void instantiateHEpubkey(const ZZ_pX a1, const ZZ_pX f, const ZZ q, const ZZ_pX e, const ZZ_pX s);

public:
	HEpubkey(LNV2011Params *params);
	HEpubkey(LNV2011Params *params, const ZZ_pX a, const ZZ_pX b, const ZZ_pX e);
	HEpubkey(LNV2011Params *params, const ZZ_pX f, const ZZ q, const ZZ_pX e, const ZZ_pX s);
	HEpubkey(LNV2011Params *params, const ZZ_pX f, const ZZ q, const ZZ_pX s, const int i);
	HEpubkey(LNV2011Params *params, const ZZ_pX a1, const ZZ_pX f, const ZZ q, const ZZ_pX e, const ZZ_pX s);
	virtual ~HEpubkey();

	const ZZ_pX& get_a0() const;
	const ZZ_pX& get_a1() const;

	const ZZ_pX& get_HomomorphismKey_e() const;
};

#endif /* LNV2011_HEPUBKEY_H_ */
