/*
 * HEcrypto.h
 *
 *  Created on: Jul 1, 2015
 *      Author: ph
 */

#ifndef LNV2011_HECRYPTO_H_
#define LNV2011_HECRYPTO_H_



#include <array>
#include <vector>

#include <common/JSONParser.h>
#include "BaseCrypto.h"
#include "HEctxt.h"
#include "HEpubkey.h"

using namespace std;
using namespace NTL;


class HEcrypto : public BaseCrypto{
private:
	// declare variables for common polynomials
	ZZ_pX _u, _g, _h;
	vector<ZZ> _modulis;

	JSONParser *_jsonParser;
	HEpubkey *_he_pubkey;
	LNV2011Params *_params;

	vector<ZZ_pX> convert_c2_intoBaseT_Polys(const ZZ_pX c2);
	ZZ_pX sum_c2i_ai(const vector<ZZ_pX> c2_baseT_polys, const vector<HEpubkey *> he_homorphismKeys);
	ZZ_pX sum_c2i_bi(const vector<ZZ_pX> c2_baseT_polys, const vector<HEpubkey *> he_homorphismKeys);
	void relinearization(ZZ_pX &c0_relin, ZZ_pX &c1_relin, const ZZ_pX c0, const ZZ_pX c1, const ZZ_pX c2, const vector<HEpubkey *> he_homorphismKeys);


public:
	HEcrypto(LNV2011Params *params, JSONParser *jsonParser);
	/*
	 * add public key to the crypto-system
	 */
	void addHEpubkey(HEpubkey *he_pubkey){
		delete _he_pubkey;
		_he_pubkey = he_pubkey;
	}

	virtual ~HEcrypto();

	void printParameters();

	void encrypt(ZZ_pX &c0, ZZ_pX &c1, const ZZ_pX m);
	void encrypt(ZZ_pX &c0, ZZ_pX &c1, const ZZ_pX m, const ZZ t, const HEpubkey * pk);
	HEctxt *encrypt(const ZZ_pX m);


	void decrypt(ZZ_pX &m, HEctxt *ctxt, const ZZ_pX s);
	void decryptWithoutModT(ZZX &m, HEctxt *ctxt, const ZZ_pX s);
	void decrypt(ZZ_pX &m, const ZZ_pX c0, const ZZ_pX c1, const ZZ_pX s);
	void decrypt(ZZ_pX &m, const ZZ_pX c0, const ZZ_pX c1, const ZZ_pX c2, const ZZ_pX s);
	void decrypt_withoutModT(ZZX &m, const ZZ_pX c0, const ZZ_pX c1, const ZZ_pX c2, const ZZ_pX s);
	HEctxt *addition(const HEctxt *m1_ctxt, const HEctxt *m2_ctxt);
	HEctxt *subtraction(const HEctxt *m1_ctxt, const HEctxt *m2_ctxt);
	HEctxt *addition_threeElements(const HEctxt *m1_ctxt, const HEctxt *m2_ctxt);
	HEctxt *subtraction_threeElements(const HEctxt *m1_ctxt, const HEctxt *m2_ctxt);
	void HE_multiplication(ZZ_pX &c0_relin, ZZ_pX &c1_relin, const ZZ_pX m, const ZZ_pX m1, const vector<HEpubkey *> he_homorphismKeys);
	void multiplication(ZZ_pX &new_c0, ZZ_pX &new_c1, ZZ_pX &new_c2, const ZZ_pX m, const ZZ_pX m1);
	void multiplication(ZZ_pX &new_c0, ZZ_pX &new_c1, ZZ_pX &new_c2, const HEctxt *m_ctxt, const HEctxt *m1_ctxt);
	void multiplication(ZZ_pX &new_c0, ZZ_pX &new_c1, const HEctxt *m_ctxt, const ZZ_pX m1);
	HEctxt *multiplication(const HEctxt *m_ctxt, const HEctxt *m1_ctxt);
	HEctxt *multiplication_with_relinearization(const HEctxt *m_ctxt, const HEctxt *m1_ctxt, const vector<HEpubkey *> he_homorphismKeys);



	vector<HEctxt *> encrypt_vector(const vector<ZZ_pX> m_vector, const HEpubkey *pk);
	vector<HEctxt *> encrypt_vector(const vector<ZZ_pX> m_vector, const vector<HEpubkey *> pk_vector);
	vector<HEctxt *> addition_vector(const vector<HEctxt *> m_ctxt_vector, const vector<HEctxt *> m1_ctxt_vector);
	vector<HEctxt *> subtraction_vector(const vector<HEctxt *> m_ctxt_vector, const vector<HEctxt *> m1_ctxt_vector);
	vector<HEctxt *> multiplication_vector(const vector<HEctxt *> m_ctxt_vector, const vector<HEctxt *> m1_ctxt_vector);
	vector<HEctxt *> multiplication_vector(const vector<HEctxt *> m_ctxt_vector, const vector<HEctxt *> m1_ctxt_vector, const vector<HEpubkey *> he_homorphismKeys);
	vector<ZZ_pX> decrypt_vector_different_t(const vector<HEctxt *> m_ctxt_vector, const ZZ_pX s);
	vector<ZZ_pX> decrypt_vector_same_t(const vector<HEctxt *> m_ctxt_vector, const ZZ_pX s);

	void delete_HEctxt_vector(vector<HEctxt *> vec);

	const ZZ_pX& get_g() const;
	const ZZ_pX& get_h() const;
	const ZZ_pX& get_u() const;

	const vector<ZZ>& getModulis() const {
		return _modulis;
	}
};


#endif /* LNV2011_HECRYPTO_H_ */
