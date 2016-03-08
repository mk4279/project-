/*
 * Fujitsu.cpp
 *
 *  Created on: Aug 22, 2015
 *      Author: ph
 */

#include <common/libntl_addon.h>
#include <common/util.h>


#include "HEcrypto.h"
#include "HEkeygen.h"
#include "HEpubkey.h"
#include "HEctxt.h"


using namespace std;
using namespace NTL;

ZZX reverse_pack(const ZZX poly, long degree){
	ZZX poly_r=ZZX(poly);
	flipPloy(poly_r, degree+1);
	poly_r *= -1;
	return poly_r;
}

void runFujitsuHammingDistanceExample(LNV2011Params *params, HEcrypto* he_crypto, HEkeygen *he_keygen, HEpubkey *he_pubkey, const ZZ_pX sk){
	ZZX phi_ZZX = to_ZZX(he_crypto->get_phi());
	ZZX identity, identity_r, pattern_a, pattern_b, pattern_b_r;

	identity = createZZXPolyOfOne(params->n);
	printZZX("identity", identity);

	identity_r = reverse_pack(identity, params->n);
	printZZX("identity_r", identity_r);

	while(pattern_a==0)
		pattern_a = generatePolyZZXWithCoeffFromRq(params->n);
	modT(pattern_a, params->t);
	printZZX("pattern_a", pattern_a);

	while(pattern_b==0)
		pattern_b = generatePolyZZXWithCoeffFromRq(params->n);
	modT(pattern_b, params->t);
	printZZX("pattern_b", pattern_b);

	pattern_b_r = reverse_pack(pattern_b, params->n);
	printZZX("pattern_b_r", pattern_b_r);

	ZZX HD_a, HD_b, HD_ab, HD_result;

	HD_a = pattern_a * identity_r;
	modF(HD_a, phi_ZZX);
#ifdef DEBUG
	printZZX("HD_a", HD_a);
#endif
	HD_b = pattern_b_r * identity;
	modF(HD_b, phi_ZZX);
#ifdef DEBUG
	printZZX("HD_b", HD_b);
#endif

	HD_ab = pattern_a * pattern_b_r;
	modF(HD_ab, phi_ZZX);
#ifdef DEBUG
	printZZX("HD_ab", HD_ab);
#endif

	HD_result = HD_a + HD_b - 2 * HD_ab;
	printZZX("result_plain", HD_result);


	/*
	 * encrypted version
	 */
	ZZ_pX identity_ZZ_pX = to_ZZ_pX(identity);
	ZZ_pX identity_r_ZZ_pX = to_ZZ_pX(identity_r);
	ZZ_pX pattern_a_ZZ_pX = to_ZZ_pX(pattern_a);
	ZZ_pX pattern_b_r_ZZ_pX = to_ZZ_pX(pattern_b_r);

	ZZ_pX new_a0, new_a1;
	HEctxt *enc_pattern_a = he_crypto->encrypt(pattern_a_ZZ_pX);
	he_crypto->multiplication(new_a0, new_a1, enc_pattern_a, identity_r_ZZ_pX);

	ZZ_pX new_b0, new_b1;
	HEctxt *enc_pattern_b_r = he_crypto->encrypt(pattern_b_r_ZZ_pX);
	he_crypto->multiplication(new_b0, new_b1, enc_pattern_b_r, identity_ZZ_pX);

	ZZ_pX new_ab0, new_ab1, new_ab2;
	he_crypto->multiplication(new_ab0, new_ab1, new_ab2, pattern_a_ZZ_pX, pattern_b_r_ZZ_pX);

	ZZ_pX new_result0, new_result1, new_result2;
	new_result0 = new_a0 + new_b0 - 2 * new_ab0;
	new_result1 = new_a1 + new_b1 - 2 * new_ab1;
	new_result2 = -2*new_ab2;

	ZZX result_dec;
	he_crypto->decrypt_withoutModT(result_dec, new_result0, new_result1, new_result2, sk);
	modT(result_dec, to_ZZ(params->n));
	printZZX("result_dec", result_dec);


	//vector<HEpubkey *> he_homomorphismKeys = he_keygen->createHomomorphismKeys(sk);

	delete enc_pattern_a;
	delete enc_pattern_b_r;
}


void runFujitsuEuclideanDistanceExample(LNV2011Params *params, HEcrypto* he_crypto, HEkeygen *he_keygen, HEpubkey *he_pubkey, const ZZ_pX sk){
	ZZX phi_ZZX = to_ZZX(he_crypto->get_phi());
	ZZX pattern_a, pattern_a_r, pattern_b, pattern_b_r;

	while(pattern_a==0)
		pattern_a = generatePolyZZXWithCoeffFromRq(params->n);
	//modT(pattern_a, params->t);
	printZZX("pattern_a", pattern_a);

	while(pattern_b==0)
		pattern_b = generatePolyZZXWithCoeffFromRq(params->n);
	//modT(pattern_b, params->t);
	printZZX("pattern_b", pattern_b);


	pattern_a_r = reverse_pack(pattern_a, params->n);
	printZZX("pattern_a_r", pattern_a_r);
	pattern_b_r = reverse_pack(pattern_b, params->n);
	printZZX("pattern_b_r", pattern_b_r);

	ZZX ED_a, ED_b, ED_ab, ED_result;

	ED_a = pattern_a * pattern_a_r;
	modF(ED_a, phi_ZZX);
#ifdef DEBUG
	printZZX("ED_a", ED_a);
#endif
	ED_b = pattern_b_r * pattern_b;
	modF(ED_b, phi_ZZX);
#ifdef DEBUG
	printZZX("ED_b", ED_b);
#endif

	ED_ab = pattern_a * pattern_b_r;
	modF(ED_ab, phi_ZZX);
#ifdef DEBUG
	printZZX("ED_ab", ED_ab);
#endif

	ED_result = ED_a + ED_b - 2 * ED_ab;
	modQ(ED_result, he_crypto->get_q());
	printZZX("result_plain (A^2+B^2-2AB)", ED_result);
	ZZX diff_ab = pattern_a - pattern_b;
	ZZX diff_ab_r = reverse_pack(diff_ab, params->n);
	ED_result = diff_ab * diff_ab_r;
	modF(ED_result, phi_ZZX);
	modQ(ED_result, he_crypto->get_q());
	printZZX("result_plain (A-B)*(A-B)", ED_result);


	ZZ_pX pattern_a_ZZ_pX = to_ZZ_pX(pattern_a);
	ZZ_pX pattern_a_r_ZZ_pX = to_ZZ_pX(pattern_a_r);
	ZZ_pX pattern_b_ZZ_pX = to_ZZ_pX(pattern_b);
	ZZ_pX pattern_b_r_ZZ_pX = to_ZZ_pX(pattern_b_r);

	ZZ_pX new_aa0, new_aa1, new_aa2;
	he_crypto->multiplication(new_aa0, new_aa1, new_aa2, pattern_a_ZZ_pX, pattern_a_r_ZZ_pX);

	ZZ_pX new_bb0, new_bb1, new_bb2;
	he_crypto->multiplication(new_bb0, new_bb1, new_bb2, pattern_b_ZZ_pX, pattern_b_r_ZZ_pX);

	ZZ_pX new_ab0, new_ab1, new_ab2;
	he_crypto->multiplication(new_ab0, new_ab1, new_ab2, pattern_a_ZZ_pX, pattern_b_r_ZZ_pX);

	ZZ_pX new_result0, new_result1, new_result2;
	new_result0 = new_aa0 + new_bb0 - 2 * new_ab0;
	new_result1 = new_aa1 + new_bb1 - 2 * new_ab1;
	new_result2 = new_aa2 + new_bb2 - 2 * new_ab2;

	ZZX result_dec;
	he_crypto->decrypt_withoutModT(result_dec, new_result0, new_result1, new_result2, sk);
	printZZ_pX("result_dec", to_ZZ_pX(result_dec));
}





