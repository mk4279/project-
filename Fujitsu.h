/*
 * Fujitsu.h
 *
 *  Created on: Aug 22, 2015
 *      Author: ph
 */

#ifndef LNV2011_FUJITSU_H_
#define LNV2011_FUJITSU_H_


void runFujitsuHammingDistanceExample(LNV2011Params *params, HEcrypto* he_crypto, HEkeygen *he_keygen, HEpubkey *he_pubkey, const ZZ_pX sk);

void runFujitsuEuclideanDistanceExample(LNV2011Params *params, HEcrypto* he_crypto, HEkeygen *he_keygen, HEpubkey *he_pubkey, const ZZ_pX sk);


#endif /* LNV2011_FUJITSU_H_ */
