/*
 * HEkeygen.h
 *
 *  Created on: Jul 1, 2015
 *      Author: ph
 */

#ifndef LNV2011_HEKEYGEN_H_
#define LNV2011_HEKEYGEN_H_

#include <vector>

#include <common/JSONParser.h>
#include "HEcrypto.h"
#include "HEpubkey.h"

using namespace std;
using namespace NTL;

class HEkeygen {
private:
	HEcrypto *_he_crypto;
	LNV2011Params *_params;
	JSONParser *_jsonParser;

	HEpubkey *extractHEpubkey(json_t *hkey_obj);
	void createHomomorphismKeyJSONArrayElt(json_t *hkey_json_obj, const int index, const HEpubkey *hKey);

public:
	HEkeygen(LNV2011Params *params, HEcrypto *he_crypto,  JSONParser *jsonParser);
	virtual ~HEkeygen();

	HEpubkey *createHEpubkey(const ZZ_pX s);
	vector<HEpubkey *> createHEpubkey(const ZZ_pX s, const vector<ZZ> t_vector);
	vector<HEpubkey *> createHomomorphismKeys(const ZZ_pX s);

	json_t *convertHomomorphismKeysToJSONObj(vector<HEpubkey *> hKeys);
	vector<HEpubkey *> reconstructHomomorphismKeysFromJSONObj(json_t *jsonData);

	void destoryhKeys(vector<HEpubkey *> hKeys);
};

#endif /* LNV2011_HEKEYGEN_H_ */
