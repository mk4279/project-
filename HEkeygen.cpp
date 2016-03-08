/*
 * HEkeygen.cpp
 *
 *  Created on: Jul 1, 2015
 *      Author: ph
 */

#include <common/distribution.h>
#include <common/libntl_addon.h>
#include <common/util.h>

#include "LNV2011.h"
#include "HEkeygen.h"


HEkeygen::HEkeygen(LNV2011Params *params, HEcrypto *he_crypto, JSONParser *jsonParser) {
	_he_crypto = he_crypto;
	_jsonParser = jsonParser;
	_params = params;
}



HEkeygen::~HEkeygen() {

}




HEpubkey* HEkeygen::createHEpubkey(const ZZ_pX s){
	if(_jsonParser->isLoadedFromFile()){
		return new HEpubkey(_params, _jsonParser->readJSONZZ_pX("a1"), _he_crypto->get_phi(), _he_crypto->get_q(), _he_crypto->get_e(), s);
	}else{
		return new HEpubkey(_params, _he_crypto->get_phi(), _he_crypto->get_q(), _he_crypto->get_e(), s);
	}
}


vector<HEpubkey *> HEkeygen::createHEpubkey(const ZZ_pX s, const vector<ZZ> t_vector){
	vector<HEpubkey *> pk_vector;
	for(int i=0;i<t_vector.size();i++){
		_params->t=t_vector[i];
		pk_vector.push_back(new HEpubkey(_params, _he_crypto->get_phi(), _he_crypto->get_q(), _he_crypto->get_e(), s));
	}
	return pk_vector;
}



vector<HEpubkey *> HEkeygen::createHomomorphismKeys(const ZZ_pX s){
	int v_size= log_baseT(_params->t, _he_crypto->get_q());
	vector<HEpubkey *> hKeys(v_size+1);

	for(int i=0; i<=v_size; ++i){
		hKeys[i] = new HEpubkey(_params, _he_crypto->get_phi(), _he_crypto->get_q(), s, i);
	}

	// the following code is for debugging
#ifdef DEBUG
	printHeader("Generate homomorphic key");
	int i=0;
	for(HEpubkey *k : hKeys){
		printINT("i", i);
		printZZ_pX("a", k->get_a0());
		printZZ_pX("b", k->get_a1());
		printZZ_pX("e", k->get_HomomorphismKey_e());
		++i;
	}
#endif

	return hKeys;
}



void HEkeygen::createHomomorphismKeyJSONArrayElt(json_t *hkey_json_obj, const int index, const HEpubkey *hKey){
	json_t *hkey_elems = json_object();

	json_t *tmp_a = _jsonParser->encodeZZ_pXToJSON(hKey->get_a0());
	json_object_set(hkey_elems, "a", tmp_a);
	json_decref(tmp_a);

	json_t *tmp_b = _jsonParser->encodeZZ_pXToJSON(hKey->get_a1());
	json_object_set(hkey_elems, "b", tmp_b);
	json_decref(tmp_b);

	json_t *tmp_e = _jsonParser->encodeZZ_pXToJSON(hKey->get_HomomorphismKey_e());
	json_object_set(hkey_elems, "e", tmp_e);
	json_decref(tmp_e);

	json_object_set(hkey_json_obj, to_string(index).c_str(), hkey_elems);
	json_decref(hkey_elems);
}


json_t *HEkeygen::convertHomomorphismKeysToJSONObj(vector<HEpubkey *> hKeys){
	json_t *hkeys_json = json_object();
	for(int i=0; i<hKeys.size(); i++){
		//json_t *hkey_json_obj = json_object();
		createHomomorphismKeyJSONArrayElt(hkeys_json, i, hKeys[i]);
		//json_array_append(hkeys_json_arr, hkey_json_obj);
	}
	json_t *tmp=json_integer(hKeys.size());
	json_object_set(hkeys_json, "size", tmp);
	json_decref(tmp);
	return hkeys_json;
}


// internal method, don't need to json_decref(), we do it in the upper level method
HEpubkey *HEkeygen::extractHEpubkey(json_t *hkey_obj){
	json_t *hkey_a = json_object_get(hkey_obj, "a");
	ZZ_pX a = _jsonParser->extractZZ_pX(hkey_a);

	json_t *hkey_b = json_object_get(hkey_obj, "b");
	ZZ_pX b = _jsonParser->extractZZ_pX(hkey_b);

	json_t *hkey_e = json_object_get(hkey_obj, "e");
	ZZ_pX e = _jsonParser->extractZZ_pX(hkey_e);

#ifdef DEBUG
	printZZ_pX("a", a);
	printZZ_pX("b", b);
	printZZ_pX("e", e);
#endif

	return new HEpubkey(_params, a, b, e);
}


vector<HEpubkey *> HEkeygen::reconstructHomomorphismKeysFromJSONObj(json_t *jsonData){
	vector<HEpubkey *> hKeys;
	const int length = json_integer_value(json_object_get(jsonData, "size"));
#ifdef DEBUG
	printHeader("Homomorphism Keys");
#endif

	for(int i=0; i<length; i++){
		json_t *hkey_obj = json_object_get(jsonData, to_string(i).c_str());
#ifdef DEBUG
		printINT("i", i);
#endif
		hKeys.push_back(extractHEpubkey(hkey_obj));
	}

	return hKeys;
}

void HEkeygen::destoryhKeys(vector<HEpubkey *> hKeys){
	for(int i=0; i<hKeys.size(); i++){
		HEpubkey *tmp = hKeys[i];
		tmp->~HEpubkey();
		delete tmp;
	}
}



