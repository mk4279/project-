/*
 * HEcrypto.cpp
 *
 *  Created on: Jul 1, 2015
 *      Author: ph
 */
#include <vector>

#include <common/distribution.h>
#include <common/libntl_addon.h>
#include <common/util.h>

#include "LNV2011ZZX.h"
#include "HEcryptoZZX.h"




HEcryptoZZX::HEcryptoZZX(LNV2011ZZXParams *params, JSONParser *jsonParser){
	// just to initialize it with empty pubkey; we set a pubkey later
	_he_pubkey=new HEpubkeyZZX(params);
	_jsonParser = jsonParser;
	_params = params;

	// generate basic parameters f, e, q
	generateParameters(params, _jsonParser);

	if(_jsonParser->isLoadedFromFile()){
		_u = _jsonParser->readJSONZZX("u");
		_g = _jsonParser->readJSONZZX("g");
		_h = _jsonParser->readJSONZZX("h");
	}else{
		_u = generatePolyZZXWithCoeffFromGaussianNoise(_params->n, _params->mean, _params->standard_deviation);
		_g = generatePolyZZXWithCoeffFromGaussianNoise(_params->n, _params->mean, _params->standard_deviation);
		_h = generatePolyZZXWithCoeffFromGaussianNoise(_params->n, _params->mean, _params->standard_deviation);
	}
}


HEcryptoZZX::~HEcryptoZZX() {

}

void HEcryptoZZX::printParameters(){
	BaseCryptoZZX::printParameters(_jsonParser);
	printZZX("u", _u);
	_jsonParser->addZZXToJSON("u", _u);
	printZZX("g", _g);
	_jsonParser->addZZXToJSON("g", _g);
	printZZX("h", _h);
	_jsonParser->addZZXToJSON("h", _h);
}


void HEcryptoZZX::encrypt(ZZX &c0, ZZX &c1, const ZZX m){
	if(deg(_he_pubkey->get_a0())==-1) {
		cout << "error! add public key to the crypto object" << endl;
		exit(0);
	}

	// compute c0
	ZZX a0_u = MulMod_f(_he_pubkey->get_a0(), _u, _f);
	ZZX t_g = MulMod_f(_params->t, _g, _f);
	c0 = (a0_u +t_g + m);
	c0 = Mod_f(c0, _f);

	// compute c1
	ZZX a1_u = MulMod_f(_he_pubkey->get_a1(), _u, _f);
	ZZX t_h = MulMod_f(_params->t, _h, _f);
	c1 = a1_u + t_h;
	c1 = Mod_f(c1, _f);
}


HEctxtZZX *HEcryptoZZX::encrypt(const ZZX m){
	ZZX c0, c1;
	encrypt(c0, c1, m);
	return new HEctxtZZX(c0, c1);
}


void HEcryptoZZX::decrypt(ZZX &m, const ZZX c0, const ZZX c1, const ZZX s){
	ZZX c1_s = MulMod_f(c1, s, _f);
	m = c0 + c1_s;
	m = Mod_f(m, _f);
	modQhalf(m, _q);
	modT(m, _params->t);

}

void HEcryptoZZX::decrypt(ZZX &m, HEctxtZZX *ctxt, const ZZX s){
	decrypt(m, ctxt->getC0(), ctxt->getC1(), s);
}

void HEcryptoZZX::decrypt(ZZX &m, const ZZX c0, const ZZX c1, const ZZX c2, const ZZX s){
	ZZX c1_s = c1 * s;   		  // c1 * s
#ifdef VERBOSE
	printZZX("c1 * s", c1_s);
#endif
	c1_s= Mod_f(c1_s, _f);   // m = c0 + c1*s
#ifdef VERBOSE
	printZZX("c1 * s", c1_s);
#endif

	ZZX ss = s * s;				  // s*s
#ifdef VERBOSE
	printZZX("s * s", ss);
#endif
	ss = Mod_f(ss, _f);
#ifdef VERBOSE
	printZZX("s * s", ss);
#endif

	ZZX c2_ss = c2 * ss;		  // c2 * ss
#ifdef VERBOSE
	printZZX("c2 * ss", c2_ss);
#endif
	c2_ss= Mod_f(c2_ss, _f);  // m = c0 + c1*s + c2*s*s
#ifdef VERBOSE
	printZZX("c2 * ss", c2_ss);
#endif

	m = c0 + c1_s + c2_ss;

#ifdef VERBOSE
	printZZX("mm", m);
#endif
	m = Mod_f(m, _f);
#ifdef VERBOSE
	printZZX("mm", m);
#endif

	modQhalf(m, _q);
	modT(m, _params->t);
}


/*
 * HE addition
 */

HEctxtZZX *HEcryptoZZX::addition(HEctxtZZX *m1_ctxt, HEctxtZZX *m2_ctxt){
	ZZX c0_sum, c1_sum;
	c0_sum = m1_ctxt->getC0() + m2_ctxt->getC0();
	c0_sum = Mod_f(c0_sum, _f);

	c1_sum = m1_ctxt->getC1() + m2_ctxt->getC1();
	c1_sum = Mod_f(c1_sum, _f);

#ifdef DEBUG
	printZZX("c0 + c0_1", c0_sum);
	printZZX("c1 + c1_1", c1_sum);
#endif

	return new HEctxtZZX(c0_sum, c1_sum);
}


vector<ZZX> HEcryptoZZX::convert_c2_intoBaseT_Polys(const ZZX c2){
	ZZ tmp_coeff;

	// convert all coefficients in c2 into baseT form
	// check correctness here: http://korn19.ch/coding/base_converter.php
	vector< vector<ZZ> > coeffs_baseT_vector;
	for(int j=0; j<=deg(c2); j++){
		tmp_coeff = conv<ZZ>(coeff(c2, j));
#ifdef VERBOSE
		printZZ("old", tmp_coeff);
#endif
		// we add q if the coefficient is less than 0;
		// TODO try with preserving - sign
		if(tmp_coeff<0) tmp_coeff += _q;

		//tmp_coeff = abs(tmp_coeff);
#ifdef VERBOSE
		printZZ("new", tmp_coeff);
#endif
		coeffs_baseT_vector.push_back (convertCoeffToBaseT(tmp_coeff, _q, _params->t));
	}

	/*
	 * now create new ceil(log_t(q)) number of polynomials, c2,i
	 */
	// e.g., t=2, q=97, v_size is 7; we are doing int log, so it round down +1 is correct
	// TODO find proper log_t function
	int v_size= log_baseT(_params->t, _q)+1;
	vector<ZZX> c2_basedT_polys(v_size);
	for(int i=0; i<v_size;i++){
		ZZX new_poly;
		for(int j=0; j<=deg(c2); j++){
			vector<ZZ> coeff_baseT =  coeffs_baseT_vector[j];
			//			if(j==0 || j==deg(c2))    // for setting the first and the last term negative
			//				SetCoeff(new_poly, j, -coeff_baseT[c2_basedT_polys.size()-1-i]);
			//			else
			SetCoeff(new_poly, j, coeff_baseT[c2_basedT_polys.size()-1-i]);
		}

#ifdef DEBUG
		cout<< "c2," << i << " : ";
		printZZX(new_poly);
#endif
		c2_basedT_polys[i] = new_poly;
	}

	return c2_basedT_polys;
}




ZZX HEcryptoZZX::sum_c2i_ai(const vector<ZZX> c2_baseT_polys, const vector<HEpubkeyZZX *> he_homorphismKeys){
	ZZX sum_c2;
	for(int i=0; i<he_homorphismKeys.size(); ++i){
		HEpubkeyZZX *he_pubkey = he_homorphismKeys[i];
		sum_c2 += MulMod_f(c2_baseT_polys[i], he_pubkey->get_a0(), _f);
		sum_c2=Mod_f(sum_c2, _f);
#ifdef VERBOSE
		printZZX("ai", he_pubkey->get_a0());
#endif
	}
	return sum_c2;
}

ZZX HEcryptoZZX::sum_c2i_bi(const vector<ZZX> c2_baseT_polys, const vector<HEpubkeyZZX *> he_homorphismKeys){
	ZZX sum_c2;
	for(int i=0; i<he_homorphismKeys.size(); ++i){
		HEpubkeyZZX *he_pubkey = he_homorphismKeys[i];
		sum_c2 += MulMod_f(c2_baseT_polys[i], he_pubkey->get_a1(), _f);
		sum_c2= Mod_f(sum_c2, _f);
#ifdef VERBOSE
		printZZX("bi", he_pubkey->get_a1());
#endif
	}
	return sum_c2;
}


/*
 * Relinearization
 */
void HEcryptoZZX::relinearization(ZZX &c0_relin, ZZX &c1_relin, const ZZX c0, const ZZX c1, const ZZX c2, const vector<HEpubkeyZZX *> he_homorphismKeys){
	vector<ZZX> c2_baseT_polys = convert_c2_intoBaseT_Polys(c2);
	ZZX _sum_c2iai = sum_c2i_ai(c2_baseT_polys, he_homorphismKeys);
	ZZX _sum_c2ibi = sum_c2i_bi(c2_baseT_polys, he_homorphismKeys);

	c0_relin = c0 + _sum_c2iai;
	c0_relin = Mod_f(c0_relin, _f);

	c1_relin = c1 + _sum_c2ibi;
	c1_relin = Mod_f(c1_relin, _f);
}




/*
 * HE multiplication
 */
void HEcryptoZZX::HE_multiplication(ZZX &c0_relin, ZZX &c1_relin, const ZZX m, const ZZX m1, const vector<HEpubkeyZZX *> he_homorphismKeys){
	ZZX new_c0, new_c1, new_c2;
	multiplication(new_c0, new_c1, new_c2, m, m1);
	relinearization(c0_relin, c1_relin, new_c0, new_c1, new_c2, he_homorphismKeys);
}


void HEcryptoZZX::multiplication(ZZX &new_c0, ZZX &new_c1, ZZX &new_c2, const ZZX m, const ZZX m1){
	ZZX c0, c1;
	encrypt(c0, c1, m);

	ZZX c0_1, c1_1;
	encrypt(c0_1, c1_1, m1);

#ifdef DEBUG
	printZZX("c0", c0);
	printZZX("c1", c1);
	printZZX("c0_1", c0_1);
	printZZX("c1_1", c1_1);
#endif

	// new c0
	new_c0 = c0 * c0_1;
#ifdef VERBOSE
	printZZX("c0 * c0_1", new_c0);
#endif
	new_c0 = Mod_f(new_c0, _f);

	// new c1
	ZZX c0_c1_1 = c0 * c1_1;  // must (mod f and q) after every operation
#ifdef VERBOSE
	printZZX("c0 * c1_1", c0_c1_1);
#endif
	ZZX c0_1_c1 = c0_1 * c1;  // otherwise, there will be error
#ifdef VERBOSE
	printZZX("c0_1 * c1", c0_1_c1);
#endif
	new_c1 = Mod_f(c0_c1_1, _f) + Mod_f(c0_1_c1, _f);
#ifdef VERBOSE
	printZZX("c0_c1_1 +  c0_1_c1", new_c1);
#endif
	new_c1 = Mod_f(new_c1, _f);

	// new c2
	new_c2 = c1 * c1_1;
#ifdef VERBOSE
	printZZX("c1 * c1_1", new_c2);
#endif
	new_c2 = Mod_f(new_c2, _f);


#ifdef DEBUG
	printZZX("c0_mult", new_c0);
	printZZX("c1_mult", new_c1);
	printZZX("c2_mult", new_c2);
#endif
}


const ZZX& HEcryptoZZX::get_g() const {
	return _g;
}

const ZZX& HEcryptoZZX::get_h() const {
	return _h;
}

const ZZX& HEcryptoZZX::get_u() const {
	return _u;
}


