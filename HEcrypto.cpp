/*
 * HEcrypto.cpp
 *
 *  Created on: Jul 1, 2015
 *      Author: ph
 */

#include <common/distribution.h>
#include <common/libntl_addon.h>
#include <common/CRT_ops.h>
#include <common/util.h>

#include "LNV2011.h"
#include "HEcrypto.h"

#include <vector>

HEcrypto::HEcrypto(LNV2011Params *params, JSONParser *jsonParser){
	// just to initialize it with empty pubkey; we set a pubkey later
	_he_pubkey=new HEpubkey(params);
	_jsonParser = jsonParser;
	_params = params;

	// generate basic parameters f, e, q
	generateParameters(params, _jsonParser);

	if(_jsonParser->isLoadedFromFile()){
		_u = _jsonParser->readJSONZZ_pX("u");
		_g = _jsonParser->readJSONZZ_pX("g");
		_h = _jsonParser->readJSONZZ_pX("h");
	}else{
		_u = generatePolyWithCoeffFromGaussianNoise(_params->n, _params->mean, _params->standard_deviation);
		_g = generatePolyWithCoeffFromGaussianNoise(_params->n, _params->mean, _params->standard_deviation);
		_h = generatePolyWithCoeffFromGaussianNoise(_params->n, _params->mean, _params->standard_deviation);
	}

	_modulis = getListOfModulis();
}


HEcrypto::~HEcrypto() {

}

void HEcrypto::printParameters(){
	BaseCrypto::printParameters(_jsonParser);
	printZZ_pX("u", _u);
	_jsonParser->addZZ_pXToJSON("u", _u);
	printZZ_pX("g", _g);
	_jsonParser->addZZ_pXToJSON("g", _g);
	printZZ_pX("h", _h);
	_jsonParser->addZZ_pXToJSON("h", _h);
}


void HEcrypto::encrypt(ZZ_pX &c0, ZZ_pX &c1, const ZZ_pX m){
	if(deg(_he_pubkey->get_a0())==-1) {
		cout << "error! add public key to the crypto object" << endl;
		exit(0);
	}

	// compute c0
	ZZ_pX a0_u = MulMod_f(_he_pubkey->get_a0(), _u, _phi);
	ZZ_pX t_g = MulMod_f(_params->t, _g, _phi);
	c0 = (a0_u +t_g + m);
	c0 = Mod_f(c0, _phi);

	// compute c1
	ZZ_pX a1_u = MulMod_f(_he_pubkey->get_a1(), _u, _phi);
	ZZ_pX t_h = MulMod_f(_params->t, _h, _phi);
	c1 = a1_u + t_h;
	c1 = Mod_f(c1, _phi);
}

void HEcrypto::encrypt(ZZ_pX &c0, ZZ_pX &c1, const ZZ_pX m, const ZZ t, const HEpubkey * pk){
	if(deg(pk->get_a0())==-1) {
		cout << "error! add public key to the crypto object" << endl;
		exit(0);
	}

	// compute c0
	ZZ_pX a0_u = MulMod_f(pk->get_a0(), _u, _phi);
	ZZ_pX t_g = MulMod_f(t, _g, _phi);
	c0 = (a0_u +t_g + m);
	c0 = Mod_f(c0, _phi);

	// compute c1
	ZZ_pX a1_u = MulMod_f(pk->get_a1(), _u, _phi);
	ZZ_pX t_h = MulMod_f(t, _h, _phi);
	c1 = a1_u + t_h;
	c1 = Mod_f(c1, _phi);
}


HEctxt *HEcrypto::encrypt(const ZZ_pX m){
	ZZ_pX c0, c1;
	encrypt(c0, c1, m);
	return new HEctxt(c0, c1);
}

vector<HEctxt *> HEcrypto::encrypt_vector(const vector<ZZ_pX> m_vector, const HEpubkey *pk){
	vector<HEctxt *> ctxt_vector;
	ZZ_pX c0, c1;
	for(int i=0;i<m_vector.size();i++){
		encrypt(c0, c1, m_vector[i], _params->t, pk);
		ctxt_vector.push_back(new HEctxt(c0, c1));
	}
	return ctxt_vector;
}


vector<HEctxt *> HEcrypto::encrypt_vector(const vector<ZZ_pX> m_vector, const vector<HEpubkey *> pk_vector){
	vector<HEctxt *> ctxt_vector;
	ZZ_pX c0, c1;
	for(int i=0;i<m_vector.size();i++){
		encrypt(c0, c1, m_vector[i], to_ZZ(_modulis[i]), pk_vector[i]);
		ctxt_vector.push_back(new HEctxt(c0, c1));
	}
	return ctxt_vector;
}


void HEcrypto::decrypt(ZZ_pX &m, const ZZ_pX c0, const ZZ_pX c1, const ZZ_pX s){
	ZZ_pX c1_s = MulMod_f(c1, s, _phi);
	m = c0 + c1_s;
	m = Mod_f(m, _phi);
	ZZX tmp_m = to_ZZX(m);
	modQhalf(tmp_m, _q);
	modT(tmp_m, _params->t);
	m=to_ZZ_pX(tmp_m);
}

void HEcrypto::decrypt(ZZ_pX &m, HEctxt *ctxt, const ZZ_pX s){
	decrypt(m, ctxt->getC0(), ctxt->getC1(), s);
}

void HEcrypto::decryptWithoutModT(ZZX &m, HEctxt *ctxt, const ZZ_pX s){
	ZZ_pX c1_s = MulMod_f(ctxt->getC1(), s, _phi);
	ZZ_pX tmp_m = ctxt->getC0() + c1_s;
	tmp_m = Mod_f(tmp_m, _phi);
	m = to_ZZX(tmp_m);
	modQhalf(m, _q);
}

void HEcrypto::decrypt(ZZ_pX &m, const ZZ_pX c0, const ZZ_pX c1, const ZZ_pX c2, const ZZ_pX s){
	ZZ_pX c1_s = c1 * s;   		  // c1 * s
#ifdef VERBOSE
	printZZ_pX("c1 * s", c1_s);
#endif
	c1_s= Mod_f(c1_s, _phi);   // m = c0 + c1*s
#ifdef VERBOSE
	printZZ_pX("c1 * s", c1_s);
#endif

	ZZ_pX ss = s * s;				  // s*s
#ifdef VERBOSE
	printZZ_pX("s * s", ss);
#endif
	ss = Mod_f(ss, _phi);
#ifdef VERBOSE
	printZZ_pX("s * s", ss);
#endif

	ZZ_pX c2_ss = c2 * ss;		  // c2 * ss
#ifdef VERBOSE
	printZZ_pX("c2 * ss", c2_ss);
#endif
	c2_ss= Mod_f(c2_ss, _phi);  // m = c0 + c1*s + c2*s*s
#ifdef VERBOSE
	printZZ_pX("c2 * ss", c2_ss);
#endif

	m = c0 + c1_s + c2_ss;

#ifdef VERBOSE
	printZZ_pX("mm", m);
#endif
	m = Mod_f(m, _phi);
#ifdef VERBOSE
	printZZ_pX("mm", m);
#endif

	ZZX tmp_m = to_ZZX(m);
	modQhalf(tmp_m, _q);
	modT(tmp_m, _params->t);
	m=to_ZZ_pX(tmp_m);
}

void HEcrypto::decrypt_withoutModT(ZZX &m, const ZZ_pX c0, const ZZ_pX c1, const ZZ_pX c2, const ZZ_pX s){
	ZZ_pX c1_s = c1 * s;   		  // c1 * s
#ifdef VERBOSE
	printZZ_pX("c1 * s", c1_s);
#endif
	c1_s= Mod_f(c1_s, _phi);   // m = c0 + c1*s
#ifdef VERBOSE
	printZZ_pX("c1 * s", c1_s);
#endif

	ZZ_pX ss = s * s;				  // s*s
#ifdef VERBOSE
	printZZ_pX("s * s", ss);
#endif
	ss = Mod_f(ss, _phi);
#ifdef VERBOSE
	printZZ_pX("s * s", ss);
#endif

	ZZ_pX c2_ss = c2 * ss;		  // c2 * ss
#ifdef VERBOSE
	printZZ_pX("c2 * ss", c2_ss);
#endif
	c2_ss= Mod_f(c2_ss, _phi);  // m = c0 + c1*s + c2*s*s
#ifdef VERBOSE
	printZZ_pX("c2 * ss", c2_ss);
#endif

	ZZ_pX tmp_m  = c0 + c1_s + c2_ss;

#ifdef VERBOSE
	printZZ_pX("mm", m);
#endif
	tmp_m = Mod_f(tmp_m, _phi);
#ifdef VERBOSE
	printZZ_pX("mm", m);
#endif

	m=to_ZZX(tmp_m);
	modQhalf(m, _q);
}


vector<ZZ_pX> HEcrypto::decrypt_vector_different_t(const vector<HEctxt *> m_ctxt_vector, const ZZ_pX s){
	vector<ZZ_pX> results;
	ZZX m_tmp;
	HEctxt *ctxt_tmp;
	for(int i=0;i<m_ctxt_vector.size();i++){
		ctxt_tmp=m_ctxt_vector[i];
		decrypt_withoutModT(m_tmp, ctxt_tmp->getC0(), ctxt_tmp->getC1(), ctxt_tmp->getC2(), s);
		modT(m_tmp, to_ZZ(_modulis[i]));
		results.push_back(to_ZZ_pX(m_tmp));
	}
	return results;
}

vector<ZZ_pX> HEcrypto::decrypt_vector_same_t(const vector<HEctxt *> m_ctxt_vector, const ZZ_pX s){
	vector<ZZ_pX> results;
	ZZX m_tmp;
	HEctxt *ctxt_tmp;
	for(int i=0;i<m_ctxt_vector.size();i++){
		ctxt_tmp=m_ctxt_vector[i];
		decrypt_withoutModT(m_tmp, ctxt_tmp->getC0(), ctxt_tmp->getC1(), ctxt_tmp->getC2(), s);
		modT(m_tmp, _params->t);
		results.push_back(to_ZZ_pX(m_tmp));
	}
	return results;
}


/*
 * HE addition
 */
HEctxt *HEcrypto::addition(const HEctxt *m1_ctxt, const HEctxt *m2_ctxt){
	ZZ_pX c0_sum, c1_sum;
	c0_sum = m1_ctxt->getC0() + m2_ctxt->getC0();
	//c0_sum = Mod_f(c0_sum, _phi);

	c1_sum = m1_ctxt->getC1() + m2_ctxt->getC1();
	//c1_sum = Mod_f(c1_sum, _phi);

#ifdef DEBUG
	printZZ_pX("c0 + c0_1", c0_sum);
	printZZ_pX("c1 + c1_1", c1_sum);
#endif

	return new HEctxt(c0_sum, c1_sum);
}

HEctxt *HEcrypto::addition_threeElements(const HEctxt *m1_ctxt, const HEctxt *m2_ctxt){
	ZZ_pX c0_sum, c1_sum, c2_sum;
	c0_sum = m1_ctxt->getC0() + m2_ctxt->getC0();
	//c0_sum = Mod_f(c0_sum, _phi);

	c1_sum = m1_ctxt->getC1() + m2_ctxt->getC1();
	//c1_sum = Mod_f(c1_sum, _phi);

	c2_sum = m1_ctxt->getC2() + m2_ctxt->getC2();

#ifdef DEBUG
	printZZ_pX("c0 + c0_1", c0_sum);
	printZZ_pX("c1 + c1_1", c1_sum);
	printZZ_pX("c2 + c2_1", c2_sum);
#endif

	return new HEctxt(c0_sum, c1_sum, c2_sum);
}

vector<HEctxt *> HEcrypto::addition_vector(const vector<HEctxt *> m_ctxt_vector, const vector<HEctxt *> m1_ctxt_vector){
	vector<HEctxt *> sum_vector;
	assert(m_ctxt_vector.size()==m1_ctxt_vector.size());

	for(int i=0;i<m_ctxt_vector.size();i++){
		sum_vector.push_back(addition_threeElements(m_ctxt_vector[i], m1_ctxt_vector[i]));
	}
	return sum_vector;
}


/*
 * HE subtraction
 */
HEctxt *HEcrypto::subtraction(const HEctxt *m1_ctxt, const HEctxt *m2_ctxt){
	ZZ_pX c0_diff, c1_diff;
	c0_diff = m1_ctxt->getC0() - m2_ctxt->getC0();
	//c0_diff = Mod_f(c0_diff, _phi);

	c1_diff = m1_ctxt->getC1() - m2_ctxt->getC1();
	//c1_diff = Mod_f(c1_diff, _phi);

#ifdef DEBUG
	printZZ_pX("c0 - c0_1", c0_diff);
	printZZ_pX("c1 - c1_1", c1_diff);
#endif

	return new HEctxt(c0_diff, c1_diff);
}

HEctxt *HEcrypto::subtraction_threeElements(const HEctxt *m1_ctxt, const HEctxt *m2_ctxt){
	ZZ_pX c0_rem, c1_rem, c2_rem;
	c0_rem = m1_ctxt->getC0() - m2_ctxt->getC0();
	//c0_sum = Mod_f(c0_sum, _phi);

	c1_rem = m1_ctxt->getC1() - m2_ctxt->getC1();
	//c1_sum = Mod_f(c1_sum, _phi);

	c2_rem = m1_ctxt->getC2() - m2_ctxt->getC2();

#ifdef DEBUG
	printZZ_pX("c0 - c0_1", c0_rem);
	printZZ_pX("c1 - c1_1", c1_rem);
	printZZ_pX("c2 - c2_1", c2_rem);
#endif

	return new HEctxt(c0_rem, c1_rem, c2_rem);
}


vector<HEctxt *> HEcrypto::subtraction_vector(const vector<HEctxt *> m_ctxt_vector, const vector<HEctxt *> m1_ctxt_vector){
	vector<HEctxt *> sum_vector;
	assert(m_ctxt_vector.size()==m1_ctxt_vector.size());

	for(int i=0;i<m_ctxt_vector.size();i++){
		sum_vector.push_back(subtraction_threeElements(m_ctxt_vector[i], m1_ctxt_vector[i]));
	}
	return sum_vector;
}




vector<ZZ_pX> HEcrypto::convert_c2_intoBaseT_Polys(const ZZ_pX c2){
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
	vector<ZZ_pX> c2_basedT_polys(v_size);
	for(int i=0; i<v_size;i++){
		ZZ_pX new_poly;
		for(int j=0; j<=deg(c2); j++){
			vector<ZZ> coeff_baseT =  coeffs_baseT_vector[j];
			//			if(j==0 || j==deg(c2))    // for setting the first and the last term negative
			//				SetCoeff(new_poly, j, -coeff_baseT[c2_basedT_polys.size()-1-i]);
			//			else
			SetCoeff(new_poly, j, conv<ZZ_p>(coeff_baseT[c2_basedT_polys.size()-1-i]));
		}

#ifdef DEBUG
		cout<< "c2," << i << " : ";
		printZZ_pX(new_poly);
#endif
		c2_basedT_polys[i] = new_poly;
	}

	return c2_basedT_polys;
}




ZZ_pX HEcrypto::sum_c2i_ai(const vector<ZZ_pX> c2_baseT_polys, const vector<HEpubkey *> he_homorphismKeys){
	ZZ_pX sum_c2;
	for(int i=0; i<he_homorphismKeys.size(); ++i){
		HEpubkey *he_pubkey = he_homorphismKeys[i];
		sum_c2 += MulMod_f(c2_baseT_polys[i], he_pubkey->get_a0(), _phi);
		sum_c2=Mod_f(sum_c2, _phi);
#ifdef VERBOSE
		printZZ_pX("ai", he_pubkey->get_a0());
#endif
	}
	return sum_c2;
}

ZZ_pX HEcrypto::sum_c2i_bi(const vector<ZZ_pX> c2_baseT_polys, const vector<HEpubkey *> he_homorphismKeys){
	ZZ_pX sum_c2;
	for(int i=0; i<he_homorphismKeys.size(); ++i){
		HEpubkey *he_pubkey = he_homorphismKeys[i];
		sum_c2 += MulMod_f(c2_baseT_polys[i], he_pubkey->get_a1(), _phi);
		sum_c2= Mod_f(sum_c2, _phi);
#ifdef VERBOSE
		printZZ_pX("bi", he_pubkey->get_a1());
#endif
	}
	return sum_c2;
}


/*
 * Relinearization
 */
void HEcrypto::relinearization(ZZ_pX &c0_relin, ZZ_pX &c1_relin, const ZZ_pX c0, const ZZ_pX c1, const ZZ_pX c2, const vector<HEpubkey *> he_homorphismKeys){
	vector<ZZ_pX> c2_baseT_polys = convert_c2_intoBaseT_Polys(c2);
	ZZ_pX _sum_c2iai = sum_c2i_ai(c2_baseT_polys, he_homorphismKeys);
	ZZ_pX _sum_c2ibi = sum_c2i_bi(c2_baseT_polys, he_homorphismKeys);

	c0_relin = c0 + _sum_c2iai;
	c0_relin = Mod_f(c0_relin, _phi);

	c1_relin = c1 + _sum_c2ibi;
	c1_relin = Mod_f(c1_relin, _phi);
}




/*
 * HE multiplication
 */
void HEcrypto::HE_multiplication(ZZ_pX &c0_relin, ZZ_pX &c1_relin, const ZZ_pX m, const ZZ_pX m1, const vector<HEpubkey *> he_homorphismKeys){
	ZZ_pX new_c0, new_c1, new_c2;
	multiplication(new_c0, new_c1, new_c2, m, m1);
	relinearization(c0_relin, c1_relin, new_c0, new_c1, new_c2, he_homorphismKeys);
}



void HEcrypto::multiplication(ZZ_pX &new_c0, ZZ_pX &new_c1, ZZ_pX &new_c2, const ZZ_pX m, const ZZ_pX m1){
	HEctxt *m_ctxt = encrypt(m);
	HEctxt *m1_ctxt = encrypt(m1);

	multiplication(new_c0, new_c1, new_c2, m_ctxt, m1_ctxt);
	delete m_ctxt;
	delete m1_ctxt;
}


void HEcrypto::multiplication(ZZ_pX &new_c0, ZZ_pX &new_c1, ZZ_pX &new_c2, const HEctxt *m_ctxt, const HEctxt *m1_ctxt){
#ifdef DEBUG
	printZZ_pX("c0", c0);
	printZZ_pX("c1", c1);
	printZZ_pX("c0_1", c0_1);
	printZZ_pX("c1_1", c1_1);
#endif

	// new c0
	new_c0 = m_ctxt->getC0() * m1_ctxt->getC0();
#ifdef VERBOSE
	printZZ_pX("c0 * c0_1", new_c0);
#endif
	new_c0 = Mod_f(new_c0, _phi);

	// new c1
	ZZ_pX c0_c1_1 = m_ctxt->getC0() * m1_ctxt->getC1();  // must (mod f and q) after every operation
#ifdef VERBOSE
	printZZ_pX("c0 * c1_1", c0_c1_1);
#endif
	ZZ_pX c0_1_c1 = m1_ctxt->getC0() * m_ctxt->getC1();  // otherwise, there will be error
#ifdef VERBOSE
	printZZ_pX("c0_1 * c1", c0_1_c1);
#endif
	new_c1 = Mod_f(c0_c1_1, _phi) + Mod_f(c0_1_c1, _phi);
#ifdef VERBOSE
	printZZ_pX("c0_c1_1 +  c0_1_c1", new_c1);
#endif
	new_c1 = Mod_f(new_c1, _phi);

	// new c2
	new_c2 = m_ctxt->getC1() * m1_ctxt->getC1();
#ifdef VERBOSE
	printZZ_pX("c1 * c1_1", new_c2);
#endif
	new_c2 = Mod_f(new_c2, _phi);


#ifdef DEBUG
	printZZ_pX("c0_mult", new_c0);
	printZZ_pX("c1_mult", new_c1);
	printZZ_pX("c2_mult", new_c2);
#endif
}

HEctxt *HEcrypto::multiplication(const HEctxt *m_ctxt, const HEctxt *m1_ctxt){
	ZZ_pX new_c0, new_c1, new_c2;
	multiplication(new_c0, new_c1, new_c2, m_ctxt, m1_ctxt);
	return new HEctxt(new_c0, new_c1, new_c2);
}

HEctxt *HEcrypto::multiplication_with_relinearization(const HEctxt *m_ctxt, const HEctxt *m1_ctxt, const vector<HEpubkey *> he_homorphismKeys){
	ZZ_pX c0_relin, c1_relin;
	ZZ_pX new_c0, new_c1, new_c2;
	multiplication(new_c0, new_c1, new_c2, m_ctxt, m1_ctxt);
	relinearization(c0_relin, c1_relin, new_c0, new_c1, new_c2, he_homorphismKeys);

	return new HEctxt(c0_relin, c1_relin);
}


void HEcrypto::multiplication(ZZ_pX &new_c0, ZZ_pX &new_c1, const HEctxt *m_ctxt, const ZZ_pX m1){
	// new c0
	new_c0 = m_ctxt->getC0() * m1;
#ifdef VERBOSE
	printZZ_pX("c0 * m1", new_c0);
#endif
	new_c0 = Mod_f(new_c0, _phi);

	new_c1 = m_ctxt->getC1() * m1;
#ifdef VERBOSE
	printZZ_pX("c1 * m1", new_c1);
#endif
	new_c1 = Mod_f(new_c1, _phi);


#ifdef DEBUG
	printZZ_pX("c0_mult", new_c0);
	printZZ_pX("c1_mult", new_c1);
#endif
}


vector<HEctxt *> HEcrypto::multiplication_vector(const vector<HEctxt *> m_ctxt_vector, const vector<HEctxt *> m1_ctxt_vector, const vector<HEpubkey *> he_homorphismKeys){
	vector<HEctxt *> prod_vector;
	assert(m_ctxt_vector.size()==m1_ctxt_vector.size());

	for(int i=0;i<m_ctxt_vector.size();i++){
		prod_vector.push_back(multiplication_with_relinearization(m_ctxt_vector[i], m1_ctxt_vector[i], he_homorphismKeys));
	}

	return prod_vector;
}

vector<HEctxt *> HEcrypto::multiplication_vector(const vector<HEctxt *> m_ctxt_vector, const vector<HEctxt *> m1_ctxt_vector){
	vector<HEctxt *> prod_vector;
	assert(m_ctxt_vector.size()==m1_ctxt_vector.size());

	for(int i=0;i<m_ctxt_vector.size();i++){
		prod_vector.push_back(multiplication(m_ctxt_vector[i], m1_ctxt_vector[i]));
	}

	return prod_vector;
}

void HEcrypto::delete_HEctxt_vector(vector<HEctxt *> vec){
	for(int i=0;i<vec.size();i++){
		HEctxt *tmp = vec[i];
		tmp->~HEctxt();
		delete tmp;
	}
}


const ZZ_pX& HEcrypto::get_g() const {
	return _g;
}

const ZZ_pX& HEcrypto::get_h() const {
	return _h;
}

const ZZ_pX& HEcrypto::get_u() const {
	return _u;
}


