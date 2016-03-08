/*
 * demo.cpp
 *
 *  Created on: Jul 25, 2015
 *      Author: ph
 */


/*
 * documentation:
 * http://www.shoup.net/ntl/doc/tour-struct.html
 * http://www.shoup.net/ntl/doc/ZZ_pX.cpp.html
 * http://www.shoup.net/ntl/doc/ZZ.cpp.html
 * http://www.cams.aub.edu.lb/docs/files/NTL/conversions.txt
 */

#include <common/distribution.h>
#include <common/JSONParser.h>
#include <common/libntl_addon.h>
#include <common/timing.h>
#include <common/profiler.h>
#include <common/util.h>

#include <fstream>
#include <iostream>
#include <ctime>
#include <math.h>


#include "HEcrypto.h"
#include "HEkeygen.h"
#include "HEpubkey.h"
#include "HEctxt.h"
#include "Fujitsu.h"
#include "app/GeoHashExample.h"

void writeToFile(const std::string content) {
	ofstream jsonfile;
	std::time_t current_timestamp = std::time(nullptr);

	const string _folder = "experiments";
	mkdir(_folder.c_str(), 0777);
	string filename = to_string(current_timestamp) + "_ctxt.txt";
	jsonfile.open(_folder + "/" + filename);

	jsonfile << content;
	jsonfile.close();
}

/*
 * load the parameters from an initial JSON file
 */
void load_parameters(LNV2011Params *params, const string params_filename){
	JSONParser *paramsParser = new JSONParser();
	paramsParser->loadJSONFromFile(params_filename);
	printHeader("Loading Parameters from "+params_filename);
	//	cout<< "test" << endl;
	//		paramsParser->dumpJSONTxt();
	//		cout<< "test" << endl;
	params->mean=paramsParser->readJSONINT("mean");
	printINT("mean", params->mean);
	params->standard_deviation=paramsParser->readJSONDOUBLE("standard_deviation");
	printDOUBLE("standard_deviation", params->standard_deviation);
	params->n=paramsParser->readJSONINT("n");
	printINT("n", params->n);
	params->q_length_in_bits=paramsParser->readJSONINT("q_length_in_bits");
	printINT("q_length_in_bits", params->q_length_in_bits);
	params->t=paramsParser->readJSONZZ("t");
	printZZ("t", params->t);
	delete paramsParser;
}

/*
 * export parameter to store it in record file
 */
void export_parameters(JSONParser *jsonParser, LNV2011Params *params){
	jsonParser->addINTToJSON("mean", params->mean);
	jsonParser->addDOUBLEToJSON("standard_deviation", params->standard_deviation);
	jsonParser->addZZToJSON("t", params->t);
	jsonParser->addINTToJSON("n", params->n);
	jsonParser->addINTToJSON("q_length_in_bits", params->q_length_in_bits);
}

void LNV2011_test_HE_operations(const string params_filename, PerformanceStats *pStats, timing *timer, bool exportParameters){
	//load the parameters
	LNV2011Params params;
	if(!params_filename.empty())
		load_parameters(&params, params_filename);

	JSONParser *jsonParser = new JSONParser();

	// instantiate a new HEcrypto object
	HEcrypto* he_crypto = new HEcrypto(&params, jsonParser);

	// pass HEcrypto object for accessing f, q, e
	HEkeygen *he_keygen = new HEkeygen(&params, he_crypto, jsonParser);

	/*
	 * sample s from GaussianNoise
	 */
	ZZ_pX param_s;
	while(param_s==0){
		timer->start();
		param_s=generatePolyWithCoeffFromGaussianNoise(params.n, params.mean, params.standard_deviation);
		timer->stop("sk_sample", false);
	}
	printZZ_pX("s", param_s);
	jsonParser->addZZ_pXToJSON("s", param_s);

	/*
	 * compute a0, a1
	 */
	timer->start();
	HEpubkey *he_pubkey = he_keygen->createHEpubkey(param_s);
	timer->stop("pubkey_gen", false);
	he_crypto->addHEpubkey(he_pubkey);

	printZZ_pX("a0", he_pubkey->get_a0());
	jsonParser->addZZ_pXToJSON("a0", he_pubkey->get_a0());
	printZZ_pX("a1", he_pubkey->get_a1());
	jsonParser->addZZ_pXToJSON("a1", he_pubkey->get_a1());


	printHeader("HEcrypto parameters");
	he_crypto->printParameters();


	printHeader("Basic parameters");
	export_parameters(jsonParser, &params);
	printZZ("t", params.t);
	printINT("n=2^k", params.n);

	ZZ_pX m, m1;
	while(m==0)
		m = generatePolyWithCoeffFromRq(params.n);
	modT(m, params.t);

	printZZ_pX("m", m);
	jsonParser->addZZ_pXToJSON("m", m);
	while(m1==0)
		m1 = generatePolyWithCoeffFromRq(params.n);
	modT(m1, params.t);
	printZZ_pX("m1", m1);
	jsonParser->addZZ_pXToJSON("m1", m1);


	/*
	 * test case for encrypt and decrypt m
	 */
	printHeader("Test encrypt/decrypt");
	// encryption
	timer->start();
	HEctxt *m_ctxt = he_crypto->encrypt(m);
	timer->stop("encryption", false);

	string text;
	text = "c0: " + returnStringFromZZ_pX(m_ctxt->getC0()) + " c1: " + returnStringFromZZ_pX(m_ctxt->getC0());
	writeToFile(text);

	printZZ_pX("c0", m_ctxt->getC0());
	jsonParser->addZZ_pXToJSON("c0", m_ctxt->getC0());
	printZZ_pX("c1", m_ctxt->getC1());
	jsonParser->addZZ_pXToJSON("c1", m_ctxt->getC1());

	// decryption
	ZZ_pX m_dec;
	timer->start();
	he_crypto->decrypt(m_dec, m_ctxt, param_s);
	timer->stop("decryption_deg1", false);
	printZZ_pX("m_dec", m_dec);
	if(m==m_dec){
		pStats->encrypt_successed++;
		pStats->success=true;
	}else{
		cerr << "encryption failed! -> " << endl;
		printZZ_pX("m", m);
		pStats->encrypt_failed++;
		pStats->success=false;
	}


	/*
	 * addition of m and m1
	 */
	printHeader("Test Addition");
	ZZ_pX add_m_m1_plain = m + m1;
	modT(add_m_m1_plain, params.t);
	printZZ_pX("m+m1", add_m_m1_plain);

	ZZ_pX add_m_m1_dec;
	HEctxt *m1_ctxt = he_crypto->encrypt(m1);
	timer->start();
	HEctxt *sum_ctxt = he_crypto->addition(m_ctxt, m1_ctxt);
	timer->stop("addition", false);
	he_crypto->decrypt(add_m_m1_dec, sum_ctxt, param_s);
	printZZ_pX("dec(enc(m)+enc(m1))", add_m_m1_dec);
	if(add_m_m1_plain == add_m_m1_dec){
		pStats->addition_successed++;
		pStats->success=true;
	}else{
		cerr << "Addition failed! -> " << endl;
		printZZ_pX("add_m_m1_plain", add_m_m1_plain);
		pStats->addition_failed++;
		pStats->success=false;
	}

	delete sum_ctxt;

	/*
	 * subtraction of m and m1
	 */
	printHeader("Test Subtraction");
	ZZX diff_m_m1_plain = to_ZZX(m) - to_ZZX(m1);
	modT(diff_m_m1_plain, params.t);
	ZZ_pX diff_m_m1_plain_zz_pX = to_ZZ_pX(diff_m_m1_plain);
	printZZ_pX("m-m1", diff_m_m1_plain_zz_pX);

	ZZ_pX diff_m_m1_dec;
	timer->start();
	HEctxt *diff_ctxt = he_crypto->subtraction(m_ctxt, m1_ctxt);
	timer->stop("subtraction", false);
	he_crypto->decrypt(diff_m_m1_dec, diff_ctxt, param_s);
	printZZ_pX("dec(enc(m)-enc(m1))", diff_m_m1_dec);
	if (diff_m_m1_plain_zz_pX == diff_m_m1_dec) {
		pStats->subtraction_succeeded++;
		pStats->success = true;
	} else {
		cerr << "Subtraction failed! -> " << endl;
		printZZ_pX("diff_m_m1_plain", diff_m_m1_plain_zz_pX);
		pStats->subtraction_failed++;
		pStats->success = false;
	}
	delete diff_ctxt;

	printHeader("Test Multiplication");
	ZZ_pX prod_m_m1_plain = m * m1;
	modF(prod_m_m1_plain, he_crypto->get_phi());
	ZZX tmp_mm = to_ZZX(prod_m_m1_plain);
	modQhalf(tmp_mm, he_crypto->get_q());
	modT(tmp_mm, params.t);
	prod_m_m1_plain = to_ZZ_pX(tmp_mm);
	printZZ_pX("m*m1", prod_m_m1_plain);

	ZZ_pX new_c0, new_c1, new_c2, prod_m_m1_dec;
	timer->start();
	he_crypto->multiplication(new_c0, new_c1, new_c2, m, m1);
	timer->stop("multiplication", false);

	timer->start();
	he_crypto->decrypt(prod_m_m1_dec, new_c0, new_c1, new_c2, param_s);
	timer->stop("mult_decryption_deg2", false);
	printZZ_pX("dec(enc(m)*enc(m1))", prod_m_m1_dec);
	if(prod_m_m1_plain == prod_m_m1_dec){
		pStats->multiplication_successed++;
		pStats->success=true;
	}else{
		cerr << "Multiplication failed! -> " << endl;
		printZZ_pX("prod_m_m1_plain", prod_m_m1_plain);
		pStats->multiplication_failed++;
		pStats->success=false;
	}


	/*
	 * generate Homomorphism Keys
	 */
	timer->start();
	vector<HEpubkey *> he_homomorphismKeys = he_keygen->createHomomorphismKeys(param_s);
	timer->stop("hKeygen", false);
	jsonParser->addHomomorphismKeysToJSON(he_keygen->convertHomomorphismKeysToJSONObj(he_homomorphismKeys));

	/*
	 * multiplication of m and m1
	 */
	printHeader("Test Multiplication with relinearization");
	ZZ_pX c0_relin, c1_relin, prod_m_m1_dec_relin;
	timer->start();
	he_crypto->HE_multiplication(c0_relin, c1_relin, m, m1, he_homomorphismKeys);
	timer->stop("multiplication_relin", false);
	timer->start();
	he_crypto->decrypt(prod_m_m1_dec_relin, c0_relin, c1_relin, param_s);
	timer->stop("multrelin_decryption_deg2", false);
	printZZ_pX("dec(enc(m)*enc(m1))", prod_m_m1_dec_relin);
	if(prod_m_m1_plain == prod_m_m1_dec_relin){
		pStats->multiplication_relin_successed++;
		pStats->success=true;
	}else{
		cerr << "Multiplication relin failed! -> " << endl;
		printZZ_pX("prod_m_m1_plain", prod_m_m1_plain);
		pStats->multiplication_relin_failed++;
		pStats->success=false;
	}

	//printINT("MEM (k)", getVirtualMemValue());

	// applications
	printHeader("Hamming Distance");
	//runFujitsuHammingDistanceExample(&params, he_crypto, he_keygen, he_pubkey, param_s);
	cout<<endl;

	printHeader("Euclidean Distance");
	//runFujitsuEuclideanDistanceExample(&params, he_crypto, he_keygen, he_pubkey, param_s);
	cout<<endl;


	printHeader("GeoHashing example");
//	runGeoHashingExample();
//	cout<<endl;

	//runGeoHashTest(&params, he_crypto, he_keygen, he_pubkey, param_s);
	//runGeoHashTest2(&params, he_crypto, he_keygen, he_pubkey, param_s, timer);
	//cout<<endl;


	printHeader("Euclidean distance example");
	//runEuclideanDistanceTest(&params, he_crypto, he_keygen, param_s, timer);
	cout<<endl;




	if(exportParameters){
		string json_file_prefix = "t"+to_string(to_long(params.t))+"n"+to_string(params.n)+"q"+ZZToString(he_crypto->get_q());
		jsonParser->dumpJSONToFile(json_file_prefix);
	}

	// delete objects
	//he_keygen->destoryhKeys(he_homomorphismKeys);  //TODO investigate what happen to this; it causes memory leaks
	delete m_ctxt;
	delete m1_ctxt;
	delete jsonParser;
	delete he_pubkey;
	delete he_keygen;
	delete he_crypto;
}




/*
 * reload parameters from the stored JSON record file
 */
void reload_parameters(JSONParser *jsonParser, LNV2011Params *params){
	params->mean=jsonParser->readJSONINT("mean");
	params->standard_deviation = jsonParser->readJSONDOUBLE("standard_deviation");
	params->t = jsonParser->readJSONZZ("t");

	params->n = jsonParser->readJSONINT("n");
	params->q_length_in_bits = jsonParser->readJSONINT("q_length_in_bits");
}


/*
 * This function is for re-testing the HE operations with the same parameters recorded in another run
 */
void LNV2011_retest_HE_operations(PerformanceStats *pStats, timing *timer, const string json_file){
	LNV2011Params params;
	printHeader("recompute "+json_file);
	JSONParser *jsonParser = new JSONParser();

	jsonParser->loadJSONFromFile(json_file);
	reload_parameters(jsonParser, &params);

	// instantiate a new HEcrypto object
	HEcrypto* he_crypto = new HEcrypto(&params, jsonParser);

	// pass HEcrypto object for accessing f, q, e
	HEkeygen *he_keygen = new HEkeygen(&params, he_crypto, jsonParser);

	/*
	 * sample s from GaussianNoise
	 */
	ZZ_pX s = jsonParser->readJSONZZ_pX("s");
	printZZ_pX("s", s);

	/*
	 * compute a0, a1
	 */
	printHeader("Compute a0, a1");
	timer->start();
	HEpubkey *he_pubkey = he_keygen->createHEpubkey(s);
	timer->stop("pubkeygen", false);
	printZZ_pX("a0", he_pubkey->get_a0());
	printZZ_pX("a1", he_pubkey->get_a1());
	he_crypto->addHEpubkey(he_pubkey);


	printHeader("HEcrypto parameters");
	he_crypto->printParameters();


	printHeader("Basic parameters");
	printZZ("t", params.t);
	printINT("n=2^k", params.n);

	ZZ_pX m = jsonParser->readJSONZZ_pX("m");
	printZZ_pX("m", m);
	ZZ_pX m1 = jsonParser->readJSONZZ_pX("m1");
	printZZ_pX("m1", m1);


	/*
	 * test case for encrypt and decrypt m
	 */
	printHeader("Test encrypt/decrypt");
	// encryption
	timer->start();
	HEctxt *m_ctxt = he_crypto->encrypt(m);
	timer->stop("encryption", false);
	printZZ_pX("c0", m_ctxt->getC0());
	printZZ_pX("c1", m_ctxt->getC1());

	// decryption
	ZZ_pX m_dec;
	timer->start();
	he_crypto->decrypt(m_dec, m_ctxt, s);
	timer->stop("decryption", false);
	printZZ_pX("m_dec", m_dec);
	if(m==m_dec){
		pStats->encrypt_successed++;
		pStats->success=true;
	}else{
		cerr << "encryption failed! -> " << endl;
		printZZ_pX("m", m);
		pStats->encrypt_failed++;
		pStats->success=false;
	}


	/*
	 * addition of m and m1
	 */
	printHeader("Test Addition");
	ZZ_pX add_m_m1_plain = m + m1;
	modT(add_m_m1_plain, params.t);
	printZZ_pX("m+m1", add_m_m1_plain);

	HEctxt *m1_ctxt = he_crypto->encrypt(m1);

	ZZ_pX add_m_m1_dec;
	timer->start();
	HEctxt *sum_ctxt = he_crypto->addition(m_ctxt, m1_ctxt);
	timer->stop("addition", false);
	he_crypto->decrypt(add_m_m1_dec, sum_ctxt, s);
	printZZ_pX("dec(enc(m)+enc(m1))", add_m_m1_dec);
	if(add_m_m1_plain == add_m_m1_dec){
		pStats->addition_successed++;
		pStats->success=true;
	}else{
		cerr << "Addition failed! -> " << endl;
		printZZ_pX("add_m_m1_plain", add_m_m1_plain);
		pStats->addition_failed++;
		pStats->success=false;
	}

	delete sum_ctxt;

	printHeader("Test Multiplication");
	ZZ_pX prod_m_m1_plain = m * m1;
	modF(prod_m_m1_plain, he_crypto->get_phi());
	ZZX tmp_mm = to_ZZX(prod_m_m1_plain);
	modQhalf(tmp_mm, he_crypto->get_q());
	modT(tmp_mm, params.t);
	prod_m_m1_plain = to_ZZ_pX(tmp_mm);
	printZZ_pX("m*m1", prod_m_m1_plain);

	ZZ_pX new_c0, new_c1, new_c2, prod_m_m1_dec;
	timer->start();
	he_crypto->multiplication(new_c0, new_c1, new_c2, m, m1);
	timer->stop("multiplication", false);
	timer->start();
	he_crypto->decrypt(prod_m_m1_dec, new_c0, new_c1, new_c2, s);
	timer->stop("mult_decryption_deg2", false);
	printZZ_pX("dec(enc(m)*enc(m1))", prod_m_m1_dec);
	if(prod_m_m1_plain == prod_m_m1_dec){
		pStats->multiplication_successed++;
		pStats->success=true;
	}else{
		cerr << "Multiplication failed! -> " << endl;
		printZZ_pX("prod_m_m1_plain", prod_m_m1_plain);
		pStats->multiplication_failed++;
		pStats->success=false;
	}


	/*
	 * get Homomorphism Keys
	 */
	vector<HEpubkey *> he_homomorphismKeys = he_keygen->reconstructHomomorphismKeysFromJSONObj(jsonParser->getHomomorphismKeysJSONObj());

	/*
	 * multiplication of m and m1
	 */
	printHeader("Test Multiplication with relinearization");
	ZZ_pX c0_relin, c1_relin, prod_m_m1_dec_relin;
	timer->start();
	he_crypto->HE_multiplication(c0_relin, c1_relin, m, m1, he_homomorphismKeys);
	timer->stop("multiplication_relin", false);
	timer->start();
	he_crypto->decrypt(prod_m_m1_dec_relin, c0_relin, c1_relin, s);
	timer->stop("multrelin_decryption_deg2", false);
	printZZ_pX("dec(enc(m)*enc(m1))", prod_m_m1_dec_relin);
	if(prod_m_m1_plain == prod_m_m1_dec_relin){
		pStats->multiplication_relin_successed++;
		pStats->success=true;
	}else{
		cerr << "Multiplication relin failed! -> " << endl;
		printZZ_pX("prod_m_m1_plain", prod_m_m1_plain);
		pStats->multiplication_relin_failed++;
		pStats->success=false;
	}

	// delete objects
	he_keygen->destoryhKeys(he_homomorphismKeys);
	delete m_ctxt;
	delete m1_ctxt;
	delete jsonParser;
	delete he_pubkey;
	delete he_crypto;
	delete he_keygen;
}

void terminal_test(const string params_filename, PerformanceStats *pStats, timing *timer, float run_count){
	float progress=0;
	for(int i=0;i<run_count;i++){
		LNV2011_test_HE_operations(params_filename, pStats, timer, false);
		progress=(float)i/(float)run_count * 100;
		if(fmod(progress, 10) ==0)
			cout<<"Progress "<< progress <<"%/100%" << endl;
	}
}


void LNV2011_run_demo(Options *opts){
	PerformanceStats pStats;
	timing timer;

	switch(opts->rType){
	case TERMINAL:
		terminal_test(opts->scheme_params_filename, &pStats, &timer, opts->terminal_run_count);
		break;
	case RERUN:
		if(fileExists(opts->rerun_json_filename)){
			LNV2011_retest_HE_operations(&pStats, &timer, opts->rerun_json_filename);
		}else{
			cerr<<"invalid file or file does not exist!"<<endl;
			abort();
		}
		break;
	case ECLIPSE:
	default:
		LNV2011_test_HE_operations(opts->scheme_params_filename, &pStats, &timer, true);
		break;
	}

	printPStats(&pStats);
	cout<<endl;
	timer.show();
}
