/*
 * demo.h
 *
 *  Created on: Jul 25, 2015
 *      Author: ph
 */

#ifndef LNV2011_DEMO_H_
#define LNV2011_DEMO_H_

#include <string>

void LNV2011_test_HE_operations(PerformanceStats *pStats, timing *timer);

void LNV2011_retest_HE_operations(PerformanceStats *pStats, timing *timer, const std::string json_file);

void LNV2011_run_demo(Options *opts);

void writeToFile(const std::string content);


#endif /* LNV2011_DEMO_H_ */
