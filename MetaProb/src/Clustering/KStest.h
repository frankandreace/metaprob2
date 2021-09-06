/*
 * KStest.h
 *
 *  Created on: 07/mar/2016
 *      Author: samuele
 */

#ifndef SRC_CLUSTERING_KSTEST_H_
#define SRC_CLUSTERING_KSTEST_H_

#include <cmath>
#include <algorithm>
#include <boost/math/distributions/normal.hpp>
#include "nr3.h"
using namespace std;

// NR-Struktur fuer KS-Verteilung
//Approximate
struct KSdist {
    // Kap 14.3.3. S.737 + Kap 6.14.56 S. 334

    double pks(double z) {
        // Return cumulative distribution function.
        if (z < 0.) throw ("bad z in KSdist");
        if (z == 0.) return 0.;
        if (z < 1.18) {
            double y = exp(-1.23370055013616983 / SQR(z));
            return 2.25675833419102515 * sqrt(-log(y))*(y + pow(y, 9) + pow(y, 25) + pow(y, 49));
        } else {
            double x = exp(-2. * SQR(z));
            return 1. - 2. * (x - pow(x, 4) + pow(x, 9));
        }
    }

    // Return complementary cumulative distribution function.

    double qks(double z) {
        if (z < 0.) throw ("bad z in KSdist");
        if (z == 0.) return 1.;
        if (z < 1.18) return 1. - pks(z);
        double x = exp(-2. * SQR(z));
        return 2. * (x - pow(x, 4) + pow(x, 9));
    }

    // more code in the book ...
};

// NR-Funkton zum KS-Test
void kstwo(vector<long double> &data1, vector<long double> &data2, long double &d, long double &prob);
void ksone(vector<long double> &data, long double func(const long double), long double &d, long double &prob);
long double probks(const long double alam); //Approximate to n sum
long double normalCDF(const long double);

#endif /* SRC_CLUSTERING_KSTEST_H_ */
