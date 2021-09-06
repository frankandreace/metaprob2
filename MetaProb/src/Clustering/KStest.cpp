/*
 * KStest.cpp
 *
 *  Created on: 07/mar/2016
 *      Author: samuele
 */

#include "KStest.h"
// NR-Funkton zum KS-Test
void kstwo(vector<long double> &data1, vector<long double> &data2, long double &d, long double &prob) {
    // Given an array data1[0..n1-1], and an array data2[0..n2-1], this routine returns the K–S
    // statistic d and the p-value prob for the null hypothesis that the data sets are drawn from the
    // same distribution. Small values of prob show that the cumulative distribution function of data1
    // is significantly different from that of data2. The arrays data1 and data2 are modified by being
    // sorted into ascending order.
    int j1 = 0, j2 = 0, n1 = data1.size(), n2 = data2.size();
    double d1, d2, dt, en1, en2, en, fn1 = 0.0, fn2 = 0.0;
    KSdist ks;
    sort(data1.begin(), data1.end());
    sort(data2.begin(), data2.end());

    en1 = n1;
    en2 = n2;
    d = 0.0;
    while (j1 < n1 && j2 < n2) // If we are not done...
    {
        if ((d1 = data1[j1]) <= (d2 = data2[j2])) // Next step is in data1.
        {
            do {
                fn1 = ++j1 / en1;
            } while (j1 < n1 && d1 == data1[j1]);
        }
        if (d2 <= d1) // Next step is in data2.
        {
            do {
                fn2 = ++j2 / en2;
            } while (j2 < n2 && d2 == data2[j2]);
        }
        if ((dt = abs(fn2 - fn1)) > d) {
            //cout << " fn2 " << fn2 << "; fn1 " << fn1 << "; d " << d << endl;
            d = dt;
        }
    }
    en = sqrt(en1 * en2 / (en1 + en2));
    //cout << "* en " << en << " aus en1 " << en1 << " und en2 " << en2 << endl;
    prob = probks((en + 0.12 + 0.11 / en) * d);//ks.qks((en + 0.12 + 0.11 / en) * d); // Compute p-value.
}

void ksone(vector<long double> &data, long double func(const long double), long double &d, long double &prob)
{
	int j;
	long double dt,en,ff,fn,fo=0.0;

	int n=data.size();
	sort(data.begin(), data.end());
	en=n;
	d=0.0;
	for (j=0;j<n;j++) {
		fn=(j+1)/en;
		ff=func(data[j]);
		dt=MAX(fabs(fo-ff),fabs(fn-ff));
		if (dt > d) d=dt;
		fo=fn;
	}
	en=sqrt(en);
	prob=probks((en+0.12+0.11/en)*d);
}

long double probks(const long double alam) //Al posto di ks.qks
{
	const long double EPS1=1.0e-6,EPS2=1.0e-16;
	int j;
	long double a2,fac=2.0,sum=0.0,term,termbf=0.0;

	a2 = -2.0*alam*alam;
	for (j=1;j<=100;j++) {
		term=fac*exp(a2*j*j);
		sum += term;
		if (fabs(term) <= EPS1*termbf || fabs(term) <= EPS2*sum) return sum;
		fac = -fac;
		termbf=fabs(term);
	}
	return 1.0;
}

long double normalCDF(const long double val)
{
	boost::math::normal_distribution<long double> normal;
	return boost::math::cdf(normal, val);
}
