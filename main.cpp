#include<stdio.h>
#include<stdlib.h>
#include<NTL/matrix.h>
#include<NTL/ZZ.h>
#include<NTL/ZZ_p.h>
#include<NTL/RR.h>
#include<NTL/vec_ZZ.h>
#include<NTL/vec_RR.h>
#include<NTL/mat_ZZ.h>
#include<NTL/mat_ZZ_p.h>
#include<NTL/mat_RR.h>
#include<NTL/LLL.h>
#include<NTL/HNF.h>
#include<iostream>
#include<time.h>
#include<math.h>
#include<fstream>
#include<time.h>
#include<cmath>
#include<iomanip>
#include <windows.h>; 
#include <stdio.h>;
#include <random>



NTL_CLIENT;
using namespace std;

RR PI = to_RR(3.14159265358979323846264338327950288419716939937510582097494);


int sp;//smoothing parameter
int t;
int precision;

ofstream out("result.txt");
ofstream out2("result2.txt");

int outputfile = 1;

int EX1 = 1; //basic

int EX2 = 1; //PDG14

int EX3 = 0; //a=1-11,b=1 

int EX4 = 1; // Revised PDG14

int EX5 = 1; //Revised MW17


RR Gfun(RR x, RR c, RR s) {
	return exp(-1 * PI* pow(abs((x - c) / s), to_RR(2)));
}

 
void init(int sp_t, int t_t, int precision_t) {
	t = t_t;
	sp = sp_t;
	precision = precision_t;
	PI = ComputePi_RR();
}


void set_pre(int p) {
	precision = p;
	PI = ComputePi_RR();
}

RR* init_Gaussian(RR s, int * ra) {

	int range;
	conv(range, ceil(s*t));
	//cout << range << endl;
	RR* gau = new RR[2 * range+1];
	memset(gau, 0, (2 * range+1) * sizeof(RR));
	RR trans;
	if (precision > 0) {
		trans.SetPrecision(precision);
	}

	RR sum = to_RR(0);
	for (int i = 0; i <= 2 * range; i++) {

		RR tmp= Gfun(to_RR(i), to_RR(range), s);	 
		trans = to_RR(tmp);
		tmp = trans;

		long tl = NumBits(tmp.x) + tmp.e;
		if (tl > (-1 * tmp.precision())) {
			sum += tmp;
			gau[i] = tmp;
		}
		else
			gau[i] = 0;
	}
 
	//cout << sum << " " << log(sum) / log(2) << endl;
	//system("pause");
	for (int i = 0; i <= 2 * range; i++) {
		/*RR tmp = Gfun(to_RR(i), to_RR(range), s);	 
		trans = to_RR(tmp);
		tmp = trans;		
		long tl =   NumBits(tmp.x) + tmp.e;*/
		//cout << tl << " 111 " << (tl >(-1 * tmp.precision()))<<endl;
		//if (tl > (-1 * tmp.precision()))
		if(gau[i]!=0)
			gau[i] = gau[i] / sum;
		else
			gau[i] = 0;
	 
	}

	ra[0] = 2 * range;
	return gau;
}



void compare(RR* gc, RR* gau, int rc) {
	//compare

	RR trans;	 
	if (precision > 0) {
		trans.SetPrecision(precision);
	}


	RR d_sd = to_RR(0), d_kl = to_RR(0), d_ml = to_RR(0), d_re = to_RR(0);

	int ren = -1;
	int mln = -1;

	int ren2 = -1;
	int mln2 = -1;


	int cnt = 0;

	for (int i = 0; i <= rc; i++) {

	 
		//out2   << i << " " << log(abs(gc[i] - gau[i])) / log(2)<<" "<< gau[i]<< " " << gc[i] <<"  "<<abs(i-rc/2)<< endl;

		if (gc[i] > 0) {
		//	cout << gau[i].x << " " << gau[i].e << " " << i << " " << gau[i] << endl;
			long tl = NumBits(gc[i].x) + gc[i].e;
		
			if (tl > (-1 * gc[i].precision())  ) {
				d_sd += abs(gc[i] - gau[i]);
			}

			if (abs(gc[i] - gau[i]) / gc[i]> d_re&&tl > (-1 * gc[i].precision())&&i<=rc/2) {
				//cout << log(abs(gc[i] - gau[i])) / log(2) << " " << log(gc[i]) / log(2) << " " << gc[i].precision() << " " << i << " " << rc <<" "<< log(abs(abs(gc[i] - gau[i]) / gc[i] - d_re)) / log(2) << endl;
				d_re = abs(gc[i] - gau[i]) / gc[i];
				ren = i; 
				ren2 = rc-i;
			}

			if (gau[i] > 0) {
				d_kl += log(gc[i] / gau[i])*gc[i];
				if (abs(log(gc[i] / gau[i])) > d_ml&&tl > (-1 * gc[i].precision()) && i <= rc / 2) {
					d_ml = abs(log(gc[i] / gau[i]));
					mln = i;
					mln2 = rc-i;
				}				 
			}
		}

	}
	d_sd=d_sd/2;
	// system("pause");
	//cout <<setprecision(20)<< d_sd << " " << d_kl << " " << d_ml << " " << d_re << " " << endl;

	gau[0].SetOutputPrecision(6);

	if (outputfile == 1)
	{
		out << "Delta_SD: " << log(d_sd) / log(2) << "  Delta_KL: " << log(abs(d_kl)) / log(2) << " " << "  Delta_RE: " << log(d_re) / log(2) << "  Delta_RE's x: " << ren - rc / 2 << "," << (ren2 - rc / 2) << "  Delta_ML: " << log(d_ml) / log(2) << "  Delta_ML's x: " << mln - rc / 2 << "," << (mln2 - rc / 2) << endl;
		out << endl;
	}
	cout  << "Delta_SD: " << log(d_sd) / log(2) << "  Delta_KL: " << log(abs(d_kl)) / log(2)  << " " << "  Delta_RE: " << log(d_re) / log(2) << "  Delta_RE's x: " << ren - rc / 2 << "," <<  (ren2 - rc / 2) << "  Delta_ML: " << log(d_ml) / log(2) << "  Delta_ML's x: " << mln - rc / 2 << "," << (mln2 - rc / 2) << endl;
	cout << endl;
}

 

RR* eres_dx(RR* d1,RR* d2,int a,int b,RR s1,RR s2,int num,int* ra) {
	RR s_new = sqrt(a*a*s1*s1 + b*b*s2*s2);
	int range ;
	conv(range, ceil(t*s_new));

	

	RR* gau;
	
	 
	gau = new RR[2 * range + 1];
	memset(gau, 0, (2 * range + 1) * sizeof(RR));
	 
	
	
	RR trans;	 
	if (precision > 0) {
		trans.SetPrecision(precision);		 
	}
	
	RR ty = to_RR(0);
	 
	for (int i = 0; i <= num; i++) 

		for (int j = 0; j <= num; j++) {
			int tem = a*i + b*j - (num*(a + b) / 2)+range;
			if (tem >= 0 && tem <= 2 * range) {
				 
				gau[tem] += d1[i] * d2[j];
	 
			}
		}
	 
	for (int i = 0; i <= 2 * range; i++) {
		trans = to_RR(gau[i]);
		gau[i] = trans;
	}
	
	ra[0] = 2 * range;
	
	return gau;
}



int main() {

	init(6, 6, 200);
	RR* g1 = NULL;
	RR* g2 = NULL;
	RR* g3 = NULL;
	int r1, r2,r3;
	double s0;
	RR s1, s2,tmp;
	long oldpr;

	RR* gc;
	int rc;
	RR* gi1, gi2;
	int ri1, ri2;
	int a, b;
	RR s_new;
	int max_pre = 200;


	// Experiment 1
	if (EX1 == 1) {
		cout << "Experiment 1: " << endl;

		s0 = 34;
		s1 = to_RR(s0);
		s2 = to_RR(s0);

		a = 1;
		b = 1;
		set_pre(max_pre);
		s_new = sqrt(a*a*s1*s1 + b*b*s2*s2);
		gc = init_Gaussian(s_new, &rc);
		gi1 = init_Gaussian(s1, &ri1);
		//gi2 = init_Gaussian(s2, &ri2);



		for (int i = 53; i <=173; i ++) {

			set_pre(i);

			g1 = init_Gaussian(s1, &r1);
			g2 = init_Gaussian(s2, &r2);
			g3 = eres_dx(g1, g2, a, b, s1, s2, r1, &r3);

			
			if (outputfile == 1) {
				out << "precision: " << tmp.precision() << endl;
				out << "input_errors:" << endl;
			}			 
			cout << "precision: " << tmp.precision() << endl;
			cout << "input_errors:" << endl;		 

			compare(gi1, g1, ri1);

			cout << "output_errors:" << endl;
			if (outputfile == 1) {
				out << "output_errors:" << endl;
			}


			compare(gc, g3, rc);

			delete[] g1;
			delete[] g2;
			delete[] g3;

			cout << endl;
			if (outputfile == 1) {
				out  << endl;
			}
		}

		delete[] gc;
		delete[] gi1;
	}

	if (EX2 == 1) {
		// Experiment 2  a=11,b=1,s0=19.53

		cout << "Experiment 2: " << endl;

		s0 = 19.53;
		s1 = to_RR(s0);
		s2 = to_RR(s0);

		a = 11;
		b = 1;
		set_pre(max_pre);
		s_new = sqrt(a*a*s1*s1 + b*b*s2*s2);
		gc = init_Gaussian(s_new, &rc);
		gi1 = init_Gaussian(s1, &ri1);

		for (int i = 53; i <=173; i ++) {

			set_pre(i);

			g1 = init_Gaussian(s1, &r1);
			g2 = init_Gaussian(s2, &r2);
			g3 = eres_dx(g1, g2, a, b, s1, s2, r1, &r3);


			if (outputfile == 1) {
				out << "precision: " << tmp.precision() << endl;
				out << "input_errors:" << endl;
			}
			cout << "precision: " << tmp.precision() << endl;
			cout << "input_errors:" << endl;

			compare(gi1, g1, ri1);

			cout << "output_errors:" << endl;
			if (outputfile == 1) {
				out << "output_errors:" << endl;
			}


			compare(gc, g3, rc);

			delete[] g1;
			delete[] g2;
			delete[] g3;

			cout << endl;
			if (outputfile == 1) {
				out << endl;
			}
		}

		delete[] gc;
		delete[] gi1;
	}

	if (EX3 == 1) {

		// Experiment 3  a=1 to 11,b=1,s0=19.53

		cout << "Experiment 3: " << endl;

		s0 = 19.53;
		s1 = to_RR(s0);
		s2 = to_RR(s0);

		b = 1;

		for (int j = 1; j <= 11; j++) {
			a = j;
			set_pre(max_pre);
			s_new = sqrt(a*a*s1*s1 + b*b*s2*s2);
			gc = init_Gaussian(s_new, &rc);
			gi1 = init_Gaussian(s1, &ri1);

			cout << "a=" <<a<< endl;
			if (outputfile == 1)
				out << "a=" << a << endl;

			for (int i = 53; i <= 173; i++) {

				set_pre(i);

				g1 = init_Gaussian(s1, &r1);
				g2 = init_Gaussian(s2, &r2);
				g3 = eres_dx(g1, g2, a, b, s1, s2, r1, &r3);


				if (outputfile == 1) {
					out << "precision: " << tmp.precision() << endl;
					out << "input_errors:" << endl;
				}
				cout << "precision: " << tmp.precision() << endl;
				cout << "input_errors:" << endl;

				compare(gi1, g1, ri1);

				cout << "output_errors:" << endl;
				if (outputfile == 1) {
					out << "output_errors:" << endl;
				}

				compare(gc, g3, rc);


				delete[] g1;
				delete[] g2;
				delete[] g3;

				cout << endl;
				if (outputfile == 1) {
					out << endl;
				}
			}

			delete[] gc;
			delete[] gi1;
		}
	}

	if (EX4 == 1) {
		// Experiment 4  a=7,b=1,s0=30.8185

		cout << "Experiment 4_1: " << endl;
		double ts =  215.73;


		s0 = ts/7.0;
		s1 = to_RR(s0);
		s2 = to_RR(s0);

		a = 7;
		b = 1;
		set_pre(max_pre);
		s_new = sqrt(a*a*s1*s1 + b*b*s2*s2);
		gc = init_Gaussian(s_new, &rc);
		gi1 = init_Gaussian(s1, &ri1);

		for (int i = 53; i <= 200 ; i++) {

			set_pre(i);

			g1 = init_Gaussian(s1, &r1);
			g2 = init_Gaussian(s2, &r2);
			g3 = eres_dx(g1, g2, a, b, s1, s2, r1, &r3);


			if (outputfile == 1) {
				out << "precision: " << tmp.precision() << endl;
				out << "input_errors:" << endl;
			}
			cout << "precision: " << tmp.precision() << endl;
			cout << "input_errors:" << endl;

			compare(gi1, g1, ri1);

			cout << "output_errors:" << endl;
			if (outputfile == 1) {
				out << "output_errors:" << endl;
			}


			compare(gc, g3, rc);

			delete[] g1;
			delete[] g2;
			delete[] g3;

			cout << endl;
			if (outputfile == 1) {
				out << endl;
			}
		}

		delete[] gc;
		delete[] gi1;


		/*cout << "Experiment 4_2: " << endl;

		s0 = ts / 6.0;
		s1 = to_RR(s0);
		s2 = to_RR(s0);

		a = 6;
		b = 1;
		set_pre(max_pre);
		s_new = sqrt(a*a*s1*s1 + b*b*s2*s2);
		gc = init_Gaussian(s_new, &rc);
		gi1 = init_Gaussian(s1, &ri1);

		for (int i = 53; i <= 200; i++) {

			set_pre(i);

			g1 = init_Gaussian(s1, &r1);
			g2 = init_Gaussian(s2, &r2);
			g3 = eres_dx(g1, g2, a, b, s1, s2, r1, &r3);


			if (outputfile == 1) {
				out << "precision: " << tmp.precision() << endl;
				out << "input_errors:" << endl;
			}
			cout << "precision: " << tmp.precision() << endl;
			cout << "input_errors:" << endl;

			compare(gi1, g1, ri1);

			cout << "output_errors:" << endl;
			if (outputfile == 1) {
				out << "output_errors:" << endl;
			}


			compare(gc, g3, rc);

			delete[] g1;
			delete[] g2;
			delete[] g3;

			cout << endl;
			if (outputfile == 1) {
				out << endl;
			}
		}

		delete[] gc;
		delete[] gi1;*/
	 


	}

	// Experiment 5 a=s_(i-1)/sqrt(2)/\eta b=max(a-1,1);
	if (EX5 == 1) {
		cout << "Experiment 5: " << endl;

	 

		s0 = 34;
		s1 = to_RR(s0);
		s2 = to_RR(s0);

		a = 4;
		b = 3;
		set_pre(max_pre);
		s_new = sqrt(a*a*s1*s1 + b*b*s2*s2);
		gc = init_Gaussian(s_new, &rc);
		gi1 = init_Gaussian(s1, &ri1);
		//gi2 = init_Gaussian(s2, &ri2);



		for (int i = 53; i <= 173; i++) {

			set_pre(i);

			g1 = init_Gaussian(s1, &r1);
			g2 = init_Gaussian(s2, &r2);
			g3 = eres_dx(g1, g2, a, b, s1, s2, r1, &r3);


			if (outputfile == 1) {
				out << "precision: " << tmp.precision() << endl;
				out << "input_errors:" << endl;
			}
			cout << "precision: " << tmp.precision() << endl;
			cout << "input_errors:" << endl;

			compare(gi1, g1, ri1);

			cout << "output_errors:" << endl;
			if (outputfile == 1) {
				out << "output_errors:" << endl;
			}


			compare(gc, g3, rc);

			delete[] g1;
			delete[] g2;
			delete[] g3;

			cout << endl;
			if (outputfile == 1) {
				out << endl;
			}
		}

		delete[] gc;
		delete[] gi1;
	}
	/*system("pause");
	for (int i = 0; i <= r1; i++) {

		cout<<setprecision(20) << i << " " << ((double *)g1)[i] << endl;

	}*/

	system("pause");
	return 0;



}