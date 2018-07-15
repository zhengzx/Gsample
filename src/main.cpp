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

RR PI ;
RR EE = to_RR(2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274274663919320030599218174135966290435729003342952605956307381323286279434907632338298807531952510190115738341879307021540891499348841
);
int outputfile = 1;
 
int precision;

int max_precision;

ofstream out("result.txt");

int Renyi_a =  2;

int EX1 = 1; //a=11,b=1,s1=s2=19.53sqrt(2pi),t=3-8, pre=53-200 (\Delta_{SD} and \Delta_{KL})

int EX2 = 0; //a=4,b=3,s1=s2=34,t=3-8,pre=53-200 (\Delta_{RE} and \Delta_{ML})

int EX3 = 0; //a=11,b=1,s1=s2=19.53sqrt(2pi),t=5.35 ,pre=53-200 ([PDG14]: pre=72, Modified [PDG14]: pre=130)

int EX4 = 0; //a=4,b=3,s1=s2=34,t=6,pre=53-200 ([MW17]: pre=60, \Delta_{ML}\le 2^{-55}  Modified [MW17]: pre=113, \Delta_{KL}\le 2^{-110})

int EX5 = 0; //a=11,b=1,s1=s2=19.53sqrt(2pi),t=3,5,7,  Renyi_a=2,200 ,pre=53-200 (  \Delta_{RD} )


RR Gfun(RR x, RR c, RR s) {
	return exp(-1 * PI* pow(abs((x - c) / s), to_RR(2)));
}

RR etatoep(RR t) {

	return 2 / (pow(EE, t*t*PI) - 2);

}


RR* init_Gaussian(RR s, int * ra, RR t) {

	int range;
	conv(range, ceil(s * t));
	//cout << range << endl;
	RR* gau = new RR[2 * range + 1];
	memset(gau, 0, (2 * range + 1) * sizeof(RR));
	RR trans;
	
	trans.SetPrecision(max_precision);
	PI = ComputePi_RR();

	RR sum = to_RR(0);
	for (int i = 0; i <= 2 * range; i++) {

		RR tmp = Gfun(to_RR(i), to_RR(range), s);
		sum += tmp;
		gau[i] = tmp;

	}

	for (int i = 0; i <= 2 * range; i++) {
		gau[i] = gau[i] / sum;		
	}

	//convert precision
	trans.SetPrecision(precision);
	sum = to_RR(0);
	for (int i = 0; i <= 2 * range; i++) {
		
		trans = to_RR(gau[i]);		
		gau[i] = trans;	
		sum += gau[i];
	}
	for (int i = 0; i <= 2 * range; i++) {
		gau[i] = gau[i] / sum;
	}
	//cout << sum << endl;
	//system("pause");

	ra[0] = 2 * range;
	return gau;
}

RR* eres_dx(RR* d1, RR* d2, int a, int b, RR s1, RR s2, int num, int* ra, RR t) {
	RR s_new = sqrt(a*a*s1*s1 + b*b*s2*s2);
	int range;
	conv(range, ceil(t * s_new));

	RR* gau;

	gau = new RR[2 * range + 1];
	memset(gau, 0, (2 * range + 1) * sizeof(RR));

	RR trans;
	trans.SetPrecision(max_precision);

	RR ty = to_RR(0);

	for (int i = 0; i <= num; i++)

		for (int j = 0; j <= num; j++) {
			int tem = a*i + b*j - (num*(a + b) / 2) + range;
			if (tem >= 0 && tem <= 2 * range) {

				gau[tem] += d1[i] * d2[j];

			}
		}
	
	
	//convert precision
	//trans.SetPrecision(precision);

	//for (int i = 0; i <= 2 * range; i++) {

	//	//cout << NumBits(gau[i].x) << " " << gau[i] << endl;
	//	trans = to_RR(gau[i]);
	//	gau[i] = trans;
	//	/*cout << NumBits(gau[i].x) << " " << gau[i] << " " << log(abs(gau[i]))/log(2) << endl;
	//	system("pause");*/
	//}

	ra[0] = 2 * range;

	return gau;
}

void compare(RR* gau, RR* gc, int rc,int r1) {
	//compare

	RR trans;
	 
	trans.SetPrecision(max_precision);

	RR d_sd = to_RR(0), d_kl = to_RR(0), d_kll = to_RR(0), d_ml = to_RR(0), d_re = to_RR(0), d_ry = to_RR(0);

	int ren = -1;
	int mln = -1;

	int ren2 = -1;
	int mln2 = -1;


	int cnt = 0;
	int ci = 0;

	RR tt1 = to_RR(0);
	RR tt2 = to_RR(0);
	for (int i = 0; i <= rc; i++) {
		ci = r1 / 2 - rc / 2 + i;
		//cout << r1 << " " << rc << " "<<ci<<endl;

		if (gau[i] > 0) {
			d_sd += abs(gc[ci] - gau[i]);			
			
			if (abs(gc[ci] - gau[i]) / gau[i]> d_re  && i <= rc / 2) {
				//cout << log(abs(gc[i] - gau[i])) / log(2) << " " << log(gc[i]) / log(2) << " " << gc[i].precision() << " " << i << " " << rc <<" "<< log(abs(abs(gc[i] - gau[i]) / gc[i] - d_re)) / log(2) << endl;
				d_re = abs(gc[ci] - gau[i]) / gau[i];
				ren = i;
				ren2 = rc - i;
			}

			if (gc[ci] > 0) {
				d_ry += power(gau[i],Renyi_a) / power(gc[ci],Renyi_a-1);
				tt1 += gau[i];
				tt2 += gc[ci];
				d_kll += log(gc[ci] / gau[i])*gc[ci];
				d_kl += log(gau[i] / gc[ci])*gau[i];
				if (abs(log(gc[ci] / gau[i])) > d_ml  && i <= rc / 2) {
					d_ml = abs(log(gc[ci] / gau[i]));
					mln = i;
					mln2 = rc - i;
				}
			}
		}

	}
	
	d_sd = d_sd / 2.0;
	gau[0].SetOutputPrecision(6);

	if (outputfile == 1)
	{
		out << "Delta_SD: " << log(d_sd) / log(2) << "  Delta_KL: " << log(abs(d_kl)) / log(2)    << "  Delta_RE: " << log(d_re) / log(2) << "  Delta_RE's x: " << ren - rc / 2 << "," << (ren2 - rc / 2) << "  Delta_ML: " << log(d_ml) / log(2) << "  Delta_ML's x: " << mln - rc / 2 << "," << (mln2 - rc / 2) << "  Delta_Renyi: "<< log(abs(to_RR(1)-pow(d_ry,to_RR(1)/(to_RR(Renyi_a-1)))))/log(to_RR(2)) << endl;
		out << endl;
	}
	cout << "Delta_SD: " << log(d_sd) / log(2) << "  Delta_KL: " << log(abs(d_kl)) / log(2) << "  Delta_KLL: " << log(abs(d_kll)) / log(2)   << "  Delta_RE: " << log(d_re) / log(2) << "  Delta_RE's x: " << ren - rc / 2 << "," << (ren2 - rc / 2) << "  Delta_ML: " << log(d_ml) / log(2) << "  Delta_ML's x: " << mln - rc / 2 << "," << (mln2 - rc / 2) << "  Delta_Renyi: " << log(abs(to_RR(1) - pow(d_ry, to_RR(1) / (to_RR(Renyi_a - 1))))) / log(to_RR(2))   << endl;

	
	cout << endl;
	// cout << log(gau[rc/2]) / log(2) << endl;
}


int main() {


	max_precision = 500;
	RR pre;
	pre.SetPrecision(max_precision);



	RR* g1 = NULL;
	RR* g2 = NULL;
	RR* g3 = NULL;
	int r1, r2, r3;
	double s0;
	RR s1, s2, tmp;
	long oldpr;

	RR* gc;
	int rc;
	RR* gi1, gi2;
	int ri1, ri2;
	int a, b;
	RR s_new;
	RR t;

	RR* gc2;
	int rc2;
	if (EX1 == 1) {
		// Experiment 1  a=11,b=1,s1=s2=19.53sqrt(2pi),t=3-8 ,pre=53-200 (\Delta_{SD} and \Delta_{KL})
		 

		
		cout << "Experiment 1: " << endl;
		if (outputfile == 1) {
			out << "Experiment 1: " << endl;
			 
		}


		pre.SetPrecision(max_precision);

		PI = ComputePi_RR();
		double po;
		conv(po, PI);

	

		a = 11;
		b = 1;

		//s_new = sqrt(a*a*s1*s1 + b*b*s2*s2);
		//s_new = sqrt(to_RR(11)*to_RR(11)*to_RR(19.53)*to_RR(19.53)*  2 * PI  + to_RR(19.53)*to_RR(19.53) * 2 * PI );
		s1 = sqrt(  PI/log(to_RR(2)) / to_RR(a*a+b*b))*to_RR(254);
		s2 = sqrt(PI / log(to_RR(2)) / to_RR(a*a + b*b))*to_RR(254);
		s_new = sqrt(PI / log(to_RR(2))  )*to_RR(254);

		for (int j = 300; j <= 800; j += 10) {

			t = to_RR(j/100.0);

			precision = max_precision;
			gc = init_Gaussian(s_new, &rc, to_RR(10));
			cout << t<<" epsilon_t: " << log( etatoep(t)) / log(2)<<endl;

			 

			if (outputfile == 1) {
				out << "t: " << t << " epsilon_t: " << log(etatoep(t)) / log(2) << endl;
				//cout << "input_errors:" << endl;
			}
			 
			//precision = max_precision;
			
			//gi1 = init_Gaussian(s1, &ri1, t);

			for (int i = 53; i <= 200; i+=1 ) {


				precision = i;

				
				g1 = init_Gaussian(s1, &r1, t);
				g2 = init_Gaussian(s2, &r2, t);

				 
				//cout << log(g1[r1 / 2]) / log(2) << endl;
				g3 = eres_dx(g1, g2, a, b, s1, s2, r1, &r3, t);


				if (outputfile == 1) {
					out << "precision: " << precision << endl;
					//out << "input_errors:" << endl;
				}
				cout << "precision: " << precision << endl;
				//cout << "input_errors:" << endl;

				//compare(gi1, g1, ri1);

				cout << "output_errors:" << endl;
				if (outputfile == 1) {
					out << "output_errors:" << endl;
				}

			   
				compare(g3, gc, r3, rc);
				 
				

				delete[] g1;
				delete[] g2;
				delete[] g3;
				

				cout << endl;
				if (outputfile == 1) {
					out << endl;
				}
			}
			delete[] gc;
			 

		}
		

	}


	if (EX2 == 1) {
		// Experiment 2  a=4,b=3,s1=s2=34,t=3-8,pre=53-200 (\Delta_{RE} and \Delta_{ML})
		 

		cout << "Experiment 2: " << endl;
		if (outputfile == 1) {
			out << "Experiment 2: " << endl;

		}
		pre.SetPrecision(max_precision);

		PI = ComputePi_RR();
		double po;
		conv(po, PI);

		s0 = 34;
		s1 = to_RR(s0);
		s2 = to_RR(s0);


		a = 4;
		b = 3;

		s_new = sqrt(to_RR(a)*to_RR(a)*s1*s1 + to_RR(b)*to_RR(b)*s2*s2);


		for (int j = 300; j <= 800; j += 10) {

			t = to_RR(j / 100.0);

			precision = max_precision;
			gc = init_Gaussian(s_new, &rc, to_RR(10));
			cout << t << " epsilon_t: " << log(etatoep(t)) / log(2) << endl;
			if (outputfile == 1) {
				out << "t: " << t << " epsilon_t: " << log(etatoep(t)) / log(2) << endl;
				//out << "input_errors:" << endl;
			}

			//precision = max_precision;

			//gi1 = init_Gaussian(s1, &ri1, t);

			for (int i = 53; i <= 200; i += 1) {


				precision = i;


				g1 = init_Gaussian(s1, &r1, t);
				g2 = init_Gaussian(s2, &r2, t);

				//cout << log(g1[r1 / 2]) / log(2) << endl;
				g3 = eres_dx(g1, g2, a, b, s1, s2, r1, &r3, t);


				if (outputfile == 1) {
					out << "precision: " << precision << endl;
					//out << "input_errors:" << endl;
				}
				cout << "precision: " << precision << endl;
				//cout << "input_errors:" << endl;

				//compare(gi1, g1, ri1);

				cout << "output_errors:" << endl;
				if (outputfile == 1) {
					out << "output_errors:" << endl;
				}


				compare(g3, gc, r3, rc);

				delete[] g1;
				delete[] g2;
				delete[] g3;


				cout << endl;
				if (outputfile == 1) {
					out << endl;
				}
			}
			delete[] gc;


		}


	}

	if (EX3 == 1) {
		// Experiment 3  a=11,b=1,s1=s2=19.53sqrt(2pi),t=5.35 ,pre=53-200 ([PDG14]: pre=72, Modified [PDG14]: pre=130)

 

		cout << "Experiment 3: " << endl;
		if (outputfile == 1) {
			out << "Experiment 3: " << endl;

		}
		pre.SetPrecision(max_precision);

		PI = ComputePi_RR();
		double po;
		conv(po, PI);

	 


		a = 11;
		b = 1;

		s1 = sqrt(PI / log(to_RR(2)) / to_RR(a*a + b*b))*to_RR(254);
		s2 = sqrt(PI / log(to_RR(2)) / to_RR(a*a + b*b))*to_RR(254);
		s_new = sqrt(PI / log(to_RR(2)))*to_RR(254);


		for (int j = 535; j <= 535; j += 1) {

			t = to_RR(j / 100.0);

			precision = max_precision;
			gc = init_Gaussian(s_new, &rc, to_RR(10));
			cout << t << " epsilon_t: " << log(etatoep(t)) / log(2) << endl;
			if (outputfile == 1) {
				out << "t: " << t << " epsilon_t: " << log(etatoep(t)) / log(2) << endl;
				//out << "input_errors:" << endl;
			}

			//precision = max_precision;

			//gi1 = init_Gaussian(s1, &ri1, t);

			for (int i = 53; i <= 200; i += 1) {


				precision = i;


				g1 = init_Gaussian(s1, &r1, t);
				g2 = init_Gaussian(s2, &r2, t);

				//cout << log(g1[r1 / 2]) / log(2) << endl;
				g3 = eres_dx(g1, g2, a, b, s1, s2, r1, &r3, t);


				if (outputfile == 1) {
					out << "precision: " << precision << endl;
					//out << "input_errors:" << endl;
				}
				cout << "precision: " << precision << endl;
				//cout << "input_errors:" << endl;

				//compare(gi1, g1, ri1);

				cout << "output_errors:" << endl;
				if (outputfile == 1) {
					out << "output_errors:" << endl;
				}


				compare(g3, gc, r3, rc);

				delete[] g1;
				delete[] g2;
				delete[] g3;


				cout << endl;
				if (outputfile == 1) {
					out << endl;
				}
			}
			delete[] gc;


		}


	}

	if (EX4 == 1) {
		// Experiment 4  a=4,b=3,s1=s2=34,t=6,pre=53-200 ([MW17]: pre=60, \Delta_{ML}\le 2^{-55}  Modified [MW17]: pre=113, \Delta_{KL}\le 2^{-110})

	 
		cout << "Experiment 4: " << endl;
		if (outputfile == 1) {
			out << "Experiment 4: " << endl;

		}
		pre.SetPrecision(max_precision);

		PI = ComputePi_RR();
		double po;
		conv(po, PI);

		s0 = 34;
		s1 = to_RR(s0);
		s2 = to_RR(s0);


		a = 4;
		b = 3;

		s_new = sqrt(to_RR(a)*to_RR(a)*s1*s1 + to_RR(b)*to_RR(b)*s2*s2);


		for (int j = 600; j <= 600; j += 10) {

			t = to_RR(j / 100.0);

			precision = max_precision;
			gc = init_Gaussian(s_new, &rc, to_RR(10));
			cout << t << " epsilon_t: " << log(etatoep(t)) / log(2) << endl;
			if (outputfile == 1) {
				out << "t: " << t << " epsilon_t: " << log(etatoep(t)) / log(2) << endl;
				//out << "input_errors:" << endl;
			}

			//precision = max_precision;

			//gi1 = init_Gaussian(s1, &ri1, t);

			for (int i = 53; i <= 200; i += 1) {


				precision = i;


				g1 = init_Gaussian(s1, &r1, t);
				g2 = init_Gaussian(s2, &r2, t);

				//cout << log(g1[r1 / 2]) / log(2) << endl;
				g3 = eres_dx(g1, g2, a, b, s1, s2, r1, &r3, t);


				if (outputfile == 1) {
					out << "precision: " << precision << endl;
					//out << "input_errors:" << endl;
				}
				cout << "precision: " << precision << endl;
				//cout << "input_errors:" << endl;

				//compare(gi1, g1, ri1);

				cout << "output_errors:" << endl;
				if (outputfile == 1) {
					out << "output_errors:" << endl;
				}


				compare(g3, gc, r3, rc);

				delete[] g1;
				delete[] g2;
				delete[] g3;


				cout << endl;
				if (outputfile == 1) {
					out << endl;
				}
			}
			delete[] gc;


		}


	}

	if (EX5 == 1) {
		// Experiment 5   a=11,b=1,s1=s2=19.53sqrt(2pi),t=3,5,7,  Renyi_a=2,200 ,pre=53-200 (  \Delta_{RD} )



		cout << "Experiment 5: " << endl;
		if (outputfile == 1) {
			out << "Experiment 5: " << endl;

		}


		pre.SetPrecision(max_precision);

		PI = ComputePi_RR();
		double po;
		conv(po, PI);



		a = 11;
		b = 1;

		//s_new = sqrt(a*a*s1*s1 + b*b*s2*s2);
		//s_new = sqrt(to_RR(11)*to_RR(11)*to_RR(19.53)*to_RR(19.53)*  2 * PI  + to_RR(19.53)*to_RR(19.53) * 2 * PI );
		s1 = sqrt(PI / log(to_RR(2)) / to_RR(a*a + b*b))*to_RR(254);
		s2 = sqrt(PI / log(to_RR(2)) / to_RR(a*a + b*b))*to_RR(254);
		s_new = sqrt(PI / log(to_RR(2)))*to_RR(254);

		for (int j = 300; j <= 800; j += 200) {

			t = to_RR(j / 100.0);

			precision = max_precision;
			gc = init_Gaussian(s_new, &rc, to_RR(10));
			cout << t << " epsilon_t: " << log(etatoep(t)) / log(2) << endl;



			if (outputfile == 1) {
				out << "t: " << t << " epsilon_t: " << log(etatoep(t)) / log(2) << endl;
				//cout << "input_errors:" << endl;
			}

			//precision = max_precision;

			//gi1 = init_Gaussian(s1, &ri1, t);

			for (int i = 53; i <= 200; i += 1) {


				precision = i;


				g1 = init_Gaussian(s1, &r1, t);
				g2 = init_Gaussian(s2, &r2, t);


				//cout << log(g1[r1 / 2]) / log(2) << endl;
				g3 = eres_dx(g1, g2, a, b, s1, s2, r1, &r3, t);


				if (outputfile == 1) {
					out << "precision: " << precision << endl;
					//out << "input_errors:" << endl;
				}
				cout << "precision: " << precision << endl;
				//cout << "input_errors:" << endl;

				//compare(gi1, g1, ri1);

				cout << "output_errors:" << endl;
				if (outputfile == 1) {
					out << "output_errors:" << endl;
				}

				for (int k = 2; k <= 200; k += 198) {

					Renyi_a = k;
					cout << "Renyi_a:" << " " << k << endl;
					compare(g3, gc, r3, rc);
				}


				delete[] g1;
				delete[] g2;
				delete[] g3;


				cout << endl;
				if (outputfile == 1) {
					out << endl;
				}
			}
			delete[] gc;


		}


	}

	system("pause");
	return 0;
}