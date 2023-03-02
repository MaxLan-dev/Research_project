#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <vector>
#include <stdio.h>
using namespace std;
double n_e(double I_6716, double I_6731) {
    double x = I_6716/I_6731;
    return pow(10, 8.448 - 15.101 * x + 14.419 * pow(x, 2) - 5.115 * pow(x, 3));
}
//double Te_method( double R, double R_2, double R_3, double x,);

double Te_method( double R, double R_2, double R_3, double n_e, double *t2, double *t3 ) {
    
    double t_3;
    double t_3new = 1;
    double C_T;
    //double y;
    double V; 
    double x_3; 
    double O_H;
    int count = 0;

    do {
        t_3 = t_3new;
        x_3 = 1.e-4 * n_e * pow(t_3, -0.5);
        V = (1. + 0.0004*x_3)/(1.+0.044*x_3);
        C_T = (8.44 - 1.09 * t_3 + 0.5 * pow(t_3, 2) - 0.08 * pow(t_3, 3)) * V;
        //y = R_3 / R;
        t_3new = 1.432 / (log10(R_3 / R) - log10(C_T));
	count++; if (count>=50) {t_3=t_3new=0;break;}
        //cout << "t_3 : " << t_3 << " t_3new = " << t_3new <<  endl;
    }
    while (fabs(t_3new-t_3)/t_3 > 1e-8 );
    
    double t_2; 
    double x_2;
    t_2 = 0.7 * t_3 + 0.3;
    x_2 = pow(10, -4) * n_e * pow(t_2, -0.5);
    double logO1_H1, logO2_H1;
    logO2_H1 = (log10(R_3) + 6.200 + (1.251/t_3new) - (0.55 * log10(t_3new)) - (0.014 * t_3new)) - 12 ;
    logO1_H1 = (log10(R_2) + 5.961 + (1.676/t_2) - (0.40 * log10(t_2)) - (0.034 * t_2) + log10(1+1.35*x_2)) - 12;
    double O1_h1, O2_h1;
    O1_h1 = pow(10, logO1_H1);
    O2_h1 = pow(10, logO2_H1);
    O_H = O1_h1 + O2_h1;
   // cout << scientific;
 //   cout << "O+/H+ : " << O1_h1 << '\n' << "O++/H+ : " << O2_h1 << '\n' << "O/H : " << O_H <<'\n'  ;
    *t2 = t_2;
    *t3 = t_3;
    return O_H ;
    
}
int main( int argc, char *argv[])
{
    double R_2;
    double t_2, t_3;
    double R_3;
    double R, I_6716, I_6731;
    const int FLMN = 6;
    double I_4959, I_5007;
    double O_H_grd, O_H;
        
    
    if (argc != 4 ) {
    cout << "./Te_method <input file> <output file>"<< endl ;
    cout << "Please enter I [OII] 3727+3729 / H b  : ";
    cin >> R_2;
    cout << "Please enter I [OIII] 4959+5007 / H b : ";
    cin >> R_3;
    cout << "Please enter I [OII] 4363 / H b : ";
    cin >> R;
    cout << "Please enter I 6716 : ";
    cin >> I_6716;
    cout << "Please enter I6731 : ";
    cin >> I_6731;
    cout << "O/H =  " << scientific << std :: setprecision(2) << Te_method(R, R_2, R_3, n_e(I_6716, I_6731) ,&t_2, &t_3) << " Te(O+) = " << t_2 * 1.e4 << "K, Te(O++) = " << t_3 * 1.e4 << "K\n" ; 
    }
    else {
        ///cout << "Input file name : " << argv[1] << endl;
        //ifstream in(argv[1]);
        FILE * pFile;
        char mystring [200],Model[20];
	pFile = fopen (argv[1] , "r");
        ofstream out(argv[2], std :: ios :: app);
        ofstream bad(argv[3], std :: ios :: app);
	long int counter = 0;
	double Metals = -0.3;
	double Dev=0.;
	while(!feof(pFile))
	 { 
	    
	    if(Metals>=0.3) Metals = -0.3; else Metals+= 0.3;
		fgets (mystring , 12 , pFile);
		printf( "%s ! counter=%li\n" ,mystring, (long int)(counter+1) );
	        fscanf(pFile,"\n");

    	       for(int j=0;j<FLMN;j++)
		{
    		double dummy;
	    	fgets (mystring , 16 , pFile);
	    	printf("%s ",mystring);
	        fscanf (pFile, "%le", &dummy);
		printf("%le\n",dummy);
		fscanf(pFile, "\n");
		if(!j) R_2=dummy;
		else if(j==1) R = dummy;
		else if(j==2) I_4959 = dummy;
		else if(j==3) I_5007 = dummy;
		else if(j==4) I_6716 = dummy;
		else if(j==5) {I_6731 = dummy; fscanf(pFile, "\n"); if(feof(pFile)) break;}
		}

		R_3 = I_4959 + I_5007; 
		O_H = Te_method(R, R_2, R_3, n_e(I_6716, I_6731) ,&t_2, &t_3);
		O_H_grd = pow(10.,(Metals-3.481146));
		O_H_grd = 12.+log10(O_H_grd);O_H = 12.+log10(O_H);
		Dev = fabs(O_H_grd - O_H);
		// cout << scientific << std :: setprecision(2) << O_H << "\t" << pow(10.,(Metals-3.481146)) <<"\n" ;
		if(counter<10) sprintf(Model,"00000000%li\t",counter);else if(counter<100) sprintf(Model,"0000000%li\t",counter);else if(counter<1000) sprintf(Model,"000000%li\t",counter);
		if(Dev<=0.5) out<< Model <<"\t" << std :: setprecision(3) << O_H_grd <<"\t" << std :: setprecision(3) 
		   <<  O_H <<"\t" << std :: setprecision(2)  <<  Dev << " dex" << endl;
		else bad<< Model <<"\t" << std :: setprecision(3) << O_H_grd <<"\t" << std :: setprecision(3) 
		   <<  O_H <<"\t" << std :: setprecision(2)  <<  Dev << " dex"  << endl;
 	   counter++;
	 }
	   fclose (pFile);
	   out.close();
	   bad.close();
	}
	         
   return 0;}



