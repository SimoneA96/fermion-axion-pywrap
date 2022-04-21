#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <string>
using namespace std;

// Parameters
long double y[5], y_axi[5], diff[5], F[5];
long double** var;
long double* riso, *riso2, *beta;
long double dpsidr, ADM_mass;
long double k1[5], k2[5], k3[5], k4[5];
const long double Pi=3.14159265358979323846L, KK=100.0, m=1.0, polind=2;
long double fa_axion;
long double eps;
long double phi0c, P0c, omega, omega2, l_1, l_2;
long double delta_r, r, R, rtmp;
long double numbos, numfer, numbos2,numfer2, r99bos, r99fer, r99tot, rbos_iso, rfer_iso, r_isotropic_max, factor, rmax;
int nbos,nfer, nmax;
long double massbos, massfer, massT, massT2;
long double* conformal, *conformal2;
bool part, isotropic;
long double* bospotential;

long double epsilon = 1e-8;
int flag = 0;

long double min(long double z_1, long double z_2){
	if (z_1<z_2) return z_1;
	else return z_2;
}
long double V(long double phi) {
    return 2.0 * pow(m,2) * pow(fa_axion,2) / 0.22 * (1.0 - sqrt(1.0 - 4 * 0.22 * pow(sin( phi/(2.0*fa_axion)),2)));
}
long double dV(long double phi) {
     if(phi!=0)
         return 2.0 * pow(m,2)  * fa_axion * sin(phi/(2.0*fa_axion)) * cos(phi/(2.0*fa_axion)) /
	    (phi * sqrt(1.0 - 4.0 * 0.22 * pow(sin(phi/(2.0*fa_axion)),2) ));
     else
	 return 0.;
}

int func(long double r,  long double h) {
	long double a=y_axi[0];
	long double alp=y_axi[1];
	long double phi=y_axi[2];
	long double dphi=y_axi[3];
	long double P=y_axi[4];
	long double rho_tmp;

        eps = pow(KK,1/polind)*pow(P,(polind-1)/polind)/(polind-1);
        rho_tmp = pow(P/KK,1/polind);

        F[0] = 0.5*a*((1-a*a)/r + 8*Pi*r*( pow(dphi,2) + omega2 * pow(a*phi/alp, 2) + a*a * V(phi)+a*a*rho_tmp*(1+eps)));
        F[1] = 0.5*alp*((a*a-1)/r + 8*Pi*r*( pow(dphi,2) + omega2 * pow(a*phi/alp, 2) - a*a * V(phi)+a*a*P));
        F[2] = dphi;
        F[3] = -(1+a*a-8*Pi*pow(r*a,2)*(V(phi)+0.5*(rho_tmp*(1+eps)-P)))*dphi/r+(dV(phi)-omega2/pow(alp,2))* a*a * phi;
        F[4] =-(rho_tmp*(1+eps)+P)*0.5*((a*a-1)/r + 8*Pi*r*( pow(dphi,2) + omega2 * pow(a*phi/alp, 2) - a*a * V(phi)+a*a*P));

	return 0;
}

int Runge_Kutta_4(long double r,  long double omega2, long double h)  {
	long double half_h=0.5*h;
	for (int i=0; i<=4; i++) y_axi[i]=y[i];
	func(r, h);
	for (int i=0; i<=4; i++) k1[i]=h*F[i];
	for (int i=0; i<=4; i++) y_axi[i]=y[i]+0.5*k1[i];
//control for P not to become negative
        if (y_axi[4]<=0) {
		y_axi[4]=0;
		k1[4]=0;
		}
	func(r+half_h, h);
	for (int i=0; i<=4; i++) k2[i]=h*F[i];
	for (int i=0; i<=4; i++) y_axi[i]=y[i]+0.5*k2[i];
        if (y_axi[4]<=0) {
		y_axi[4]=0;
		k2[4]=0;
		}
	func(r+half_h, h);
	for (int i=0; i<=4; i++) k3[i]=h*F[i];
	for (int i=0; i<=4; i++) y_axi[i]=y[i]+k3[i];
        if (y_axi[4]<=0) {
		y_axi[4]=0;
		k3[4]=0;
		}
	func(r+h, h);
	for (int i=0; i<=4; i++) k4[i]=h*F[i];
	for (int i=0; i<=4; i++) diff[i]=(k1[i]+k4[i])/6 + (k2[i]+k3[i])/3;
	return 0;
}

int main(int argc, char *argv[]){

	
	if (argc<3){
	  cout<<"You should run './Solve data' "<<'\n';
	  cout<<"where data is the folder where you want to put your output"<<endl;
	  cout<<"and fa_exp is fa_axion = pow(10,fa_exp) for the axion potential"<<endl;
	  return 0;
	}

	int n, num, ntmp, input_num=1;
	ifstream  input;
	string folder;
	folder = argv[1];
	fa_axion = pow(10,atof(argv[2]));
	string tmp;
	string fa_str;
	string phi_str;
	string rho_str;
	input.open(folder + "/Output.txt", ios::in );
	input >> input_num >> num;

	//dynamical allocation
	var = new long double*[5];
	for (int i = 0;i<5;i++)
		var[i] = new long double[num];

	riso = new long double[num];
	riso2 = new long double[num];
	beta = new long double[num];
	conformal = new long double[num];
	conformal2 = new long double[num];
	bospotential = new long double[num];

	ofstream rfile;
	ofstream solfile;
	ofstream solfile1;
	ofstream solfile2;
	ofstream solfile3;
	ofstream solfile4;
	ofstream solfile5;
	ofstream alpsolfile;
	ofstream asolfile;
        ofstream evolutionfile;
	ofstream massfile;
	ofstream conffile;
	ofstream infofile;
	ofstream bospotential_file;
	ofstream numbosfer;
	ofstream existence_plot;

	fa_str = to_string(fa_axion);

	numbosfer.open(folder + "/numbosfer.txt");
        existence_plot.open(folder + "/existence_plot.txt");

	for(int j=1; j<=input_num; j++) {
		y[0]=1;
		y[1]=1;
		input >> phi0c >> P0c >> omega >> R;
		cout<< phi0c<<" "<<P0c<< " "<<omega<<endl;
		y[2]=phi0c;
		y[3]=0;
		y[4]=P0c;

		n=0;
		omega2=pow(omega,2);
		delta_r=R/(num);
		cout << setprecision(9);
		r=delta_r;
		part = 0;
		numbos = 0.0;
		numfer = 0.0;
		while (r<=R  +delta_r) {
			for (int i=0; i<=4; i++) {
				var[i][n]=y[i];
			}
	    		n++;
	    		Runge_Kutta_4(r, omega2, delta_r);
	    		for (int i=0; i<=4; i++) {
	    			y[i]+=diff[i];
	    		}
//control for P not to become negative
			if(y[4]<=0) 
				y[4]=0;

			if ((y[2] < 0) || diff[2] >0) {
				if(flag==0){
					cout<<r<<endl;
					flag=1;
					rmax = r;
					nmax = n;
				}
// I need to put a value for phi, because of the derivative of axion potential has phi at the denominator! 		

	        		y[2]=0;
				y[3]=0;
	    	    		
			}

		//Setting to 0 the pressure when it starts to do steps (for low values of P)
			if (y[4] == var[4][n-1] || var[4][n-1]==0) {

				y[4]=0;
	    	    		
			}


			r+=delta_r;
		}
//number of bosons and fermions
	        tmp = to_string(j);
		r=delta_r;
		n=0;
		numbos=0;
		numfer=0;
		while (r<=R) {
			numbos = numbos + 4.0*M_PI*var[0][n]*omega*pow(var[2][n]*r,2.0)/var[1][n]*delta_r;
			numfer = numfer + 4.0*M_PI*var[0][n]*pow(var[4][n]/KK,1/polind)*pow(r,2.0)*delta_r;
	    		n++;
			r+=delta_r;
		}

		r=delta_r;
		n=0;
		numbos2=0;
		numfer2=0;

//R99 or R95 for boson and fermions and total

		while (r<R) {
			numbos2 = numbos2 + 4.0*M_PI*var[0][n]*omega*pow(var[2][n]*r,2.0)/var[1][n]*delta_r;
			numfer2 = numfer2 + 4.0*M_PI*var[0][n]*pow(var[4][n]/KK,1/polind)*pow(r,2.0)*delta_r;
	    		n++;
			if(numbos2/numbos<0.99){
				r99bos=r;
				nbos = n;
				}
//			if(var[4][n] == 0 && var[4][n-1]!=0){
			if(numfer2/numfer<0.99){
				r99fer=r;
				nfer = n;
				}

			r+=delta_r;
		}

//Change of R99 to isotropic coordinates

		r_isotropic_max = pow( (1+sqrt(var[0][num-1]))/2,2)*R/var[0][num-1];


//Isotropic radius
	isotropic=0;
	if(isotropic){	
	r=delta_r;
	n=0;
	while(r<R){
		rtmp=r;
		factor=0;
		ntmp=n;
		while(rtmp<R){

			factor = factor + var[0][ntmp]*delta_r/rtmp;
			ntmp++;
			rtmp+=delta_r;		
		
		}
		riso[n]=r_isotropic_max*exp(-factor);
		r+=100*delta_r;
		n+=100;

	}

//using beta-factor improves the change to isotropic coordinates (see arXiv:gr-qc/9707045 Sec. IIIA)
	r=delta_r;
	n=0;
	while(r<R){
		rtmp=r;
		factor=0;
		ntmp=n;
		while(rtmp<R){

			factor = factor + (var[0][ntmp]-1)*delta_r/rtmp;
			ntmp++;
			rtmp+=delta_r;		
		
		}
		beta[n] =r_isotropic_max/R * exp(-factor);
		riso2[n]=r*beta[n];
		r+=100*delta_r;
		n+=100;

	}

	r=delta_r;
	n=0;

	while(r<R){
	
		conformal[n] = sqrt(r/riso[n]);
		conformal2[n] = sqrt(1.0/beta[n]);
		r+=100*delta_r;
		n+=100;

	}

	dpsidr = (conformal2[n-400]-8.0*conformal2[n-300]+ 8.0*conformal2[n-100]-conformal2[n])* 1.0 / (12.0*delta_r*100);
	ADM_mass = -2.0*riso2[n-200]*riso2[n-200]*dpsidr;
	cout<<"dpsidr= "<<dpsidr<<" ADM mass = "<<ADM_mass<<endl;
	}

//Boson and fermion mass as spatial volume integrals
// sqdetg = a*r^2
//rho_HD = (rho*(1+eps) + P) - P = rho + rho*eps = pow(P/K,1/polind) + P/(polind-1)
		r=delta_r;
		n=0;
		massbos=0;
		massfer=0;
		
		while (r<=R) {
			massfer = massfer + 4.0*M_PI*var[0][n]*pow(r,2)*( pow(var[4][n]/KK,1/polind) + var[4][n]/(polind-1) ) *delta_r;
			massbos = massbos + 4.0*M_PI*var[0][n]*pow(r,2.0)*(0.5*( omega2/var[1][n]/var[1][n]*var[2][n]*var[2][n] + 
				1/var[0][n]*var[3][n]*var[3][n] + V(var[2][n]) ))*delta_r;
			bospotential[n] =  4.0*M_PI*var[0][n]*pow(r,2.0)*(0.5*( omega2/var[1][n]/var[1][n]*var[2][n]*var[2][n] + 
				1/var[0][n]*var[3][n]*var[3][n] + V(var[2][n]) ))*delta_r;
	    		n++;

			r+=delta_r;
		}	
	

		massT=((R)/2)* (1-1/pow(var[0][num-1],2));
		massT2=0;
		r=delta_r;
		n=0;

		while (r<R) {
			massT2 = ((r)/2)* (1-1/pow(var[0][n],2));

			if(massT2>=0.99*massT) {
				r99tot=r;
				break;
				}
	    		n++;

			r+=delta_r;
		}


		//Normalization of the lapse
		r=delta_r;
		n=0;
		omega = omega/(var[0][num-1]*var[1][num-1]);
		while (r<R) {
			var[1][n] = var[1][n]/(var[0][num-1]*var[1][num-1]);
			n++;
			r+=delta_r;
		}

	        tmp = to_string(j);

		phi_str = to_string(phi0c);
		rho_str = to_string(pow(P0c/KK,1/polind));
		
	        rfile.open(folder + "/solr"+tmp+".txt");
	        solfile.open(folder + "/solphi"+tmp+".txt");
	        solfile1.open(folder + "/solP"+tmp+".txt");
	        solfile2.open(folder + "/solrho"+tmp+".txt");
	        solfile3.open(folder + "/solrhoiso"+tmp+".txt");
	        solfile4.open(folder + "/solphiiso"+tmp+".txt");
		solfile5.open(folder + "/solrhobosoniso"+tmp+".txt");
	        alpsolfile.open(folder + "/solalp"+tmp+".txt");
	        asolfile.open(folder + "/sola"+tmp+".txt");
	        massfile.open(folder + "/solmass"+tmp+".txt");
	        conffile.open(folder + "/conformal"+tmp+".txt");
	        infofile.open(folder + "/info"+tmp+".txt");
	        evolutionfile.open(folder + "/Mixed_"+phi_str +"_"+ rho_str +"_"+ fa_str+".dat");
		evolutionfile<<delta_r*10<<" "<<omega<<endl;

		n=0;
		r=delta_r;
		while (r<=R) {


			 if (n%10 == 0) {

				evolutionfile << r <<" "<<var[2][n]<<" "<<pow(var[1][n],2.)<<" "<<pow(var[0][n],2.)
					<<"\n"; //<<" "<<var[4][n]<<" "<<	pow(var[4][n]/KK,1/polind)<<"\n";	

			
			}


			 if (n%10 == 0) {
				rfile << r  << "\n";

				solfile << r<<" "<<var[2][n]  << "\n";

				solfile1 << r<<" "<<var[4][n]  << "\n";

				solfile2 << r <<" "<< pow(var[4][n]/KK,1/polind)  << "\n";
	
				alpsolfile << r <<" "<< var[1][n] << "\n";

				asolfile << (r)<<" "<< var[0][n] <<"\n";

				massfile << r <<" "<< ((r)/2)* \
				(1-1/pow(var[0][n],2)) << "\n";	

				solfile3 <<riso[n]<<" "<< pow(var[4][n]/KK,1/polind)  << "\n";

				solfile4 <<riso[n]<<" "<< var[2][n]  << "\n";

				solfile5 <<riso[n]<<" "<<0.5*( omega2/var[1][n]/var[1][n]*var[2][n]*var[2][n] + 
					1/var[0][n]/var[0][n]*var[3][n]*var[3][n] + V(var[2][n])) <<"\n";

				conffile <<riso[n]<<" "<< conformal[n]  << " "<<riso2[n]<<" "<<conformal2[n]<< "\n";
	
			
			}
	
			n++;
			r+=delta_r;
		}

		existence_plot<<pow(P0c/KK,1/polind)<<" "<<var[2][0]<<" "<<omega<<" "<<((R)/2)* (1-1/pow(var[0][num-1],2))<<endl;
		numbosfer<< var[2][0]<<" "<<pow(P0c/KK,1/polind)<<" "<<numbos<<" "<<numfer<<" "<<r99bos<<" "<<r99fer<<" "<<endl;
		numbosfer<< var[2][0]<<" "<<pow(P0c/KK,1/polind)<<" "<<numbos<<" "<<numfer<<" "<<r99bos<<" "<<r99fer<<" "<<endl;
		rfile << endl;
		solfile << endl;
		solfile1 << endl;
		solfile2 << endl;
		solfile3 << endl;
		solfile4 << endl;
		solfile5 <<endl;
		alpsolfile << endl;
		massfile << endl;
		conffile << endl;
		asolfile << endl;
		evolutionfile<<endl;

	       cout<<"grr="<<var[0][num-1]<<endl;
               cout<< "mass= "<<((R)/2)* (1-1/pow(var[0][num-1],2))<<" numbos: "<<numbos<< " numfer: "<<numfer<<" omega: "<<omega<<" r99bos: "<<r99bos<<" r99fer: "<<r99fer<<" r99tot: "<<r99tot<<endl;
		cout<<" Nb/Nf: "<<numbos/numfer<< " Eb/Ef: "<<massbos/massfer<<endl;
	       cout<< "mass_bos= "<<massbos<<" mass_fer= "<<massfer<<endl;
		cout<<"compactness= "<<((R)/2)* (1-1/pow(var[0][num-1],2))/r99tot<<endl;
		cout<<"compactness_boson= "<<massbos/r99bos<<endl;
	//	cout<< "mass_corr= "<<((R)/2)* (1-1/pow(var[0][num-1],2))<<endl;
		cout<< "E_binding= "<< ((R)/2)* (1-1/pow(var[0][num-1],2)) - numbos - numfer <<endl;

	if(isotropic)
		 cout<<"r99bos_iso: "<<riso[nbos]<<" r99fer_iso: "<<riso[nfer]<<endl;


               infofile<< "mass= "<<((R)/2)* (1-1/pow(var[0][num-1],2))<< "numbos: "<<numbos<< " numfer: "<<numfer<<" omega: "<<omega<<" r99bos: "<<r99bos<<" r99fer: "<<r99fer<<" r99tot: "<<r99tot<<endl;
		infofile<<" Nb/Nf: "<<numbos/numfer<< " Eb/Ef: "<<massbos/massfer<<endl;
	       infofile<< "mass_bos= "<<massbos<<" mass_fer= "<<massfer<<endl;
		infofile<<"compactness= "<<((R)/2)* (1-1/pow(var[0][num-1],2))/r99tot<<endl;
		infofile<<"compactness_boson= "<<massbos/r99bos<<endl;
	//	infofile<< "mass_corr= "<<((rmax)/2)* (1-1/pow(var[0][nmax],2))<<endl;
	       infofile<< "E_binding= "<< ((R)/2)* (1-1/pow(var[0][num-1],2)) - numbos - numfer <<endl;

	if(isotropic)
		 infofile<<"r99bos_iso: "<<riso[nbos]<<" r99fer_iso: "<<riso[nfer]<<endl;

	       rfile.close();
	       solfile.close();
	       solfile1.close();
	       solfile2.close();
	       solfile3.close();
	       solfile4.close();
	       solfile5.close();
	       alpsolfile.close();
	       massfile.close();
	       conffile.close();
	       asolfile.close();
	       evolutionfile.close();
	       infofile.close();

	}

	numbosfer.close();

	for(int i = 0; i < 5; i++) {
   	  delete [] var[i];
	}
	delete [] var;
	delete [] riso;
	delete [] riso2;
	delete [] beta;
	delete [] conformal;
	delete [] conformal2;
	delete [] bospotential;

	input.close();

	return 0;
}
