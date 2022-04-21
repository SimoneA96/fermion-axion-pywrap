#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <string>
#include <sys/stat.h>
using namespace std;
int num;
long double y[5], y_axi[5], diff[5], F[5];
long double** var;
long double k1[5], k2[5], k3[5], k4[5];
const long double Pi=3.14159265358979323846L, m=1.0, KK=100, polind = 2.;
long double fa_axion;
long double R, Rmax;
long double eps;
long double Cutphi, CutP;
long double omega2;
long double delta_r, r;
int sign (long double z) {
    if (z>0) return 1;
    if (z==0) return 0;
    if (z<0) return -1;
}
long double V(long double phi) {
    return 2.0 * pow(m,2) * pow(fa_axion,2) / 0.22 * (1.0 - sqrt(1.0 - 4 * 0.22 * pow(sin( phi/(2.0*fa_axion)),2)));
}
long double dV(long double phi) {
     if(phi!=0)
         return 2.0 * pow(m,2)  * fa_axion * sin(phi/(2.0*fa_axion)) * cos(phi/(2.0*fa_axion)) /
	    (phi * sqrt(1.0 - 4.0 * 0.22 * pow(sin(phi/(2.0*fa_axion)),2) ));
     else
         return 0;
}
int func(long double r,  long double omega2, long double h) {
    long double a=y_axi[0];
    long double alp=y_axi[1];
    long double phi=y_axi[2];
    long double dphi=y_axi[3];
    long double P=y_axi[4];
    long double rho_tmp;
//OSS: eps=K*rho --> eps=10*sqrt(P) eps = pow(KK,1/polind)*pow(P,(polind-1)/polind)/(polind-1)
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
    func(r, omega2, h);
    for (int i=0; i<=4; i++) k1[i]=h*F[i];
    for (int i=0; i<=4; i++) y_axi[i]=y[i]+0.5*k1[i];
//control for P not to become negative
        if (y_axi[4]<=0) {
		y_axi[4]=0;
		k1[4]=0;
		}
    func(r+half_h, omega2, h);
    for (int i=0; i<=4; i++) k2[i]=h*F[i];
    for (int i=0; i<=4; i++) y_axi[i]=y[i]+0.5*k2[i];
        if (y_axi[4]<=0) {
		y_axi[4]=0;
		k2[4]=0;
		}
    func(r+half_h, omega2, h);
    for (int i=0; i<=4; i++) k3[i]=h*F[i];
    for (int i=0; i<=4; i++) y_axi[i]=y[i]+k3[i];
        if (y_axi[4]<=0) {
		y_axi[4]=0;
		k3[4]=0;
		}
    func(r+h, omega2, h);
    for (int i=0; i<=4; i++) k4[i]=h*F[i];
    for (int i=0; i<=4; i++) diff[i]=(k1[i]+k4[i])/6 + (k2[i]+k3[i])/3;
    return 0;
}

int ode_int(long double phi0c, long double rho0c, long double omega, long double* a){
    int n=0;
    y[0]=1;
    y[1]=1;
    y[2]=phi0c;
    y[3]=0;
    y[4]=KK*pow(rho0c,polind);
    delta_r=R/(num);
    omega2=pow(omega,2);
    r=delta_r;
    while (r<=R + delta_r) {
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

	//Setting to 0 the pressure when it starts to do steps (for low values of P)
        if (y[4] == var[4][n-1] || var[4][n-1]==0) y[4]=0;

        if (y[2] > Cutphi || y[2] < 0 || y[4]> CutP || diff[2]>0) {
		Rmax = r;
		break;
		}
        r+=delta_r;
    }
    *a = y[0];
    return sign(y[2]);
}

int main(int argc, char *argv[]){
	int input_num;
   	long double order=-19, phi0c, rho0c;
	long double grr, mass;
	if (argc<3){
	  cout<<"You should run './Shoot data fa_exp' "<<'\n';
	  cout<<"where data is the folder where you want to put your output"<<endl;
	  cout<<"and fa_exp is fa_axion = pow(10,fa_exp) for the axion potential"<<endl;
	  return 0;
	}

   	ifstream  input;
   	ofstream  output;
	string folder;
	folder = argv[1];
	char* folder2;
	
	folder2 = argv[1];
	fa_axion = pow(10,atof(argv[2]));
	mkdir( folder2, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0;
   	input.open("Input.txt", ios::in);
   	output.open(folder + "/Output.txt", ios::out);
   	input >> order >> num >> input_num;

	var = new long double*[5];
	for (int i = 0;i<5;i++)
		var[i] = new long double[num];

   	output << input_num <<" "<< num <<endl;
   	output << setprecision(2-order);
   	cout << setprecision(2-order); 


   	if (order<0) {int j=(-1)*order;
      		order=1;
      		for (int i=0; i<j; i++) {
			order=order*0.1;
      		}
      		cout << order << "\n";
   	}
   	for (int j=0; j<input_num; j++) {
     		input >> phi0c >> rho0c >> R;
     		cout << j <<" "<<phi0c<<" "<<rho0c<< endl;
		Rmax = R;
		//I set here the cuts for phi and P..if the solution goes over these values it breaks. To avoid having nans.
		Cutphi = 1.1*phi0c;
		CutP = 1.1* KK*pow(rho0c,polind);
     		long double omega_1=1.0*m,omega_2=2.0*m,omega_mid;
     		long double f_1=ode_int(phi0c, rho0c, omega_1, &grr);
     		long double f_2=ode_int(phi0c, rho0c, omega_2, &grr);
     		long double f_mid;
     		while (f_1*f_2>0) {
        		omega_2=omega_2*1.1;
        		omega_1=omega_1/1.1;
        		f_1=ode_int(phi0c, rho0c, omega_1, &grr);
        		f_2=ode_int(phi0c, rho0c, omega_2, &grr);
			//control for nans
			if (f_1*f_2 != 1 && f_1*f_2 != -1 && f_1*f_2!=0){
				cout<<"nan detected;"<<'\n';
				cout<<"omega_1= "<<omega_1 <<" omega_2= "<<omega_2<<endl; 
				cout<<"f_1= "<<f_1 <<" f_2= "<<f_2<<endl; 
			}
			// end control
     		}
     		while((omega_2-omega_1)/omega_1>2*order) {
        		omega_mid=(omega_1+omega_2)/2;
        		cout << omega_mid << "  ";
        		f_mid=ode_int(phi0c, rho0c, omega_mid, &grr);
			
			//control for nans
			if (f_mid != 1 && f_mid != -1 && f_1*f_2!=0){
				cout<<"nan detected;"<<'\n';
				cout<<"omega_mid= "<<omega_mid <<endl; 
				cout<<"f_mid= "<<f_mid<<endl; 
			}
			// end control
        		if (f_mid>0){
				if(f_1<0){
					omega_2=omega_mid;
					f_2=f_mid;
				}
				else {
					omega_1=omega_mid;
					f_1=f_mid;
				}
			}
			else{
				if(f_1<0){
					omega_1=omega_mid;
					f_1=f_mid;
				}
				else {
					omega_2=omega_mid;
					f_2=f_mid;
				}
			}

		}
		cout<<endl;

		mass = (Rmax/2.) * (1-1/pow(grr,2));


		output <<  phi0c << " " << KK*pow(rho0c,polind) << " " << omega_mid <<" "<<R<< endl;
		cout<<"grr= "<<grr<<" mass= "<<mass<<endl;
	}

	for (int i=0;i<5;i++)
		delete [] var[i];
	delete [] var;
	input.close();
	output.close();
	return 0;
}
