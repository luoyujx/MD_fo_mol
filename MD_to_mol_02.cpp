#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <direct.h>
using namespace std;

//------------------------units&constants------------------------//
// time: ps
// distance: nm
// energy: amu (nm^2/ps^2)
// CONST=6.02*10^23*e^2/(4pi*e0) 
#define MAX_NUM 200 // max number of particles
#define MAX_LABEL 20 // max number of lables
#define PI   3.14159282      // Pi
#define E_0  8.854187817e-12 // epsilon_0 [F/m]
#define AMU  1.6605402e-27   // atomic mass unit [kg]
#define K_B  1.380658e-23    // Boltzmann constant[J/K]
#define UNIT_E 1.602177e-19   // elemental charge[C]
#define EPSILON 166.4    // epsilon/k_B for Kr [K]
#define SIGMA 0.365      // LJ parameter sigma [nm] for Kr 
//#define EPSILON 120.0      // LJ parameter epsilon/k_B for Ar
//#define SIGMA 0.341      // LJ parameter sigma [nm] for Ar 
//#define EPSILON 36.7     // LJ parameter epsilon/k_B for Ne  
//#define SIGMA 0.279      // LJ parameter sigma [nm] for Ne 
#define  ALPHA_CGS  2.46e-30       // polarizabiliy [m^3] for Kr, CGS unit,  4*pi*e_0*alpha_CGS[cm^3] =10^6 * alpha_SI[C m^2 / V] 
//#define  ALPHA_CGS  1.62e-30       // polarizabiliy [m^3] for Ar, CGS unit
//#define  ALPHA_CGS  0.390e-30      // polarizabiliy [m^3] for Ne, CGS unit

//---------------parameters of a charged particle---------------//
class denka{
	private:
		int PN;		                        // number of particles (charged + neutral)
		int CN;			                    // number of charge
		long double position[MAX_NUM][3];   // coordinates of i-th atom [nm]
		long double velocity[MAX_NUM][3];   // velocity of i-th atom    [nm/ps]
		long double mass[MAX_NUM];          // mass of i-th atom        [amu]
		long double momentum[MAX_NUM][3];   // momentum of i-th atom    [amu nm/ps]
		int	charge_init[MAX_NUM];		    // initial charge of i-th atom
		long double	charge[MAX_NUM];		// charge of i-th atom (including charge build-up)
		int	Q_limit[MAX_NUM];		        // maxium charge of i-th atom
		long double t_charge_bu;			// charge build-up time
		int random_charge[MAX_NUM];         // target for charge distribution
		int num_random_atom;				// number of target for charge distribution
		int num_fixed_atom;				    // number of particles with constant charge
		string label[MAX_NUM];              // label for identifying particles
		int num_label;				        // number of independent label 
		string valid_label[MAX_NUM];        // list of independent label 
		long double t,t0,tn;                // time 
		int n_step;			                // number of divisions of time
		long double interval_t;				// duration of each division of time = (tn-t0)/n_step
		long double q[MAX_NUM][3][2];       // inaccuracy due to approximation
		long double px, py, pz;				// momentums of the system
		long double refmomentum;			// for accuracy evaluation. ? 
		long double minmomentum;			// for accuracy evaluation. ?
		long double refenergy;				// for accuracy evaluation. ?
		int snap_num;						// how many steps of snapshot to take

		ofstream foutinformation;			// name of configuration file, filename[i] + label
		char InfoFile[32];					// string of the configuration file
		char filename[MAX_LABEL*2][32];		// name of the exported file
		ofstream outfile[MAX_LABEL*2];		// string of the exported file
		char filename_2body[MAX_LABEL][MAX_LABEL][32]; // name of the exported file
		ofstream outfile_2body[MAX_LABEL][MAX_LABEL];  // name of the exported file
		char SnapFile[32];			        // name of the snapshot file
		ofstream fsnapshot;					// file stream of the snapshot files

		char pot_flag;			// interaction model to use default: C, C: Coulomb, L: Lenard-Jones, M: Coulomb+LJ, N: Coulomb + LJ + Induced dipole
		char pot_label[8];		// label for interaction models C: Coulomb, L: Lenard-Jones, M: Coulomb+LJ, N: Coulomb + LJ + Induced dipole

	public:
		int trial;              // how many times to repeat the simulation
		void keisan(int n_th);
		long double kyori(
			long double,long double,long double ,long double,long double,long double
		);
		int dice(int faces);
		void rkg();

		long double force(int i,int j); //???q?????v?Z?B???q i ?? dv/dt ??E?????   j ????W??????
		long double potential_energy(); //?S?|?e???V?????G?l???M?[??v?Z [amu nm^2/ps^2]
		long double kinetic_energy();   //?S?^???G?l???M?[??v?Z [amu nm^2/ps^2]
		long double PE_plus_KE();       //?S?G?l???M?[?iPE+KE?j??v?Z [amu nm^2/ps^2]
		void cal_momentum();    //?????_????d????z????
		void dist_charge();    //?????_????d????z????

		void mk_out_dir(char *makingdirectory);
		void open_output_files(char *makingdirectory);
		void close_output_files(char *makingdirectory);
		int get_data(char *input_file);
		void snap_output(long double snap_time);
		void message();
};

//-------------------------help message-------------------------//
void denka::message(){
	cout << "program.exe parameterfile" << endl; 
	cout << " " << endl; 
	cout << "?p?????[?^?t?@?C??????? " << endl; 
	cout << "?P?s??@?@?J?n????[ps]  ?I??????[ps]  MD?X?e?b?v??[????]  ?X?i?b?v?V???b?g??" << endl; 
	cout << "?@?@      ????????X?y?[?X?????B" << endl; 
	cout << "?@?@      ?X?i?b?v?V???b?g??=0????X?i?b?v?V???b?g??o???????B" << endl; 
	cout << " " << endl; 
	cout << "?Q?s??@?@?S???q???@?S?d????  ?????p?w??  ?d??build-up????[ps]" << endl; 
	cout << "?@?@?@?@?@?????p: default=C: Coulomb, L: LJ, M: C+LJ, N: C+LJ+Induced dipole " << endl; 
	cout << " " << endl; 
	cout << "?R?s??@?@?V?~?????[?V????????s??" << endl; 
	cout << " " << endl; 
	cout << "?S?s???~????????q?f?[?^??? " << endl; 
	cout << "???? " << endl; 
	cout << "???x??  ????[amu]  x[nm] y[nm] z[nm] vx[nm/ps] vy[nm/ps] vz[nm/ps] ?????d?? ???d??" << endl; 
	cout << "???x????T???????B?????]????????????B " << endl; 
	cout << "?????d???u????V?~?????[?V?????B?????d???O???????" << endl; 
	cout << "???d??????q??????_????d???z????B " << endl; 
	exit(0);
}

void denka::mk_out_dir(char *makingdirectory){
	char NewDir[32];
	int i;
	// create a directory according to the name of the configuration file
    i=0;
	do {
		NewDir[i]=makingdirectory[i];
        i=i+1;
        if(i>31) {cout<<"Filename is too long. Be fewer 32 letters."<<endl;exit (1);}
	} 
	while(makingdirectory[i]!='.');
    NewDir[i]=NULL;
	cout<<"simulation result -> "<<NewDir<<endl;
    _mkdir(NewDir); // create the directory
}

// create the output directory
void denka::open_output_files(char *makingdirectory){
	char DirName[32];
	int i,j;
	i=0;
	
	do {
		DirName[i]=makingdirectory[i];
		i=i+1;
		if(i>31){cout<<"Filename must be fewer than 32 letters."<<endl;exit (1);}
		}
	while(makingdirectory[i]!='.');	
	
	DirName[i]=NULL;
	// make an restructured version of the configuration file
	strcpy_s(InfoFile, DirName);
	strcat_s(InfoFile, "\\Infomation.txt");
	foutinformation.open(InfoFile);
	cout << InfoFile << endl;
	if(!foutinformation){cout<<"cannot make Infomation FILE. \n";}
	foutinformation<<"START TIME T0 "<<t0<<" ps"<<endl;
	foutinformation<<"STOP TIME Tn "<<tn<<" ps"<<endl;
	foutinformation<<"kizami n"<< n_step <<endl;
	foutinformation<<"trial number "<<trial<<endl;
	foutinformation<<"Number of simulated particles is "<<PN<<endl;
	foutinformation<<"Number of charge is "<<CN<<endl;
	foutinformation<<"Interaction : "<< pot_label <<endl;
	foutinformation<<"New directory name is "<<DirName<<endl;
	foutinformation<<"N"<<DirName<<endl;
	foutinformation<<endl;
	foutinformation<<endl;
	foutinformation<<endl;
	// open the restructured version of the configuration file
	if(MAX_LABEL < num_label){cout << "too much label!" <<endl;}
	
	for(i=0; i<num_label;i++){
		strcpy_s(filename[2*i]  ,DirName);			// label i-th ion
		strcpy_s(filename[2*i+1],DirName);			// label i-th neutral 
		strcat_s(filename[2*i]  ,"\\");
		strcat_s(filename[2*i+1],"\\");
		strcat_s(filename[2*i]  ,valid_label[i].c_str());
		strcat_s(filename[2*i+1],valid_label[i].c_str());
		strcat_s(filename[2*i]  ,"_P.dat");
		strcat_s(filename[2*i+1],"_N.dat");
		// write the file
		outfile[2*i].open(  filename[2*i]);
		if(!outfile[2*i]) {cout << "cannot make \n ";}
		outfile[2*i+1].open(filename[2*i+1]);
		if(!outfile[2*i+1]) {cout << "cannot make \n ";}
	}

	for(i=0;i<num_label;i++){
		for(j=i;j<num_label;j++){
			strcpy_s(filename_2body[i][j]  ,DirName);			// label i-th ion
			strcat_s(filename_2body[i][j]  ,"\\");
			strcat_s(filename_2body[i][j]  ,"2body_");
			strcat_s(filename_2body[i][j]  ,valid_label[i].c_str());
			strcat_s(filename_2body[i][j]  ,valid_label[j].c_str());
			strcat_s(filename_2body[i][j]  ,".dat");
			// write the file
			outfile_2body[i][j].open(filename_2body[i][j]);
			if(!outfile_2body[i][j]){cout << "cannot make \n";}
		}
	}

	if(snap_num>0){
		strcpy_s(SnapFile ,DirName);		
		strcat_s(SnapFile,"\\Snapshot.dat");
		fsnapshot.open(SnapFile);
		cout << SnapFile << endl;
		if(!fsnapshot){cout<<" cannot make Infomation FILE. \n";}
	}
}

void denka::close_output_files(char *makingdirectory){
	int i;
	for(i=0; i < num_label; i++){
		outfile[2*i].close(); outfile[2*i+1].close();
	}
	foutinformation.close();
	if(snap_num > 0){
		fsnapshot.close();
	}
}

void denka::keisan(int n_th) {
	int i,j,k,temp_i;
	int suffix_1,suffix_2;
	long double firstP[3] = {0, 0, 0}; // return the sum of the initial momentum 
	long double firstenergy = 0;
	long double finalenergy = 0;
	double Vx1,Vy1,Vz1;
	double Vx2,Vy2,Vz2;
	double InnerProduct;
	long double ev[MAX_NUM]; 	       // return the kinetic energy of each particle
	int	file_num;

    cout<<"you entered... "<<endl;
	cout<<"START TIME T0 "<<t0<<" ps"<<endl;
	cout<<"STOP TIME Tn "<<tn<<" ps"<<endl;
	cout<<"kizami n "<< n_step <<endl;
	cout<<"trial number "<< trial <<endl;
	cout<<"Number of simulated particles is "<<PN<<endl;
	cout<<"Number of total charge is "<<CN<<endl;

	cal_momentum();
	cout<<"initial total momentum of system : [amu nm/ps]"<<endl<<
	"px="<<px<<" py="<<py<<" pz="<<pz<<endl;
	// calculate the sum of the initial energy
	// the output is converted from [amu nm^2/ps^2] to [eV]
	cout<<"first potential energy  "<< (AMU*1e6/UNIT_E)*potential_energy() <<"  eV"<<endl;  
	cout<<"first kinetic energy  "<< (AMU*1e6/UNIT_E)*kinetic_energy() <<"  eV"<<endl;
	firstenergy=PE_plus_KE();  ////?P??? [amu nm^2/ps^2]
	cout<<"first total energy  "<< (AMU*1e6/UNIT_E)*firstenergy<<"  eV"<<endl;		 
    // loop of Runge Kutta Gill for every step
	for(i = 0;i <= n_step; i++){
		t = t0 + i*interval_t;				
		for(i=0;i<PN;i++){
			charge[i]=charge_init[i]*exp(-t/t_charge_bu);
		}
        rkg();
		if((snap_num!=0)&&(i%snap_num==0)){
			snap_output(t);
		}
		// calculate the energy
   		for(k=0;k<PN;k++){
			ev[k]=0.5*(AMU*1e6/UNIT_E)*mass[k]*(velocity[k][0]*velocity[k][0] + velocity[k][1]*velocity[k][1]+velocity[k][2]*velocity[k][2]);
		}; // the output is converted from [amu nm^2/ps^2] to [eV]
	}
	// export the snapshot after the last step finishes
	if((snap_num != 0) && (n_step%snap_num != 0) ){
		snap_output(t);
	} 
	cal_momentum();
	cout<<"total momentum / [amu nm/ps]"<<endl<<
	"px="<<px<<" py="<<py<<" pz="<<pz<<endl;
	// calculate the sum of the initial kinetic energy
	cout<<"final potential energy  "<< (AMU*1e6/UNIT_E)*potential_energy() <<"  eV"<<endl; 
	cout<<"final kinetic energy  "<< (AMU*1e6/UNIT_E)*kinetic_energy() <<"  eV"<<endl;  
	finalenergy=PE_plus_KE();  // unit: [amu nm^2/ps^2]
	cout<<"final total energy  "<< (AMU*1e6/UNIT_E)*finalenergy<<"  eV"<<endl;
	if(refenergy<(fabs((finalenergy-firstenergy)/finalenergy))) {
		refenergy = (fabs((finalenergy-firstenergy)/finalenergy));
	}
	if(refmomentum<fabs(kyori(px,firstP[0],py,firstP[1],pz,firstP[2])/minmomentum)){
		refmomentum=fabs(kyori(px,firstP[0],py,firstP[1],pz,firstP[2])/minmomentum);
	}
    cout<<"seido (first total momentum-final total momentum)/(minmomentum) "<<refmomentum<<endl;
	cout<<"seido max(first-final/final)energy = "<<refenergy<<endl;
    foutinformation<<endl;
	foutinformation<<"momentum seido "<<refmomentum<<endl;
	foutinformation<<"energy seido "<<refenergy<<endl;
	// write the results to the exported file
	// the results are sorted by the lables
	for(i=0; i<PN; i++){
        if(charge_init[i]==0){
			for(j=0; j<num_label; j++){
				if(label[i] == valid_label[j]){file_num = 2*j + 1;};
			}
		}
        if(charge_init[i]>0){
			for(j=0; j<num_label; j++){
				if(label[i]==valid_label[j]){file_num=2*j;};
			}
		}
		outfile[file_num] <<n_th<< " " << charge_init[i]<< " " << ev[i]<< " " << momentum[i][0] << " " << momentum[i][1] << " " << momentum[i][2]<< " "<< "0"<< endl;
		foutinformation  <<n_th<< " "  <<label[i]<< " " << charge_init[i]<< " " << ev[i]<< " " << momentum[i][0] << " " << momentum[i][1] << " " << momentum[i][2]<< " "<< "0"<< endl;
		if(charge_init[i]>0){
			for(j=i+1;j<PN;j++){
				if(charge_init[j]>0) {
					for(k=0;k<num_label;k++){
						if(label[i]==valid_label[k]){suffix_1=k;};
						if(label[j]==valid_label[k]){suffix_2=k;};
					}
					if(suffix_1>suffix_2){temp_i= suffix_1; suffix_1=suffix_2; suffix_2=temp_i;}
					Vx1= momentum[i][0]; Vy1= momentum[i][1]; Vz1= momentum[i][2]; 
					Vx2= momentum[j][0]; Vy2= momentum[j][1]; Vz2= momentum[j][2]; 
					InnerProduct = (Vx1*Vx2 + Vy1*Vy2 + Vz1*Vz2)/sqrt((Vx1*Vx1 + Vy1*Vy1 + Vz1*Vz1)*(Vx2*Vx2 + Vy2*Vy2 + Vz2*Vz2));
					outfile_2body[suffix_1][suffix_2]  <<n_th<< " "  << label[i]  << " "  << charge_init[i]  << " "  << label[j]  << " "  << charge_init[j]  << " "  << InnerProduct  << " "  << ev[i]  << " " << ev[j]  << " " << endl;
				}
			}
		}
	}
}

// calculated the momentums
void denka::cal_momentum(){
	int i,j;
	// miltiply the mass with the velocity to get the momentum
	for(i=0;i<PN;i++){
		for(j=0;j<=2;j++){momentum[i][j]=mass[i]*velocity[i][j];};
	}
	px=0; py=0; pz=0;
	// sum of the momentums
	for(i=0;i<PN;i++) {
		px=px+momentum[i][0]; py=py+momentum[i][1]; pz=pz+momentum[i][2];
		// evaluate the accuracy of the momentums
		if(i==0){
			minmomentum=kyori(momentum[i][0],0,momentum[i][1],0,momentum[i][2],0);
		}
		else if(i){
			if(minmomentum>kyori(momentum[i][0],0,momentum[i][1],0,momentum[i][2],0)){
				minmomentum=kyori(momentum[i][0],0,momentum[i][1],0,momentum[i][2],0);
			}
		}
	}
}

void denka::snap_output(long double snap_time){
	int i;
	fsnapshot << " simlation time : " << snap_time << " [ps]" << endl;
	for(i=0;i<PN;i++){
		fsnapshot	<< label[i] << "     "  << mass[i]  
					<< "     "  << position[i][0] << "     "  << position[i][1] << "     "  << position[i][2]
					<< "     "  << velocity[i][0] << "     "  << velocity[i][1] << "     "  << velocity[i][2]
					<< "     "  << charge[i] << endl;
	};
	fsnapshot	<< endl;  

}

long double denka::kyori(long double a1,long double b1,long double a2,long double b2,long double a3,long double b3){
	return sqrt((a1-b1)*(a1-b1)+(a2-b2)*(a2-b2)+(a3-b3)*(a3-b3));
}

// distance between two points
void denka::rkg(){
    static long double p1=0.2928932188134524, p2=1.707106781186547; // p1 = 1 - sqrt(1/2), p2 = 1 + sqrt(1/2) 
	long double ck0[MAX_NUM][3][2];
	long double q1[MAX_NUM][3][2];
	long double r1[MAX_NUM][3][2];
	long double y1[MAX_NUM][3][2];
    int i, j, k;
	
	for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				ck0[i][j][k]=0; q1[i][j][k]=0;  r1[i][j][k]=0;y1[i][j][k]=0;
			};
		};
	}
    // i for ions, j for coordinations
	// start the 4-step-calculation
	// step 1
	for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				if(k){
					ck0[i][j][k]=force(i,j)*interval_t/2;
				}
	   			else{
					ck0[i][j][k]=velocity[i][j]*interval_t/2;
				};
	   		};
		};
	}
    for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				if(k){
					y1[i][j][k]=velocity[i][j]+ck0[i][j][k]-q[i][j][k];
				}
				else{
					y1[i][j][k]=position[i][j]+ck0[i][j][k]-q[i][j][k];
				}
			};
		};
	}
	for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				if(k){
					r1[i][j][k]=y1[i][j][k]-velocity[i][j];
				}
	   			else{
					r1[i][j][k]=y1[i][j][k]-position[i][j];
				}
			};
		};
	}
	for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				q[i][j][k]=q[i][j][k]+r1[i][j][k]*3-ck0[i][j][k];
			};
		};
	}

    for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				if(k){
					velocity[i][j]=y1[i][j][k];
				}
	   			else{
					position[i][j]=y1[i][j][k];
				}
			};
		};
	}
	// step 2
    for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				if(k){
					ck0[i][j][k]=force(i,j)*interval_t;
				}
	   			else{
					ck0[i][j][k]=velocity[i][j]*interval_t;
				}
	   		};
		};
	}
    for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				if(k){
					y1[i][j][k]=velocity[i][j]+p1*(ck0[i][j][k]-q[i][j][k]);
				}
	   			else{
					y1[i][j][k]=position[i][j]+p1*(ck0[i][j][k]-q[i][j][k]);
				}
			};
		};
	}
    for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				if(k){
					r1[i][j][k]=y1[i][j][k]-velocity[i][j];
				}
				else{
					r1[i][j][k]=y1[i][j][k]-position[i][j];
				}
			};
		};
	}
	for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				q[i][j][k]=q[i][j][k]+r1[i][j][k]*3-ck0[i][j][k]*p1;
			};
		};
	}
	for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				if(k){
					velocity[i][j]=y1[i][j][k];
				}
				else{
					position[i][j]=y1[i][j][k];
				}
			};
		};
	}

	// step 3
	for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				if(k){
					ck0[i][j][k]=force(i,j)*interval_t;
				}
				else{
					ck0[i][j][k]=velocity[i][j]*interval_t;
				}
			};
		};
	}       
	for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				if(k){
					y1[i][j][k]=velocity[i][j]+p2*(ck0[i][j][k]-q[i][j][k]);
				}
				else{
					y1[i][j][k]=position[i][j]+p2*(ck0[i][j][k]-q[i][j][k]);
				}
			};
		};
	}
	for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				if(k){
				r1[i][j][k]=y1[i][j][k]-velocity[i][j];
				}
				else{
					r1[i][j][k]=y1[i][j][k]-position[i][j];
				}
			};
		};
	}
	for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				q[i][j][k]=q[i][j][k]+r1[i][j][k]*3-ck0[i][j][k]*p2;
			};
		};
	}
	for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				if(k){
					velocity[i][j]=y1[i][j][k];
				}
				else{
					position[i][j]=y1[i][j][k];
				}
			};
		};
	}

    // step 4
	// xx = t + interval_t / 2; 
	for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				if(k){
					ck0[i][j][k]=force(i,j)*interval_t/2;
				}
				else{
					ck0[i][j][k]=velocity[i][j]*interval_t/2;
				}
			};
		};
	}       
	for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				if(k){
					y1[i][j][k]=velocity[i][j]+(ck0[i][j][k]-q[i][j][k])/3;
				}
				else{
					y1[i][j][k]=position[i][j]+(ck0[i][j][k]-q[i][j][k])/3;
				}
			};
		};
	}
	for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				if(k){
					r1[i][j][k]=y1[i][j][k]-velocity[i][j];
				}
				else{
					r1[i][j][k]=y1[i][j][k]-position[i][j];
				};
			};
		};
	}
	for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				q[i][j][k]=q[i][j][k]+r1[i][j][k]*3-ck0[i][j][k];
			};
		};
	}
	for(i=0;i<=PN-1;i++){
		for(j=0;j<=2;j++){
			for(k=0;k<=1;k++){
				if(k){
					velocity[i][j]=y1[i][j][k];
				}
				else{
					position[i][j]=y1[i][j][k];
				};
			};
		};
	}
}


//Runge Kutta Gill calculation
int denka::dice(int faces){ 	
	return(1 + rand()%faces);
}  
// n-face dice return random interger from 1 to n

void denka::dist_charge() {
	int i,num,max_charge,fixed_charge;         
		srand( (unsigned)time(NULL) );
		max_charge=0; fixed_charge=0;
		for(i=0;i<PN;i++){max_charge=max_charge+Q_limit[i];
				fixed_charge=fixed_charge+charge_init[i];
		}
		if(CN > max_charge){
			cout << " charge number mismatch" << endl; exit(0);
		}
		for(i=0;i<(CN - fixed_charge);i++){
			num = dice(PN)-1;
			while( charge_init[num]>=Q_limit[num]){
				num = dice(PN)-1;
			}
			charge_init[num]=charge_init[num]+1;
		}	// randomly distribute the charge

}

// calculate the force, calculate the dv/dt of ith particle   
// j: x, y, z coordination energy: [amu nm^2/ps^2] force: [amu nm/ps^2]
long double denka::force(int i ,int j){
	int seisu=0;
	long double nagasa=0;
	long double uhen=0;
	switch(pot_flag){
		case 'C'	: 
			for (seisu=0;seisu<PN;seisu++)  {
				nagasa=kyori(position[i][0],position[seisu][0],position[i][1],position[seisu][1],position[i][2],position[seisu][2]);
				if(seisu-i) {uhen=uhen+(UNIT_E*UNIT_E/(4*PI*E_0)*1e3)*charge[i]*charge[seisu]*(position[i][j]-position[seisu][j])/pow(nagasa,3) /(AMU*mass[i]);}
			};  // couloumb interaction
			break;
		case 'L'	:
			for (seisu=0;seisu<PN;seisu++)  {
				nagasa=kyori(position[i][0],position[seisu][0],position[i][1],position[seisu][1],position[i][2],position[seisu][2]);
					if(seisu-i) {
					uhen=uhen+(4*(EPSILON * K_B *1e-6)*(-6*pow(SIGMA/nagasa,6)+12*pow(SIGMA/nagasa,12)))*(position[i][j]-position[seisu][j])/pow(nagasa,2) /(AMU*mass[i]);
				};
			};  //  LJ interaction between neutral ions, couloumb interaction ignored 
			break;
		case 'M' : 
			for (seisu=0;seisu<PN;seisu++)  {
				nagasa=kyori(position[i][0],position[seisu][0],position[i][1],position[seisu][1],position[i][2],position[seisu][2]);
					if(seisu-i) {
					uhen=uhen+((UNIT_E*UNIT_E/(4*PI*E_0)*1e3)*charge[i]*charge[seisu]/nagasa + 4*(EPSILON * K_B *1e-6)*(-6*pow(SIGMA/nagasa,6)+12*pow(SIGMA/nagasa,12)))
							*(position[i][j]-position[seisu][j])/pow(nagasa,2) /(AMU*mass[i]);
				};
			};  // couloumb + LJ interactions
			break;
		case 'N'	:
			for (seisu=0;seisu<PN;seisu++)  {
				nagasa=kyori(position[i][0],position[seisu][0],position[i][1],position[seisu][1],position[i][2],position[seisu][2]);
					if(seisu-i) {
						uhen=uhen+ ( (UNIT_E*UNIT_E/(4*PI*E_0)*1e3)*charge[i]*charge[seisu]/nagasa 
								+ 4*(EPSILON * K_B *1e-6)*(-6*pow(SIGMA/nagasa,6)+12*pow(SIGMA/nagasa,12))
								- ( ALPHA_CGS * UNIT_E*UNIT_E/(8*PI*E_0) *1e30 )*charge[i]/pow(nagasa,4) )
								* (position[i][j]-position[seisu][j])/pow(nagasa,2) /(AMU*mass[i]);
					};
			};  // couloumn + LJ + ion exitation dipole interactions
			break;
		default	: 
			for (seisu=0;seisu<PN;seisu++)  {
				nagasa=kyori(position[i][0],position[seisu][0],position[i][1],position[seisu][1],position[i][2],position[seisu][2]);
					if(seisu-i) {
					uhen=uhen+(UNIT_E*UNIT_E/(4*PI*E_0)*1e3)*charge[i]*charge[seisu]*(position[i][j]-position[seisu][j])/pow(nagasa,3) /(AMU*mass[i]);
					};
			};  // dafult: couloumb interaction only
			break;
	}
    return uhen;
}

// total energy
long double denka::PE_plus_KE(){
	return potential_energy() + kinetic_energy();
} // total energy
// total potential energy

long double denka::potential_energy(){
	int i,j;
	long double nagasa=0;
	long double t_energy=0;
	switch(pot_flag){
		case 'C'	: 
			for(i=0;i<PN;i++) {
				for(j=i+1;j<PN;j++){
					nagasa=kyori(position[i][0],position[j][0],position[i][1],position[j][1],position[i][2],position[j][2]);
					t_energy=t_energy+(UNIT_E*UNIT_E/(4*PI*E_0*AMU)*1e3)*charge[i]*charge[j]/nagasa;
				};	
			};  
			break;
		case 'L'	: 
			for(i=0;i<PN;i++) {
				for(j=i+1;j<PN;j++){
					nagasa=kyori(position[i][0],position[j][0],position[i][1],position[j][1],position[i][2],position[j][2]);
					t_energy=t_energy+4* (EPSILON * K_B *1e-6/AMU)*(pow(SIGMA/nagasa,12)-pow(SIGMA/nagasa,6)) ;
				};	
			};  
			break;
		case 'M'	: 
			for(i=0;i<PN;i++) {
				for(j=i+1;j<PN;j++){
					nagasa=kyori(position[i][0],position[j][0],position[i][1],position[j][1],position[i][2],position[j][2]);
					t_energy=t_energy+(UNIT_E*UNIT_E/(4*PI*E_0*AMU)*1e3)*charge[i]*charge[j]/nagasa
									+4*(EPSILON * K_B *1e-6/AMU)*(pow(SIGMA/nagasa,12)-pow(SIGMA/nagasa,6)) ;
				};	
			};  
			break;
		case 'N'	: 
			for(i=0;i<PN;i++) {
				for(j=i+1;j<PN;j++){
					nagasa=kyori(position[i][0],position[j][0],position[i][1],position[j][1],position[i][2],position[j][2]);
					t_energy=t_energy+(UNIT_E*UNIT_E/(4*PI*E_0*AMU)*1e3)*charge[i]*charge[j]/nagasa
									+4*(EPSILON * K_B *1e-6/AMU)*(pow(SIGMA/nagasa,12)-pow(SIGMA/nagasa,6)) 
									- (ALPHA_CGS * UNIT_E * UNIT_E * 1e30/(8*PI*E_0*AMU) )*charge[i]/pow(nagasa,4) ;
				};	
			};  
			break;
		default	: 
			for(i=0;i<PN;i++) {
				for(j=i+1;j<PN;j++){
					nagasa=kyori(position[i][0],position[j][0],position[i][1],position[j][1],position[i][2],position[j][2]);
					t_energy=t_energy+(UNIT_E*UNIT_E/(4*PI*E_0*AMU)*1e3)*charge[i]*charge[j]/nagasa;
				};	
			};  
			break;
	}
	return t_energy;
} 
					
// total kinetic energy
long double denka::kinetic_energy(){
	int i;
	long double t_energy=0;
	for(i=0;i<PN;i++){
		t_energy=t_energy+0.5*mass[i]*(velocity[i][0]*velocity[i][0]+velocity[i][1]*velocity[i][1]+velocity[i][2]*velocity[i][2]);
	}
	return t_energy;
} 

int denka::get_data(char *input_file){

	char buf[512];
	int  i,j,k,LINE_NUM,flag;

	ifstream fin(input_file);

	if(!fin){cout << "cannot open \n ";
			return 1;
    }

	LINE_NUM=0;

		while( fin.getline( buf, sizeof(buf) )){
			LINE_NUM=LINE_NUM+1;
			istringstream is( buf );
			if(LINE_NUM==1){ is >> t0 >> tn >> n_step>> snap_num ;}
			if(LINE_NUM==2){ is >> PN >> CN >> pot_flag >> t_charge_bu;}
			if(LINE_NUM==3){ is >> trial ;}
			if(LINE_NUM >3){ is >> label[LINE_NUM-4] >> mass[LINE_NUM-4] 
								>> position[LINE_NUM-4][0] >> position[LINE_NUM-4][1] >> position[LINE_NUM-4][2] 
								>> velocity[LINE_NUM-4][0] >> velocity[LINE_NUM-4][1] >> velocity[LINE_NUM-4][2] 
								>> charge[LINE_NUM-4] >> Q_limit[LINE_NUM-4] ;
			}
		}
	fin.close();

	interval_t=(tn-t0)/n_step;	

	switch(pot_flag){
		case 'C'	: strcpy_s(pot_label, "Coulomb");	break;
		case 'L'	: strcpy_s(pot_label, "LJ");	break;
		case 'M'	: strcpy_s(pot_label, "C+LJ");	break;
		case 'N'	: strcpy_s(pot_label, "C+LJ+ID");	break;
		default		: strcpy_s(pot_label, "Coulomb");	break;
	}

	for(j=0;j<3;j++) {for(i=0;i<PN;i++) {momentum[i][j]=0;} ;} // reset
	for(i=0;i<PN;i++) {for(j=0;j<3;j++)  {for(k=0;k<2;k++) {q[i][j][k]=0;};};} // reset the inaccuracy

	valid_label[0] = label[0];
	num_label=1;
	for(i=1;i<PN;i++){
		flag=0;
		for(j=0;j<num_label;j++){ if(label[i] == valid_label[j]){ flag=1;} ;}
		if(flag==0){valid_label[num_label] = label[i]; num_label=num_label+1 ;}
	} // check the total number of labels, the number of non-duplicated lables

	cout << LINE_NUM-3 << " particle's data. " << PN << " particles are valid." << endl;
	num_random_atom=0;
	num_fixed_atom=0;
	for(i=0;i<PN;i++){
		if(charge[i]==-1){random_charge[num_random_atom]=i;num_random_atom=num_random_atom+1;};
		if(charge[i]!=-1){num_fixed_atom=num_fixed_atom+1;};
		}
	return 0;
}
int main(int argc,char *argv[]){
	int i;
	cout<<"********************************PROGRAM START********************************"<<endl;
	cout<<"length:nm,time:ps,energy:eV"<<endl;    
	denka sim_data;				// initialized the class
	if(argv[1]==NULL){          // display the help messange if there is no configuration file assigned
		sim_data.message();
	}
	sim_data.mk_out_dir(argv[1]); // make the output directory
	sim_data.get_data(argv[1]);	  // get the initial files
	sim_data.open_output_files(argv[1]); // make the output files

	for(i=0;i<sim_data.trial;i++){
		cout<<"******"<<i+1<<"th simulation"<< "****************************************************"<<endl;
		sim_data.dist_charge();       // distribute the charge randomly
		sim_data.keisan(i+1);         // n-step simulation for Runge-Kutta-Gill calculation
	}
	sim_data.close_output_files(argv[1]);
}

