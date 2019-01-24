// MD_fo_mol.cpp : �R���\�[�� �A�v���P�[�V�����̃G���g�� �|�C���g���`���܂��B
//

#include "stdafx.h"


/*���ϗ���
2008/1/24	��؃v���O������RKG���g���āA�Ƃ肠���������B
			�ϐ����𐮗����āA�\������i�߂�B�p�����[�^�t�@�C���̏�����ύX�B
			���q��]��TOF�APIPICO�쐬�@�\�͍폜�B�i�����炭���̃v���O�����ɕ����������ǂ��B�j
			���ݍ�p��Coulomb�ALJ�AID���܂ށB
2008/1/25	snapshot�o�͋@�\��t���B
			�ϐ��̌^�𐮗��B�i���ʂ�double�^��int�ɕύX�j
			�G�l���M�[�ۑ��Ɖ^���ʕۑ��̊m�F���[�`���͍ă`�F�b�N���邱�ƁI
2012/10/~   �����C�I�������������B
			���́A�o�̓t�@�C���t�H�[�}�b�g�ύX�B�S�Ă̌��ʂ�information.txt�ɓ���邱�Ƃɂ���B
*/ //���ϗ����I���


//CoulombExplosionMomentum.cpp
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

//#define CONST 138.6 
//���Ԃ̒P�ʂ� ps �����̒P�ʂ� nm �G�l���M�[�̒P�ʂ� amu nm^2/ps^2
//CONST=6.02*10^23*e^2/(4pi*e0)        ��ؒ萔

#define MAX_NUM 200 
//�ő嗱�q��

#define MAX_LABEL 20 
//�ő�label��

#define PI   3.14159282      //Pi
#define E_0  8.854187817e-12 //epsilon_0 [F/m]
#define AMU  1.6605402e-27   //atomic mass unit [kg]
#define K_B  1.380658e-23    //Boltzmann constant[J/K]
#define UNIT_E 1.602177e-19   //elemental charge[C]
#define EPSILON 166.4    //epsilon/k_B for Kr [K]
#define SIGMA 0.365      //LJ parameter sigma [nm] for Kr 
//#define EPSILON 120.0      //LJ parameter epsilon/k_B for Ar
//#define SIGMA 0.341      //LJ parameter sigma [nm] for Ar 
//#define EPSILON 36.7     //LJ parameter epsilon/k_B for Ne  
//#define SIGMA 0.279      //LJ parameter sigma [nm] for Ne 
#define  ALPHA_CGS  2.46e-30       //polarizabiliy [m^3] for Kr, CGS unit,  4*pi*e_0*alpha_CGS[cm^3] =10^6 * alpha_SI[C m^2 / V] 
//#define  ALPHA_CGS  1.62e-30       //polarizabiliy [m^3] for Ar, CGS unit
//#define  ALPHA_CGS  0.390e-30      //polarizabiliy [m^3] for Ne, CGS unit


//�d�ׂ̌��
class denka
{
	private:
	int PN;		                        // number of simulated atom (ion + neutral)
	int CN;			                    // number of charge
	long double position[MAX_NUM][3];   // coordinates of i-th atom    [nm]
	long double velocity[MAX_NUM][3];   // velocity of i-th atom      [nm/ps]
	long double mass[MAX_NUM];          // mass of i-th atom          [amu]
	long double momentum[MAX_NUM][3];   // momentum of i-th atom      [amu nm/ps]
  	int	charge[MAX_NUM];		        // charge of i-th atom
  	int	Q_limit[MAX_NUM];		        // �ő�charge of i-th atom
  	int random_charge[MAX_NUM];         // �d�ׂ������_���ɕ��z���錴�q�̔ԍ������Ă���
	int num_random_atom;				// �d�ׂ������_���ɕ��z���錴�q�̐�
	int num_fixed_atom;				    // �����f�[�^�œd�ׂ��m�肳��Ă��錴�q�̐�
	string label[MAX_NUM];              // label for identifying particles
	int num_label;				        // number of independent label 
	string valid_label[MAX_NUM];        // list of independent label 
	long double t,t0,tn;                // time 
	int n_step;			                //���Ԃ̍��ݐ�
	long double interval_t;				//���Ԃ̍��݋� = (tn-t0)/n_step
	long double q[MAX_NUM][3][2];       //�ۂ߂̌덷
	long double px,py,pz;				//�n�S�̂̉^���ʂ̐����B 
	long double refmomentum;			//���x�]���Ŏg�p�B�ő�̑S�^���ʂ̘a�̐����B 
	long double minmomentum;			//���x�]���Ŏg�p�B�ŏ��̉^���ʐ����B�@
	long double refenergy;				// ���x�]���Ŏg�p
	int snap_num;						// snapshot is took every snap_num step

	ofstream foutinformation;			//���t�@�C�����@fimename[i]��label�ԍ�i�̃t�@�C�����@
    char InfoFile[32];					//���t�@�C���X�g���[���@
    char filename[MAX_LABEL*2][32];		//�o�̓t�@�C�����@fimename[i]��label�ԍ�i�̃t�@�C�����@
	ofstream outfile[MAX_LABEL*2];		//�o�̓t�@�C���X�g���[�� outfile[i]��label�ԍ�i�̃t�@�C��
    char filename_2body[MAX_LABEL][MAX_LABEL][32];		//�o�̓t�@�C�����@fimename[i]��label�ԍ�i�̃t�@�C�����@
	ofstream outfile_2body[MAX_LABEL][MAX_LABEL];		//�o�̓t�@�C���X�g���[�� outfile[i]��label�ԍ�i�̃t�@�C��
	char SnapFile[32];			        //�X�i�b�v�V���b�g�t�@�C�����@
    ofstream fsnapshot;					//�X�i�b�v�V���b�g�p�t�@�C���X�g���[���@

	char pot_flag;			//���ݍ�p�@default: C, C: Coulomb, L: Lenard-Jones, M: Coulomb+LJ, N: Coulomb + LJ + Induced dipole
	char pot_label[8];		//���ݍ�p���ʃ��x���@C: Coulomb, L: Lenard-Jones, M: Coulomb+LJ, N: Coulomb + LJ + Induced dipole

public:
	int trial;              // ���s��
	void keisan(int n_th);
	long double kyori(long double,long double,long double ,long double,long double,long double);
	int dice(int faces);
	void rkg();

	long double force(int i,int j); //���q�Ԃ̗͂��v�Z�B���q i �� dv/dt �̉E�ӂ�Ԃ�   j �͍��W������
	long double potential_energy(); //�S�|�e���V�����G�l���M�[���v�Z [amu nm^2/ps^2]
	long double kinetic_energy();   //�S�^���G�l���M�[���v�Z [amu nm^2/ps^2]
	long double PE_plus_KE();       //�S�G�l���M�[�iPE+KE�j���v�Z [amu nm^2/ps^2]
	void cal_momentum();    //�����_���ɓd�ו��z����
	void dist_charge();    //�����_���ɓd�ו��z����

    void mk_out_dir(char *makingdirectory);
    void open_output_files(char *makingdirectory);
    void close_output_files(char *makingdirectory);
	int get_data(char *input_file);
	void snap_output(long double snap_time);
	void message();

};




void denka::message(){

cout << "usage: program.exe  parameterfile " << endl; 
cout << " " << endl; 
cout << "�p�����[�^�t�@�C���̏��� " << endl; 
cout << "�P�s�ځ@�@�J�n����[ps]  �I������[ps]  MD�X�e�b�v��[����]  �X�i�b�v�V���b�g��" << endl; 
cout << "�@�@      �����̊Ԃ̓X�y�[�X�ŋ�؂�B" << endl; 
cout << "�@�@      �X�i�b�v�V���b�g��=0�Ȃ�΃X�i�b�v�V���b�g���o�͂��Ȃ��B" << endl; 
cout << " " << endl; 
cout << "�Q�s�ځ@�@�S���q���@�S�d����  ���ݍ�p�w��" << endl; 
cout << "�@�@�@�@�@�S���q�� > �S�d�����B" << endl; 
cout << "�@�@�@�@�@���ݍ�p: default=C: Coulomb,���݂�Coulomb�̂ݗL�� " << endl; 
cout << " " << endl; 
cout << "�R�s�ځ@�@�V�~�����[�V�����̎��s��" << endl; 
cout << " " << endl; 
cout << "�S�s�ڈȍ~�ɏ������q�f�[�^������ " << endl; 
cout << "���� " << endl; 
cout << "���x��  ����[amu]  x[nm] y[nm] z[nm] vx[nm/ps] vy[nm/ps] vz[nm/ps] �����d�� �ő�d��" << endl; 
cout << "���x���͂T�����ȓ��B����ɏ]���Č��ʂ𕪗ނ���B " << endl; 
cout << "�����d�ׂ�u���ăV�~�����[�V�����B�����d�ׂ͂O�ȏ�̐���" << endl; 
cout << "�ő�d�ׂ܂ŗ��q�Ƀ����_���ɓd�ׂ𕪔z����B " << endl; 
exit(0);
}

void denka::mk_out_dir(char *makingdirectory){

char NewDir[32];
int i;

//���s�t�@�C�����Ɠ����i�g���q�����j�f�B���N�g�������쐬�@�E�E.exe ���s�t�@�C���@�p�����^�t�@�C��
    i=0;
	do {
		NewDir[i]=makingdirectory[i];
        i=i+1;
        if(i>31) {cout<<"Filename is too long. Be fewer 32 letters."<<endl;exit (1);}
	} while(makingdirectory[i]!='.');
    NewDir[i]=NULL;
	cout<<"simulation result -> "<<NewDir<<endl;
    _mkdir(NewDir); //�f�B���N�g���쐬����


}
//�o�̓f�B���N�g���쐬

void denka::open_output_files(char *makingdirectory){
char DirName[32];
int i,j;

    i=0;
	do {
		DirName[i]=makingdirectory[i];
        i=i+1;
        if(i>31) {cout<<"Filename is too long. Be fewer 32 letters."<<endl;exit (1);}
	} while(makingdirectory[i]!='.');
    DirName[i]=NULL;


         //�V�~�����[�V���������̏o�̓t�@�C���̍쐬
		   strcpy_s(InfoFile,DirName);
		   strcat_s(InfoFile,"\\Infomation.txt");
			foutinformation.open(InfoFile);
				cout << InfoFile << endl;
			if(!foutinformation) {cout<<" cannot make Infomation FILE. \n";}
		  foutinformation<<"START TIME T0 "<<t0<<" ps"<<endl;
	      foutinformation<<"STOP TIME Tn "<<tn<<" ps"<<endl;
	      foutinformation<<"kizami n "<< n_step <<endl;
		  foutinformation<<"trial number "<<trial<<endl;
	      foutinformation<<"Number of simulated particles is "<<PN<<endl;
		  foutinformation<<"Number of charge is "<<CN<<endl;
		  foutinformation<<"Interaction : "<< pot_label <<endl;
		  foutinformation<<"New directory name is "<<DirName<<endl;
		  foutinformation<<"N"<<DirName<<endl;
		  foutinformation<<endl;
		  foutinformation<<endl;
		  foutinformation<<endl;


			//�t�@�C�����J���A"label"P.dat (ion) labelN.dat (neutral)  �ȂǂȂ�

			if(MAX_LABEL<num_label){cout << "too much label! " <<endl; }

			for(i=0;i<num_label;i++){
				strcpy_s(filename[2*i]  ,DirName);			// label i-th ion
				strcpy_s(filename[2*i+1],DirName);			// label i-th neutral 
				strcat_s(filename[2*i]  ,"\\");
				strcat_s(filename[2*i+1],"\\");
				strcat_s(filename[2*i]  ,valid_label[i].c_str());
				strcat_s(filename[2*i+1],valid_label[i].c_str());
				strcat_s(filename[2*i]  ,"_P.dat");
				strcat_s(filename[2*i+1],"_N.dat");
              //�t�@�C���쐬
				outfile[2*i].open(  filename[2*i]);
//				cout << filename[2*i] << endl;
                if(!outfile[2*i]) {cout << "cannot make \n ";}
				outfile[2*i+1].open(filename[2*i+1]);
//				cout << filename[2*i+1] << endl;
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
              //�t�@�C���쐬
				outfile_2body[i][j].open(  filename_2body[i][j]);
//				cout << filename[2*i] << endl;
                if(!outfile_2body[i][j]) {cout << "cannot make \n ";}
			}
			}

			if(snap_num>0){
				strcpy_s(SnapFile ,DirName);		
				strcat_s(SnapFile,"\\Snapshot.dat");
				fsnapshot.open(SnapFile);
				cout << SnapFile << endl;
				if(!fsnapshot) {cout<<" cannot make Infomation FILE. \n";}
			}
}

void denka::close_output_files(char *makingdirectory){
int i;

			for(i=0;i<num_label;i++){outfile[2*i].close();outfile[2*i+1].close();}
			foutinformation.close();
			if(snap_num>0){	fsnapshot.close();}

}

void denka::keisan(int n_th) {
	int i,j,k,temp_i;
	int suffix_1,suffix_2;
	long double firstP[3]={0,0,0};   //�����̑S�^���ʌv�Z
	long double firstenergy=0;
	long double finalenergy=0;
	double Vx1,Vy1,Vz1;
	double Vx2,Vy2,Vz2;
	double InnerProduct;
	long double ev[MAX_NUM]; //�e���q�̉^���G�l���M�[
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


//		  //�����̑S�G�l���M�[�v�Z
			cout<<"first potential energy  "<< (AMU*1e6/UNIT_E)*potential_energy() <<"  eV"<<endl;  //�Ō�� [amu nm^2/ps^2] --> [eV]�@�֕ϊ�
			cout<<"first kinetic energy  "<< (AMU*1e6/UNIT_E)*kinetic_energy() <<"  eV"<<endl;  //�Ō�� [amu nm^2/ps^2] --> [eV]�@�֕ϊ�
			firstenergy=PE_plus_KE();  ////�P�ʂ� [amu nm^2/ps^2]
			cout<<"first total energy  "<< (AMU*1e6/UNIT_E)*firstenergy<<"  eV"<<endl;  //�Ō�� [amu nm^2/ps^2] --> [eV]�@�֕ϊ�

			 
            //Runge Kutta Gill �Ōv�Z�@n �͎��Ԃ̍��݂̐��B���̐����� for �Ōv�Z
			for(i=0;i<=n_step;i++) {
				t=t0+i*interval_t;
                rkg();

				if( (snap_num!=0) && (i%snap_num==0) ){snap_output(t);}

           //�G�l���M�[�v�Z
   				for(k=0;k<PN;k++)  {
						ev[k]=0.5*(AMU*1e6/UNIT_E)*mass[k]*(velocity[k][0]*velocity[k][0] + velocity[k][1]*velocity[k][1]+velocity[k][2]*velocity[k][2]);
			  	};    // �G�l���M�[�P�� [amu nm^2/ps^2] --> [eV]�@�֕ϊ�
			}

			if( (snap_num!=0) && (n_step%snap_num!=0) ){snap_output(t);} //�Ō��step�ł�sanpshot���o�͂���B


			cal_momentum();
			cout<<"total momentum / [amu nm/ps]"<<endl<<
				"px="<<px<<" py="<<py<<" pz="<<pz<<endl;

//		  //�����̑S�G�l���M�[�v�Z
			cout<<"final potential energy  "<< (AMU*1e6/UNIT_E)*potential_energy() <<"  eV"<<endl;  //�Ō�� [amu nm^2/ps^2] --> [eV]�@�֕ϊ�
			cout<<"final kinetic energy  "<< (AMU*1e6/UNIT_E)*kinetic_energy() <<"  eV"<<endl;  //�Ō�� [amu nm^2/ps^2] --> [eV]�@�֕ϊ�
			finalenergy=PE_plus_KE();  ////�P�ʂ� [amu nm^2/ps^2]
			cout<<"final total energy  "<< (AMU*1e6/UNIT_E)*finalenergy<<"  eV"<<endl;  //�Ō�� [amu nm^2/ps^2] --> [eV]�@�֕ϊ�

			if(refenergy<(fabs((finalenergy-firstenergy)/finalenergy))) {
				refenergy = (fabs((finalenergy-firstenergy)/finalenergy));}

			if(refmomentum<fabs(kyori(px,firstP[0],py,firstP[1],pz,firstP[2])/minmomentum))  //�傫������refmomentum�ɍ̗p
			{refmomentum=fabs(kyori(px,firstP[0],py,firstP[1],pz,firstP[2])/minmomentum);}
       
           
		   cout<<" seido (first total momentum-final total momentum)/(minmomentum) "<<refmomentum<<endl;
		   cout<<" seido max(first-final/final)energy = "<<refenergy<<endl;
        

           foutinformation<<endl;
		   foutinformation<<"momentum seido "<<refmomentum<<endl;
		   foutinformation<<"energy seido   "<<refenergy<<endl;



		   //label�ɏ]���Č��ʂ��t�@�C���ɏo�͂���
			for(i=0;i<PN;i++)  {
               if(charge[i]==0) {
				   for(j=0;j<num_label;j++){if(label[i]==valid_label[j]){file_num=2*j+1;};}
//						cout << charge[i] << " " << label[i] << " " << filename[file_num] << endl;
			   }
               if(charge[i]>0) {
				   for(j=0;j<num_label;j++){if(label[i]==valid_label[j]){file_num=2*j;};}
//						cout << charge[i] << " " << label[i] << " " << filename[file_num] << endl;
			   }
			outfile[file_num] <<n_th<< " " << charge[i]<< " " << ev[i]<< " " << momentum[i][0] << " " << momentum[i][1] << " " << momentum[i][2]<< " "<< "0"<< endl;
			foutinformation  <<n_th<< " "  <<label[i]<< " " << charge[i]<< " " << ev[i]<< " " << momentum[i][0] << " " << momentum[i][1] << " " << momentum[i][2]<< " "<< "0"<< endl;


				if(charge[i]>0) {
					for(j=i+1;j<PN;j++){
						if(charge[j]>0) {
							for(k=0;k<num_label;k++){
								if(label[i]==valid_label[k]){suffix_1=k;};
								if(label[j]==valid_label[k]){suffix_2=k;};
							}
							if(suffix_1>suffix_2){temp_i= suffix_1; suffix_1=suffix_2; suffix_2=temp_i; }
							Vx1= momentum[i][0]; Vy1= momentum[i][1]; Vz1= momentum[i][2]; 
							Vx2= momentum[j][0]; Vy2= momentum[j][1]; Vz2= momentum[j][2]; 
							InnerProduct = (Vx1*Vx2 + Vy1*Vy2 + Vz1*Vz2)/sqrt((Vx1*Vx1 + Vy1*Vy1 + Vz1*Vz1)*(Vx2*Vx2 + Vy2*Vy2 + Vz2*Vz2));
							outfile_2body[suffix_1][suffix_2]  <<n_th<< " "  << label[i]  << " "  << charge[i]  << " "  << label[j]  << " "  << charge[j]  << " "  << InnerProduct  << " "  << ev[i]  << " " << ev[j]  << " " << endl;
						}
					}
				}

			}

			

}

//�^���ʌv�Z
void denka::cal_momentum(){
int i,j;

//�^���ʂ����߂�
for(i=0;i<PN;i++) {for(j=0;j<=2;j++)  { momentum[i][j]=mass[i]*velocity[i][j];};}

px=0;	py=0;	pz=0;

//�S�^���ʂ̘a
	for(i=0;i<PN;i++) {
		px=px+momentum[i][0];        py=py+momentum[i][1];		pz=pz+momentum[i][2];
		//�^���ʐ��x�]���@�^���ʃx�N�g���̐�Βl�ŕ]������ �悸�ł��������C�I���̉^���ʃx�N�g����T��
		if(i==0) {minmomentum=kyori(momentum[i][0],0,momentum[i][1],0,momentum[i][2],0);
		}else  if(i) {
			if(minmomentum>kyori(momentum[i][0],0,momentum[i][1],0,momentum[i][2],0))   
			  {minmomentum=kyori(momentum[i][0],0,momentum[i][1],0,momentum[i][2],0);}
		}
		}
}			//�S�^���ʂ̘a�I��


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

long double denka::kyori(long double a1,long double b1,long double a2,long double b2,long double a3,long double b3)  {
return sqrt((a1-b1)*(a1-b1)+(a2-b2)*(a2-b2)+(a3-b3)*(a3-b3));
}



//�Q�_�Ԃ̋������v�Z

void denka::rkg()
      
      {
       static long double p1=0.2928932188134524, p2=1.707106781186547;   //  p1 = 1 - sqrt(1/2),  p2 = 1 + sqrt(1/2) 
	   long double ck0[MAX_NUM][3][2];
	   long double q1[MAX_NUM][3][2];
	   long double r1[MAX_NUM][3][2];
	   long double y1[MAX_NUM][3][2];
       int i,j,k;
   
	   for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {ck0[i][j][k]=0; q1[i][j][k]=0;  r1[i][j][k]=0;y1[i][j][k]=0;};};}

      
	    //i ���q��\�� j ���W��\��
	   // �S�X�e�b�v�̌v�Z�J�n
	   //  step 1
	   for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {ck0[i][j][k]=force(i,j)*interval_t/2;}
	   else    {ck0[i][j][k]=velocity[i][j]*interval_t/2;}
	   ;};}  ;}
        
       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {y1[i][j][k]=velocity[i][j]+ck0[i][j][k]-q[i][j][k];}
	   else    {y1[i][j][k]=position[i][j]+ck0[i][j][k]-q[i][j][k];}
					  };};}


       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {r1[i][j][k]=y1[i][j][k]-velocity[i][j];}
	   else    {r1[i][j][k]=y1[i][j][k]-position[i][j];}
					  };};}

	  
	   for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {q[i][j][k]=q[i][j][k]+r1[i][j][k]*3-ck0[i][j][k];
	   };};}

       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {velocity[i][j]=y1[i][j][k];}
	   else    {position[i][j]=y1[i][j][k];}
					  };};}
       


       //  step 2
//       xx=t+interval_t*0.5;

       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {ck0[i][j][k]=force(i,j)*interval_t;}
	   else    {ck0[i][j][k]=velocity[i][j]*interval_t;}
	   };};}


        for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {y1[i][j][k]=velocity[i][j]+p1*(ck0[i][j][k]-q[i][j][k]);}
	   else    {y1[i][j][k]=position[i][j]+p1*(ck0[i][j][k]-q[i][j][k]);}
					  };};}


       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {r1[i][j][k]=y1[i][j][k]-velocity[i][j];}
	   else    {r1[i][j][k]=y1[i][j][k]-position[i][j];}
					  };};}

      
	    for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {q[i][j][k]=q[i][j][k]+r1[i][j][k]*3-ck0[i][j][k]*p1;
	   };};}
      
	   
	   for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {velocity[i][j]=y1[i][j][k];}
	   else    {position[i][j]=y1[i][j][k];}
					  };};}
    
	  

       //  step 3
       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {ck0[i][j][k]=force(i,j)*interval_t;}
	   else    {ck0[i][j][k]=velocity[i][j]*interval_t;}
	   };};}       

        for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {y1[i][j][k]=velocity[i][j]+p2*(ck0[i][j][k]-q[i][j][k]);}
	   else    {y1[i][j][k]=position[i][j]+p2*(ck0[i][j][k]-q[i][j][k]);}
					  };};}


       
	   for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {r1[i][j][k]=y1[i][j][k]-velocity[i][j];}
	   else    {r1[i][j][k]=y1[i][j][k]-position[i][j];}
					  };};}

       
        for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {q[i][j][k]=q[i][j][k]+r1[i][j][k]*3-ck0[i][j][k]*p2;
	   };};}


	   for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {velocity[i][j]=y1[i][j][k];}
	   else    {position[i][j]=y1[i][j][k];}
					  };};}

	  
       
       //  step 4
//     xx=t+interval_t/2;
      
       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {ck0[i][j][k]=force(i,j)*interval_t/2;}
	   else    {ck0[i][j][k]=velocity[i][j]*interval_t/2;}
	   };};}       
    

       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++) 
	   {if(k)  {y1[i][j][k]=velocity[i][j]+(ck0[i][j][k]-q[i][j][k])/3;}
	   else    {y1[i][j][k]=position[i][j]+(ck0[i][j][k]-q[i][j][k])/3;}
	   };};}      


       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++)
	   {if(k)  {r1[i][j][k]=y1[i][j][k]-velocity[i][j];}
	   else    {r1[i][j][k]=y1[i][j][k]-position[i][j];};};};}


       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++)
	   {q[i][j][k]=q[i][j][k]+r1[i][j][k]*3-ck0[i][j][k];};};}
	   

       for(i=0;i<=PN-1;i++) {for(j=0;j<=2;j++) {for(k=0;k<=1;k++)
	   {if(k)  {velocity[i][j]=y1[i][j][k];}
	   else    {position[i][j]=y1[i][j][k];};};};}

          
    
 	  }



//Runge Kutta Gill calculation
int denka::dice(int faces) { 	return (1 + rand()%faces);}  
// n-face dice  �Ԃ�l1�`n

void denka::dist_charge() {
int i,num,max_charge,fixed_charge;         

     srand( (unsigned)time(NULL) );

	max_charge=0; fixed_charge=0;
	for(i=0;i<PN;i++){max_charge=max_charge+Q_limit[i];
			fixed_charge=fixed_charge+charge[i];
		}

	if(CN > max_charge){ cout << " charge number mismatch" << endl; exit(0);}

	for(i=0;i<(CN - fixed_charge);i++){
			num = dice(PN)-1;
	        while( charge[num]>=Q_limit[num]){num = dice(PN)-1; }
				charge[num]=charge[num]+1;
		}						//  �d��(CN)�������_�����z�B

}


	  

//�͂̌v�Z   ���q i �� dv/dt �̉E�ӂ�Ԃ�   j �͍��W������  �����̃G�l���M�[�P�ʂ� [amu nm^2/ps^2]�@�͂� [amu nm/ps^2]
long double denka::force(int i ,int j)
{
  
   int seisu=0;
   long double nagasa=0;
   long double uhen=0;
	     
		   switch(pot_flag){
				case 'C'	: 
				   for (seisu=0;seisu<PN;seisu++)  {
						nagasa=kyori(position[i][0],position[seisu][0],position[i][1],position[seisu][1],position[i][2],position[seisu][2]);
						if(seisu-i) {uhen=uhen+(UNIT_E*UNIT_E/(4*PI*E_0)*1e3)*charge[i]*charge[seisu]*(position[i][j]-position[seisu][j])/pow(nagasa,3) /(AMU*mass[i]);}
				   };  //�ꉿ�C�I���Ԃ̃N�[�������ݍ�p
					break;
				case 'L'	:
				   for (seisu=0;seisu<PN;seisu++)  {
						nagasa=kyori(position[i][0],position[seisu][0],position[i][1],position[seisu][1],position[i][2],position[seisu][2]);
						   if(seisu-i) {
							uhen=uhen+(4*(EPSILON * K_B *1e-6)*(-6*pow(SIGMA/nagasa,6)+12*pow(SIGMA/nagasa,12)))*(position[i][j]-position[seisu][j])/pow(nagasa,2) /(AMU*mass[i]);
					   };
					};  //�������q�Ԃ�LJ���ݍ�p�@�C�I���ԃN�[�����͖���
					break;
				case 'M' : 
				   for (seisu=0;seisu<PN;seisu++)  {
						nagasa=kyori(position[i][0],position[seisu][0],position[i][1],position[seisu][1],position[i][2],position[seisu][2]);
						   if(seisu-i) {
							uhen=uhen+((UNIT_E*UNIT_E/(4*PI*E_0)*1e3)*charge[i]*charge[seisu]/nagasa + 4*(EPSILON * K_B *1e-6)*(-6*pow(SIGMA/nagasa,6)+12*pow(SIGMA/nagasa,12)))
						   		 *(position[i][j]-position[seisu][j])/pow(nagasa,2) /(AMU*mass[i]);
					   };
					};  //�N�[���� + LJ���ݍ�p
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
				   };  //�N�[���� + LJ + �C�I���|�U�N�o�Ɏq���ݍ�p
					break;
				default	: 
				   for (seisu=0;seisu<PN;seisu++)  {
						nagasa=kyori(position[i][0],position[seisu][0],position[i][1],position[seisu][1],position[i][2],position[seisu][2]);
						   if(seisu-i) {
							uhen=uhen+(UNIT_E*UNIT_E/(4*PI*E_0)*1e3)*charge[i]*charge[seisu]*(position[i][j]-position[seisu][j])/pow(nagasa,3) /(AMU*mass[i]);
						   };
				   };  //�f�t�H���g�͈ꉿ�C�I���Ԃ̃N�[�������ݍ�p�̂݁B
					break;
			}
    return uhen;
}


//�S�G�l���M�[�v�Z
long double denka::PE_plus_KE(){

	return potential_energy() + kinetic_energy();

} //�S�G�l���M�[�iPE+KE�j���v�Z [amu nm^2/ps^2]
//�S�|�e���V�����G�l���M�[�v�Z
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
		};  //�ꉿ�C�I���Ԃ̃N�[�������ݍ�p�B�G�l���M�[�P�ʂ� [amu nm^2/ps^2]
		break;
	case 'L'	: 
		for(i=0;i<PN;i++) {
			for(j=i+1;j<PN;j++){
				nagasa=kyori(position[i][0],position[j][0],position[i][1],position[j][1],position[i][2],position[j][2]);
				t_energy=t_energy+4* (EPSILON * K_B *1e-6/AMU)*(pow(SIGMA/nagasa,12)-pow(SIGMA/nagasa,6)) ;
			};	
		};  //�������q�Ԃ�LJ���ݍ�p�@�C�I���ԃN�[�����͖����B�G�l���M�[�P�ʂ� [amu nm^2/ps^2]
		break;
	case 'M'	: 
		for(i=0;i<PN;i++) {
			for(j=i+1;j<PN;j++){
				nagasa=kyori(position[i][0],position[j][0],position[i][1],position[j][1],position[i][2],position[j][2]);
				t_energy=t_energy+(UNIT_E*UNIT_E/(4*PI*E_0*AMU)*1e3)*charge[i]*charge[j]/nagasa
								 +4*(EPSILON * K_B *1e-6/AMU)*(pow(SIGMA/nagasa,12)-pow(SIGMA/nagasa,6)) ;
			};	
		};  ////�N�[���� + LJ���ݍ�p�B�G�l���M�[�P�ʂ� [amu nm^2/ps^2]
		break;
	case 'N'	: 
		for(i=0;i<PN;i++) {
			for(j=i+1;j<PN;j++){
				nagasa=kyori(position[i][0],position[j][0],position[i][1],position[j][1],position[i][2],position[j][2]);
				t_energy=t_energy+(UNIT_E*UNIT_E/(4*PI*E_0*AMU)*1e3)*charge[i]*charge[j]/nagasa
								 +4*(EPSILON * K_B *1e-6/AMU)*(pow(SIGMA/nagasa,12)-pow(SIGMA/nagasa,6)) 
								 - (ALPHA_CGS * UNIT_E * UNIT_E * 1e30/(8*PI*E_0*AMU) )*charge[i]/pow(nagasa,4) ;
			};	
		};  ////�N�[���� + LJ���ݍ�p+ �C�I���|�U�N�o�Ɏq���ݍ�p�B�G�l���M�[�P�ʂ� [amu nm^2/ps^2]
		break;
	default	: 
		for(i=0;i<PN;i++) {
			for(j=i+1;j<PN;j++){
				nagasa=kyori(position[i][0],position[j][0],position[i][1],position[j][1],position[i][2],position[j][2]);
				t_energy=t_energy+(UNIT_E*UNIT_E/(4*PI*E_0*AMU)*1e3)*charge[i]*charge[j]/nagasa;
			};	
		};  //�ꉿ�C�I���Ԃ̃N�[�������ݍ�p�B�G�l���M�[�P�ʂ� [amu nm^2/ps^2]
		break;
}
return t_energy;
} //�S���q��potential�G�l���M�[���v�Z
					
//�S�^���G�l���M�[�v�Z
long double denka::kinetic_energy(){
int i;
long double t_energy=0;

for(i=0;i<PN;i++) {	t_energy=t_energy+0.5*mass[i]*(velocity[i][0]*velocity[i][0]+velocity[i][1]*velocity[i][1] + velocity[i][2]*velocity[i][2]);}
return t_energy;
} //�S���q�̉^���G�l���M�[���v�Z


int denka::get_data(char *input_file){

  char buf[512];
  int  i,j,k,LINE_NUM,flag;

 ifstream fin(input_file);

 if(!fin) {cout << "cannot open \n ";
          return 1;}


 LINE_NUM=0;

    while( fin.getline( buf, sizeof(buf) )){
		LINE_NUM=LINE_NUM+1;
//		cout << buf << endl;
		istringstream is( buf );
		if(LINE_NUM==1){ is >> t0 >> tn >> n_step>> snap_num ;}
		if(LINE_NUM==2){ is >> PN >> CN >> pot_flag ;}
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

for(j=0;j<3;j++) {for(i=0;i<PN;i++) {momentum[i][j]=0;} ;}//������
for(i=0;i<PN;i++) {for(j=0;j<3;j++)  {for(k=0;k<2;k++) {q[i][j][k]=0;};};}   //�ۂ߂̌덷������

valid_label[0] = label[0];
num_label=1;
for(i=1;i<PN;i++){
	flag=0;
	for(j=0;j<num_label;j++){ if(label[i] == valid_label[j]){ flag=1;} ;}
	if(flag==0){valid_label[num_label] = label[i]; num_label=num_label+1 ;}
}			//  label�̐��ƁA�d�����Ȃ����x�����`�F�b�N

cout << LINE_NUM-3 << " particle's data. " << PN << " particles are valid." << endl;


num_random_atom=0;
num_fixed_atom=0;
for(i=0;i<PN;i++){
	if(charge[i]==-1){random_charge[num_random_atom]=i;num_random_atom=num_random_atom+1;};
	if(charge[i]!=-1){num_fixed_atom=num_fixed_atom+1;};
	}

//cout << num_random_atom << " ���q�ɓd�ׂ𕪔z����B " << endl;
//cout << num_fixed_atom << " ���q�͓d�ׂ��m��B " << endl;
//cout << CN << " �d�א��B " << endl;

//for(i=0;i<num_label;i++){ 	cout << label[i] << " " << valid_label[i] << endl;}
//cout << num_label << " label���B " << endl;


return 0;
}


//argv[1] ���f�[�^�t�@�C�� ���������@�ŏI�����@���ԍ��ݐ��@���q���@���ʁ@�d�ׁ@�����ʒuxyz�@�����xxyz �̏��@
//argy[2] ���p�����[�^�t�@�C���@�C�I�����@�m�茴�q�ԍ��@0 �����_�����q�ԍ� 0

//int _tmain(int argc, _TCHAR* argv[])
//{
//	return 0;
//}
//
//
int main(int argc,char *argv[])
{
int i;


cout<<"********************************PROGRAM START********************************"<<endl;
cout<<"length:nm,time:ps,energy:eV"<<endl;    

  denka sim_data;				//class ����
  if(argv[1]==NULL){sim_data.message();}

  sim_data.mk_out_dir(argv[1]); //�o�̓f�B���N�g���쐬

	sim_data.get_data(argv[1]);	//�����f�[�^�擾 �t�@�C���쐬�݂̂Ɏg�p����
	sim_data.open_output_files(argv[1]);//�o�̓t�@�C���쐬

for(i=0;i<sim_data.trial;i++) {
	sim_data.get_data(argv[1]);	//�����f�[�^�擾

	cout<<"******"<<i+1<<"th simulation"<< "****************************************************"<<endl;
	sim_data.dist_charge();       //�d�ׂ������_���ɕ��z�B
	sim_data.keisan(i+1);       //��step�̃V�~�����[�V������Runge-Kutta-Gill�ŉ���


 }
	sim_data.close_output_files(argv[1]);


}

