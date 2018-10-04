#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <time.h>
#define N 16 //numero de dados na rede é 4*(N^3)
#define J1 0.034 //J1 em eV/kb
#define J2 -0.156 //J2 em eV/kb
#define S 3.5   // spin do EuTe
#define Bi 8.5  //Campo inicial T
#define Bf 8.5  //Campo final T
#define Ti 2.01  //temperatura inicial K
#define Tf 22.01  //temperatura final K
#define NPTS 21  //numero de pontos entre tf e ti ou Bi e Bf
#define PDEV 20 //número de dados para a estatística
#define REP 15 //número de vezes que repete a estatística
#define MCI 50  //numero de iteracoes de monte carlo de cada levantamento de dados
#define Jxf 0.095
#define KB 8.61733e-5 //constante de boltzmann eV/K
#define muB 5.7883818e-5//magneton de bohr eV/T
#define gs 2.0023192 //fator gyromagnetico
#define PI 3.14159265
#define mef 0.3*9.1e-31 //mass efetiva
#define eps 8.3*8.85e-12 //condutividade
#define hbar 1.05e-34 //constnte de plank
#define Joule 6.24e18 //conversao Joule->eV
#define q 1.6e-19 //carga do eletron
#define ar 6.6e-10


void inicializa2(float M[N][N][N][3], int n);
long rodarM(float M[N][N][N][3], float T, float BB, int n, long idum);
double eff(float M[N][N][N][3], float V[3], float BB, int n, int x, int y, int z);
double exf(float M[N][N][N][3], float V[3], int x, int y, int z);
double ezeeman(float M[N][N][N][3], float V[3], float BB, int n, int x, int y, int z);
double cinetica();
double coulomb();
float magnetizacao(float M[N][N][N][3], double Mag[4], float BB, int n);
void MonteCarlo(float T, float BB, double media[10]);
void Dist(int i,int j,int k,int n,float R[3]);
float ran1(long *idum);
void RaioVar(float B);
void RaioNum(float B);

float M0[N][N][N][3],M1[N][N][N][3],M2[N][N][N][3],M3[N][N][N][3];
double U, Rb, RQ;

float B[3];
FILE *ofp;

int main(){
   int i,j,k;
   float T,BB;
   double bsf, Ap, Rp;
   double media[10];
   char ofname[15];

   srand(time(NULL));
   
   /*####### Raio de Bohr #######*/
   //Rb=1.5;
   //RQ=9.00;

   printf("digite o nome do arquivo de saida:\n");
   scanf("%s", &ofname);
   ofp=fopen(ofname, "w");
   fprintf(ofp,"# N: %d, dados para estatistica: %d, iteracoes: %d Jfx: %f J2: %f Ti: %f Bi: %f\n",N,PDEV,MCI,Jxf,J2,Ti,Bi);
   fprintf(ofp,"B T Rb Eff NDEff Exf NDExf Ez NDEz Ek Ev\n");
   printf("concluido:\n  0%c\n",'%');

   for(i=0;i<NPTS;i++){
      T=Ti+i*(Tf-Ti)/(NPTS-1);
      BB=Bi+i*(Bf-Bi)/(NPTS-1);
      
      //RaioNum(BB);
      RaioVar(BB);

      for(k=0;k<8;k++)
         media[k]=0;
      
      for(j=0;j<REP;j++){
         MonteCarlo(T,BB,media);
      }
      for(k=0;k<6;k++)
         media[k]=media[k]/REP;
      media[6]=cinetica();
      media[7]=coulomb();
      fprintf(ofp,"%f %f %f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",BB,T, Rb,media[0],media[1]*2.355,media[2],media[3]*2.355,media[4],media[5]*2.355,media[6],media[7]);
      printf("\033[F\033[J%2d%c\n",100*i/(NPTS-1),'%');
   }

   printf("\n");
   fclose(ofp);
   return 0;
}

void MonteCarlo(float T, float BB, double media[10]){
   int mci, n, i;
   long int idum;
   double Mag[5], d;
   double sum0=0, sumquad0=0, sum1=0, sumquad1=0, sum2=0, sumquad2=0, sum3=0, sumquad3=0, sum4=0, sumquad4=0;
   d=sqrt(3);
   B[0]=1.0/d;B[1]=1.0/d;B[2]=1.0/d;
   
   for(n=0;n<PDEV;n++){
       inicializa2(M0,0);
       inicializa2(M1,1);
       inicializa2(M2,2);
       inicializa2(M3,3);
       for(mci=MCI;mci>0;mci--){ //numero de iteracoes de monte carlo
          idum=rodarM(M0,T,BB,0,idum);
          idum=rodarM(M1,T,BB,1,idum);
          idum=rodarM(M2,T,BB,2,idum);
          idum=rodarM(M3,T,BB,3,idum);
       }
       for(i=0;i<4;i++)
          Mag[i]=0;
       magnetizacao(M0,Mag,BB,0); //calcula a energia
       magnetizacao(M1,Mag,BB,1);
       magnetizacao(M2,Mag,BB,2);
       magnetizacao(M3,Mag,BB,3);
       
       sum0=sum0+Mag[0];
       sumquad0=sumquad0+Mag[0]*Mag[0];
       sum1=sum1+Mag[1];
       sumquad1=sumquad1+Mag[1]*Mag[1];
       sum2=sum2+Mag[2];
       sumquad2=sumquad2+Mag[2]*Mag[2];
   }
   media[0]=media[0]+sum0/PDEV;
   media[1]=media[1]+sqrt((sumquad0-sum0*sum0/PDEV)/(PDEV-1));
   media[2]=media[2]+sum1/PDEV;
   media[3]=media[3]+sqrt((sumquad1-sum1*sum1/PDEV)/(PDEV-1));
   media[4]=media[4]+sum2/PDEV;
   media[5]=media[5]+sqrt((sumquad2-sum2*sum2/PDEV)/(PDEV-1));
}

long rodarM(float M[N][N][N][3], float T, float BB, int n, long idum){
   int i,j,k,l;
   float r,V[3];
   double P,H;
   double theta, phi, st;
   for(i=0;i<N;i++) //cada um dos spins da rede
      for(j=0;j<N;j++)
         for(k=0;k<N;k++){ //troca M[i][j][k][] por um vetor V[] de modulo 1
            theta=PI*ran1(&idum);
            phi=2.0*PI*ran1(&idum);
            V[2]=S*cos(theta);
            st=sin(theta);
            V[0]=S*st*cos(phi);
            V[1]=S*st*sin(phi);
            H=2*eff(M,V,BB,n,i,j,k)+ezeeman(M,V,BB,n,i,j,k)+exf(M,V,i,j,k); //calcula a dif de energia entre o antigo e novo estado
            if(H<=0) //realiza a troca
               for(l=0;l<3;l++)
                  M[i][j][k][l]=V[l];
            else{
               P=expl(-H/T); //calcula a probabilidade de troca
               r=ran1(&idum);
               if(r<P) //realiza a troca
                  for(l=0;l<3;l++)
                     M[i][j][k][l]=V[l];
            }
         }
   return idum;
}

double eff(float M[N][N][N][3], float V[3], float BB, int n, int x, int y, int z){
   int i,j,k,l;
   float NN=0, NN2=0, NNN=0, NNN2=0, H=0, H2=0, ri;
   //calcula energia de exchange de NN
/**/
   if(n==0)
      for(i=-1;i<1;i++)
         for(j=-1;j<1;j++)
            for(k=0;k<3;k++){
               NN=NN+M1[(x+N+i)%N][(y+N+j)%N][z][k]*M0[x][y][z][k]; //xy
               NN=NN+M2[(x+N+i)%N][y][(z+N+j)%N][k]*M0[x][y][z][k]; //xz
               NN=NN+M3[x][(y+N+i)%N][(z+N+j)%N][k]*M0[x][y][z][k]; //yz
               NN2=NN2+M1[(x+N+i)%N][(y+N+j)%N][z][k]*V[k]; //xy
               NN2=NN2+M2[(x+N+i)%N][y][(z+N+j)%N][k]*V[k]; //xz
               NN2=NN2+M3[x][(y+N+i)%N][(z+N+j)%N][k]*V[k]; //yz
            }
   if(n==1)
      for(i=-1;i<1;i++)
         for(j=-1;j<1;j++)
            for(k=0;k<3;k++){
               NN=NN+M0[(x+N+i+1)%N][(y+N+j+1)%N][z][k]*M1[x][y][z][k]; //xy
               NN=NN+M3[(x+N+i+1)%N][y][(z+N+j)%N][k]*M1[x][y][z][k]; //xz
               NN=NN+M2[x][(y+N+i+1)%N][(z+N+j)%N][k]*M1[x][y][z][k]; //yz
               NN2=NN2+M0[(x+N+i+1)%N][(y+N+j+1)%N][z][k]*V[k]; //xy
               NN2=NN2+M3[(x+N+i+1)%N][y][(z+N+j)%N][k]*V[k]; //xz
               NN2=NN2+M2[x][(y+N+i+1)%N][(z+N+j)%N][k]*V[k]; //yz
            }
   if(n==2)
      for(i=-1;i<1;i++)
         for(j=-1;j<1;j++)
            for(k=0;k<3;k++){
               NN=NN+M3[(x+N+i+1)%N][(y+N+j)%N][z][k]*M2[x][y][z][k]; //xy
               NN=NN+M0[(x+N+i+1)%N][y][(z+N+j+1)%N][k]*M2[x][y][z][k]; //xz
               NN=NN+M1[x][(y+N+i)%N][(z+N+j+1)%N][k]*M2[x][y][z][k]; //yz
               NN2=NN2+M3[(x+N+i+1)%N][(y+N+j)%N][z][k]*V[k]; //xy
               NN2=NN2+M0[(x+N+i+1)%N][y][(z+N+j+1)%N][k]*V[k]; //xz
               NN2=NN2+M1[x][(y+N+i)%N][(z+N+j+1)%N][k]*V[k]; //yz
            }
   if(n==3)
      for(i=-1;i<1;i++)
         for(j=-1;j<1;j++)
            for(k=0;k<3;k++){
               NN=NN+M2[(x+N+i)%N][(y+N+j+1)%N][z][k]*M3[x][y][z][k]; //xy
               NN=NN+M1[(x+N+i)%N][y][(z+N+j+1)%N][k]*M3[x][y][z][k]; //xz
               NN=NN+M0[x][(y+N+i+1)%N][(z+N+j+1)%N][k]*M3[x][y][z][k]; //yz
               NN2=NN2+M2[(x+N+i)%N][(y+N+j+1)%N][z][k]*V[k]; //xy
               NN2=NN2+M1[(x+N+i)%N][y][(z+N+j+1)%N][k]*V[k]; //xz
               NN2=NN2+M0[x][(y+N+i+1)%N][(z+N+j+1)%N][k]*V[k]; //yz
            }
/**/
   //calcula energia de exchange de NNN
/**/
   for(i=-1;i<2;i=i+2){
      for(k=0;k<3;k++){
         NNN=NNN+M[(x+N+i)%N][y][z][k]*M[x][y][z][k];
         NNN=NNN+M[x][(y+N+i)%N][z][k]*M[x][y][z][k];
         NNN=NNN+M[x][y][(z+N+i)%N][k]*M[x][y][z][k];
         NNN2=NNN2+M[(x+N+i)%N][y][z][k]*V[k];
         NNN2=NNN2+M[x][(y+N+i)%N][z][k]*V[k];
         NNN2=NNN2+M[x][y][(z+N+i)%N][k]*V[k];
      }
   }
/**/
   //calcula a energia final
   H=(-J1*NN-J2*NNN);
   H2=(-J1*NN2-J2*NNN2);
   return H2-H; //retorna a diferença de energia
}

double ezeeman(float M[N][N][N][3], float V[3], float BB, int n, int x, int y, int z){
   float H=0, H2=0, ri;
   //energia zeeman devido ao campo externo
   H=-gs*muB/KB*BB*(M[x][y][z][0]*B[0]+M[x][y][z][1]*B[1]+M[x][y][z][2]*B[2]);
   H2=-gs*muB/KB*BB*(V[0]*B[0]+V[1]*B[1]+V[2]*B[2]);

   return H2-H; //retorna a diferença de energia
}

double exf(float M[N][N][N][3], float V[3], int x, int y, int z){
   float H=0, H2=0, ri;
   //energia zeeman devido ao campo do elétron
   return 0;
   ri=(x-N/2)*(x-N/2)+(y-N/2)*(y-N/2)+(z-N/2)*(z-N/2);
   ri=sqrt(ri);
   H=-Jxf/(4*PI*Rb*Rb*Rb)*(M[x][y][z][0]*B[0]+M[x][y][z][1]*B[1]+M[x][y][z][2]*B[2])*exp(-2*ri/Rb)/KB;
   H2=-Jxf/(4*PI*Rb*Rb*Rb)*(V[0]*B[0]+V[1]*B[1]+V[2]*B[2])*exp(-2*ri/Rb)/KB;
   return H2-H; //retorna a diferença de energia
}

double cinetica(){
   float ek;
   ek=(hbar*hbar*Joule)/(2*mef*Rb*ar*Rb*ar);
   return ek;
}

double coulomb(){
   float ev;
   ev=-(q*q*Joule)/(4*PI*eps*Rb*ar);
   return ev;
}

void inicializa2(float M[N][N][N][3], int n){//inicializa com uma determinada orientação
   int i,j,k,l;
   float d;
   d=sqrt(6);
   for(i=0;i<N;i++)
      for(j=0;j<N;j++)
         for(k=0;k<N;k++){ //cada um dos spins da rede
            if(((i+j+k)%2==0 && n!=0) || ((i+j+k)%2==1 && n==0)){
               M[i][j][k][0]=1.0/d*3.5;
               M[i][j][k][1]=1.0/d*3.5;
               M[i][j][k][2]=-2.0/d*3.5;
            }
            else{
               M[i][j][k][0]=-1.0/d*3.5;
               M[i][j][k][1]=-1.0/d*3.5;
               M[i][j][k][2]=2.0/d*3.5;
            }
         }
}

float magnetizacao(float M[N][N][N][3], double Mag[5], float BB, int n){
   int i,j,k,l,norm;
   float V[3];
   double a=0, b=0, c=0, d=0, e=0, ri;
   /**energia**/
   V[0]=0;
   V[1]=0;
   V[2]=0;
   l=0;
   for(i=0;i<N;i++)
      for(j=0;j<N;j++)
         for(k=0;k<N;k++) //cada um dos spins da rede
            {//if(((i*i-i*N+N*N/4)+(j*j-j*N+N*N/4)+(k*k-k*N+N*N/4))<RQ){
               a=a+eff(M,V,BB,n,i,j,k)*KB;
               b=b+exf(M,V,i,j,k)*KB;
               c=c+ezeeman(M,V,BB,n,i,j,k)*KB;
               l=l++;
            }
   Mag[0]=Mag[0]-a;
   Mag[1]=Mag[1]-b;
   Mag[2]=Mag[2]-c;
   /**/
}

void Dist(int i,int j,int k,int n,float R[3]){
   if(n==0){
      R[0]=i;R[1]=j;R[2]=k;
   }
   if(n==1){
      R[0]=(i+0.5);R[1]=(j+0.5);R[2]=k;
   }
   if(n==2){
      R[0]=(i+0.5);R[1]=j;R[2]=(k+0.5);
   }
   if(n==3){
      R[0]=i;R[1]=(j+0.5);R[2]=(k+0.5);
   }
}
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum){
   int j, r;
   long k;
   static long iy=0;
   static long iv[NTAB];
   float temp;
   
   /**/
   if(*idum <= 0 || !iy){
      if(-(*idum)<1) *idum=1;
      else *idum=-(*idum);
      for (j=NTAB+7;j>=0;j--){
         k=(*idum)/IQ;
         *idum=IA*(*idum-k*IQ)-IR*k;
         if(*idum<0) *idum+=IM;
         if(j<NTAB) iv[j] = *idum;
      }
      iy=iv[0];
   }
   k=(*idum)/IQ;
   *idum=IA*(*idum-k*IQ)-IR*k;
   if (*idum < 0) *idum += IM;
   j=iy/NDIV;
   iy=iv[j];
   iv[j]= *idum;
   if((temp=AM*iy) > RNMX) return RNMX;
   /**/
   //temp=1.0*rand()/(RAND_MAX+1.0);
   return temp;
}

void RaioVar(float B){
   if(B>=0)  Rb=1.369;
   if(B>0.3) Rb=1.373;
   if(B>0.4) Rb=1.378;
   if(B>0.5) Rb=1.382;
   if(B>0.6) Rb=1.387;
   if(B>0.7) Rb=1.391;
   if(B>0.8) Rb=1.396;
   if(B>0.9) Rb=1.400;
   if(B>1.0) Rb=1.405;
   if(B>1.1) Rb=1.409;
   if(B>1.2) Rb=1.414;
   if(B>1.3) Rb=1.418;
   if(B>1.4) Rb=1.423;
   if(B>1.5) Rb=1.427;
   if(B>1.6) Rb=1.432;
   if(B>1.7) Rb=1.441;
   if(B>1.8) Rb=1.445;
   if(B>1.9) Rb=1.450;
   if(B>2.0) Rb=1.454;
   if(B>2.1) Rb=1.459;
   if(B>2.2) Rb=1.463;
   if(B>2.3) Rb=1.468;
   if(B>2.4) Rb=1.477;
   if(B>2.5) Rb=1.481;
   if(B>2.6) Rb=1.486;
   if(B>2.7) Rb=1.490;
   if(B>2.8) Rb=1.499;
   if(B>2.9) Rb=1.504;
   if(B>3.0) Rb=1.508;
   if(B>3.1) Rb=1.517;
   if(B>3.2) Rb=1.522;
   if(B>3.3) Rb=1.531;
   if(B>3.4) Rb=1.535;
   if(B>3.5) Rb=1.544;
   if(B>3.6) Rb=1.549;
   if(B>3.7) Rb=1.553;
   if(B>3.8) Rb=1.562;
   if(B>3.9) Rb=1.571;
   if(B>4.0) Rb=1.576;
   if(B>4.1) Rb=1.585;
   if(B>4.2) Rb=1.589;
   if(B>4.3) Rb=1.598;
   if(B>4.4) Rb=1.607;
   if(B>4.5) Rb=1.616;
   if(B>4.6) Rb=1.621;
   if(B>4.7) Rb=1.630;
   if(B>4.8) Rb=1.639;
   if(B>4.9) Rb=1.648;
   if(B>5.0) Rb=1.657;
   if(B>5.1) Rb=1.666;
   if(B>5.2) Rb=1.675;
   if(B>5.3) Rb=1.688;
   if(B>5.4) Rb=1.697;
   if(B>5.5) Rb=1.706;
   if(B>5.6) Rb=1.720;
   if(B>5.7) Rb=1.733;
   if(B>5.8) Rb=1.747;
   if(B>5.9) Rb=1.760;
   if(B>6.0) Rb=1.783;
   if(B>6.1) Rb=1.805;
   if(B>6.2) Rb=1.850;
   if(B>6.3) Rb=1.954;
   if(B>6.4) Rb=2.206;
}
void RaioNum(float B){
   if(B>=0.00) Rb=1.675;
   if(B>0.634) Rb=1.678;
   if(B>1.268) Rb=1.683;
   if(B>1.902) Rb=1.692;
   if(B>2.536) Rb=1.705;
   if(B>3.170) Rb=1.725;
   if(B>3.804) Rb=1.754;
   if(B>4.438) Rb=1.796;
   if(B>5.072) Rb=1.862;
   if(B>5.706) Rb=1.973;
   if(B>6.340) Rb=2.211;
}
