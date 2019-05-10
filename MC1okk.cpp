/*排除了头对头构型*/
/*改变了温度更新函数*/
#include<iostream>
#include<cmath>
#include<fstream>
#include <stdlib.h>
#include<time.h>
using namespace std;
int main()
{
double min(double x,double y);
const int L=200,N=2000,M=1000;
int A[L][L],A1[L][L];
   for(int i=0;i<L;i++)
   {
        for(int j=0;j<=L;j++)
            {
                A[i][j]=0;
                A1[i][j]=0;
            }
   };

int B[6][2];
    B[0][0]=-1,B[0][1]=0,B[1][0]=-1,B[1][1]=1,B[2][0]=0,B[2][1]=1,
    B[3][0]=1,B[3][1]=0,B[4][0]=1,B[4][1]=-1,B[5][0]=0,B[5][1]=-1;

/*起始分子1位置*/
int P1[N][2],D1[N],a,b,c,i=0;/*P1：分子1位置，D1：分子方向，D11：与分子手性有关*/
int fk,fk1,fkk,fk11;
int A2[3*N][3];
   for(int i=0;i<3*N;i++)
   {
        for(int j=0;j<=3;j++)
            {
                A2[i][j]=0;
            }
   };


srand((int)time(0));
while(i<N)
{
    i=i+1;
    a=rand()%L;
    b=rand()%L;
    c=rand()%6;
    fk=(a-B[c][0]+L)%L;
    fk1=(b-B[c][1])%L;
    if(A[a][b]==0&&A[fk][fk1]==0)
    {
        A[a][b]=i,A[fk][fk1]=-i;
        P1[i-1][0]=a,P1[i-1][1]=b;
        D1[i-1]=c;
        A2[1+3*(i-1)+B[D1[i-1]][0]][1+B[D1[i-1]][1]]=1;
        fk=1+3*(i-1)+B[(D1[i-1]+1)%6][0],fk1=1+B[(D1[i-1]+1)%6][1];
        A2[fk][fk1]=2;
        fk=1+3*(i-1)+B[(D1[i-1]+5)%6][0],fk1=1+B[(D1[i-1]+5)%6][1];
        A2[fk][fk1]=3;
    }
    else
    {
        i=i-1;
    };
};

/*for(int i=0;i<3*N;i++)
{
    for(int j=0;j<3;j++)
    {
        cout<<A2[i][j]<<" ";
    }
    cout<<"\n";
}*/


/*起始分子2位置*/
int P2[M][2],D2[M][2],m,n=0;
int d,e,f,g;
i=0;
while(i<M)
{
    i=i+1;
    a=rand()%L,b=rand()%L;
    m=rand()%6,n=(m+2)%6;
    c=B[m][0];d=B[m][1];e=B[n][0];f=B[n][1];
    fk=(a+c+L)%L;fk1=(b+d+L)%L;
    fkk=(a+e+L)%L;fk11=(b+f+L)%L;
    if(A[a][b]==0&&A[fk][fk1]==0&&A[fkk][fk11]==0)
    {
        A[a][b]=-i-N,A[fk][fk1]=-i-N,A[fkk][fk11]=-i-N;
        P2[i-1][0]=a,P2[i-1][1]=b;
        D2[i-1][0]=m,D2[i-1][1]=n;
        A1[a][b]=-i-N;
    }
    else
    {
        i=i-1;
    }

};

/*for(int i=0;i<L;i++)
{
    for(int j=0;j<L;j++)
    {
        cout<<A[L-1-i][j]<<" ";
    }
    cout<<"\n";
}*/
/*上面ok了*/

/*分子开始运动*/
int k=9900000,k1=0;
double T=1,F1=0,F2=0,Eo=0,Etol=0,En=0,E1=0,E2=0,E3=0,E4=0;
int Insum11=0,Insum12=0,in11=0,in12=0;
int a1,b1,c1,d1,m1,D,DD1,DD2,Dx1,Dx2,Dk,wtf,wtf1;
k1=0;

while(T>0.01)
{
    if(k1%10000==0)
    {
        T=T*0.999;
        cout<<T<<"  "<<Etol<<'\n';
    };
    /*先动的是金属原子 */

    for(int k2=0;k2<N+M;k2++)
    {
    int rd1=rand()%(N+M);
    if(rd1<N)
    {
        i=rd1;
        Eo=0,En=0;
        in11=0,in12=0;
        a=rand()%L,b=rand()%L;
        D=rand()%3;
        D=(D1[i]+D+5)%6;
        if(A[a][b]==0&&A[(a-B[D][0]+L)%L][(b-B[D][1]+L)%L]==0)
        {
            for(int sz=-1;sz<2;sz++)
            {
                c=(sz+D1[i]+6)%6;
                /*cout<<A2[1+3*i+B[c][0]][1+B[c][1]]<<"   "<<sz<<'\n';*/
                if(A[(P1[i][0]+B[c][0]+L)%L][(P1[i][1]+B[c][1]+L)%L]>0&&sz==0)
                {
                    d=A[(P1[i][0]+B[c][0]+L)%L][(P1[i][1]+B[c][1]+L)%L];
                    if(A2[1+3*(d-1)-B[c][0]][1-B[c][1]]==2)
                    {
                        Eo=Eo+0.5;
                        in11=in11-1;
                    }
                }

                if(A[(P1[i][0]+B[c][0]+L)%L][(P1[i][1]+B[c][1]+L)%L]>0&&sz==1)
                {
                    d=A[(P1[i][0]+B[c][0]+L)%L][(P1[i][1]+B[c][1]+L)%L];
                    if(A2[1+3*(d-1)-B[c][0]][1-B[c][1]]==1)
                    {
                        Eo=Eo+0.5;
                        in11=in11-1;
                    }
                }









                if(A[(P1[i][0]+B[c][0]+L)%L][(P1[i][1]+B[c][1]+L)%L]<-N&&sz==-1)
                {
                    fk=(P1[i][0]+2*B[c][0]+L)%L;
                    fkk=(P1[i][1]+2*B[c][1]+L)%L;
                    fk1=(P1[i][0]+B[c][0]+L)%L;
                    fk11=(P1[i][1]+B[c][1]+L)%L;
                    if(A1[fk][fkk]==A[fk1][fk11])
                    {
                        Eo=Eo+1;
                        in12=in12-1;
                    }
                }
            }

            A[P1[i][0]][P1[i][1]]=0;
            A[(P1[i][0]-B[D1[i]][0]+L)%L][(P1[i][1]-B[D1[i]][1]+L)%L]=0;

            for(int sz=-1;sz<2;sz++)
            {
                c=(sz+D+6)%6;
                if(A[(a+B[c][0]+L)%L][(b+B[c][1]+L)%L]>0&&sz==0)
                {
                    d=A[(a+B[c][0]+L)%L][(b+B[c][1]+L)%L];
                    if(A2[1+3*(d-1)-B[c][0]][1-B[c][1]]==2)
                    {
                        En=En+0.5;
                        in11=in11+1;
                    }
                }


                if(A[(a+B[c][0]+L)%L][(b+B[c][1]+L)%L]>0&&sz==1)
                {
                    d=A[(a+B[c][0]+L)%L][(b+B[c][1]+L)%L];
                    if(A2[1+3*(d-1)-B[c][0]][1-B[c][1]]==1)
                    {
                        En=En+0.5;
                        in11=in11+1;
                    }
                }





                if(A[(a+B[c][0]+L)%L][(b+B[c][1]+L)%L]<-N&&sz==-1)
                {
                    fk=(a+2*B[c][0]+L)%L;
                    fkk=(b+2*B[c][1]+L)%L;
                    fk1=(a+B[c][0]+L)%L;
                    fk11=(b+B[c][1]+L)%L;
                    if(A1[fk][fkk]==A[fk1][fk11])
                    {
                        En=En+1;
                        in12=in12+1;
                    }
                }
            }

            A[P1[i][0]][P1[i][1]]=i+1;
            A[(P1[i][0]-B[D1[i]][0]+L)%L][(P1[i][1]-B[D1[i]][1]+L)%L]=-i-1;

            F1=min(1,exp((En-Eo)/T));
            F2=(rand()+1)/(RAND_MAX+1.0);
            if(F2<F1)
            {
                A2[1+3*i+B[D1[i]][0]][1+B[D1[i]][1]]=0;
                A2[1+3*i+B[(D1[i]+1+6)%6][0]][1+B[(D1[i]+1+6)%6][1]]=0;
                A2[1+3*i+B[(D1[i]-1+6)%6][0]][1+B[(D1[i]-1+6)%6][1]]=0;
                A[P1[i][0]][P1[i][1]]=0,A[(P1[i][0]-B[D1[i]][0]+L)%L][(P1[i][1]-B[D1[i]][1]+L)%L]=0;
                A[a][b]=i+1,A[(a-B[D][0]+L)%L][(b-B[D][1]+L)%L]=-i-1;
                P1[i][0]=a,P1[i][1]=b,D1[i]=D,Etol=Etol+Eo-En;
                A2[1+3*i+B[D1[i]][0]][1+B[D1[i]][1]]=1;
                A2[1+3*i+B[(D1[i]+1+6)%6][0]][1+B[(D1[i]+1+6)%6][1]]=2;
                A2[1+3*i+B[(D1[i]-1+6)%6][0]][1+B[(D1[i]-1+6)%6][1]]=3;
                E1=E1+Eo;E2=E2+En;
                Insum11=Insum11+in11;
                Insum12=Insum12+in12;
            }
        }
    }
    /*再动分子2*/

    if(rd1>=N)
    {
        i=rd1-N;
        Eo=0,En=0;
        in11=0,in12=0;
        a=rand()%L,b=rand()%L;
        Dk=(rand()%3)-1;
        Dx1=D2[i][0];
        Dx2=D2[i][1];
        DD1=(D2[i][0]+Dk+6)%6;
        DD2=(D2[i][1]+Dk+6)%6;
        a1=P2[i][0],b1=P2[i][1];
        if(A[a][b]==0&&A[(a+B[DD1][0]+L)%L][(b+B[DD1][1]+L)%L]==0&&A[(a+B[DD2][0]+L)%L][(b+B[DD2][1]+L)%L]==0)
        {
            if(A[(a1+2*B[Dx1][0]+L)%L][(b1+2*B[Dx1][1]+L)%L]>0)
            {
                d=A[(a1+2*B[Dx1][0]+L)%L][(b1+2*B[Dx1][1]+L)%L];
                fk=1+3*(d-1)-B[Dx1][0];
                fk1=1-B[Dx1][1];
                if(A2[fk][fk1]==3)
                {
                    Eo=Eo+1;
                    in12=in12-1;
                }
            }

            if(A[(a1+2*B[Dx2][0]+L)%L][(b1+2*B[Dx2][1]+L)%L]>0)
            {
                d=A[(a1+2*B[Dx2][0]+L)%L][(b1+2*B[Dx2][1]+L)%L];
                fk=1+3*(d-1)-B[Dx2][0];
                fk1=1-B[Dx2][1];
                if(A2[fk][fk1]==3)
                {
                    Eo=Eo+1;
                    in12=in12-1;
                }
            }



            A[a1][b1]=0;
            A1[a1][b1]=0;
            A[(a1+B[Dx1][0]+L)%L][(b1+B[Dx1][1]+L)%L]=0;
            A[(a1+B[Dx2][0]+L)%L][(b1+B[Dx2][1]+L)%L]=0;
            if(A[(a+2*B[DD1][0]+L)%L][(b+2*B[DD1][1]+L)%L]>0)
            {
                d=A[(a+2*B[DD1][0]+L)%L][(b+2*B[DD1][1]+L)%L];
                fk=1+3*(d-1)-B[DD1][0];
                fk1=1-B[DD1][1];
                if(A2[fk][fk1]==3)
                {
                    En=En+1;
                    in12=in12+1;
                }
            }

            if(A[(a+2*B[DD2][0]+L)%L][(b+2*B[DD2][1]+L)%L]>0)
            {
                d=A[(a+2*B[DD2][0]+L)%L][(b+2*B[DD2][1]+L)%L];
                fk=1+3*(d-1)-B[DD2][0];
                fk1=1-B[DD2][1];
                if(A2[fk][fk1]==3)
                {
                    En=En+1;
                    in12=in12+1;
                }
            }
            A[a1][b1]=-N-i-1;
            A1[a1][b1]=-N-i-1;
            A[(a1+B[Dx1][0]+L)%L][(b1+B[Dx1][1]+L)%L]=-N-i-1;
            A[(a1+B[Dx2][0]+L)%L][(b1+B[Dx2][1]+L)%L]=-N-i-1;

            F1=min(1,exp((En-Eo)/T));
            F2=(rand()+1)/(RAND_MAX+1.0);
            if(F2<F1)
            {
                A[P2[i][0]][P2[i][1]]=0;
                A1[P2[i][0]][P2[i][1]]=0;
                A[(P2[i][0]+B[D2[i][0]][0]+L)%L][(P2[i][1]+B[D2[i][0]][1]+L)%L]=0;
                A[(P2[i][0]+B[D2[i][1]][0]+L)%L][(P2[i][1]+B[D2[i][1]][1]+L)%L]=0;
                P2[i][0]=a,P2[i][1]=b;
                D2[i][0]=(D2[i][0]+Dk+6)%6;D2[i][1]=(D2[i][1]+Dk+6)%6;
                Etol=Etol+Eo-En;
                A[P2[i][0]][P2[i][1]]=-N-i-1;
                A1[P2[i][0]][P2[i][1]]=-N-i-1;
                A[(P2[i][0]+B[D2[i][0]][0]+L)%L][(P2[i][1]+B[D2[i][0]][1]+L)%L]=-N-i-1;
                A[(P2[i][0]+B[D2[i][1]][0]+L)%L][(P2[i][1]+B[D2[i][1]][1]+L)%L]=-N-i-1;
                E3=E3+Eo;E4=E4+En;
                Insum11=Insum11+in11;
                Insum12=Insum12+in12;
            }
        }
    }
    }
    k1=k1+1;


};

/*cout<<Etol<<"\n";
cout<<E1<<"   "<<E2<<"   "<<E3<<"   "<<E4<<"   "<<'\n';
cout<<Insum11<<"   "<<Insum12<<'\n';
for(int i=0;i<L;i++)
{
    for(int j=0;j<L;j++)
    {
        cout<<A[L-1-i][j]<<" ";
    }
    cout<<"\n";
}

for(int i=0;i<3*N;i++)
{
    for(int j=0;j<3;j++)
    {
        cout<<A2[i][j]<<" ";
    }
    cout<<"\n";
}

for(int i=0;i<L;i++)
{
    for(int j=0;j<L;j++)
    {
        cout<<A1[L-1-i][j]<<" ";
    }
    cout<<"\n";
}*/





ofstream outfile;
outfile.open("ndata.txt");
for(int i=0;i<N;i++)
{
    outfile<<P1[i][0]<<endl;
    outfile<<P1[i][1]<<endl;
}

outfile.close();

outfile.open("ndata1.txt");
for(int i=0;i<M;i++)
{
    outfile<<P2[i][0]<<endl;
    outfile<<P2[i][1]<<endl;
}

outfile.close();

outfile.open("ndata2.txt");
for(int i=0;i<N;i++)
{
    outfile<<D1[i]<<endl;
}

outfile.close();

outfile.open("ndata3.txt");
for(int i=0;i<M;i++)
{
    outfile<<D2[i][0]<<endl;
    outfile<<D2[i][1]<<endl;
}

outfile.close();
}

double min(double x,double y)
{
    double z;
    if(x<y)
    {
        z=x;
    };
    if(x>=y)
    {
        z=y;
    }
    return z;
}
