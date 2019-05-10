/*120度分子与金属原子模拟*/
/*上次得到的数据结果抖动严重，十分怀疑是代码问题*/
/*抖动问题已解决，盖由随机数选择引起。但无法得到高阶三角形*/
#include<iostream>
#include<cmath>
#include<fstream>
#include <stdlib.h>
#include<time.h>
using namespace std;
int main()
{
    double min(double x,double y);
    double  rand1();
    const int L=200,N=1200,M=1800;
    srand((int)time(0));

    double S0[1100],S1[1100],S2[1100],S3[1100];
    int S=0;
    for(int i=0;i<1100;i++)
    {
        S0[i]=0;S1[i]=0;S2[i]=0;S3[i]=0;
    }
    int A[L][L],A1[L][L];/*初始化矩阵*/
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            A[i][j]=0;
            A1[i][j]=0;
        }
    }

    int B[6][2];/*方向矩阵*/
    B[0][0]=-1,B[0][1]=0,B[1][0]=-1,B[1][1]=1,B[2][0]=0,B[2][1]=1,
    B[3][0]=1,B[3][1]=0,B[4][0]=1,B[4][1]=-1,B[5][0]=0,B[5][1]=-1;

    /*放置金属原子*/
    int P1[N][2],a,b,fk0,fk1,x,y,z;
    int i=0;
    while(i<N)
    {
        a=rand()%L;b=rand()%L;
        if(A[a][b]==0)
        {
            x=0;
            for(int j=0;j<6;j++)
            {
                fk0=B[j][0],fk1=B[j][1];
                if(A[(a+fk0+L)%L][(b+fk1+L)%L]==0)
                {
                    x=x+1;
                }
            }
            if(x==6)
            {

                A[a][b]=i+1;
                P1[i][0]=a,P1[i][1]=b;
                i=i+1;
            }
        }
    }

    /*放置连接分子*/
    int P2[M][2],D[M],c,fk2,fk3;
    i=0;
    while(i<M)
    {
        a=rand()%L,b=rand()%L,c=rand()%6;
        fk0=B[c][0],fk1=B[c][1],fk2=B[(c+2)%6][0],fk3=B[(c+2)%6][1];
        if(A[a][b]==0&&A[(a+fk0+L)%L][(b+fk1+L)%L]==0&&A[(a+fk2+L)%L][(b+fk3+L)%L]==0)
        {
            x=0;
            for(int j=0;j<6;j++)
            {
                fk2=B[j][0],fk3=B[j][1];
                if(A[(a+fk2+L)%L][(b+fk3+L)%L]==0)
                {
                    x=x+1;
                }
            }
            fk0=B[c][0],fk1=B[c][1];
            for(int j=0;j<6;j++)
            {
                if(j!=c)
                {
                    fk2=B[j][0],fk3=B[j][1];
                    if(A[(a+fk2+fk0+L)%L][(b+fk3+fk1+L)%L]==0)
                    {
                        x=x+1;
                    }
                }
            }
            fk0=B[(c+2)%6][0],fk1=B[(c+2)%6][1];
            for(int j=0;j<6;j++)
            {

                if(j!=(c+2)%6)
                {
                    fk2=B[j][0],fk3=B[j][1];
                    if(A[(a+fk2+fk0+L)%L][(b+fk3+fk1+L)%L]==0)
                    {
                        x=x+1;
                    }
                }
            }
            y=0;
            if(x==16)
            {
                fk0=B[c][0],fk1=B[c][1];
                if(A[(a+2*fk0+L)%L][(b+2*fk1+L)%L]>0)
                {
                    z=A[(a+2*fk0+L)%L][(b+2*fk1+L)%L];
                    fk2=B[(c+4)%6][0],fk3=B[(c+4)%6][1];
                    if(A[(a+2*fk0+L+fk2)%L][(b+2*fk1+L+fk3)%L]==0)
                    {
                        fk2=B[(c+2)%6][0],fk3=B[(c+2)%6][1];
                        if(A[(a+2*fk0+L+fk2)%L][(b+2*fk1+L+fk3)%L]==0)
                        {
                            y=1;
                        }
                    }
                }
                if(A[(a+2*fk0+L)%L][(b+2*fk1+L)%L]==0)
                {
                    y=1;
                }
                fk0=B[(c+2)%6][0],fk1=B[(c+2)%6][1];
                if(A[(a+2*fk0+L)%L][(b+2*fk1+L)%L]>0)
                {
                    z=A[(a+2*fk0+L)%L][(b+2*fk1+L)%L];
                    fk2=B[(c+4)%6][0],fk3=B[(c+4)%6][1];
                    if(A[(a+2*fk0+L+fk2)%L][(b+2*fk1+L+fk3)%L]==0)
                    {
                        fk2=B[(c+6)%6][0],fk3=B[(c+6)%6][1];
                        if(A[(a+2*fk0+L+fk2)%L][(b+2*fk1+L+fk3)%L]==0)
                        {
                            y=y+1;
                        }
                    }
                }
                if(A[(a+2*fk0+L)%L][(b+2*fk1+L)%L]==0)
                {
                    y=y+1;
                }
                if(y==2)
                {
                    A[a][b]=-i-1;
                    fk0=B[c][0],fk1=B[c][1],fk2=B[(c+2)%6][0],fk3=B[(c+2)%6][1];
                    A[(a+fk0+L)%L][(b+fk1+L)%L]=-i-1,A[(a+fk2+L)%L][(b+fk3+L)%L]=-i-1;
                    P2[i][0]=a,P2[i][1]=b,D[i]=c;
                    A1[a][b]=-1;
                    i=i+1;
                    cout<<i<<'\n';
                }
            }


        }
    }


    /*分子运动*/
    int k=10580000,k1=0,rd1=0;
    double T=1,T1=T,F1=0,F2=0,Etol=0,E1=0,E2=0;
    while(T>=0.01)
    {
        if(k1%1000==0)
        {
            T=T-0.001;
            cout<<T<<" "<<Etol<<'\n';
        }
        for(int k2=0;k2<(N+M);k2++)
        {
            rd1=rand()%(N+M);
            if(rd1<N)
            {
                int i=rd1;
                E1=0,E2=0;
                a=P1[i][0],b=P1[i][1];
                for(int j=0;j<6;j++)
                {
                    fk0=B[j][0],fk1=B[j][1];
                    if(A[(a+fk0+L)%L][(b+fk1+L)%L]<0)
                    {
                        E1=E1+1;
                    }
                }
                a=rand()%L,b=rand()%L;
                if(A[a][b]==0)
                {
                    x=0,y=0,z=0,E2=0;
                    for(int j=0;j<6;j++)
                    {
                        fk0=B[j][0],fk1=B[j][1];
                        if(A[(a+fk0+L)%L][(b+fk1+L)%L]>0)
                        {
                            x=x+1;
                        }
                        if(A[(a+fk0+L)%L][(b+fk1+L)%L]<0)
                        {
                            if(A[(a+fk0+L)%L][(b+fk1+L)%L]==A[(a+2*fk0+L)%L][(b+2*fk1+L)%L])
                            {
                                if(A1[(a+fk0+L)%L][(b+fk1+L)%L]==0)
                                {
                                    E2=E2+1;
                                }
                                else
                                {
                                    y=y+1;
                                }
                            }
                            else
                            {
                                y=y+1;
                            }
                        }
                        fk1=B[(j+5)%6][0],fk2=B[(j+5)%6][1];
                        if(A[(a+fk0+L)%L][(b+fk1+L)%L]!=0&&A[(a+fk1+L)%L][(b+fk2+L)%L]!=0)
                        {
                            z=z+1;
                        }
                    }



                   /* int Switch=0,f=0,f1=0,fold=0,fk00,fk11,fk22,fk33;
                    for(int l=0;l<6;l++)
                    {
                        fk00=B[l][0],fk11=B[l][1];
                        if(A[(a+fk00+L)%L][(b+fk11+L)%L]<0)
                        {
                            fold=fold+1;f=l;
                        }
                    }

                    if(fold==3)
                    {
                        f=f%2;
                        int Newf=0;
                        for(int l=0;l<3;l++)
                        {
                            Newf=f+2*l;
                            fk00=B[Newf][0],fk11=B[Newf][1];
                            fk22=B[(Newf+1)%6][0],fk33=B[(Newf+1)%6][1];
                            if(A[(a+fk00+L)%L][(b+fk11+L)%L]==A[(a+2*fk00+fk22+L)%L][(b+2*fk11+fk33+L)%L])
                            {
                                f1=f1+1;
                            }
                        }
                        if(f1==0||f1==3)
                        {
                            Switch=1;
                        }
                        f1=0;
                    }*/





                    if(x==0&&y==0&&z==0)
                    {
                        int Switch=0,f=0,f1=0,fk00,fk11,fk22,fk33;
                        if(E2==3)
                        {
                            fk00=B[1][0],fk11=B[1][1];
                            if(A[(a+fk00+L)%L][(b+fk11+L)%L]<0)
                            {
                                f=1;
                            }
                            for(int l=0;l<3;l++)
                            {
                                f=f+2*l;
                                fk00=B[f][0],fk11=B[f][1],fk22=B[(f+1)%6][0],fk33=B[(f+1)%6][1];
                                if(A[(a+fk00+L)%L][(b+fk11+L)%L]==A[(a+2*fk00+fk22+L)%L][(b+2*fk11+fk33+L)%L])
                                {
                                    f1=f1+1;
                                }
                            }
                            cout<<"f1= "<<f1<<endl;
                            if(f1==0||f1==3)
                            {
                                Switch=1;
                            }
                        }

                        F1=min(1,exp(-(E1-E2)/T));
                        F2=rand1();
                        if(F2<F1&&Switch==0)
                        {
                            /*cout<<"这是金属原子"<<"  "<<E1<<"  "<<E2<<"  "<<F1<<" "<<F2<<"  "<<"改变"<<"  "<<Etol<<endl;*/
                            fk0=P1[i][0],fk1=P1[i][1];A[fk0][fk1]=0,A[a][b]=i+1;
                            P1[i][0]=a,P1[i][1]=b;Etol=Etol+E1-E2;
                            cout<<"金属原子   "<<Etol<<"  "<<E1<<"  "<<E2<<endl;
                        }
                    }
                }
            }
            if(rd1>=N)
            {
                rd1=rd1-N;
                int i=rd1;
                E1=0,E2=0;
                a=P2[i][0],b=P2[i][1],c=D[i];
                fk0=B[c][0],fk1=B[c][1];
                if(A[(a+2*fk0+L)%L][(b+2*fk1+L)%L]>0)
                {
                    E1=E1+1;
                }
                fk0=B[(c+2)%6][0],fk1=B[(c+2)%6][1];
                if(A[(a+2*fk0+L)%L][(b+2*fk1+L)%L]>0)
                {
                    E1=E1+1;
                }
                a=rand()%L,b=rand()%L,c=(c-1+rand()%3+6)%6;
                fk0=B[c][0],fk1=B[c][1],fk2=B[(c+2)%6][0],fk3=B[(c+2)%6][1];
                if(A[a][b]==0&&A[(a+fk0+L)%L][(b+fk1+L)%L]==0&&A[(a+fk2+L)%L][(b+fk3+L)%L]==0)
                {
                    x=0;
                    for(int j=0;j<6;j++)
                    {
                        fk2=B[j][0],fk3=B[j][1];
                        if(A[(a+fk2+L)%L][(b+fk3+L)%L]==0)
                        {
                            x=x+1;
                        }
                    }
                    fk0=B[c][0],fk1=B[c][1];
                    for(int j=0;j<6;j++)
                    {
                        fk2=B[j][0],fk3=B[j][1];
                        if(j!=c)
                        {
                            if(A[(a+fk2+fk0+L)%L][(b+fk3+fk1+L)%L]==0)
                            {
                                x=x+1;
                            }
                        }
                    }
                    fk0=B[(c+2)%6][0],fk1=B[(c+2)%6][1];
                    for(int j=0;j<6;j++)
                    {
                        fk2=B[j][0],fk3=B[j][1];
                        if(j!=(c+2)%6)
                        {
                            if(A[(a+fk2+fk0+L)%L][(b+fk3+fk1+L)%L]==0)
                            {
                                x=x+1;
                            }
                        }
                    }
                    y=0;








                    fk0=B[c][0],fk1=B[c][1],fk2=B[(c+2)%6][0],fk3=B[(c+2)%6][1];
                    A[a][b]=-i-1,A[(a+fk0+L)%L][(b+fk1+L)%L]=-i-1,A[(a+fk2+L)%L][(b+fk3+L)%L]=-i-1;


                    int Switch=0;
                    if(x==16)
                    {
                        fk0=B[c][0],fk1=B[c][1];
                        if(A[(a+2*fk0+L)%L][(b+2*fk1+L)%L]>0)
                        {
                            z=A[(a+2*fk0+L)%L][(b+2*fk1+L)%L];
                            fk2=B[(c+4)%6][0],fk3=B[(c+4)%6][1];
                            if(A[(a+2*fk0+L+fk2)%L][(b+2*fk1+L+fk3)%L]==0)
                            {
                                fk2=B[(c+2)%6][0],fk3=B[(c+2)%6][1];
                                if(A[(a+2*fk0+L+fk2)%L][(b+2*fk1+L+fk3)%L]==0)
                                {
                                    y=y+1;E2=E2+1;
                                    /*int f=0,f1=0,fold=0,fk00,fk11,fk22,fk33;
                                    for(int l=0;l<6;l++)
                                    {
                                        fk00=B[l][0],fk11=B[l][1];
                                        if(A[(a+2*fk0+fk00+L)%L][(b+2*fk1+fk11+L)%L]<0)
                                        {
                                            fold=fold+1;f=l;
                                        }
                                    }
                                    if(fold==3)
                                    {
                                        f=f%2;
                                        int Newf=0;
                                        for(int l=0;l<3;l++)
                                        {
                                            Newf=f+2*l;
                                            fk00=B[Newf][0],fk11=B[Newf][1];
                                            fk22=B[(Newf+1)%6][0],fk33=B[(Newf+1)%6][1];
                                            if(A[(a+2*fk0+2*fk00+L)%L][(b+2*fk1+2*fk11+L)%L]==A[(a+2*fk0+2*fk00+2*fk22+L)%L][(b+2*fk1+2*fk11+2*fk33+L)%L])
                                            {
                                                f1=f1+1;
                                            }
                                        }
                                        if(f1==0||f1==3)
                                        {
                                            Switch=1;
                                        }
                                        f1=0;

                                    }*/
                                }
                            }
                        }
                        else if(A[(a+2*fk0+L)%L][(b+2*fk1+L)%L]==0)
                        {
                            y=y+1;
                        }
                        fk0=B[(c+2)%6][0],fk1=B[(c+2)%6][1];
                        if(A[(a+2*fk0+L)%L][(b+2*fk1+L)%L]>0)
                        {
                            z=A[(a+2*fk0+L)%L][(b+2*fk1+L)%L];
                            fk2=B[(c+4)%6][0],fk3=B[(c+4)%6][1];
                            if(A[(a+2*fk0+L+fk2)%L][(b+2*fk1+L+fk3)%L]==0)
                            {
                                fk2=B[(c+6)%6][0],fk3=B[(c+6)%6][1];
                                if(A[(a+2*fk0+L+fk2)%L][(b+2*fk1+L+fk3)%L]==0)
                                {
                                    y=y+1;E2=E2+1;
                                    /*int f=0,f1=0,fold=0,fk00,fk11,fk22,fk33;
                                    for(int l=0;l<6;l++)
                                    {
                                        fk00=B[l][0],fk11=B[l][1];
                                        if(A[(a+2*fk0+fk00+L)%L][(b+2*fk1+fk11+L)%L]<0)
                                        {
                                            fold=fold+1;f=l;
                                        }
                                    }
                                    if(fold==3)
                                    {
                                        f=f%2;
                                        int Newf=0;
                                        for(int l=0;l<3;l++)
                                        {
                                            Newf=f+2*l;
                                            fk00=B[Newf][0],fk11=B[Newf][1];
                                            fk22=B[(Newf+1)%6][0],fk33=B[(Newf+1)%6][1];
                                            if(A[(a+2*fk0+2*fk00+L)%L][(b+2*fk1+2*fk11+L)%L]==A[(a+2*fk0+2*fk00+2*fk22+L)%L][(b+2*fk1+2*fk11+2*fk33+L)%L])
                                            {
                                                f1=f1+1;
                                            }
                                        }
                                        if(f1==0||f1==3)
                                        {
                                            Switch=1;
                                        }
                                        f1=0;
                                    }*/
                                }
                            }
                        }
                        else if(A[(a+2*fk0+L)%L][(b+2*fk1+L)%L]==0)
                        {
                            y=y+1;
                        }


                        fk0=B[c][0],fk1=B[c][1],fk2=B[(c+2)%6][0],fk3=B[(c+2)%6][1];
                        A[a][b]=0,A[(a+fk0+L)%L][(b+fk1+L)%L]=0,A[(a+fk2+L)%L][(b+fk3+L)%L]=0;

                        if(y==2)
                        {
                            F1=min(1,exp(-(E1-E2)/T));
                            F2=rand1();
                            if(F2<F1)
                            {

                                fk0=P2[i][0],fk1=P2[i][1];
                                fk2=B[D[i]][0],fk3=B[D[i]][1];
                                A1[fk0][fk1]=0;
                                A[fk0][fk1]=0,A[(fk0+fk2+L)%L][(fk1+fk3+L)%L]=0;
                                fk2=B[(D[i]+2)%6][0],fk3=B[(D[i]+2)%6][1];
                                A[(fk0+fk2+L)%L][(fk1+fk3+L)%L]=0;
                                P2[i][0]=a,P2[i][1]=b,D[i]=c;
                                fk0=P2[i][0],fk1=P2[i][1];
                                fk2=B[D[i]][0],fk3=B[D[i]][1];
                                A1[fk0][fk1]=-1;
                                A[fk0][fk1]=-i-1,A[(fk0+fk2+L)%L][(fk1+fk3+L)%L]=-i-1;
                                fk2=B[(D[i]+2)%6][0],fk3=B[(D[i]+2)%6][1];
                                A[(fk0+fk2+L)%L][(fk1+fk3+L)%L]=-i-1;
                                Etol=Etol+E1-E2;
                                cout<<Etol<<"  "<<E1<<"  "<<E2<<endl;
                                /*cout<<"这是有机分子"<<"  "<<E1<<"  "<<E2<<"  "<<F1<<" "<<F2<<"  "<<"改变"<<"  "<<Etol<<endl;*/
                            }
                        }
                    }
                }
            }

        }
        k1=k1+1;
        if(T1!=T)
        {
            T1=T;
            double num0=0,num1=0,num2=0,num3=0,zj=0;
            for(int wh=0;wh<N;wh++)
            {
                int zj2=0;
                for(int wh1=0;wh1<6;wh1++)
                {
                    zj=A[(P1[wh][0]+B[wh1][0]+L)%L][(P1[wh][1]+B[wh1][1]+L)%L];
                    if(zj<0)
                    {
                        zj2=zj2+1;
                    }
                }
                if(zj2==0)
                {
                    num0=num0+1;
                }
                if(zj2==1)
                {
                    num1=num1+1;
                }
                if(zj2==2)
                {
                    num2=num2+1;
                }
                if(zj2==3)
                {
                    num3=num3+1;
                }
            }
            S0[S]=num0/N;S1[S]=num1/N;S2[S]=num2/N;S3[S]=num3/N;
            S=S+1;
        }
        if(k1==9000000)
        {
            ofstream outfile;
            outfile.open("dataA.txt");
            for(int i=0;i<N;i++)
            {
                 outfile<<P1[i][0]<<endl;outfile<<P1[i][1]<<endl;
            }
            outfile.close();
            outfile.open("data1A.txt");
            for(int i=0;i<M;i++)
            {
                 outfile<<P2[i][0]<<endl;outfile<<P2[i][1]<<endl;
            }
            outfile.close();
            outfile.open("data2A.txt");
            for(int i=0;i<M;i++)
            {
                 outfile<<D[i]<<endl;
            }
            outfile.close();
        }
        if(k1==9010000)
        {
            ofstream outfile;
            outfile.open("dataB.txt");
            for(int i=0;i<N;i++)
            {
                 outfile<<P1[i][0]<<endl;outfile<<P1[i][1]<<endl;
            }
            outfile.close();
            outfile.open("data1B.txt");
            for(int i=0;i<M;i++)
            {
                 outfile<<P2[i][0]<<endl;outfile<<P2[i][1]<<endl;
            }
            outfile.close();
            outfile.open("data2B.txt");
            for(int i=0;i<M;i++)
            {
                 outfile<<D[i]<<endl;
            }
            outfile.close();
        }
    }



    /*用来输出矩阵*/
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            cout<<A[i][j]<<' ';
        }
        cout<<'\n';
    }

ofstream outfile;
outfile.open("data.txt");

for(int i=0;i<N;i++)
{
    outfile<<P1[i][0]<<endl;outfile<<P1[i][1]<<endl;
}
outfile.close();

outfile.open("data1.txt");
for(int i=0;i<M;i++)
{
    outfile<<P2[i][0]<<endl;outfile<<P2[i][1]<<endl;
}
outfile.close();

outfile.open("data2.txt");
for(int i=0;i<M;i++)
{
    outfile<<D[i]<<endl;
}
outfile.close();

outfile.open("dataNew0.txt");
for(int i=0;i<1000;i++)
{
    outfile<<S0[i]<<endl;
}
outfile.close();

outfile.open("dataNew1.txt");
for(int i=0;i<1000;i++)
{
    outfile<<S1[i]<<endl;
}
outfile.close();

outfile.open("dataNew2.txt");
for(int i=0;i<1000;i++)
{
    outfile<<S2[i]<<endl;
}
outfile.close();

outfile.open("dataNew3.txt");
for(int i=0;i<1000;i++)
{
    outfile<<S3[i]<<endl;
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
double  rand1()
{
    double result=0,k=1;
    for(int i=0;i<9;i++)
    {
        k=k*10;
        result=result+(rand()%10)/k;
    }
    return result;
}
