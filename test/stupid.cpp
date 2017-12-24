#include<bits/stdc++.h>

using namespace std;

int main()
{
    int  N, i,cc;
    scanf("%d",&cc);
    while(cc--){
    cin>>N;

    long long X[N+1], Y[N+1];
    for(i=0;i<N;++i)
    {
        cin>>X[i]>>Y[i];
    }
    X[N]=X[0];
    Y[N]=Y[0];

    int K;
    scanf("%d", &K);
    while(K--){
    int Xs, Ys;
    scanf("%d%d", &Xs, &Ys);
    double s=0;
    
    bool ok = 0;
    for(i=0;i<N;++i)
    {
        long long scalar=(X[i]-Xs)*(X[i+1]-Xs)+(Y[i]-Ys)*(Y[i+1]-Ys);
        long long vector=(X[i]-Xs)*(Y[i+1]-Ys)-(Y[i]-Ys)*(X[i+1]-Xs);
        
        if( vector == 0 && scalar <= 0 )
        {
            cout<<"BORDER"<<endl;
            ok = 1;
            break;
        }

        s+=atan2(vector, scalar);
    }
    
    if(!ok){
    if( s>-1 && s<1)
        cout<<"OUTSIDE"<<endl;
    else
        cout<<"INSIDE"<<endl;
        }
	}
}    
    return 0;
}
