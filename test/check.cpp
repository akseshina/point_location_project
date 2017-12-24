#include<bits/stdc++.h>
using namespace std;

#ifdef WIN32
	#define lld "%I64d"
	#define llu "%I64u"
#else
	#define lld "%lld"
	#define llu "%llu"
#endif

typedef unsigned int uint;
typedef long long LL;
typedef unsigned long long ULL;
typedef double dbl;
typedef long double LD;
typedef vector<int> vi;
typedef pair<int,int> pii;
#define pb push_back
#define mp make_pair
#define fst first
#define snd second
#define sz(x) (int)(x.size())

#define forn(i,n) for(int i=0;i<(n);++i)
#define fornr(i,n) for(int i=(n)-1;i>=0;--i)
#define forab(i,a,b) for(int i=(a);i<(b);++i)
#define forba(i,a,b) for(int i=(b)-1;i>=(a);--i)
#define forit(i,A) for(__typeof((A).begin()) i=(A).begin();i!=(A).end();++i)
#define all(A) (A).begin(),(A).end()

char S[100];

vi f(const char *fname){
	assert(freopen(fname,"r",stdin));
	vi v;
	while(gets(S)){
		if(!strcmp(S, ""))
			break;
		else if(!strcmp(S, "INSIDE"))
			v.pb(1);
		else if(!strcmp(S, "OUTSIDE"))
			v.pb(-1);
		else if(!strcmp(S, "BORDER"))
			v.pb(0);
		else
			assert(0);
	}
	return v;
}

int main(){
	vi v1 = f("test.out"), v2 = f("test2.out");

	assert(v1.size() == v2.size());
	forn(i, sz(v1))
		if(v2[i] != 0)
			assert(v1[i] == v2[i]);

	return 0;
}
