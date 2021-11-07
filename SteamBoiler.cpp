//*Author : Franklin*
#include<bits/stdc++.h>
using namespace std;
#define sim template < class c
#define ris return * this
#define dor > debug & operator <<
#define eni(x) sim > typename enable_if<sizeof dud<c>(0) x 1, debug&>::type operator<<(c i)
sim > struct rge { c b, e; };
sim > rge<c> range(c i, c j) { return rge<c>{i, j}; }
sim > auto dud(c* x) -> decltype(cerr << *x, 0);
sim > char dud(...);
struct debug {
#ifndef LOCAL
    ~debug() { cerr << endl; }
    eni(!=) { cerr << boolalpha << i; ris; }
    eni(==) { ris << range(begin(i), end(i)); }
    sim, class b dor(pair < b, c > d) {
        ris << "(" << d.first << ", " << d.second << ")";
    }
    sim dor(rge<c> d) {
        *this << "[";
        for (auto it = d.b; it != d.e; ++it) *this << &", " [ 2 * (it == d.b)] << *it;
        ris << "]";
    }
#else
    sim dor(const c&) { ris; }
#endif
};
#define imie(...) "" << #__VA_ARGS__ " = " << (__VA_ARGS__) << ""
#define arr array
#define eps 1e-4
#define deb(x) cout << #x << " = " << x << "\n"

using ll = long long;
const int INF = 1e9 + 5;
const int N = 1e6 + 5;
const double EPS = 1e-6;

struct Solver {
    vector<double> p_list = {0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
    vector<double> vf_list = {0.0011, 0.0011, 0.0011, 0.0011, 0.0011, 0.0011, 0.0011};
    vector<double> vg_list = {0.273, 0.255, 0.24, 0.229, 0.215, 0.204, 0.194};
    vector<double> Uf_list = {696.33, 708.475, 720.02, 731.065, 741.61, 751.755, 761.5};
    vector<double> Ug_list = {2570.9, 2573.75, 2576, 2575.35, 2578.5, 2580.2, 2582};
    vector<double> Hg_list = {2762, 2765, 2768, 2770, 2772, 2774, 2776};
    vector<double> slope_list;

    double interpolate(double y,double x1,double x2,double y1,double y2) {
        return (((x1 - x2)*(y - y2))/(y1 - y2)) + x2;
    }
    // double get_value(int x) {


    // }
    //Custom bound function - O(log(N))
    double fs1(int s,int e,double x) {
        if(s == e) return p_list[s] <= x ? s : -1;
        int mid = (s + e) >> 1;
        if(x < p_list[mid]) return fs1(s,mid,x);
        double ret = fs1(mid+1,e,x);
        return ret == -1 ? mid : ret;
    }

    double getjustlarge(double p) {
        double x = upper_bound(p_list.begin(),p_list.end(),p) - p_list.begin();
        return p_list[x];
    }

    double getjustsmall(double p) {
        double x = fs1(0,p_list.size(),p);
        return p_list[x];
    }

    double get_vf(double p) {
        double ans = 0;
        if(p < 0.7 || p > 1.00) ans = -1;
        if(p >= 0.7 && p <= 1.00) {
            auto it = find(p_list.begin(),p_list.end(),p);
            if(it != p_list.end()) {
                int idx = distance(p_list.begin(),it);
                ans = vf_list[idx];
            }else {
                double p1 = getjustlarge(p);
                double p2 = getjustsmall(p);

                int idx1 = find(p_list.begin(),p_list.end(),p1)-p_list.begin();
                int idx2 = find(p_list.begin(),p_list.end(),p2)-p_list.begin();

                double vf1 = vf_list[idx1];
                double vf2 = vf_list[idx2];
                ans = interpolate(p,vf1,vf2,p1,p2);
            }
        }
        return ans;
    }

    double get_vg(double p) {
        double ans = 0;
        if(p < 0.7 || p > 1.00) ans = -1;
        if(p >= 0.7 && p <= 1.00) {
            auto it = find(p_list.begin(),p_list.end(),p);
            if(it != p_list.end()) {
                int idx = distance(p_list.begin(),it);
                ans = vg_list[idx];
            }else {
                double p1 = getjustlarge(p);
                double p2 = getjustsmall(p);

                int idx1 = find(p_list.begin(),p_list.end(),p1)-p_list.begin();
                int idx2 = find(p_list.begin(),p_list.end(),p2)-p_list.begin();

                double vg1 = vg_list[idx1];
                double vg2 = vg_list[idx2];
                ans = interpolate(p,vg1,vg2,p1,p2);
            }
        }
        return ans;
    }

    double get_Uf(double p) {
        double ans = 0;
        if(p < 0.7 || p > 1.00) ans = -1;
        if(p >= 0.7 && p <= 1.00) {
            auto it = find(p_list.begin(),p_list.end(),p);
            if(it != p_list.end()) {
                int idx = distance(p_list.begin(),it);
                ans = Uf_list[idx];
            }else {
                double p1 = getjustlarge(p);
                double p2 = getjustsmall(p);

                int idx1 = find(p_list.begin(),p_list.end(),p1)-p_list.begin();
                int idx2 = find(p_list.begin(),p_list.end(),p2)-p_list.begin();

                double Uf1 = Uf_list[idx1];
                double Uf2 = Uf_list[idx2];
                ans = interpolate(p,Uf1,Uf2,p1,p2);
            }
        }
        return ans;
    }

    double get_Ug(double p) {
        double ans = 0;
        if(p < 0.7 || p > 1.00) ans = -1;
        if(p >= 0.7 && p <= 1.00) {
            auto it = find(p_list.begin(),p_list.end(),p);
            if(it != p_list.end()) {
                int idx = distance(p_list.begin(),it);
                ans = Ug_list[idx];
            }else {
                double p1 = getjustlarge(p);
                double p2 = getjustsmall(p);

                int idx1 = find(p_list.begin(),p_list.end(),p1)-p_list.begin();
                int idx2 = find(p_list.begin(),p_list.end(),p2)-p_list.begin();

                double Ug1 = Ug_list[idx1];
                double Ug2 = Ug_list[idx2];
                ans = interpolate(p,Ug1,Ug2,p1,p2);
            }
        }
        return ans;
    }

    double get_Hg(double p) {
        double ans = 0;
        if(p < 0.7 || p > 1.00) ans = -1;
        if(p >= 0.7 && p <= 1.00) {
            auto it = find(p_list.begin(),p_list.end(),p);
            if(it != p_list.end()) {
                int idx = distance(p_list.begin(),it);
                ans = Hg_list[idx];
            }else {
                double p1 = getjustlarge(p);
                double p2 = getjustsmall(p);

                int idx1 = find(p_list.begin(),p_list.end(),p1)-p_list.begin();
                int idx2 = find(p_list.begin(),p_list.end(),p2)-p_list.begin();

                double Hg1 = Hg_list[idx1];
                double Hg2 = Hg_list[idx2];
                ans = interpolate(p,Hg1,Hg2,p1,p2);
            }
        }
        return ans;
    }

    double g1(double p,double x) {
        double ans = 0;
        if(p < 0.7 || p > 1.00) ans = -1;
        if(p >= 0.7 && p <= 1.00) {
            double value1 = x/get_vf(p);
            double value2 = (1-x)/get_vg(p);
            ans = value1+value2;
        }
        return ans;
    }

    double g2(double p,double x) {
        double ans = 0;
        if(p < 0.7 || p > 1.00) ans = -1;
        if(p >= 0.7 && p <= 1.00) {
            double value1 = (x/get_vf(p)) * get_Uf(p);
            double value2 = ((1-x)/get_vg(p)) * get_Ug(p);
            ans = value1+value2;
        }
        return ans;
    }

    double solve(double pf,double pi,double xi,double xf) {
        double res = 0;
        if(pi < 0.7 || pi > 1.00) res = -1;
        if(pi >= 0.7 && pi <= 1.00) {
            res = (g2(pf,xf) - g2(pi,xi))/(g1(pf,xf) - g1(pi,xi));
        }
        return res;
    }

    double fs2(int s,int e,double x) {
        if(s == e) return slope_list[s] <= x ? s : -1;
        int mid = (s + e) >> 1;
        if(x < slope_list[mid]) return fs2(s,mid,x);
        double ret = fs2(mid+1,e,x);
        return ret == -1 ? mid : ret;
    }

    double getjustlargeslope(double trueslope) {
        double x = upper_bound(slope_list.begin(),slope_list.end(),trueslope) - slope_list.begin();
        return slope_list[x];
    }

    double getjustsmallslope(double trueslope) {
        double x = fs2(0,slope_list.size(),trueslope);
        return slope_list[x];
    }

    double get_Pf(double s,double pi,double xi,double xf) {
        double ans = -1;

        if(xi > xf) reverse(slope_list.begin(),slope_list.end());
        double s1 = getjustlargeslope(s);
        double s2 = getjustsmallslope(s);

        //deb(s1); deb(s2);

        auto it1 = find(slope_list.begin(),slope_list.end(),s1);
        auto it2 = find(slope_list.begin(),slope_list.end(),s2);

        if(it1 != slope_list.end() && it2 != slope_list.end()) {
            int idx1 = it1-slope_list.begin();
            int idx2 = it2-slope_list.begin();

            double p1 = (xi < xf)?p_list[idx1]:p_list[6-idx1];
            double p2 = (xi < xf)?p_list[idx2]:p_list[6-idx2];

            //deb(p1); deb(p2);

            while(abs(p1-p2) >= EPS) {
                double pMid = (p1 + p2)/2;
                assert(solve(pMid,pi,xi,xf) != -1);
                if(solve(pMid,pi,xi,xf) > s) {
                    p1 = pMid;
                }else {
                    p2 = pMid;
                }
            }
            ans = p1;
        }
        return ans;
    }
};

int main() {
    ios::sync_with_stdio(0); cin.tie(0);

#ifndef JUDGE
    //freopen("input.txt", "r", stdin);
    //freopen("output.txt", "w", stdout);
#endif

    double P1,P2,m1dot,m2dot,t,xi,xf,Pi;
    cin >> P1 >> P2 >> m1dot >> m2dot >> t >> Pi >> xi >> xf;
    P1 /= 1000; P2 /= 1000; Pi /= 1000;

    Solver s;
    for(auto x : s.p_list) s.slope_list.push_back(s.solve(x,Pi,xi,xf));

    double a = (m1dot - m2dot) * t;
    double b = (m1dot*s.get_Hg(P1) - m2dot*s.get_Hg(P2)) * t;
    //true slope
    double ts = b/a;
    double final_pressure = s.get_Pf(ts,Pi,xi,xf);
    double volume = a/(s.g1(final_pressure,xf) - s.g1(Pi,xi));

    //debug() << imie(s.slope_list);
    //debug() << imie(s.p_list);

    double check = (m1dot - m2dot) * (xi - xf);
    // main condition to handle all cases at once

    if(final_pressure != -1 && final_pressure > P2 && final_pressure < P1 && volume > 0 && check < 0) {
        assert(final_pressure > P2 && final_pressure < P1);
        final_pressure *= 1000;
        deb(final_pressure);
        deb(volume);
    }else {
        cout << "No Valid solution" << "\n";
    }
    return 0;
}
//971 711 8000 7000 0.4 800 0.94 0.95
