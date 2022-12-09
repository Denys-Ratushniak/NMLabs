#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
typedef long double ld;
const ld eps = 0.5e-6;

ld fx1(ld x){
    return  (10 * (x * x + 3.5)) / (0.65 + 0.17 * x * x);
}


ld fx2(ld x){
    return  (x * x) * sqrt(1 + x) / sqrt(3 * x + 13.5);
}


ld S(ll type, ld left, ld right, ld n){
    ld h = (right - left) / ld(n);
    ld sumy = 0;
    ld y1yn = 0;
    ld x = left;
    ld now = 0;
    for(int i = 1; i <= n; ++ i){
        now = 0;
        if(type == 1) now = fx1(x);
        else now = fx2(x);

        if(i == 1 || i == n) y1yn += now;
        else sumy += now;

        x += h;
    }
    return h * (y1yn / 2.0 + sumy);
}

int main()
{
    ll type;
    cout << "Enter 1 for integral of (10 * (x * x + 3.5)) / (0.65 + 0.17 * x * x) from 2 to 6\n";
    cout << "Enter 2 for integral of (x * x) * sqrt(1 + x) / sqrt(3 * x + 13.5) from 1 to 2\n";
    cin >> type;
    if(type != 1 && type != 2){
        cout << "Error, enter 1 or 2\n";
        return 0;
    }
    ll n0;
    if(type == 1) cout << "Enter number of points n0 for first step h = (b - a)/n0 to calculate first integral\n";
    else cout << "Enter number of points n0 for first step h = (b - a)/n0 to calculate second integral\n";

    cin >> n0;
    if(n0 <= 1) {
        cout << "Number of points must be greater than 1\n";
        return 0;
    }
    ld t1 = clock();
    ld left;
    ld right;
    ld h;
    if(type == 1) {
        left = 2;
        right = 6;
    }
    else {
        left = 1;
        right = 3;
    }
    ld RESULT = 0;
    ll step = 1;
    cout << "Step/number of points/ number of points * 2, result for n, result for 2n, error" << endl;
    while(1){
        ll n2 = n0 * 2;
        ld ans1 = S(type, left, right, n0);
        ld ans2 = S(type, left, right, n2);

        cout << fixed << setprecision(10) << "STEP " << step << ": " << n0 << " " << n2 << " " << ans1 << " " << ans2 << " " << fabs(ans1 - ans2) << endl;
        n0 = n2;
        step++;
        if(fabs(ans1 - ans2) < eps){
            RESULT = ans2;
            break;
        }

    }
    ld t2 = clock();
    cout << fixed << setprecision(10) << RESULT << " Time: " << (t2 - t1)/CLOCKS_PER_SEC << " seconds" << endl;
    return 0;
}
