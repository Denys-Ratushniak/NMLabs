#include <bits/stdc++.h>
using namespace std;

typedef long double ld;
typedef long long ll;
typedef vector< vector<ld > > matrix;

ld eps;

int main()
{
    ld x,y;
    cout << "Input first guess x and y\n";
    cin >> x >> y;
    vector<ld> X,Y;
    X.push_back(x);
    Y.push_back(y);
    eps = 1e-2;
    ll steps = 0;
    while(1){
        ld oldx = X.back();
        ld oldy = Y.back();
        steps++;
        x = 3 + cos(oldy);
        y = sqrt(-0.8 - cos(oldx - 1));

        X.push_back(x);
        Y.push_back(y);

        if(steps > 2 && steps % 2 == 1 && fabs(x - X[X.size()-3]) < eps) break;
        if(steps > 2 && steps % 2 == 0 && fabs(y - Y[Y.size()-3]) < eps) break;
    }
    if(ll(X.size()) % 2 == 0){
        x = X.back();
        y = acos(x - 3);
    }
    else {
        y = Y.back();
        x = acos(-0.8 - y * y)+ 1;
    }
    cout << "SEIDEL METHOD TO SOLVE NONLINEAR SOE - > x1 = " << fixed << setprecision(5) << x << " y1 = " << y << "\n";
    cout << "x2 = " << fixed << setprecision(5) << x << " y2 = " << -y << "\n";
    cout << "First f(x1, y1) = 0 -> " << fixed << setprecision(5) << (cos(x - 1) + y * y + 0.8) << "\n";
    cout << "Second f(x1, y1) = 0 -> " << fixed << setprecision(5) << (x - cos(y) - 3) << "\n";
    cout << "First f(x2, y2) = 0 -> " << fixed << setprecision(5) << (cos(x - 1) + y * y + 0.8) << "\n";
    cout << "Second f(x2, y2) = 0 -> " << fixed << setprecision(5) << (x - cos(-y) - 3) << "\n";
    return 0;
}
