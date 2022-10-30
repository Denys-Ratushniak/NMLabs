#include <bits/stdc++.h>
#include <windows.h>
using namespace std;

typedef long double ld;
typedef long long ll;
typedef vector< vector<ld > > matrix;
ld eps;

ld k1,k2,k3,k4,k5;
ll type;
ld f(ld x)
{
    if(type == 1) return log(x) + (x + 1) * (x + 1) * (x + 1);
    return k1 * log(x) + k2 * pow((k3 * x + k4), k5);
}

ld derivative(ld x)
{
    if(type == 1) return (1.0 / x) + 3 * (x + 1) * (x + 1);
    return k1 / x + k2 * k3 * k5 * pow((k3 * x + k4), k5 - 1);
}

ld F(ld x, ld m)
{
    return (x - m * f(x));
}

ll totalans = 0;


void find_ans(ld a, ld b){
    ld m = 1.0 / max(derivative(a), derivative(b));
    ld x0 = (a + b) / 2.0;
    vector<ld> x;
    x.emplace_back(x0);

    ld y;
    ll steps = 0;

    while(steps <= 10000000)
    {
        y = F(x.back(), m);
        steps ++;
        if(fabs(y - x.back()) < eps) break;
        x.emplace_back(y);
    }
    cout << fixed << setprecision(10) << "Found in range [" << a << ", " << b << " ] answer with " << steps << " steps x = ";
    cout << fixed << setprecision(30) << x.back() << "\nf(x) = " << f(x.back()) << "\n";
    totalans++;
}

int main()
{
    cout << "Input 1 for f(x) = ln(x) + (x + 1) ^ 3 and other number for f(x) = k1 * ln(x) + k2(k3 * x + k4) ^ k5\n";
    cin >> type;
    if(type == 1)
    {
        /// f(x) = ln(x) + (x + 1) ^ 3
        ld a,b,x0;
        cout << "Here you can input multiple test to solve equation ln(x) + (x + 1)^3 = 0\n";
        while(1)
        {
            cout << "Input range [a, b], then starting x0 (a <= x0 <= b) and eps\n";
            cin >> a >> b >> x0 >> eps;
            if(a > b)
            {
                cout  << "a cannot be greater than b, try again\n";
                continue;
            }

            if(a <= 0)
            {
                cout  << "a cannot be less or equal than 0 due to ln(x), try again \n";
                continue;
            }

            if(x0 < a || b < x0)
            {
                cout  << "x0 must be in range [a, b] , try again\n";
                continue;
            }
            if(f(a) * f(b) > 0)
            {
                cout << "As function y = ln(x) + (x + 1)^3 is increasing and f(a) * f(b) > 0, there is no root in given range [a, b] , try again\n";
                continue;
            }

            ld m = 1.0 / max(derivative(a), derivative(b));

            vector<ld> x;
            x.emplace_back(x0);

            ld y;
            ll steps = 0;

            while(1)
            {
                y = F(x.back(), m);
                steps ++;
                if(fabs(y - x.back()) < eps) break;
                x.emplace_back(y);
            }
            cout << "Found answer with " << steps << " steps x = ";
            cout << fixed << setprecision(30) << x.back() << "\nf(x) = " << f(x.back()) << "\n";
        }


        return 0;
    }

    ///f(x) = k1 * ln(x) + k2(k3 * x + k4) ^ k5

    ld a,b,x0;
    ll total_seg = 0;
    cout << "Here you can input multiple test to solve equation k1 * ln(x) + k2(k3 * x + k4) ^ k5 = 0\n";
    while(1)
    {
        cout << "Input not zero k1, k2, k3, k4, k5\n";
        cin >> k1 >> k2 >> k3 >> k4 >> k5;
        totalans = 0;
        if(k1*k2*k3*k4*k5 == 0) {

            cout << "Not all k is zero, try again\n";
            continue;
        }

        cout << "\nInput range [a, b] where to find root(s), eps, and total segments for finding roots\n";
        cin >> a >> b >> eps >> total_seg;
        if(total_seg <= 0){
            cout  << "total segments must be greater than 0\n";
            continue;
        }
        if(a > b)
        {
            cout  << "a cannot be greater than b, try again\n";
            continue;
        }
        if(a <= 0)
        {
            cout  << "a cannot be less or equal than 0 due to ln(x), try again \n";
            continue;
        }

        ld len = (b - a) / total_seg;
        ld a1 = a;
        ld b1 = a1 + len;
        for(ll k = 0; k < total_seg; k++){
            if(f(a1) * f(b1) < 0) find_ans(a1, b1);
            a1 += len;
            b1 += len;
        }
        if(totalans == 0){
            cout << "No root found\n";
        }
    }

    /// 0.1 3 1 0.0001

    /// -100 1 1 1 3
    /// 0.1 5 0.00000000001 100

    return 0;
}
