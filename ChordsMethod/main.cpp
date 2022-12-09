#include <bits/stdc++.h>
using namespace std;

typedef long double ld;
typedef long long ll;
typedef vector< vector<ld > > matrix;

ld eps;
ld a,b,total_seg;

ld f(ld x)
{
    return pow(x, -x) - x * x;
}


ld deriv(ld x, ld eps)
{
    ld h = 0.5;
    ld d0 = 0;
    ld d = (f(x + h) - f(x - h)) / (2 * h);
    while(fabs(d - d0) > eps)
    {
        d0 = d;
        h *= 0.5;
        d = (f(x + h) - f(x - h)) / (2 * h);
    }
    return d;
}

ld secondderiv(ld x, ld eps)
{
    ld h = 0.5;
    ld d0 = 0;
    ld d = (f(x + h) - 2 * f(x) + f(x - h)) / (h * h);
    while(fabs(d - d0) > eps)
    {
        d0 = d;
        h *= 0.5;
        d = (f(x + h) - 2 * f(x) + f(x - h)) / (h * h);
    }
    return d;
}


void chords(ld a, ld b)
{
    ll steps = 0;
	ld x = 0;

	if(deriv(a, eps) * secondderiv(b, eps) < 0) swap(a, b);

	while(1)
	{
		x = a - (f(a) * (b - a)) / (f(b) - f(a));
		steps++;
		if (abs(a - x) < eps)
			break;
		a = x;
	}

	cout << fixed << setprecision(10) << "CHORDS METHOD ->  x = " << x << " Iterations: " << steps << "\n";
	return;
}


void tangents(ld a, ld b){

    ll steps = 0;
    ld x;
    if(deriv(a, eps) * secondderiv(a, eps) > 0) x = a;
    else x= b;
    do {
        x = x - f(x) / deriv(x, eps);
        steps += 1;
    }
    while(fabs(f(x)) >= eps);

    cout << fixed << setprecision(10) << "TANGENTS METHOD ->  x = " << x << " Iterations: " << steps << "\n";
}



bool input()
{
    cin >> a >> b >> total_seg;
    eps = 1e-5;
    if(total_seg <= 0)
    {
        cout  << "total segments must be greater than 0\n";
        return 0;
    }
    if(a > b)
    {
        cout  << "a cannot be greater than b, try again\n";
        return 0;
    }
    if(a <= 0)
    {
        cout  << "a cannot be less or equal than 0 due to x^(-x), try again \n";
        return 0;
    }
    return 1;
}

int main()
{
    while(!input());

    ld len = (b - a) / total_seg;
    ld a1 = a;
    ld b1 = a1 + len;
    ll total = 0;
    for(ll k = 0; k < total_seg; k++)
    {
        if(f(a1) * f(b1) < 0) {
            chords(a1, b1);
            tangents(a1, b1);
            total++;
        }
        a1 += len;
        b1 += len;
    }
    if(total == 0) cout << "No roots were found\n";

    return 0;
}
