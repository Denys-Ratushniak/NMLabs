#include <bits/stdc++.h>
#include <windows.h>
using namespace std;
typedef long double ld;
typedef long long ll;
typedef vector< vector<ld > > matrix;
const ld eps = 1e-6;

ifstream tests("input.txt");
ofstream answers("output.txt");

void print(vector<ld> &A, ll cnt = 3)
{
    cout << fixed << setprecision(cnt);
    for(int i = 0; i < A.size(); ++i) cout << A[i] << " ";
    cout << endl;
    cout << fixed << setprecision(0);
}

void print(vector<ld> &A, ofstream &answers, ll cnt = 3)
{
    answers << fixed << setprecision(cnt);
    for(int i = 0; i < A.size(); ++i) answers << A[i] << " ";
    answers << "\n\n";
    answers << fixed << setprecision(0);
}

void print(matrix &A, ofstream &answers, ll cnt=3)
{
    ll n = A.size();
    answers << fixed << setprecision(cnt);
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < A[i].size(); ++j) answers << A[i][j] << " ";
        answers << "\n";
    }
    answers << fixed << setprecision(0);
    answers << "\n";
}

void print(matrix &A, ll cnt=3)
{
    ll n = A.size();
    cout << fixed << setprecision(cnt);
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < A[i].size(); ++j) cout << A[i][j] << " ";
        cout << endl;
    }
    cout << fixed << setprecision(0);
    cout << endl;
}

matrix operator * (const matrix &A, const matrix &B)
{
    ll n = A.size();
    matrix res(n, vector<ld>(n, 0.0));
    for(int k = 0; k < n; ++k)
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < n; ++j)
                res[i][j] += A[i][k] * B[k][j];
    return res;
}

void matrix_mult_vector(matrix &A, vector<ld> &B, vector<ld> &res)
{
    ll n = A.size();
    res.resize(n, 0);
    for(int i = 0; i < n; ++ i)
    {
        res[i] = 0;
        for(int j = 0; j < n; ++ j)
            res[i] += A[i][j] * B[j];
    }

}


matrix get_random_matrix(ll size_)
{
    ///randomize due to clock
    srand(clock());
    matrix res;
    res.resize(size_);
    for(int i = 0; i < size_; ++i)
    {
        res[i].resize(size_);
        for(int j = 0; j < size_; ++j) res[i][j] = rand() - rand()/2;
    }
    return res;
}


void check_ans(matrix &A, vector<ld> &X, vector<ld> &B, ll cnt=0)
{
    answers << "HERE IS WHAT B SHOULD BE " << "\n";
    print(B, answers, cnt);
    answers << "HERE IS CALCULATED B" << "\n";
    ll n = A.size();
    vector<ld> ans(n,0);
    vector<ld> diff(n);
    ld euNorm = 0.0;
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            ans[i] += A[i][j] * X[j];
        }
        diff[i] = B[i] - ans[i];
        euNorm += diff[i] * diff[i];
    }
    print(ans, answers, cnt);
    answers << "HERE is B - B(calculated)" << "\n";
    print(diff, answers, cnt);
    answers << "Norm of this vector = " << fixed << setprecision(cnt) << sqrt(euNorm) << "\n";
    answers << fixed << setprecision(0);
}


ld eu_norm_vector(vector<ld> &A)
{
    ld res = 0;
    for(auto to: A) res += to * to;
    return sqrt(res);
}


ld eu_norm_matrix(matrix &A)
{
    ld res = 0;
    for(int i = 0; i < A.size(); ++ i)
        for(int j = 0; j < A[i].size(); ++j)
            res += A[i][j] * A[i][j];
    return sqrt(res);
}


bool is_determinant_zero(matrix &A)
{
    ll n = A.size();
    matrix AE;
    AE.resize(n);

    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++ j)
            AE[i].emplace_back(A[i][j]);

    for(int k = 0; k < n; ++k)
    {
        int row_max_el = k;

        /// finding maximum element
        for(int i = k + 1; i < n; ++i)
            if(fabs(AE[i][k]) > fabs(AE[row_max_el][k])) row_max_el = i;

        if(fabs(AE[row_max_el][k]) < eps) return 0;


        /// change of 2 rows for optimization of division by the maximum element
        if(k != row_max_el)
            for(int j = k; j < n; ++j) swap(AE[k][j], AE[row_max_el][j]);


        /// forming zeros under main diagonal
        for(int i = k + 1; i < n; ++i)
        {
            ld mik = -AE[i][k]/AE[k][k];

            AE[i][k] = 0;
            for(int j = k + 1; j < n; ++j) AE[i][j] += mik * AE[k][j];
        }
    }
    for(int i = 0; i < n; ++i)
    {
        if(fabs(AE[i][i]) < eps) return 0;
    }
    return 1;

}

void solve(matrix &A)
{
    ll n = A.size();
    answers << "Size = " << A.size() << "\n";
    answers << "Matrix A:\n";
    print(A, answers);

    vector<ld> x[2], y[2], lambda[2];
    ld y_norm;
    x[0].resize(n);
    x[1].resize(n);
    y[0].resize(n);
    y[1].resize(n);
    lambda[0].resize(n);
    lambda[1].resize(n);
    vector<bool> was[2];
    was[0].resize(n);
    was[1].resize(n);

    for(int i = 0; i < n; ++ i) y[0][i] = i + 1;
    y_norm = eu_norm_vector(y[0]);

    answers << "Vector y0:\n";
    print(y[0], answers, 7);

    answers << fixed << setprecision(10) << "Norm of y0 = " << y_norm << "\n";

    for(int i = 0; i < n; ++ i)
    {
        x[0][i] = y[0][i] / y_norm;
        if(fabs(x[0][i]) > eps) was[0][i] = 1;
        else was[0][i] = 0;
    }

    ld ANS = 0;

    answers << "\n";

    vector<ld> ansV(n);
    ll step;
    for(step = 1; ; ++step)
    {
        matrix_mult_vector(A, x[!(step&1)], y[step&1]);

        answers << "STEP = " << step << " \n";


        y_norm = eu_norm_vector(y[step&1]);

        answers << "Vector y" << step << ":\n";
        print(y[step&1], answers, 7);
        answers << fixed << setprecision(10) << "Norm of y" << step << " = " << y_norm << "\n";

        for(int i = 0; i < n; ++ i)
        {
            x[step&1][i] = y[step&1][i] / y_norm;
            if(step == 1)
            {
                if(fabs(x[step&1][i]) > eps) was[1][i] = 1;
                else was[1][0] = 0;
            }
        }



        answers << "Vector x" << step << ":\n" ;
        print(x[step&1], answers, 7);



        for(int i = 0; i < n; ++ i)
        {
            if(fabs(x[!(step&1)][i]) > eps ) lambda[step&1][i] = y[step&1][i] / x[!(step&1)][i];
        }

        answers << "Lambda" << step << "Vector:\n" ;
        print(lambda[step&1], answers, 7);
        answers << "\n";

        if(step > 1)
        {
            ld sum = 0;
            ll total = 0;
            bool good = true;
            for(int i = 0; i < n; ++ i)
            {
                if(was[0][i] & was[1][i])
                {
                    total ++;
                    sum += lambda[step&1][i];
                    if(fabs(lambda[0][i] - lambda[1][i]) >= eps) good = false;
                }
            }

            if(good)
            {
                ANS = sum / ld(total);
                for(int i = 0; i < n; ++ i) ansV[i] = x[step&1][i];
                break;
            }

            for(int i = 0; i < n; ++ i) was[0][i] = was[1][i];
            for(int i = 0; i < n; ++ i)
            {
                if(fabs(x[step&1][i]) > eps) was[1][i] = 1;
                else was[1][0] = 0;

            }
        }
    }

    answers << fixed << setprecision(10) << "Biggest Eigen Value = " << ANS << "\n";
    answers << "Eigen Vector:\n";
    print(ansV, answers, 7);

    for(int i = 0; i < n; ++ i) A[i][i] -= ANS;
    answers << "TOTAL STEPS = " << step << "\n";
    if(is_determinant_zero(A)) answers << "Determinant of (A - lambda * E) matrix -> 0 so found Eigen Value is found correctly\n";
    //cout << step << "\n";

    answers << "---------------------------------------------------------\n";
}
int main()
{
    ll n;
    while(tests >> n)
    {
        matrix A(n, vector<ld>(n, 0));
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < n; ++j)
                tests >> A[i][j];

        solve(A);
    }
    ll sz = 500;
    matrix A(sz, vector<ld>(sz, 0));
    A = get_random_matrix(sz);
    solve(A);
    return 0;
}
