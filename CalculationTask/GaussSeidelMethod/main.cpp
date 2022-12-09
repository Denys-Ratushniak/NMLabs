#include <bits/stdc++.h>
using namespace std;

typedef long double ld;
typedef long long ll;
typedef vector< vector<ld > > matrix;

ld eps;
void print(vector<ld> &A, ll cnt = 5)
{
    cout << fixed << setprecision(cnt);
    for(int i = 0; i < A.size(); ++i) cout << A[i] << " ";
    cout << "\n";
    cout << fixed << setprecision(0);
}

void print(matrix &A, ll cnt=5)
{
    ll n = A.size();
    cout << fixed << setprecision(cnt);
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < A[i].size(); ++j) cout << A[i][j] << " ";
        cout << "\n";
    }
    cout << fixed << setprecision(0);
    cout << "\n";
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

ld euclidnorm(vector<ld> &guess, vector<ld> &ans){
    ld norm = 0;
    for(int i = 0; i < guess.size(); ++i)
        norm += (guess[i] - ans[i]) * (guess[i] - ans[i]);

    return sqrt(norm);
}


void seidel(matrix &A, vector<ld> &B)
{
    ll n = A.size();
    vector<ld> Y(n, 0);
    vector<ld> x(n, 0);
    x[0] = 1;
    x[1] = 1;
    x[2] = 1;
    vector<ld> y(n, 0);
    vector<ld> ybefor(n, 0);
    ll steps = 0;

    while(1)
    {
        for(int i = 0; i < n; ++i) ybefor[i] = y[i];
        for(int i = 0; i < n; i++)
        {
            y[i] = (B[i] / A[i][i]);
            for(int j = 0; j < n; j++)
            {
                if(j == i) continue;
                y[i] = y[i] - ((A[i][j] / A[i][i]) * x[j]);
                x[i] = y[i];
            }

        }

        ld norm = euclidnorm(y, ybefor);
        steps ++;

        if(steps % 1000000 <= 2) {
            cout << "Step: " << steps << setprecision(5) << " norm = "<< norm << " current answer: \n";
            print(y);
        }
        //cout << "\n";

        if(norm < eps) break;
    }
    cout << "SEIDEL METHOD ->  x = "; print(y);
    cout << "Iterations: " << steps << "\n";
}

bool is_determinant_zero(matrix A)
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

        /// if maximum element goes to 0 -> inverse matrix does not exist
        if(fabs(AE[row_max_el][k]) < eps) return true;


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
        if(fabs(AE[i][i]) < eps) return true;
    }
    return false;

}

bool check_valid(matrix &A, ll &num)
{
    matrix B;
    ll n = A.size();
    B.resize(1);
    B[0].resize(1);
    B[0][0] = A[0][0];
    vector<ld> add(1);
    num = 0;
    if(fabs(A[0][0]) < eps) return false;
    for(int i = 1; i < n; ++i)
    {
        for(int num = 0; num < i; ++num) B[num].emplace_back(A[num][i]);
        for(int j = 0; j < i; ++j) add[j] = A[i][j];
        add.push_back(A[i][i]);
        B.push_back(add);
        if(is_determinant_zero(B)) {
            num = i;
            return false;
        }
    }
    return true;
}

void decompose_LU(matrix &A, matrix &L, matrix &U)
{
    ll n = A.size();
    vector<ld> null(n, 0);
    L.resize(n, vector<ld>(n, 0));
    U.resize(n, vector<ld>(n, 0));
    for(int s = 0; s < n; s++)
    {
        for(int j = s; j < n; j++)
        {
            U[s][j] = A[s][j];
            for(int k = 0; k < s; k++) U[s][j] -= L[s][k] * U[k][j];
        }
        for(int i = s; i < n; i++)
        {
            L[i][s] = A[i][s];
            for(int k = 0; k < s; k++) L[i][s] -= (L[i][k] * U[k][s]);
            L[i][s] /= U[s][s];
        }
    }
}


void LY_eq_B(matrix &L, vector<ld> &Y, vector<ld> &B)
{
    ll n = L.size();
    for(int i = 0; i < n; ++ i)
    {
        ld h = 0;
        for(int j = 0; j < i; ++j)
        {
            h += L[i][j] * Y[j];
        }
        Y[i] = B[i] - h;
    }
}

void UX_eq_Y(matrix &U, vector<ld> &X, vector<ld> &Y)
{
    ll n = U.size();
    for(int i = n - 1; i >= 0; -- i)
    {
        ld h = 0;
        for(int j = n - 1; j > i; --j)
        {
            h += U[i][j] * X[j];
        }
        X[i] = (Y[i] - h)/U[i][i];
    }
}


void LU(matrix &A, vector<ld> &B)
{
    ll n = A.size();
    ll num;
    if(!check_valid(A, num))
    {
        cout << "Cannot do LU decomposition, because " << num << "th leading principal minors is equal to 0\n";
        return;
    }

    matrix L,U;
    decompose_LU(A, L, U);
    cout << "L:\n" ;print(L);
    cout << "U:\n" ;print(U);
    vector<ld> Y(n);
    vector<ld> X(n);
    LY_eq_B(L, Y, B);
    UX_eq_Y(U, X, Y);
    cout << "LU METHOD -> x = ";
    print(X);

}

bool check_DD(matrix &A){
    ll n = A.size();
    bool if_one = 0;
    for(int i = 0; i < n; ++ i){
        ld sum = 0;
        for(int j = 0; j < n; ++ j){
            sum += fabs(A[i][j]);
        }
        sum -= fabs(A[i][i]);
        if(fabs(A[i][i]) > sum) if_one = 1;
        if(fabs(A[i][i]) < sum) return false;
    }
    return if_one;
}

int main()
{
    eps = 1e-2;
    matrix A;
    vector<ld> B,X;
    //A = {{5, 2, 4}, {-1, 3, -6}, {7, -1, -3}};
    //B = {9, 8, 7};

    A = {{10, 2, 2}, {1, 15, 1}, {1, 1, 23}};
    B = {1, 1, 1};

    cout << "A:\n";print(A);
    cout << "B:\n";print(A);

    LU(A, B);
    if(check_DD(A)) seidel(A, B);
    else {
        cout << "Matrix isnt diagonally dominant so Seidel Method may not coverge" << "\n";
        return 0;
    }
    return 0;
}

/*
 5  2  4  9
-1  3 -6  8
 7 -1 -3  7

 12 1 7 16
-1  3 -6  8
 7 -1 -3  7

 12  1   7  16
-1   3  -6   8
 0 20  -45   63*/



