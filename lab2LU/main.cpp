#include <bits/stdc++.h>
using namespace std;
typedef long double ld;
typedef long long ll;
typedef vector< vector<ld > > matrix;
const ld eps = 1e-10;

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


void solve(matrix &A, vector<ld> &B)
{
    answers << "Size = " << A.size() << "\n";
    answers << "Matrix A:\n";
    print(A, answers);
    answers << "B:\n";
    print(B, answers);
    ll num;
    if(!check_valid(A, num))
    {
        answers << "Cannot do LU decomposition, because " << num << "th leading principal minors is equal to 0\n";
        answers << "---------------------------------------------------------\n";
        return;
    }

    matrix L,U;
    decompose_LU(A, L, U);
    answers << "Matrix L:\n";
    print(L, answers, 5);
    answers << "Matrix U:\n";
    print(U, answers, 5);
    matrix mult = L * U;
    answers << "Matrix mult(L * U):\n";
    print(mult, answers, 5);
    ld euNorm = 0;
    answers << "Matrix A - mult(L * U):\n";
    answers << fixed << setprecision(22);
    ll n = A.size();
    for(int i = 0; i < n; ++ i)
    {
        for(int j = 0; j < n; ++j)
        {
            euNorm += (A[i][j] - mult[i][j]) * (A[i][j] - mult[i][j]);
            answers << A[i][j] - mult[i][j] << " ";
        }
        answers << "\n";
    }
    answers << fixed << setprecision(22) << "\nNorm of matrix A - mult = " << sqrt(euNorm) << "\n";
    answers << fixed << setprecision(0);
    vector<ld> Y(n);
    vector<ld> X(n);
    LY_eq_B(L, Y, B);
    UX_eq_Y(U, X, Y);
    answers << "Y: ";
    print(Y, answers, 22);
    answers << "X: ";
    print(X, answers, 22);
    check_ans(A, X, B, 22);
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
        vector<ld> B(n);
        for(int i = 0 ; i < n; ++i) tests >> B[i];
        solve(A, B);
    }
    ll sz = 300;
    matrix A(sz, vector<ld>(sz, 0));
    A = get_random_matrix(sz);
    vector<ld> B(sz);
    for(int i = 0; i < sz; ++i) B[i] = rand() - rand() / 2;
    //solve(A, B);
    return 0;
}
