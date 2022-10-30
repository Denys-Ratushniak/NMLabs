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

matrix get_random_matrix(ll size_, ll mxVal)
{
    ///randomize due to clock
    srand(clock());
    matrix res;
    res.resize(size_);
    for(int i = 0; i < size_; ++i)
    {
        res[i].resize(size_);
        for(int j = 0; j < size_; ++j) res[i][j] = (rand() - rand()/2)%mxVal;
    }
    return res;
}

ld linf_norm(matrix &A){
    ll n = A.size();
    ld norm = 0;
    for(int i = 0; i < n; ++ i){
        ld sum = 0;
        for(int j = 0; j < n; ++ j){
            sum += fabs(A[i][j]);
        }
        norm = max(norm, sum);
    }
    return norm;
}

ld l1_norm(matrix &A){
    ll n = A.size();
    ld norm = 0;
    for(int j = 0; j < n; ++ j){
        ld sum = 0;
        for(int i = 0; i < n; ++ i){
            sum += fabs(A[i][j]);
        }
        norm = max(norm, sum);
    }
    return norm;
}


ld eu_norm(matrix &A){
    ll n = A.size();
    ld sum = 0;
    for(int i = 0; i < n; ++ i){
        for(int j = 0; j < n; ++ j){
            sum += fabs(A[i][j]) * fabs(A[i][j]);
        }
    }
    return sqrt(sum);
}


void check_ans(matrix &A, vector<ld> &X, vector<ld> &B, ll cnt=0)
{
    answers << "\nHERE IS ANSWER(X): ";
    print(X, answers, cnt);
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


void solve(matrix &A, vector<ld> &B, ld allowed_error = 0.0000000001)
{
    answers << "Size = " << A.size() << "\n";
    answers << "Matrix A:\n";
    print(A, answers);
    answers << "B:\n";
    print(B, answers);

    ll n = A.size();
    matrix D(n, vector<ld>(n, 0));
    matrix LplusR(n, vector<ld>(n, 0));
    for(int i = 0; i < n; ++ i){
        for(int j = 0; j < n; ++ j){
            if(i == j) D[i][j] = -1.0/A[i][j];
            else LplusR[i][j] = A[i][j];
        }
    }
    matrix C(n, vector<ld>(n, 0));
    C = D * LplusR;
    ld n1 = l1_norm(C);
    ld n2 = linf_norm(C);
    ld n3 = eu_norm(C);
    answers << "Matrix C(-D^-1 * (L+R)):\n";
    print(C, answers, 5);
    answers << "\n Norms of C(-D^-1 * (L+R)):\n";
    answers << fixed << setprecision(15) << n1 << " " << n2 << " " << n3 << "\n";

    if(min({n1,n2,n3}) - 1 > eps){
        answers << "All norms of matrix C(-D^-1 * (L+R)) is > 1, so Jacobi method won`t converge\n";
        answers << fixed << setprecision(0) << "---------------------------------------------------------\n";
        return;
    }

    if(!check_DD(A)){
        answers << "This matrix is not diagonally dominant, so it is not garanteed for jacobi method to converge\n";
        answers << fixed << setprecision(0) << "---------------------------------------------------------\n";
        return;
    }


    answers << fixed << setprecision(15) << "CURRENT EPS:" << eps << "\n";
    answers << fixed << setprecision(0);
    vector<ld> x(n, 0.0);
    vector<ld> tmp_x(n, 0.0);

    ld error;
    ll step = 1;
    do
    {
        error = 0;
        tmp_x = B;
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < n; ++j)
                if (i != j) tmp_x[i] -= A[i][j] * x[j];

        answers << fixed << setprecision(10) << "Step: " << step << " Current X: ";

        for(int i = 0; i < n; ++ i)
        {
            ld x_upd = tmp_x[i] / A[i][i];
            ld e = fabs(x[i] - x_upd);
            answers << x_upd << " ";
            x[i] = x_upd;
            error = max(error, e);
        }
        answers << "\n";
        step++;
    }
    while (error > allowed_error);

    check_ans(A, x, B, 13);

    answers << fixed << setprecision(0) << "---------------------------------------------------------\n";
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
    A = get_random_matrix(sz, 1000000);
    vector<ld> B(sz);
    for(int i = 0; i < sz; ++i) B[i] = rand() - rand() / 2;
    for(int i = 0; i < sz; ++i) A[i][i] = 1000000000 + 2*i;
    //solve(A, B);
    return 0;
}
