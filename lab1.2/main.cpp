#include <bits/stdc++.h>
using namespace std;
typedef long double ld;
typedef long long ll;
typedef vector< vector<ld > > matrix;
const ld eps = 1e-10;
ifstream tests("input.txt");
ofstream answers("output.txt");
void print(matrix &A)
{
    ll n = A.size();
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < A[i].size(); ++j) answers << A[i][j] << " ";
        answers << "\n";
    }
    answers << "\n";
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

pair<bool, matrix> get_inverse_matrix(matrix A)
{
    ll n = A.size();
    matrix AE;
    AE.resize(n);
    /// forming matrix A|E
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++ j) AE[i].emplace_back(A[i][j]);

        for(int j = 0; j < n; ++ j) AE[i].emplace_back(i == j);
    }
    for(int k = 0; k < n; ++k)
    {
        int row_max_el = k;

        /// finding maximum element
        for(int i = k + 1; i < n; ++i)
            if(fabs(AE[i][k]) > fabs(AE[row_max_el][k])) row_max_el = i;

        /// if maximum element goes to 0 -> inverse matrix does not exist
        if(fabs(AE[row_max_el][k]) < eps) return {false, A};


        /// change of 2 rows for optimization of division by the maximum element
        if(k != row_max_el)
            for(int j = k; j < 2 * n; ++j) swap(AE[k][j], AE[row_max_el][j]);


        /// forming zeros under main diagonal
        for(int i = k + 1; i < n; ++i){
            ld mik = -AE[i][k]/AE[k][k];

            AE[i][k] = 0;
            for(int j = k + 1; j < 2 * n; ++j) AE[i][j] += mik * AE[k][j];
        }
    }

    for(int k = n - 1; k >= 0; --k){

        /// making 1 on the main diagonal
        for(int j = n; j < 2 * n; ++ j) AE[k][j] /= AE[k][k];
        AE[k][k] = 1;

        /// making 0 above main diagonal
        for(int i = k - 1; i >= 0; --i){
            ld mn = -AE[i][k];
            AE[i][k] = 0;
            for(int j = n; j < 2 * n; ++j)
            {
                AE[i][j] += AE[k][j] * mn;
            }
        }
    }

    ///forming inverse matrix
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            A[i][j] = AE[i][j + n];

    return {true, A};
}

matrix get_random_matrix(ll size_)
{
    ///randomize due to clock
    srand(clock());
    matrix res;
    res.resize(size_);
    for(int i = 0; i < size_; ++i){
        res[i].resize(size_);
        for(int j = 0; j < size_; ++j) res[i][j] = rand() - rand()/2;
    }
    return res;
}

void solve(matrix &A)
{
    ll n = A.size();
    pair<bool, matrix> inv_A = get_inverse_matrix(A);
    answers << "Our matrix:\n";
    print(A);
    if(inv_A.first)
    {
        answers << "Inverse matrix:\n";
        print(inv_A.second);

        answers << "The product of the given matrix and its inverse:\n";
        matrix E = A * inv_A.second;
        print(E);

        ld sum = 0.0;
        /// Calculating deviation - Frobenius norm
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < n; ++j)
                sum += (E[i][j] - (i == j)) * (E[i][j] - (i == j));

        answers << "Deviation from the unit matrix = Frobenius norm = " << sqrt(sum) << "\n";

        ld norm_A = 0, norm_invA = 0;
        ld sum_A, sum_invA;
        /// Calculating norms and cond
        for(int i = 0; i < n; ++i)
        {
            sum_A = 0;
            sum_invA = 0;

            for(int j = 0; j < n; ++j)
            {
                sum_A += fabs(A[i][j]);
                sum_invA += fabs(inv_A.second[i][j]);
            }

            norm_A = max(norm_A, sum_A);
            norm_invA = max(norm_invA, sum_invA);
        }

        answers << "Norm of the given matrix = " << norm_A << "\n";
        answers << "Norm of the inverse matrix = " << norm_invA << "\n";
        answers << "Cond of the given matrix = " << norm_A * norm_invA << "\n";
    }
    else answers << "The inverse matrix cannot be found, since the determinant of this matrix is zero\n";
    answers << "-------------------------------------------------------\n\n";
}

int main()
{
    ll n;
    while(tests >> n){
        matrix A(n, vector<ld>(n, 0));
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < n; ++j)
                tests >> A[i][j];
        solve(A);
    }
    ll cnt_of_random_matrix = 1;
    while(cnt_of_random_matrix--){
        matrix A = get_random_matrix(300);
        solve(A);
    }
    return 0;
}
