#include <bits/stdc++.h>
#include <sys/time.h>
#include <omp.h>
#define Vf vector<vector<double>>
#define epsilon 0.25

#define MAX_THREAD 128

using namespace std;

int num = 32;
template <typename T>
double sgn(T val)
{
        return (val > T(0)) - (val < T(0));
}

void print_array(int *I1, int *I2, int rows, int columns)
{
        printf("Array I1");
        cout << endl;
        for (int i = 0; i < rows; i++)
        {
                cout << I1[i] << " ";
        }
        cout << endl;
        cout << "Array I2" << endl;
        for (int j = 0; j < columns; j++)
        {
                cout << I2[j] << " ";
        }
        cout << endl;
}
void checking_input_parameters(int argc, int rows, int columns)
{
        if (argc < 2)
        {
                cout << "Invalid input" << endl;
                exit(0);
        }
        if (rows != columns)
        {
                cout << "invalid matrix" << endl;
                exit(0);
        }
}
void fill_ones_in_diagonal(vector<vector<double>> &arr)
{
        for (int i = 0; i < arr.size(); i++)
        {
                for (int j = 0; j < arr[i].size(); j++)
                {
                        if (i == j)
                        {
                                arr[i][j] = 1;
                        }
                        else
                        {
                                arr[i][j] = 0;
                        }
                }
        }
}
void assiging_sin_cos_operations_to_matrix(Vf &matrix_U, Vf &matrix_V, double cos_theta, double sin_theta, double tau, int i, int j, int columns)
{
        int k = 0 ;
        while(k<columns)
        {       
                double tau2;
                tau = matrix_V[i][k];
                tau2 = tau;
                matrix_V[i][k] = cos_theta * tau2 - sin_theta * matrix_V[j][k];
                matrix_V[j][k] = sin_theta * tau2 + cos_theta * matrix_V[j][k];

                double tau1 ;
                tau1 = tau;
                tau1 = matrix_U[i][k];
                matrix_U[i][k] = tau1*cos_theta - matrix_U[j][k]*sin_theta;
                matrix_U[j][k] = tau1*sin_theta + matrix_U[j][k]*cos_theta;


                k++;
        }
}

void calculate_alpha_beta_gamma(double &alpha, double &beta, double &gamma, Vf &matrix_U, int columns, int i, int j)
{
        for (int k = 0; k < columns; k++)
        {
                alpha += (matrix_U[i][k] * matrix_U[i][k]);
                gamma += (matrix_U[i][k] * matrix_U[j][k]);
                beta += (matrix_U[j][k] * matrix_U[j][k]);
        }
}
pair<double, double> schur(vector<vector<double>>& A, int p, int q)
{
    pair<double, double> cs ;
    
    double app = 0.0 ;
    double aqq = 0.0 ;
    double apq = 0.0 ;

    double t , tau ;

    for(int k = 0 ; k < A.size() ; k++)
    {
        app+= A[p][k] * A[p][k] ;
        aqq+= A[q][k] * A[q][k] ;
        apq+= A[p][k] * A[q][k] ;
    }

    if(apq != 0)
    {
        tau = (app - aqq) / (2* apq) ;
        if(tau >= 0)
            t = 1/(tau + (sqrt(1 + pow(tau,2)))) ;
        else
            t = -1/(tau + (sqrt(1 + pow(tau,2)))) ;

        double c = 1/(1+pow(t,2)) ;
        double s = c * t ;

        cs = make_pair(c,s) ;
       

        
    }

    else
    {
        cs = make_pair(1.0,0.0) ;
    }

    return cs ;
}

void parallel_svd(int rows, int columns, vector<vector<double>> &A, vector<vector<double>> &U, vector<double> &S, vector<vector<double>> &V)
{

        vector<vector<double>> matrix_U(columns + 1, vector<double>(columns, 0));
        vector<vector<double>> matrix_V(columns + 1, vector<double>(columns, 0));

        double alpha;
        double beta;
        double gamma;
        double cos_theta;
        double zeta;
        double tau;
        double sin_theta;
        double sub_zeta;
        double converge = 1.0;
        int acum = 0;
        double temp1, temp2;
        double time_1, time_2;
        vector<int> I1(A.size() + 1, 0);
        vector<int> I2(A.size() + 1, 0);
        timeval start, end, end2;
        vector<double> C(MAX_THREAD);

        matrix_U = A;
        fill_ones_in_diagonal(matrix_V);

        double conv;
        while (true)
        {
                if (epsilon >= converge || acum == 20)
                {
                        break;
                }
                else
                {
                        acum += 1;
                }
                converge = 0.0;
                for (int l = 1; l < A.size(); l++)
                {
                        int r1 = 0;
                        int r2 = 0;
                        for (int i = 0; i + l < A.size(); i++)
                        {
                                if (i % (2 * l) < l)
                                {
                                        r1 += 1;
                                        I1[r1] = i;
                                }
                                else
                                {
                                        r2 += 1;
                                        I2[r2] = i;
                                }
                        }

                        for (int k = 0; k < num; k++)
                        {
                                C[k] = converge;
                        }

#pragma omp parallel for num_threads(num)
                        for (int p = 1; p <= r1; p++)
                        {
                                int k = omp_get_thread_num();
                                int i = I1[p];
                                int j = i + l;
                                double alpha = 0;
                                double beta = 0;
                                double gamma = 0;
                                double zeta, tau, cos_theta, sin_theta;
                                calculate_alpha_beta_gamma(alpha, beta, gamma, matrix_U, A.size() ,  i,  j);
                                double temp = max(C[k], abs(gamma) / sqrt(alpha * beta));
                                C[k] = temp;
                                zeta = (beta - alpha) / (2.0 * gamma);
                                tau = sgn(zeta) / (abs(zeta) + sqrt(1.0 + (zeta * zeta)));
                                cos_theta = 1.0 / (sqrt(1.0 + (tau * tau)));
                                sin_theta = cos_theta * tau;
                                assiging_sin_cos_operations_to_matrix(matrix_U, matrix_V, cos_theta, sin_theta, tau, i, j, columns);
                        }
#pragma omp parallel for num_threads(num)
                        for (int p = 1; p <= r2; p++)
                        {
                                int k = omp_get_thread_num();
                                int i = I2[p], j = i + l;
                                double alpha = 0, beta = 0, gamma = 0;
                                double zeta, tau, cos_theta, sin_theta;
                                calculate_alpha_beta_gamma(alpha, beta, gamma, matrix_U, A.size() ,  i,  j);
                                double temp = max(C[k], abs(gamma) / sqrt(alpha * beta));
                                C[k] = temp;
                                zeta = (beta - alpha) / (2.0 * gamma);
                                tau = sgn(zeta) / (abs(zeta) + sqrt(1.0 + (zeta * zeta)));

                                cos_theta = 1.0 / (sqrt(1.0 + (tau * tau)));
                                sin_theta = cos_theta * tau;
                                assiging_sin_cos_operations_to_matrix(matrix_U, matrix_V, cos_theta, sin_theta, tau, i, j, columns);
                        }

                        for (int k = 0; k < num; k++)
                        {
                                double temp2 = max(converge,C[k]);
                                converge = temp2;
                        }
                               
                }
        }

        for (int i = 0; i < rows; i++)
        {

                tau = 0;
                for (int j = 0; j < columns; j++)
                {
                        tau = tau + pow(matrix_U[i][j], 2);
                }
                tau = sqrt(tau);

                for (int j = 0; j < columns; j++)
                {
                        matrix_U[i][j] = matrix_U[i][j] / tau;
                        if (i == j)
                        {
                                S[i] = tau;
                        }
                }
        }

        for (int i = 0; i < rows; i++)
        {

                for (int j = 0; j < columns; j++)
                {

                        U[i][j] = matrix_U[j][i];
                        V[i][j] = matrix_V[j][i];
                }
        }
}

Vf transpose_matrix(Vf a)
{
        Vf b(a[0].size(), vector<double>(a.size(), 0));
        for (int i = 0; i < a.size(); i++)
        {
                for (int j = 0; j < a[i].size(); j++)
                {
                        b[j][i] = a[i][j];
                }
        }
        return b;
}

Vf dense_matrix_mul(Vf &a, Vf &b)
{
        int m = a.size();
        int n = a[0].size();
        int l = b.size();
        int k = b[0].size();

        if (n != l)
        {
                cout << "Invalid matrix multiplication" << endl;
                exit(1);
        }
        else
        {
                Vf result;
                vector<double> temp(k, 0);
                for (int i = 0; i < m; i++)
                        result.push_back(temp);

                for (int i = 0; i < m; i++)
                        for (int j = 0; j < k; j++)
                                for (int r = 0; r < n; r++)
                                {
                                        result[i][j] += a[i][r] * b[r][j];
                                }
                return result;
        }
}

vector<vector<double>> sort_column_wise(std::vector<std::vector<double>> v, std::vector<double> sin_theta, int k)
{
        vector<vector<double>> res(v.size(), vector<double>());
        vector<double> copy_s = sin_theta;
        vector<int> index;
        int in = 0;
        while (in < k)
        {
                int max_i = 0;
                double max_val = 0;
                for (int i = 0; i < copy_s.size(); i++)
                {
                        if (copy_s[i] > max_val)
                        {
                                max_val = copy_s[i];
                                max_i = i;
                        }
                }

                index.push_back(max_i);
                copy_s[max_i] = 0;
                in++;
        }
        for (int i = 0; i < index.size(); i++)
        {
                for (int j = 0; j < v.size(); j++)
                {
                        res[j].push_back(v[j][index[i]]);
                }
        }
        return res;
}

void calculate_time_spend(double time_1, timeval start, timeval end)
{
        time_1 = (end.tv_sec - start.tv_sec) * 1000.0;
        time_1 = time_1 + (end.tv_usec - start.tv_usec) / 1000.0;
        cout << "Total time spent: " << time_1 << " milli seconds." << endl
             << endl;
}
void filling_matrix_from_stream(Vf&A , string filename , int rows , int cols)
{
        ifstream filestream(filename);
        for (  int i = 0 ; i  < rows ; i++)
                for (int j = 0 ; j < cols ; j++)
                filestream>>A[i][j];

}
int main(int argc, char *argv[])
{

        int rows, columns;
        rows = atoi(argv[1]);
        columns = atoi(argv[2]);
        string s ;
        num = int(atoi(argv[3]));
        s = atoi(argv[4]);
        int percentage_to_compress ;
        percentage_to_compress = int(atoi(argv[5]));
        int vectors_to_take = rows  - int(percentage_to_compress*rows/100);

        double time_1;
        double time_2;
        timeval start;
        timeval end;
        timeval end2;

        checking_input_parameters(argc, rows, columns);
        vector<vector<double>> A(columns, vector<double>(columns, 0));
        vector<vector<double>> U(columns, vector<double>(columns, 0));
        vector<vector<double>> V(columns, vector<double>(columns, 0));
        vector<double> S(columns, 0);
        
        filling_matrix_from_stream(A,s,rows,columns);

        gettimeofday(&start, NULL);
        parallel_svd(rows, columns, A, U, S, V);
        gettimeofday(&end, NULL);

        Vf projection_matrix = sort_column_wise(V, S, vectors_to_take);
        cout << "SHape of projection matrix" << endl;
        cout << projection_matrix.size() << endl;
        cout << projection_matrix[0].size() << endl;
        cout << endl;

        Vf compressed_matrix = dense_matrix_mul(A, projection_matrix);
        cout << "SHape of compressed matrix" << endl;
        cout << compressed_matrix.size() << endl;
        cout << compressed_matrix[0].size() << endl;
        cout << endl;

        cout << "Recoonstructing the iamge " << endl;

        Vf projection_matrix_tranpose = transpose_matrix(projection_matrix);
        cout << projection_matrix_tranpose.size() << endl;
        cout << projection_matrix_tranpose[0].size() << endl;
        Vf reconstructed_matrix = dense_matrix_mul(compressed_matrix, projection_matrix_tranpose);
        cout << "SHape of reconstructed matrix" << endl;
        cout << reconstructed_matrix.size() << endl;
        cout << reconstructed_matrix[0].size() << endl;
        cout << endl;

        ofstream matrixfile_output("reconstructed_matrix.txt");

        for (int i = 0; i < reconstructed_matrix.size(); i++)
        {
                for (int j = 0; j < reconstructed_matrix[i].size(); j++)
                {
                        matrixfile_output << reconstructed_matrix[i][j];
                        matrixfile_output << " ";
                }
                matrixfile_output << endl;
        }

        calculate_time_spend(time_1, start, end);

        return 0;
}
