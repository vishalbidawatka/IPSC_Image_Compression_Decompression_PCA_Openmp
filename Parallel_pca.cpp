#include <bits/stdc++.h>
#include <sys/time.h>
#include <omp.h>
#define Vf vector< vector < double > > 
#define epsilon 0.5
#define num 4
#define MAX_THREAD 100

using namespace std;

template <typename T>
double sgn(T val)
{
        return (val > T(0)) - (val < T(0));
}

void print_array(int *I1, int *I2, int M, int N)
{
        printf("Array I1");
        cout << endl;
        for (int i = 0; i < M; i++)
        {
                cout << I1[i] << " ";
        }
        cout << endl;
        cout << "Array I2" << endl;
        for (int j = 0; j < N; j++)
        {
                cout << I2[j] << " ";
        }
        cout << endl;
}
void checking_input_parameters(int argc, int M, int N)
{
        if (argc < 2)
        {
                cout << "Invalid input" << endl;
                exit(0);
        }
        if (M != N)
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



void parallel_svd(vector<vector<double>>& A,vector<vector<double>>& U,vector<double>& S,vector<vector<double>>& V)
{       
       
        double c;
        double t;
        double s;
        double converge = 1.0;
        int acum = 0;
        double temp1, temp2;
        int n = A.size() ;
        double elapsedTime, elapsedTime2;
        vector<int> I1( A.size() + 1,0);
        vector<int> I2( A.size() + 1,0);
        timeval start, end, end2;
        vector<double> C(MAX_THREAD);

        fill_ones_in_diagonal(V);

   

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

                            pair<double, double> cs = schur(A,i,j) ;

                            double c = cs.first ;
                            double s = cs.second ;

                            for(int k=0; k<n; k++){
                            double t = A[i][k];
                            A[i][k] = c*t - s*A[j][k];
                            A[j][k] = s*t + c*A[j][k];

                            t = V[i][k];
                            V[i][k] = c*t - s*V[j][k];
                            V[j][k] = s*t + c*V[j][k];

                          }

						#pragma omp parallel for num_threads(num)
                        for (int p = 1; p <= r2; p++)
                        {
                                                        {
                            int k = omp_get_thread_num();
                                int i = I2[p];
                                int j = i + l;

                            pair<double, double> cs = schur(A,i,j) ;

                            double c = cs.first ;
                            double s = cs.second ;

                            for(int k=0; k<n; k++){
                            double t = A[i][k];
                            A[i][k] = c*t - s*A[j][k];
                            A[j][k] = s*t + c*A[j][k];

                            t = V[i][k];
                            V[i][k] = c*t - s*V[j][k];
                            V[j][k] = s*t + c*V[j][k];

                          }
                      }
                        }
                

                for (int k = 0; k < num; k++)
                        converge = max(converge, C[k]);
                }
        }

        for (int i = 0; i < n; i++)
                {

                    double t = 0;
                    for (int j = 0; j < n; j++)
                    {
                            t = t + pow(A[i][j], 2);
                    }
                    t = sqrt(t);

                    for (int j = 0; j < n; j++)
                    {
                        if (i == j)
                            S[i] = t;
                            
                    }
                }


}


}




Vf transpose_matrix(Vf a)
{       
        Vf b(a[0].size()  , vector<double>(a.size(),0));
        for(int i = 0; i < a.size() ; i++)
        {
                for(int j = 0 ; j < a[i].size() ; j++)
                {
                        b[j][i] = a[i][j];
                }
        }
        return b;
}


Vf dense_matrix_mul(Vf &a , Vf &b)
{
    int m = a.size();
    int n = a[0].size();
    int l = b.size();
    int k = b[0].size();

    if( n != l)
    {
        cout<<"Invalid matrix multiplication"<<endl;
        exit(1) ;
    }
    else
    {
        Vf result;
        vector<double> temp(k,0);
        for(int i = 0 ; i < m ; i++)
            result.push_back(temp);
        
        for(int i = 0 ; i < m ; i++)
            for(int j = 0 ; j < k ; j++)
                for(int r = 0 ; r < n ; r++)
                {
                    result[i][j] += a[i][r]*b[r][j] ;
                }
    return result;
    }
    
}


vector<vector<double>> sort_column_wise (std::vector<std::vector<double >> v,std::vector<double> s,int k)
{       
	vector< vector<double> > res(v.size(),vector<double> ());
	vector<double> copy_s = s;
	vector<int>index;
	int in=0;
	while(in<k){
		int max_i=0;
		double max_val=0;
		for(int i=0;i<copy_s.size();i++){
			if(copy_s[i]>max_val){
				max_val=copy_s[i];
				max_i=i;
			}
		}
            
		index.push_back(max_i);
		copy_s[max_i]=0;
		in++;
	}
	for(int i=0;i<index.size();i++){
		for(int j=0;j<v.size();j++){
			res[j].push_back(v[j][index[i]]);
		}
	}
	return res;
}







int main(int argc, char *argv[])
{

        int M, N;

        string T, P, Db;
        M = atoi(argv[1]);
        N = atoi(argv[2]);

        double elapsedTime, elapsedTime2;
        timeval start, end, end2;
        checking_input_parameters(argc, M, N);
       
        vector<vector<double>> A(N, vector<double>(N, 0));
        vector<vector<double>> U(N, vector<double>(N, 0));
        vector<vector<double>> V(N, vector<double>(N, 0));
        vector<double> S(N, 0);


        ifstream matrixfile("0.txt");
        for (int i = 0; i < M; i++)
        {
                for (int j = 0; j < N; j++)
                {

                        matrixfile >> A[i][j];
                }
        } 
 
 		gettimeofday(&start, NULL);
        parallel_svd(A,U,S,V);
        gettimeofday(&end, NULL);
      
      

        
        Vf projection_matrix = sort_column_wise(V,S,400);
        cout<<"SHape of projection matrix"<<endl;
        cout<<projection_matrix.size()<<endl;
        cout<<projection_matrix[0].size()<<endl;
        cout<<endl;

        Vf compressed_matrix = dense_matrix_mul(A,projection_matrix);
        cout<<"SHape of compressed matrix"<<endl;
        cout<<compressed_matrix.size()<<endl;
        cout<<compressed_matrix[0].size()<<endl;
        cout<<endl;


        cout<<"Recoonstructing the iamge "<<endl;

        Vf projection_matrix_tranpose = transpose_matrix(projection_matrix);
        cout<<projection_matrix_tranpose.size()<<endl;
        cout<<projection_matrix_tranpose[0].size()<<endl;
        Vf reconstructed_matrix = dense_matrix_mul(compressed_matrix,projection_matrix_tranpose);
        cout<<"SHape of reconstructed matrix"<<endl;
        cout<<reconstructed_matrix.size()<<endl;
        cout<<reconstructed_matrix[0].size()<<endl;
        cout<<endl;

        ofstream matrixfile_output("reconstructed_matrix.txt");

        for(int i = 0 ; i < reconstructed_matrix.size() ; i++)
        {
                for (int j = 0 ; j < reconstructed_matrix[i].size() ; j++)
                {
                        matrixfile_output << reconstructed_matrix[i][j];
                        matrixfile_output<<" ";
                }
                matrixfile_output<<endl;
        }


    elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0;
    elapsedTime += (end.tv_usec - start.tv_usec) / 1000.0;
    cout<<"Time: "<<elapsedTime<<" ms."<<endl<<endl;



	return 0 ;
}
