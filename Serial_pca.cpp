#include <bits/stdc++.h>
#include <sys/time.h>
#define epsilon 1e-08
#define Vf vector< vector < double > > 
using namespace std;




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


template <typename T>
double sgn(T val)
{
        return (val > T(0)) - (val < T(0));
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



pair<double, double> schur(vector<vector<double>>& A, int p, int q, double &offA)
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
        offA = max(offA,abs(c)/sqrt(app * aqq)) ;

        
    }

    else
    {
        cs = make_pair(1.0,0.0) ;
        offA = max(offA,1.0 / sqrt(app * aqq)) ;
    }

    return cs ;
}


void serial_svd(int n, vector<vector<double>>& A,vector<vector<double>>& U,vector<double>& S,vector<vector<double>>& V)
{
    double offA = 0.0 ;
    fill_ones_in_diagonal(V);
    

   
    while(offA > epsilon)
    {
      

        offA = 0.0 ;
        for(int i = 0 ; i < n ; i++)
        {
            for(int j = i+1 ; j < n ; j++)
            {
                pair<double, double> cs = schur(A,i,j,offA) ;

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
        string s ;

        s = argv[3];
        int percentage_to_compress ;
        percentage_to_compress = int(atoi(argv[4]));
        int vectors_to_take = M  - int(percentage_to_compress*M/100);
        double elapsedTime, elapsedTime2;
        timeval start, end, end2;
        checking_input_parameters(argc, M, N);

        vector<vector<double>> A(N, vector<double>(N, 0));
        vector<vector<double>> U(N, vector<double>(N, 0));
        vector<vector<double>> V(N, vector<double>(N, 0));
        vector<double> S(N, 0);


        filling_matrix_from_stream(A,s,M,N);
        // Calculate A.TA

        // Vf AT = transpose_matrix(A) ;
        // Vf ATA = dense_matrix_mul(AT,A) ;


 
 		gettimeofday(&start, NULL);
        for(int i = 0 ; i < 8 ; i++)
            serial_svd(N,A,U,S,V);
        gettimeofday(&end, NULL);
      
      

        
        Vf projection_matrix = sort_column_wise(V,S,vectors_to_take);
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
