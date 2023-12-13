#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "iostream"
#include <ctime>
#include <fstream>
#include <string>
using namespace std;

const double PI = 3.141592653589793238463;
const double  k = 0.59;
const double  c = 1.65;
const double  R = 7000;
const double  T = 200;
const double const S = 0.04;
const double  L = 2 * PI * R;
//const int I_array = 16000;
//const int K_array = 250;
//const int I_array[] = { 8000, 16000, 50000 };
//const int K_array[] = { 10 ,  50,   190 };
const int I_array[] = { 50000 };
const int K_array[] = { 190 };
const double a = sqrt(k / c);
const string extension = ".txt";

double phi(double x) {
    if (-PI * R / 2 <= x && x <= PI * R / 2) {
        return 1.0;
    }
    return 0.0;

}


void writeArray(string filename, double* my_array, int size) {
    ofstream outputFile(filename);
    if (outputFile.is_open()) {
        for (int i = 0; i < size; i++) {
            outputFile << my_array[i] << " ";
        }
        outputFile.close();
        std::cout << "Массив успешно записан в файл " << filename << std::endl;
    }
    else {
        std::cerr << "Не удалось открыть файл для записи." << std::endl;
    }
}


void print_Array(string name_array, double* my_array, int size) {
    cout << name_array << "---->>> " << endl;
    for (int i = 0; i < size; i++) {
        cout << my_array[i] << " ";
    }
    cout << endl;
}



string areArraysEqual(double** u_1, double** u_2, int K, int I) {
    for (int i = 0; i <= K; i++) {
        for (int j = 0; j <= I; j++) {
            if (u_1[i][j] != u_2[i][j]) return "Не равны!!!";
        };
    }
    return "Равны";
}
int main(int argc, char** argv) {
    int length = sizeof(I_array) / sizeof(I_array[0]);
    for (int count = 0; count < length; count++) {
        int I = I_array[count] - 1;
        int K = K_array[count] - 1;

        double step_x = (2 * PI * R) / (I_array[count] - 1);
        double step_t = T / (K_array[count] - 1);
        double gamma = (k * step_t) / (c * pow(step_x, 2));

        MPI_Init(&argc, &argv);






        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);





        cout << I << " " << K << endl;

        int size_for_proc = (I + 1) / size;
        int start_in_proc = rank * size_for_proc;
        int end_in_proc = start_in_proc + size_for_proc - 1;

        //cout << "Proc " << rank << " start = " << start_in_proc << " end = " << end_in_proc << " size = " << size_for_proc << endl;


        int rank_per_gran;
        int rank_per_period;
        int tag = (rank % 2 == 0) ? 1 : 0;;
        if (rank == 0) {
            rank_per_gran = size_for_proc - 1;
            rank_per_period = 1;
        }
        else {
            rank_per_gran = 0;
            rank_per_period = size_for_proc - 2;
        }
        cout << "Proc " << rank << " start = " << start_in_proc << " end = " << end_in_proc << " size = " << size_for_proc << " rank_per_gran = " << rank_per_gran << " rank_per_period = " << rank_per_period << endl;

        double* x_array_loc = new double[size_for_proc];
        for (int i = 0; i < size_for_proc; i++) {
            x_array_loc[i] = -PI * R + (i + start_in_proc) * step_x;
        }

        double** u_loc = new double* [2];
        for (int j = 0; j < 2; j++) {
            u_loc[j] = new double[size_for_proc];
            for (int i = 0; i < size_for_proc; i++) {
                if (j == 0) {
                    u_loc[j][i] = phi(x_array_loc[i]);
                }
                else {
                    u_loc[j][i] = 0.0;
                }

            }
        }


        
        int nubmer_new = 1;
        int nubmer_old = 0;
        int temp = 0;
        double u_gran = 0.0;
        double u_period = 0.0;
        clock_t start_time_parallel = clock();
        for (int k = 0; k < K; k++) {
            
            /*nubmer_new = (k % 2 == 0) ? 1 : 0;
            nubmer_old = (k % 2 != 0) ? 1 : 0;*/

            MPI_Send(&u_loc[nubmer_old][rank_per_gran], 1, MPI_DOUBLE, tag, 0, MPI_COMM_WORLD);
            MPI_Recv(&u_gran, 1, MPI_DOUBLE, tag, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Send(&u_loc[nubmer_old][rank_per_period], 1, MPI_DOUBLE, tag, 0, MPI_COMM_WORLD);
            MPI_Recv(&u_period, 1, MPI_DOUBLE, tag, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


            if (rank == 0) {
                u_loc[nubmer_new][size_for_proc - 1] = u_loc[nubmer_old][size_for_proc - 1] + gamma * (u_gran - 2 * u_loc[nubmer_old][size_for_proc - 1] + u_loc[nubmer_old][size_for_proc - 1 - 1]);
                u_loc[nubmer_new][0] = (u_loc[nubmer_old][1] + u_period) / 2;
            }
            else {
                u_loc[nubmer_new][0] = u_loc[nubmer_old][0] + gamma * (u_loc[nubmer_old][1] - 2 * u_loc[nubmer_old][0] + u_gran);
                u_loc[nubmer_new][size_for_proc - 1] = (u_period + u_loc[nubmer_old][size_for_proc - 2]) / 2;
            }
            for (int i = 1; i < size_for_proc - 1; i++) {
                u_loc[nubmer_new][i] = u_loc[nubmer_old][i] + gamma * (u_loc[nubmer_old][i + 1] - 2 * u_loc[nubmer_old][i] + u_loc[nubmer_old][i - 1]);
            }
            temp = nubmer_new;
            nubmer_new = nubmer_old;
            nubmer_old = temp;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) {
            clock_t end_time_parallel = clock();
            double duration_parallel = static_cast<double>(end_time_parallel - start_time_parallel) / CLOCKS_PER_SEC;
            cout << " Time - >" << duration_parallel << endl;
        }
        string name_file_x_array = "E:\\study\\Мага\\Параллельные\\arrays\\x_array_MPI_"+to_string(rank) + extension;
        string name_file_u_loc = "E:\\study\\Мага\\Параллельные\\arrays\\u_t_array_MPI_"+to_string(rank) + extension;



        writeArray(name_file_x_array, x_array_loc, size_for_proc);
        writeArray(name_file_u_loc, u_loc[nubmer_new], size_for_proc);
        delete[] x_array_loc;
        for (int j = 0; j < 2; j++) {
            delete[] u_loc[j];
        }
        delete[] u_loc;


        MPI_Finalize();
    }

    return 0;
}
