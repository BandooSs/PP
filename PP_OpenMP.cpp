#include <iostream>
#include <omp.h>
#include <cmath>
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
//const int I_array[] = { 100,1000,2000,4000,8000,10000 };
//const int K_array[] = { 10,200,800,2200,4400,7000 };

//const int I_array[] = { 100 ,1000,2000,4000,8000,16000 };
//const int K_array[] = { 10,25,100,400,1000,2000 };
// 
//const int I_array[] = { 100 ,1000, 2000, 4000,8000, 16000, 20000  };
//const int K_array[] = { 5,  10,   20,    30,  60 ,  150,   250   }; 

const int I_array[] = { 8000, 16000, 50000 };
const int K_array[] = { 10 ,  50,   190 };
// 
// 
// 
//const int I_array[] = { 6  };
//const int K_array[] = { 6 };
const double a = sqrt(k / c);
const string extension = ".txt";


double phi(double x) {
    if (-PI * R / 2 < x && x < PI * R / 2) {
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

int main()
{
    setlocale(LC_ALL, "Russian");
    int num_threads = 2;  // Замените на желаемое количество потоков
    int length = sizeof(I_array) / sizeof(I_array[0]);
    //omp_set_num_threads(num_threads);
    cout << "Всего потоков - > " << num_threads << endl;
    double* massive_time = new double[length];
    for (int step = 0; step < length; step++) {
        massive_time[step] = 0.0;
    }
    for (int step = 0; step < 12; step++) {
        
        for (int count = 0; count < length; count++) {

            double step_x = (2 * PI * R) / (I_array[count] - 1);
            double step_t = T / (K_array[count] - 1);
            double* x_array = new double[I_array[count]];
            double* t_array = new double[K_array[count]];
            
            for (int i = 0; i < I_array[count]; i++) {
                x_array[i] = -PI * R + i * step_x;
            }

            for (int i = 0; i < K_array[count]; i++) {
                t_array[i] = 0 + i * step_t;
            }

            int I = I_array[count] - 1;
            int K = K_array[count] - 1;
            cout << I_array[count] << " " << K_array[count] << endl;;
            
            omp_set_num_threads(num_threads);
            
            /*for (int k = 0; k < K; k++) {
                #pragma omp parallel for shared(u_new,u_old)
                for (int i = 1; i < I; i++) {
                    u_new[i] = u_old[i] + gamma * (u_old[i + 1] - 2 * u_old[i] + u_old[i - 1]);
                }
                u_new[I] = (u_old[1] + u_old[I - 1]) / 2;
                u_new[0] = u_new[I];
                #pragma omp parallel for shared(u_new,u_old)
                for (int i = 0; i <= I; i++) {
                    u_old[i] = u_new[i];
                }

            }*/
            double** u = new double* [2];
            for (int j = 0; j < 2; j++) {
                u[j] = new double[I + 1];
            }
            for (int j = 0; j < 2; j++) {
                for (int i = 0; i <= I; i++) {
                    if (j == 0) {
                        u[j][i] = phi(x_array[i]);
                    }
                    else {
                        u[j][i] = 0.0;
                    }

                }
            }
            double gamma = (k * step_t) / (c * pow(step_x, 2));
            int nubmer_new = 0;
            int nubmer_old = 0;

            clock_t start_time_parallel = clock();
            for (int k = 0; k < K; k++) {
                nubmer_new = (k % 2 == 0) ? 1 : 0;
                nubmer_old = (k % 2 != 0) ? 1 : 0;
                //cout << " u_new-> " << nubmer_new << " u_old-> " << nubmer_old << endl;
                #pragma omp parallel for shared(u)
                for (int i = 1; i < I; i++) {
                    u[nubmer_new][i] = u[nubmer_old][i] + gamma * (u[nubmer_old][i + 1] - 2 * u[nubmer_old][i] + u[nubmer_old][i - 1]);
                }
                #pragma omp sections
                {
                    #pragma omp section
                    u[nubmer_new][I] = (u[nubmer_old][1] + u[nubmer_old][I - 1]) / 2;
                    #pragma omp section
                    u[nubmer_new][0] = u[nubmer_new][I];
                }

            }
            clock_t end_time_parallel = clock();
            
            
            // Рассчитываем разницу между начальным и конечным временем
            double duration_parallel = static_cast<double>(end_time_parallel - start_time_parallel) / CLOCKS_PER_SEC;

            cout << " Time - >" << duration_parallel << endl;
            massive_time[count] += duration_parallel;

            /*string name_file_x_array = "E:\\study\\Мага\\Параллельные\\arrays\\x_array_OpenMP_" + to_string(count) + extension;
            string name_file_u_t_array = "E:\\study\\Мага\\Параллельные\\arrays\\u_t_array_OpenMP_" + to_string(count) + extension;
            string name_file_t_array = "E:\\study\\Мага\\Параллельные\\arrays\\t_array_OpenMP_" + to_string(count) + extension;

            writeArray(name_file_x_array, x_array, I_array[count]);
            writeArray(name_file_u_t_array, u[nubmer_new], I_array[count]);
            writeArray(name_file_t_array, t_array, K_array[count]);*/

            delete[] x_array;
            delete[] t_array;
            for (int j = 0; j < 2; j++) {
                delete[]u[j];
            }
            delete[]u;
            
        }
        cout << endl << endl;
        
    }
    cout << "Среднее время:\n";
    for (int step = 0; step < length; step++) {
        cout<<I_array[step]<<" "<< K_array[step]<<" "<<massive_time[step]/12<<endl;
    }
    
    return 0;
}



