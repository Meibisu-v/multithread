#include <iostream>
#include <thread>
#include <vector>
#include <algorithm>
#include <fstream>
#include <time.h>
#define MATRIX_SIZE 1500
#define MAX_THREAD_NUM 1500
#define BLOCK_SIZE 64
#define EPS 0.00001

void create_matrix(double **matrix, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        matrix[i] = new double[n];
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            matrix[i][j] = rand()%3;
        }
    }

}

void print_matrix(double **matrix, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            std::cout << matrix[i][j] << ' ';
        }
        std::cout << '\n';
    }
}

void delete_matrix(double **matrix, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        delete[] matrix[i];
    }
    delete[]matrix;
}

void simple_mult(double **first, double **second, double **result, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            for (size_t k = 0; k < n; ++k) {
                result[i][j] += first[i][k] * second[k][j];
            }
        }
    }
}

void simple_mult_with_treads(double **first, double **second, double **result, size_t n, size_t thread_num){
    std::vector<std::thread> threads;
    threads.reserve(thread_num);
    auto row_mult = [](double **first, double **second, double **result, size_t row, size_t n, size_t row_count){
        for (size_t shift = 0; shift < row_count && row + shift < n; ++shift) {
            for (size_t i = 0; i < n; ++i) {
                double temp = 0;
                for (size_t j = 0; j < n; ++j) {
                    temp += first[row + shift][j] * second[j][i];
                }
                result[row + shift][i] = temp;
            }
        }
    };
    for (size_t i = 0; i < n; i += n / thread_num) {
        threads.push_back(std::thread(row_mult, first, second, result, i , n, n / thread_num));
    }
    //count rest
    if (n % thread_num != 0)
        threads.emplace_back(std::thread(row_mult, first, second, result, n - n % thread_num, n, n % thread_num));
    for (auto&& thread : threads) {
        thread.join();
    }
}

void block_mult(double **first, double **second, double **result, size_t size, size_t block_size) {
    for (size_t i = 0; i < size; i += block_size) {
        for (size_t j = 0; j < size; j += block_size) {
            for (size_t k = 0; k < size; ++k) {
                for (size_t i_block = i; i_block < std::min(i + block_size, size); ++i_block) {
                    double temp = 0;
                    for (size_t j_block = j; j_block < std::min(j + block_size, size); ++j_block){
                        temp += first[k][j_block] * second[j_block][i_block];
                    }
                    result[k][i_block] += temp;
                }
            }
        }
    }
}

void block_mult_with_threads(double **first, double **second, double **result, size_t size, size_t block_size) {
    std::vector<std::thread> threads;
    auto block_mult = [](double **first, double **second, double **result, size_t i , size_t j, size_t size, size_t block_size) {
        for (size_t k = 0; k < size; ++k) {
            for (size_t i_block = i; i_block < std::min(i + block_size, size); ++i_block) {
                double temp = 0;
                for (size_t j_block = j; j_block < std::min(j + block_size, size); ++j_block){
                    temp += first[k][j_block] * second[j_block][i_block];
                }
                result[k][i_block] += temp;
            }
        }
    };
    for (size_t i = 0; i < size; i += block_size) {
        for (size_t j = 0; j < size; j += block_size) {
            threads.emplace_back(std::thread(block_mult, first, second, result, i, j, size, block_size));
        }
    }    
    // std::cout << block_size << ' ' << threads.size() << '\n';
    for (auto&& thread : threads) {
        thread.join();
    }
}

void reset_matrix(double **matrix, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
        matrix[i][j] = 0;
        }
    }
}

bool is_equal(double **a, double ** b, size_t n){
    bool is_equal = true;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (abs(a[i][j] - b[i][j]) > EPS) {
                is_equal = false;
                return is_equal;
            }
        }
    }
    return is_equal;
}

void test_mult(double **first, double **second, double **result, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            for (size_t k = 0; k < n; ++k) {
                result[i][j] += first[i][k] * second[k][j];
            }
        }
    }
}

int main(int argc, char *argv[]) {
    size_t n = MATRIX_SIZE;
    size_t block_size = BLOCK_SIZE;
    double **first_matrix = new double*[n];
    double **sec_matrix = new double*[n];
    double **result_matrix = new double*[n];
    // double **test_result = new double*[n];
    create_matrix(first_matrix, n);
    create_matrix(sec_matrix, n);
    create_matrix(result_matrix, n);
    // create_matrix(test_result, n);

    reset_matrix(result_matrix, n);
    // reset_matrix(test_result, n);
    std::vector<size_t> size = {8, 20, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1150, 1200}; 
    std::ofstream out;

    out.open ("simple_mult_one_thread.txt");    
    for (size_t n : size) {
        clock_t tStart = clock();
        simple_mult(first_matrix, sec_matrix, result_matrix, n);
        out << n << ' ' << (double)(clock() - tStart)/CLOCKS_PER_SEC << std::endl;
        reset_matrix(result_matrix, n); 
    }
    out.close();

    out.open("block_mult_one_thread.txt");
    for (size_t n : size) {
        clock_t tStart = clock();
        block_mult(first_matrix, sec_matrix, result_matrix, n, block_size);
        out << n << ' ' << (double)(clock() - tStart)/CLOCKS_PER_SEC << std::endl;
        reset_matrix(result_matrix, n);
    }
    out.close();

    out.open("simple_mult_multi_thread.txt");
    for (size_t thread_num = 1; thread_num < MATRIX_SIZE; thread_num += 100) {
        clock_t tStart = clock();
        simple_mult_with_treads(first_matrix, sec_matrix, result_matrix, MATRIX_SIZE, thread_num);
        out << thread_num << ' ' << (double)(clock() - tStart)/CLOCKS_PER_SEC << std::endl;
        reset_matrix(result_matrix, MATRIX_SIZE);
    }
    out.close();

    out.open("block_mult_multi_thread.txt");
    for(size_t thread_num = 1; thread_num * thread_num <= MAX_THREAD_NUM; ++thread_num) {
        clock_t tStart = clock();
        block_mult_with_threads(first_matrix, sec_matrix, result_matrix, MATRIX_SIZE, MATRIX_SIZE / thread_num);
        out << thread_num * thread_num << ' ' << (double)(clock() - tStart)/CLOCKS_PER_SEC << std::endl;
        reset_matrix(result_matrix, MATRIX_SIZE);
    }
    out.close();
    // delete_matrix(test_result, MATRIX_SIZE);
    delete_matrix(first_matrix, MATRIX_SIZE);
    delete_matrix(sec_matrix, MATRIX_SIZE);
    delete_matrix(result_matrix, MATRIX_SIZE);
}