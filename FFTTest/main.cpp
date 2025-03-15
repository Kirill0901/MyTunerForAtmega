#include <complex>
#include <vector>
#include <iostream>
#include <fstream>
#define MAXN 256
#define LOGN 8

using namespace std;

typedef complex<double> base;

const double PI = acos(-1);
const base I(0.0, 1.0);

// X - входные значения
void fft(base X[], int n, int chet){
    if (n == 1){
        return;
    }
    base X0[n / 2];
    base X1[n / 2];
    for (int i = 0, p = 0; i < n; i += 2, ++p){
        X0[p] = X[i];
        X1[p] = X[i + 1];
    }
    fft(X0, n / 2, chet + chet);
    fft(X1, n / 2, chet);

    for (int i = 0; i < n / 2; ++i){
        base w = exp(-2 * M_PI * I * double(i)/ double(n));
        X[i] = X0[i] + X1[i] * w;
        X[i + n / 2] = X0[i] - X1[i] * w;
    }
}

// function to reverse bits of a number
unsigned int reverseBits(unsigned int n, unsigned int bits)
{
    unsigned int rev = 0;

    // traversing bits of 'n' from the right
    while (bits > 0) {
        // bitwise left shift
        // 'rev' by 1
        rev <<= 1;

        // if current bit is '1'
        if ((n & 1) == 1)
            rev ^= 1;

        // bitwise right shift
        // 'n' by 1
        n >>= 1;
        bits--;
    }

    // required number
    return rev;
}

void fft_less_memory(base X[], base C[], int n, int log2n){
    // Устанавливаем нужный порядок для нерекурсивного FFT
    for (int i = 0; i < n; ++i){
        unsigned int ind = reverseBits(i, log2n);
        X[i] = {X[i].real(), X[ind].real()};
    }
    for (int i = 0; i < n; ++i){
        X[i] = {X[i].imag(), 0};
    }

    for (int k = 2; k <= n; k <<= 1){
        int k_half = k >> 1;
        int ndk = n / k;
        for (int i = 0; i < n; i += k){
            for (int j = i, t = 0; j < i + k_half; ++j, t += ndk){
                base w = C[t];
                base first = X[j];
                X[j] = X[j] + X[j + k_half] * w;
                X[j + k_half] = first - X[j + k_half] * w;
            }
        }
    }
}

int main(){

    int n = MAXN, log2n = LOGN;

    base X[n], X1[n];
    for (int i = 0; i < n; ++i){
        double x = 255 / 4;
        cin >> x;
        x *= 0.53836 - 0.46164 * cos(2 * PI * i / (MAXN - 1));
        X[i] = x;
        X1[i] = x;
    }

    base C[n / 2];
    for (int i = 0; i < n / 2; ++i){
        C[i] = exp(-2 * M_PI * I * double(i)/ double(n));
    }

    fft(X, n, 1);
    fft_less_memory(X1, C, n, log2n);
    for (int i = 0; i < n; ++i){
        //cout << X1[i] << "\n";
        cout << X1[i].real() * X1[i].real() + X1[i].imag() * X1[i].imag() << "\n";
    }
}
