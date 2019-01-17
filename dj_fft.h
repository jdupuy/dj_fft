/* dj_fft.h - public domain FFT library
by Jonathan Dupuy

   INTERFACING

   define DJ_ASSERT(x) to avoid using assert.h.

   QUICK NOTES

*/

#ifndef DJ_INCLUDE_FFT_H
#define DJ_INCLUDE_FFT_H

#include <complex> // std::complex
#include <vector>  // std::vector

namespace dj {
namespace fft {
    // FFT argument: std::vector<std::complex>
    template <typename T> using arg = std::vector<std::complex<T>>;

    // FFT direction specifier
    enum class e_dir {DIR_FWD = +1, DIR_BWD = -1};

    // FFT routines
    template <typename T> arg<T> eval_1d(const arg<T> &xi, const e_dir &dir);
    template <typename T> arg<T> eval_2d(const arg<T> &xi, const e_dir &dir);
    template <typename T> arg<T> eval_3d(const arg<T> &xi, const e_dir &dir);
} // namespace fft
} // namespace dj

//
//
//// end header file ///////////////////////////////////////////////////////////
#endif // DJ_INCLUDE_FFT_H

#include <cmath>
#include <cstdint>
#include <cstring> // std::memcpy

#ifndef DJ_ASSERT
#   include <cassert>
#   define DJ_ASSERT(x) assert(x)
#endif

namespace dj {
namespace fft {

constexpr auto Pi = 3.141592653589793238462643383279502884;

/*
 * Returns offset to most significant bit
 * NOTE: only works for positive power of 2s
 * examples:
 * 1b      -> 0d
 * 100b    -> 2d
 * 100000b -> 5d
 */
int findMSB(int x)
{
    DJ_ASSERT(x > 0 && "invalid input");
    int p = 0;

    while (x > 1) {
        x>>= 1;
        ++p;
    }

    return p;
}


/*
 *  Bit reverse an integer given a word of nb bits
 *  NOTE: Only works for 32-bit words max
 *  examples:
 *  10b      -> 01b
 *  101b     -> 101b
 *  1011b    -> 1101b
 *  0111001b -> 1001110b
 */
int bitr(uint32_t x, int nb)
{
    DJ_ASSERT(nb > 0 && 32 > nb && "invalid bit count");
    x = ( x               << 16) | ( x               >> 16);
    x = ((x & 0x00FF00FF) <<  8) | ((x & 0xFF00FF00) >>  8);
    x = ((x & 0x0F0F0F0F) <<  4) | ((x & 0xF0F0F0F0) >>  4);
    x = ((x & 0x33333333) <<  2) | ((x & 0xCCCCCCCC) >>  2);
    x = ((x & 0x55555555) <<  1) | ((x & 0xAAAAAAAA) >>  1);

    return ((x >> (32 - nb)) & (0xFFFFFFFF >> (32 - nb)));
}


/*
 * Computes a Fourier transform, i.e.,
 * xo[k] = 1/sqrt(N) sum(j=0 -> N-1) xi[j] exp(i 2pi j k / N)
 * with O(N log N) complexity using the butterfly technique
 *
 * NOTE: Only works for arrays whose size is a power-of-two
 */
template <typename T> arg<T> eval_1d(const arg<T> &xi, const e_dir &dir)
{
    DJ_ASSERT((xi.size() & (xi.size() - 1)) == 0 && "invalid input size");
    int cnt = (int)xi.size();
    int msb = findMSB(cnt);
    T nrm = T(1) / std::sqrt(cnt);
    arg<T> xo(cnt);

    // pre-process the input data
    for (int j = 0; j < cnt; ++j)
        xo[j] = nrm * xi[bitr(j, msb)];

    // fft passes
    for (int i = 0; i < msb; ++i) {
        int bm = 1 << i; // butterfly mask
        int bw = 2 << i; // butterfly width
        T ang = T(dir) * Pi / T(bm); // precomputation

        // fft butterflies
        for (int j = 0; j < (cnt/2); ++j) {
            int i1 = ((j >> i) << (i + 1)) + j % bm; // left wing
            int i2 = i1 ^ bm;                        // right wing
            T a1 = ang * (i1 ^ bw); // left wing rotation
            T a2 = ang * (i2 ^ bw); // right wing rotation
            std::complex<T> tmp = xo[i1];

            xo[i1]+= std::polar(T(1), a1) * xo[i2];
            xo[i2] = tmp + std::polar(T(1), a2) * xo[i2];
        }
    }

    return xo;
}

/*
 * Computes a 2D Fourier transform
 * with O(N^2 log N) complexity using the butterfly technique
 *
 * NOTE: the input must be a square matrix whose size is a power-of-two
 */
template <typename T> arg<T> eval_2d(const arg<T> &xi, const e_dir &dir)
{
    DJ_ASSERT((xi.size() & (xi.size() - 1)) == 0 && "invalid input size");
    int cnt2 = (int)xi.size();   // NxN
    int msb = findMSB(cnt2) / 2; // lg2(N) = lg2(sqrt(NxN))
    int cnt = 1 << msb;          // N = 2^lg2(N)
    T nrm = T(1) / T(cnt);
    arg<T> xo(cnt2);

    // pre-process the input data
    for (int j2 = 0; j2 < cnt; ++j2)
    for (int j1 = 0; j1 < cnt; ++j1) {
        int k2 = bitr(j2, msb);
        int k1 = bitr(j1, msb);

        xo[j1 + cnt * j2] = nrm * xi[k1 + cnt * k2];
    }

    // fft passes
    for (int i = 0; i < msb; ++i) {
        int bm = 1 << i; // butterfly mask
        int bw = 2 << i; // butterfly width
        float ang = T(dir) * Pi / T(bm); // precomputation

        // fft butterflies
        for (int j2 = 0; j2 < (cnt/2); ++j2)
        for (int j1 = 0; j1 < (cnt/2); ++j1) {
            int i11 = ((j1 >> i) << (i + 1)) + j1 % bm; // xmin wing
            int i21 = ((j2 >> i) << (i + 1)) + j2 % bm; // ymin wing
            int i12 = i11 ^ bm;                         // xmax wing
            int i22 = i21 ^ bm;                         // ymax wing
            T a11 = ang * (i11 ^ bw); // upper left rotation
            T a12 = ang * (i12 ^ bw); // upper right rotation
            T a21 = ang * (i21 ^ bw); // lower left rotation
            T a22 = ang * (i22 ^ bw); // lower right rotation
            int k11 = i11 + i21 * cnt; // array offset
            int k12 = i12 + i21 * cnt; // array offset
            int k21 = i11 + i22 * cnt; // array offset
            int k22 = i12 + i22 * cnt; // array offset

            // FFT-X
            std::complex<T> tmp1 = xo[k11];
            std::complex<T> tmp2 = xo[k21];
            xo[k11]+= std::polar(T(1), a11) * xo[k12];
            xo[k12] = tmp1 + std::polar(T(1), a12) * xo[k12];
            xo[k21]+= std::polar(T(1), a11) * xo[k22];
            xo[k22] = tmp2 + std::polar(T(1), a12) * xo[k22];

            // FFT-Y
            std::complex<T> tmp3 = xo[k11];
            std::complex<T> tmp4 = xo[k12];
            xo[k11]+= std::polar(T(1), a21) * xo[k21];
            xo[k21] = tmp3 + std::polar(T(1), a22) * xo[k21];
            xo[k12]+= std::polar(T(1), a21) * xo[k22];
            xo[k22] = tmp4 + std::polar(T(1), a22) * xo[k22];
        }
    }

    return xo;
}



} // namespace fft
} // namespace dj
