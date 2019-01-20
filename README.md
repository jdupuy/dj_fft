# dj_fft: Header-only FFT library

[![Build Status](https://travis-ci.org/jdupuy/dj_fft.svg?branch=master)](https://travis-ci.org/jdupuy/dj_fft)
[![Build status](https://ci.appveyor.com/api/projects/status/nwcgmc1l74h8sudk?svg=true)](https://ci.appveyor.com/project/jdupuy/dj-fft)

### Details
This repository provides a header-only library to compute fourier transforms in 1D, 2D, and 3D. Its goal is to provide a fast and easy-to-use fast fourier transform algorithm. 

### Usage
The 1D, 2D, and 3D FFT routines return an `std::vector<std::complex<T>>`, given another `std::vector<std::complex<T>>` as input, which specifies the data that must be transformed, as well as an enum class `dj::fft::e_dir`, which specifies in which direction the FFT must be computed (specify `dj::fft::e_dir::DIR_FWD` for the forward direction and `dj::fft::e_dir::DIR_BWD` for the backward direction).

Note that the input vector is expected to be of size `N` for 1D FFT, `NxN` for a 2D FFT, and `NxNxN` for a 3D FFT, where `N` *must* be a power of two. Note that the 2D and 3D vectors are expected to be arranged in a flat row-major fashion, i.e., the 2D and 3D elements `(i, j)` and `(i, j, k)` are respectively located at index `i + N * j` and `i + N * (j + N * k)` in memory.

Below is a C++ pseudocode for computing a 2D FFT in forward direction:
```c++
int N = size_of_your_input; // input size
auto myData = std::vector<std::complex<T>>(N * N); // input data

// prepare data
for (int j = 0; j < N; ++j) {
  for (int i = 0; i < N; ++i) {
    myData[i + N * j] = some_value; // set element (i, j)
  }
}

// compute forward 2D FFT
auto fftData = dj::fft::eval_2d(myData, dj::fft::e_dir::FFT_FWD);

// print the data
for (int j = 0; j < N; ++j) {
  for (int i = 0; i < N; ++i) {
    printf("{%f, %f} ", fftData[i + N * j].real(), fftData[i + N * j].imag());
  }
  printf("\n");
}
```

To see examples that compile, see the examples/ directory. 

### Future Extensions
Note that at some point I'll also add a GPU port of this code -- stay tuned ;)

### License
The code from this repository is released in public domain.
