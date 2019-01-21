# dj_fft: Header-only FFT library

[![Build Status](https://travis-ci.org/jdupuy/dj_fft.svg?branch=master)](https://travis-ci.org/jdupuy/dj_fft)
[![Build status](https://ci.appveyor.com/api/projects/status/nwcgmc1l74h8sudk?svg=true)](https://ci.appveyor.com/project/jdupuy/dj-fft)

### Details
This repository provides a header-only library to compute fourier transforms in 1D, 2D, and 3D. Its goal is to provide a fast and easy-to-use fast fourier transform algorithm. 

### Cloning

Clone the repository and all its submodules using the following command:
```sh
git clone --recursive git@github.com:jdupuy/dj_fft.git
```

If you accidentally omitted the `--recursive` flag when cloning the repository you can retrieve the submodules like so:
```sh
git submodule update --init --recursive
```

### Usage
The 1D, 2D, and 3D FFT routines return an `std::vector<std::complex<T>>`, given another `std::vector<std::complex<T>>` as input, which specifies the data that must be transformed, as well as an enum class `dj::fft::e_dir`, which specifies in which direction the FFT must be computed (specify `dj::fft::e_dir::DIR_FWD` for the forward direction and `dj::fft::e_dir::DIR_BWD` for the backward direction).

Note that the input vector is expected to be of size `N` for 1D FFT, `NxN` for a 2D FFT, and `NxNxN` for a 3D FFT, where `N` *must* be a power of two. Note that the 2D and 3D vectors are expected to be arranged in a flat row-major fashion, i.e., the 2D and 3D elements `(i, j)` and `(i, j, k)` are respectively located at index `i + N * j` and `i + N * (j + N * k)` in memory.

Below is a C++ pseudocode for computing a 2D FFT in forward direction:
```c++
#define DJ_FFT_IMPLEMENTATION // define this in exactly *one* .cpp file
#include "dj_fft.h"

some_function()
{
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
}
```

To see examples that compile, see the examples/ directory. 

### GPU Acceleration
Additionally, the library provides GPU accelerated 1D, 2D, and 3D FFTs for `std::vector<std::complex<float>>` inputs. GPU acceleration is especially relevant for large 2D and 3D datasets. For instance:
- for an input of size 4096x4096, a regular 2D FFT completes in roughly 18 seconds on an intel i7-8086k, and 0.9 seconds on an NVidia RTX 2080
- for an input of size 512x512x512, a regular 3D FFT completes in roughly 131 seconds on an intel i7-8086k, and 8.2 seconds on an NVidia RTX 2080

GPU acceleration recquires the user to create an OpenGL4.5 context. This can be achieved with a cross-platform windowing library such as GLFW (https://www.glfw.org/), and an OpenGL function loader such as glad (https://glad.dav1d.de/). To enable the GPU accelerated routines, the user must define the GPU acceleration flag by defining a ``#define DJ_FFT_ENABLE_GPU`` before including the ``dj_fft.h`` file. Below is a C++ pseudocode for computing a 1D FFT in backward direction on the GPU:

```c++
#define DJ_FFT_ENABLE_GPU // enables GPU acceleration
#define DJ_FFT_IMPLEMENTATION // define this in exactly *one* .cpp file
#include "dj_fft.h"

some_function()
{
  int N = size_of_your_input; // input size
  auto myData = std::vector<std::complex<float>>(N); // input data

  // prepare data
  for (int i = 0; i < N; ++i) {
    myData[i] = some_float_value; // set element (i)
  }

  // compute backward 1D FFT
  auto fftData = dj::fft::eval_1d_gpu(myData, dj::fft::e_dir::FFT_BWD);

  // print the data
  for (int i = 0; i < N; ++i) {
    printf("{%f, %f}\n", fftData[i].real(), fftData[i].imag());
  }
}
```
Note that the return values of a GPU FFT may differ slightly from that of a regular FFT, due to the way floating point arithmetic is implemented.

For a complete example that compiles, see the examples/ directory. 

### License
The code from this repository is released in public domain.
