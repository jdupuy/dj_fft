
#include "glad/glad.h"
#include "GLFW/glfw3.h"
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <chrono>

#define DJ_FFT_ENABLE_GPU // enable GPU FFT

#define DJ_FFT_IMPLEMENTATION
#include "dj_fft.h"

#define LOG(fmt, ...)  fprintf(stdout, fmt, ##__VA_ARGS__); fflush(stdout);

float rng()
{
    return (double)rand() / (double)RAND_MAX;
}

int main(int argc, char **argv)
{
    // create an OpenGL context
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);

    GLFWwindow* window = glfwCreateWindow(128, 128, "dj_fft", NULL, NULL);
    if (window == NULL) {
        printf("OpenGL window creation failed!\n");
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    // Load OpenGL functions
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        printf("window creation failed!\n");
        return -1;
    }

#if 1 // perform a 1D FFT
    {
        srand(1234);
        int cnt = 8192;
        dj::fft::arg<float> xi;

        for (int i = 0; i < cnt; ++i)
            xi.push_back(std::complex<float>((2.0f * rng() - 1.0f) * 64.f, 0.0f));

        auto s1 = std::chrono::high_resolution_clock::now();
        dj::fft::arg<float> cpu = dj::fft::eval_1d(xi, dj::fft::e_dir::DIR_FWD);
        auto s2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> d1 = s2 - s1;
        printf("=> [%i] CPU: %f s\n", cnt, d1.count()); fflush(stdout);

        auto s3 = std::chrono::high_resolution_clock::now();
        dj::fft::arg<float> gpu = dj::fft::eval_1d_gpu(xi, dj::fft::e_dir::DIR_FWD);
        auto s4 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> d2 = s4 - s3;
        printf("=> [%i] GPU: %f s\n", cnt, d2.count()); fflush(stdout);

#if 0 // display content
    for (int i = 0; i < cpu.size(); ++i) {
        printf("{%f %f} vs {%f %f}\n",
               cpu[i].real(), cpu[i].imag(),
               gpu[i].real(), gpu[i].imag());
    }
#endif
    }
#endif

#if 1 // perform a 2D FFT
    {
        srand(3234);
        int cnt = 4096;
        dj::fft::arg<float> xi;

        for (int i = 0; i < cnt * cnt; ++i)
            xi.push_back(std::complex<float>((2.0f * rng() - 1.0f) * 64.f, 0.0f));

        auto s1 = std::chrono::high_resolution_clock::now();
        dj::fft::arg<float> cpu = dj::fft::eval_2d(xi, dj::fft::e_dir::DIR_FWD);
        auto s2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> d1 = s2 - s1;
        printf("=> [%i^2] CPU: %f s\n", cnt, d1.count()); fflush(stdout);

        auto s3 = std::chrono::high_resolution_clock::now();
        dj::fft::arg<float> gpu = dj::fft::eval_2d_gpu(xi, dj::fft::e_dir::DIR_FWD);
        auto s4 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> d2 = s4 - s3;
        printf("=> [%i^2] GPU: %f s\n", cnt, d2.count()); fflush(stdout);

#if 0 // display content
        for (int j = 0; j < cnt; ++j) {
            for (int i = 0; i < cnt; ++i) {
                printf("{%f %f} vs {%f %f}\n",
                       cpu[i+cnt*j].real(), cpu[i+cnt*j].imag(),
                       gpu[i+cnt*j].real(), gpu[i+cnt*j].imag());
             }
            printf("\n");
        }
#endif
    }
#endif

#if 1 // perform a 3D FFT
    {
        srand(3234);
        int cnt = 512;
        dj::fft::arg<float> xi;

        for (int i = 0; i < (cnt * cnt * cnt); ++i)
            xi.push_back(std::complex<float>((2.0f * rng() - 1.0f) * 64.f, 0.0f));

        auto s1 = std::chrono::high_resolution_clock::now();
        dj::fft::arg<float> cpu = dj::fft::eval_3d(xi, dj::fft::e_dir::DIR_FWD);
        auto s2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> d1 = s2 - s1;
        printf("=> [%i^3] CPU: %f s\n", cnt, d1.count()); fflush(stdout);

        auto s3 = std::chrono::high_resolution_clock::now();
        dj::fft::arg<float> gpu = dj::fft::eval_3d_gpu(xi, dj::fft::e_dir::DIR_FWD);
        auto s4 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> d2 = s4 - s3;
        printf("=> [%i^3] GPU: %f s\n", cnt, d2.count()); fflush(stdout);

#if 0 // display content
        for (int k = 0; k < cnt; ++k) {
            for (int i = 0; i < cnt; ++i) {
                for (int j = 0; j < cnt; ++j) {
                    printf("{%f %f} vs {%f %f}\n",
                           cpu[i+cnt*(j+cnt*k)].real(), cpu[i+cnt*(j+cnt*k)].imag(),
                           gpu[i+cnt*(j+cnt*k)].real(), gpu[i+cnt*(j+cnt*k)].imag());
                 }
                printf("\n");
            }
            printf("\n\n");
        }
#endif
    }
#endif

    glfwTerminate();

    return 0;
}

