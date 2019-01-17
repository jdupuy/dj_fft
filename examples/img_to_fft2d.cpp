
#include <cstdio>
#include <cstdint>
#include <cstdlib>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "dj_fft.h"

static void usage(const char *appName)
{
    printf("%s -- FFT for raster image data\n\n", appName);
    printf("usage: %s path_to_image_file\n", appName);
    printf("options:\n"
           "    --help\n"
           "        print help\n"
           "\n"
           "    --imageFormat can be either of:\n"
           "        u8  (--u8)  for uint8 image data file (default)\n"
           "        u16 (--u16) for uint16 image data file\n"
           "        hdr (--hdr) for fp32 image data file\n"
           "\n"
           "[NOTE: for square power-of-two resolution images only!]\n"
           );
}

bool is_power_of_two(int x)
{
    return (x & (x - 1)) == 0;
}

bool check_image(int x, int y)
{
    return (x == y) && is_power_of_two(x) && is_power_of_two(y);
}

int main(int argc, char **argv)
{
    enum {FORMAT_UINT8, FORMAT_UINT16, FORMAT_HDR};
    int format = FORMAT_UINT8;
    const char *imgFile = nullptr;
    dj::fft::arg<float> imgData, imgDataFFT;
    std::vector<float> fftNrm, fftArg;

    for (int i = 1; i < argc; ++i) {
        if (!strcmp("--hdr", argv[i])) {
            format = FORMAT_HDR;
        } else if (!strcmp("--u16", argv[i])) {
            format = FORMAT_UINT16;
        } else if (!strcmp("--u8", argv[i])) {
            format = FORMAT_UINT8;
        } else if (!strcmp("--help", argv[i])) {
            usage(argv[0]);
            return 0;
        } else {
            imgFile = argv[i];
        }
    }

    if (!imgFile) {
        usage(argv[0]);
        return 0;
    }

    // build input data
    if (format == FORMAT_HDR) {
        int x, y, comp;
        const float *texels = stbi_loadf(imgFile, &x, &y, &comp, 1);

        if (!check_image(x, y)) {
            printf("invalid image resolution\n");
            return 0;
        }

        imgData.resize(x * y);
        for (int j = 0; j < y; ++j)
        for (int i = 0; i < x; ++i)
            imgData[i + x * j] = texels[i + x * j];

        stbi_image_free((unsigned char *)texels);
    } else if (format == FORMAT_UINT16) {
        int x, y, comp;
        const uint16_t *texels = stbi_load_16(imgFile, &x, &y, &comp, 1);

        if (!check_image(x, y)) {
            printf("invalid image resolution\n");
            return 0;
        }

        imgData.resize(x * y);
        for (int j = 0; j < y; ++j)
        for (int i = 0; i < x; ++i)
            imgData[i + x * j] = float(texels[i + x * j]) / 64535.f;

        stbi_image_free((unsigned char *)texels);
    } else {
        int x, y, comp;
        const uint8_t *texels = stbi_load(imgFile, &x, &y, &comp, 1);

        if (!check_image(x, y)) {
            printf("invalid image resolution\n");
            return 0;
        }

        imgData.resize(x * y);
        for (int j = 0; j < y; ++j)
        for (int i = 0; i < x; ++i)
            imgData[i + x * j] = float(texels[i + x * j]) / 255.f;

        stbi_image_free((unsigned char *)texels);
    }

    // compute FFT:
    imgDataFFT = dj::fft::eval_2d(imgData, dj::fft::e_dir::DIR_FWD);

    // build norm array and phase array
    fftNrm.resize(imgDataFFT.size());
    fftArg.resize(imgDataFFT.size());

    for (int i = 0; i < (int)imgDataFFT.size(); ++i) {
        fftNrm[i] = std::abs(imgDataFFT[i]);
        fftArg[i] = std::arg(imgDataFFT[i]);
        if (fftArg[i] < 0.0f) fftArg[i]+= /*Pi*/std::acos(-1.f);
    }

    // write to .hdr file
    int x = sqrt(fftNrm.size());
    stbi_write_hdr("fft_nrm.hdr", x, x, 1, &fftNrm[0]);
    stbi_write_hdr("fft_arg.hdr", x, x, 1, &fftArg[0]);

    return 0;
}
