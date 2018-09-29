#include "image.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

/**
 * Image
 **/
Image::Image (int width_, int height_){

    assert(width_ > 0);
    assert(height_ > 0);

    width           = width_;
    height          = height_;
    num_pixels      = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;
    
    data.raw = new uint8_t[num_pixels*4];
    int b = 0; //which byte to write to
    for (int j = 0; j < height; j++){
        for (int i = 0; i < width; i++){
            data.raw[b++] = 0;
            data.raw[b++] = 0;
            data.raw[b++] = 0;
            data.raw[b++] = 0;
        }
    }

    assert(data.raw != NULL);
}

Image::Image (const Image& src){
    
    width           = src.width;
    height          = src.height;
    num_pixels      = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;
    
    data.raw = new uint8_t[num_pixels*4];
    
    //memcpy(data.raw, src.data.raw, num_pixels);
    *data.raw = *src.data.raw;
}

Image::Image (char* fname){

    int numComponents; //(e.g., Y, YA, RGB, or RGBA)
    data.raw = stbi_load(fname, &width, &height, &numComponents, 4);
    
    if (data.raw == NULL){
        printf("Error loading image: %s", fname);
        exit(-1);
    }
    

    num_pixels = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;
    
}

Image::~Image (){
    delete data.raw;
    data.raw = NULL;
}

void Image::Write(char* fname){
    
    int lastc = strlen(fname);

    switch (fname[lastc-1]){
       case 'g': //jpeg (or jpg) or png
         if (fname[lastc-2] == 'p' || fname[lastc-2] == 'e') //jpeg or jpg
            stbi_write_jpg(fname, width, height, 4, data.raw, 95);  //95% jpeg quality
         else //png
            stbi_write_png(fname, width, height, 4, data.raw, width*4);
         break;
       case 'a': //tga (targa)
         stbi_write_tga(fname, width, height, 4, data.raw);
         break;
       case 'p': //bmp
       default:
         stbi_write_bmp(fname, width, height, 4, data.raw);
    }
}

void Image::AddNoise (double factor)
{
    // Creating random noise that increases based on factor
    int x, y;
    for (x = 0; x < Width(); x++) {
        for (y = 0; y < Height(); y++) {
            Pixel p = GetPixel(x, y);
            p.r = rand() % 2 > 0 ? fmin(p.r + factor, 255) : fmax(0, fmin(p.r - factor, 255));
            p.g = rand() % 2 > 0 ? fmin(p.g + factor, 255) : fmax(0, fmin(p.g - factor, 255));
            p.b = rand() % 2 > 0 ? fmin(p.b + factor, 255) : fmax(0, fmin(p.b - factor, 255));
            GetPixel(x, y) = p;
        }
    }
}

void Image::Brighten (double factor)
{
    int x,y;
    for (x = 0; x < Width(); x++)
    {
        for (y = 0; y < Height(); y++)
        {
            Pixel p = GetPixel(x, y);
            Pixel scaled_p = p*factor;
            GetPixel(x,y) = scaled_p;
        }
    }
}


void Image::ChangeContrast (double factor)
{
    int r = 0;
    int g = 0;
    int b = 0;

    // Get the global greyscale
    for (int x = 0; x < Width(); x++) {
        for (int y = 0; y < Height(); y++) {
            Pixel p = GetPixel(x, y);
            r += p.r;
            g += p.g;
            b += p.b;
        }
    }

    r /= (Width() * Height());
    g /= (Width() * Height());
    b /= (Width() * Height());

    // Calculate global luminance
    int luminance = (int)(0.30*r + 0.59*g + 0.11*b);

    // Calculate new contrast based on global luminance
    for (int x = 0; x < Width(); x++) {
        for (int y = 0; y < Height(); y++) {
            Pixel p = GetPixel(x, y);

            int diff_r = (int)(p.r - luminance);
            int diff_g = (int)(p.g - luminance);
            int diff_b = (int)(p.b - luminance);

            p.r = fmax(0, fmin(255, (int)((p.r - diff_r) + (factor * diff_r))));
            p.g = fmax(0, fmin(255, (int)((p.g - diff_g) + (factor * diff_g))));
            p.b = fmax(0, fmin(255, (int)((p.b - diff_b) + (factor * diff_b))));

            GetPixel(x, y) = p;
        }
    }
}


void Image::ChangeSaturation(double factor)
{
    // Change saturation based on local greyscale
    for (int x = 0; x < Width(); x++) {
        for (int y = 0; y < Height(); y++) {
            Pixel p = GetPixel(x, y);
            
            int avg = (p.r + p.g + p.b) / 3.0;

            int diff_r = (int)(p.r - avg);
            int diff_g = (int)(p.g - avg);
            int diff_b = (int)(p.b - avg);

            p.r = fmax(0, fmin(255, (int)((p.r - diff_r) + (factor * diff_r))));
            p.g = fmax(0, fmin(255, (int)((p.g - diff_g) + (factor * diff_g))));
            p.b = fmax(0, fmin(255, (int)((p.b - diff_b) + (factor * diff_b))));

            GetPixel(x, y) = p;
        }
    }
}


Image* Image::Crop(int x, int y, int w, int h)
{
    if (x + w > Width() || y + h > Height() | x < 0 | y < 0 | w < 0 | h < 0) {
        return NULL;
    }

    Image* img = new Image(w, h);

    int max_x = x + w;
    int max_y = y + h;

    // Copy pixels over from same spot in original image
    int i = 0;
    for (x = max_x - w; x < max_x; x++) {
        int j = 0;
        for (y = max_y - h; y < max_y; y++) {
            Pixel p = GetPixel(x, y);
            img->GetPixel(i, j) = p;
            j++;
        }
        i++;
    }
    return img;
}


void Image::ExtractChannel(int channel)
{
    Pixel* p_ptr = new Pixel();

    // Switch between channels based on input
    switch (channel) {
    case 0:
        p_ptr = new Pixel(1, 0, 0, 1);
        break;
    case 1:
        p_ptr = new Pixel(0, 1, 0, 1);
        break;
    case 2:
        p_ptr = new Pixel(0, 0, 1, 1);
        break;
    default:
        p_ptr = new Pixel(1, 1, 1, 1);
        break;
    }

    for (int i = 0; i < Width(); i++) {
        for (int j = 0; j < Height(); j++) {
            Pixel p = GetPixel(i, j);
            Pixel new_p = Pixel(p.r * p_ptr->r, p.g * p_ptr->g, p.b * p_ptr->b, p.a * p_ptr->a);
            GetPixel(i, j) = new_p;
        }
    }

    delete p_ptr;
}


void Image::Quantize (int nbits)
{
    for (int x = 0; x < Width(); x++) {
        for (int y = 0; y < Height(); y++) {
            Pixel p = GetPixel(x, y);

            Pixel newP = PixelQuant(p, nbits);

            GetPixel(x, y) = newP;
        }
    }
}

void Image::RandomDither (int nbits)
{
    int factor = rand() % 8 + 1;
    AddNoise(4.0 * factor);

    Quantize(nbits);
}


static int Bayer4[4][4] =
{
    {15,  7, 13,  5},
    { 3, 11,  1,  9},
    {12,  4, 14,  6},
    { 0,  8,  2, 10}
};


void Image::OrderedDither(int nbits)
{
    /* WORK HERE */
}

/* Error-diffusion parameters */
const double
    ALPHA = 7.0 / 16.0,
    BETA  = 3.0 / 16.0,
    GAMMA = 5.0 / 16.0,
    DELTA = 1.0 / 16.0;

void Image::FloydSteinbergDither(int nbits)
{
    // Run through pixels of image
    for (int x = 0; x < Width(); x++) {
        for (int y = 0; y < Height(); y++) {
            Pixel p = GetPixel(x, y);

            Pixel newP = PixelQuant(p, nbits);

            // If we are not at a boundary, perform floydsteinbergdither
            if (x != 0 && x != Width() - 1 && y != Height() - 1) {
                // Get diff between quantitisized pixel and original
                int diffR = p.r - newP.r;
                int diffG = p.g - newP.g;
                int diffB = p.b - newP.b;

                // dither information to neighbors
                Pixel r = GetPixel(x + 1, y);
                Pixel rd = GetPixel(x + 1, y + 1);
                Pixel d = GetPixel(x, y + 1);
                Pixel ld = GetPixel(x - 1, y + 1);

                r.r += ALPHA * diffR;
                r.g += ALPHA * diffG;
                r.b += ALPHA * diffB;

                rd.r += DELTA * diffR;
                rd.g += DELTA * diffG;
                rd.b += DELTA * diffB;

                d.r += GAMMA * diffR;
                d.g += GAMMA * diffG;
                d.b += GAMMA * diffB;

                ld.r += BETA * diffR;
                ld.g += BETA * diffG;
                ld.b += BETA * diffB;


                GetPixel(x + 1, y) = r;
                GetPixel(x + 1, y + 1) = rd;
                GetPixel(x, y + 1) = d;
                GetPixel(x - 1, y + 1) = ld;
            }

            GetPixel(x, y) = newP;
        }
    }
}

void Image::Blur(int n)
{
    // Create our 2D array of pixels
    int size = 2 * n + 1;
    Pixel** pixels = new Pixel*[size];

    for (int i = 0; i < size; i++) {
        pixels[i] = new Pixel[size];
    }

    for (int x = 0; x < Width(); x++) {
        for (int y = 0; y < Height(); y++) {
            //Pixel p = GetPixel(x, y);
            //fprintf(stderr, "pr: %d, pg: %d, pb: %d\n", p.r, p.g, p.b);

            // Using mirror technique on edges
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    int r = x;
                    if (x > Width() - n - 1) {
                        r = Width() - n - 1 - abs(x - Width() - n - 1);
                    }

                    int c = y;
                    if (y > Height() - n - 1) {
                        c = Height() - n - 1 - abs(y - Height() - n - 1);
                    }

                    pixels[i][j] = GetPixel(abs(r - n + i), abs(c - n + j));
                }
            }

            int r = 0;
            int g = 0;
            int b = 0;
            double total = 0;

            // Dynamic sizing with size of square, all adds up to one
            for (int i = 0; i < size; i++) {
                int d_i = abs(i - n);
                for (int j = 0; j < size; j++) {
                    int d_j = abs(j - n);

                    r += pow(0.5, d_i + d_j) * pixels[i][j].r;
                    g += pow(0.5, d_i + d_j) * pixels[i][j].g;
                    b += pow(0.5, d_i + d_j) * pixels[i][j].b;

                    total += pow(0.5, d_i + d_j);
                }
            }

            r /= total;
            g /= total;
            b /= total;

            GetPixel(x, y) = Pixel(
                fmax(0, fmin(255, r)),
                fmax(0, fmin(255, g)),
                fmax(0, fmin(255, b)));
        }
    }

    for (int i = 0; i < size; i++) {
        delete pixels[i];
    }
    delete pixels;
}

void Image::Sharpen(int n)
{
    // Create new image that is a clone of the input
    Image *img = new Image(Width(), Height());

    for (int x = 0; x < Width(); x++) {
        for (int y = 0; y < Height(); y++) {
            img->GetPixel(x, y) = GetPixel(x, y);
        }
    }

    // Blur our input
    Blur(n);

    // Get the difference between the original image and the blur
    for (int x = 0; x < Width(); x++) {
        for (int y = 0; y < Height(); y++) {
            Pixel old_p = img->GetPixel(x, y);
            Pixel p = GetPixel(x, y);

            int rDiff = old_p.r - p.r;
            int gDiff = old_p.r - p.r;
            int bDiff = old_p.b - p.b;

            old_p.r = fmax(0, fmin(255, (int)(old_p.r + rDiff)));
            old_p.g = fmax(0, fmin(255, (int)(old_p.g + gDiff)));
            old_p.b = fmax(0, fmin(255, (int)(old_p.b + bDiff)));

            GetPixel(x, y) = old_p;
        }
    }
}

void Image::EdgeDetect()
{
    // Create our 2D array of pixels
    int extend = 1;
    int size = 2 * extend + 1;

    Pixel** pixels = new Pixel*[size];

    for (int i = 0; i < size; i++) {
        pixels[i] = new Pixel[size];
    }

    // Run through image with 3x3 filter
    for (int x = 0; x < Width(); x++) {
        for (int y = 0; y < Height(); y++) {

            // Set the pixels of our 3x3 array using mirror technique on edge
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {

                    int c = x;
                    if (x > Width() - extend - 1) {
                        c = Width() - extend - 1 - abs(x - Width() - extend - 1);
                    }

                    int r = y;
                    if (y > Height() - extend - 1) {
                        r = Height() - extend - 1 - abs(y - Height() - extend - 1);
                    }

                    pixels[i][j] = GetPixel(c + abs(i - extend), r + abs(j - extend));
                }
            }

            int r = 0;
            int g = 0;
            int b = 0;

            // Multiply current pixel by 8 and all neighbors by -1
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {

                    if (i == extend && j == extend) {
                        r += 8 * pixels[i][j].r;
                        g += 8 * pixels[i][j].g;
                        b += 8 * pixels[i][j].b;
                        continue;
                    }

                    r += -1 * pixels[i][j].r;
                    g += -1 * pixels[i][j].g;
                    b += -1 * pixels[i][j].b;
                }
            }

            // Set pixel to sum of edge detect convolution
            GetPixel(x, y) = Pixel(
                fmax(0, fmin(255, r)),
                fmax(0, fmin(255, g)),
                fmax(0, fmin(255, b)));
        }
    }

    for (int i = 0; i < size; i++) {
        delete pixels[i];
    }
    delete pixels;
}

Image* Image::Scale(double sx, double sy)
{
    /* WORK HERE */
    return NULL;
}

void Image::NonPhotorealism()
{
    char name[] = "sample_images/orangehaze.png";
    Image* img = new Image(name);

    for (int x = 0; x < Width(); x++) {
        for (int y = 0; y < Height(); y++) {
            Pixel p = GetPixel(x, y);

            int val = (int)((p.r + p.g + p.b) / 3.0);
            GetPixel(x, y) = img->GetPixel(val / 256.0 * img->Width(), 0);
        }
    }
    delete img;
}

Image* Image::Rotate(double angle)
{
    /* WORK HERE */
    return NULL;
}

void Image::Fun()
{
    /* WORK HERE */
}

Image* Image::Compression() {
    int scale = 2;

    Image* img = new Image(Width() / scale, Height() / scale);

    for (int x = 0; x+scale-1 < Width(); x+= scale) {
        for (int y = 0; y+scale-1 < Height(); y+= scale) {
            Pixel p = GetPixel(x, y);
            Pixel r = GetPixel(x + 1, y);
            Pixel d = GetPixel(x + 1, y + 1);
            Pixel l = GetPixel(x, y + 1);

            int red = (int)((1.0 / 4.0) * (p.r + r.r + d.r + l.r));
            int gre = (int)((1.0 / 4.0) * (p.g + r.g + d.g + l.g));
            int blu = (int)((1.0 / 4.0) * (p.b + r.b + d.b + l.b));

            img->GetPixel(x / scale, y / scale) = Pixel(
                fmax(0, fmin(255, red)),
                fmax(0, fmin(255, gre)),
                fmax(0, fmin(255, blu)));
        }
    }

    return img;
}

/**
 * Image Sample
 **/
void Image::SetSamplingMethod(int method)
{
    assert((method >= 0) && (method < IMAGE_N_SAMPLING_METHODS));
    sampling_method = method;
}


Pixel Image::Sample (double u, double v){
    /* WORK HERE */
    return Pixel();
}