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

void Image::CopyPixels(Image* img) {
    for (int i = 0; i < Width(); i++) {
        for (int j = 0; j < Height(); j++) {
            GetPixel(i, j) = img->GetPixel(i, j);
        }
    }
}

void Image::AddNoise (double factor)
{
    // Creating random noise that increases based on factor
    int x, y;
    for (x = 0; x < Width(); x++) {
        for (y = 0; y < Height(); y++) {
            Pixel p = GetPixel(x, y);
            p.r = rand() % 2 > 0 ? fmin(p.r + (rand() % (int)factor), 255) : fmax(0, fmin(p.r - (rand() % (int)factor), 255));
            p.g = rand() % 2 > 0 ? fmin(p.g + (rand() % (int)factor), 255) : fmax(0, fmin(p.g - (rand() % (int)factor), 255));
            p.b = rand() % 2 > 0 ? fmin(p.b + (rand() % (int)factor), 255) : fmax(0, fmin(p.b - (rand() % (int)factor), 255));
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
    int factor = rand() % 4 + 5;
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

void Image::FloydSteinbergDither(int nbits) {
    // Run through pixels of image
    for (int x = 0; x < Width(); x++) {
        for (int y = 0; y < Height(); y++) {
            Pixel p = GetPixel(x, y);

            Pixel newP = PixelQuant(p, nbits);

            // If we are not at a boundary, perform floydsteinbergdither
            if (x != 0 && x != Width() - 1 && y != Height() - 1) {
                // Get diff between quantisized pixel and original
                int diffR = p.r - newP.r;
                int diffG = p.g - newP.g;
                int diffB = p.b - newP.b;

                // dither information to neighbors
                Pixel r = GetPixel(x + 1, y);
                Pixel rd = GetPixel(x + 1, y + 1);
                Pixel d = GetPixel(x, y + 1);
                Pixel ld = GetPixel(x - 1, y + 1);

                r.r = fmax(0, fmin(255, r.r + ALPHA * diffR));
                r.g = fmax(0, fmin(255, r.g + ALPHA * diffG));
                r.b = fmax(0, fmin(255, r.b + ALPHA * diffB));

                rd.r = fmax(0, fmin(255, rd.r + DELTA * diffR));
                rd.g = fmax(0, fmin(255, rd.g + DELTA * diffG));
                rd.b = fmax(0, fmin(255, rd.b + DELTA * diffB));

                d.r = fmax(0, fmin(255, d.r + GAMMA * diffR));
                d.g = fmax(0, fmin(255, d.g + GAMMA * diffG));
                d.b = fmax(0, fmin(255, d.b + GAMMA * diffB));

                ld.r = fmax(0, fmin(255, ld.r + BETA * diffR));
                ld.g = fmax(0, fmin(255, ld.g + BETA * diffG));
                ld.b = fmax(0, fmin(255, ld.b + BETA * diffB));


                GetPixel(x + 1, y) = r;
                GetPixel(x + 1, y + 1) = rd;
                GetPixel(x, y + 1) = d;
                GetPixel(x - 1, y + 1) = ld;
            }

            GetPixel(x, y) = newP;
        }
    }
}

void Image::Blur(int n) {
    // Create our 2D array of pixels
    int size = 2 * n + 1;
    Pixel** pixels = new Pixel*[size];

    Image* img = new Image(Width(), Height());

    for (int i = 0; i < size; i++) {
        pixels[i] = new Pixel[size];
    }

    for (int x = 0; x < Width(); x++) {
        for (int y = 0; y < Height(); y++) {

            // Using mirror technique on edges
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {

                    int r = abs(x);
                    if (x > Width() - n - 1) {
                        r = abs(Width() - n - 1 - abs(x - Width() - n - 1));
                    }

                    int c = abs(y);
                    if (y > Height() - n - 1) {
                        c = abs(Height() - n - 1 - abs(y - Height() - n - 1));
                    }

                    pixels[i][j] = GetPixel(r + abs(i - n), c + abs(j - n));
                }
            }

            double r = 0;
            double g = 0;
            double b = 0;
            double total = 0;

            // Dynamic sizing with size of square, all adds up to one
            for (int i = 0; i < size; i++) {
                int d_i = abs(i - n);
                for (int j = 0; j < size; j++) {
                    int d_j = abs(j - n);

                    r += pow(1.0 / n, pow(d_i, 2) + pow(d_j, 2)) * pixels[i][j].r;
                    g += pow(1.0 / n, pow(d_i, 2) + pow(d_j, 2)) * pixels[i][j].g;
                    b += pow(1.0 / n, pow(d_i, 2) + pow(d_j, 2)) * pixels[i][j].b;

                    total += pow(1.0 / n, pow(d_i, 2) + pow(d_j, 2));
                }
            }

            r /= total;
            g /= total;
            b /= total;

            img->GetPixel(x, y) = Pixel(
                fmax(0, fmin(255, (int)r)),
                fmax(0, fmin(255, (int)g)),
                fmax(0, fmin(255, (int)b)));
        }
    }

    CopyPixels(img);

    delete img;

    for (int i = 0; i < size; i++) {
        delete pixels[i];
    }
    delete pixels;
}

void Image::Sharpen(int n) {
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

void Image::EdgeDetect() {
    // Create our 2D array of pixels
    int extend = 1;
    int size = 2 * extend + 1;

    Image* img = new Image(Width(), Height());

    Pixel** pixels = new Pixel*[size];

    for (int i = 0; i < size; i++) {
        pixels[i] = new Pixel[size];
    }

    ChangeSaturation(0);
    ChangeContrast(1.5);

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
            img->GetPixel(x, y) = Pixel(
                fmax(0, fmin(255, r)),
                fmax(0, fmin(255, g)),
                fmax(0, fmin(255, b)));
        }
    }

    CopyPixels(img);

    delete img;

    for (int i = 0; i < size; i++) {
        delete pixels[i];
    }
    delete pixels;
}

Image* Image::Scale(double sx, double sy) {
    int img_width = (int)(Width() * sx);
    int img_height = (int)(Height() * sy);

    Image* img = new Image(img_width, img_height);
    
    for (int x = 0; x < img_width; x++) {
        for (int y = 0; y < img_height; y++) {
            img->GetPixel(x, y) = Sample(x / sx, y / sy);
        }
    }
    return img;
}

void Image::NonPhotorealism() {
    char name[] = "sample_images/darkpinkhaze.png";
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

Image* Image::Rotate(double angle) {
    int new_width =
        fmax((Width() - 1) * cos(-angle) - (Height() - 1) * sin(-angle), (Width() - 1) * cos(-angle) - 0 * sin(-angle)) -
        fmin(0 * cos(-angle) - 0 * sin(-angle), 0 * cos(-angle) - (Height() - 1) * sin(-angle));

    int new_height =
        fmax(0 * sin(-angle) + (Height() - 1) * cos(-angle), (Width() - 1) * sin(-angle) + (Height() - 1) * cos(-angle)) -
        fmin(0 * sin(-angle) + 0 * cos(-angle), (Width() - 1) * sin(-angle) + 0 * cos(-angle));

    Image* img = new Image(new_width, new_height);

    for (int x = 0; x < new_width; x++) {
        for (int y = 0; y < new_height; y++) {
            float u = (Width() / 2.0) + (x - new_width / 2.0) * cos(-angle) - (y - new_height / 2.0) * sin(-angle);
            float v = (Height() / 2.0) + (x - new_width / 2.0) * sin(-angle) + (y - new_height / 2.0) * cos(-angle);

            if (u >= 0 && v >= 0 && u < Width() && v < Height())
                img->GetPixel(x, y) = Sample(u, v);
            else
                img->GetPixel(x, y) = Pixel(0, 0, 0);
        }
    }

    return img;
}

void Image::Halftone() {
    Image* img = new Image(Width(), Height());

    int i = -1;

    for (int x = 0; x < Width() - 1; x += 2) {
        for (int y = 0; y < Height() - 1; y += 2) {
            i++;

            if (i % 2 == 0) {
                img->GetPixel(x, y) = Pixel(255, 255, 255);
                img->GetPixel(x + 1, y) = Pixel(255, 255, 255);
                img->GetPixel(x + 1, y + 1) = Pixel(255, 255, 255);
                img->GetPixel(x, y + 1) = Pixel(255, 255, 255);

                continue;
            }

            Pixel p = GetPixel(x, y);
            Pixel r = GetPixel(x + 1, y);
            Pixel rd = GetPixel(x + 1, y + 1);
            Pixel d = GetPixel(x, y + 1);

            double red = (p.r + r.r + rd.r + d.r) / 4.0;
            double gre = (p.g + r.g + rd.g + d.g) / 4.0;
            double blu = (p.b + r.b + rd.b + d.b) / 4.0;

            double avg = (red + gre + blu) / 3.0;

            if (avg / 255.0 < 0.2) {
                img->GetPixel(x, y) = Pixel(255, 255, 255);
                img->GetPixel(x + 1, y) = Pixel(255, 255, 255);
                img->GetPixel(x + 1, y + 1) = Pixel(255, 255, 255);
                img->GetPixel(x, y + 1) = Pixel(255, 255, 255);
            }
            else if (avg / 255.0 < 0.4) {
                img->GetPixel(x, y) = Pixel(0, 0, 0);
                img->GetPixel(x + 1, y) = Pixel(255, 255, 255);
                img->GetPixel(x + 1, y + 1) = Pixel(255, 255, 255);
                img->GetPixel(x, y + 1) = Pixel(255, 255, 255);
            }
            else if (avg / 255.0 < 0.6) {
                img->GetPixel(x, y) = Pixel(0, 0, 0);
                img->GetPixel(x + 1, y) = Pixel(255, 255, 255);
                img->GetPixel(x + 1, y + 1) = Pixel(0, 0, 0);
                img->GetPixel(x, y + 1) = Pixel(255, 255, 255);
            }
            else if (avg / 255.0 < 0.8) {
                img->GetPixel(x, y) = Pixel(0, 0, 0);
                img->GetPixel(x + 1, y) = Pixel(0, 0, 0);
                img->GetPixel(x + 1, y + 1) = Pixel(0, 0, 0);
                img->GetPixel(x, y + 1) = Pixel(255, 255, 255);
            }
            else {
                img->GetPixel(x, y) = Pixel(0, 0, 0);
                img->GetPixel(x + 1, y) = Pixel(0, 0, 0);
                img->GetPixel(x + 1, y + 1) = Pixel(0, 0, 0);
                img->GetPixel(x, y + 1) = Pixel(0, 0, 0);
            }
        }
    }

    CopyPixels(img);

    delete img;
}

void Image::Fun() {

    Image* img = new Image(Width(), Height());
    double r = fmin(Width(), Height()) / 2.0;

    for (int x = 0; x < Width(); x++) {
        for (int y = 0; y < Height(); y++) {
            double center_x = x - (Width() / 2.0);
            double center_y = y - (Height() / 2.0);

            double theta = atan2(center_x, center_y) - M_PI / 2.0;

            double d = Distance(0, 0, center_x, center_y) / r;
            double pol_x = r * pow(d, 2) * cos(-theta);
            double pol_y = r * pow(d, 2) * sin(-theta);

            float u = pol_x + Width() / 2.0;
            float v = pol_y + Height() / 2.0;

            if (pol_x * pol_x + pol_y * pol_y < r * r)
                img->GetPixel(x, y) = Sample(u, v);
            else
               img->GetPixel(x, y) = Pixel(0, 0, 0);
        }
    }

    CopyPixels(img);

    delete img;
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


Pixel Image::Sample (double u, double v) {

    int r = 0;
    int g = 0;
    int b = 0;

    if (u < 0)
        return Pixel(0, 0, 0);
    if (v < 0)
        return Pixel(0, 0, 0);
    if (u > Width() - 1)
        return Pixel(0, 0, 0);
    if (v > Height() - 1)
        return Pixel(0, 0, 0);

    switch (sampling_method) {

        // Closest point
    case IMAGE_SAMPLING_POINT:

        return GetPixel((int)u, (int)v);

        // Bilinear (4 closest)
    case IMAGE_SAMPLING_BILINEAR: {

        double d1 = 1.0 / Distance(u, v, floor(u), floor(v)) < INFINITY ? 1.0 / Distance(u, v, floor(u), floor(v)) : 1.0;
        double d2 = 1.0 / Distance(u, v, floor(u), ceil(v)) < INFINITY ? 1.0 / Distance(u, v, floor(u), ceil(v)) : 1.0;
        double d3 = 1.0 / Distance(u, v, ceil(u), floor(v)) < INFINITY ? 1.0 / Distance(u, v, ceil(u), floor(v)) : 1.0;
        double d4 = 1.0 / Distance(u, v, ceil(u), ceil(v)) < INFINITY ? 1.0 / Distance(u, v, ceil(u), ceil(v)) : 1.0;

        if (ceil(u) < Width() && ceil(v) < Height()) {
            double sum = d1 + d2 + d3 + d4;

            //fprintf(stderr, "sum: %f\n", sum);
            //fprintf(stderr, "d1: %f, d2: %f, d3: %f, d4: %f\n", d1 / sum, d2 / sum, d3 / sum, d4 / sum);
            //fprintf(stderr, "flooru: %f, ceilu: %f, floorv: %f, ceilb: %f\n", floor(u), ceil(u), floor(v), ceil(v));

            Pixel p = 
                GetPixel(floor(u), floor(v)) * (d1 / sum) +
                GetPixel(floor(u), ceil(v)) * (d2 / sum) +
                GetPixel(ceil(u), floor(v)) * (d3 / sum) +
                GetPixel(ceil(u), ceil(v)) * (d4 / sum);

            //fprintf(stderr, "p.r: %f, p.g: %f, p.b: %f\n", p.r, p.g, p.b);
            return p;
        }
        else if (ceil(u) == Width() && ceil(v) != Height()) {
            double sum = d1 + d2;

            return
                GetPixel(floor(u), floor(v)) * (d1 / sum) +
                GetPixel(floor(u), ceil(v)) * (d2 / sum);
        }
        else if (ceil(u) != Width() && ceil(v) == Height()) {
            double sum = d1 + d3;

            return
                GetPixel(floor(u), floor(v)) * (d1 / sum) +
                GetPixel(ceil(u), floor(v)) * (d3 / sum);
        }
        else {
            return GetPixel(floor(u), floor(v));
        }
    }

        // Gaussian
    case IMAGE_SAMPLING_GAUSSIAN: {

        double red = 0;
        double gre = 0;
        double blu = 0;

        double weightTot = 0;

        for (int x = (int)u - 2; x <= (int)u + 2; x++) {
            for (int y = (int)v - 2; y <= (int)v + 2; y++) {

                int r = abs(x);
                if (x >= Width())
                    r = abs(Width() - 1 - abs(x - Width() - 1));

                int c = abs(y);
                if (y >= Height())
                    c = abs(Height() - 1 - abs(y - Height() - 1));

                double weight = (1.0 / (2 * M_PI)) * pow(M_E, -pow(Distance(u, v, r, c), 2) / 2.0);

                red += weight * GetPixel(r, c).r;
                gre += weight * GetPixel(r, c).g;
                blu += weight * GetPixel(r, c).b;

                weightTot += weight;
            }
        }

        return Pixel(
            fmax(0, fmin(255, (int)(red / weightTot))),
            fmax(0, fmin(255, (int)(gre / weightTot))),
            fmax(0, fmin(255, (int)(blu / weightTot))));

    }

        // Shouldn't reach here based on assert in setmethod
    default:
        return Pixel();
    }
}

double Image::Distance(double x, double y, double x1, double y1) {
    return sqrt(pow(x - x1, 2) + pow(y - y1, 2));
}