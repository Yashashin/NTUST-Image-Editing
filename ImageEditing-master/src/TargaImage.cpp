///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#include "Globals.h"
#include "TargaImage.h"
#include "libtarga.h"
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include<ctime>

using namespace std;

// constants
const int           RED             = 0;                // red channel
const int           GREEN           = 1;                // green channel
const int           BLUE            = 2;                // blue channel
const unsigned char BACKGROUND[3]   = { 0, 0, 0 };      // background color


// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1 ; i <= s ; i++)
        res = (n - i + 1) * res / i ;

    return res;
}// Binomial


///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
   data = new unsigned char[width * height * 4];
   ClearToBlack();
}// TargaImage



///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char *d)
{
    int i;

    width = w;
    height = h;
    data = new unsigned char[width * height * 4];

    for (i = 0; i < width * height * 4; i++)
	    data[i] = d[i];
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image) 
{
   width = image.width;
   height = image.height;
   data = NULL; 
   if (image.data != NULL) {
      data = new unsigned char[width * height * 4];
      memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
    if (data)
        delete[] data;
}// ~TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
    unsigned char   *rgb = new unsigned char[width * height * 3];
    int		    i, j;

    if (! data)
	    return NULL;

    // Divide out the alpha
    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = i * width * 4;
	    int out_offset = i * width * 3;

	    for (j = 0 ; j < width ; j++)
        {
	        RGBA_To_RGB(data + (in_offset + j*4), rgb + (out_offset + j*3));
	    }
    }

    return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char *filename)
{
    TargaImage	*out_image = Reverse_Rows();

    if (! out_image)
	    return false;

    if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
    {
	    cout << "TGA Save Error: %s\n", tga_error_string(tga_get_last_error());
	    return false;
    }

    delete out_image;

    return true;
}// Save_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char *filename)
{
    unsigned char   *temp_data;
    TargaImage	    *temp_image;
    TargaImage	    *result;
    int		        width, height;

    if (!filename)
    {
        cout << "No filename given." << endl;
        return NULL;
    }// if

    temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
    if (!temp_data)
    {
        cout << "TGA Error: %s\n", tga_error_string(tga_get_last_error());
	    width = height = 0;
	    return NULL;
    }
    temp_image = new TargaImage(width, height, temp_data);
    free(temp_data);

    result = temp_image->Reverse_Rows();

    delete temp_image;

    return result;
}// Load_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
   
       for (int i = 0; i < this->height; i++)
       {
           for (int j = 0; j < this->width; j++)
           {
               double r, g, b;
               r = this->data[(i * this->width + j) * 4];
               g = this->data[(i * this->width + j) * 4 + 1];
               b = this->data[(i * this->width + j) * 4 + 2];
               r = r * 0.30;
               g = g * 0.59;
               b = b * 0.11;
               unsigned char tmp = r + g + b;
               this->data[(i * this->width + j) * 4] = tmp;
               this->data[(i * this->width + j) * 4 + 1] = tmp;
               this->data[(i * this->width + j) * 4 + 2] = tmp;

           }
       }
	return true;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
    const double r1 = 256 / 8;
	const double r2 = 256 / 4;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
            double r, g, b;
            r = this->data[(i * this->width + j) * 4];
            g = this->data[(i * this->width + j) * 4 + 1];
            b = this->data[(i * this->width + j) * 4 + 2];
          
            r = r / r1;
            g = g / r1;
            b = b / r2;
           
            r = (int)r * 32.0;
            g = (int)g * 32.0;
            b = (int)b * 64.0;
            this->data[(i * this->width + j) * 4] = (char)r;
            this->data[(i * this->width + j) * 4 + 1] = (char)g;
            this->data[(i * this->width + j) * 4 + 2] = (char)b;
		}
	}
    return true;
}// Quant_Uniform

int  computeDistance(double r1, double g1, double b1, double r2, double g2, double b2)
{
    int result = (r1 - r2) * (r1 - r2) + (g1 - g2) * (g1 - g2) + (b1 - b2) * (b1 - b2);
    return result;
}


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity()
{
    this->Quant_Uniform();
    vector<vector<int>> record;
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            bool flag = false;
            for (int k = 0; k < record.size(); k++)
            {
                if (record[k][0] == data[(width * i + j) * 4] && record[k][1] == data[(width * i + j) * 4 + 1] && record[k][2] == data[(width * i + j) * 4+2]) //same color
                {
                    record[k][3]++;
                    flag = true;
                    break;
                }
            }
            if (flag == false)
            {
                vector<int>arr;
                arr.push_back(data[(width * i + j) * 4]);
                arr.push_back(data[(width * i + j) * 4+1]);
                arr.push_back(data[(width * i + j) * 4+2]);
                arr.push_back(1);
                record.push_back(arr);
            }
        }
    }
    if (record.size() < 256)
    {
        cout << "less than 256" << endl;
        return false;
    }

    for (int i = 1; i < record.size(); i++)
    {
        int keyR = record[i][0];
        int keyG = record[i][1];
        int keyB = record[i][2];
        int key = record[i][3];
        int j = i - 1;
        while (j >= 0 && key > record[j][3])
        {
            record[j + 1][0] = record[j][0];
            record[j + 1][1] = record[j][1];
            record[j + 1][2] = record[j][2];
            record[j + 1][3] = record[j][3];
            j--;
        }
        record[j + 1][0] = keyR;
        record[j + 1][1] = keyG;
        record[j + 1][2] = keyB;
        record[j + 1][3] = key;
    }
    
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int max = 10000000;
            int index = -1;
            for (int k = 0; k < 256; k++)
            {
                int dis = computeDistance((double)data[(width * i + j) * 4], (double)data[(width * i + j) * 4 + 1], (double)data[(width * i + j) * 4 + 2], (double)record[k][0], (double)record[k][1], (double)record[k][2]);
                if (max >= dis)
                {
                    max = dis;
                    index = k;
                }
            }
            data[(width * i + j) * 4] = record[index][0];
            data[(width * i + j) * 4+1] = record[index][1];
            data[(width * i + j) * 4+2] = record[index][2];
        }
    }
    return true;
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
    this->To_Grayscale();
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            if (data[(i * width + j) * 4] / 255.0 > 0.5)
            {
                data[(i * width + j) * 4] = 255;
                data[(i * width + j) * 4+1] = 255;
                data[(i * width + j) * 4 + 2] = 255;
            }
            else
            {
                data[(i * width + j) * 4] = 0;
                data[(i * width + j) * 4 + 1] = 0;
                data[(i * width + j) * 4 + 2] = 0;
            }
        }
    }
    return true;
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
    srand(time(NULL));
    this->To_Grayscale();
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width;j++)
        {
            int random = rand() % 20;
            double deviation = -0.2 + 0.02 * random;
            if (data[(i * width + j) * 4] / 255.0 + deviation > 0.5)
            {
                data[(i * width + j) * 4] = 255;
                data[(i * width + j) * 4 + 1] = 255;
                data[(i * width + j) * 4 + 2] = 255;
            }
            else
            {
                data[(i * width + j) * 4] = 0;
                data[(i * width + j) * 4 + 1] = 0;
                data[(i * width + j) * 4 + 2] = 0;
            }
        }
    }
    return true;
}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_FS()
{
    //double a = 255;
    //data[0] = a;
    //if (data[0] - 256.0 <0)
    //{
    //    cout << "jhaha";
    //}
	this->To_Grayscale();
	for (int i = 0; i < height; i++)
	{
		if (i % 2 == 0)
		{
			for (int j = 0; j < width; j++)
			{
				double error;
				if (data[(i * width + j) * 4] / 255.0 >= 0.5)
				{
					error = (data[(i * width + j) * 4] - 255);
					data[(i * width + j) * 4] = 255;
					data[(i * width + j) * 4 + 1] = 255;
					data[(i * width + j) * 4 + 2] = 255;
				}
				else
				{
					error = data[(i * width + j) * 4];
					data[(i * width + j) * 4] = 0;
					data[(i * width + j) * 4 + 1] = 0;
					data[(i * width + j) * 4 + 2] = 0;
				}
    
				if (j + 1 < width)
				{
                    if (data[(i * width + j + 1) * 4] + error * 7 / 16.0 > 255)
                    {
                        data[(i * width + j + 1) * 4] = 255;
                    }
                    else if (data[(i * width + j + 1) * 4] + error * 7 / 16.0 < 0)
                    {
                        data[(i * width + j + 1) * 4] = 0;
                    }
                    else
                    {
                        data[(i * width + j + 1) * 4] += error * 7 / 16.0;
                    }
                   
				}
				if (j - 1 >= 0 && i + 1 < height)
				{
                    if (data[((i + 1) * width + j - 1) * 4] + error *3 / 16.0 > 255)
                    {
                        data[((i + 1) * width + j - 1) * 4] = 255;
                    }
                    else if (data[((i + 1) * width + j - 1) * 4] + error * 3 / 16.0 < 0)
                    {
                        data[((i + 1) * width + j - 1) * 4] = 0;
                    }
                    else
                    {
                        data[((i + 1) * width + j - 1) * 4] += error * 3 / 16.0;
                    }
				}
				if (i + 1 < height)
				{
                    if (data[((i + 1) * width + j) * 4] + error * 5 / 16.0 > 255)
                    {
                        data[((i + 1) * width + j) * 4] = 255;
                    }
                    else if (data[((i + 1) * width + j) * 4] + error *5 / 16.0 < 0)
                    {
                        data[((i + 1) * width + j) * 4] = 0;
                    }
                    else
                    {
                        data[((i + 1) * width + j) * 4] += error * 5 / 16.0;
                    }
				}
				if (j + 1 < width && i + 1 < height)
				{
                    if (data[((i + 1) * width + j + 1) * 4] + error * 1 / 16.0 > 255)
                    {
                        data[((i + 1) * width + j + 1) * 4] = 255;
                    }
                    else if (data[((i + 1) * width + j + 1) * 4] + error * 1 / 16.0 < 0)
                    {
                        data[((i + 1) * width + j + 1) * 4] = 0;
                    }
                    else
                    {
                        data[((i + 1) * width + j + 1) * 4] += error * 1 / 16.0;
                    }
				}
			}
		}
		else
		{
            for (int j = width - 1; j >= 0; j--)
            {
                double error;
                if (data[(i * width + j) * 4] / 255.0 >= 0.5)
                {
                    error = (data[(i * width + j) * 4] - 255);
                    data[(i * width + j) * 4] = 255;
                    data[(i * width + j) * 4 + 1] = 255;
                    data[(i * width + j) * 4 + 2] = 255;
                }
                else
                {
                    error = data[(i * width + j) * 4];
                    data[(i * width + j) * 4] = 0;
                    data[(i * width + j) * 4 + 1] = 0;
                    data[(i * width + j) * 4 + 2] = 0;
                }
           
                if (j - 1 >= 0)
                {
                   
                    if (data[(i * width + j - 1) * 4] + error*7/ 16.0> 255)
                    {
                        data[(i * width + j - 1) * 4] = 255;
                    }
                    else if (data[(i * width + j - 1) * 4] + error * 7 / 16.0 < 0)
                    {
                        data[(i * width + j - 1) * 4] = 0;
                    }
                    else
                    {
                        data[(i * width + j - 1) * 4] += error * 7 / 16.0;
                    }
                }
                if (j + 1 <width && i + 1 < height)
                {
                    if (data[((i + 1) * width + j + 1) * 4] + error * 3 / 16.0 > 255)
                    {
                        data[((i + 1) * width + j + 1) * 4] = 255;
                    }
                    else if (data[((i + 1) * width + j + 1) * 4] + error * 3 / 16.0 < 0)
                    {
                        data[((i + 1) * width + j + 1) * 4] = 0;
                    }
                    else
                    {
                        data[((i + 1) * width + j + 1) * 4] += error * 3 / 16.0;
                    }
                }
                if (i + 1 < height)
                {
                    if (data[((i + 1) * width + j) * 4] + error * 5 / 16.0 > 255)
                    {
                        data[((i + 1) * width + j) * 4] = 255;
                    }
                    else if (data[((i + 1) * width + j) * 4] + error * 5 / 16.0 < 0)
                    {
                        data[((i + 1) * width + j) * 4] = 0;
                    }
                    else
                    {
                        data[((i + 1) * width + j) * 4] += error * 5 / 16.0;
                    }
                }
                if (j -1 >=0 && i + 1 < height)
                {
                    if (data[((i + 1) * width + j - 1) * 4] + error * 1 / 16.0 > 255)
                    {
                        data[((i + 1) * width + j - 1) * 4] = 255;
                    }
                    else if (data[((i + 1) * width + j - 1) * 4] + error * 1 / 16.0< 0)
                    {
                        data[((i + 1) * width + j - 1) * 4] = 0;
                    }
                    else
                    {
                        data[((i + 1) * width + j - 1) * 4] += error * 1 / 16.0;
                    }
                }
            }
		}
	}
return true;
}// Dither_FS



void insertSort(double *arr,int size)
{
    for (int i = 1; i < size; i++)
    {
        double key = arr[i * 2];
        double key2 = arr[i * 2+1];
        int j = i - 1;
        while (j >= 0 && key <=arr[j * 2])
        {
            arr[(j+1) * 2] = arr[j * 2];
            arr[(j+1) * 2+1] = arr[j * 2+1];
            j--;
        }
        arr[(j + 1) * 2] = key;
        arr[(j + 1) * 2 + 1] = key2;
    }
}

///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
    this->To_Grayscale();
    double sum = 0;
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            sum += (data[(i * width + j) * 4])/255.0;
        }
    }
    double threshold = sum / (width * height);

    double *arr = new  double[width * height*2];
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
           arr[(i * width + j)*2] = (double)data[(i * width + j) * 4];
           arr[(i * width + j)*2+1] = (double)(i*width+j);
        }
    }

    insertSort(arr,height * width);

    for (int i = 0; i < (1-threshold)*width*height; i++)
    {
        data[(int)(arr[i * 2 + 1] * 4)] = 0;
        data[(int)(arr[i * 2 + 1] * 4 + 1)] = 0;
        data[(int)(arr[i * 2 + 1] * 4 + 2)] = 0;
    }
    for (int i = (1-threshold)*width*height; i < width*height; i++)
    {
        data[(int)(arr[i * 2 + 1] * 4)] = 255;
        data[(int)(arr[i * 2 + 1] * 4 + 1)] = 255;
        data[(int)(arr[i * 2 + 1] * 4 + 2)] = 255;
    }
    return true;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
    const double threshold[4][4] = {
        {0.7059,0.3529,0.5882,0.2353},
        {0.0588,0.9412,0.8325,0.4118},
        {0.4706,0.7647,0.8824,0.1176},
        {0.1765,0.5294,0.2941,0.6471}
    };
    this->To_Grayscale();
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            char value;
            if ((double)data[(i * width + j)*4] >= threshold[i % 4][j % 4]*256)
            {
               value = (char)255;
            }
            else
            {
                value = (char)0;
            }
            this->data[(i * this->width + j) * 4] = value;
            this->data[(i * this->width + j) * 4 + 1] = value;
            this->data[(i * this->width + j) * 4 + 2] = value;
        }
    }
    return true;
}// Dither_Cluster


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
    for (int i = 0; i < height; i++)
    {
        if (i % 2 == 0)
        {
            for (int j = 0; j < width; j++)
            {
				double errorR, errorG, errorB;
				double newR = data[(i * width + j) * 4] / 32.0;
                newR = (int)newR;
                newR *= 32;
				double newG = data[(i * width + j) * 4+1] / 32.0;
                newG = (int)newG;
                newG *= 32;
				double newB = data[(i * width + j) * 4+2] / 64.0;
				newB =(int)newB;
                newB *= 64;
				errorR = data[(i * width + j) * 4] - newR;
				errorG = data[(i * width + j) * 4 + 1] - newG;
				errorB = data[(i * width + j) * 4 + 2] - newB;

                data[(i * width + j) * 4]= newR;
                data[(i * width + j) * 4 + 1] = newG;
                data[(i * width + j) * 4 + 2] = newB;

                if (j + 1 < width)
                {
                    if (data[(i * width + j + 1) * 4] + errorR * (7 / 16.0) > 255)
                    {
                        data[(i * width + j + 1) * 4] = 255;
                    }
                    else if (data[(i * width + j + 1) * 4] + errorR * (7 / 16.0) < 0)
                    {
                        data[(i * width + j + 1) * 4] = 0;
                    }
                    else
                    {
                        data[(i * width + j + 1) * 4] += errorR * (7 / 16.0);
                    }
                    if (data[(i * width + j + 1) * 4 + 1] +errorG * (7 / 16.0) > 255)
                    {
                        data[(i * width + j + 1) * 4+1] = 255;
                    }
                    else if (data[(i * width + j + 1) * 4 + 1] + errorG * (7 / 16.0) < 0)
                    {
                        data[(i * width + j + 1) * 4 + 1] = 0;
                    }
                    else
                    {
                        data[(i * width + j + 1) * 4 + 1] += errorG * (7 / 16.0);
                    }
                    if (data[(i * width + j + 1) * 4 + 2] +errorB * (7 / 16.0) > 255)
                    {
                        data[(i * width + j + 1) * 4 + 2] = 255;
                    }
                    else if (data[(i * width + j + 1) * 4+2] + errorB * (7 / 16.0) < 0)
                    {
                        data[(i * width + j + 1) * 4+2] = 0;
                    }
                    else
                    {
                        data[(i * width + j + 1) * 4 + 2] += errorB * (7 / 16.0);
                    }
                }
                if (j - 1 >= 0 && i + 1 < height)
                {
                    if (data[((i + 1) * width + j - 1) * 4] +errorR * (3 / 16.0) > 255)
                    {
                        data[((i + 1) * width + j - 1) * 4] = 255;
                    }
                    else if (data[((i + 1) * width + j - 1) * 4] + errorR * (3 / 16.0) < 0)
                    {
                        data[((i + 1) * width + j - 1) * 4] = 0;
                    }
                    else
                    {
                        data[((i + 1) * width + j - 1) * 4] += errorR * (3 / 16.0);
                    }
                    if (data[((i + 1) * width + j - 1) * 4+1] + errorG * (3 / 16.0) > 255)
                    {
                        data[((i + 1) * width + j - 1) * 4+1] = 255;
                    }
                    else if (data[((i + 1) * width + j - 1) * 4 + 1] + errorG * (3 / 16.0) < 0)
                    {
                        data[((i + 1) * width + j - 1) * 4 + 1] = 0;
                    }
                    else
                    {
                        data[((i + 1) * width + j - 1) * 4 + 1] += errorG * (3 / 16.0);
                    }
                    if (data[((i + 1) * width + j - 1) * 4+2] + errorB * (3 / 16.0) > 255)
                    {
                        data[((i + 1) * width + j - 1) * 4+2] = 255;
                    }
                    else if (data[((i + 1) * width + j - 1) * 4+2] + errorB * (3 / 16.0) < 0)
                    {
                        data[((i + 1) * width + j - 1) * 4+2] = 0;
                    }
                    else
                    {
                        data[((i + 1) * width + j - 1) * 4 + 2] += errorB * (3 / 16.0);
                    }
                }
                if (i + 1 < height)
				{
					if (data[((i + 1) * width + j) * 4] + errorR * (5 / 16.0) > 255)
					{
						data[((i + 1) * width + j) * 4] = 255;
					}
					else if (data[((i + 1) * width + j) * 4] + errorR * (5 / 16.0) < 0)
					{
						data[((i + 1) * width + j) * 4] = 0;
					}
                    else
                    {
                        data[((i + 1) * width + j) * 4] += errorR * (5 / 16.0);
                    }
					if (data[((i + 1) * width + j) * 4 + 1] + errorG * (5 / 16.0) > 255)
					{
						data[((i + 1) * width + j) * 4 + 1] = 255;
					}
					else if (data[((i + 1) * width + j) * 4 + 1] + errorG * (5 / 16.0) < 0)
					{
						data[((i + 1) * width + j) * 4 + 1] = 0;
					}
                    else
                    {
                        data[((i + 1) * width + j) * 4 + 1] += errorG * (5 / 16.0);
                    }
					if (data[((i + 1) * width + j) * 4 + 2] +errorB * (5 / 16.0) > 255)
					{
						data[((i + 1) * width + j) * 4 + 2] = 255;
					}
					else if (data[((i + 1) * width + j) * 4 + 2] + errorB * (5 / 16.0) < 0)
					{
						data[((i + 1) * width + j) * 4 + 2] = 0;
					}
                    else
                    {
                        data[((i + 1) * width + j) * 4 + 2] += errorB * (5 / 16.0);
                    }
				}
                if (j + 1 < width && i + 1 < height)
                {
                    if (data[((i + 1) * width + j + 1) * 4] + errorR * (1 / 16.0) > 255)
                    {
                        data[((i + 1) * width + j + 1) * 4] = 255;
                    }
                    else if (data[((i + 1) * width + j + 1) * 4] +errorR * (1 / 16.0) < 0)
                    {
                        data[((i + 1) * width + j + 1) * 4] = 0;
                    }
                    else
                    {
                        data[((i + 1) * width + j + 1) * 4] += errorR * (1 / 16.0);
                    }
                    if (data[((i + 1) * width + j + 1) * 4+1] +errorG * (1 / 16.0) > 255)
                    {
                        data[((i + 1) * width + j + 1) * 4+1] = 255;
                    }
                    else if (data[((i + 1) * width + j + 1) * 4 + 1] + errorG * (1 / 16.0) < 0)
                    {
                        data[((i + 1) * width + j + 1) * 4 + 1] = 0;
                    }
                    else
                    {
                        data[((i + 1) * width + j + 1) * 4 + 1] += errorG * (1 / 16.0);
                    }
                    if (data[((i + 1) * width + j + 1) * 4 + 2] +errorB * (1 / 16.0) > 255)
                    {
                        data[((i + 1) * width + j + 1) * 4 + 2] = 255;
                    }
                    else if (data[((i + 1) * width + j + 1) * 4+2] +errorB * (1 / 16.0) < 0)
                    {
                        data[((i + 1) * width + j + 1) * 4+2] = 0;
                    }
                    else
                    {
                        data[((i + 1) * width + j + 1) * 4 + 2] += errorB * (1 / 16.0);
                    }
                }
            }
        }
        else
        {
            for (int j = width - 1; j >= 0; j--)
            {
                double errorR, errorG, errorB;
                double newR = data[(i * width + j) * 4] / 32.0;
                newR = (int)newR;
                newR *= 32;
                double newG = data[(i * width + j) * 4 + 1] / 32.0;
                newG = (int)newG;
                newG *= 32;
                double newB = data[(i * width + j) * 4 + 2] / 64.0;
                newB = (int)newB;
                newB *= 64;
                errorR = data[(i * width + j) * 4] - newR;
                errorG = data[(i * width + j) * 4 + 1] - newG;
                errorB = data[(i * width + j) * 4 + 2] - newB;

                data[(i * width + j) * 4] = newR;
                data[(i * width + j) * 4 + 1] = newG;
                data[(i * width + j) * 4 + 2] = newB;

                if (j - 1 >= 0)
                {
                    if (data[(i * width + j - 1) * 4] + errorR * (7 / 16.0) > 255)
                    {
                        data[(i * width + j - 1) * 4] = 255;
                    }
                    else if (data[(i * width + j - 1) * 4] + errorR * (7 / 16.0) < 0)
                    {
                        data[(i * width + j - 1) * 4] = 0;
                    }
                    else
                    {
                        data[(i * width + j - 1) * 4] += errorR * (7 / 16.0);
                    }
                    if (data[(i * width + j - 1) * 4+1] +errorG * (7 / 16.0) > 255)
                    {
                        data[(i * width + j - 1) * 4+1] = 255;
                    }
                    else if (data[(i * width + j - 1) * 4 + 1] + errorG * (7 / 16.0) < 0)
                    {
                        data[(i * width + j - 1) * 4 + 1] = 0;
                    }
                    else
                    {
                        data[(i * width + j - 1) * 4 + 1] += errorG * (7 / 16.0);
                    }
                    if (data[(i * width + j - 1) * 4+2] + errorB * (7 / 16.0) > 255)
                    {
                        data[(i * width + j - 1) * 4+2] = 255;
                    }
                    else if (data[(i * width + j - 1) * 4+2] + errorB * (7 / 16.0) < 0)
                    {
                        data[(i * width + j - 1) * 4+2] = 0;
                    }
                    else
                    {
                        data[(i * width + j - 1) * 4 + 2] += errorB * (7 / 16.0);
                    }
                }
                if (j + 1 < width && i + 1 < height)
                {
                    if (data[((i + 1) * width + j + 1) * 4] + errorR * (3 / 16.0) > 255)
                    {
                        data[((i + 1) * width + j + 1) * 4] = 255;
                    }
                    else if (data[((i + 1) * width + j + 1) * 4] + errorR * (3 / 16.0) < 0)
                    {
                        data[((i + 1) * width + j + 1) * 4] = 0;
                    }
                    else
                    {
                        data[((i + 1) * width + j + 1) * 4] += errorR * (3 / 16.0);
                    }
                    if (data[((i + 1) * width + j + 1) * 4+1] + errorG * (3 / 16.0) > 255)
                    {
                        data[((i + 1) * width + j + 1) * 4+1] = 255;
                    }
                    else if (data[((i + 1) * width + j + 1) * 4 + 1] + errorG * (3 / 16.0) < 0)
                    {
                        data[((i + 1) * width + j + 1) * 4 + 1] = 0;
                    }
                    else
                    {
                        data[((i + 1) * width + j + 1) * 4 + 1] += errorG * (3 / 16.0);
                    }
                    if (data[((i + 1) * width + j + 1) * 4+2] + errorB * (3 / 16.0) > 255)
                    {
                        data[((i + 1) * width + j + 1) * 4l+2] = 255;
                    }
                    else if (data[((i + 1) * width + j + 1) * 4+2] + errorB * (3 / 16.0) < 0)
                    {
                        data[((i + 1) * width + j + 1) * 4+2] = 0;
                    }
                    else
                    {
                        data[((i + 1) * width + j + 1) * 4 + 2] += errorB * (3 / 16.0);
                    }
                }
                if (i + 1 < height)
                {
                    data[((i + 1) * width + j) * 4] ;
                    data[((i + 1) * width + j) * 4+1] ;
                    data[((i + 1) * width + j) * 4+2] ;
                    if (data[((i + 1) * width + j) * 4] + errorR * (5 / 16.0) > 255)
                    {
                        data[((i + 1) * width + j) * 4] = 255;
                    }
                    else if (data[((i + 1) * width + j) * 4] + errorR * (5 / 16.0) < 0)
                    {
                        data[((i + 1) * width + j) * 4] = 0;
                    }
                    else
                    {
                        data[((i + 1) * width + j) * 4] += errorR * (5 / 16.0);
                    }
                    if (data[((i + 1) * width + j) * 4+1] + errorG * (5 / 16.0) > 255)
                    {
                        data[((i + 1) * width + j) * 4+1] = 255;
                    }
                    else if (data[((i + 1) * width + j) * 4 + 1] + errorG * (5 / 16.0) < 0)
                    {
                        data[((i + 1) * width + j) * 4 + 1] = 0;
                    }
                    else
                    {
                        data[((i + 1) * width + j) * 4 + 1] += errorG * (5 / 16.0);
                    }
                    if (data[((i + 1) * width + j) * 4+2] + errorB * (5 / 16.0) > 255)
                    {
                        data[((i + 1) * width + j) * 4+2] = 255;
                    }
                    else if (data[((i + 1) * width + j) * 4+2] + errorB * (5 / 16.0) < 0)
                    {
                        data[((i + 1) * width + j) * 4+2] = 0;
                    }
                    else
                    {
                        data[((i + 1) * width + j) * 4 + 2] += errorB * (5 / 16.0);
                    }
                }
                if (j - 1 >= 0 && i + 1 < height)
                {
                    data[((i + 1) * width + j - 1) * 4];
                    data[((i + 1) * width + j - 1) * 4+1];
                    data[((i + 1) * width + j - 1) * 4+2] ;
                    if (data[((i + 1) * width + j - 1) * 4] + errorR * (1 / 16.0) > 255)
                    {
                        data[((i + 1) * width + j - 1) * 4] = 255;
                    }
                    else if (data[((i + 1) * width + j - 1) * 4] + errorR * (1 / 16.0) < 0)
                    {
                        data[((i + 1) * width + j - 1) * 4] = 0;
                    }
                    else
                    {
                        data[((i + 1) * width + j - 1) * 4] += errorR * (1 / 16.0);
                    }
                   if (data[((i + 1) * width + j - 1) * 4+1] + errorG * (1 / 16.0) > 255)
                    {
                        data[((i + 1) * width + j - 1) * 4+1] = 255;
                    }
                   else if (data[((i + 1) * width + j - 1) * 4 + 1] + errorG * (1 / 16.0) < 0)
                   {
                       data[((i + 1) * width + j - 1) * 4 + 1] = 0;
                   }
                   else
                   {
                       data[((i + 1) * width + j - 1) * 4 + 1] += errorG * (1 / 16.0);
                   }
                    if (data[((i + 1) * width + j - 1) * 4+2] + errorB * (1 / 16.0) > 255)
                    {
                        data[((i + 1) * width + j - 1) * 4+2] = 255;
                    }
                    else if (data[((i + 1) * width + j - 1) * 4+2] + errorB * (1 / 16.0) < 0)
                    {
                        data[((i + 1) * width + j - 1) * 4+2] = 0;
                    }
                    else
                    {
                        data[((i + 1) * width + j - 1) * 4 + 2] += errorB * (1 / 16.0);
                    }
                }
            }
        }
    }
    return true;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout <<  "Comp_Over: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_In: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Out: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Atop: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Xor: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Xor


///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this imag and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
    if (!pImage)
        return false;

    if (width != pImage->width || height != pImage->height)
    {
        cout << "Difference: Images not the same size\n";
        return false;
    }// if

    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i+1] = abs(rgb1[1] - rgb2[1]);
        data[i+2] = abs(rgb1[2] - rgb2[2]);
        data[i+3] = 255;
    }

    return true;
}// Difference


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box() //method 3
{
    unsigned char* newArr = new unsigned char[width * height*4];
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            newArr[(i * width + j) * 4] = data[(i * width + j) * 4];
            newArr[(i * width + j) * 4+1] = data[(i * width + j) * 4+1];
            newArr[(i * width + j) * 4+2] = data[(i * width + j) * 4+2];
            newArr[(i * width + j) * 4 + 3] = data[(i * width + j) * 4 + 3];
        }
    }
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            double sum1 = 0;
            double sum2 = 0;
            double sum3 = 0;
            for (int m = -2; m < 3; m++)
            {
                for (int n = -2; n < 3; n++)
                {
                    int p=i+m, r=j+n;
                    if (p < 0 || p>=height)
                    {
                        p = i - m;
                    }
                    if (r < 0 || r >= width)
                    {
                        r = j - n;
                    }
                        sum1 += newArr[(p * width + r) * 4];
                        sum2 += newArr[(p * width + r) * 4 + 1];
                        sum3 += newArr[(p * width + r) * 4 + 2];
                }
            }
            data[(i * width + j) * 4] = sum1 / 25;
            data[(i * width + j) * 4+1] = sum2 / 25;
            data[(i * width + j) * 4+2] = sum3 / 25;
        }
    }
    return true;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
    unsigned char* newArr = new unsigned char[width * height * 4];
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            newArr[(i * width + j) * 4] = data[(i * width + j) * 4];
            newArr[(i * width + j) * 4 + 1] = data[(i * width + j) * 4 + 1];
            newArr[(i * width + j) * 4 + 2] = data[(i * width + j) * 4 + 2];
            newArr[(i * width + j) * 4 + 3] = data[(i * width + j) * 4 + 3];
        }
    }
    double table[5][5] = { 0.0 };
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            table[i][j] = table[i][4 - j] = table[4 - i][4 - j] = table[4 - i][j] = ((j + 1) * (i + 1))/81.0;
        }
    }
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            double sum1 = 0;
            double sum2 = 0;
            double sum3 = 0;
            for (int m = -2; m < 3; m++)
            {
                for (int n = -2; n < 3; n++)
                {
                    int p = i + m, r = j + n;
                    if (p < 0 || p >= height)
                    {
                        p = i - m;
                    }
                    if (r < 0 || r >= width)
                    {
                        r = j - n;
                    }
                    sum1 += newArr[(p * width + r) * 4]*table[m+2][n+2];
                    sum2 += newArr[(p * width + r) * 4 + 1]*table[m+2][n+2];
                    sum3 += newArr[(p * width + r) * 4 + 2]*table[m+2][n+2];
                }
            }
            data[(i * width + j) * 4] = sum1 ;
            data[(i * width + j) * 4 + 1] = sum2;
            data[(i * width + j) * 4 + 2] = sum3;
        }
    }
    return true;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
    int vec1[5] = { 1,4,6,4,1 };
    double divide = 256.0;
    double** table;
    table = new double* [5];
    for (int i = 0; i < 5; i++)
    {
        table[i] = new double[5];
    }
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            table[i][j] = vec1[i] * vec1[j] / divide;
        }
    }
    unsigned char* newArr = new unsigned char[width * height * 4];
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            newArr[(i * width + j) * 4] = data[(i * width + j) * 4];
            newArr[(i * width + j) * 4 + 1] = data[(i * width + j) * 4 + 1];
            newArr[(i * width + j) * 4 + 2] = data[(i * width + j) * 4 + 2];
            newArr[(i * width + j) * 4 + 3] = data[(i * width + j) * 4 + 3];
        }
    }
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            double sum1 = 0;
            double sum2 = 0;
            double sum3 = 0;
            for (int m = (int)(5 / 2) * -1; m <= (int)(5 / 2); m++)
            {
                for (int n = (int)(5 / 2) * -1; n <= (int)(5 / 2); n++)
                {
                    int p = i + m, r = j + n;
                    if (p < 0 || p >= height)
                    {
                        p = i - m;
                    }
                    if (r < 0 || r >= width)
                    {
                        r = j - n;
                    }
                    sum1 += newArr[(p * width + r) * 4] * table[m + (int)(5 / 2)][n + (int)(5 / 2)];
                    sum2 += newArr[(p * width + r) * 4 + 1] * table[m + (int)(5 / 2)][n + (int)(5 / 2)];
                    sum3 += newArr[(p * width + r) * 4 + 2] * table[m + (int)(5 / 2)][n + (int)(5 / 2)];
                }
            }
            data[(i * width + j) * 4] = sum1;
            data[(i * width + j) * 4 + 1] = sum2;
            data[(i * width + j) * 4 + 2] = sum3;
        }
    }
    return true;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N( unsigned int N )
{
    vector<int> vec1;
    vector<int> vec2;
    for (int i = 1; i <= N-1; i++)
    {
        vec2.push_back(1);
        for (int j = 1; j <i; j++)
        {
            vec2.push_back(vec1[j - 1] + vec1[j]);
        }
        vec2.push_back(1);
        vec1.clear();
        vec1 = vec2;
        vec2.clear();
    }
    double divide = 0.0;
    double** table;
    table = new double*[N];
    for (int i = 0; i < N; i++)
    {
        divide += vec1[i];
        table[i] = new double[N];
    }
    divide = divide * divide;
    for (int i = 0; i < N; i++)
    {  
        for (int j = 0; j < N; j++)
        {
            table[i][j] = vec1[i] * vec1[j]/divide;
        }
    }
    unsigned char* newArr = new unsigned char[width * height * 4];
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            newArr[(i * width + j) * 4] = data[(i * width + j) * 4];
            newArr[(i * width + j) * 4 + 1] = data[(i * width + j) * 4 + 1];
            newArr[(i * width + j) * 4 + 2] = data[(i * width + j) * 4 + 2];
            newArr[(i * width + j) * 4 + 3] = data[(i * width + j) * 4 + 3];
        }
    }
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            double sum1 = 0;
            double sum2 = 0;
            double sum3 = 0;
            for (int m = (int)(N/2)*-1;m<=(int)(N / 2); m++)
            {
                for (int n = (int)(N / 2) * -1; n <= (int)(N / 2); n++)
                {
                    int p = i + m, r = j + n;
                    if (p < 0 || p >= height)
                    {
                        p = i - m;
                    }
                    if (r < 0 || r >= width)
                    {
                        r = j - n;
                    }
                    sum1 += newArr[(p * width + r) * 4] * table[m + (int)(N/2)][n + (int)(N / 2)];
                    sum2 += newArr[(p * width + r) * 4 + 1] * table[m + (int)(N / 2)][n + (int)(N / 2)];
                    sum3 += newArr[(p * width + r) * 4 + 2] * table[m + (int)(N / 2)][n + (int)(N / 2)];
                }
            }
            data[(i * width + j) * 4] = sum1;
            data[(i * width + j) * 4 + 1] = sum2;
            data[(i * width + j) * 4 + 2] = sum3;
        }
    }
   return true;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
    unsigned char* newArr = new unsigned char[width * height * 4];
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            newArr[(i * width + j) * 4] = data[(i * width + j) * 4];
            newArr[(i * width + j) * 4 + 1] = data[(i * width + j) * 4 + 1];
            newArr[(i * width + j) * 4 + 2] = data[(i * width + j) * 4 + 2];
            newArr[(i * width + j) * 4 + 3] = data[(i * width + j) * 4 + 3];
        }
    }
    double table[5][5] = { 0.0 };
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            table[i][j] = table[i][4 - j] = table[4 - i][4 - j] = table[4 - i][j] = ((j + 1) * (i + 1)) / 81.0;
        }
    }
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            if (i == 2 && j == 2)
            {
                table[i][j] = 1 - table[i][j];
            }
            else
            {
                table[i][j] *= -1.0;
            }
        }
    }
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            double sum1 = 0;
            double sum2 = 0;
            double sum3 = 0;
            for (int m = -2; m < 3; m++)
            {
                for (int n = -2; n < 3; n++)
                {
                    int p = i + m, r = j + n;
                    if (p < 0 || p >= height)
                    {
                        p = i - m;
                    }
                    if (r < 0 || r >= width)
                    {
                        r = j - n;
                    }
                    sum1 += newArr[(p * width + r) * 4] * table[m + 2][n + 2];
                    sum2 += newArr[(p * width + r) * 4 + 1] * table[m + 2][n + 2];
                    sum3 += newArr[(p * width + r) * 4 + 2] * table[m + 2][n + 2];
                }
            }
            if (sum1 > 255)
            {
                sum1 = 255;
            }
            else if (sum1 < 0)
            {
                sum1 = 0;
            }
            if (sum2 > 255)
            {
                sum2 = 255;
            }
            else if (sum2 < 0)
            {
                sum2 = 0;
            }
            if (sum3 > 255)
            {
                sum3 = 255;
            }
            else if (sum3 < 0)
            {
                sum3 = 0;
            }
            data[(i * width + j) * 4] = sum1;
            data[(i * width + j) * 4 + 1] = sum2;
            data[(i * width + j) * 4 + 2] = sum3;
        }
    }
    return true;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
    unsigned char* newArr = new unsigned char[width * height * 4];
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            newArr[(i * width + j) * 4] = data[(i * width + j) * 4];
            newArr[(i * width + j) * 4 + 1] = data[(i * width + j) * 4 + 1];
            newArr[(i * width + j) * 4 + 2] = data[(i * width + j) * 4 + 2];
            newArr[(i * width + j) * 4 + 3] = data[(i * width + j) * 4 + 3];
        }
    }
    double table[5][5] = { 0.0 };
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            table[i][j] = table[i][4 - j] = table[4 - i][4 - j] = table[4 - i][j] = ((j + 1) * (i + 1)) / 81.0;
        }
    }
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            if (i == 2 && j == 2)
            {
                table[i][j] = 2 - table[i][j];
            }
            else
            {
                table[i][j] *= -1.0;
            }
        }
    }
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            double sum1 = 0;
            double sum2 = 0;
            double sum3 = 0;
            for (int m = -2; m < 3; m++)
            {
                for (int n = -2; n < 3; n++)
                {
                    int p = i + m, r = j + n;
                    if (p < 0 || p >= height)
                    {
                        p = i - m;
                    }
                    if (r < 0 || r >= width)
                    {
                        r = j - n;
                    }
                    sum1 += newArr[(p * width + r) * 4] * table[m + 2][n + 2];
                    sum2 += newArr[(p * width + r) * 4 + 1] * table[m + 2][n + 2];
                    sum3 += newArr[(p * width + r) * 4 + 2] * table[m + 2][n + 2];
                }
            }
            if (sum1 > 255)
            {
                sum1 = 255;
            }
            else if (sum1 < 0)
            {
                sum1 = 0;
            }
            if (sum2 > 255)
            {
                sum2 = 255;
            }
            else if (sum2 < 0)
            {
                sum2 = 0;
            }
            if (sum3 > 255)
            {
                sum3 = 255;
            }
            else if (sum3 < 0)
            {
                sum3 = 0;
            }
            data[(i * width + j) * 4] = sum1;
            data[(i * width + j) * 4 + 1] = sum2;
            data[(i * width + j) * 4 + 2] = sum3;
        }
    }
    return true;
}// Filter_Enhance


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{
    ClearToBlack();
    return false;
}



///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
    int w2 = width / 2;
    int h2 = height / 2;
    unsigned char *newArr = new unsigned char[w2 *h2*4];
 
    double table[3][3] = { 1/16.0,1/8.0,1/16.0,1/8.0,1/4.0,1/8.0,1/16.0,1/8.0,1/16.0 };
    //for (int i = 0; i < 1; i++)
    //{
    //    for (int j = 0; j < 1; j++)
    //    {
    //        table[i][j] = table[i][2 - j] = table[2 - i][2 - j] = table[2 - i][j] = ((j + 1) * (i + 1)) / 16.0;
    //    }
    //}
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            double sum1 = 0;
            double sum2 = 0;
            double sum3 = 0;
            for (int m = -1; m < 2; m++)
            {
                for (int n = -1; n < 2; n++)
                {
                    int p = i + m, r = j + n;
                    if (p < 0 || p >= height)
                    {
                        p = i - m;
                    }
                    if (r < 0 || r >= width)
                    {
                        r = j - n;
                    }
                    sum1 += data[(p * width + r) * 4] * table[m + 1][n + 1];
                    sum2 += data[(p * width + r) * 4 + 1] * table[m + 1][n + 1];
                    sum3 += data[(p * width + r) * 4 + 2] * table[m + 1][n + 1];
                }
            }
            if (sum1 > 255)
            {
                sum1 = 255;
            }
            else if (sum1 < 0)
            {
                sum1 = 0;
            }
            if (sum2 > 255)
            {
                sum2 = 255;
            }
            else if (sum2 < 0)
            {
                sum2 = 0;
            }
            if (sum3 > 255)
            {
                sum3 = 255;
            }
            else if (sum3 < 0)
            {
                sum3 = 0;
            }
            newArr[((i/2) * w2 + j/2) * 4] = sum1;
            newArr[((i/2) * w2 + j/2) * 4 + 1] = sum2;
            newArr[((i/2) * w2 + j/2) * 4 + 2] = sum3;
            newArr[((i/2) * w2 + j/2) * 4 + 3] = data[(i * width + j) * 4 + 3];
        }
    }
    data = newArr;
    width /= 2;
    height /= 2;
    return true;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
    ClearToBlack();
    return false;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
    ClearToBlack();
    return false;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
    ClearToBlack();
    return false;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
    const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

    unsigned char  alpha = rgba[3];

    if (alpha == 0)
    {
        rgb[0] = BACKGROUND[0];
        rgb[1] = BACKGROUND[1];
        rgb[2] = BACKGROUND[2];
    }
    else
    {
	    float	alpha_scale = (float)255 / (float)alpha;
	    int	val;
	    int	i;

	    for (i = 0 ; i < 3 ; i++)
	    {
	        val = (int)floor(rgba[i] * alpha_scale);
	        if (val < 0)
		    rgb[i] = 0;
	        else if (val > 255)
		    rgb[i] = 255;
	        else
		    rgb[i] = val;
	    }
    }
}// RGA_To_RGB


///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char   *dest = new unsigned char[width * height * 4];
    TargaImage	    *result;
    int 	        i, j;

    if (! data)
    	return NULL;

    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = (height - i - 1) * width * 4;
	    int out_offset = i * width * 4;

	    for (j = 0 ; j < width ; j++)
        {
	        dest[out_offset + j * 4] = data[in_offset + j * 4];
	        dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
	        dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
	        dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
        }
    }

    result = new TargaImage(width, height, dest);
    delete[] dest;
    return result;
}// Reverse_Rows


///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
    memset(data, 0, width * height * 4);
}// ClearToBlack


///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
   int radius_squared = (int)s.radius * (int)s.radius;
   for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
      for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
         int x_loc = (int)s.x + x_off;
         int y_loc = (int)s.y + y_off;
         // are we inside the circle, and inside the image?
         if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
            int dist_squared = x_off * x_off + y_off * y_off;
            if (dist_squared <= radius_squared) {
               data[(y_loc * width + x_loc) * 4 + 0] = s.r;
               data[(y_loc * width + x_loc) * 4 + 1] = s.g;
               data[(y_loc * width + x_loc) * 4 + 2] = s.b;
               data[(y_loc * width + x_loc) * 4 + 3] = s.a;
            } else if (dist_squared == radius_squared + 1) {
               data[(y_loc * width + x_loc) * 4 + 0] = 
                  (data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
               data[(y_loc * width + x_loc) * 4 + 1] = 
                  (data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
               data[(y_loc * width + x_loc) * 4 + 2] = 
                  (data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
               data[(y_loc * width + x_loc) * 4 + 3] = 
                  (data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
            }
         }
      }
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
               unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
   radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}

