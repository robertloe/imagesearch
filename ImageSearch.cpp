// build:
//
// rc imagesearch.rc
// cl -EHsc -DWIN32 -D_WINDOWS -O2 imagesearch.cpp imagesearch.res

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tchar.h>
#include <math.h>


#include <iostream>
using namespace std;
#include <string>
#include <map>
#include <algorithm>
#include <vector>


#pragma comment(lib,"ws2_32.lib")
#pragma comment(lib,"kernel32.lib")
#pragma comment(lib,"user32.lib")
#pragma comment(lib,"gdi32.lib")
#pragma comment(lib,"winspool.lib")
#pragma comment(lib,"comdlg32.lib")
#pragma comment(lib,"advapi32.lib")
#pragma comment(lib,"shell32.lib")
#pragma comment(lib,"ole32.lib")
#pragma comment(lib,"oleaut32.lib")
#pragma comment(lib,"uuid.lib")
#pragma comment(lib,"odbc32.lib")
#pragma comment(lib,"odbccp32.lib")
#pragma comment(lib, "user32.lib")


#define IDD_DIALOG5 100
#define IDC_EDIT1 1000

char m_parameter[256];


int DialogType = 0;

LRESULT CALLBACK SetParameters(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
    UINT wmId,wmEvent;

    switch (message)
    {

	case WM_INITDIALOG:
		if (DialogType == 0)
			SetWindowText(hDlg,"Enter horizontal scale 0.0 to 1.0");
		if (DialogType == 1)
			SetWindowText(hDlg,"Enter threshold 0 to 255");
		SetDlgItemText(hDlg,IDC_EDIT1,m_parameter);
		ShowWindow(hDlg,SW_SHOWNORMAL);
		SetFocus(GetDlgItem(hDlg,IDC_EDIT1));
		return (FALSE);

	case WM_COMMAND:

		wmId = LOWORD(wParam);
		wmEvent = HIWORD(wParam);
		switch (wmId)
		{
		case IDOK:
			GetDlgItemText(hDlg,IDC_EDIT1,m_parameter,sizeof(m_parameter)-1);
			EndDialog(hDlg,TRUE);
			return(TRUE);

		case IDCANCEL:
            EndDialog(hDlg,FALSE);
			return(TRUE);
	    }
		break;
    }
	return (FALSE);
}



////////

#define TRUE 1
#define FALSE 0

int XSIZE,YSIZE;
int STORE;

typedef struct
{
	double real;
	double imag;
} COMPLEX;

COMPLEX *data,*filter,*result,*change,*pm,*savefilter;
double *power;
char *mark;

typedef pair<double,string> ListPair;
typedef map<double,string> List;

List save;
List best;


////////


#define ID_FILE_OPEN	1001
#define ID_FILE_OPEN2	1002
#define ID_FILE_PATTERN	1003
#define ID_FILE_PROCESS	1004
#define ID_FILE_ENERGYCENTROID	1005
#define ID_EXIT	1006


// Global variables

// The main window class name.
static TCHAR szWindowClass[] = _T("win32app");

// The string that appears in the application's title bar.
static TCHAR szTitle[] = _T("Simple");


char m_filterfile[256],m_bmpfile[256];

int xClient,yClient;
double xAdjust = 1.0;

OPENFILENAME ofn;
char szDirName[256];
char szFile[256]={"\0"},szFileTitle[256]={"\0"};
char szInputFile[256],szOutputFile[256];
char szReadFilter[81]=
{"BMP files (*.BMP)\0*.BMP\0All files (*.*)\0*.*\0\0"};
char szFilterFilter[81]=
{"FLT files (*.FLT)\0*.FLT\0All files (*.*)\0*.*\0\0"};

HBITMAP m_bitmap = NULL;
int m_width,m_height;

HBITMAP m_processed = NULL;

HGLOBAL hMem = NULL;

HINSTANCE hInst;

// Forward declarations of functions included in this code module:
LRESULT CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM);
void CreateAMenu(HWND hWnd);


//////////


/*
	This function writes out a 24-bit Windows bitmap file that is readable by Microsoft Paint.
   The image data	is a 1D array of (r, g, b) triples, where individual (r, g, b) values can
	each take on values between 0 and 255, inclusive.

   The input to the function is:
	   char *filename:						A string representing the filename that will be written
		unsigned int width:					The width, in pixels, of the bitmap
		unsigned int height:					The height, in pixels, of the bitmap
		unsigned char *image:				The image data, where each pixel is 3 unsigned chars (r, g, b)


   Written by Greg Slabaugh (greg_slabaugh@hotmail.com), 10/19/00
*/
int write24BitBmpFile(char *filename, unsigned int width, unsigned int height, unsigned char *image)
{
	BITMAPINFOHEADER bmpInfoHeader;
	BITMAPFILEHEADER bmpFileHeader;
	FILE *filep;
	unsigned int row, column;
	unsigned int extrabytes, bytesize;
	unsigned char *paddedImage = NULL, *paddedImagePtr, *imagePtr;


	/* The .bmp format requires that the image data is aligned on a 4 byte boundary.  For 24 - bit bitmaps,
	   this means that the width of the bitmap  * 3 must be a multiple of 4. This code determines
	   the extra padding needed to meet this requirement. */
   extrabytes = (4 - (width * 3) % 4) % 4;


	// This is the size of the padded bitmap
	bytesize = (width * 3 + extrabytes) * height;


	// Fill the bitmap file header structure
	bmpFileHeader.bfType = 'MB';   // Bitmap header
	bmpFileHeader.bfSize = 0;      // This can be 0 for BI_RGB bitmaps
	bmpFileHeader.bfReserved1 = 0;
	bmpFileHeader.bfReserved2 = 0;
	bmpFileHeader.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);


	// Fill the bitmap info structure
	bmpInfoHeader.biSize = sizeof(BITMAPINFOHEADER);
	bmpInfoHeader.biWidth = width;
	bmpInfoHeader.biHeight = height;
	bmpInfoHeader.biPlanes = 1;
	bmpInfoHeader.biBitCount = 24;            // 8 - bit bitmap
	bmpInfoHeader.biCompression = BI_RGB;
	bmpInfoHeader.biSizeImage = bytesize;     // includes padding for 4 byte alignment
	bmpInfoHeader.biXPelsPerMeter = 0;
	bmpInfoHeader.biYPelsPerMeter = 0;
	bmpInfoHeader.biClrUsed = 0;
	bmpInfoHeader.biClrImportant = 0;




	// Open file
	if ((filep = fopen(filename, "wb")) == NULL) {
		printf("Error opening file %s\n", filename);
		return FALSE;
	}


	// Write bmp file header
	if (fwrite(&bmpFileHeader, 1, sizeof(BITMAPFILEHEADER), filep) < sizeof(BITMAPFILEHEADER)) {
		printf("Error writing bitmap file header\n");
		fclose(filep);
		return FALSE;
	}


	// Write bmp info header
	if (fwrite(&bmpInfoHeader, 1, sizeof(BITMAPINFOHEADER), filep) < sizeof(BITMAPINFOHEADER)) {
		printf("Error writing bitmap info header\n");
		fclose(filep);
		return FALSE;
	}


	// Allocate memory for some temporary storage
	paddedImage = (unsigned char *)calloc(sizeof(unsigned char), bytesize);
	if (paddedImage == NULL) {
		printf("Error allocating memory \n");
		fclose(filep);
		return FALSE;
	}


	/* This code does three things.  First, it flips the image data upside down, as the .bmp
	   format requires an upside down image.  Second, it pads the image data with extrabytes
		number of bytes so that the width in bytes of the image data that is written to the
		file is a multiple of 4.  Finally, it swaps (r, g, b) for (b, g, r).  This is another
		quirk of the .bmp file format. */

	for (row = 0; row < height; row++) {
//LOE		imagePtr = image + (height - 1 - row) * width * 3;
		imagePtr = image + (height - 1 - row) * width * 4;
		paddedImagePtr = paddedImage + row * (width * 3 + extrabytes);
		for (column = 0; column < width; column++) {
//LOE			*paddedImagePtr = *(imagePtr + 2);
//LOE			*(paddedImagePtr + 1) = *(imagePtr + 1);
//LOE			*(paddedImagePtr + 2) = *imagePtr;
			*paddedImagePtr = *(imagePtr + 0);
			*(paddedImagePtr + 1) = *(imagePtr + 1);
			*(paddedImagePtr + 2) = *(imagePtr + 2);
//			imagePtr += 3;
			imagePtr += 4;
			paddedImagePtr += 3;
		}
	}


	// Write bmp data
	if (fwrite(paddedImage, 1, bytesize, filep) < bytesize) {
		printf("Error writing bitmap data\n");
		free(paddedImage);
		fclose(filep);
		return FALSE;
	}


	// Close file
	fclose(filep);
	free(paddedImage);
	return TRUE;
}


int Powerof2(int n,int *m,int *twopm)
{
	if (n <= 1) {
		*m = 0;
		*twopm = 1;
		return(FALSE);
	}

   *m = 1;
   *twopm = 2;
   do {
      (*m)++;
      (*twopm) *= 2;
   } while (2*(*twopm) <= n);

   if (*twopm != n)
		return(FALSE);
	else
		return(TRUE);
}



//	http://local.wasp.uwa.edu.au/~pbourke/miscellaneous/dft/
//	http://paulbourke.net/miscellaneous/dft/

/*-------------------------------------------------------------------------
   This computes an in-place complex-to-complex FFT
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform

     Formula: forward
                  N-1
                  ---
              1   \          - j k 2 pi n / N
      X(n) = ---   >   x(k) e                    = forward transform
              N   /                                n=0..N-1
                  ---
                  k=0

      Formula: reverse
                  N-1
                  ---
                  \          j k 2 pi n / N
      X(n) =       >   x(k) e                    = forward transform
                  /                                n=0..N-1
                  ---
                  k=0
*/
int FFT(int dir,int m,double *x,double *y)
{
   long nn,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   nn = 1;
   for (i=0;i<m;i++)
      nn *= 2;

   /* Do the bit reversal */
   i2 = nn >> 1;
   j = 0;
   for (i=0;i<nn-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0;
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0;
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<nn;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1;
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1)
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<nn;i++) {
         x[i] /= (double)nn;
         y[i] /= (double)nn;
      }
   }

   return(TRUE);
}



//	http://local.wasp.uwa.edu.au/~pbourke/miscellaneous/dft/
//	http://paulbourke.net/miscellaneous/dft/

/*-------------------------------------------------------------------------
   Perform a 2D FFT inplace given a complex 2D array
   The direction dir, 1 for forward, -1 for reverse
   The size of the array (nx,ny)
   Return false if there are memory problems or
      the dimensions are not powers of 2
*/
int FFT2D(COMPLEX *c,int nx,int ny,int dir)
{
   int i,j;
   int m,twopm;
   double *real,*imag;

   /* Transform the rows */
   real = (double *)malloc(nx * sizeof(double));
   imag = (double *)malloc(nx * sizeof(double));
   if (real == NULL || imag == NULL)
      return(FALSE);
   if (!Powerof2(nx,&m,&twopm) || twopm != nx)
      return(FALSE);
   for (j=0;j<ny;j++) {
      for (i=0;i<nx;i++) {
         real[i] = c[i*ny+j].real;
         imag[i] = c[i*ny+j].imag;
      }
      FFT(dir,m,real,imag);
      for (i=0;i<nx;i++) {
         c[i*ny+j].real = real[i];
         c[i*ny+j].imag = imag[i];
      }
   }
   free(real);
   free(imag);

   /* Transform the columns */
   real = (double *)malloc(ny * sizeof(double));
   imag = (double *)malloc(ny * sizeof(double));
   if (real == NULL || imag == NULL)
      return(FALSE);
   if (!Powerof2(ny,&m,&twopm) || twopm != ny)
      return(FALSE);
   for (i=0;i<nx;i++) {
      for (j=0;j<ny;j++) {
         real[j] = c[i*ny+j].real;
         imag[j] = c[i*ny+j].imag;
      }
      FFT(dir,m,real,imag);
      for (j=0;j<ny;j++) {
         c[i*ny+j].real = real[j];
         c[i*ny+j].imag = imag[j];
      }
   }
   free(real);
   free(imag);

   return(TRUE);
}



void ProcessImage(char *infile, char *filterfile, string &results)
{
	int i,j,adjust,loc,choice,loc2;
	char outfile[512],buffer[256],aa[256],bb[256],cc[256];
	int r,g,b,rowval,colval;
	double val,maxval;
	double lowpower = 0.0;
	FILE *in;

	results = "";

	DWORD newarea = m_height * m_width * 4;

	BYTE *hpbits = new BYTE[newarea];
	printf("Image: Read %d byes; allocated %d bytes\n",GetBitmapBits(m_bitmap,newarea,hpbits),newarea);


//	calculate total size to avoid circular convolution

	int totalWidth,totalHeight;
	int imageWidth,imageHeight,patternWidth,patternHeight;
	int newwidth,newheight;

	imageWidth = m_width;
	imageHeight = m_height;
	patternWidth = 0;
	patternHeight = 0;

	totalWidth = imageWidth + patternWidth;
	totalHeight = imageHeight + patternHeight;

	adjust = (int)(log((double)totalWidth) / log(2.0));
	newwidth = (int)pow(2.0,(double)adjust);
	if (newwidth < totalWidth) adjust += 1;
	newwidth = (int)pow(2.0,(double)adjust);
	sprintf(buffer,"width %d %d %d %d\r\n",imageWidth,patternWidth,totalWidth,newwidth);
	results += buffer;

	adjust = (int)(log((double)totalHeight) / log(2.0));
	newheight = (int)pow(2.0,(double)adjust);
	if (newheight < totalHeight) adjust += 1;
	newheight = (int)pow(2.0,(double)adjust);
	sprintf(buffer,"height %d %d %d %d\r\n",imageHeight,patternWidth,totalHeight,newheight);
	results += buffer;


	XSIZE = newwidth;
	YSIZE = newheight;
	STORE = XSIZE * YSIZE;


	data = new COMPLEX[STORE];
	filter = new COMPLEX[STORE];
	result = new COMPLEX[STORE];
	change = new COMPLEX[STORE];
	power = new double[STORE];



	for (i=0; i<YSIZE; ++i)
	{
		for (j=0; j<XSIZE; ++j)
		{
			data[i*XSIZE+j].real = 0;
			data[i*XSIZE+j].imag = 0;

			filter[i*XSIZE+j].real = 0;
			filter[i*XSIZE+j].imag = 0;
		}
	}



//	convert rgb to grayscale


	loc = 0;
	for (i=0; i<m_height && i < YSIZE; ++i)
	{
		for (j=0; j<m_width && j < XSIZE; ++j)
		{
			r = hpbits[loc++];
			g = hpbits[loc++];
			b = hpbits[loc++];
			loc++;

//			data[i*XSIZE + j].real = (r + g + b) / 3 ;

            data[i*XSIZE + j].real = (((r * 30.0) / 100.0) + ((g * 59.0) / 100.0) + ((b * 11.0) / 100.0));
		}
	}

	delete [] hpbits ;




	in = fopen(filterfile,"r");
	if (in == NULL)
	{
		printf("Error opening %s for input.\n",filterfile);
		return ;
	}

	i = 0;
	fgets(buffer,255,in);
	while (feof(in) == 0 && i < 3)
	{
		sscanf(buffer,"%s %s %s",aa,bb,cc);
		filter[i*XSIZE+0].real =  atof(aa);
		filter[i*XSIZE+1].real =  atof(bb);
		filter[i*XSIZE+2].real =  atof(cc);
		++i;

		fgets(buffer,255,in);
	}
	fclose(in);

	sprintf(buffer,"%4.2f %4.2f %4.2f\r\n",filter[0*XSIZE+0].real,filter[0*XSIZE+1].real,filter[0*XSIZE+2].real);
	results += buffer;
	sprintf(buffer,"%4.2f %4.2f %4.2f\r\n",filter[1*XSIZE+0].real,filter[1*XSIZE+1].real,filter[1*XSIZE+2].real);
	results += buffer;
	sprintf(buffer,"%4.2f %4.2f %4.2f\r\n",filter[2*XSIZE+0].real,filter[2*XSIZE+1].real,filter[2*XSIZE+2].real);
	results += buffer;




//	perform correlation / filtering

	FFT2D(data,XSIZE,YSIZE,1);
	FFT2D(filter,XSIZE,YSIZE,1);

	for (i=0; i<YSIZE; ++i)
	{
		for (j=0; j<XSIZE; ++j)
		{
			loc = i*XSIZE+j;
			result[loc].real = data[loc].real * filter[loc].real - data[loc].imag * filter[loc].imag ;
			result[loc].imag = data[loc].imag * filter[loc].real + data[loc].real * filter[loc].imag ;
		}
	}

	FFT2D(result,XSIZE,YSIZE,-1);



	maxval = 0;
	rowval = 0;
	colval = 0;
	for (i=0; i<YSIZE; ++i)
	{
		for (j=0; j<XSIZE; ++j)
		{
	 		val = sqrt(result[i*XSIZE+j].real*result[i*XSIZE+j].real+result[i*XSIZE+j].imag*result[i*XSIZE+j].imag);
	 		power[i*XSIZE+j] = val;

	 		if (val > maxval)
	 		{
				maxval = val;
				rowval = i;
				colval = j;
			}
		}
	}



//	generate RESULTS.BMP that shows the best match location and center

	newarea = m_height * m_width * 4;

	hpbits = new BYTE[newarea];

	BYTE *lpMem;

	lpMem = (BYTE *)GlobalLock(hMem);


	loc = 0;
	loc2 = 0;
	for (i=0; i<m_height; ++i)
	{
		for (j=0; j<m_width; ++j)
		{
	 		val = power[i*XSIZE+j];
			r = (int)(val / maxval * 255.0);
			g = r;
			b = r;

			hpbits[loc++] = r;
			hpbits[loc++] = g;
			hpbits[loc++] = b;
			hpbits[loc++] = 0;

			lpMem[loc2++] = r;
			lpMem[loc2++] = g;
			lpMem[loc2++] = b;
			lpMem[loc2++] = 0;
		}
	}


	if (m_processed) DeleteObject(m_processed);

	m_processed = CreateBitmap(m_width, m_height, 1, 32, lpMem);
	GlobalUnlock(hMem);

	strcpy(outfile,infile);
	outfile[strcspn(outfile,".")] = 0;
	strcat(outfile,"_");

	i = strlen(filterfile) - 1;
	while (i > 0 && filterfile[i] != '\\') --i;
	++i;

	strcat(outfile,&filterfile[i]);
	outfile[strcspn(outfile,".")] = 0;
	strcat(outfile,"_results.bmp");
	write24BitBmpFile(outfile, m_width,m_height,hpbits);

	delete [] hpbits ;


	sprintf(buffer,"Output is %s\n",outfile);
	results += buffer;

	delete [] data;
	delete [] filter;
	delete [] result;
	delete [] change;
	delete [] power;

}




void PatternMatch(BYTE *imageBits, int imageWidth, int imageHeight,
	BYTE *patternBits, int patternWidth, int patternHeight,
	double azfov, double elfov,
	int addNoise, int levelTranslate,
	double matchCutoff, string &results)
{
	int FullAnalysis = 0;

	int i,j,adjust,loc,newwidth,newheight;

	int r,g,b,rowval,colval;
	double val,maxval,bestpossible;
	double lowpower = 0.0;

	HBITMAP m_bitmap;
    BITMAP bmImage,bmPattern;
	char patfile[256],outfile[256],buffer[256];

	List::iterator it;
	List::reverse_iterator rit;


	results = "";

//	calculate total size to avoid circular convolution

	int totalWidth,totalHeight;

	totalWidth = imageWidth + patternWidth;
	totalHeight = imageHeight + patternHeight;

	adjust = (int)(log((double)totalWidth) / log(2.0));
	newwidth = (int)pow(2.0,(double)adjust);
	if (newwidth < totalWidth) adjust += 1;
	newwidth = (int)pow(2.0,(double)adjust);
//	printf("width %d %d %d %d\n",imageWidth,patternWidth,totalWidth,newwidth);

	adjust = (int)(log((double)totalHeight) / log(2.0));
	newheight = (int)pow(2.0,(double)adjust);
	if (newheight < totalHeight) adjust += 1;
	newheight = (int)pow(2.0,(double)adjust);
//	printf("height %d %d %d %d\n",imageHeight,patternWidth,totalHeight,newheight);


	XSIZE = newwidth;
	YSIZE = newheight;
	STORE = XSIZE * YSIZE;

	data = new COMPLEX[STORE];
	filter = new COMPLEX[STORE];
	result = new COMPLEX[STORE];
	change = new COMPLEX[STORE];
	pm = new COMPLEX[STORE];
	savefilter = new COMPLEX[STORE];
	power = new double[STORE];
	mark = new char[STORE];


	for (i=0; i<YSIZE; ++i)
	{
		for (j=0; j<XSIZE; ++j)
		{
			data[i*XSIZE+j].real = 0;
			data[i*XSIZE+j].imag = 0;

			pm[i*XSIZE+j].real = 0;
			pm[i*XSIZE+j].imag = 0;

			filter[i*XSIZE+j].real = 0;
			filter[i*XSIZE+j].imag = 0;

			savefilter[i*XSIZE+j].real = 0;
			savefilter[i*XSIZE+j].imag = 0;
		}
	}





//	convert rgb to grayscale and store in data array

	double avgPixel = 0.0;
	int lowPixel = 255;
	int highPixel = 0;
	double totalBits = 0.0;
	int totalCount = 0;
	int k;

	loc = 0;
	for (i=0; i<imageHeight && i < YSIZE; ++i)
	{
		for (j=0; j<imageWidth && j < XSIZE; ++j)
		{
			r = imageBits[loc++];
			g = imageBits[loc++];
			b = imageBits[loc++];
			loc++;

			k = i*XSIZE + j;
			data[k].real = (r + g + b) / 3 ;

			if (data[k].real < lowPixel) lowPixel = data[k].real;
			if (data[k].real > highPixel) highPixel = data[k].real;

			totalBits = totalBits + ((double)(r + g + b) / 3.0);
			++totalCount;
		}
	}


	if (totalCount > 0)
		avgPixel = totalBits / (double)totalCount;
	sprintf(buffer,"Average image pixel %f\r\n",avgPixel);
	results += buffer;
	sprintf(buffer,"Low image pixel %d\r\n",lowPixel);
	results += buffer;
	sprintf(buffer,"High image pixel %d\r\n",highPixel);
	results += buffer;







//	load pattern into pm array and do pattern match with the same pattern

	loc = 0;
	for (i=0; i<patternHeight && i < YSIZE; ++i)
	{
		for (j=0; j<patternWidth && j < XSIZE; ++j)
		{
			r = patternBits[loc++];
			g = patternBits[loc++];
			b = patternBits[loc++];
			loc++;

			filter[i*XSIZE + j].real = (r + g + b) / 3 ;
			savefilter[i*XSIZE + j].real = (r + g + b) / 3 ;
			pm[(i+25)*XSIZE + (j+25)].real = (r + g + b) / 3 ;
		}
	}


//	flip horizontal and vertical

	int actualWidth,actualHeight;

	actualWidth = patternWidth;
	if (actualWidth > XSIZE) actualWidth = XSIZE;
	actualHeight = patternHeight;
	if (actualHeight > YSIZE) actualHeight = YSIZE;

	for (i=0; i<actualHeight; ++i)
	{
		for (j=0; j<actualWidth; ++j)
		{
			filter[i*XSIZE + j].real = savefilter[ (actualHeight-i-1)*XSIZE + (actualWidth-j-1) ].real;
		}
	}
	for (i=0; i<actualHeight; ++i)
	{
		for (j=0; j<actualWidth; ++j)
		{
			savefilter[i*XSIZE + j].real = filter[i*XSIZE + j].real;
		}
	}




	FFT2D(pm,XSIZE,YSIZE,1);
	FFT2D(filter,XSIZE,YSIZE,1);

	for (i=0; i<YSIZE; ++i)
	{
		for (j=0; j<XSIZE; ++j)
		{
			loc = i*XSIZE+j;
			result[loc].real = pm[loc].real * filter[loc].real - pm[loc].imag * filter[loc].imag ;
			result[loc].imag = pm[loc].imag * filter[loc].real + pm[loc].real * filter[loc].imag ;
		}
	}

	FFT2D(result,XSIZE,YSIZE,-1);

//	find the maxval of the resulting power;  this will be used as the denominator
//	in the percentage calculation of the actual pattern match

	maxval = 0;
	rowval = 0;
	colval = 0;
	for (i=0; i<YSIZE; ++i)
	{
		for (j=0; j<XSIZE; ++j)
		{
			val = sqrt(result[i*XSIZE+j].real*result[i*XSIZE+j].real+result[i*XSIZE+j].imag*result[i*XSIZE+j].imag);

			if (val > maxval)
			{
				maxval = val;
				rowval = i;
				colval = j;
			}
		}
	}

	bestpossible = maxval;

	sprintf(buffer,"Pattern: %f\r\n",bestpossible);
	results += buffer;


//	restore the filter

	for (i=0; i<YSIZE; ++i)
	{
		for (j=0; j<XSIZE; ++j)
		{
			filter[i*XSIZE+j].real = 0;
			filter[i*XSIZE+j].imag = 0;
		}
	}

	for (i=0; i<patternHeight; ++i)
	{
		for (j=0; j<patternWidth; ++j)
		{
			filter[i*XSIZE + j].real = savefilter[i*XSIZE + j].real;
		}
	}



	if (addNoise == 1)
	{
		for (i=0; i<YSIZE; ++i)
		{
			for (j=0; j<XSIZE; ++j)
			{
				adjust = (rand() % 40) - 20;
				data[i*XSIZE+j].real += adjust;
			}
		}
	}


//	level translate

	if (levelTranslate == 1)
	{
		int diff = highPixel - lowPixel;

		if (diff <= 0) diff = 255;
		double work;

		for (i=0; i<imageHeight && i < YSIZE; ++i)
		{
			for (j=0; j<imageWidth && j < XSIZE; ++j)
			{
				loc = i*XSIZE + j;

				data[loc].real -= (int)lowPixel;

				if (data[loc].real < 0) data[loc].real = 0;
				work = data[loc].real;
				work = work * 255.0 / (double)diff;
				data[loc].real = (int)work;
			}
		}
	}


//	generate results.bmp to show results of adding noise

	maxval = 0;
	rowval = 0;
	colval = 0;
	for (i=0; i<YSIZE; ++i)
	{
		for (j=0; j<XSIZE; ++j)
		{
	 		val = sqrt(data[i*XSIZE+j].real*data[i*XSIZE+j].real+data[i*XSIZE+j].imag*data[i*XSIZE+j].imag);
	 		if (val > maxval)
	 		{
				maxval = val;
				rowval = i;
				colval = j;
			}
		}
	}

	DWORD newarea = imageHeight * imageWidth * 4;

	BYTE *hpbits = new BYTE[newarea];

	loc = 0;
	for (i=0; i<imageHeight; ++i)
	{
		for (j=0; j<imageWidth; ++j)
		{
	 		val = sqrt(data[i*XSIZE+j].real*data[i*XSIZE+j].real+data[i*XSIZE+j].imag*data[i*XSIZE+j].imag);
			r = (int)(val / maxval * 255.0);
			g = r;
			b = r;
			hpbits[loc++] = r;
			hpbits[loc++] = g;
			hpbits[loc++] = b;
			hpbits[loc++] = 0;
		}
	}

	strcpy(outfile,"results.bmp");
	write24BitBmpFile(outfile, imageWidth,imageHeight,hpbits);

	delete [] hpbits;






//	perform correlation (pattern match)

	FFT2D(data,XSIZE,YSIZE,1);
	FFT2D(filter,XSIZE,YSIZE,1);

	for (i=0; i<YSIZE; ++i)
	{
		for (j=0; j<XSIZE; ++j)
		{
			loc = i*XSIZE+j;
			result[loc].real = data[loc].real * filter[loc].real - data[loc].imag * filter[loc].imag ;
			result[loc].imag = data[loc].imag * filter[loc].real + data[loc].real * filter[loc].imag ;
		}
	}

	FFT2D(result,XSIZE,YSIZE,-1);






//	search through the results and find the max power

	maxval = 0;
	rowval = 0;
	colval = 0;

	if (FullAnalysis == 1)
	{
		for (i=0; i<YSIZE; ++i)
		{
			for (j=0; j<XSIZE; ++j)
			{
				val = sqrt(result[i*XSIZE+j].real*result[i*XSIZE+j].real+result[i*XSIZE+j].imag*result[i*XSIZE+j].imag);
				power[i*XSIZE+j] = val;
				mark[i*XSIZE+j] = 0;


				sprintf(buffer,"%d %d",i,j);
				save.insert(ListPair(val,buffer));

				if (val > maxval)
				{
					maxval = val;
					rowval = i;
					colval = j;
				}
			}
		}
	}
	else
	{
		for (i=0; i<YSIZE; ++i)
		{
			for (j=0; j<XSIZE; ++j)
			{
				val = sqrt(result[i*XSIZE+j].real*result[i*XSIZE+j].real+result[i*XSIZE+j].imag*result[i*XSIZE+j].imag);
				power[i*XSIZE+j] = val;
				mark[i*XSIZE+j] = 0;

				if (val > maxval)
				{
					maxval = val;
					rowval = i;
					colval = j;
				}
			}
		}


		sprintf(buffer,"%d %d",rowval,colval);
		best.insert(ListPair(maxval,buffer));

	}


//	process saved power values and find local maximums
//	use to mark other potential matches (if image has multiple copies of pattern)

	FILE *fp;
	int islocalmax;
	int step,row,col,value,count,max;
	double detection;

	it = save.begin();
	lowpower = it->first;

	for (it=save.begin(); it!=save.end(); ++it)
	{
		strcpy(buffer,it->second.c_str());
		sscanf(buffer,"%d %d",&rowval,&colval);

		islocalmax = 1;
		for (value=3; value <= 11 && islocalmax; value += 2)
		{
			step = value / 2;
			for (col=-step; col<=step && islocalmax; ++col)
			{
				row = -step;
				if (rowval+row >= 0 && rowval+row < YSIZE && colval+col >= 0 && colval+col < XSIZE)
				{
					if (power[(rowval+row)*XSIZE+(colval+col)] >= power[rowval*XSIZE+colval]) islocalmax = 0;
				}
				row = step;
				if (rowval+row >= 0 && rowval+row < YSIZE && colval+col >= 0 && colval+col < XSIZE)
				{
					if (power[(rowval+row)*XSIZE+(colval+col)] >= power[rowval*XSIZE+colval]) islocalmax = 0;
				}
			}
			for (row=-step+1; row<=step-1 && islocalmax; ++row)
			{
				col = -step;
				if (rowval+row >= 0 && rowval+row < YSIZE && colval+col >= 0 && colval+col < XSIZE)
				{
					if (power[(rowval+row)*XSIZE+(colval+col)] >= power[rowval*XSIZE+colval]) islocalmax = 0;
				}
				col = step;
				if (rowval+row >= 0 && rowval+row < YSIZE && colval+col >= 0 && colval+col < XSIZE)
				{
					if (power[(rowval+row)*XSIZE+(colval+col)] >= power[rowval*XSIZE+colval]) islocalmax = 0;
				}
			}
		}


		detection = (it->first-lowpower)/(bestpossible);

		if (islocalmax && detection >= matchCutoff)
		{
			best.insert(ListPair(it->first,buffer));
		}
	}



//	display match results


	sprintf(buffer,"Max value %f @ col,row = %d,%d\r\n",maxval,colval,rowval);
	results += buffer;

	sprintf(buffer,"Center position = %f,%f\r\n",
		(double)colval-(double)patternWidth/2.0,
		(double)rowval-(double)patternHeight/2.0);
	results += buffer;


	double y1,y2,y3,dx,dy;


//	find true maximum power location by interpolating powers around the detected maximum

//	interpolate x

	if (rowval >= 0 && rowval < YSIZE && colval-1 >= 0 && colval-1 < XSIZE)
	{
		y1 = power[rowval*XSIZE+(colval-1)];
	}
	else
	{
		y1 = 0.0;
	}


	y2 = maxval;

	if (rowval >= 0 && rowval < YSIZE && colval+1 >= 0 && colval+1 < XSIZE)
	{
		y3 = power[rowval*XSIZE+(colval+1)];
	}
	else
	{
		y3 = 0.0;
	}


	dx = (y3-y1) / (2.0 * (2.0 * y2 - y1 - y3)) ;

	sprintf(buffer,"X: %f %f %f = %f, %f\r\n",y1,y2,y3,dx,(double)colval+dx);
	results += buffer;


//	interpolate y

	if (rowval-1 >= 0 && rowval-1 < YSIZE && colval >= 0 && colval < XSIZE)
	{
		y1 = power[(rowval-1)*XSIZE+(colval)];
	}
	else
	{
		y1 = 0.0;
	}

	y2 = maxval;

	if (rowval+1 >= 0 && rowval+1 < YSIZE && colval >= 0 && colval < XSIZE)
	{
		y3 = power[(rowval+1)*XSIZE+(colval)];
	}
	else
	{
		y3 = 0.0;
	}

	dy = (y3-y1) / (2.0 * (2.0 * y2 - y1 - y3)) ;

	sprintf(buffer,"Y: %f %f %f = %f, %f\r\n",y1,y2,y3,dy,(double)rowval+dy);
	results += buffer;


	double posx,posy,centerx,centery;

	centerx = (double)colval-(double)patternWidth/2.0+dx;
	centery = (double)rowval-(double)patternHeight/2.0+dy;
	sprintf(buffer,"Adjusted center position = %f,%f\r\n",centerx,centery);
	results += buffer;


	posx = (centerx - ((double)(imageWidth-1)/2.0)) * azfov / (double)imageWidth;
	posy = (((double)(imageHeight-1)/2.0) - centery) * elfov / (double)imageHeight;

	sprintf(buffer,"Adjusted center position = %f mrad,%f mrad\r\n",posx,posy);
	results += buffer;

	val = (maxval-lowpower)/(bestpossible-lowpower);
	sprintf(buffer,"Match %.3f\r\n",val);
	results += buffer;


//	open results.bmp and mark center of pattern with red X

	HBITMAP hOldBitmap;
	HDC hDC,hMemDC;
	HPEN hPen,hOldPen;
	int y,x,fontorient,fontsize,fontweight;
	HFONT hfont=NULL,holdfont=NULL;

	y = rowval-patternHeight/2;
	x = colval-patternWidth/2;

	strcpy(outfile,"results.bmp");
	m_bitmap = (HBITMAP)LoadImage(NULL,outfile,IMAGE_BITMAP,0,0,LR_LOADFROMFILE);
    GetObject(m_bitmap,sizeof(BITMAP),&bmImage);

	hPen = (HPEN)CreatePen(PS_SOLID,(int)2,RGB(255,0,0));
	hDC = GetDC(NULL);
	hMemDC = CreateCompatibleDC(hDC);
	SetBkMode(hMemDC,TRANSPARENT);

	hOldBitmap = (HBITMAP)SelectObject(hMemDC,m_bitmap);
	hOldPen = (HPEN)SelectObject(hMemDC,hPen);

	SetTextColor(hMemDC,RGB(255,0,0));

	fontsize = 9;
	fontorient = 0;
	fontweight = 500;
	hfont = CreateFont(0-fontsize,0,fontorient,fontorient,fontweight,0,0,0,ANSI_CHARSET,
		OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,PROOF_QUALITY, //DEFAULT_QUALITY,
		FF_DONTCARE, "Arial");

	holdfont = (HFONT)SelectObject(hMemDC,hfont);


	MoveToEx(hMemDC,x-6,y-6,NULL);
	LineTo(hMemDC,x+6,y+6);
	MoveToEx(hMemDC,x-6,y+6,NULL);
	LineTo(hMemDC,x+6,y-6);

	sprintf(buffer,"low %.1f",lowpower);
	TextOut(hMemDC,10,10,buffer,strlen(buffer));

	i = 0;

	for (rit=best.rbegin(); rit != best.rend() && i <10; ++rit)
	{
		strcpy(buffer,rit->second.c_str());
		sscanf(buffer,"%d %d",&rowval,&colval);


		val = (rit->first-lowpower)/(bestpossible-lowpower);
		if (val > 0.1)  // || 1)
		{
			Rectangle(hMemDC,colval-2,rowval-2,colval+2,rowval+2);
			sprintf(buffer,"%.3f",val);
			TextOut(hMemDC,colval,rowval+6,buffer,strlen(buffer));
		}
		++i;
	}
	newarea = imageHeight * imageWidth * 4;

	hpbits = new BYTE[newarea];
	printf("%d %d\n",GetBitmapBits(m_bitmap,newarea,hpbits),newarea);

	SelectObject(hMemDC,hOldBitmap);
	SelectObject(hMemDC,hOldPen);
	SelectObject(hMemDC,holdfont);

	DeleteDC(hMemDC);

	ReleaseDC(NULL,hDC);

	DeleteObject(hPen);
	DeleteObject(m_bitmap);




	BYTE *lpMem;

	lpMem = (BYTE *)GlobalLock(hMem);


	loc = 0;
	for (i=0; i<m_height; ++i)
	{
		for (j=0; j<m_width; ++j)
		{
			lpMem[loc] = hpbits[loc++];
			lpMem[loc] = hpbits[loc++];
			lpMem[loc] = hpbits[loc++];
			lpMem[loc] = hpbits[loc++];
		}
	}


	if (m_processed) DeleteObject(m_processed);

	m_processed = CreateBitmap(m_width, m_height, 1, 32, lpMem);
	GlobalUnlock(hMem);


	strcpy(outfile,"results.bmp");
	write24BitBmpFile(outfile, imageWidth,imageHeight,hpbits);

	delete [] hpbits;



	sprintf(buffer,"Output file is results.bmp\n");
	results += buffer;



	delete [] data;
	delete [] filter;
	delete [] result;
	delete [] change;
	delete [] pm;
	delete [] savefilter;
	delete [] power;
	delete [] mark;


}




void EnergyCentroid(char *infile, int level, string &results)
{
	int i,j,adjust,loc,choice,loc2;
	char outfile[512],buffer[256],aa[256],bb[256],cc[256];
	int r,g,b,rowval,colval;
	double val,maxval;
	double lowpower = 0.0;
	FILE *in;
	int gray;
	results = "";

	DWORD newarea = m_height * m_width * 4;

	BYTE *hpbits = new BYTE[newarea];
	GetBitmapBits(m_bitmap,newarea,hpbits);


//	calculate total size to avoid circular convolution

	int totalWidth,totalHeight;
	int imageWidth,imageHeight,patternWidth,patternHeight;
	int newwidth,newheight;

	imageWidth = m_width;
	imageHeight = m_height;
	patternWidth = 0;
	patternHeight = 0;

	totalWidth = imageWidth + patternWidth;
	totalHeight = imageHeight + patternHeight;

	adjust = (int)(log((double)totalWidth) / log(2.0));
	newwidth = (int)pow(2.0,(double)adjust);
	if (newwidth < totalWidth) adjust += 1;
	newwidth = (int)pow(2.0,(double)adjust);
	printf("width %d %d %d %d\n",imageWidth,patternWidth,totalWidth,newwidth);

	adjust = (int)(log((double)totalHeight) / log(2.0));
	newheight = (int)pow(2.0,(double)adjust);
	if (newheight < totalHeight) adjust += 1;
	newheight = (int)pow(2.0,(double)adjust);
	printf("height %d %d %d %d\n",imageHeight,patternWidth,totalHeight,newheight);


	XSIZE = newwidth;
	YSIZE = newheight;
	STORE = XSIZE * YSIZE;


	data = new COMPLEX[STORE];


	for (i=0; i<YSIZE; ++i)
	{
		for (j=0; j<XSIZE; ++j)
		{
			data[i*XSIZE+j].real = 0;
			data[i*XSIZE+j].imag = 0;
		}
	}


//	convert rgb to grayscale

	loc = 0;
	for (i=0; i<m_height; ++i)
	{
		for (j=0; j<m_width; ++j)
		{
			r = hpbits[loc++];
			g = hpbits[loc++];
			b = hpbits[loc++];
			loc++;
			gray = (r + g + b) / 3;

			if (gray >= level)
				data[i*m_width + j].real = 0xff;
			else
				data[i*m_width + j].real = 0x00;

		}
	}


//	do centroid

	unsigned long x,y,count;
	x = 0;
	y = 0;
	count = 0;

	loc = 0;
	for (i=0; i<m_height; ++i)
	{
		for (j=0; j<m_width; ++j)
		{
	 		val = data[i*m_width + j].real;
			if (val >= 0x80)
			{
				x += j;
				y += i;
				++count;
			}
		}
	}

	double fx=0,fy=0;

	if (count > 0)
	{
		fx = (double)x / (double)count;
		fy = (double)y / (double)count;
	}
	sprintf(buffer,"Centroid at %f,%f\r\n",fx,fy);
	results += buffer;

//	generate RESULTS.BMP that shows the best match location and center



	BYTE *lpMem;

	lpMem = (BYTE *)GlobalLock(hMem);


	loc = 0;
	loc2 = 0;
	for (i=0; i<m_height; ++i)
	{
		for (j=0; j<m_width; ++j)
		{
	 		val = data[i*m_width + j].real;
	 		r = val;
			g = r;
			b = r;

			hpbits[loc++] = r;
			hpbits[loc++] = g;
			hpbits[loc++] = b;
			hpbits[loc++] = 0;

			lpMem[loc2++] = r;
			lpMem[loc2++] = g;
			lpMem[loc2++] = b;
			lpMem[loc2++] = 0;
		}
	}


	if (m_processed) DeleteObject(m_processed);

	m_processed = CreateBitmap(m_width, m_height, 1, 32, lpMem);
	GlobalUnlock(hMem);





//	open RESULTS.BMP and mark center of pattern with red X

	HBITMAP hOldBitmap;
	HDC hDC,hMemDC;
	HPEN hPen,hOldPen;
	int fontorient,fontsize,fontweight;
	HFONT hfont=NULL,holdfont=NULL;

	y = fy;
	x = fx;

	hPen = (HPEN)CreatePen(PS_SOLID,(int)2,RGB(255,0,0));
	hDC = GetDC(NULL);
	hMemDC = CreateCompatibleDC(hDC);
	SetBkMode(hMemDC,TRANSPARENT);

	hOldBitmap = (HBITMAP)SelectObject(hMemDC,m_processed);
	hOldPen = (HPEN)SelectObject(hMemDC,hPen);

	SetTextColor(hMemDC,RGB(255,0,0));

	fontsize = 9;
	fontorient = 0;
	fontweight = 500;
	hfont = CreateFont(0-fontsize,0,fontorient,fontorient,fontweight,0,0,0,ANSI_CHARSET,
		OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,PROOF_QUALITY, //DEFAULT_QUALITY,
		FF_DONTCARE, "Arial");

	holdfont = (HFONT)SelectObject(hMemDC,hfont);


	MoveToEx(hMemDC,x-6,y-6,NULL);
	LineTo(hMemDC,x+6,y+6);
	MoveToEx(hMemDC,x-6,y+6,NULL);
	LineTo(hMemDC,x+6,y-6);

	sprintf(buffer,"%f,%f",fx,fy);
	TextOut(hMemDC,10,10,buffer,strlen(buffer));

	GetBitmapBits(m_processed,newarea,hpbits);

	SelectObject(hMemDC,hOldBitmap);
	SelectObject(hMemDC,hOldPen);
	SelectObject(hMemDC,holdfont);

	DeleteDC(hMemDC);

	ReleaseDC(NULL,hDC);

	DeleteObject(hPen);

	strcpy(outfile,infile);
	outfile[strcspn(outfile,".")] = 0;
	strcat(outfile,"_results.bmp");
	write24BitBmpFile(outfile, m_width, m_height, hpbits);

	delete [] hpbits;

	sprintf(buffer,"Output is %s\r\n",outfile);
	results += buffer;

	delete [] data;
}




//////////

int WINAPI WinMain(HINSTANCE hInstance,
                   HINSTANCE hPrevInstance,
                   LPSTR lpCmdLine,
                   int nCmdShow)
{
    WNDCLASSEX wcex;

    wcex.cbSize = sizeof(WNDCLASSEX);
    wcex.style          = CS_HREDRAW | CS_VREDRAW;
    wcex.lpfnWndProc    = WndProc;
    wcex.cbClsExtra     = 0;
    wcex.cbWndExtra     = 0;
    wcex.hInstance      = hInstance;
    wcex.hIcon          = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_APPLICATION));
    wcex.hCursor        = LoadCursor(NULL, IDC_ARROW);
    wcex.hbrBackground  = (HBRUSH)(COLOR_WINDOW+1);
    wcex.lpszMenuName   = NULL;
    wcex.lpszClassName  = szWindowClass;
    wcex.hIconSm        = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_APPLICATION));

    if (!RegisterClassEx(&wcex))
    {
        MessageBox(NULL,
            _T("Call to RegisterClassEx failed!"),
            _T("Win32 Guided Tour"),
            NULL);

        return 1;
    }

    hInst = hInstance; // Store instance handle in our global variable

    // The parameters to CreateWindow explained:
    // szWindowClass: the name of the application
    // szTitle: the text that appears in the title bar
    // WS_OVERLAPPEDWINDOW: the type of window to create
    // CW_USEDEFAULT, CW_USEDEFAULT: initial position (x, y)
    // 500, 100: initial size (width, length)
    // NULL: the parent of this window
    // NULL: this application does not have a menu bar
    // hInstance: the first parameter from WinMain
    // NULL: not used in this application
    HWND hWnd = CreateWindow(
        szWindowClass,
        szTitle,
        WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT,
        800, 600,
        NULL,
        NULL,
        hInstance,
        NULL
    );

    if (!hWnd)
    {
        MessageBox(NULL,
            _T("Call to CreateWindow failed!"),
            _T("Win32 Guided Tour"),
            NULL);

        return 1;
    }

    // The parameters to ShowWindow explained:
    // hWnd: the value returned from CreateWindow
    // nCmdShow: the fourth parameter from WinMain
    ShowWindow(hWnd,
        nCmdShow);
    UpdateWindow(hWnd);

    // Main message loop:
    MSG msg;
    while (GetMessage(&msg, NULL, 0, 0))
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }

    return (int) msg.wParam;
}

//
//  FUNCTION: WndProc(HWND, UINT, WPARAM, LPARAM)
//
//  PURPOSE:  Processes messages for the main window.
//
//  WM_PAINT    - Paint the main window
//  WM_DESTROY  - post a quit message and return
//
//
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    PAINTSTRUCT ps;
    HDC hdc,hMemDC;
	HBITMAP hOldBitmap;
	BITMAP bm,bmImage,bmPattern;
	double scale,yscale;
	BYTE *lpMem;
	char buffer[256];
	HBITMAP patternbitmap = NULL;
	FILE *out;
	int level;

    switch (message)
    {
	case WM_CREATE:
		CreateAMenu(hWnd);
		szFile[0] = 0;
		memset(&ofn,0,sizeof(OPENFILENAME));
		ofn.lStructSize = sizeof(OPENFILENAME);
		ofn.hwndOwner = hWnd;
		ofn.hInstance = (HINSTANCE)hInst;
		ofn.lpstrFilter = szReadFilter;
		ofn.nFilterIndex = 1;
		ofn.lpstrFile = szFile;
		ofn.nMaxFile = sizeof(szFile);
		ofn.lpstrFileTitle = szFileTitle;
		ofn.nMaxFileTitle = sizeof(szFileTitle);
		GetCurrentDirectory(sizeof(szDirName)-1,szDirName);
		ofn.lpstrInitialDir = szDirName;
		ofn.lCustData = 0L;

		strcpy(m_filterfile,"kirsch_8.txt");
		break;

	case WM_COMMAND:
		if (LOWORD(wParam) == ID_FILE_OPEN || LOWORD(wParam) == ID_FILE_OPEN2)
		{
			DialogType = 0;
			xAdjust = 1.0;
			if (LOWORD(wParam) == ID_FILE_OPEN2)
			{
				if (DialogBox(hInst,MAKEINTRESOURCE(IDD_DIALOG5),hWnd,(DLGPROC)SetParameters) == TRUE)
				{
					xAdjust = atof(m_parameter);
					if (xAdjust < 0.0001) xAdjust = 1.0;
				}
			}

			ofn.Flags = OFN_PATHMUSTEXIST | OFN_EXPLORER ;
			ofn.lpstrFilter = szReadFilter;

			if (GetOpenFileName(&ofn))
			{
				if (m_bitmap != NULL) DeleteObject(m_bitmap);
				m_bitmap = (HBITMAP)LoadImage(NULL,szFile,IMAGE_BITMAP,0,0,LR_LOADFROMFILE);
				strcpy(m_bmpfile,szFile);

				GetObject(m_bitmap,sizeof(BITMAP),&bm);
				m_width = bm.bmWidth;
				m_height = bm.bmHeight;

				if (hMem) GlobalFree(hMem);

				DWORD size = m_width * m_height * 4; // <2Gb
				hMem = GlobalAlloc(0, (SIZE_T)size);
				lpMem = (BYTE *)GlobalLock(hMem);

				sprintf(buffer,"%d",GetBitmapBits(m_bitmap,size,lpMem));
				m_processed = CreateBitmap(m_width, m_height, 1, 32, lpMem);
				GlobalUnlock(hMem);

				InvalidateRect(hWnd,NULL,TRUE);
			}
		}
		if (LOWORD(wParam) == ID_FILE_PATTERN)
		{
			ofn.Flags = OFN_PATHMUSTEXIST | OFN_EXPLORER ;
			ofn.lpstrFilter = szReadFilter;

			if (GetOpenFileName(&ofn))
			{
				if (patternbitmap != NULL) DeleteObject(patternbitmap);
				patternbitmap = (HBITMAP)LoadImage(NULL,szFile,IMAGE_BITMAP,0,0,LR_LOADFROMFILE);

			    GetObject(m_bitmap,sizeof(BITMAP),&bmImage);
				DWORD newarea = bmImage.bmHeight * bmImage.bmWidth * 4;
				BYTE *hpbits = new BYTE[newarea];
				GetBitmapBits(m_bitmap,newarea,hpbits);

				GetObject(patternbitmap,sizeof(BITMAP),&bmPattern);
				DWORD patarea = bmPattern.bmHeight * bmPattern.bmWidth * 4;
				BYTE *patbits = new BYTE[patarea];
				GetBitmapBits(patternbitmap,patarea,patbits);

				double azfov = 14.0;
				double elfov = 10.5;
				int addnoise = 1;
				int leveltranslate = 0;

				SetCursor(LoadCursor(NULL,IDC_WAIT));
				ShowCursor(TRUE);

				string results;
				PatternMatch(hpbits,bmImage.bmWidth,bmImage.bmHeight,
					patbits,bmPattern.bmWidth,bmPattern.bmHeight,
					azfov,elfov,
					addnoise,leveltranslate,
					0.5,
					results);

				delete [] hpbits;
				delete [] patbits;

				DeleteObject(patternbitmap);

				ShowCursor(FALSE);
				SetCursor(LoadCursor(NULL,IDC_ARROW));

				InvalidateRect(hWnd,NULL,TRUE);

				MessageBox(NULL,results.c_str(),"Info",MB_OK);

				out = fopen("results.txt","w");
				if (out != NULL)
				{
					fprintf(out,"%s",results.c_str());
					fclose(out);
				}

			}
		}
		if (LOWORD(wParam) == ID_FILE_PROCESS)
		{
			if (m_bitmap)
			{
				ofn.Flags = OFN_PATHMUSTEXIST | OFN_EXPLORER ;
				ofn.lpstrFilter = szFilterFilter;

				if (GetOpenFileName(&ofn))
				{
					strcpy(m_filterfile,szFile);
					string results;
					ProcessImage(m_bmpfile,m_filterfile,results);
					InvalidateRect(hWnd,NULL,TRUE);

					MessageBox(NULL,results.c_str(),"Info",MB_OK);

					out = fopen("results.txt","w");
					if (out != NULL)
					{
						fprintf(out,"%s",results.c_str());
						fclose(out);
					}
				}
			}
		}
		if (LOWORD(wParam) == ID_FILE_ENERGYCENTROID)
		{
			if (m_bitmap)
			{
				DialogType = 1;
				level = 128;
				if (DialogBox(hInst,MAKEINTRESOURCE(IDD_DIALOG5),hWnd,(DLGPROC)SetParameters) == TRUE)
				{
					level = atoi(m_parameter);
				}
				string results;
				EnergyCentroid(m_bmpfile,level,results);
				InvalidateRect(hWnd,NULL,TRUE);

				MessageBox(NULL,results.c_str(),"Info",MB_OK);

				out = fopen("results.txt","w");
				if (out != NULL)
				{
					fprintf(out,"%s",results.c_str());
					fclose(out);
				}
			}
		}
		if (LOWORD(wParam) == ID_EXIT)
		{
			PostQuitMessage(0);
		}
		break;
	case WM_SIZE:
		xClient = LOWORD(lParam);
		yClient = HIWORD(lParam);
		InvalidateRect(hWnd,NULL,TRUE);
		break;
    case WM_PAINT:

		if (m_bitmap != NULL)
		{
			scale = (double)xClient / (double)(m_width * 2.0 * xAdjust);
			yscale = (double)yClient / (double)m_height;
			if (yscale < scale) scale = yscale;

//			scale = 1.0;

			hdc = BeginPaint(hWnd, &ps);
			hMemDC = CreateCompatibleDC(hdc);
			SetStretchBltMode(hdc,STRETCH_ANDSCANS);
			hOldBitmap = (HBITMAP)SelectObject(hMemDC,m_bitmap);
			StretchBlt(hdc,0,0,(int)((double)m_width*scale*xAdjust),(int)((double)m_height*scale),
				hMemDC,0,0,(int)(m_width),(int)(m_height),SRCCOPY);
			SelectObject(hMemDC,m_processed);
			StretchBlt(hdc,m_width*scale*xAdjust,0,(int)((double)m_width*scale*xAdjust),(int)((double)m_height*scale),
				hMemDC,0,0,(int)(m_width),(int)(m_height),SRCCOPY);
			SelectObject(hMemDC,hOldBitmap);
			DeleteDC(hMemDC);
			EndPaint(hWnd, &ps);
		}
        break;
    case WM_DESTROY:
    	if (m_bitmap) DeleteObject(m_bitmap);
    	if (m_processed) DeleteObject(m_processed);
		if (hMem) GlobalFree(hMem);

        PostQuitMessage(0);
        break;
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
        break;
    }

    return 0;
}







void CreateAMenu(HWND hWnd)
{
    HMENU hMenu = CreateMenu();
    HMENU hSubMenu = CreatePopupMenu();

    AppendMenu(hSubMenu, MF_STRING, ID_FILE_OPEN, "&File open");
    AppendMenu(hSubMenu, MF_STRING, ID_FILE_OPEN2, "File open scale");
    AppendMenu(hSubMenu, MF_STRING, ID_FILE_PATTERN, "Pattern Match");
    AppendMenu(hSubMenu, MF_STRING, ID_FILE_PROCESS, "Filter");
    AppendMenu(hSubMenu, MF_STRING, ID_FILE_ENERGYCENTROID, "Energy centroid");
    AppendMenu(hSubMenu, MF_STRING, ID_EXIT, "&Exit");
    AppendMenu(hMenu, MF_STRING | MF_POPUP, (UINT)hSubMenu, "&File");

//    hSubMenu = CreatePopupMenu();
//    AppendMenu(hSubMenu, MF_STRING, ID_SHOW_ALL_ITEM, "Show &All Data");
//    AppendMenu(hSubMenu, MF_STRING, ID_SELECT_REPORT_ITEM, "S&eelect report");
//    AppendMenu(hMenu, MF_STRING | MF_POPUP, (UINT)hSubMenu, "&Reports");

    SetMenu(hWnd, hMenu);
}