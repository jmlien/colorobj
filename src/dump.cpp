# include "dump.h"
#include <iostream>
using namespace std;

#define PNG_DEBUG 3
#include <png.h>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

/****************************************************************************
Pack the 3 color components (red, green and blue) into one single number.
Since the raster represents each pixel as a single number, packcolor
needs to be used to convert your 3 components to a form imagesetpixel()
can use.
*****************************************************************************/
uint packcolors(uchar r, uchar g, uchar b)
{
    return(r << 24 | g << 16 | b << 8);
}

uchar unpack_red(uint rgb)
{
    return(0x000000ff & (rgb >> 24));
}
uchar unpack_green(uint rgb)
{
    return(0x000000ff & (rgb >> 16));
}
uchar unpack_blue(uint rgb)
{
    return(0x000000ff & (rgb >> 8));
}

/****************************************************************************
Initialize an image, takes the desired width and height and allocates
space for it.  If you don't use this routine to initialize your
image, imagesetpixel() and imagewritepixel() will return without
doing anything.  imagereadpixel() automatically initializes the
image itself.
*****************************************************************************/
bool imagemake(int width, int height, Image *image)
{
    image->comc = 0;
    image->comv = NULL;
    image->width = width;
    image->height = height;
    image->samplesperprimary = 255;
    image->pixels = (uint*)malloc(sizeof(uint)*image->width*image->height);
    memset(image->pixels,0,sizeof(uint)*image->width*image->height);
    return(true);
}

/****************************************************************************
This function will set a single pixel at coords x,y to the value of
color.  If the image hasn't been initialized, this routine will not
do anything.
*****************************************************************************/
void imagesetpixel(uint color, int x, int y, Image *image)
{
    if(image->pixels != NULL)
        image->pixels[y * image->width + x] = color;
}

/****************************************************************************
This function will read a ppm file into an image.  It will return
with true or false, depending on whether it successfully read the
ppm file.
*****************************************************************************/
bool imageread(const char *filename, Image *image)
{
    int     i;
    char        str[100];
    uchar       *buf;
    uint        size;
    FILE        *fp;

    image->comc = 0;
    image->comv = NULL;
    if(filename == NULL)
        fp = stdin;
    else
    {
        fp = fopen(filename,"r");
        if(fp == NULL)
        {
            fprintf(stderr,"imageread: %s - no such file or directory\n",filename);
            return(false);
        }
    }
    fgets(str,100,fp);
    if(strcmp(str,"P6\n") != 0)
    {
        fprintf(stderr,"imageread: input not in ppm format\n");
        return(false);
    }
    do
    {
        fgets(str,100,fp);
        if(str[0] == '#')
        {
            if(image->comv == NULL)
                image->comv = (char**)malloc(sizeof(char*));
            else image->comv = (char**)realloc((char**)image->comv,(int)sizeof(char*) *
                image->comc);
            image->comv[image->comc] = (char*)malloc(sizeof(char) * (int)strlen(str));
            strcpy(image->comv[image->comc],str);
            image->comc++;
        }
    } while(str[0] == '#');
    sscanf(str,"%i %i\n",&image->width,&image->height);
    size = image->width * image->height;
    fgets(str,100,fp);
    sscanf(str,"%i\n",&image->samplesperprimary);
    buf = (uchar *)malloc(sizeof(uchar) * size * 3);
    image->pixels = (uint*)malloc(sizeof(uint) * size);
    if(fread((uchar *)buf,sizeof(uchar),size * 3,fp) != size * 3)
    {
        fprintf(stderr,"imageread: file too short\n");
        return(false);
    }
    if(fp != stdin)
        fclose(fp);
    for(i=0; i < (int)size; i++)
        image->pixels[i] = packcolors(buf[i*3],buf[i*3+1],buf[i*3+2]);
    free((uchar*)buf);
    return(true);
}

/****************************************************************************
This function will write an image to a ppm file.  It will return
with true or false, depending on whether it successfully wrote the
ppm file.
*****************************************************************************/
bool imagewrite(const char *filename, Image image)
{
    int     i,j, index, indexp;
    uchar   *buf;
    int     size;
    FILE    *fp;

    if(image.pixels == NULL)
        return(false);
    fp = fopen(filename,"wb");
    if(fp == NULL)
        fp = stdout;

    size = image.width * image.height;
    buf = (uchar *)malloc(sizeof(uchar) * size * 3);

    fprintf(fp,"P6\n");
    fprintf(fp,"%i %i\n",image.width,image.height);
    fprintf(fp,"%i\n",image.samplesperprimary);

    for(i=0; i < (int)image.width; i++)
    {
        for(j=0; j < (int)image.height; j++)
        {
            index  = (i+j*image.width)*3;
            indexp = i+(image.height-j-1)*image.width;

            buf[index]   = image.pixels[indexp] >> 0  & 0xff; //G
            buf[index+1] = image.pixels[indexp] >> 8  & 0xff; //B
            buf[index+2] = image.pixels[indexp] >> 16 & 0xff; //R
        }
    }

    if(fwrite((uchar *)buf,sizeof(uchar),size * 3,fp) != (uint)size * 3)
    {
        fprintf(stderr,"imagewriteppm: couldn't write to output file\n");
        return(false);
    }


    if(fp != stdout)
        fclose(fp);
    free((uchar*)buf);

    return(true);
}


bool dump(const char * filename, int w, int h)
{
   //Guess filename
   const char * index=strrchr(filename, '.');
   char ppmFilename[256]="";
   strncpy(ppmFilename,filename,index-filename+1);
   strcat(ppmFilename,"ppm");

   Image new_image;
   imagemake(w, h, &new_image);
   glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, new_image.pixels);

   if( imagewrite(ppmFilename, new_image)==false ) return false;
   free(new_image.pixels);

   return true;
}



bool save_buffer_png(const char *filename, unsigned char * pixels, int w, int h)
{
    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);//
    if (!png){
       return false;
     }

    png_infop info = png_create_info_struct(png);//

    if (!info) {
        png_destroy_write_struct(&png, &info);//
        return false;
    }

    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        png_destroy_write_struct(&png, &info);
        return false;
    }

    // if (setjmp(png_jmpbuf(png)))
    // {
    //    /* If we get here, we had a problem with the file */
    //    fclose(fp);
    //    png_destroy_write_struct(&png, &info);
    //    cout << "png_destroy_write_struct failed" << endl;
    //    return false;
    // }


    png_init_io(png, fp);
    png_set_IHDR(png, info, w, h, 8 , PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    png_colorp palette = (png_colorp)png_malloc(png, PNG_MAX_PALETTE_LENGTH * sizeof(png_color));
    if (!palette) {
        fclose(fp);
        png_destroy_write_struct(&png, &info);
        return false;
    }

    png_set_PLTE(png, info, palette, PNG_MAX_PALETTE_LENGTH);
    png_write_info(png, info);
    png_set_packing(png);

    png_bytepp rows = (png_bytepp)png_malloc(png, h * sizeof(png_bytep));
    for (int i = 0; i < h; ++i)
        rows[i] = (png_bytep)(pixels + (h - i - 1) * w * 3);

    png_write_image(png, rows);

    png_write_end(png, info);
    png_free(png, palette);
    png_destroy_write_struct(&png, &info);

    fclose(fp);
    free(rows);
    return true;
}

bool dump_png(const char * filename, int w, int h)
{
   //Guess filename
   const char * index=strrchr(filename, '.');
   char pngFilename[256]="";
   strncpy(pngFilename,filename,index-filename+1);
   strcat(pngFilename,"png");

   unsigned char *pixels=new unsigned char[w*h*3];
   glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, pixels);

   if( save_buffer_png(pngFilename, pixels, w, h)==false )
   {
      cerr<<"! Error: Failed to save "<<pngFilename<<endl;
      delete [] pixels;
      return false;
    }

   delete [] pixels;

   return true;
}
