//------------------------------------------------------------------------------
//  Copyright 2007-2014 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "draw.h"
#include "dump.h"
#include <GLFW/glfw3.h>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

void save_buffer(GLFWwindow* window);

extern string output_filename;

//-----------------------------------------------------------------------------
//
// keyboard/mouse/etc... callback functions
//
//-----------------------------------------------------------------------------

static void error_callback(int error, const char* description)
{
    cerr << description << endl;
}

void draw(masc::model * m3d, flattener::flattened_mesh * m2d, int width, int height)
{
     GLFWwindow* window;
     glfwSetErrorCallback(error_callback);
     if (!glfwInit())
         exit(EXIT_FAILURE);
     window = glfwCreateWindow(width, height, "color obj", NULL, NULL);
     cout<<"create width="<<width<<" height="<<height<<endl;

     if (!window)
     {
         glfwTerminate();
         exit(EXIT_FAILURE);
     }

     glfwMakeContextCurrent(window);

     //while (!glfwWindowShouldClose(window))
     {
         int width, height;
         glfwGetFramebufferSize(window, &width, &height);
         glViewport(0, 0, width, height);

         glClear(GL_COLOR_BUFFER_BIT);

         glMatrixMode(GL_PROJECTION);
         glLoadIdentity();
         glOrtho(0, 1.f, 0, 1.f, 1.f, -1.f);

         glMatrixMode(GL_MODELVIEW);
         glLoadIdentity();

         glBegin(GL_TRIANGLES);
         for(int tid=0; tid<m3d->t_size; tid++)
         {
           auto& t3d = m3d->tris[tid];
           auto& t2d = m2d->get_triangle(tid);

           for(int i=0;i<3;i++)
           {
             glColor3dv(m3d->vertices[t3d.v[i]].c.get());
             //glNormal3d(0,0,1);
             glVertex2dv(t2d.cage[i].get());
           }
         }
         glEnd();

         save_buffer(window);
     }
     glfwDestroyWindow(window);
     glfwTerminate();
}

void save_buffer(GLFWwindow* window)
{
  int width, height;
  glfwGetFramebufferSize(window, &width, &height);
  cout<<"width="<<width<<" height="<<height<<endl;
  string file=output_filename+"png";
  dump_png(file.c_str(), width, height);
  cout<<"- Save to "<<file<<endl;
}
