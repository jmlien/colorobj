//------------------------------------------------------------------------------
//  Copyright 2007-2017 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#pragma once


//#include "gettime.h"
#include "model.h"

/*
#ifdef _WIN32
extern "C"{
#include "triangulate.h"
}
#else
#include "triangulate.h"
#endif
*/

//-----------------------------------------------------------------------------
// INPUTS
string output_filename;
string input_filename;
string param_method;

//-----------------------------------------------------------------------------
// Intermediate data
//mksum g_M;
//mksum& getM(){ return g_M; } //get singleton

static int mk_count = 0;

//-----------------------------------------------------------------------------
//read M+ from file
void save2file(FILE * fp);
void printUsage(char * name);

//-----------------------------------------------------------------------------
bool parseArg(int argc, char ** argv)
{
  double scale=1;
  Vector3d rot_angles;
	string env_filename;

  for(int i=1;i<argc;i++)
  {
        if(argv[i][0]=='-')
        {
            if(strcmp(argv[i],"-scale")==0) scale=atof(argv[++i]);
            else if(strcmp(argv[i],"-output")==0) output_filename=argv[++i];
            else if(strcmp(argv[i],"-input")==0) input_filename=argv[++i];
            else if(strcmp(argv[i],"-cgal")==0){ param_method="cgal"; }
            else if(strcmp(argv[i],"-libigl")==0){ param_method="libigl"; }
            else if(strcmp(argv[i],"-lp")==0){ param_method="lp"; }
        }
        else{
          input_filename = argv[i];
        }
  }

	if (input_filename.empty())
	{
		cerr << "! Error: No .obj file is given." << endl;
		printUsage(argv[0]);
		return false;
	}

  if(output_filename.empty())
  {
    //get output from input
    const char * filename=input_filename.c_str();
    const char * index=strrchr(filename, '.');
    char output[256]="";
    strncpy(output,filename,index-filename);
    output_filename=output;
    output_filename+="_tex.";
  }

  if(param_method.empty()) param_method="cgal";

  return true;
}

void printUsage(char * name)
{
    int offset=20;
    cerr<<"Usage: "<<name<<" [options] *.obj \n";

    offset=50;
    cerr<<"\n-- Report bugs to: Jyh-Ming Lien jmlien@cs.gmu.edu";
    cerr<<endl; //done
}

//-----------------------------------------------------------------------------

void save2file(ofstream& fout)
{

}

void save2file(FILE * fp)
{

}
