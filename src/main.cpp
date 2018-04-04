//------------------------------------------------------------------------------
//  Copyright 2007-2017 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "main.h"
#include "gettime.h"
#include "angular-map.h"
#include "draw.h"

void save_obj(masc::model * m3d, flattener::flattened_mesh * m2d);

//
//
//
//
//   The MAIN function
//
//
//
//

int main(int argc, char ** argv)
{
    //
    //cout.precision(10);
	//srand48(0);
    masc::model m;

    //
    double start=getTime();
    if(!parseArg(argc,argv))
    {
        return 1;
    }
    double end=getTime();

    if (!m.build(input_filename)){
      cerr<<"! Error: Failed to create model from file ("<<input_filename<<")"<<endl;
      return 1;
    }

    cout<<"- Reading model takes "<<end-start<<" ms"<<endl;

    //unfold here
    flattener::flattened_mesh fm(m);
    flattener::ANGULAR_MAP_CREATOR * amap_creator=NULL;

    if(param_method=="cgal")
    {
      amap_creator = new flattener::CGAL_AMAP_Creator(flattener::CGAL_AMAP_Creator::Floater, flattener::CGAL_AMAP_Creator::Square_Arc_Length);
      //amap_creator = new flattener::CGAL_AMAP_Creator(flattener::CGAL_AMAP_Creator::Tutte, flattener::CGAL_AMAP_Creator::Square_Arc_Length);
    }

		//unwrap
    start=getTime();
    amap_creator->create(&m,fm);
    end=getTime();
    cout<<"- Flattening takes "<<end-start<<" ms"<<endl;

    //init openGL and draw with color
    draw(&m,&fm,1024,1024);

    //save the framebutter as texture
    save_obj(&m,&fm);

    return 0;
}

void save_obj(masc::model * m3d, flattener::flattened_mesh * m2d)
{
  string obj_filename = output_filename+"obj";
  string mlt_filename = obj_filename+".mtl";

  //save obj file
  {
    ofstream fout(obj_filename);
    if(!fout.good())
    {
      cerr<<"! Error: Failed to open "<<obj_filename<<endl;
      return;
    }

    fout<<"# Input: "<<input_filename<<"; Output: "<<obj_filename<<"\n";
    fout<<"# Created by Jyh-Ming Lien, jmlien@cs.gmu.edu, https://cs.gmu.edu/~jmlien"<<endl;
    fout<<"mtllib "<<mlt_filename<<"\n";

    //save vertices
    for(int vid=0; vid<m3d->v_size; vid++)
    {
      auto& v3d = m3d->vertices[vid];
      auto& v2d = m2d->get_vertex(vid);
      fout<<"v "<<v3d.p<<"\n";
      fout<<"vt "<<v2d.pos<<"\n";
    }

    //save faces
    for(int tid=0; tid<m3d->t_size; tid++)
    {
      auto& t3d = m3d->tris[tid];
      fout<<"f";
      for(int i=0;i<3;i++) fout<<" "<<t3d.v[i]+1<<"/"<<t3d.v[i]+1;
      fout<<"\n";
    }

    fout.close();
    cout<<"- Saved to "<<obj_filename<<endl;
  }

  //save mlt file
  {
    ofstream fout(mlt_filename);
    if(!fout.good())
    {
      cerr<<"! Error: Failed to open "<<mlt_filename<<endl;
      return;
    }


    fout<<"newmtl material_0\n"
        <<"Ka 0.200000 0.200000 0.200000\n"
        <<"Kd 1.000000 1.000000 1.000000\n"
        <<"Ks 1.000000 1.000000 1.000000\n"
        <<"Tr 1.000000\nillum 2\nNs 0.000000\n"
        <<"map_Kd "<<output_filename<<"png";
    fout.close();
    cout<<"- Saved to "<<mlt_filename<<endl;
  }
}
