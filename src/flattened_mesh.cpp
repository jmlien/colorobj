//------------------------------------------------------------------------------
//  Copyright 2010-2017 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "flattened_mesh.h"
#include "Quaternion.h"
#include "intersection.h"
#include "simple_svg_1.0.0.hpp"

namespace flattener
{
  inline void backup(Point2d _from[3], Point2d _to[3])
  {_from[0]=_to[0]; _from[1]=_to[1]; _from[2]=_to[2];}

  inline float area(const Vector2d& a, const Vector2d& b, const Vector2d& c)
  {
      Vector2d v=b-a;
      Vector2d u=c-a;
      return fabs(v[0]*u[1]-v[1]*u[0])/2;
  }

  inline Point2d center(const Point2d& a, const Point2d& b, const Point2d& c)
  {
      return Point2d( (a[0]+b[0]+c[0])/3, (a[1]+b[1]+c[1])/3 );
  }

  //is c on the left of ab
  inline bool isLeft(const Vector2d& a, const Vector2d& b, const Vector2d& c)
  {
      const auto v=b-a;
      const auto u=c-a;
      return v[0]*u[1]-v[1]*u[0]>0;
  }

  //check if triangle xyz encloses point a
  inline bool isEnclosed(const Vector2d& a, const Vector2d& x, const Vector2d& y, const Vector2d& z)
  {
    if(isLeft(x,y,a)==false) return false;
    if(isLeft(y,z,a)==false) return false;
    if(isLeft(z,x,a)==false) return false;
    return true;
  }

  //check if triangle xyz encloses triangle abc
  inline bool isEnclosed(const Vector2d& a, const Vector2d& b, const Vector2d& c,
                         const Vector2d& x, const Vector2d& y, const Vector2d& z)
  {
    if(isEnclosed(a,x,y,z)==false) return false;
    if(isEnclosed(b,x,y,z)==false) return false;
    if(isEnclosed(c,x,y,z)==false) return false;
    return true;
  }

  inline void scale(Vector2d& x, Vector2d& y, Vector2d& z, float factor)
  {
    x=x*factor;
    y=y*factor;
    z=z*factor;
  }

  //rotate p around o for randian
  inline void rotate(float radian, const Vector2d& o, Vector2d& p)
  {

    Vector3d vec(p[0]-o[0],p[1]-o[1],0);
    Quaternion<double> q = Quaternion<double>::get(radian, Vector3d(0,0,1));
    vec=q.rotate(vec);
    p.set(vec[0]+o[0],vec[1]+o[1]);
  }

  //rotate triangle abc around o for randian
  inline void rotate_triangle(float radian, const Vector2d& o,
                               Vector2d& a, Vector2d& b, Vector2d& c)
  {
    rotate(radian,o,a);
    rotate(radian,o,b);
    rotate(radian,o,c);
  }

  inline void rotate_triangle(float radian, const Point2d& o,
                               Point2d& a, Point2d& b, Point2d& c)
  {
    const Vector2d _o(o.get());
    Vector2d _a(a.get());
    Vector2d _b(b.get());
    Vector2d _c(c.get());
    rotate(radian,_o,_a);
    rotate(radian,_o,_b);
    rotate(radian,_o,_c);
    a.set(_a.get());
    b.set(_b.get());
    c.set(_c.get());
  }

  //move edge xy so they touch vertex a or b
  inline void translate_edges(const Vector2d& a, const Vector2d& b,
                                    Vector2d& x,       Vector2d& y)
  {
    Vector2d v=y-x;
    Vector2d n(-v[1],v[0]);
    n=n.normalize();

    double d1=n*(a-x);
    double d2=n*(b-x);
    if(d1*d2<0)
    {
      cout<<"d1="<<d1<<" d2="<<d2<<endl;
      abort();
    }

    auto d=min(fabs(d1),fabs(d2));
    x=x+n*d;
    y=y+n*d;
    //cout<<" D="<<d;
  }

  //move edges of triangle xyz inward so they touch vertices a,b,c
  //assuming that triangle xyz contains triangle abc
  inline void translate_edges(const Vector2d& a, const Vector2d& b, const Vector2d& c,
                                    Vector2d& x,       Vector2d& y,       Vector2d& z)
  {
    Vector2d x1=x,x2=x;
    Vector2d y1=y,y2=y;
    Vector2d z1=z,z2=z;

    translate_edges(a,b,x1,y1);
    translate_edges(b,c,y2,z1);
    translate_edges(c,a,z2,x2);

    double p[2],q[2],r[2];
    char r1=LineLineInt(x1.get(),y1.get(),z2.get(),x2.get(),p); //this is new x
    char r2=LineLineInt(x1.get(),y1.get(),y2.get(),z1.get(),q); //this is new y
    char r3=LineLineInt(y2.get(),z1.get(),z2.get(),x2.get(),r); //this is new z

    //cout<<" r1="<<r1<<" r2="<<r2<<" r3="<<r3<<endl;

    if(r1=='0' || r2=='0' || r3=='0'){
       cout<<" r1="<<r1<<" r2="<<r2<<" r3="<<r3<<endl;
       abort();
     }

    x.set(p);
    y.set(q);
    z.set(r);
  }

  //compute angle between x and a (from x to a)
  inline float angle(const Vector2d& x, const Vector2d& a)
  {
    float radian = acos(a.normalize()*x.normalize());
    if(a[1]*x[0]-x[1]*a[0]<0) radian=-radian;
    return radian;
  }

  inline float align_triangle(const Vector2d& a, const Vector2d& b, const Vector2d& c,
                                    Vector2d& x,       Vector2d& y,       Vector2d& z)
  {
    float rxa = angle(x,a);
    float ryb = angle(y,b);
    float rzc = angle(z,c);

    Vector2d x_bkup=x, y_bkup=y, z_bkup=z;

    const Vector2d o(0,0);
    float r=0;

    if( rxa>0&&ryb>0&&rzc>0 )
    {
      r=min(rxa, min(ryb,rzc));
      rotate_triangle(r,o,x,y,z);
    }
    else if( rxa<0&&ryb<0&&rzc<0 )
    {
      r=max(rxa, max(ryb,rzc));
      rotate_triangle(r,o,x,y,z);
    }

    if(isEnclosed(a,b,c,x,y,z)==false) //failed...
    {
      x=x_bkup;
      y=y_bkup;
      z=z_bkup;
      return 0;
    }

    return r;
  }

  //rotate triangle xyz so that oa is aligned with ox, where o is world origin
  inline void align(const Vector2d& a, const Vector2d& b, const Vector2d& c,
                          Vector2d& x,       Vector2d& y,       Vector2d& z)
  {
    float radian = angle(x,a);
    const Vector2d o(0,0);
    rotate_triangle(radian,o,x,y,z);
  }

  //scale triangle xyz to enclose triangle abc
  inline void enclose(const Vector2d& a, const Vector2d& b, const Vector2d& c,
                            Vector2d& x,       Vector2d& y,       Vector2d& z)
  {
    //scale up
    while( isEnclosed(a,b,c,x,y,z)==false )
    {
      scale(x,y,z,2);
    }

    //cout<<"scaled up"<<endl;
    //cout<<"tri area="<<area(a,b,c)<<endl;
    //cout<<"before area="<<area(x,y,z)<<endl;

    //now, scale down
    int iteration=0;
    while(true)
    {
      cout<<"["<<iteration++<<"] ";
      translate_edges(a,b,c,x,y,z);
      // cout<<"before a=("<<a[0]<<","<<a[1]<<") b=("<<b[0]<<","<<b[1]<<") c=("<<c[0]<<","<<c[1]
      //     <<")\nx=("<<x[0]<<","<<x[1]<<") y=("<<y[0]<<","<<y[1]<<") z=("<<z[0]<<","<<z[1]<<")"<<endl;
      // cout<<endl;

      float amount=align_triangle(a,b,c,x,y,z);
      //cout<<" R="<<amount<<endl;

      if(fabs(amount)<1e-10) break;
    }

    //cout<<"after a=("<<a[0]<<","<<a[1]<<") b=("<<b[0]<<","<<b[1]<<") c=("<<c[0]<<","<<c[1]
    //    <<")\nx=("<<x[0]<<","<<x[1]<<") y=("<<y[0]<<","<<y[1]<<") z=("<<z[0]<<","<<z[1]<<")"<<endl;
    //cout<<"after area="<<area(x,y,z)<<endl;
    cout<<endl;

  }

  inline void print(string label, const Vector2d& a, const Vector2d& b, const Vector2d& c)
  {
    cout<<label<<": ("<<a[0]<<","<<a[1]<<") ("<<b[0]<<","<<b[1]<<") ("<<c[0]<<","<<c[1]<<")"<<endl;
  }


  inline void print(string label, const Point2d& a, const Point2d& b, const Point2d& c)
  {
    cout<<label<<": ("<<a[0]<<","<<a[1]<<") ("<<b[0]<<","<<b[1]<<") ("<<c[0]<<","<<c[1]<<")"<<endl;
  }

  float triangle::find_smallest_cage() //find the smallest enclosing triangle similar to cage
  {
    cout<<"----------------------------"<<endl;
    cout<<"find_smallest_cage"<<endl;
    //
    Point2d t_center = center(vertex[0],vertex[1],vertex[2]);
    Point2d c_center = center(cage[0],cage[1],cage[2]);
    Vector2d vertex_centered[3]={ vertex[0]-t_center, vertex[1]-t_center, vertex[2]-t_center};
    //

    float min_area=FLT_MAX;
    for(int i=0;i<3;i++)
    {
      //--------------------------------------------------------------------------
      Vector2d cage_centered[3]={ cage[0]-c_center, cage[1]-c_center, cage[2]-c_center};
      //
      align(vertex_centered[i],vertex_centered[(i+1)%3],vertex_centered[(i+2)%3],cage_centered[i],cage_centered[(i+1)%3],cage_centered[(i+2)%3]);
      enclose(vertex_centered[0],vertex_centered[1],vertex_centered[2],cage_centered[0],cage_centered[1],cage_centered[2]);
      //
      float cage_area=area(cage_centered[0],cage_centered[1],cage_centered[2]);
      if(cage_area<min_area)
      {
        min_area=cage_area;
        for(int j=0;j<3;j++) this->cage_min[j].set(cage_centered[j].get());
      }
      //--------------------------------------------------------------------------
    }

    //
    float input_area=area(vertex_centered[0],vertex_centered[1],vertex_centered[2]);
    cout<<"input area="<<input_area<<endl;
    cout<<"cage  area="<<min_area<<"(X"<<min_area/input_area<<")"<<endl;
    print("input",vertex_centered[0],vertex_centered[1],vertex_centered[2]);
    print("cage",cage_min[0],cage_min[1],cage_min[2]);

    return (cage_min[0]-cage_min[1]).norm()/(cage[0]-cage[1]).norm(); //min_area/input_area;
  }

  void triangle::scale_cage(float s)
  {
    for(int i=0;i<3;i++) for(int j=0;j<2;j++) cage[i][j]*=s;
  }

  void triangle::move_to_cage()
  {
    //find transform from min_cage to cage

    //find the center of cage, min_cage, and triangle
    Point2d t_center  = center(vertex[0],vertex[1],vertex[2]);
    Point2d c_center  = center(cage[0],cage[1],cage[2]);
    Point2d cm_center = center(cage_min[0],cage_min[1],cage_min[2]);
    Vector2d vertex_new[3]={ vertex[0]-t_center, vertex[1]-t_center, vertex[2]-t_center};
    //known: cage_min contains vertex_centered


    //translation: cm_center to c_center
    Vector2d T = c_center-cm_center;

    //rotation (cage_min[0]-cm_center) to (cage[0]-c_center)
    Vector2d u=(cage_min[0]-cm_center);
    Vector2d v=(cage[0]-c_center);
    float r=angle(u,v);

    //move this triangle using the given transform
    const Point2d O(0,0);
    rotate_triangle(r,cm_center-O, vertex_new[0],vertex_new[1],vertex_new[2]);
    for(int i=0;i<3;i++) vertex[i]=O+vertex_new[i]+T;

    //Let us move the cage_min as well
    Vector2d cage_new[3]={cage_min[0]-O, cage_min[1]-O, cage_min[2]-O};
    rotate_triangle(r,cm_center-O, cage_new[0],cage_new[1],cage_new[2]);
    for(int i=0;i<3;i++) cage_min[i]=cage_new[i]+T;
  }


  //------------------------------------------------------------------
  //put the triangles in cage
  //------------------------------------------------------------------
  void flattened_mesh::place_in_cage()
  {
    //use the angular map to find the intial placement
    float max_scale=0;
    for(auto& f : triangles)
    {
      float scale=f.find_smallest_cage();
      if(scale>max_scale) max_scale=scale;
    }

    cout<<"max_scale="<<max_scale<<endl;

    //not scale all cages by max_scale
    for(auto& f : triangles)
    {
      f.scale_cage(max_scale);
      f.move_to_cage();
    }
  }

  //------------------------------------------------------------------
  //build edges
  //------------------------------------------------------------------

  void flattened_mesh::build_edges()
  {
    //for each edge in masc::model, we build an edge
    for(uint ie=0;ie<pm->e_size;ie++)
    {
      //
      auto& pm_e=pm->edges[ie];
      //
      edge& e = this->edges[ie];
      e.eid=ie;
      e.fid=pm_e.fid;

      if(e.fid.size()>2)
      {
        cerr<<"! Error: Edge "<<ie<<" is non-manifold, incident to "<<e.fid.size()<<" faces"<<endl;
        abort();
      }


      if(e.fid.size()==2 && e.fid.front()==e.fid.back())
      {
        e.fid.pop_back();
        pm_e.fid.pop_back();
      }


      //
      for(auto fid: e.fid)
      {
        auto& f = pm->tris[fid];
        for(int i=0;i<3;i++)
        {
          int ni=(i+1)%3;
          if( (f.v[i]==pm_e.vid[0] && f.v[ni]==pm_e.vid[1]) || (f.v[i]==pm_e.vid[1] && f.v[ni]==pm_e.vid[0]))
          {
            e.vid.push_back(i);
            e.vid.push_back(ni);
            break;
          }//
        }//end for i
      }//end for fid

    }//end for ie

    for(uint ie=0;ie<pm->e_size;ie++)
    {
       if(isValidEdge(ie)==false)
       {
         cout<<" edge "<<ie<<" is invalid"<<endl;
       }
    }
  }

  void flattened_mesh::build_vertices()
  {
    for(uint iv=0;iv<pm->v_size;iv++)
    {
      //vertices[iv].eid=pm->vertices[iv].m_e;
      auto& edges=pm->vertices[iv].m_e;
      vertices[iv].eid=vector<uint>(edges.begin(), edges.end());
      for(auto eid:edges)
      {
        if(isBorderEdge(eid)){ vertices[iv].is_interior=false; break;}
      }//end for eid
    }//end for iv
  }


  void flattened_mesh::lapaician_smoothing(int max_iteration, float step_size)
  {
    //for each face, move to the average of the its neighboring positions
    int success_condense_count=0;
    int iteration=0;
    uint fsize=this->triangles.size();
    do{

      success_condense_count=0;

      //for each face contract to root
      for(uint fid=0;fid<fsize;fid++ )
        if(lapaician_smoothing(fid,step_size)) success_condense_count++;

      // for(uint fid=0;fid<fsize;fid++ )
      //   if(update_face_orientation(fid,PI/180)) success_condense_count++;

      if(iteration%100==0)
        export_to_svg("flatten.svg");

      cout<<"-[Lapaician_smoothing] iteration "<<iteration++<<endl;
    }
    while(success_condense_count>0 && iteration<max_iteration);//stop until no face can be condesnsed
  }

  //------------------------------------------------------------------
  // condensation methods
  //------------------------------------------------------------------

  //iterative condensation to incident vertices
  void flattened_mesh::condense_to_vertices(int max_iteration, float step_size)
  {
    int success_condense_count=0;
    int iteration=0;
    uint fsize=this->triangles.size();
    uint vsize=this->vertices.size();

    do{

      success_condense_count=0;
      //update vertex locations
      for(uint vid=0;vid<vsize;vid++ )
      {
        if(this->vertices[vid].is_interior)
          update_vertex_position(vid);
      }

      //for each face contract to incident vertex
      for(uint fid=0;fid<fsize;fid++ )
        if(condense_to_vertices(fid,step_size)) success_condense_count++;

      for(uint fid=0;fid<fsize;fid++ )
       if(update_face_orientation(fid,PI/180)) success_condense_count++;

      if(iteration%100==0) export_to_svg("flatten.svg");
      cout<<"-[Condense by collapsing to vertices] iteration "<<iteration++<<endl;
    }
    while(success_condense_count>0 && iteration<max_iteration);//stop until no face can be condesnsed
  }

  void flattened_mesh::condense_to_root_with_edge_collapsing(int max_iteration, float step_size)
  {
    int success_condense_count=0;
    int iteration=0;
    uint fsize=this->triangles.size();

    //find root id
    vector< pair<float, uint> > sorted_fids;
    for(auto& f: triangles)
    {
      float len=0;
      for(int i=0;i<3;i++) len+=getEdgeLength(f.eid[i]);
      sorted_fids.push_back(make_pair(len,f.fid));
    }
    sort(sorted_fids.begin(),sorted_fids.end());

    //identify rootid
    uint rootid=sorted_fids.front().second;

    //
    do{

      success_condense_count=0;

      //for each face contract to root
      for(uint fid=0;fid<fsize;fid++ )
        if(condense_by_edge_collapsing(fid,step_size)) success_condense_count++;

      for(uint fid=0;fid<fsize;fid++ )
        if(condense_to_face(fid,rootid, step_size)) success_condense_count++;

      for(uint fid=0;fid<fsize;fid++ )
       if(update_face_orientation(fid,PI/180)) success_condense_count++;

      if(iteration%100==0) export_to_svg("flatten.svg");
      cout<<"-[Condense to root with edge collapsing] iteration "<<iteration++<<endl;
    }
    while(success_condense_count>0 && iteration<max_iteration);//stop until no face can be condesnsed
  }


  //------------------------------------------------------------------
  //directly iterative condensation
  //------------------------------------------------------------------

  void flattened_mesh::condense_to_root(int max_iteration, float step_size)
  {
    int success_condense_count=0;
    int iteration=0;

    //sort all faces
    vector< pair<float, uint> > sorted_fids;
    for(auto& f: triangles)
    {
      float len=0;
      for(int i=0;i<3;i++) len+=getEdgeLength(f.eid[i]);
      sorted_fids.push_back(make_pair(len,f.fid));
    }
    sort(sorted_fids.begin(),sorted_fids.end());

    //identify rootid
    uint rootid=sorted_fids.front().second;

    do
    {
      success_condense_count=0;

      //for each face contract to root
      for(auto it=sorted_fids.rbegin();it!=sorted_fids.rend();it++)
        if(condense_to_face(it->second,rootid, step_size)) success_condense_count++;

      for(auto it=sorted_fids.rbegin();it!=sorted_fids.rend();it++)
        if(update_face_orientation(it->second,PI/180)) success_condense_count++;

      if(iteration%100==0) export_to_svg("flatten.svg");
      cout<<"-[Condense to root] iteration "<<iteration++<<endl;
    }
    while(success_condense_count>0 && iteration<max_iteration);//stop until no face can be condesnsed
  }

  //------------------------------------------------------------------
  // condense by contracting edges
  //------------------------------------------------------------------
  void flattened_mesh::condense_by_edge_collapsing(int max_iteration, float step_size)
  {
    int success_condense_count=0;
    int iteration=0;
    uint fsize=this->triangles.size();
    do{

      success_condense_count=0;

      //for each face contract to root
      for(uint fid=0;fid<fsize;fid++ )
        if(condense_by_edge_collapsing(fid,step_size)) success_condense_count++;

      for(uint fid=0;fid<fsize;fid++ )
        if(update_face_orientation(fid,PI/180)) success_condense_count++;

      if(iteration%100==0) export_to_svg("flatten.svg");
      cout<<"-[Condense by edge collapsing] iteration "<<iteration++<<endl;
    }
    while(success_condense_count>0 && iteration<max_iteration);//stop until no face can be condesnsed
  }

  //------------------------------------------------------------------
  //move a given face to the center of its neighbors
  //return true if successful
  //------------------------------------------------------------------
  bool flattened_mesh::lapaician_smoothing(uint fid, const float step_size)
  {

    auto& f=triangles[fid];
    Vector2d displacement;
    const float min_len=20;
    Point2d f_center=center(f.cage_min[0],f.cage_min[1],f.cage_min[2]);
    Point2d nei_center(0,0);

    for(uint eid : f.eid)
    {
      if( isBorderEdge(eid) ) return true; //do not need to move the border triangles

      edge& e=edges[eid];
      uint ofid = e.otherf(fid);
      if(ofid==UINT_MAX) continue;

      auto& of=triangles[ofid];
      Point2d of_center=center(of.cage_min[0],of.cage_min[1],of.cage_min[2]);
      nei_center[0] += of_center[0];
      nei_center[1] += of_center[1];
    }
    nei_center[0]/=3;
    nei_center[1]/=3;

    displacement=(nei_center-f_center)*step_size;
    return attempt_translate(displacement,f);
  }

  //------------------------------------------------------------------
  //condense a given face to its neighbors
  //return true if successful
  //------------------------------------------------------------------
  bool flattened_mesh::condense_by_edge_collapsing(uint fid, const float step_size)
  {
    auto& f=triangles[fid];
    Vector2d displacement;
    const float min_len=20;

    for(uint eid : f.eid)
    {
      if( isBorderEdge(eid) ) continue;

      vector<Point2d> geo;
      getEdgeGeometry(eid, geo);
      edge& e=edges[eid];

      Vector2d v1=geo[3]-geo[0];
      Vector2d v2=geo[2]-geo[1];

      if(e.fid.front()!=fid) //geo[0]&geo[1] are not vertices of fid
      {
        v1=-v1;
        v2=-v2;
      }

      if(v1.norm()>min_len) displacement=displacement+v1;
      if(v2.norm()>min_len) displacement=displacement+v2;
    }

    displacement=displacement*step_size;
    return attempt_translate(displacement,f);
  }

  //------------------------------------------------------------------
  //condense a given face to a given face
  //return true if successful
  //------------------------------------------------------------------
  bool flattened_mesh::condense_to_face(uint fid, uint rootid, const float step_size)
  {
    if(fid==rootid) return true;

    Point2d root_center, f_center;
    list<Vector2d>  displacements;

    //set<uint> affected_eids(f.eid.begin(),f.eid.end());
    auto& root=triangles[rootid];
    auto& f=triangles[fid];

    root_center=center(root.cage_min[0],root.cage_min[1],root.cage_min[2]);
    f_center=center(f.cage_min[0],f.cage_min[1],f.cage_min[2]);
    Vector2d vec_to_root=(root_center-f_center).normalize();

    //find edges that would pull this triangle to root
    //and identify displacement
    for(uint eid : f.eid)
    {
      if( isBorderEdge(eid) ) continue;

      Vector2d dir;
      vector<Point2d> geo;
      getEdgeGeometry(eid, geo);
      edge& e=edges[eid];

      Vector2d v1=geo[3]-geo[0];
      Vector2d v2=geo[2]-geo[1];
      if(e.fid.front()!=fid) //geo[0]&geo[1] are vertices of fid
      {
        v1=-v1;
        v2=-v2;
      }

      double d1=v1*vec_to_root;
      double d2=v2*vec_to_root;
      if(d1>0 && d2>0) dir=(v1+v2)/2;
      else if(d1>0) dir=d1;
      else if(d2>0) dir=d2;
      //if(dir.norm()>20)
      //displacement=displacement+dir;
      if(dir.normsqr()>0) displacements.push_back(dir);
    }

    //mix the displacements
    if(displacements.size()==2)
    {
        Vector2d v=displacements.front()+displacements.back();
        displacements.push_front(v); //try this first
    }
    else if(displacements.size()==3)
    {
        const Vector2d& v1=displacements.front();
        const Vector2d& v2=*(++displacements.begin());
        const Vector2d& v3=displacements.back();
        Vector2d v4=v1+v2+v3;
        Vector2d v5=v1+v2;
        Vector2d v6=v2+v3;
        Vector2d v7=v1+v3;
        displacements.push_front(v7);
        displacements.push_front(v6);
        displacements.push_front(v5);
        displacements.push_front(v4); //try this first
    }

    for(auto& displacement:displacements)
    {
      if(displacement.norm()<20) continue; //step too small
      displacement=displacement*step_size;
      if(attempt_translate(displacement,f));
        return true;
    }
    return false;
  }

  //------------------------------------------------------------------
  //condense a given face to a given face
  //return true if successful
  //------------------------------------------------------------------
  bool flattened_mesh::condense_to_vertices(uint fid, const float step_size)
  {
    Vector2d displacement;
    auto& f=triangles[fid];
    Point2d f_center=center(f.cage_min[0],f.cage_min[1],f.cage_min[2]);

    // for(uint eid : f.eid)
    // {
    //   if( isBorderEdge(eid) ) return true; //do not need to move the border triangles
    // }

    //find edges that would pull this triangle to root
    //and identify displacement
    float areasum=0;
    for(uint vid : f.vid)
    {
      auto& v=this->vertices[vid];
      if( !v.is_interior ) continue; //for interior vertex only
      displacement=displacement+(v.pos-f_center); //*v.area;
      areasum+=v.area;
    }

    //displacement=displacement/areasum; //normalize

    //if(displacement.norm()<20) return true; //step too small
    displacement=displacement*step_size;

    return attempt_translate(displacement,f);
  }

  //rotate the given face so that its edge modules are perpendicular to the edges
  bool flattened_mesh::update_face_orientation(uint fid, float step_size)
  {
    auto& f=triangles[fid];

    //
    double radian_sum=0;

    int count=0;
    for(uint eid : f.eid)
    {
      if( isBorderEdge(eid) ) continue;

      vector<Point2d> geo;
      getEdgeGeometry(eid, geo);
      edge& e=edges[eid];

      Vector2d v1=geo[3]-geo[0];
      Vector2d v2=geo[2]-geo[1];
      Vector2d u =geo[1]-geo[0];

      if(e.fid.front()!=fid) //geo[0]&geo[1] are not vertices of fid
      {
        v1=-v1;
        v2=-v2;
        u=geo[3]-geo[2];
      }
      u=u.normalize();
      double d1=v1.norm();
      double d2=v2.norm();

      //
      if(d1>0)
      {
        v1=v1/d1; //normalize
        double radian=angle(u,v1);
        radian+=((radian>0)?(-PI/2):(PI/2));
        radian_sum+=radian;
        count++;
      }

      if(d2>0)
      {
        v2=v2/d2; //normalize
        double radian=angle(u,v2);
        radian+=((radian>0)?(-PI/2):(PI/2));
        radian_sum+=radian;
        count++;
      }
    }

    auto radian=radian_sum/count;//(radian_sum>0)?step_size:-step_size;
    return attempt_rotate(radian,f);
  }

  //------------------------------------------------------------------
  //reduce edge length level by level
  //------------------------------------------------------------------
  void flattened_mesh::condense_level()
  {
    int rootid=findRoot();
    map<uint,list<uint> > levelsets;
    uint maxlevel=buildLevelSets(rootid, levelsets);

    for(uint level=maxlevel;level>0;level--)
    {
      cout<<"condense level "<<level<<" with "<<levelsets[level].size()<<" triangles"<<endl;
      int iteration=0;
      //while(condense_level(levelsets[level]))
      {
        cout<<"iteration="<<iteration++<<endl;
      }
    }//end for
  }

  //condense this level to its parent level
  bool flattened_mesh::condense_level(list<uint>& levelset)
  {
    float step_size=0.1;

    map<uint,Vector2d> displacement;
    set<uint> affected_eids;

    for(uint fid : levelset)
    {
      list<uint> low_level_eids;
      getLowLeveLNeighbors(fid,low_level_eids);
      Vector2d dir;
      cout<<"low_level_eids size="<<low_level_eids.size()<<endl;
      for(uint eid:low_level_eids)
      {
        vector<Point2d> geo;
        getEdgeGeometry(eid, geo);
        edge& e=edges[eid];
        if(e.fid.front()==fid) //geo[0]&geo[1] are vertices of fid
        {
          dir=dir+(geo[3]-geo[0])+(geo[2]-geo[1]);
        }
        else //geo[2]&geo[3] are vertices of fid
        {
          dir=dir+(geo[0]-geo[3])+(geo[1]-geo[2]);
        }
      }//loop through
      //
      displacement[fid]=dir.normalize()*step_size;
      cout<<"displacement["<<fid<<"]="<<displacement[fid]<<endl;
    }//end for fid

    //move the triangles
    for(uint fid : levelset)
    {
      auto& f=triangles[fid];
      affected_eids.insert(f.eid.begin(),f.eid.end());
      for(int i=0;i<3;i++)
      {
        f.vertex[i]=f.vertex[i]+displacement[fid];
        f.cage_min[i]=f.cage_min[i]+displacement[fid];
      }
    }

    //check if the move is valid
    bool OK=true;
    for(uint eid : affected_eids)
    {
      if(isValidEdge(eid)==false)
      {
        OK=false;
        break;
      }
    }

    if(OK==false) //failed, move back
    {
      for(uint fid : levelset)
      {
        auto& f=triangles[fid];
        for(int i=0;i<3;i++)
        {
          f.vertex[i]=f.vertex[i]-displacement[fid];
          f.cage_min[i]=f.cage_min[i]-displacement[fid];
        }
      }//end for fid
    }//end if OK

    return OK;
  }

  //root is a face with min edge length and 3 neighbors
  uint flattened_mesh::findRoot()
  {
    uint root=UINT_MAX;
    float min_edge_length=FLT_MAX;

    for(auto& f: triangles)
    {
      float len=0;
      for(int i=0;i<3;i++) len+=getEdgeLength(f.eid[i]);
      if(len<min_edge_length)
      {
        min_edge_length=len;
        root=f.fid;
      }
    }

    return root;
  }

  //
  //Note: This does not work as intended as the density of these triangles
  //      differs a lot in different regions
  //
  //order all facets into levels
  //return max level
  uint flattened_mesh::buildLevelSets(uint rootid, map<uint,list<uint> >& levelsets)
  {
    list<uint> open;
    open.push_back(rootid);
    triangles[rootid].level=0;
    uint max_level=0;

    while(!open.empty())
    {
      uint fid=open.front();
      open.pop_front();
      auto& f=triangles[fid];
      //put f into level set
      levelsets[f.level].push_back(fid);
      max_level=f.level;

      //propagate
      for(int i=0;i<3;i++)
      {
        auto& e=edges[f.eid[i]];
        uint of = e.otherf(fid);
        if(of==UINT_MAX) continue; //not valid
        if(triangles[of].level<=f.level+1) continue;
        //
        triangles[of].level=f.level+1;
        open.push_back(of);
      }//end for i

    }//end while

    return max_level;
  }

  //check if a given edge is valid
  bool flattened_mesh::isValidEdge(uint eid) const
  {
    //border edge has not geometry and therefore is always intersection free
    if(isBorderEdge(eid)) return true;

    vector<Point2d> geo;
    list<uint> incident_fids, incident_eids;
    getEdgeGeometry(eid, geo);

    incidentObjects(eid, incident_fids, incident_eids);

    //check if geo touches the triangles incident to this edge
    //(1) make sure geo[0]geo[3] and geo[1]geo[2] does not cross each other
    if( SegSegInt(geo[0].get(),geo[3].get(),geo[1].get(),geo[2].get()) ){
      //cout<<"A";
      // cout<<" geo size="<<geo.size()<<" ";
      // cout<<" geo[0]="<< geo[0]<<" ";
      // cout<<" geo[1]="<< geo[1]<<" ";
      // cout<<" geo[2]="<< geo[2]<<" ";
      // cout<<" geo[3]="<< geo[3]<<" ";
      return false;
    }

    //move geo[0] & geo[3] closer to each other by a small amount
    float small=1e-8;
    {
      Vector2d vec= (geo[3]-geo[0]).normalize()*small;
      geo[0] = geo[0]+vec;
      geo[3] = geo[3]-vec;
    }
    //move geo[1] & geo[2] closer to each other by a small amount
    {
      Vector2d vec= (geo[2]-geo[1]).normalize()*small;
      geo[1] = geo[1]+vec;
      geo[2] = geo[2]-vec;
    }

    //(2) make sure geo[0]geo[3] and geo[1]geo[2] does not cross incident trangles
    for(uint ifid:incident_fids)
    {
      auto& t=triangles[ifid];
      for(int i=0;i<3;i++)
      {
        int ni=(i+1)%3;
        auto& p1=t.cage_min[i];
        auto& p2=t.cage_min[ni];
        if( SegSegInt(p1.get(),p2.get(),geo[1].get(),geo[2].get()) ){ return false;}
        if( SegSegInt(p1.get(),p2.get(),geo[0].get(),geo[3].get()) ){ return false; }
      }
    }

    //(3) make sure geo[0]geo[3] and geo[1]geo[2] does not cross other edges
    for(uint ieid:incident_eids)
    {
      if(isBorderEdge(ieid)) continue;

      vector<Point2d> geo2;
      getEdgeGeometry(ieid, geo2);
      if( SegSegInt(geo[0].get(),geo[3].get(),geo2[0].get(),geo2[3].get()) ){ return false;}
      if( SegSegInt(geo[0].get(),geo[3].get(),geo2[1].get(),geo2[2].get()) ){ return false;}
      if( SegSegInt(geo[1].get(),geo[2].get(),geo2[0].get(),geo2[3].get()) ){ return false;}
      if( SegSegInt(geo[1].get(),geo[2].get(),geo2[1].get(),geo2[2].get()) ){ return false;}
    }

    return true;
  }

  bool flattened_mesh::attempt_translate(const Vector2d& displacement, triangle& f)
  {
    //move the triangle
    for(int i=0;i<3;i++)
    {
      f.vertex[i]=f.vertex[i]+displacement;
      f.cage_min[i]=f.cage_min[i]+displacement;
    }

    //check if the move is valid
    bool OK=true;
    for(uint eid : f.eid)
    {
      if(isValidEdge(eid)==false)
      {
        OK=false;
        break;
      }
    }

    if(OK==false) //failed, move back
    {
      for(int i=0;i<3;i++)
      {
        f.vertex[i]=f.vertex[i]-displacement;
        f.cage_min[i]=f.cage_min[i]-displacement;
      }
    }//end if OK

    //if(OK) cout<<"OK"<<endl; else  cout<<"NO"<<endl;
    return OK;
  }

  bool flattened_mesh::attempt_rotate(float radian, triangle& f)
  {
    //rotate this triangle
    Point2d f_center=center(f.cage_min[0],f.cage_min[1],f.cage_min[2]);
    rotate_triangle(radian,f_center,f.vertex[0],f.vertex[1],f.vertex[2]);
    rotate_triangle(radian,f_center,f.cage_min[0],f.cage_min[1],f.cage_min[2]);

    //check if the move is valid
    bool OK=true;
    for(uint eid : f.eid)
    {
      if( isBorderEdge(eid) ) continue;
      if(isValidEdge(eid)==false)
      {
        OK=false;
        break;
      }
    }

    if(OK==false) //failed, rotate back
    {
      rotate_triangle(-radian,f_center,f.vertex[0],f.vertex[1],f.vertex[2]);
      rotate_triangle(-radian,f_center,f.cage_min[0],f.cage_min[1],f.cage_min[2]);
    }//end if OK

    return OK;
  }

//list<Point2d> samples40;

  //compute vertex location by avergage its copies in triangles
  void flattened_mesh::update_vertex_position(uint vid)
  {
    //if(vid==40) samples40.clear();

    //
    vertex& v=this->vertices[vid];
    Point2d pos;
    int count=0;
    float sample_dist=10;

    list< list<Point2d> > sample_lists;
    //
    for(auto eid : v.eid)
    {
      list<Point2d> samples;
      if(isBorderEdge(eid))
      {
        cerr<<"! Error: flattened_mesh::update_vertex_position can be called for interior vertex only"<<endl;
        exit(1);
      }

      vector<Point2d> geo;
      getEdgeGeometry(eid, geo);
      edge& e = this->edges[eid];

      uint fvids[4]={triangles[e.fid[0]].vid[e.vid[0]], triangles[e.fid[0]].vid[e.vid[1]],
                     triangles[e.fid[1]].vid[e.vid[2]], triangles[e.fid[1]].vid[e.vid[3]]};
      //
      if( fvids[0]==vid || fvids[3]==vid )
      {
        Vector2d vec=geo[3]-geo[0];
        int size=vec.norm()/sample_dist;
        vec=vec/size;
        for(int i=0;i<size;i++)
        {
          Point2d p=geo[0]+vec*i;
          samples.push_back(p);
          //if(vid==40) samples40.push_back(p);
          pos[0]+=p[0];
          pos[1]+=p[1];
        }
        count+=size;
      }
      else if( fvids[2]==vid || fvids[1]==vid )
      {
        Vector2d vec=geo[1]-geo[2];
        int size=vec.norm()/sample_dist;
        vec=vec/size;
        for(int i=0;i<size;i++)
        {
          Point2d p=geo[2]+vec*i;
          samples.push_back(p);
          //if(vid==40) samples40.push_back(p);
          pos[0]+=p[0];
          pos[1]+=p[1];
        }
        count+=size;
      }
      else{
        cerr<<"! Error: flattened_mesh::update_vertex_position vid="<<vid<<" is not incident to edge "<<eid<<endl;
        exit(1);
      }

      sample_lists.push_back(samples);
    }//end for eid

    pos[0]/=count;
    pos[1]/=count;

    v.pos=pos; //center
    //compute areas
    // v.area=0;
    // for(auto& samples: sample_lists)
    // {
    //   for(auto it=samples.begin();it!=samples.end();it++)
    //   {
    //     auto nit=it; nit++;
    //     if(nit==samples.end()) break;
    //     //compute area: pos, (*it) (*nit)
    //     Vector2d v1=(*it)-pos;
    //     Vector2d v2=(*nit)-pos;
    //     v.area+=fabs(v1[0]*v2[1]-v1[1]*v2[0])/2;
    //   }//end if
    // }//end samples
  }

  //------------------------------------------------------------------
  //export the flattened_mesh to svg file
  //------------------------------------------------------------------
  void flattened_mesh::export_to_svg(const string& svg_filename)
  {

    //create a svg file
    //convert the mesh to trianglesß

    //------------------------------------------------------------------
    //find dimension
    float min_x=FLT_MAX, min_y=FLT_MAX, max_x=-FLT_MAX, max_y=-FLT_MAX;
    //for (auto pHalfedge = mesh.halfedges_begin(); pHalfedge != mesh.halfedges_end(); pHalfedge++)
    for(auto& tri: triangles)
    {
      for(int i=0;i<3;i++)
      {
        float x = tri.cage[i][0];
        float y = tri.cage[i][1];
        if(x<min_x) min_x=x;
        if(x>max_x) max_x=x;
        if(y<min_y) min_y=y;
        if(y>max_y) max_y=y;
      }
    }


    float scale=1.0f;
    int width =(int)ceil(max_x-min_x);
    int height=(int)ceil(max_y-min_y);

    cout<<"width="<<width<<" height="<<height<<endl;

    const int desired_width=10000;

    //if(width<desired_width)
    {
      scale = desired_width*1.0f/width;
      width=(int)ceil(scale*width);
      height=(int)ceil(scale*height);
      min_x*=scale;
      max_x*=scale;
      min_y*=scale;
      max_y*=scale;
    }

    int padding=(int)ceil(min(width,height)*1.0f/20);
    width+=padding;
    height+=padding;
    svg::Dimensions dimensions(width,height);
    svg::Document doc(svg_filename, svg::Layout(dimensions, svg::Layout::TopLeft));

    //------------------------------------------------------------------

    for(auto& tri: triangles)
    {
      //draw cage
      {
        auto& v1 = tri.cage[0];
        auto& v2 = tri.cage[1];
        auto& v3 = tri.cage[2];

        svg::Polygon poly(svg::Stroke(0.1, svg::Color::Black));
        poly << svg::Point(scale * v1[0]+padding-min_x,scale * v1[1]+padding-min_y);
        poly << svg::Point(scale * v2[0]+padding-min_x,scale * v2[1]+padding-min_y);
        poly << svg::Point(scale * v3[0]+padding-min_x,scale * v3[1]+padding-min_y);
        poly.endBoundary();
        doc << poly;
      }

      //draw min cage
      {
        //cout<<"is triangle="<<f->is_triangle()<<endl;
        auto& v1 = tri.cage_min[0];
        auto& v2 = tri.cage_min[1];
        auto& v3 = tri.cage_min[2];

        svg::Polygon poly(svg::Fill(svg::Color::Yellow), svg::Stroke(0.1, svg::Color::Black));
        poly << svg::Point(scale * v1[0]+padding-min_x,scale * v1[1]+padding-min_y);
        poly << svg::Point(scale * v2[0]+padding-min_x,scale * v2[1]+padding-min_y);
        poly << svg::Point(scale * v3[0]+padding-min_x,scale * v3[1]+padding-min_y);
        poly.endBoundary();
        doc << poly;
      }

      //draw triangle
      {
        //cout<<"is triangle="<<f->is_triangle()<<endl;
        auto& v1 = tri.vertex[0];
        auto& v2 = tri.vertex[1];
        auto& v3 = tri.vertex[2];

        svg::Polygon poly(svg::Fill(svg::Color::Red), svg::Stroke(0.2, svg::Color::Black));
        poly << svg::Point(scale * v1[0]+padding-min_x,scale * v1[1]+padding-min_y);
        poly << svg::Point(scale * v2[0]+padding-min_x,scale * v2[1]+padding-min_y);
        poly << svg::Point(scale * v3[0]+padding-min_x,scale * v3[1]+padding-min_y);
        poly.endBoundary();
        doc << poly;
      }
    }//end for(auto& tri: triangles)

    //draw edges
    cout<<"draw edges"<<endl;
    for(auto& edge:edges)
    {
      //cout<<"edge.fid.size()="<<edge.fid.size()<<endl;
      vector<Point2d> geo;
      getEdgeGeometry(edge, geo);

      svg::Color color=(isValidEdge(edge.eid)?svg::Color::Blue:svg::Color::Magenta);
      svg::Polygon poly(svg::Fill(color), svg::Stroke(0.2, svg::Color::Black));
      for(Point2d& p:geo) poly << svg::Point(scale * p[0]+padding-min_x,scale * p[1]+padding-min_y);
      poly.endBoundary();
      doc << poly;
    }

    //draw vertices
    for(auto& v: vertices)
    {
      const auto& p=v.pos;
      svg::Circle circ(svg::Point(scale * p[0]+padding-min_x,scale * p[1]+padding-min_y), 10, svg::Fill(svg::Color::Red), svg::Stroke(0.2, svg::Color::Black));
      doc << circ;
    }
    //draw triangle id
    char id[32]={};
    for(auto& tri: triangles)
    {
      //cout<<"is triangle="<<f->is_triangle()<<endl;
      auto o=center(tri.vertex[0],tri.vertex[1],tri.vertex[2]);
      sprintf(id,"%d",(int)tri.fid);
      svg::Text text(svg::Point(scale *o[0]+padding-min_x,scale * o[1]+padding-min_y),
                     id, svg::Fill(svg::Color::Lime),svg::Font(5));
      doc << text;
    }

    for(int iv=0;iv<vertices.size();iv++)
    {
      auto& v=vertices[iv];
      //cout<<"is triangle="<<f->is_triangle()<<endl;
      auto & o=v.pos;
      sprintf(id,"%d",(int)iv);
      svg::Text text(svg::Point(scale *o[0]+padding-min_x,scale * o[1]+padding-min_y),
                     id, svg::Fill(svg::Color::Yellow),svg::Font(5));
      doc << text;
    }
    //------------------------------------------------------------------

    doc.save();
    cout << "- Saved " << svg_filename << endl;
  }


}//end namespace flattener
