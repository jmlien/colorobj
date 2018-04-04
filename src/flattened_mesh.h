//------------------------------------------------------------------------------
//  Copyright 2010-2017 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#pragma once

#include "model.h"
#include <cfloat>

namespace flattener
{


  struct triangle
  {
    triangle()
    {
      fid=UINT_MAX;
      level=UINT_MAX;
    }

    Point2d vertex[3];
    Point2d cage[3];
    Point2d cage_min[3]; //local smallest cage similar to cage that can enclose this triangle

    //find the smallest enclosing triangle similar to cage
    //return scaling factor from cage to cage_min
    float find_smallest_cage();

    //scale cage by factor s
    void scale_cage(float s);

    //fit this triangle to cage
    //assumption: this triangle can already fit to cage_min
    //the goal is to find transform between cage_min and cage
    //and move this triangle inside cage
    void move_to_cage();

    vector<uint> eid;
    vector<uint> vid;
    uint fid; //triangle id to the original triangle in the mesh
    uint level;
  };

  ostream & operator << (ostream & out, const triangle& tri);

  struct edge
  {
    edge(){ eid=UINT_MAX; }
    uint otherf(uint id)
    {
      for(uint f:fid) if(f!=id) return f;
      return UINT_MAX;
    }

    uint eid;
    vector<uint> vid; //vertex ids in adjacent triangle (2 or 4 elements)
    vector<uint> fid; //one or two incident triangles
  };

  struct vertex
  {
    vertex(){is_interior=true;area=0;}
    bool is_interior;

    vector<uint> eid;

    Point2d pos;
    float area;
  };

  ostream & operator << (ostream & out, const edge& e);

  class flattened_mesh
  {
  public:
    flattened_mesh(masc::model& m):pm(&m)
    {
      triangles=vector<triangle>(m.t_size,triangle());
      edges=vector<edge>(m.e_size,edge());
      vertices=vector<vertex>(m.v_size,vertex());

      //initialize flattened_mesh
      for(uint fid=0;fid<pm->t_size;fid++)
      {
        triangles[fid].fid=fid;
        triangles[fid].eid=vector<uint>(pm->tris[fid].e, pm->tris[fid].e+3);
        triangles[fid].vid=vector<uint>(pm->tris[fid].v, pm->tris[fid].v+3);
        project( pm, pm->tris[fid], triangles[fid]);
      }
    }

    void init()
    {
        build_edges();
        build_vertices();
    }

    //place triangles in cage
    void place_in_cage();

    void lapaician_smoothing(int max_iteration, float step_size);

    //iterative condensation to incident vertices
    void condense_to_vertices(int max_iteration, float step_size);

    void condense_to_root_with_edge_collapsing(int max_iteration, float step_size);

    void condense_by_edge_collapsing(int max_iteration, float step_size);

    //iterative condensation to a root face
    void condense_to_root(int max_iteration, float step_size);

    //reduce edge length level by level
    void condense_level();


    void export_to_svg(const string& filename);

    //access functions
    triangle& get_triangle(uint id) { return triangles[id]; }
    vertex& get_vertex(uint id) { return vertices[id]; }

  protected:

    bool condense_to_vertices(uint fid, const float step_size);

    //compute vertex location by avergage its copies in triangles
    void update_vertex_position(uint vid);

    //moving the face to the center of neighbors
    bool lapaician_smoothing(uint fid, const float step_size);

    //condesne by moving the face to each neighboring triangle
    bool condense_by_edge_collapsing(uint fid, const float step_size);

    //build edges
    void build_edges();

    //build vertices
    //Note: this must be called after build_edges
    void build_vertices();

    //condense a given face
    //return true if successful
    bool condense_to_face(uint fid, uint rootid, const float step_size);

    //condense the given level
    //return true if successful
    bool condense_level(list<uint>& levelset);

    uint findRoot(); //root is a face with min edge length and 3 neighbors

    //order all facets into levels, return max level
    uint buildLevelSets(uint rootid, map<uint,list<uint> >& levelsets);

    //check if a given edge is valid
    bool isValidEdge(uint eid) const;
    bool isBorderEdge(uint eid) const { return edges[eid].fid.size()<2;}

    void getEdgeGeometry(uint eid, vector<Point2d>& geo) const
    {
      getEdgeGeometry(edges[eid], geo);
    }

    void getEdgeGeometry(const edge& e, vector<Point2d>& geo) const
    {
      for(int i=0;i<e.fid.size();i++)
      {
        uint fid=e.fid[i];
        uint vid1=e.vid[i*2];
        uint vid2=e.vid[i*2+1];

        auto& v1 = triangles[fid].cage_min[vid1];
        auto& v2 = triangles[fid].cage_min[vid2];
        geo.push_back(v1);
        geo.push_back(v2);
      }
    }

    float getEdgeLength(uint eid) const
    {
      vector<Point2d> geo;
      getEdgeGeometry(edges[eid], geo);
      if(geo.size()<4) return FLT_MAX/2; //not a regular edge module
      return (geo[1]-geo[2]).norm()+(geo[0]-geo[3]).norm();
    }

    void project(masc::model * m, const masc::triangle& t3d, triangle& t2d )
    {
      const Point3d& p1=m->vertices[t3d.v[0]].p;
      const Point3d& p2=m->vertices[t3d.v[1]].p;
      const Point3d& p3=m->vertices[t3d.v[2]].p;

      Vector3d u=(p2-p1).normalize();
      Vector3d v=t3d.n%u;
      Vector3d w=(p3-p1);
      t2d.vertex[0].set(0,0);
      t2d.vertex[1].set((p2-p1)*u, 0);
      t2d.vertex[2].set(w*u, w*v);
    }

    //collect edges of given fid that connect to the lower level faces
    void getLowLeveLNeighbors(uint fid,list<uint>& low_level_eids)
    {
      auto& t=triangles[fid];
      for(int i=0;i<3;i++)
      {
        uint ofid=edges[t.eid[i]].otherf(fid);
        if(ofid==UINT_MAX) continue;
        auto& ot=triangles[ofid];

        if(ot.level>=t.level) continue;
        low_level_eids.push_back(t.eid[i]);
      }
    }

    //collect faces and edges incident to this given edge
    void incidentObjects(uint eid, list<uint>& fids, list<uint>& eids) const
    {
      masc::edge & e=pm->edges[eid];
      set<uint> fset;
      for(short d=0;d<2;d++)
      {
        auto & v=pm->vertices[e.vid[d]];
        for(auto& veid : v.m_e)
        {
           if(veid!=eid) eids.push_back(veid);
           fset.insert(pm->edges[veid].fid.begin(),pm->edges[veid].fid.end());
        }//end e
      }//end d

      fids.insert(fids.end(),fset.begin(),fset.end());
    }

    //rotate the given face so that its edge modules are perpendicular to the edges
    bool update_face_orientation(uint fid, float step_size);

    //attempts to translate and rotate
    //if failed, the triangle will be moved/rotated back and return false
    //if suceessed, return true and the triangle will be moved/rotated
    bool attempt_translate(const Vector2d& T, triangle& t);
    bool attempt_rotate(float radian, triangle& t);

  private:
    masc::model * pm;
    vector<vertex> vertices;
    vector<triangle> triangles;
    vector<edge> edges;
  };

  ostream & operator << (ostream & out, const flattened_mesh& fm);

};//end name space flattener
