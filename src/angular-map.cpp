#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/IO/Polyhedron_iostream.h>
//#include <CGAL/Surface_mesh_parameterization/Parameterization_SurfaceMesh_3.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>
//#include <CGAL/Surface_mesh_parameterization/Parameterization_mesh_patch_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Barycentric_mapping_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Mean_value_coordinates_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/LSCM_parameterizer_3.h>
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <boost/foreach.hpp>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "angular-map.h"
#include "simple_svg_1.0.0.hpp"

namespace flattener
{

  struct editable_mesh
  {
    struct V
    {
      uint old_vid;
      uint vid;
    };

    struct E{
      uint v[2];
      uint f[2];
    };

    struct F
    {
      uint v[3];
      uint e[3];
      uint nei[3];
    };

    uint add_v(uint vid)
    {
      V v;
      v.old_vid=vid;
      v.vid=last_vid++;
      vs_[v.vid]=v;
      return v.vid;
    }

    uint add_e(uint fid)
    {
      //V v;
      //v.old_vid=vid;
      //v.vid=last_vid++;
      //vs_[v.vid]=v;
      return 0;//v.vid;
    }

    void add_f(masc::triangle& t)
    {
      F f;
      for(int i=0;i<3;i++) f.v[i]=add_v(t.v[i]);
      //for(int i=0;i<3;i++) f.e[i]=add_e(f,i);
    }

    vector<F>   fs_;
    map<uint,V> vs_;
    map<uint,E> es_;
    uint last_vid;
  };

  bool ANGULAR_MAP_CREATOR::cut(masc::model * m)
  {
    return cut(m, 0);
  }

  bool ANGULAR_MAP_CREATOR::cut(masc::model * m, uint seed)
  {
    //set<uint> vids;
    set<uint> fids;
    vector< vector<uint> > neighbors(m->t_size, vector<uint>(3,UINT_MAX)); //one for each triangle

    //construct an mesh by iteratively add facets of m into it
    //starting from the seed facet
    list<uint> open;
    open.push_back(seed);
    fids.insert(seed);
    //for(short i=0;i<3;i++) vids.insert(m->tris[seed].v[i]);
    uint vertex_count=3;
    uint edge_count=3;

    while(open.empty()==false)
    {
      uint fid=open.front();
      open.pop_front();
      auto& tri=m->tris[fid];

      //check neighbors
      // for(short i=0;i<3;i++)
      // {
      //   //get the other triangle
      //   auto& e=tri.edges[i];
      //   if(e.fid.size()>2){ std::cerr<<"! Error: the input mesh is non-manifold"<<std::endl; return false; }
      //   uint otherf = e.other(fid);
      //   if(otherf==UINT_MAX) continue; //no adj triangle via this edge (i.e. a hole edge)
      //   //options that can be apply
      //   if(fids.find(otherf)==fids.end()) //this facet is new
      //   {
      //     //check if we can add the face and edge edge
      //     edge_count+=2;
      //     vertex_count++;
      //   }
      //   else //this face has been added
      //   {
      //     if(neighbors[fid][i]!=otherf) //not a neighbor yet....
      //     {
      //       //check if we can add this face as the neighbor
      //       edge_count--;
      //       vertex_count--; //this is not the only case...
      //     }//end if not a neighbor yet
      //   }//end if this facet is new
      // }
    }

    return true;
  }


  // ----------------------------------------------------------------------------
  // Private types
  // ----------------------------------------------------------------------------
  typedef CGAL::Simple_cartesian<double>       Kernel;
  typedef Kernel::Point_3                      Point_3;
  typedef Kernel::Point_2                      Point_2;
  typedef CGAL::Surface_mesh<Kernel::Point_3>  SurfaceMesh;

  typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor    vertex_descriptor;
  typedef boost::graph_traits<SurfaceMesh>::face_descriptor      face_descriptor;

  typedef SurfaceMesh::Property_map<vertex_descriptor, Point_2>   UV_pmap;
  namespace SMP = CGAL::Surface_mesh_parameterization;
  //typedef CGAL::Parameterization_SurfaceMesh_3<Polyhedron> SurfaceMesh;


  //converting masc model to CGAL mesh
  void m2m( masc::model * m_, SurfaceMesh& m)
  {
    assert(m_); //make sure m is at least not null

    //init
    typedef SurfaceMesh::Vertex_index VID;
    typedef SurfaceMesh::Face_index   FID;
    vector<VID> vids;
    vids.reserve(m_->v_size);

    // add the polyhedron vertices
    for( uint i=0; i<m_->v_size; i++ ){
      auto& pos=m_->vertices[i].p;
      VID vid=m.add_vertex( Point_3( pos[0], pos[1], pos[2] ) );
      vids.push_back(vid);
    }

    // add the polyhedron triangles
    for( uint i=0; i<m_->t_size; i++ ){
     auto& tri = m_->tris[i];
     m.add_face( vids[tri.v[0]], vids[tri.v[1]], vids[tri.v[2]] );
    }
  }

  //------------------------------------------------------------------
  //export the parameterized mesh to svg file
  //------------------------------------------------------------------
  void export_to_svg(SurfaceMesh & mesh, UV_pmap uv_map)
  {

    //create a svg file
    //convert the mesh to trianglesß
    char svg_filename[256];
    sprintf(svg_filename, "%s.svg", "cgal-param");

    //------------------------------------------------------------------
    //find dimension
    float min_x=FLT_MAX, min_y=FLT_MAX, max_x=-FLT_MAX, max_y=-FLT_MAX;
    //for (auto pHalfedge = mesh.vertex(); pHalfedge != mesh.halfedges_end(); pHalfedge++)
    //BOOST_FOREACH( vertex_descriptor vd, mesh.vertices())
    for(auto vd : mesh.vertices())
    {
        Point_2 uv=uv_map[vd];
        float x =  uv.x();
        float y =  uv.y();
        if(x<min_x) min_x=x;
        if(x>max_x) max_x=x;
        if(y<min_y) min_y=y;
        if(y>max_y) max_y=y;
    }
    float scale=1.0f;
    int width =(int)ceil(max_x-min_x);
    int height=(int)ceil(max_y-min_y);
    const int desired_width=800;
    if(width<desired_width)
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
    //draw the external boundary
    //for (auto f = mesh.facets_begin(); f != mesh.facets_end(); f++)
    for(auto fd : mesh.faces())
    {
        vector<Point_2> uvs;
        for(auto vd : vertices_around_face(mesh.halfedge(fd), mesh))
        {
          uvs.push_back(uv_map[vd]);
        }

        assert(uvs.size()==3);

        svg::Polygon poly(svg::Fill(svg::Color::Yellow), svg::Stroke(0.5, svg::Color::Black));
        poly << svg::Point(scale * uvs[0].x()+padding-min_x,scale * uvs[0].y()+padding-min_y);
        poly << svg::Point(scale * uvs[1].x()+padding-min_x,scale * uvs[1].y()+padding-min_y);
        poly << svg::Point(scale * uvs[2].x()+padding-min_x,scale * uvs[2].y()+padding-min_y);
        poly.endBoundary();
        doc << poly;
    }

    //------------------------------------------------------------------


    doc.save();
    cout << "- Saved " << svg_filename << endl;
  }

  bool CGAL_AMAP_Creator::create(masc::model * m, flattened_mesh& fm)
  {

    //create mesh
    SurfaceMesh mesh;

    // build a CGAL surface mesh from model *m
    m2m(m,mesh);

    //***************************************
    // Parameterization
    //***************************************
    UV_pmap uv_map = mesh.add_property_map<vertex_descriptor, Point_2>("h:uv").first;

    // Border parameterizer
    typedef SMP::Square_border_arc_length_parameterizer_3<SurfaceMesh> s_arc_bd;
    typedef SMP::Square_border_uniform_parameterizer_3<SurfaceMesh> s_uni_bd;
    typedef SMP::Circular_border_uniform_parameterizer_3<SurfaceMesh> c_uni_bd;
    typedef SMP::Circular_border_arc_length_parameterizer_3<SurfaceMesh> c_arc_bd;

    // Discrete Authalic Parameterization (square border)
    // with Eigen solver
    typedef SMP::Discrete_authalic_parameterizer_3<SurfaceMesh, s_arc_bd> authalic_s_arc_bd;
    typedef SMP::Discrete_authalic_parameterizer_3<SurfaceMesh, s_uni_bd> authalic_s_uni_bd;
    typedef SMP::Discrete_authalic_parameterizer_3<SurfaceMesh, c_arc_bd> authalic_c_arc_bd;
    typedef SMP::Discrete_authalic_parameterizer_3<SurfaceMesh, c_uni_bd> authalic_c_uni_bd;

    typedef SMP::Discrete_conformal_map_parameterizer_3<SurfaceMesh, s_arc_bd> conformal_s_arc_bd;
    typedef SMP::Discrete_conformal_map_parameterizer_3<SurfaceMesh, s_uni_bd> conformal_s_uni_bd;
    typedef SMP::Discrete_conformal_map_parameterizer_3<SurfaceMesh, c_arc_bd> conformal_c_arc_bd;
    typedef SMP::Discrete_conformal_map_parameterizer_3<SurfaceMesh, c_uni_bd> conformal_c_uni_bd;

    typedef SMP::Barycentric_mapping_parameterizer_3<SurfaceMesh, s_arc_bd> barycentric_s_arc_bd;
    typedef SMP::Barycentric_mapping_parameterizer_3<SurfaceMesh, s_uni_bd> barycentric_s_uni_bd;
    typedef SMP::Barycentric_mapping_parameterizer_3<SurfaceMesh, c_arc_bd> barycentric_c_arc_bd;
    typedef SMP::Barycentric_mapping_parameterizer_3<SurfaceMesh, c_uni_bd> barycentric_c_uni_bd;

    typedef SMP::Mean_value_coordinates_parameterizer_3<SurfaceMesh, s_arc_bd> mv_s_arc_bd;
    typedef SMP::Mean_value_coordinates_parameterizer_3<SurfaceMesh, s_uni_bd> mv_s_uni_bd;
    typedef SMP::Mean_value_coordinates_parameterizer_3<SurfaceMesh, c_arc_bd> mv_c_arc_bd;
    typedef SMP::Mean_value_coordinates_parameterizer_3<SurfaceMesh, c_uni_bd> mv_c_uni_bd;

    typedef SMP::LSCM_parameterizer_3<SurfaceMesh, s_arc_bd> lscm_s_arc_bd;
    typedef SMP::LSCM_parameterizer_3<SurfaceMesh, s_uni_bd> lscm_s_uni_bd;
    typedef SMP::LSCM_parameterizer_3<SurfaceMesh, c_arc_bd> lscm_c_arc_bd;
    typedef SMP::LSCM_parameterizer_3<SurfaceMesh, c_uni_bd> lscm_c_uni_bd;

    SMP::Error_code err = SMP::ERROR_EMPTY_MESH;

    // A halfedge on the border
    halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(mesh).first;

    if(param_method_==Tutte)
    {
      switch(bd_type_)
      {
        case Square_Arc_Length:   err = SMP::parameterize(mesh, barycentric_s_arc_bd(), bhd, uv_map); break;
        case Square_Uniform:      err = SMP::parameterize(mesh, barycentric_s_uni_bd(), bhd, uv_map); break;
        case Circular_Arc_Length: err = SMP::parameterize(mesh, barycentric_c_arc_bd(), bhd, uv_map); break;
        case Circular_Uniform:    err = SMP::parameterize(mesh, barycentric_c_uni_bd(), bhd, uv_map); break;
        default: break;
      }
    }
    else if(param_method_==Conformal_Map)
    {
      switch(bd_type_)
      {
        case Square_Arc_Length:   err = SMP::parameterize(mesh, conformal_s_arc_bd(), bhd, uv_map); break;
        case Square_Uniform:      err = SMP::parameterize(mesh, conformal_s_uni_bd(), bhd, uv_map); break;
        case Circular_Arc_Length: err = SMP::parameterize(mesh, conformal_c_arc_bd(), bhd, uv_map); break;
        case Circular_Uniform:    err = SMP::parameterize(mesh, conformal_c_uni_bd(), bhd, uv_map); break;
        default: break;
      }
    }
    else if(param_method_==Floater)
    {
      switch(bd_type_)
      {
        case Square_Arc_Length:   err = SMP::parameterize(mesh, mv_s_arc_bd(), bhd, uv_map); break;
        case Square_Uniform:      err = SMP::parameterize(mesh, mv_s_uni_bd(), bhd, uv_map); break;
        case Circular_Arc_Length: err = SMP::parameterize(mesh, mv_c_arc_bd(), bhd, uv_map); break;
        case Circular_Uniform:    err = SMP::parameterize(mesh, mv_c_uni_bd(), bhd, uv_map); break;
        default: break;
      }
    }
    else if(param_method_==Authalic)
    {
      switch(bd_type_)
      {
        case Square_Arc_Length:   err = SMP::parameterize(mesh, authalic_s_arc_bd(), bhd, uv_map); break;
        case Square_Uniform:      err = SMP::parameterize(mesh, authalic_s_uni_bd(), bhd, uv_map); break;
        case Circular_Arc_Length: err = SMP::parameterize(mesh, authalic_c_arc_bd(), bhd, uv_map); break;
        case Circular_Uniform:    err = SMP::parameterize(mesh, authalic_c_uni_bd(), bhd, uv_map); break;
        default: break;
      }
    }
    else if(param_method_==LeastSquares)
    {
      switch(bd_type_)
      {
        case Square_Arc_Length:   err = SMP::parameterize(mesh, lscm_s_arc_bd(), bhd, uv_map); break;
        case Square_Uniform:      err = SMP::parameterize(mesh, lscm_s_uni_bd(), bhd, uv_map); break;
        case Circular_Arc_Length: err = SMP::parameterize(mesh, lscm_c_arc_bd(), bhd, uv_map); break;
        case Circular_Uniform:    err = SMP::parameterize(mesh, lscm_c_uni_bd(), bhd, uv_map); break;
        default: break;
      }
    }

    switch(err) {
      case SMP::OK: // Success
          break;
      case SMP::ERROR_EMPTY_MESH: // Input mesh not supported
      case SMP::ERROR_NON_TRIANGULAR_MESH:
      case SMP::ERROR_NO_TOPOLOGICAL_DISC:
      case SMP::ERROR_BORDER_TOO_SHORT:
          std::cerr << "! Error: Onput mesh not supported: " << SMP::get_error_message(err) << std::endl;
          return false;
          break;
      default: // Error
          std::cerr << "! Error: " << SMP::get_error_message(err) << std::endl;
          return false;
          break;
    };

    //compute the mesh angles in parameterized domain
    uint fid=0;
    // for (auto f = mesh.facets_begin(); f != mesh.facets_end(); f++, fid++)
    // {
    for(auto fd : mesh.faces())
    {
        vector<Point_2> uvs;
        for(auto vd : vertices_around_face(mesh.halfedge(fd), mesh))
        {
          uvs.push_back(uv_map[vd]);
        }

        assert(uvs.size()==3);

        //cout<<"is triangle="<<f->is_triangle()<<endl;
        auto& ff = fm.get_triangle(fid);
        ff.cage[0].set(uvs[0].x(),uvs[0].y());
        ff.cage[1].set(uvs[1].x(),uvs[1].y());
        ff.cage[2].set(uvs[2].x(),uvs[2].y());
        fid++;
    }

    uint vid=0;
    for(auto vd : mesh.vertices())
    {
      auto& vv = fm.get_vertex(vid);
      Point_2 uv=uv_map[vd];
      vv.pos.set(uv.x(), uv.y());
      vid++;
    }

    //export to svg
    export_to_svg(mesh, uv_map);

    return true;
  }

  bool LP_AMAP_Creator::create(masc::model * m, flattened_mesh& fm)
  {
    return true;
  }

}//end namespace flattener
