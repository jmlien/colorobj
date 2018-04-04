//------------------------------------------------------------------------------
//  Copyright 2007-2009 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _BF_MODEL_H_
#define _BF_MODEL_H_

#include <mathtool/Point.h>
#include <mathtool/Vector.h>
#include <mathtool/Matrix.h>
#include <mathtool/Quaternion.h>
using namespace mathtool;

#include <string>
#include <cassert>
#include <set>
#include <map>
#include <algorithm>
using namespace std;

#include "objReader.h"

//shorten some notations
typedef unsigned int uint;

//defined in cd.h
//struct cd_model;
namespace masc
{
//
class Mesh_sphere; //defined in spheres.h

struct object3D
{
	object3D()
	{
		for (int i = 0; i < 3; i++)
		{
			mat_color[i] = 1;
			mat_specular[i] = 1;
			mat_emission[i] = 0;
		}

		mat_shininess = 128;
	}

	virtual ~object3D()
	{

	}

	//
	string name;

	//material properties
	mathtool::Vector3d mat_color;      //diffuse color
	mathtool::Vector3d mat_specular;   //specular color
	mathtool::Vector3d mat_emission;   //emission color
	float    mat_shininess;  //how shiny the specular is

	//texture filenames
	string color_texture_filename;
	string normalmap_texture_filename;
};

inline ofstream & operator<<(ofstream & out, const object3D& obj)
{
	if (obj.name.empty() == false) out << "name=" << obj.name;
	out << " shininess=" << obj.mat_shininess << " cr=" << obj.mat_color[0] << " cg=" << obj.mat_color[1] << " cb=" << obj.mat_color[2];
	out << " sr=" << obj.mat_specular[0] << " sg=" << obj.mat_specular[1] << " sb=" << obj.mat_specular[2];
	if (obj.mat_emission.normsqr() != 0) out << " er=" << obj.mat_emission[0] << " eg=" << obj.mat_emission[1] << " eb=" << obj.mat_emission[2];
	if (obj.color_texture_filename.empty() == false) out << " text_c=" << obj.color_texture_filename;
	if (obj.normalmap_texture_filename.empty() == false) out << " text_n=" << obj.normalmap_texture_filename;

	return out;
}

//a sphere
struct sphere : public object3D
{
	sphere(){ radius = 1; refractive_index = 1; }

	mathtool::Point3d center;
	REAL radius;

	//
	mathtool::Vector3d transparency;  //how transparent this object is
	REAL refractive_index; //we assume that the refractive index of air is 1
	//diamond is 2.419
	//amber is about 1.55
	//water is 1.333 and ice is about 1.31

	mathtool::Vector3d reflectiveness; //how mirrow like this object is
};


//a triangle of the model
struct triangle
{
	triangle(){ bounding_sphere = NULL;  }

	uint v[3]; //vertex id
	uint e[3]; //edge id

	mathtool::Vector3d n; //normal
	mathtool::Vector3d t; //tangent
	mathtool::Vector3d b; //bi-tangent

	Mesh_sphere * bounding_sphere;

	//backups
	mathtool::Vector3d bk_n;

};

//a vertex of the model
struct vertex
{
	vertex(){ concave = false; visited_GID = 0; }
	mathtool::Point3d p;  //position
  mathtool::Vector3d c; //color
	//
	// per vertex TBN and UV coord
	// note: in theory, a vertex may have multiple of these values but
	//       in this implmentation only one copy per vertex. You should
	//       modify this code if this assumption is not true
	//
	mathtool::Vector3d n; //normal
	mathtool::Vector3d t; //tangent
	mathtool::Vector3d b; //binormal
	mathtool::Vector2d uv; //texture uv

	list<uint> m_f;
	list<uint> m_e; //a list of edges

	//used for search
	unsigned long visited_GID;
	uint from;

	//backups
	mathtool::Point3d bk_p;

	//if concave, set to true
	bool concave;
};

//an edge of the model
struct edge
{
	edge(){ type = 'x'; vid[0] = vid[1] = UINT_MAX; }
	uint vid[2];
	vector<uint> fid;

	mathtool::Vector3d v;       //parallel vector
	mathtool::Vector3d in_n[2]; //inface normals

	//backups
	mathtool::Vector3d bk_v;       //parallel vector
	mathtool::Vector3d bk_in_n[2]; //inface normals

	//type, c-convex, r-reflex, p-plane, b-border
	char type;
};

class model : public object3D
{
public:

	//initialization
	model()
	{
		v_size = e_size = t_size = 0;
		vertices = NULL;
		edges = NULL;
		tris = NULL;
		m_R = 0;
		current_scale = 1.0;
		//m_CD = NULL;
		cast_shadow = false;
		m_negated = false;
	}

	~model(){}

	void destroy();

	void copy(model * m)
	{
		v_size = m->v_size;
		e_size = m->e_size;
		t_size = m->t_size;

		tris = new triangle[t_size];
		edges = new edge[e_size];
		vertices = new vertex[v_size];

		assert(tris || edges || vertices);

		memcpy(tris, m->tris, sizeof(triangle)*t_size);
		memcpy(edges, m->edges, sizeof(edge)*e_size);
		memcpy(vertices, m->vertices, sizeof(vertex)*v_size);

		m_COM = m->m_COM;
		m_R = m->m_R;
		m_filename = m->m_filename;
	}

	//build this model
	bool build(const string & name);

	//build this model
	bool build(objreader::objModel<REAL>& data);

	//negate point/facets ...
	void negate();

	//reverse facets ...
	void reverse();

	//perturb/unperturb normals and vector but not vertex positions
	//void perturb(double noise);
	//void unperturb();

	//
	// set current transform matrix
	// note: these functions do not change the coordinates of the vertices/edges/faces
	//       unless the function "transform()" is called to apply the current transform
	//
	void setCurrentTransform(const mathtool::Vector3d& T, const Quaternion<REAL>& Q, const REAL S);
	const mathtool::Matrix4x4 &  getCurrentTransform() const {
		return current_transform;
	}

	const mathtool::Vector3d &  getCurrentTranslation() const {
		return current_translate;
	}

	const Quaternion<REAL> &  getCurrentOrientation() const {
		return current_orientation;
	}

	REAL getCurrentScale() const {
		return current_scale;
	}

	bool negated() const
	{
		return m_negated;
	}

	//use current transform to transform the give point and vector
	Point3d X(const Point3d& pos) const;
	Vector3d X(const Vector3d& vec) const;


	//move COM to (0,0,0), scale R to 1
	void normalize();

	const mathtool::Point3d& getCOM() const;
	const REAL getR() const;

	// filename
	const string& getFilename() const;

	//collision detection model
	void buid_collision_detection();

	//cd_model * getCD() const { return m_CD; }

	//data
	vertex   * vertices;  //vertices
	triangle * tris;      //triangles
	edge     * edges;     //edges
	uint v_size, e_size, t_size;
	bool cast_shadow;     //does this model used to cast shadow? initially this value is false

	//used only when this model is constructed by mksum
	//store a list of triangle ids that are created from the ROI feature seed in minkowski.cpp
	list<uint> roi_seed_tri_ids;

protected:

	//current transform of this model
	mathtool::Vector3d   current_translate;      //current position defined
	REAL                 current_scale;          //current scale defined
	Quaternion<REAL>     current_orientation;    //current orientation defined
	mathtool::Matrix4x4  current_transform;



	//
	//
	//
	// (private) Transformation operations
	//
	// these should be be called from outside
	// the coodinates of the vertices and normal/tangent/bitangent changes
	// the external users should use "setCurrentTransform"
	//
	// Transform
	// Rotation
	// Scale
	// Negation (public, above)
	// Reverse (turn inside out) (public, above)
	//

	//
	//note: these functions change the coodinates of vertice/faces/edges...
	//

	virtual void transform();

	virtual void restore_from_backup();

	//transform this model
	void transform(const mathtool::Vector3d& T, const mathtool::Matrix3x3& M, REAL S);

	//rotate points
	void rotate(const mathtool::Matrix2x2& m);
	void rotate(const mathtool::Matrix3x3& M);

	//scale the model
	void scale(REAL s);

	//

	//cd_model * m_CD;

	// compute COM and R
	void computeCOM_R();

	mathtool::Point3d m_COM;
	REAL m_R;

	string m_filename;

	bool m_negated;
};

//output
inline ofstream & operator<<(ofstream & out, const model& m)
{
	const mathtool::Matrix4x4 & M = m.getCurrentTransform();
	out << m.getFilename() << " T4x4 = ";
	for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) out << M[i][j] << " ";
	out << static_cast<const object3D&>(m);

	return out;
}


}//end of namespace masc

#endif //_BF_MODEL_H_
