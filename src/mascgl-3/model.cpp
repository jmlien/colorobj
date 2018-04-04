//------------------------------------------------------------------------------
//  Copyright 2007-2014 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "model.h"
#include "modelgraph/ModelGraph.h"
//#include "cd.h"

#ifdef WIN32
#pragma warning(disable:4244)
#endif

namespace masc
{
//
// build the model from a file
// the file will contain an obj model
// once obj file is read, vertices and facets will be built
// other information associated with facets and vertices,
// as well as edges will be build also using model graph
//

bool model::build(const string & name)
{
	this->m_filename = name;
	//if(pts!=NULL) return false; //built

	//build/copy model
	objreader::objReader<REAL> reader(name);
	if (!reader.Read())
		return false;
	auto& data = reader.getModel();

	if (data.textcoords.empty())
		cerr << "! Warning: Model does not have texture coordinates. Texturing will not work." << endl;

	if(data.colors.empty())
	{
		cerr << "! Warning: Model does not have vertex color. Output texture will be black." << endl;
	}

	return build(data);
}

bool model::build(objreader::objModel<REAL>& data)
{
	using namespace objreader;
	{
		//allocate memory
		v_size = (uint)data.pts.size();
		t_size = (uint)data.polys.size();
		vertices = new vertex[v_size];   //
		tris = new triangle[t_size];     //
		assert(vertices&&tris);        //make sure enough memory

		//copy vertices
		for (uint i = 0; i < v_size; i++){
			vertices[i].p.set(&data.pts[i].x);
			vertices[i].bk_p = vertices[i].p;  //backup for modification
			if(!data.colors.empty())
			{
				vertices[i].c.set(&data.colors[i].x);
			}
		}

		//copy triangles
		int tid = 0;
		for (list<polygon>::iterator i = data.polys.begin(); i != data.polys.end(); i++)
		{
			list<int>& ids = i->pts;
			//check if triangle
			if (ids.size() != 3){
				cerr << "! Error: polygon " << tid << " is not a triangle." << endl;
				return false;
			}

			int vid = 0;
			list<int>::iterator t = i->textures.begin();
			list<int>::iterator n = i->normals.begin();
			for (list<int>::iterator j = ids.begin(); j != ids.end(); j++, vid++)
			{
				tris[tid].v[vid] = *j;
				vertices[*j].m_f.push_back(tid);

				//handle vertex texture coordinates
				if (t != i->textures.end())
				{
					if (*t < data.textcoords.size())
					{
						Vector2d uv = Vector2d(&data.textcoords[*t].x); //texture coordinates
						vertices[*j].uv = uv;
					}
					t++;
				}
				//done with vertex texture coordinates

				//handle vertex normal
				Vector3d vn = Vector3d(&data.normals[*n].x); //vertex normal
				vertices[*j].n = vn;
				n++;
				//done with vertex normal
			}

			tid++; //triangle id
		}

    }//end build/copy model


    {//build model graph and collect informations
        CModelGraph<REAL> G;
        G.doInit(this);

        //create edges
        e_size=G.getEdgeSize();
        auto * ptrE=G.getEdges();
        edges=new edge[e_size];
        assert(edges);
        for(uint i=0;i<e_size;i++,ptrE=ptrE->getNext()){
            int v1=edges[i].vid[0]=ptrE->getStartPt();
            int v2=edges[i].vid[1]=ptrE->getEndPt();

            const vector<int>&tmp_fid=ptrE->getFacets();
            edges[i].fid.insert(edges[i].fid.end(),tmp_fid.begin(),tmp_fid.end());

			if( tmp_fid.size()<2 ){ //check if boundary edge
				edges[i].fid.push_back(tmp_fid.back());
				edges[i].type='b'; //bd
			}

            //compute parallel vector
            edges[i].v=ptrE->getV();

			//compute edge normals
			edges[i].in_n[0]=ptrE->getInNormal(0);

			if (edges[i].type == 'b')
				edges[i].in_n[1] = edges[i].in_n[0];
			else
				edges[i].in_n[1] = ptrE->getInNormal(1);

			//backup edge normals
			edges[i].bk_v = edges[i].v;
			edges[i].bk_in_n[0] = edges[i].in_n[0];
			edges[i].bk_in_n[1] = edges[i].in_n[1];

            vertices[v1].m_e.push_back(i);
            vertices[v2].m_e.push_back(i);
        }//end i

        //facets
        auto& facets=G.getFacets();
        for(uint i=0;i<t_size;i++)
		{
			tris[i].bk_n=tris[i].n = facets[i].n;
            for(int j=0;j<3;j++){
                tris[i].e[j]=facets[i].m_Edge[j]->getID();
            }//end j
        }//end i

        //edge type
        ptrE=G.getEdges();
        for(uint i=0;i<e_size;i++,ptrE=ptrE->getNext())
		{
			if(edges[i].type=='b') continue; //already know
            Vector3d& n1=tris[edges[i].fid[0]].n;
            Vector3d& n2=tris[edges[i].fid[1]].n;
			//
			REAL d=n1*n2;
            if(fabs(1-d)<SMALLNUMBER)
				edges[i].type='p'; //plane
            else{
                Vector3d vec=(n1%n2).normalize();
                if(vec*edges[i].v>0) edges[i].type='c'; //convex
                else edges[i].type='r'; //reflex
            }
        }

		//vertex type
		typedef list<uint>::iterator IT;
		for(uint i=0;i<v_size;i++)
		{
			vertex& v=vertices[i];

			Vector3d evec;
			for(IT ie=v.m_e.begin();ie!=v.m_e.end();ie++){
				edge& e=edges[*ie];
				Vector3d vec=e.v;
				if(e.vid[1]==i) vec=-vec;
				evec=evec+vec;
			}//end ir

			Vector3d fvec;
			for(IT f=v.m_f.begin();f!=v.m_f.end();f++){
				triangle& t=tris[*f];
				fvec=fvec+t.n;
			}//end ir

			if(evec*fvec>0) v.concave=true;
		}//end i
    }

    return true;
}


void model::buid_collision_detection()
{
	assert(false);
}


void model::destroy()
{
	delete[] tris;     tris = NULL;
	delete[] edges;    edges = NULL;
	delete[] vertices; vertices = NULL;
	v_size = e_size = t_size = 0;
	//might need to relase CD info here
	cast_shadow = false;
	//delete m_CD;
}

void model::setCurrentTransform(const Vector3d& T, const Quaternion<REAL>& Q, const REAL S)
{
	Matrix3x3 R = Q.getMatrix();
    Matrix4x4 m(R[0][0], R[0][1],   R[0][2],   T[0],
                R[1][0], R[1][1],   R[1][2],   T[1],
                R[2][0], R[2][1],   R[2][2],   T[2],
                0,       0,         0,         1);

    Matrix4x4 s(S, 0, 0, 0,
                0, S, 0, 0,
                0, 0, S, 0,
                0, 0, 0, 1);

	current_translate=T;            //current position defined
	current_scale=S;          //current scale defined
	current_orientation=Q;    //current orientation defined
	current_transform=m*s;
}

//use current transform to transform the give point
Point3d model::X(const Point3d& pos) const
{
	Vector<REAL, 4> tmp(pos[0], pos[1], pos[2], 1);
	tmp = current_transform*tmp;
	Vector3d result(tmp[0] / tmp[3], tmp[1] / tmp[3], tmp[2] / tmp[3]);
	if (m_negated) result = -result;
	return result;
}

//use current transform to transform the give vector
Vector3d model::X(const Vector3d& vec) const
{
	Vector3d v=current_orientation.rotate(vec);
	if (m_negated) v = -v;
	return v;
}

const Point3d& model::getCOM() const
{
    return this->m_COM;
}

const REAL model::getR() const
{
    return this->m_R;
}

void model::computeCOM_R()
{
    Vector3d com;
    for(uint i=0;i<this->v_size;++i)
        com = com + Vector3d(this->vertices[i].p.get());

    com = com * (1.0 / this->v_size);

    REAL r=0;
    for(uint i=0;i<this->v_size;++i)
        r = std::max(r, (this->vertices[i].p - com).normsqr());
    r = sqrt(r);

    this->m_COM.set(com.get());
    this->m_R = r;
}

void model::normalize()
{
    this->computeCOM_R();

    Vector3d T = -Vector3d(this->getCOM().get());
    REAL S = 1.0f / this->getR();

    // translate first
    this->transform(T, Matrix3x3(), 1);

    // override back position
    for(uint i=0;i<v_size;++i)
        vertices[i].bk_p = vertices[i].p;

    // scale
    this->transform(Vector3d(), Matrix3x3(), S);

    // override back position
    for(uint i=0;i<v_size;++i)
        vertices[i].bk_p = vertices[i].p;

    this->computeCOM_R();
}

const string& model::getFilename() const
{
    return m_filename;
}


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
// Negation (public)
// Reverse (turn inside out) (public)
//

void model::transform()
{
	transform(getCurrentTranslation(), getCurrentOrientation().getMatrix(), getCurrentScale());
}

void model::restore_from_backup()
{
	//restore vertices
	for (uint i = 0; i<v_size; i++)
	{
		vertices[i].p = vertices[i].bk_p;
	}

	//restore edges
	for (uint i = 0; i<e_size; i++){
		edges[i].in_n[0] = edges[i].bk_in_n[0];
		edges[i].in_n[1] = edges[i].bk_in_n[1];
		edges[i].v = edges[i].bk_v;
	}

	//restore facets
	for (uint i = 0; i < t_size; i++)
	{
		tris[i].n = tris[i].bk_n;
	}
}


// Important Note: The following functions do not correctly handle TBN of the verteices and faces...
// TODO: fix this problem...
void model::transform(const Vector3d& T, const Matrix3x3& M, REAL S)
{
	Vector3d tmp;

	//move vertices
	for (uint i = 0; i<v_size; i++)
	{
		tmp.set(vertices[i].bk_p.get());
		vertices[i].p = M*(tmp*S) + T;
		//vertices[i].n = (M*vertices[i].n).normalize(); //we should have a back up value for this.
		//vertices[i].t = (M*vertices[i].t).normalize(); //we should have a back up value for this.
		//vertices[i].b = (M*vertices[i].b).normalize(); //we should have a back up value for this.
	}

	//move edges
	for (uint i = 0; i<e_size; i++){
		edges[i].in_n[0] = (M*edges[i].bk_in_n[0]).normalize();
		edges[i].in_n[1] = (M*edges[i].bk_in_n[1]).normalize();
		edges[i].v = (M*edges[i].bk_v).normalize();
	}

	//move facets
	for (uint i = 0; i < t_size; i++)
	{
		tris[i].n = (M*tris[i].bk_n).normalize();
	}
}

//
//
// Important Note: The following functions do not correctly handle TBN of the verteices and faces...
// TODO: fix this problem...
//
//
void model::rotate(const Matrix2x2& M)
{
	Matrix3x3 M3(M[0][0], M[0][1], 0, M[1][0], M[1][1], 0, 0, 0, 1);
	Vector2d tmp;
	rotate(M3);
}


void model::rotate(const Matrix3x3& M)
{
	Vector3d tmp;

	//rotate vertices
	for (uint i = 0; i<v_size; i++){
		tmp.set(vertices[i].bk_p.get());
		vertices[i].p = M*tmp;
	}

	//rotate edges
	for (uint i = 0; i<e_size; i++){
		edges[i].in_n[0] = (M*edges[i].bk_in_n[0]).normalize();
		edges[i].in_n[1] = (M*edges[i].bk_in_n[1]).normalize();
		edges[i].v = (M*edges[i].v).normalize();
	}
	//rotate facets
	for (uint i = 0; i<t_size; i++)
		tris[i].n = (M*tris[i].bk_n).normalize();
}


void model::scale(REAL s)
{
	//scale vertices
	Point3d tmp;
	for (uint i = 0; i<v_size; i++){
		tmp.set(vertices[i].bk_p.get());
		vertices[i].p.set(tmp[0] * s, tmp[1] * s, tmp[2] * s);
	}
}


void model::negate()
{

	//for (uint i = 0; i<v_size; i++){
	//	Point3d& pt = vertices[i].p;
	//	pt.set(-pt[0], -pt[1], -pt[2]);
	//}

	//for (uint i = 0; i<t_size; i++){
	//	tris[i].n = -tris[i].n;
	//	swap(tris[i].v[1], tris[i].v[2]);
	//}

	//for (uint i = 0; i<e_size; i++){
	//	edge& e = edges[i];
	//	e.v = -e.v;
	//	e.in_n[0] = -e.in_n[0];
	//	e.in_n[1] = -e.in_n[1];
	//}

	m_negated = !m_negated;
}


void model::reverse()
{
	for (uint i = 0; i < v_size; i++)
	{
		vertices[i].concave = !vertices[i].concave;
		vertices[i].n = -vertices[i].n;
	}

	for (uint i = 0; i<e_size; i++)
	{
		if		(edges[i].type == 'r') edges[i].type = 'c';
		else if (edges[i].type == 'c') edges[i].type = 'r';

		edges[i].v = -edges[i].v;
		swap(edges[i].vid[0], edges[i].vid[1]);
		swap(edges[i].fid[0], edges[i].fid[1]); //assume that it only has two faces
	}

	for (uint i = 0; i<t_size; i++)
	{
		tris[i].n = -tris[i].n;
		swap(tris[i].v[1], tris[i].v[2]);
		swap(tris[i].e[0], tris[i].e[1]);
	}
}


}//end of namespace masc

//
//void model::perturb(double noise)
//{
//
//	//build the matrix
//	Vector3d r(noise*drand48(), noise*drand48(), noise*drand48());
//	for (int i = 0; i<3; i++){
//		if (drand48()>0.5) r[i] = -r[i];
//	}
//
//	Quaternion<REAL> q(r.get());
//	Matrix3x3 M = q.getMatrix();
//
//	//rotate edges
//	for (uint i = 0; i<e_size; i++){
//		edges[i].in_n[0] = (M*edges[i].in_n[0]).normalize();
//		edges[i].in_n[1] = (M*edges[i].in_n[1]).normalize();
//		edges[i].v = (M*edges[i].v).normalize();
//	}
//	//rotate facets
//	for (uint i = 0; i < t_size; i++)
//	{
//		tris[i].n = (M*tris[i].n).normalize();
//		tris[i].t = (M*tris[i].t).normalize();
//		tris[i].b = (M*tris[i].b).normalize();
//	}
//}
//
//void model::unperturb()
//{
//	//build the matrx
//	Matrix3x3 M = this->current_orientation.getMatrix();
//
//	//rotate edges
//	for (uint i = 0; i<e_size; i++){
//		edges[i].in_n[0] = (M*edges[i].bk_in_n[0]).normalize();
//		edges[i].in_n[1] = (M*edges[i].bk_in_n[1]).normalize();
//		edges[i].v = (M*edges[i].bk_v).normalize();
//	}
//	//rotate facets
//	for (uint i = 0; i<t_size; i++)
//		tris[i].n = (M*tris[i].bk_n).normalize();
//}
