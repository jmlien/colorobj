//------------------------------------------------------------------------------
//  Copyright 2007-2014 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------


#ifndef _MODEL_GRAPH_H_
#define _MODEL_GRAPH_H_

#include "ModelNode.h"
#include "ModelEdge.h"
#include "model.h"

namespace masc
{

template<typename T>
struct CModelFacet
{
    CModelFacet()
	{
		m_Edge[0] = m_Edge[1] = m_Edge[2] = NULL;
	}

    int vid[3];
    Vector<T,3> n;
    CModelEdge<T> * m_Edge[3];
};

template<typename T>
class CModelGraph
{
	typedef Vector<T, 3> Vector3d;
	typedef Point<T, 3>  Point3d;

public:

    //////////////////////////////////////////////////////////////////////
    // Constructor/Destructors
	CModelGraph()
	{
		m_pEdge = NULL;
		m_EdgeSize = 0;
	}

	~CModelGraph()
	{
		//Free edges
		while (m_pEdge != NULL){
			CModelEdge<T> * tmp = m_pEdge->getNext();
			delete m_pEdge;
			m_pEdge = tmp;
		}
		//Free Nodes
		m_Nodes.clear();
	}

    //////////////////////////////////////////////////////////////////////
    // Core Function
    //convert from tri to graph
	bool doInit(model* m)
	{
		m_Tail = NULL; //tail of m_pEdge

		//create nodes
		m_Nodes.reserve(m->v_size);
		for (uint iP = 0; iP<m->v_size; iP++)
			m_Nodes.push_back(CModelNode<T>(iP));

		//create edge from triangles
		unsigned int eid = 0;
		triangle * tri = m->tris;
		m_Facets.reserve(m->t_size);
		for (uint iT = 0; iT<m->t_size; iT++){
			CModelFacet<T> F;
			for (int iV = 0; iV<3; iV++) F.vid[iV] = tri[iT].v[iV];
			Point3d& v1 = m->vertices[F.vid[0]].p;
			Point3d& v2 = m->vertices[F.vid[1]].p;
			Point3d& v3 = m->vertices[F.vid[2]].p;
			F.n = ((v2 - v1) % (v3 - v1)).normalize();

			for (int iE = 0; iE<3; iE++){ //for each edge
				int a = tri[iT].v[iE]; int b = tri[iT].v[(iE + 1) % 3]; //int c=tri[iT][(iE+2)%3];
				CModelEdge<T> * e = m_Nodes[a].getEdge(m_Nodes[b]);
				if (e == NULL)
				{
					if ((e = new CModelEdge<T>(a, b, eid++)) == NULL) return false;
					m_Nodes[a].addEdge(e); m_Nodes[b].addEdge(e);
					//push new edge to end of list
					if (m_pEdge == NULL) m_pEdge = e;
					else m_Tail->setNext(e);
					m_Tail = e; m_EdgeSize++;
					e->addFacet(iT);
					e->getV() = (m->vertices[b].p - m->vertices[a].p).normalize();
					e->addInNormal((F.n%e->getV()).normalize());
				} //end of if
				else{
					e->addFacet(iT);
					e->addInNormal((F.n % (-e->getV())).normalize());
				}

				F.m_Edge[iE] = e;
			}//end of for(iE)

			m_Facets.push_back(F);

		}//end of for(iT)
		m_Tail->setNext(NULL);
		m_Model = m;

		return true;
	}

    //////////////////////////////////////////////////////////////////////
    // Access Function

    //edges
    int getEdgeSize() const { return m_EdgeSize; }
    CModelEdge<T> * getEdges(){ return m_pEdge; }

    //facets
    vector< CModelFacet<T> >& getFacets(){ return m_Facets; }

    //////////////////////////////////////////////////////////////////////
    // Private stuff
private:
    int m_EdgeSize;
    CModelEdge<T> * m_pEdge;           //a list of edges
	CModelEdge<T> * m_Tail;            //point to the end
	vector<CModelNode<T> >  m_Nodes;    //an array of nodes
	vector<CModelFacet<T> > m_Facets;   //an array of facets
    model* m_Model;
};

}

#endif
