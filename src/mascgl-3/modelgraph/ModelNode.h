//------------------------------------------------------------------------------
//  Copyright 2007-2014 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------


#ifndef _MODELGRAPH_NODEL_H_
#define _MODELGRAPH_NODEL_H_

#include <list>
using namespace std;

//////////////////////////////////////////////////////////////////////
// Include ModelGraph Headers.
#include "ModelEdge.h"
namespace masc
{
template<typename T>
class CModelNode
{
public:
    //////////////////////////////////////////////////////////////////////
    // Constructor/Destructors
    CModelNode(int Key)
	{
		m_Key = Key;
	}

    //////////////////////////////////////////////////////////////////////
    // Access Function
    CModelEdge<T> * getEdge(const CModelNode& nb) const
    {
        int Key=nb.m_Key;
        //linear search
        for(auto i=m_Edges.begin();i!=m_Edges.end();i++)
            if( (*i)->isEndPt(Key)==true ) return *i;
        return NULL;
    }

	const list<CModelEdge<T> *>& getEdges() const
    {
        return m_Edges;
    }

	void addEdge(CModelEdge<T> * e){
        if(e==NULL) return;
        m_Edges.push_back(e);
    }

    int getKey() const { return m_Key; }

    //////////////////////////////////////////////////////////////////////
    // Private Stuff
private:

    int m_Key;
	list<CModelEdge<T>*> m_Edges;
};

}

#endif // _MODELGRAPH_NODEL_H_
