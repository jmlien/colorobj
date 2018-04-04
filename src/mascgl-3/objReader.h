//------------------------------------------------------------------------------
//  Copyright 2007-2014 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _OBJ_READER_H_
#define _OBJ_READER_H_

#include <string>
#include <fstream>
#include <iostream>
#include <list>
#include <vector>
#include <ctype.h>
#include <math.h>

using namespace std;

#include <mathtool/Point.h>
#include <mathtool/Quaternion.h>
#include <mathtool/Vector.h>
using namespace mathtool;

typedef double REAL;

// Typedef common used vector type
namespace mathtool
{
          typedef Point<REAL, 2>   Point2d;
          typedef Point<REAL, 3>   Point3d;
          typedef Vector<REAL, 2> Vector2d;
          typedef Vector<REAL, 3> Vector3d;
          typedef Vector<REAL, 2> Vector2d;
          typedef Matrix<REAL, 2> Matrix2x2;
          typedef Matrix<REAL, 3> Matrix3x3;
          typedef Matrix<REAL, 4> Matrix4x4;
}

namespace objreader
{

  //convert a string to a list of tokens
  inline list<string> tokenize(char * tmp, const char * ignore)
  {
  	list<string> tokens;
  	char * tok = strtok(tmp, ignore);
  	while (tok != NULL) {
  		tokens.push_back(tok);
  		tok = strtok(NULL, ignore);
  	}
  	return tokens;
  }

  inline void getAllLabels(istream& in, list<list<string> >& tokens)
  {
  	const int size = 1024;
  	char * tmp = new char[size];

  	while (!in.eof()) {
  		in.getline(tmp, size);
  		//check for termination
  		if (string(tmp).find("}") != string::npos) //found "}"
  			break;

  		list<string> tok = tokenize(tmp, " \t[]()<>,={");
  		if (tok.empty())
  			continue;
  		string label = tok.front();
  		if (label[0] == '#')
  			continue; //comment

  		tokens.push_back(tok);
  	} //end while

  	delete[] tmp;
  }

	template<typename T>
	class Vpt
	{
	public:
		Vpt(){ x = y = z = 0; }
		Vpt(const T * v){ x = v[0]; y = v[1]; z = v[2]; }
		T x, y, z;
	};

	template<typename T>
	class V
	{
	public:
		V(){ x = y = z = 0; }
		V(const T * v){ x = v[0]; y = v[1]; z = v[2]; }

		void normalize()
		{
			T norm = (T)sqrt(x*x + y*y + z*z);
			x /= norm;
			y /= norm;
			z /= norm;
		}

		T x, y, z;
	};

	class polygon
	{
	public:
		list<int> pts;
		list<int> textures;
		list<int> normals;
	};

	template<typename T>
	class objModel
	{
	public:
		objModel() { }

		void compute_v_normal()
		{
			//check if normal information is valid from the obj file
			if (normals.empty()){ //compute normal

				normals = vector< V<T> >(pts.size(), V<T>());

				for (auto i = polys.begin(); i != polys.end(); i++)
				{
					//get 3 points, compute normal and assign to all vertices
					auto pi = i->pts.begin();
					vector< Point<T, 3> > v3;
					for (; pi != i->pts.end(); pi++){
						Vpt<T>& pt = pts[*pi];
						Point3d pos(pt.x, pt.y, pt.z);
						v3.push_back(pos);
						if (v3.size() == 3) break; //we've collected 3 points
					}
					//compute normal
					Vector<T, 3> n = ((v3[1] - v3[0]) % (v3[2] - v3[0]));

					//copy normals
					pi = i->pts.begin();
					for (; pi != i->pts.end(); pi++){

						normals[*pi].x += n[0];
						normals[*pi].y += n[1];
						normals[*pi].z += n[2];

						//normal index is the same as the vertex index
						i->normals.push_back(*pi);

					}//end copying normals

				}//end looping polygons
			}
			else{ // use the information provided
				//do nothing
			}

			//normalize
			for (auto i = normals.begin(); i != normals.end(); i++)
			{
				i->normalize();
			}
		}

		void compute_v_textcoord()
		{
			if (textcoords.empty() == false) return; //nothing to do...already setup
			//get a list of polygons, setup the text coordinate to be all zeros...
			textcoords = vector< V<T> >(1, V<T>());
			for (list<polygon>::iterator i = polys.begin(); i != polys.end(); i++)
			{
				//all texture coord point to this dummy textcoords
				for (list<int>::iterator pi = i->pts.begin(); pi != i->pts.end(); pi++)
				{
					i->textures.push_back(0);
				}
			}//end for i
		}

		vector< Vpt<T> > pts;
    vector< Vpt<T> > colors;
		vector< V<T> > normals;
		vector< V<T> > textcoords;
		list<polygon> polys;
	};

	template<typename T>
	class objReader
	{
	public:
		objReader(const string& name){ m_filename = name; }

		bool Read(){
			ifstream in(m_filename.c_str());
			if (!in.good()){
				cerr << "Error: Can't open file " << m_filename << endl;
				return false;
			}
			bool r = Read(in);
			in.close();
			return r;
		}

		const objModel<T>& getModel() const { return data; }
		objModel<T>& getModel() { return data; }

	private:

		bool Read(istream& in)
		{
      list<list<string> > tokens;
      getAllLabels(in, tokens);

			//read pts
			for( list<string> & token : tokens)
			{
        string label=token.front();
        token.pop_front();
				if (label == "v"){
          vector<T> values;
          for(string& tmp : token) values.push_back(atof(tmp.c_str()));
					Vpt<T> pt;
          pt.x=values[0];
          pt.y=values[1];
          pt.z=values[2];
					data.pts.push_back(pt);
          if(values.size()>3)
          {
            Vpt<T> pt;
            pt.x=values[3];
            pt.y=values[4];
            pt.z=values[5];
            data.colors.push_back(pt);
          }
				}
				else if (label == "vn"){
					V<T> pt;
					in >> pt.x >> pt.y >> pt.z;
					data.normals.push_back(pt);
				}
				else if (label == "vt"){
					V<T> pt;
					in >> pt.x >> pt.y;
					data.textcoords.push_back(pt);
				}
        else if(label == "usemtl")
        {
              //TODO
        }
        else if (label == "mtllib")
        {
            string mtl_path=token.front();
            // cout << "mtllib " << mtl_path << std::endl;
            // this->m_mtl_paths.push_back(mtl_path);
            // // assuming obj and mtl are in the same dir.
            // string full_path = masc::path::DirPath(this->m_filename) + "/" + mtl_path;
            // this->m_mtl_reader.Read(full_path);
        }
        else if (label == "f"){
          polygon poly;

          for(string& tmp : token) //each token value
          {
            int pos1 = (int)tmp.find('/');
  					int pos2 = (int)tmp.rfind('/');

  					string field1 = tmp.substr(0, pos1);
  					string field2 = tmp.substr(pos1 + 1, pos2 - pos1 - 1);
  					string field3 = tmp.substr(pos2 + 1);

  					if (pos1 < 0 || pos2 < 0) //has no "/"
  					{
  						field2.clear();
  						field3.clear();
  					}
  					else if (pos1 == pos2 && pos1 >= 0) //has only on "/"
  					{
  						field3.clear();
  					}

  					int id_v = atoi(field1.c_str()) - 1;
  					poly.pts.push_back(id_v);


  					if (field2.empty() == false)
  					{
  						int id_t = atoi(field2.c_str()) - 1;
  						if (id_t >= 0) poly.textures.push_back(id_t);
  					}

  					if (field3.empty() == false)
  					{
  						int id_n = atoi(field3.c_str()) - 1;
  						if (id_n >= 0) poly.normals.push_back(id_n);
  					}
          }//end parsing face

          if (!poly.pts.empty()) data.polys.push_back(poly);
        }
				else{}
			}

			data.compute_v_normal();
			data.compute_v_textcoord();

			return true;
		}

		string m_filename;
		objModel<T> data;
	};
}

#endif //_OBJ_READER_H_
