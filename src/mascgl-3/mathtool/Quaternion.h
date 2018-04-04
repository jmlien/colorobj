//------------------------------------------------------------------------------
//  Copyright 2007-2014 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#if !defined( _QUATERNION_H_NPRM_ )
#define _QUATERNION_H_NPRM_

#include "mathtool/Matrix.h"
#include "mathtool/Vector.h"
#include "mathtool/Point.h"

namespace mathtool
{
	template<class T>
    class Quaternion {
    public:

		//---------------------------------------------------------------------
		typedef Vector<T, 3> Vector3d;
		typedef Point<T, 3>  Point3d;
		typedef Matrix<T, 3> Matrix3x3;
		typedef Matrix<T, 4> Matrix4x4;
		//---------------------------------------------------------------------

        Quaternion(){ m_s=1; m_v.set(0,0,0); }
		Quaternion(T w, T i, T j, T k){ m_s = w; m_v.set(i, j, k); }
        Quaternion( T s, Vector3d v ){ m_v=v; m_s=s; }

		static Quaternion get(T radian, const Vector3d& axis)
		{ 
			Quaternion q;
			T cos_r = cos(radian/2);
			T sin_r = sin(radian/2);
			q.m_v = axis*sin_r; 
			q.m_s = cos_r;

			return q;
		}

        Quaternion( const Quaternion & q){ *this=q; }
        Quaternion( T r[3]){  //compute rotational quaternion
            //convert rot to quaternion
            T sin_r=(T)sin(r[0]/2);
            T cos_r=(T)cos(r[0]/2);
            Quaternion qx(cos_r,sin_r*Vector3d(1,0,0));
            sin_r=(T)sin(r[1]/2);
            cos_r=(T)cos(r[1]/2);
            Quaternion qy(cos_r,sin_r*Vector3d(0,1,0));
            sin_r=(T)sin(r[2]/2);
            cos_r=(T)cos(r[2]/2);
            Quaternion qz(cos_r,sin_r*Vector3d(0,0,1));
            *this=(qz*qy*qx).normalize();
        }

		/// Construct from a transformation matrix m4
		Quaternion(const Matrix4x4& m4)
		{
			Matrix3x3 m3(m4[0][0], m4[0][1], m4[0][2],
						 m4[1][0], m4[1][1], m4[1][2], 
						 m4[2][0], m4[2][1], m4[2][2]);
			set(m3);
		}

		/// Construct from a rotation matrix m3
		Quaternion(const Matrix3x3& m3)
		{
			set(m3);
		}

        ////////////////////////////////////////////////////////////////////////
        // Operations for Quaternion
        Quaternion operator*(const Quaternion & q) const {
            T s=q.m_s*m_s-q.m_v*m_v;
            Vector3d v=q.m_v*m_s+m_v*q.m_s+m_v%q.m_v;
            return Quaternion(s,v);
        }
        Quaternion operator*(const Vector3d & v) const { return *this*Quaternion(0,v); }
        Quaternion operator/(T s) const { return Quaternion(m_s/s,m_v/s); }
        Quaternion & operator=(const Quaternion & q){ set(q.m_s,q.m_v); return *this; }
        Quaternion operator+(const Quaternion & q) const { return Quaternion(m_s+q.m_s,m_v+q.m_v); }
        Quaternion operator-(const Quaternion & q) const { return Quaternion(m_s-q.m_s,m_v-q.m_v); }
        Quaternion operator-() const { return Quaternion(m_s,-m_v); }
        friend Quaternion operator*(const Vector3d & v, const Quaternion & q);
        friend istream& operator>>(istream & in, Quaternion & q );
        friend ostream& operator<<(ostream & out, const Quaternion & q );

        //////////////////////////////////////////////////////////////////////////
        //Normalization
        Quaternion normalize() const { 
            Quaternion q(*this);
            T l=q.norm();
            q=q/l;
            return q;
        }

        T norm() const { return sqrt(normsqr()); }
        T normsqr() const { return m_v.normsqr()+sqr(m_s); }

        //////////////////////////////////////////////////////////////////////////
        //Access

		Matrix3x3 getMatrix() const
		{
            T x_2=2*sqr(m_v[0]); T y_2=2*sqr(m_v[1]); T z_2=2*sqr(m_v[2]);
            T xy=2*m_v[0]*m_v[1]; T yz=2*m_v[1]*m_v[2]; T zx=2*m_v[2]*m_v[0]; 
            T sx=2*m_s*m_v[0]; T sy=2*m_s*m_v[1]; T sz=2*m_s*m_v[2]; 
			return Matrix3x3(1 - y_2 - z_2, xy - sz, zx + sy,
                             xy+sz, 1-x_2-z_2, yz-sx,
                             zx-sy, yz+sx, 1-x_2-y_2);
        }
    

        void set(T s,const Vector3d & v){ m_v=v; m_s=s; }
        void set(T q1, T q2, T q3, T q4){ m_s=q1; m_v.set(q2,q3,q4); }
		/// initialize from a rotation matrix m3
		void set(const Matrix3x3& m3)
		{
			// Compute trace of matrix 't'
			T Tt = m3.trace()+1;

			T S, X, Y, Z, W;

			if (Tt > 0.00000001f) // to avoid large distortions!
			{
				S = sqrt(Tt) * 2.f;
				X = (m3[1][2] - m3[2][1]) / S;
				Y = (m3[2][0] - m3[0][2]) / S;
				Z = (m3[0][1] - m3[1][0]) / S;
				W = 0.25f * S;
			}
			else
			{
				if (m3[0][0] > m3[1][1] && m3[0][0] > m3[2][2])
				{
					// Column 0 :
					S = sqrt(1.0f + m3[0][0] - m3[1][1] - m3[2][2]) * 2.f;
					X = 0.25f * S;
					Y = (m3[0][1] + m3[1][0]) / S;
					Z = (m3[2][0] + m3[0][2]) / S;
					W = (m3[1][2] - m3[2][1]) / S;
				}
				else if (m3[1][1] > m3[2][2])
				{
					// Column 1 :
					S = sqrt(1.0f + m3[1][1] - m3[0][0] - m3[2][2]) * 2.f;
					X = (m3[0][1] + m3[1][0]) / S;
					Y = 0.25f * S;
					Z = (m3[1][2] + m3[2][1]) / S;
					W = (m3[2][0] - m3[0][2]) / S;
				}
				else
				{   // Column 2 :
					S = sqrt(1.0f + m3[2][2] - m3[0][0] - m3[1][1]) * 2.f;
					X = (m3[2][0] + m3[0][2]) / S;
					Y = (m3[1][2] + m3[2][1]) / S;
					Z = 0.25f * S;
					W = (m3[0][1] - m3[1][0]) / S;
				}
			}

			m_s = W; 
			m_v.set( - X, -Y, -Z );
		}

        const Vector3d& getComplex() const { return m_v; }
        T getReal() const { return m_s; }

		//
		/// Do the rotation of vector 'v' with the quaternion
		Vector3d rotate(const Vector3d& v) const
		{
			return ((*this)*v*(-(*this))).getComplex();

			//const Vector3d& q_vec = m_v;
			//return v + (q_vec*2.f).cross(q_vec.cross(v) + v*coeff[0]);
			//return v + (q_vec*2.f) % (q_vec%v + v*m_s);
		}

		/// Do the rotation of point 'p' with the quaternion
		Point3d rotate(const Point3d& p) const
		{
			Vector3d v = rotate(Vector3d(p.get()));
			return Point3d(v.get());
		}


    private:
        Vector3d m_v;
        T m_s;
    };

	template<typename T> istream& operator>>(istream & in, Quaternion<T> & q)
	{
		return in;
	}

	template<typename T> ostream& operator<<(ostream & out, const Quaternion<T> & q)
	{
		return out;
	}

	template<typename T> Quaternion<T> operator*(const Vector<T,3> & v, const Quaternion<T> & q2)
	{
		Quaternion<T> q1(0, v);
		return q1*q2;
	}

}//emd of nprmlib

#endif //#if !defined( _QUATERNION_H_NPRM_ ) 

