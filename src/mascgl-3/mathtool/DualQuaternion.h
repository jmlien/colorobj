
//------------------------------------------------------------------------------
//  Copyright 2007-2014 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

//
// This code is heavily adopted from Rodolphe Vaillant's implmentation of
// Dual Quaternions skinning tutorial and C++ codes
// http://rodolphe-vaillant.fr/?e=29
//

#pragma once

#include "mathtool/Quaternion.h"
#include "mathtool/Point.h"

namespace mathtool{

/** @class DualQuaternion
    @brief Representation of a dual quaternion to express rotation and translation

    A dual quaternion (DQ) is based on the algebra of dual numbers. Dual numbers
    are somewhat close to complex numbers (a + ib) as they are writen :
    nd + eps * d where nd is the non-dual part and d the dual part and
    (eps)^2=0.

    Dual quaternion are represented as follow : q0 + eps * qe where q0
    is the non-dual part (a quaternion) and qe the dual part (another quaternion)

    With dual quaternion we can express a rotation and a translation. This
    enables us to substitute rigid transformations matrices to dual quaternions
    and transform a point with the method 'transform()'

    As a dual quaternions is the sum of two quaternions we have to store eight
    coefficients corresponding to the two quaternions.

    To move a point with a rigid transformation (i.e. solely composed
    of a translation and a rotation) you just need to construct the DQ with a
    quaternion wich express the rotation and a translation vector. You can
    now translate and rotate the point at the same time with 'transform()'.

    Linear blending of dual quaternions (DLB) is possible (dq0*w0 + dq1*w1 ...)
    where w0, w1 ... wn are scalar weights whose sum equal one. The weights
    defines the influence of each transformations expressed by the dual
    quaternions dq0, dq1 ... dqn.
    N.B : this is often used to compute mesh deformation for animation systems.

    You can compute DLB with the overloaded operators (+) and (*) and use
    the method transform() of the resulting dual quaternion to deform a point
    according to the DLB.

    @note Article: "Geometric skinning with approximate dual quaternion blending"
 */
template<typename T>
class DualQuaternion 
{
    public:

	//---------------------------------------------------------------------
	typedef Vector<T, 3> Vector3d;
	typedef Point<T, 3>  Point3d;
	typedef Matrix<T, 3> Matrix3x3;
	typedef Matrix<T, 4> Matrix4x4;
	typedef Quaternion<T> Quat;
	//---------------------------------------------------------------------

    // -------------------------------------------------------------------------
    /// @name Constructors
    // -------------------------------------------------------------------------

    /// Default constructor generates a dual quaternion with no translation
    /// and no rotation either
    DualQuaternion()
    {
		DualQuaternion res = dual_quat_from(Quat(), Vector3d(0.f, 0.f, 0.f));
        *this = res;
    }


    /// Fill directly the dual quaternion with two quaternion for the non-dual
    /// and dual part
	DualQuaternion(const Quat& q0, const Quat& qe)
    {
        _quat_0 = q0;
        _quat_e = qe;
    }

    /// Construct a dual quaternion with a quaternion 'q' which express the
    /// rotation and a translation vector
	DualQuaternion(const Quat& q, const Vector3d& t)
    {
		DualQuaternion res = dual_quat_from(q, t);
        *this = res;
    }

    /// Construct from rigid transformation 't'
	DualQuaternion(const Matrix4x4& t)
    {
		Quat q(t);
		Vector3d translation(t[0][3], t[1][3], t[2][3]);
		DualQuaternion res = dual_quat_from(q, translation);
        *this = res;
    }


    // -------------------------------------------------------------------------
    /// @name Methods
    // -------------------------------------------------------------------------

    void normalize()
    {
        T norm = _quat_0.norm();
        _quat_0 = _quat_0 / norm;
        _quat_e = _quat_e / norm;
    }

    /// Transformation of point p with the dual quaternion
	Point3d transform(const Point3d& p) const
    {
        // As the dual quaternions may be the results from a
        // linear blending we have to normalize it :
        T norm = _quat_0.norm();
		Quat qblend_0 = _quat_0 / norm;
		Quat qblend_e = _quat_e / norm;

        // Translation from the normalized dual quaternion equals :
        // 2.f * qblend_e * conjugate(qblend_0)
		Vector3d v0 = qblend_0.getComplex();
		Vector3d ve = qblend_e.getComplex();
		Vector3d trans = (ve*qblend_0.getReal - v0*qblend_e.getReal() + v0%ve) * 2.f;

        // Rotate
        return qblend_0.rotate(p) + trans;
    }

    /// Rotate a vector with the dual quaternion
	Vector3d rotate(const Vector3d& v) const
    {
        Quat tmp = _quat_0;
        tmp=tmp.normalize();
        return tmp.rotate(v);
    }

	DualQuaternion dual_quat_from(const Quat& q, const Vector3d& t) const
    {
		const Vector3d & qv = q.getComplex();
		T qw = q.getReal();
		T qi = qv[0];
		T qj = qv[1];
		T qk = qv[2];

		T w = -0.5f*(t[0] * qi + t[1] * qj + t[2] * qk);
		T i = 0.5f*(t[0] * qw + t[1] * qk - t[2] * qj);
		T j = 0.5f*(-t[0] * qk + t[1] * qw + t[2] * qi);
		T k = 0.5f*(t[0] * qj - t[1] * qi + t[2] * qw);

        return DualQuaternion(q, Quaternion(w, i, j, k) );
    }

    /// Convert the dual quaternion to a homogenous matrix
    /// N.B: Dual quaternion is normalized before conversion
    Matrix4x4 to_transformation()
    {
		Vector3d t;
        T norm = _quat_0.norm();

        // Rotation matrix from non-dual quaternion part
		Matrix3x3 m = (_quat_0 / norm).getMatrix();

        // translation vector from dual quaternion part:
		const Vector3d & qev = _quat_e.getComplex();
		T qew = _quat_e.getReal();
		T qei = qev[0];
		T qej = qev[1];
		T qek = qev[2];

		const Vector3d & q0v = _quat_0.getComplex();
		T q0w = _quat_0.getReal();
		T q0i = q0v[0];
		T q0j = q0v[1];
		T q0k = q0v[2];

		t[0] = 2.f*(-qew*q0i + qei*q0w - qej*q0k + qek*q0j) / norm;
		t[1] = 2.f*(-qew*q0j + qei*q0k + qej*q0w - qek*q0i) / norm;
		t[2] = 2.f*(-qew*q0k - qei*q0j + qej*q0i + qek*q0w) / norm;

		return Matrix4x4(m, t);
    }

    // -------------------------------------------------------------------------
    /// @name Operators
    // -------------------------------------------------------------------------

	DualQuaternion operator+(const DualQuaternion& dq) const
    {
		return DualQuaternion(_quat_0 + dq._quat_0, _quat_e + dq._quat_e);
    }

	DualQuaternion operator*(T scalar) const
    {
		return DualQuaternion(_quat_0 * scalar, _quat_e * scalar);
    }

    /// Return a dual quaternion with no translation and no rotation
    static DualQuaternion identity()
    {
		return DualQuaternion(Quat(1.f, 0.f, 0.f, 0.f), Vector3d(0.f, 0.f, 0.f));
    }

    // -------------------------------------------------------------------------
    /// @name Getters
    // -------------------------------------------------------------------------

	Quat get_dual_part() const { return _quat_e; }

	Quat get_non_dual_part() const { return _quat_0; }

	Quat translation() const { return _quat_e; }

	Quat rotation() const { return _quat_0; }

	void set_rotation(const Quat& q){ _quat_0 = q; }

    // -------------------------------------------------------------------------
    /// @name Attributes
    // -------------------------------------------------------------------------

private:
    /// Non-dual part of the dual quaternion. It also represent the rotation.
    /// @warning If you want to compute the rotation with this don't forget
    /// to normalize the quaternion as it might be the result of a
    /// dual quaternion linear blending
    /// (when overloaded operators like '+' or '*' are used)
	Quat _quat_0;

    /// Dual part of the dual quaternion which represent the translation.
    /// translation can be extracted by computing
    /// 2.f * _quat_e * conjugate(_quat_0)
    /// @warning don't forget to normalize quat_0 and quat_e :
    /// quat_0 = quat_0 / || quat_0 || and quat_e = quat_e / || quat_0 ||
	Quat _quat_e;
};

}// end of namespace

