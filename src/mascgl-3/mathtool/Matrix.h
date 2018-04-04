//------------------------------------------------------------------------------
//  Copyright 2007-2014 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------


#ifndef _H_Matrix_H_
#define _H_Matrix_H_

#include "mathtool/Vector.h"

namespace mathtool{

    //
    // The internal storage form of transforms.  The matrices in row-major
    // order (ie  mat[row][column] )
    //
	template<typename T, int D>
	class Matrix {
	protected:

		Vector<T,D> row[D];

	public:

		Matrix() { identity(); }

		Matrix(T a11, T a12, T a21, T a22)
		{
			set(a11, a12, a21, a22);
		}

		Matrix(T a11, T a12, T a13,
			   T a21, T a22, T a23,
			   T a31, T a32, T a33)
		{
			set(a11, a12, a13, a21, a22, a23, a31, a32, a33);
		}

		Matrix(T a11, T a12, T a13, T a14,
		  	   T a21, T a22, T a23, T a24,
			   T a31, T a32, T a33, T a34,
			   T a41, T a42, T a43, T a44)
		{
			set(a11, a12, a13, a14, a21, a22, a23, a24,
				a31, a32, a33, a34, a41, a42, a43, a44);
		}

		//
		// return row i of matrix as a vector
		// Routines returning an lvalue: i.e. M[i] returns ref to row i
		//
		Vector<T,D>& operator[](int i)
		{
			if (i < 0 || i >= D){
				cerr << "Matrix ("<<D<<"x"<<D<<") row index " << i << " out of bounds" << endl;
				exit(1);
			}
			return (row[i]);
		}

		const Vector<T, D>& operator[](int i) const
		{
			if (i < 0 || i >=D){
				cerr << "Matrix (" << D << "x" << D << ") row index " << i << " out of bounds" << endl;
				exit(1);
			}
			return (row[i]);
		}

		void print(int w = 7, int p = 3) const // print with width and precision
		{
			for (int i = 0; i < D; i++)
			{
				for (int j = 0; j < D; j++)
				{
					cout << setw(w) << setprecision(p) << row[i][j] << " ";
				}
				cout << endl;
			}
		}

		void set(T a11, T a12, T a21, T a22)
		{
			assert(D > 1);
			row[0].set(a11, a12);
			row[1].set(a21, a22);
		}

		void set(T a11, T a12, T a13,
			     T a21, T a22, T a23,
			     T a31, T a32, T a33)
		{
			assert(D > 2);
			row[0].set(a11, a12, a13);
			row[1].set(a21, a22, a23);
			row[2].set(a31, a32, a33);
		}

		void set(T a11, T a12, T a13, T a14,
			     T a21, T a22, T a23, T a24,
			     T a31, T a32, T a33, T a34,
			     T a41, T a42, T a43, T a44)
		{
			assert(D > 3);
			row[0].set(a11, a12, a13, a14);
			row[1].set(a21, a22, a23, a24);
			row[2].set(a31, a32, a33, a34);
			row[3].set(a41, a42, a43, a44);
		}

		void set(T value[D][D]){
			for (int i = 0; i<D; i++) for (int j = 0; j<D; j++) row[i][j] = value[i][j];
		}


		void get(T value[D][D]) const
		{
			for (int i = 0; i<D; i++)
			for (int j = 0; j<D; j++)
				value[i][j] = row[i][j];
		}

		void identity()
		{
			for (int i = 0; i < D; i++)
			{
				for (int j = 0; j < D; j++) 
				{
					if (i == j) row[i][j] = 1;
					else row[i][j] = 0;
				}
			}
		}

		Matrix<T,D> transpose() const
		{
			Matrix<T, D> transM;

			for (int i = 0; i < D; i++)
			{
				for (int j = 0; j < D; j++) 
				{
					transM[j][i] = row[i][j];
				}
			}

			return transM;
		}

		Matrix<T, D> inv() const
		{
			if (D == 2) return inv2d();
			if (D == 3) return inv3d();
			return inv_general();
		}

		T trace() const
		{
			T sum = 0;
			for (int i = 0; i < D; i++)
			{
				sum += row[i][i];
			}

			return sum;
		}

		Matrix<T, D> operator*(T a) const
		{
			Matrix<T, D> result;

			for (int i = 0; i < D; i++)
			{
				for (int j = 0; j < D; j++)
				{
					result.row[i][j] = a * row[i][j];
				}
			}

			return result;
		}

		friend Matrix<T, D> operator+(const Matrix<T, D>& m1, const Matrix<T, D>& m2)
		{
			int i, j;
			Matrix<T, D> result;

			for (i = 0; i < D; i++)
			for (j = 0; j < D; j++)
				result.row[i][j] = m1.row[i][j] + m2.row[i][j];

			return result;
		}

		friend Matrix<T, D> operator-(const Matrix<T, D>& m1, const Matrix<T, D>& m2)
		{
			int i, j;
			Matrix<T, D> result;

			for (i = 0; i < D; i++)
			for (j = 0; j < D; j++)
				result.row[i][j] = m1.row[i][j] - m2.row[i][j];

			return result;
		}

		friend Matrix<T, D> operator*(const Matrix<T, D>& m1, const Matrix<T, D>& m2)
		{
			int i, j, rc;
			Matrix<T, D> result;

			for (i = 0; i < D; i++)
			for (j = 0; j < D; j++){
				result.row[i][j] = 0;
				for (rc = 0; rc < D; rc++)
					result.row[i][j] += m1.row[i][rc] * m2[rc][j];
			}

			return result;
		}

		friend Matrix<T, D> operator*(T a, const Matrix<T, D>& m)
		{
			Matrix<T, D> result;

			int i, j;
			for (i = 0; i < D; i++)
			for (j = 0; j < D; j++)
				result.row[i][j] = a * m.row[i][j];

			return result;
		}


		// mat times vector
		friend Vector<T, D> operator*(const Matrix<T, D>& m, const Vector<T, D>& v)
		{
			int i, j;
			T sum;
			Vector<T, D> result;

			for (i = 0; i < D; i++){
				sum = 0;
				for (j = 0; j < D; j++)
					sum += m.row[i][j] * v[j];
				result[i] = sum;
			}

			return result;
		}

		// vector times mat
		friend Vector<T, D> operator*(const Vector<T, D>& v, const Matrix<T, D>& m)
		{
			int i, j;
			T sum;
			Vector<T, D> result;

			for (j = 0; j < D; j++){
				sum = 0;
				for (i = 0; i < D; i++)
					sum += v[i] * m.row[i][j];
				result[j] = sum;
			}

			return result;
		}

		// outer product
		friend Matrix<T, D> operator&(const Vector<T, D>& v1, const Vector<T, D>& v2)
		{
			Matrix<T, D> product;

			for (int i = 0; i < D; i++)
			for (int j = 0; j < D; j++)
				product.row[i][j] = v1[i] * v2[j];

			return product;
		}


	protected:

		Matrix<T, D> inv2d() const
		{
			assert(D == 2);
			Matrix<T, D> invM;
			T d;

			d = row[0][0] * row[1][1] - row[0][1] * row[1][0];

			if (d == 0.0)
				cerr << "inverse of singular M2x2" << endl;

			invM[0][0] = row[1][1] / d;
			invM[0][1] = -row[0][1] / d;
			invM[1][0] = -row[1][0] / d;
			invM[1][1] = row[0][0] / d;

			return invM;
		}

		Matrix<T, D> inv3d() const 
		{
			assert(D == 3);
			Matrix<T,D> invM;
			T d;

			d = row[0][0] * row[1][1] * row[2][2] + row[0][1] * row[1][2] * row[2][0] +
				row[0][2] * row[2][1] * row[1][0] - row[0][2] * row[1][1] * row[2][0] -
				row[0][1] * row[1][0] * row[2][2] - row[0][0] * row[2][1] * row[1][2];

			if (d == 0.0)
				cerr << "inverse of singular M3x3" << endl;

			invM[0][0] = (row[1][1] * row[2][2] - row[1][2] * row[2][1]) / d;
			invM[0][1] = (row[0][2] * row[2][1] - row[0][1] * row[2][2]) / d;
			invM[0][2] = (row[0][1] * row[1][2] - row[0][2] * row[1][1]) / d;
			invM[1][0] = (row[1][2] * row[2][0] - row[1][0] * row[2][2]) / d;
			invM[1][1] = (row[0][0] * row[2][2] - row[0][2] * row[2][0]) / d;
			invM[1][2] = (row[0][2] * row[1][0] - row[0][0] * row[1][2]) / d;
			invM[2][0] = (row[1][0] * row[2][1] - row[1][1] * row[2][0]) / d;
			invM[2][1] = (row[0][1] * row[2][0] - row[0][0] * row[2][1]) / d;
			invM[2][2] = (row[0][0] * row[1][1] - row[0][1] * row[1][0]) / d;

			return invM;
		}

		Matrix<T, D> inv_general() const
		{
			Matrix<T, D> LU_M, invM;
			int i, j, indx[D];
			T col[D];

			LU_M = LU_Decompose(*this, indx);

			for (j = 0; j < D; j++){
				for (i = 0; i < D; i++)
					col[i] = 0.0;
				col[j] = 1.0;
				LU_back_substitution(LU_M, indx, col);
				for (i = 0; i < D; i++)
					invM[i][j] = col[i];
			}

			return invM;
		}

		//////////////////////////////////////////////////////////////////////////
		// the following matrix operations are used to find the inverse of an 
		// NxN matrix. Adapted from Numerical Recipes by (Frank) Sebastian Grassia
		//////////////////////////////////////////////////////////////////////////


		Matrix<T, D> LU_Decompose(const Matrix<T, D>& M, int *indx) const
		{
			int  i, imax, j, k;
			T big, dum, sum, temp;
			T vv[4];
			Matrix<T, D> LU_M = M;
			int N = D;

			for (i = 0; i < N; i++){
				big = 0.0;
				for (j = 0; j < N; j++)
				if ((temp = fabs(M.row[i][j])) > big)
					big = temp;
				if (big == 0.0)
					cerr << "inverse of singular M4x4" << endl;

				vv[i] = 1.0f / big;
			}

			for (j = 0; j < N; j++){
				for (i = 0; i < j; i++){
					sum = LU_M[i][j];
					for (k = 0; k < i; k++)
						sum -= LU_M[i][k] * LU_M[k][j];
					LU_M[i][j] = sum;
				}
				big = 0.0;
				for (i = j; i < N; i++){
					sum = LU_M[i][j];
					for (k = 0; k < j; k++)
						sum -= LU_M[i][k] * LU_M[k][j];
					LU_M[i][j] = sum;
					if ((dum = vv[i] * fabs(sum)) >= big){
						big = dum;
						imax = i;
					}
				}
				if (j != imax){
					for (k = 0; k < N; k++){
						dum = LU_M[imax][k];
						LU_M[imax][k] = LU_M[j][k];
						LU_M[j][k] = dum;
					}
					vv[imax] = vv[j];
				}
				indx[j] = imax;
				if (j < N - 1){
					dum = 1.0f / LU_M[j][j];
					for (i = j + 1; i < N; i++)
						LU_M[i][j] *= dum;
				}
			}

			return LU_M;
		}

		void LU_back_substitution(const Matrix<T, D>& M, int *indx, T col[]) const
		{
			int i, ii = -1, ip, j;
			T sum;
			int N = D;

			for (i = 0; i < N; i++){
				ip = indx[i];
				sum = col[ip];
				col[ip] = col[i];
				if (ii >= 0)
				for (j = ii; j < i; j++)
					sum -= M.row[i][j] * col[j];
				else if (sum)
					ii = i;
				col[i] = sum;
			}

			for (i = N - 1; i >= 0; i--){
				sum = col[i];
				for (j = i + 1; j < N; j++)
					sum -= M.row[i][j] * col[j];
				col[i] = sum / M.row[i][i];
			}
		}
	};

	template<typename T>
    class M2x2 {

		typedef Vector<T, 2> Vector2d;

    protected:

		Vector2d row[2];

    public:

		M2x2(T a11 = 1, T a12 = 0,
			T a21 = 0, T a22 = 1)
        {
            set(a11, a12, a21, a22);
        }

        //
        // return row i of matrix as a vector
        // Routines returning an lvalue: i.e. M[i] returns ref to row i
        //
        Vector2d& operator[](int i)
        {
            if(i < 0 || i > 1){
                cerr << "M2x2 row index " << i << " out of bounds" << endl;
                exit(1);
            }
            return (row[i]);
        }


        const Vector2d& operator[](int i) const
        {
            if(i < 0 || i > 1){
                cerr << "M2x2 row index " << i << " out of bounds" << endl;
                exit(1);
            }
            return (row[i]);
        }

        void print(int w = 7, int p = 3) const // print with width and precision
        {
            for(int i = 0; i < 2; i++){
                cout << setw(w) << setprecision(p) << round(row[i][0], p) << " ";
                cout << setw(w) << setprecision(p) << round(row[i][1], p);
                cout << endl;
            }
        }

		void set(T a11 = 0, T a12 = 0, T a21 = 0, T a22 = 0)
        {
            row[0].set(a11, a12);
            row[1].set(a21, a22);
        }

        template<class Type>
        void set(Type value[2][2]){
            for(int i=0;i<2;i++)
                for(int j=0;j<2;j++)
                    row[i][j]=value[i][j];
        }


		void get(T value[2][2]) const
        {
            for(int i=0;i<2;i++)
                for(int j=0;j<2;j++)
                    value[i][j]=row[i][j];
        }

        void identity()
        {
            set(1, 0, 0, 1);
        }

		M2x2<T> transpose() const
        {
			M2x2<T> transM;

            transM[0][0] = row[0][0];
            transM[1][0] = row[0][1];
            transM[0][1] = row[1][0];
            transM[1][1] = row[1][1];

            return transM;
        }

        M2x2<T> inv() const
        {
			M2x2<T> invM;
			T d;

            d = row[0][0]*row[1][1] - row[0][1]*row[1][0];

            if(d == 0.0)
                cerr << "inverse of singular M2x2" << endl;

            invM[0][0] = row[1][1] / d;
            invM[0][1] = -row[0][1] / d;
            invM[1][0] = -row[1][0] / d;
            invM[1][1] = row[0][0] / d;

            return invM;
        }

        T trace() const
        {
            return row[0][0]+row[1][1];
        }

        M2x2<T> operator*(T a) const
        {
            M2x2<T> result;

            for(int i = 0; i < 2; i++){
                result.row[i][0] = a * row[i][0];
                result.row[i][1] = a * row[i][1];
            }

            return result;
        }

		friend M2x2<T> operator+(const M2x2<T>& m1, const M2x2& m2);
		friend M2x2<T> operator-(const M2x2<T>& m1, const M2x2& m2);
		friend M2x2<T> operator*(const M2x2<T>& m1, const M2x2& m2);
		friend M2x2<T> operator*(T a, const M2x2<T>& m);


        // mat times vector
		friend Vector2d operator*(const M2x2<T>& m, const Vector2d& v);

        // vector times mat
		friend Vector2d operator*(const Vector2d& v, const M2x2& m);

        // outer product
		friend M2x2<T> operator&(const Vector2d& v1, const Vector2d& v2);
    };

	template<typename T>
    class M3x3 {

		typedef Vector<T, 3> Vector3d;

    protected:

        Vector3d row[3];

    public:

		M3x3(T a11 = 1, T a12 = 0, T a13 = 0,
			T a21 = 0, T a22 = 1, T a23 = 0,
			T a31 = 0, T a32 = 0, T a33 = 1)
        {
            set(a11, a12, a13, a21, a22, a23, a31, a32, a33);
        }

        Vector3d& operator[](int i){
            if(i < 0 || i > 2){
                cerr << "M3x3 row index " << i << " out of bounds" << endl;
                exit(1);
            }

            return (row[i]);
        }

        const Vector3d& operator[](int i) const{
            if(i < 0 || i > 2){
                cerr << "M3x3 row index " << i << " out of bounds" << endl;
                exit(1);
            }

            return (row[i]);
        }

        void print(int w = 7, int p = 3) const{  // print with width and precision

            for(int i = 0; i < 3; i++){
                cout << setw(w) << setprecision(p) << round(row[i][0], p) << " ";
                cout << setw(w) << setprecision(p) << round(row[i][1], p) << " ";
                cout << setw(w) << setprecision(p) << round(row[i][2], p);
                cout << endl;
            }
        }

		void set(T a11, T a12, T a13,
			T a21, T a22, T a23,
			T a31, T a32, T a33)
        {
            row[0].set(a11, a12, a13);
            row[1].set(a21, a22, a23);
            row[2].set(a31, a32, a33);
        }

        template<class Type>
        void set(Type value[3][3]){
            for(int i=0;i<3;i++)
                for(int j=0;j<3;j++)
                    row[i][j]=value[i][j];
        }

        template<class Type>
        void get(Type value[3][3]) const
        {
            for(int i=0;i<3;i++)
                for(int j=0;j<3;j++)
                    value[i][j]=row[i][j];
        }

        void identity(){
            set(1, 0, 0, 0, 1, 0, 0, 0, 1);
        }

        M3x3<T> transpose() const {
            M3x3<T> transM;

            transM[0][0] = row[0][0];
            transM[1][0] = row[0][1];
            transM[2][0] = row[0][2];
            transM[0][1] = row[1][0];
            transM[1][1] = row[1][1];
            transM[2][1] = row[1][2];
            transM[0][2] = row[2][0];
            transM[1][2] = row[2][1];
            transM[2][2] = row[2][2];

            return transM;
        }

        M3x3<T> inv() const {
            M3x3<T> invM;
			T d;

            d = row[0][0]*row[1][1]*row[2][2] + row[0][1]*row[1][2]*row[2][0] +
                row[0][2]*row[2][1]*row[1][0] - row[0][2]*row[1][1]*row[2][0] - 
                row[0][1]*row[1][0]*row[2][2] - row[0][0]*row[2][1]*row[1][2];

            if(d == 0.0)
                cerr << "inverse of singular M3x3" << endl;

            invM[0][0] = (row[1][1]*row[2][2] - row[1][2]*row[2][1]) / d;
            invM[0][1] = (row[0][2]*row[2][1] - row[0][1]*row[2][2]) / d;
            invM[0][2] = (row[0][1]*row[1][2] - row[0][2]*row[1][1]) / d;
            invM[1][0] = (row[1][2]*row[2][0] - row[1][0]*row[2][2]) / d;
            invM[1][1] = (row[0][0]*row[2][2] - row[0][2]*row[2][0]) / d;
            invM[1][2] = (row[0][2]*row[1][0] - row[0][0]*row[1][2]) / d;
            invM[2][0] = (row[1][0]*row[2][1] - row[1][1]*row[2][0]) / d;
            invM[2][1] = (row[0][1]*row[2][0] - row[0][0]*row[2][1]) / d;
            invM[2][2] = (row[0][0]*row[1][1] - row[0][1]*row[1][0]) / d;

            return invM;
        }

		T trace() const {
            return row[0][0]+row[1][1]+row[2][2];
        }

        M3x3<T> operator*(T a) const
        {
            M3x3<T> result;
            int i, j;

            for(i = 0; i < 3; i++)
                for(j = 0; j < 3; j++)
                    result[i][j] = a * row[i][j];

            return result;
        }

		friend M3x3<T> operator+(const M3x3<T>& m1, const M3x3<T>& m2);
		friend M3x3<T> operator-(const M3x3<T>& m1, const M3x3<T>& m2);
		friend M3x3<T> operator*(const M3x3<T>& m1, const M3x3<T>& m2);
		friend M3x3<T> operator*(T a, const M3x3<T>& m);


        // mat times vector
		friend Vector3d operator*(const M3x3<T>& m, const Vector3d& v);
        // vector times mat
		friend Vector3d operator*(const Vector3d& v, const M3x3<T>& m);
        // outer product
        friend M3x3<T> operator&(const Vector3d& v1, const Vector3d& v2);
    };

	template<typename T>
    class M4x4 
	{
    protected:
		typedef Vector<T, 4> Vector4d;
		Vector4d row[4];

    public:

		M4x4(T a11 = 1, T a12 = 0, T a13 = 0, T a14 = 0,
			T a21 = 0, T a22 = 1, T a23 = 0, T a24 = 0,
			T a31 = 0, T a32 = 0, T a33 = 1, T a34 = 0,
			T a41 = 0, T a42 = 0, T a43 = 0, T a44 = 1)
        {
            set(a11, a12, a13, a14, a21, a22, a23, a24,
                a31, a32, a33, a34, a41, a42, a43, a44);
        }

        //operator Matrix();

        Vector4d& operator[](int i){
            if(i < 0 || i > 3){
                cerr << "M4x4 row index " << i << " out of bounds" << endl;
                exit(1);
            }

            return (row[i]);
        }

        const Vector4d& operator[](int i) const
        {
            if(i < 0 || i > 3){
                cerr << "M4x4 row index " << i << " out of bounds" << endl;
                exit(1);
            }

            return (row[i]);
        }

        void print(int w = 7, int p = 3) const {  // print with width and precision

            for(int i = 0; i < 4; i++){
                cout << setw(w) << setprecision(p) << round(row[i][0], p) << " ";
                cout << setw(w) << setprecision(p) << round(row[i][1], p) << " ";
                cout << setw(w) << setprecision(p) << round(row[i][2], p) << " ";
                cout << setw(w) << setprecision(p) << round(row[i][3], p);
                cout << endl;
            }
        }

		void set(T a11, T a12, T a13, T a14,
			T a21, T a22, T a23, T a24,
			T a31, T a32, T a33, T a34,
			T a41, T a42, T a43, T a44)
        {
            row[0].set(a11, a12, a13, a14);
            row[1].set(a21, a22, a23, a24);
            row[2].set(a31, a32, a33, a34);
            row[3].set(a41, a42, a43, a44);
        }

        void set(T value[4][4]){
            for(int i=0;i<4;i++)
                for(int j=0;j<4;j++)
                    row[i][j]=value[i][j];
        }

        void get(T value[4][4]) const
        {
            for(int i=0;i<4;i++)
                for(int j=0;j<4;j++)
                    value[i][j]=row[i][j];
        }

        void identity()
        {
            set(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
        }

        M4x4<T> transpose() const
        {
			M4x4<T> transM;

            transM[0][0] = row[0][0];
            transM[1][0] = row[0][1];
            transM[2][0] = row[0][2];
            transM[3][0] = row[0][3];
            transM[0][1] = row[1][0];
            transM[1][1] = row[1][1];
            transM[2][1] = row[1][2];
            transM[3][1] = row[1][3];
            transM[0][2] = row[2][0];
            transM[1][2] = row[2][1];
            transM[2][2] = row[2][2];
            transM[3][2] = row[2][3];
            transM[0][3] = row[3][0];
            transM[1][3] = row[3][1];
            transM[2][3] = row[3][2];
            transM[3][3] = row[3][3];

            return transM;
        }

		M4x4<T> inv() const
        {
			M4x4<T> LU_M, invM;
            int i, j, indx[4];
			T col[4];

            LU_M = LU_Decompose(*this, indx);

            for(j = 0; j < 4; j++){
                for(i = 0; i < 4; i++)
                    col[i] = 0.0;
                col[j] = 1.0;
                LU_back_substitution(LU_M, indx, col);
                for(i = 0; i < 4; i++)
                    invM[i][j] = col[i];
            }

            return invM;
        }

        T trace() const
        {
            return row[0][0]+row[1][1]+row[2][2]+row[3][3];
        }

		M4x4<T> operator*(T a) const
        {

			M4x4<T> result;

            int i, j;
            for(i = 0; i < 4; i++)
                for(j = 0; j < 4; j++)
                    result[i][j] = a * row[i][j];

            return result;
        }

		friend M4x4<T> operator+(const M4x4<T>& m1, const M4x4<T>& m2);
		friend M4x4<T> operator-(const M4x4<T>& m1, const M4x4<T>& m2);
		friend M4x4<T> operator*(const M4x4<T>& m1, const M4x4<T>& m2);
		friend M4x4<T> operator*(T a, const M4x4<T>& m);


        // mat times vector
		friend Vector4d operator*(const M4x4<T>& m, const Vector4d& v);
        // vector times mat
		friend Vector4d operator*(const Vector4d& v, const M4x4<T>& m);
        // outer product
		friend M4x4 operator&(const Vector4d& v1, const Vector4d& v2);

    protected:

        //////////////////////////////////////////////////////////////////////////
        // the following matrix operations are used to find the inverse of an 
        // NxN matrix. Adapted from Numerical Recipes by (Frank) Sebastian Grassia
        //////////////////////////////////////////////////////////////////////////


		M4x4<T> LU_Decompose(const M4x4<T>& M, int *indx) const
        {
            int  i, imax, j, k;
            T big, dum, sum, temp;
            T vv[4];
			M4x4<T> LU_M = M;
            int N = 4;

            for(i = 0; i < N; i++){
                big = 0.0;
                for(j = 0; j < N; j++)
                    if((temp = fabs(M.row[i][j])) > big)
                        big = temp;
                if(big == 0.0)
                    cerr << "inverse of singular M4x4" << endl;

                vv[i] = 1.0f / big;
            }

            for(j = 0; j < N; j++){
                for(i = 0; i < j; i++){
                    sum = LU_M[i][j];
                    for(k = 0; k < i; k++)
                        sum -= LU_M[i][k] * LU_M[k][j];
                    LU_M[i][j] = sum;
                }
                big = 0.0;
                for(i = j; i < N; i++){
                    sum = LU_M[i][j];
                    for(k = 0; k < j; k++)
                        sum -= LU_M[i][k] * LU_M[k][j];
                    LU_M[i][j] = sum;
                    if((dum = vv[i]*fabs(sum)) >= big){
                        big = dum;
                        imax = i;
                    }
                }
                if(j != imax){
                    for(k = 0; k < N; k++){
                        dum = LU_M[imax][k];
                        LU_M[imax][k] = LU_M[j][k];
                        LU_M[j][k] = dum;
                    }
                    vv[imax] = vv[j];
                }
                indx[j] = imax;
                if(j < N - 1){
                    dum = 1.0f / LU_M[j][j];
                    for(i = j + 1; i < N; i++)
                        LU_M[i][j] *= dum;
                }
            }

            return LU_M;
        }

		void LU_back_substitution(const M4x4<T>& M, int *indx, T col[]) const
        {
            int i, ii = -1, ip,j;
            T sum;
            int N = 4;

            for(i = 0; i < N; i++){
                ip = indx[i];
                sum = col[ip];
                col[ip] = col[i];
                if(ii >= 0)
                    for(j = ii; j < i; j++)
                        sum -= M.row[i][j] * col[j];
                else if(sum)
                    ii = i;
                col[i] = sum;
            }

            for(i = N - 1; i >= 0; i--){
                sum = col[i];
                for(j = i + 1; j < N; j++)
                    sum -= M.row[i][j] * col[j];
                col[i] = sum / M.row[i][i];
            }
        }

    };


    /////////////////////////////////////////////////////////////////////////

    // Matrix op Matrix operations

	template<typename T> M2x2<T> operator+(const M2x2<T>& m1, const M2x2<T>& m2)
    {
        int i;
		M2x2<T> result;
        
        for(i = 0; i < 2; i++){
            result.row[i][0] = m1.row[i][0] + m2.row[i][0];
            result.row[i][1] = m1.row[i][1] + m2.row[i][1];
        }
        
        return result;
    }

	template<typename T> M2x2<T> operator-(const M2x2<T>& m1, const M2x2<T>& m2)
    {
        int i;
        M2x2<T> result;
        
        for(i = 0; i < 2; i++){
            result.row[i][0] = m1.row[i][0] - m2.row[i][0];
            result.row[i][1] = m1.row[i][1] - m2.row[i][1];
        }
        
        return result;
    }

	template<typename T> M2x2<T> operator*(T a, const M2x2<T>& m)
    {
        int i;
		M2x2<T> result;
        
        for(i = 0; i < 2; i++){
            result.row[i][0] = a * m.row[i][0];
            result.row[i][1] = a * m.row[i][1];
        }
        
        return result;
    }

	template<typename T> M2x2<T> operator*(const M2x2<T>& m1, const M2x2<T>& m2)
    {
        int i, j, rc;
		M2x2<T> result;
        
        for(i = 0; i < 2; i++)
            for(j = 0; j < 2; j++){
                result.row[i][j] = 0;
                for(rc = 0; rc < 2; rc++){
                    result.row[i][j] += m1.row[i][rc] * m2.row[rc][j];
                }
            }
            
            return result;
    }

	template<typename T> M3x3<T> operator+(const M3x3<T>& m1, const M3x3<T>& m2)
    {
        int i, j;
		M3x3<T> result;
        
        for(i = 0; i < 3; i++)
            for(j = 0; j < 3; j++)
                result.row[i][j] = m1.row[i][j] + m2.row[i][j];
            
            return result;
    }

	template<typename T> M3x3<T> operator-(const M3x3<T>& m1, const M3x3<T>& m2)
    {
        int i, j;
		M3x3<T> result;
        
        for(i = 0; i < 3; i++)
            for(j = 0; j < 3; j++)
                result.row[i][j] = m1.row[i][j] - m2.row[i][j];
            
            return result;
    }

	template<typename T> M3x3<T> operator*(T a, const M3x3<T>& m)
    {
		M3x3<T> result;
        int i, j;
        
        for(i = 0; i < 3; i++)
            for(j = 0; j < 3; j++)
                result.row[i][j] = a * m.row[i][j];
            
            return result;
    }

	template<typename T> M3x3<T> operator*(const M3x3<T>& m1, const M3x3<T>& m2)
    {
        T h=m1.row[0][0];

        int i, j, rc;
		M3x3<T> result;
        
        h=m1.row[0][0];
        for(i = 0; i < 3; i++){
            for(j = 0; j < 3; j++){
                result.row[i][j] = 0;
                for(rc = 0; rc < 3; rc++){
                    result.row[i][j] += m1.row[i][rc] * m2[rc][j];
                }
            }
        }
            
        return result;
    }

	template<typename T> M4x4<T> operator+(const M4x4<T>& m1, const M4x4<T>& m2)
    {
        int i, j;
		M4x4<T> result;
        
        for(i = 0; i < 4; i++)
            for(j = 0; j < 4; j++)
                result.row[i][j] = m1.row[i][j] + m2.row[i][j];
            
            return result;
    }

	template<typename T> M4x4<T> operator-(const M4x4<T>& m1, const M4x4<T>& m2)
    {
        int i, j;
		M4x4<T> result;
        
        for(i = 0; i < 4; i++)
            for(j = 0; j < 4; j++)
                result.row[i][j] = m1.row[i][j] - m2.row[i][j];
            
            return result;
    }

	template<typename T> M4x4<T> operator*(T a, const M4x4<T>& m)
    {
		M4x4<T> result;
        
        int i, j;
        for(i = 0; i < 4; i++)
            for(j = 0; j < 4; j++)
                result.row[i][j] = a * m.row[i][j];
            
            return result;
    }

	template<typename T> M4x4<T> operator*(const M4x4<T>& m1, const M4x4<T>& m2)
    {
        int i, j, rc;
		M4x4<T> result;
        
        for(i = 0; i < 4; i++)
            for(j = 0; j < 4; j++){
                result.row[i][j] = 0;
                for(rc = 0; rc < 4; rc++)
                    result.row[i][j] += m1.row[i][rc] * m2[rc][j];
            }
            
        return result;
    }

    /* Matrix-Vector Operations */

	template<typename T> Vector<T, 2> operator*(const M2x2<T>& m, const Vector<T, 2>& v)
    {
        int i, j;
        T sum;
		Vector<T, 2> result;
        
        for(i = 0; i < 2; i++){
            sum = 0;
            for(j = 0; j < 2; j++)
                sum += m.row[i][j] * v[j];
            result[i] = sum;
        }
        
        return result;
    }

	template<typename T> Vector<T, 3> operator*(const M3x3<T>& m, const Vector<T, 3>& v)
    {
        int i, j;
        T sum;
		Vector<T, 3> result;
        
        for(i = 0; i < 3; i++){
            sum = 0;
            for(j = 0; j < 3; j++)
                sum += m.row[i][j] * v[j];
            result[i] = sum;
        }
        
        return result;
    }

	template<typename T> Vector<T, 4> operator*(const M4x4<T>& m, const Vector<T, 4>& v)
    {
        int i, j;
        T sum;
		Vector<T, 4> result;
        
        for(i = 0; i < 4; i++){
            sum = 0;
            for(j = 0; j < 4; j++)
                sum += m.row[i][j] * v[j];
            result[i] = sum;
        }
        
        return result;
    }

	template<typename T> Vector<T, 2> operator*(const Vector<T, 2>& v, const M2x2<T>& m)
    {
        int i, j;
        T sum;
		Vector<T, 2> result;
        
        for(j = 0; j < 2; j++){
            sum = 0;
            for(i = 0; i < 2; i++)
                sum += v[i] * m.row[i][j];
            result[j] = sum;
        }
        
        return result;
    }

	template<typename T> Vector<T, 3> operator*(const Vector<T, 3>& v, const M3x3<T>& m)
    {
        int i, j;
        T sum;
		Vector<T, 3> result;
        
        for(j = 0; j < 3; j++){
            sum = 0;
            for(i = 0; i < 3; i++)
                sum += v[i] * m.row[i][j];
            result[j] = sum;
        }
        
        return result;
    }

	template<typename T> Vector<T, 4> operator*(const Vector<T, 4>& v, const M4x4<T>& m)
    {
        int i, j;
        T sum;
		Vector<T, 4> result;
        
        for(j = 0; j < 4; j++){
            sum = 0;
            for(i = 0; i < 4; i++)
                sum += v[i] * m.row[i][j];
            result[j] = sum;
        }
        
        return result;
    }

    // Outer product of v1, v2 (i.e. v1 times v2 transpose)
	template<typename T> M2x2<T> operator&(const Vector<T, 2>& v1, const Vector<T, 2>& v2)
    {
        int i, j;
		M2x2<T> product;
        
        for(i = 0; i < 2; i++)
            for(j = 0; j < 2; j++)
                product.row[i][j] = v1[i] * v2[j];
            
        return product;
    }

	template<typename T> M3x3<T> operator&(const Vector<T, 3>& v1, const Vector<T, 3>& v2)
    {
		M3x3<T> product;
        
        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
                product.row[i][j] = v1[i] * v2[j];
            
        return product;
    }

	// outer product
	template<typename T> M4x4<T> operator&(const Vector<T, 4>& v1, const Vector<T, 4>& v2)
	{
		M4x4<T> product;

		for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			product.row[i][j] = v1[i] * v2[j];

		return product;
	}

} //end of mathtool namespace

#endif //_H_Matrix_H_
