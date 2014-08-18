/******************************************************************************
 * The Sven Project - All Rights Reserved
 *
 * 15 August 2014
 * ~ ry
 *
 ******************************************************************************/

#include "sven/CoreMath.hxx"
#include "gtest/gtest.h"

using namespace sven;

//Vector Tests ----------------------------------------------------------------

/*
TEST(Vector, SubscriptAndSize)
{
  Vector x(8); 
  //EXPECT_EQ(ObjectState::Materializing, x.state());

  x = Vector::Zero(8); 
  //EXPECT_EQ(ObjectState::SolidState, x.state());
  EXPECT_EQ(x.n(), 8UL);

  for(size_t i=0; i<x.n(); ++i) { EXPECT_DOUBLE_EQ(x(i), 0); }

  x(4) = 7.0;
  EXPECT_DOUBLE_EQ(x(4), 7.0);
}
*/

TEST(Vector, Add)
{
  Vector x{1,3,5,7,9}, y{0,2,4,6,8};

  //EXPECT_EQ(ObjectState::SolidState, x.state());
  //EXPECT_EQ(ObjectState::SolidState, y.state());

  EXPECT_DOUBLE_EQ(7, x(3));
  EXPECT_DOUBLE_EQ(2, y(1));

  Vector z = x + y;

  //Something like this may or may not succeed dep on thread timing
  //usleep(1);
  //EXPECT_EQ(ObjectState::Materializing, z.state());

  z.wait();
  EXPECT_DOUBLE_EQ(9, z(2));
  //EXPECT_EQ(ObjectState::SolidState, z.state());

  Vector a = !x;
  a += y;

  EXPECT_TRUE(a == z);

}

TEST(Vector, Sub)
{
  Vector x{2,4,6,8,10}, y{1,2,3,4,5};

  Vector z = x - y;

  EXPECT_TRUE(z == y);

  Vector a = x;
  a -= y;

  EXPECT_TRUE(a == z);
}

TEST(Vector, DivScale)
{
  Vector x{2,4,6,8,10};
  Scalar s{2};

  Vector y = x / s;

  y.wait();
  EXPECT_DOUBLE_EQ(1, y(0));
  EXPECT_DOUBLE_EQ(2, y(1));
  EXPECT_DOUBLE_EQ(3, y(2));
  EXPECT_DOUBLE_EQ(4, y(3));
  EXPECT_DOUBLE_EQ(5, y(4));
  
  Vector z = x;
  z /= s;

  EXPECT_TRUE(z == y);
}

TEST(Vector, MulScale)
{
  Vector x{2,4,6,8,10};
  Scalar s{2};

  Vector y = x * s;

  Vector xx{4,8,12,16,20};

  EXPECT_TRUE(y == xx);

  Vector z = x;
  z *= s;
  EXPECT_TRUE(z == y);
}

TEST(Vector, Dot)
{
  Vector x{1,3,5,7,9}, y{0,2,4,6,8};

  Scalar z = x * y;

  z.wait();
  EXPECT_EQ(1*0+2*3+4*5+6*7+8*9, z()); 
  
}

TEST(Vector, Norm)
{
  Vector x{1,3,5,7,9};

  Scalar nx = norm(x);
  Scalar _nx_{12.845232578665129};

  EXPECT_TRUE(_nx_ == nx);
}

//Scalar Tests ----------------------------------------------------------------

TEST(Scalar, Add)
{
  Scalar a{2}, b{2};

  Scalar c = a + b;

  c.wait();
  EXPECT_DOUBLE_EQ(4, c());

  Scalar d = a;
  d += b;

  EXPECT_TRUE(c == d);
}

TEST(Scalar, Sub)
{
  Scalar a{4}, b{2};
  
  Scalar c = a - b;

  c.wait();
  EXPECT_DOUBLE_EQ(2, c());

  Scalar d = a;
  d -= b;

  EXPECT_TRUE(c == d);
}

TEST(Scalar, Mul)
{
  Scalar a{4}, b{2};
  
  Scalar c = a * b;

  c.wait();
  EXPECT_DOUBLE_EQ(8, c());
  
  Scalar d = a;
  d *= b;

  EXPECT_TRUE(c == d);
}

TEST(Scalar, Div)
{
  Scalar a{5}, b{2};
  
  Scalar c = a / b;

  c.wait();
  EXPECT_DOUBLE_EQ(2.5, c());
  
  Scalar d = a;
  d /= b;

  EXPECT_TRUE(c == d);
}

//Matrix Tests ----------------------------------------------------------------

TEST(Matrix, MulVecBasic)
{
  Matrix A = Matrix::Identity(5,5);
  Vector x{1,2,3,4,5};

  Vector Ax = A * x;

  EXPECT_TRUE(Ax == x);
}

TEST(Matrix, MulVec)
{
  Matrix A = Matrix(5,5,
      {98,3,40,4,5,
       2,52,15,6,45,
       87,63,5,89,83,
       32,79,82,98,74,
       9,28,73,4,7});

  Vector x{3, 6, 90, 23, 64};

  Vector Ax = A * x;

  Vector _Ax_{4324,  4686,  8448, 14940,  7305};

  EXPECT_TRUE(_Ax_ == Ax);
}

TEST(Matrix, MulId)
{
  Matrix A = Matrix::Identity(5,5);
  
  Matrix B = Matrix(5,5,
      {98,3,40,4,5,
       2,52,15,6,45,
       87,63,5,89,83,
       32,79,82,98,74,
       9,28,73,4,7});

  Matrix AB = A * B;

  EXPECT_TRUE(AB == B);
}

TEST(Matrix, Mul)
{
  Matrix A = Matrix(5,5,
      {98,3,40,4,5,
       2,52,15,6,45,
       87,63,5,89,83,
       32,79,82,98,74,
       9,28,73,4,7});
  
  Matrix B = Matrix(5,5,
      {98,3,40,4,5,
       2,52,15,6,45,
       87,63,5,89,83,
       32,79,82,98,74,
       9,28,73,4,7});

  Matrix AB = A * B;

  Matrix _AB_(5,5,
    {13263,  3426,  4858,  4382,  4276,
     2202,  5389,  4712,  2423,  4354,
     12682, 13207, 17807, 10225, 10852,
     14230, 19184, 16313, 17800, 18291,
     7480,  6594,  1984,  7121,  7709});

  EXPECT_TRUE(AB == _AB_);

  EXPECT_DOUBLE_EQ(16313, _AB_(3, 2));

  Matrix Z = !A;
  Z *= B;

  EXPECT_TRUE(Z == AB);
}

TEST(Matrix, ColRef)
{
  Matrix A = Matrix::Identity(5,5);    
  Vector x1{0,0,1,0,0},
         x2{2,2,2,2,2};
  
  Column c2 = A.C(2);
  EXPECT_TRUE(c2 == x1);

  A.C(2) = Vector{2,2,2,2,2};

  EXPECT_TRUE(c2 == x2);

  Vector x3 = A.C(2);

  EXPECT_TRUE(x3 == x2);
}

//Sparse Matrix Tests ---------------------------------------------------------

TEST(SparseMatrix, Basics)
{
  SparseMatrix A = SparseMatrix::Identity(5,5,3);

  EXPECT_EQ(A.m(), 5UL);
  EXPECT_EQ(A.n(), 5UL);
  EXPECT_EQ(A.z(), 3UL);

  EXPECT_DOUBLE_EQ(1, A(4, 4));
}

TEST(SparseMatrix, SparseMatrixVecMul)
{
  SparseMatrix A = SparseMatrix::Identity(5,5,3);
  Vector x{1,2,3,4,5};

  Vector Ax = A * x;

  EXPECT_TRUE(Ax == x);
}

TEST(SparseMatrix, SparseMatrixVecMul2)
{
  SparseMatrix A(5, 5, 3,
      {1,   2,       2,       3,           2      },
      {0,   1,  4,   2,  4,   1,  3,  4,   0,  4  },
      {1.0, 1.0,7.0, 1.0,3.3, 2.2,1.0,1.1, 0.4,1.0});
  
  Vector x{1,2,3,4,5};

  Vector Ax = A * x;

  Vector _Ax_{1, 37, 19.5, 13.9, 5.4};

  EXPECT_TRUE(Ax == _Ax_);
}
