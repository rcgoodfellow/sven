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

TEST(Vector, SubscriptAndSize)
{
  Vector x(8); 
  EXPECT_EQ(ObjectState::Materializing, x.state());

  x = Vector::Zero(8); 
  EXPECT_EQ(ObjectState::SolidState, x.state());
  EXPECT_EQ(x.n(), 8);

  for(size_t i=0; i<x.n(); ++i) { EXPECT_DOUBLE_EQ(x(i), 0); }

  x(4) = 7.0;
  EXPECT_DOUBLE_EQ(x(4), 7.0);
}

TEST(Vector, Add)
{
  Vector x{1,3,5,7,9}, y{0,2,4,6,8};

  EXPECT_EQ(ObjectState::SolidState, x.state());
  EXPECT_EQ(ObjectState::SolidState, y.state());

  EXPECT_DOUBLE_EQ(7, x(3));
  EXPECT_DOUBLE_EQ(2, y(1));

  Vector z = x + y;

  //Something like this may or may not succeed dep on thread timing
  //usleep(1);
  //EXPECT_EQ(ObjectState::Materializing, z.state());

  EXPECT_DOUBLE_EQ(9, z(2));
  EXPECT_EQ(ObjectState::SolidState, z.state());

}

TEST(Vector, Sub)
{
  Vector x{2,4,6,8,10}, y{1,2,3,4,5};

  Vector z = x - y;

  EXPECT_TRUE(z == y);
}

TEST(Vector, DivScale)
{
  Vector x{2,4,6,8,10};
  Scalar s{2};

  Vector y = x / s;

  EXPECT_DOUBLE_EQ(1, y(0));
  EXPECT_DOUBLE_EQ(2, y(1));
  EXPECT_DOUBLE_EQ(3, y(2));
  EXPECT_DOUBLE_EQ(4, y(3));
  EXPECT_DOUBLE_EQ(5, y(4));
}

TEST(Vector, MulScale)
{
  Vector x{2,4,6,8,10};
  Scalar s{2};

  Vector y = x * s;

  Vector xx{4,8,12,16,20};

  EXPECT_TRUE(y == xx);
}

TEST(Vector, Dot)
{
  Vector x{1,3,5,7,9}, y{0,2,4,6,8};

  Scalar z = x * y;

  EXPECT_EQ(1*0+2*3+4*5+6*7+8*9, z()); 
}

//Scalar Tests ----------------------------------------------------------------

TEST(Scalar, Add)
{
  Scalar a{2}, b{2};

  Scalar c = a + b;

  EXPECT_DOUBLE_EQ(4, c());
}

TEST(Scalar, Sub)
{
  Scalar a{4}, b{2};
  
  Scalar c = a - b;

  EXPECT_DOUBLE_EQ(2, c());
}

TEST(Scalar, Mul)
{
  Scalar a{4}, b{2};
  
  Scalar c = a * b;

  EXPECT_DOUBLE_EQ(8, c());
}

TEST(Scalar, Div)
{
  Scalar a{5}, b{2};
  
  Scalar c = a / b;

  EXPECT_DOUBLE_EQ(2.5, c());
}

//Matrix Tests ----------------------------------------------------------------

TEST(Matrix, MulVec)
{
  Matrix A = Matrix::Identity(5,5);
  Vector x{1,2,3,4,5};

  Vector Ax = A * x;

  EXPECT_TRUE(Ax == x);
}
