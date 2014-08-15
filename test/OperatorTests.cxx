#include "sven/CoreMath.hxx"
#include "gtest/gtest.h"

using namespace sven;

TEST(VectorVector, SubscriptAndSize)
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

TEST(VectorVector, Add)
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

TEST(VectorVector, Dot)
{
  Vector x{1,3,5,7,9}, y{0,2,4,6,8};

  Scalar z = x * y;

  EXPECT_EQ(1*0+2*3+4*5+6*7+8*9, z()); 
}
