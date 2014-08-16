/******************************************************************************
 * The Sven Project - All Rights Reserved
 *
 * 16 August 2014
 * ~ ry
 *
 ******************************************************************************/

#include "sven/Krylov.hxx"
#include "gtest/gtest.h"

using namespace sven;

TEST(Arnoldi, Small)
{
  SparseMatrix A(5,5,4,
      {2,        3,            3,           4,                 2      },
      {0,  1,    0,   3,  2,   3,  2,  1,   2,  1,   4,  3,    3,  4  },
      {0.4,0.24, 0.74,0.4,0.3, 0.5,0.9,0.7, 0.5,0.83,0.7,0.65, 0.7,0.7});

  Vector b{0.47, 0.32, 0.34, 0.41, 0.28};
  Vector x0{0,0,0,0,0};

  Arnoldi arnoldi(5, A, x0, b);
  arnoldi();
}

