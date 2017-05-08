#include "index.h"
#include "gtest/gtest.h"

TEST(IndexTest, StartEndAndToString) {
  Index index;
  EXPECT_TRUE(index.to_string() == "");
  index.start();
  index.start();
  EXPECT_TRUE(index.to_string() == "");
  index.start();
  EXPECT_TRUE(index.to_string() == "1. ");
  index.start();
  EXPECT_TRUE(index.to_string() == "1.1. ");
  index.end();
  index.start();
  EXPECT_TRUE(index.to_string() == "1.2. ");
}