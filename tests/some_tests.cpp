#include <gtest/gtest.h>

class SomeTest : public ::testing::Test
{
};

TEST_F( SomeTest, aPassingTest )
{
    EXPECT_EQ( 2, 1+1 );
}

/*
TEST_F( SomeTest, aFailingTest )
{
    EXPECT_EQ( 3, 1+1 );
}
*/
