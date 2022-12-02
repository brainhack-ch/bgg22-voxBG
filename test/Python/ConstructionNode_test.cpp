//
// Created by guibertf on 12/2/22.
//

#include "Python/ConstructionNode.h"
#include <gtest/gtest.h>
#include "gmock/gmock.h"

using ::testing::Return;

class ConstructionNodeTest : public ::testing::Test {
protected:
    ConstructionNodeTest() {

    }

    ~ConstructionNodeTest() override {

    }

    void SetUp() override {

    }

    void TearDown() override {

    }
};

/**
 * CONSTRUCTOR TESTS
 */

TEST_F(ConstructionNodeTest, initializationWithCenterGivesCorrectSubCenters) {
    Eigen::Vector3f center;
    center << 13, -2, 4;
    float side_length = 50;

    Eigen::MatrixX3f expected_centers(8, 3);
    float shift_size = side_length / 4;
    float x = center.x();
    float y = center.y();
    float z = center.z();
    expected_centers << x-shift_size, y-shift_size, z-shift_size,
                        x - shift_size,y+shift_size, z-shift_size,
                        x+shift_size,y-shift_size, z-shift_size,
                        x+shift_size, y+shift_size, z-shift_size,
                        x-shift_size, y-shift_size, z+shift_size,
                        x - shift_size, y+shift_size, z+shift_size,
                        x+shift_size, y-shift_size,  z+shift_size,
                        x+shift_size, y+shift_size, z+shift_size;
    
    for(int i=0; i < 8; ++i){
        ConstructionNode node(center, i, side_length);
        EXPECT_FLOAT_EQ(node.side_length, 25);
        EXPECT_TRUE(node.center.isApprox(expected_centers.row(i).transpose()));
        //std::cout << expected_centers.row(i) << std::endl;
        //std::cout << node.center << std::endl;
    }
}