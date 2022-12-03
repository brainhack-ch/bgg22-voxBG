//
// Created by guibertf on 12/2/22.
//

#include "Python/ConstructionNode.h"
#include "utils.h"
#include <gtest/gtest.h>
#include "gmock/gmock.h"
#include <iostream>
#include <string>
#include <fstream>

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

TEST_F(ConstructionNodeTest, initializerWithZeroSideLengthRaisesException){
    Eigen::Vector3f center;
    center << 13, -2, 4;
    EXPECT_THROW(ConstructionNode node(center, 0, 0.0f), std::invalid_argument);
}

TEST_F(ConstructionNodeTest, initializerWithNegativeSideLengthRaisesException){
    Eigen::Vector3f center;
    center << 13, -2, 4;
    EXPECT_THROW(ConstructionNode node(center, 0, -4), std::invalid_argument);
}

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
    }
}

TEST_F(ConstructionNodeTest, initForDifferentSideLengthsReturnCorrectSubcenters){
    Eigen::Vector3f center;
    center << 13, -2, 4;
    //float side_length = 50;
    float side_lengths[5] = {3, 5, 4, 10, 100};
    for(auto &side_length: side_lengths){
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
            EXPECT_FLOAT_EQ(node.side_length, side_length/2);
            EXPECT_TRUE(node.center.isApprox(expected_centers.row(i).transpose()));
        }
    }
}


/**
* INTERSECTION FUNCTIONS TESTS
*/

/**
 * EDGE TO TRIANGLE INTERSECTION TESTS
 */

TEST_F(ConstructionNodeTest, edgeNormalToTriangleWithCorrectLengthIntersects){
    Eigen::Vector3f center;
    center << 5, 5, 5;
    float parent_length = 20.0f;
    ConstructionNode node(center, 0, parent_length);

    Eigen::MatrixX3f vertices(3, 3);
    vertices << 1.0, 0.0, 0.0,
            0., 1, 0,
            0, 0, 1;
    Eigen::MatrixX3i triangles(1, 3);
    triangles << 0, 0, 1;

    //node.insertTriangles(vertices, triangles, std::vector({0}));

    Eigen::Vector3f p1(0.0, 0.0, 0.0);
    Eigen::Vector3f p2(0.0, 0.0, 1.0);

   EXPECT_TRUE(node.checkTriangleIntersect(p1, p2, vertices));
}

TEST_F(ConstructionNodeTest, whenEdgeJustTooSmallNoIntersect){
    Eigen::Vector3f center;
    center << 5, 5, 5;
    float parent_length = 20.0f;
    ConstructionNode node(center, 0, parent_length);

    Eigen::MatrixX3f vertices(3, 3);
    vertices << 1.0, 0.0, 0.0,
            0., 1, 0,
            0, 0, 1;
    Eigen::MatrixX3i triangles(1, 3);
    triangles << 0, 0, 1;

    //node.insertTriangles(vertices, triangles, std::vector({0}));

    Eigen::Vector3f p1(0.0, 0.0, 0.0);
    Eigen::Vector3f p2(0.0, 0.0, 0.9);

    EXPECT_FALSE(node.checkTriangleIntersect(p1, p2, vertices));
}

TEST_F(ConstructionNodeTest, edgeParallelButWithDifferentElevationNoIntersect){
    Eigen::Vector3f center;
    center << 5, 5, 5;
    float parent_length = 20.0f;
    ConstructionNode node(center, 0, parent_length);

    Eigen::MatrixX3f vertices(3, 3);
    vertices << 1.0, 0.0, 0.0,
            3., 10, 0,
            2, 4, 0;
    Eigen::MatrixX3i triangles(1, 3);
    triangles << 0, 0, 1;

    //node.insertTriangles(vertices, triangles, std::vector({0}));

    Eigen::Vector3f p1(0.5, 3.0, 2.0);
    Eigen::Vector3f p2(0.2, 5.0, 2.0);

    EXPECT_FALSE(node.checkTriangleIntersect(p1, p2, vertices));
}


/**
 * EDGE TO BOX INTERSECTION TESTS
 */

TEST_F(ConstructionNodeTest, edgeOutsideBoxDoesNotIntersect){
    Eigen::Vector3f center;
    center << 5, 5, 5;
    float parent_length = 20.0f;
    ConstructionNode node(center, 0, parent_length);

    // Case 1: edge is outside by being above in Z axis
    Eigen::Vector3f p1;
    p1 << 0.5, 10, 10;
    Eigen::Vector3f p2;
    p2 << 0.3, 10, 10;
    EXPECT_FALSE(node.checkBoxIntersect(p1, p2));
}

TEST_F(ConstructionNodeTest, edgeFullyInsideIntersects){
    Eigen::Vector3f center;
    center << 5, 5, 5;
    float parent_length = 20.0f;
    ConstructionNode node(center, 0, parent_length);

    Eigen::Vector3f p1;
    p1 << 0.5, 0, 0.5;
    Eigen::Vector3f p2;
    p2 << 0.3, 0.2, 4;
    EXPECT_TRUE(node.checkBoxIntersect(p1, p2));
}

TEST_F(ConstructionNodeTest, edgePartiallyInsideIntersects){
    Eigen::Vector3f center;
    center << 5, 5, 5;
    float parent_length = 20.0f;
    ConstructionNode node(center, 0, parent_length);

    Eigen::Vector3f p1;
    p1 << 0.5, 0, 0.5;
    Eigen::Vector3f p2;
    p2 << 0.3, 20, 4;
    EXPECT_TRUE(node.checkBoxIntersect(p1, p2));
}

/**
 * TRIANGLE TO BOX INTERSECTION TESTS
 */
TEST_F(ConstructionNodeTest, triangleOutsideBoundindBoxTest){
    Eigen::Vector3f center;
    center << 5, 5, 5;
    float parent_length = 20.0f;
    ConstructionNode node(center, 0, parent_length);

    // Case 1: triangle is outside by being above in Z axis
    Eigen::Vector3f p1;
    p1 << 0.5, 10, 10;
    Eigen::Vector3f p2;
    p2 << 0.3, 10, 10;
    Eigen::Vector3f p3;
    p3 << 0.5, 5, 10;

    std::cout << node.center.transpose() << std::endl;
    std::cout << node.side_length << std::endl;
    EXPECT_FALSE(node.checkTriangleInBox(p1, p2, p3));

    // Case 2: triangle is outside by being below in Z axis
    p1 << 0.5, 10, -10;
    p2 << 0.3, 10, -10;
    p3 << 0.5, 5, -10;

    EXPECT_FALSE(node.checkTriangleInBox(p1, p2, p3));
    // Case 3: triangle outside by being above in X axis
    p1 << 10, 4, 0.5;
    p2 << 10, 4, 0.5;
    p3 << 10, 5,  0.5;

    EXPECT_FALSE(node.checkTriangleInBox(p1, p2, p3));

    // Case 4: triangle outside by being below in X axis
    p1 << -10, 4, 0.5;
    p2 << -10, 4, 0.5;
    p3 << -10, 5,  0.5;

    EXPECT_FALSE(node.checkTriangleInBox(p1, p2, p3));


    // Case 5: triangle outside by being above in Y axis
    p1 << 3, 10, 0.5;
    p2 << 3, 10, 0.5;
    p3 << 3, 10,  0.5;

    EXPECT_FALSE(node.checkTriangleInBox(p1, p2, p3));


    // Case 6: triangle outside by being below in Y axis
    p1 << 3, -10, 0.5;
    p2 << 3, -10, 0.5;
    p3 << 3, -10,  0.5;

    EXPECT_FALSE(node.checkTriangleInBox(p1, p2, p3));
}

TEST_F(ConstructionNodeTest, triangleCompletelyContainedTest){
    Eigen::Vector3f center;
    center << 5, 5, 5;
    float parent_length = 20.0f;
    ConstructionNode node(center, 0, parent_length);

    Eigen::Vector3f p1;
    p1 << 0.5, -3, -2;
    Eigen::Vector3f p2;
    p2 << 3, 0.5, 2.6;
    Eigen::Vector3f p3;
    p3 << -4, 3.4, 0;
    EXPECT_TRUE(node.checkTriangleInBox(p1, p2, p3));
}


TEST_F(ConstructionNodeTest, trianglePartiallyContainedTest){
    Eigen::Vector3f center;
    center << 5, 5, 5;
    float parent_length = 20.0f;
    ConstructionNode node(center, 0, parent_length);

    Eigen::Vector3f p1;
    p1 << 0.5, -3, -2;
    Eigen::Vector3f p2;
    p2 << 20, 0.5, 2.6;
    Eigen::Vector3f p3;
    p3 << -4, 3.4, 30;
    EXPECT_TRUE(node.checkTriangleInBox(p1, p2, p3));
}

TEST_F(ConstructionNodeTest, triangleIntersectingButOutsideTest){
    Eigen::Vector3f center;
    center << 5, 5, 5;
    float parent_length = 20.0f;
    ConstructionNode node(center, 0, parent_length);

    Eigen::Vector3f p1;
    p1 << 0, -4, 0;
    Eigen::Vector3f p2;
    p2 << 20, 4, 2.6;
    Eigen::Vector3f p3;
    p3 << -4, 3.4, 30;
    EXPECT_TRUE(node.checkTriangleInBox(p1, p2, p3));
}



TEST_F(ConstructionNodeTest, detectionOfTriangleInCaseOfSeveralTriangles){
    std::vector<Eigen::Vector3f> vertices;
    std::vector<Eigen::Vector3i> faces;

    read_and_populate_faces_from_obj_file(&vertices, &faces,
                                          "/home/guibertf/CLionProjects/graph_analysis/bgg22-voxBG/test/Python/subcube_config.obj");

    // Now allocate our dear matrix!
    Eigen::MatrixX3f mat_vertices(vertices.size(), 3);
    Eigen::MatrixX3i triangle_indices(faces.size(), 3);

    for(int i=0;  i < vertices.size(); ++i){
        mat_vertices.row(i) = vertices.at(i);
    }

    for(int i=0;  i < faces.size(); ++i){
        triangle_indices.row(i) = faces.at(i);
    }

    // Now let's create a sub node. From Blender, we can have some expectation as to how many triangles will fall within!
    Eigen::Vector3f center;
    center << 0, 0, 0;
    float parent_length = 2.0f;
    std::vector<int> indices_of_interest;
    for(int i=0; i < faces.size();++i){
        indices_of_interest.push_back(i);
    }

    int expected_triangles_per_sub_cube[8] = {2,  2,  2, 4, 0, 0, 1, 0};

    for(int i=0; i < 8; ++i){
        ConstructionNode node(center, i, parent_length);
        node.insertTriangles(mat_vertices, triangle_indices, indices_of_interest);
        EXPECT_EQ(node.triangle_ids.size(), expected_triangles_per_sub_cube[i]);
    }
}