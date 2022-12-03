//
// Created by guibertf on 12/2/22.
//

#include "Python/Octree.h"
#include <gtest/gtest.h>
#include "gmock/gmock.h"
#include "utils.h"

using ::testing::Return;

class OctreeTest : public ::testing::Test {
protected:
    OctreeTest() {

    }

    ~OctreeTest() override {

    }

    void SetUp() override {

    }

    void TearDown() override {

    }
};

/**
 * CONSTRUCTOR TESTS
 */

TEST_F(OctreeTest, initializingWithSingleTriangleAllowsToRecoverTriangle) {
    Eigen::Matrix3Xf vertices(3, 3);
    vertices << 1.0, 0.0, 0.0,
                 0., 1, 0,
                 0, 0, 1;
    Eigen::MatrixX3i triangles(1, 3);
    triangles << 0,
                 1,
                 2;

    Eigen::Vector3f center;
    center << 0.0, 0.0, 0.0;

    Octree octree(10, 3, center, vertices, triangles, 0);
    EXPECT_EQ(octree.triangleNumbers(), 1);
    EXPECT_TRUE(octree.getTriangles().isApprox(triangles));
}


TEST_F(OctreeTest, initializingWithSeveralTrianglesAllowsToRecoverTriangles) {
    Eigen::MatrixX3f vertices(4, 3);
    vertices << 1.0, 0.0, 0.0,
            0., 1, 0,
            0, 0, 1,
            0, 0, 0;
    Eigen::MatrixX3i triangles(3, 3);
    triangles << 0, 0, 1,
            1, 2, 0,
            2, 3, 3;

    Eigen::Vector3f center;
    center << 0.0, 0.0, 0.0;

    Octree octree(10, 3, center, vertices, triangles, 0);
    EXPECT_EQ(octree.triangleNumbers(), 3);
    EXPECT_TRUE(octree.getTriangles().isApprox(triangles));
}


TEST_F(OctreeTest, zeroDepthTreeConsidersRootNodeAsLeaf) {
    Eigen::MatrixX3f vertices(4, 3);
    vertices << 1.0, 0.0, 0.0,
            0., 1, 0,
            0, 0, 1,
            0, 0, 0;
    Eigen::MatrixX3i triangles(3, 3);
    triangles << 0, 0, 1,
            1, 2, 0,
            2, 3, 3;

    Eigen::Vector3f center;
    center << 0.0, 0.0, 0.0;

    Octree octree(10, 4, center, vertices, triangles, 0);
    EXPECT_TRUE(octree.isRootLeaf());
}


TEST_F(OctreeTest, zeroDepthTreeHasNoNodes) {
    Eigen::MatrixX3f vertices(4, 3);
    vertices << 1.0, 0.0, 0.0,
            0., 1, 0,
            0, 0, 1,
            0, 0, 0;
    Eigen::MatrixX3i triangles(3, 3);
    triangles << 0, 0, 1,
            1, 2, 0,
            2, 3, 3;

    Eigen::Vector3f center;
    center << 0.0, 0.0, 0.0;

    Octree octree(10, 4, center, vertices, triangles, 0);
    EXPECT_EQ(octree.getNodeNumber(), 0);
}

TEST_F(OctreeTest, oneDepthTreeProducesEightChildrenNodes) {
    Eigen::MatrixX3f vertices(4, 3);
    vertices << 1.0, 0.0, 0.0,
            0., 1, 0,
            0, 0, 1,
            0, 0, 0;
    Eigen::MatrixX3i triangles(3, 3);
    triangles << 0, 1, 2,
            0, 2, 3,
            2, 1, 3;

    Eigen::Vector3f center;
    center << 0.5, 0.5, 0.5;

    Octree octree(2, 2, center, vertices, triangles, 1);
    EXPECT_FALSE(octree.isRootLeaf());
    // The root does not count as a node
    EXPECT_EQ(octree.getNodeNumber(), 8);
}

TEST_F(OctreeTest, splittingOnRealCaseTrianglesWorksAsExpected) {
    std::vector<Eigen::Vector3f> vertices;
    std::vector<Eigen::Vector3i> faces;

    read_and_populate_from_obj_file(&vertices, &faces,
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
    /*std::vector<int> indices_of_interest;
    for(int i=0; i < faces.size();++i){
        indices_of_interest.push_back(i);
    }*/


    Octree octree(2, 4, center, mat_vertices, triangle_indices, 8);
    EXPECT_FALSE(octree.isRootLeaf());
    // The root does not count as a node
    EXPECT_EQ(octree.getNodeNumber(), 8);
}



/**
* INTERSECTION FUNCTIONS TESTS
*/


/**
* TREE TRAVERSAL TESTS
*/

