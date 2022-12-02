//
// Created by guibertf on 12/2/22.
//

#include "Python/Octree.h"
#include <gtest/gtest.h>
#include "gmock/gmock.h"

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

    Octree octree(10, 3, center, vertices, triangles);
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

    Octree octree(10, 3, center, vertices, triangles);
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

    Octree octree(10, 0, center, vertices, triangles);
    EXPECT_TRUE(octree.isRootLeaf());
}

TEST_F(OctreeTest, oneDepthTreeProducesEightChildrenNodes) {
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

    Octree octree(10, 1, center, vertices, triangles);
    EXPECT_FALSE(octree.isRootLeaf());
}

TEST_F(OctreeTest, oneDepthTreeStoresTrianglesInChildren) {
    EXPECT_FALSE(true);
}


TEST_F(OctreeTest, twoDepthTreeSubdividesChildrenProperly) {
    EXPECT_FALSE(true);
}


/**
* INTERSECTION FUNCTIONS TESTS
*/


/**
* TREE TRAVERSAL TESTS
*/

