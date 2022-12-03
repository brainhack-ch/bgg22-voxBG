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
    /*std::vector<int> indices_of_interest;
    for(int i=0; i < faces.size();++i){
        indices_of_interest.push_back(i);
    }*/


    Octree octree(2, 4, center, mat_vertices, triangle_indices, 8);
    EXPECT_FALSE(octree.isRootLeaf());
    // The root does not count as a node
    EXPECT_EQ(octree.getNodeNumber(), 8);
}


TEST_F(OctreeTest, intersectionOfEdgesWithBoxReturnsExpectedIntersects) {
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

    // Now load the edges
    std::vector<Eigen::Vector3f> edge_vertices;
    std::vector<Eigen::Vector2i> edge_assignments;
    read_and_populate_edges_from_obj_file(&edge_vertices, &edge_assignments,
                                          "/home/guibertf/CLionProjects/graph_analysis/bgg22-voxBG/test/Python/edges.obj");

    Eigen::MatrixX3f mat_vertices_edges(edge_vertices.size(), 3);
    Eigen::MatrixX2i edges_indices(edge_assignments.size(), 2);

    for(int i=0;  i < edge_vertices.size(); ++i){
        mat_vertices_edges.row(i) = edge_vertices.at(i);
    }

    for(int i=0;  i < edge_assignments.size(); ++i){
        edges_indices.row(i) = edge_assignments.at(i);
    }


    // Now let's create a sub node. From Blender, we can have some expectation as to how many triangles will fall within!
    Eigen::Vector3f center;
    center << 0, 0, 0;
    float parent_length = 2.0f;

    Octree octree(parent_length, 4, center, mat_vertices, triangle_indices, 0);
    EXPECT_TRUE(octree.isRootLeaf());
    // The root does not count as a node
    EXPECT_EQ(octree.getNodeNumber(), 0);

    bool expected_intersections[6] = {false, false, false, true, true, true};

    for(int i=0; i < edges_indices.rows(); ++i){
        Eigen::Vector3f edge_origin = mat_vertices_edges.row(edges_indices(i, 0));
        Eigen::Vector3f edge_end = mat_vertices_edges.row(edges_indices(i, 1));
        EXPECT_EQ(octree.isEdgeIntersecting(edge_origin, edge_end), expected_intersections[i]);
    }
}

TEST_F(OctreeTest, intersectionOfEdgesWithSubdividedTreeWorks) {
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

    // Now load the edges
    std::vector<Eigen::Vector3f> edge_vertices;
    std::vector<Eigen::Vector2i> edge_assignments;
    read_and_populate_edges_from_obj_file(&edge_vertices, &edge_assignments,
                                          "/home/guibertf/CLionProjects/graph_analysis/bgg22-voxBG/test/Python/edges.obj");

    Eigen::MatrixX3f mat_vertices_edges(edge_vertices.size(), 3);
    Eigen::MatrixX2i edges_indices(edge_assignments.size(), 2);

    for(int i=0;  i < edge_vertices.size(); ++i){
        mat_vertices_edges.row(i) = edge_vertices.at(i);
    }

    for(int i=0;  i < edge_assignments.size(); ++i){
        edges_indices.row(i) = edge_assignments.at(i);
    }


    // Now let's create a sub node. From Blender, we can have some expectation as to how many triangles will fall within!
    Eigen::Vector3f center;
    center << 0, 0, 0;
    float parent_length = 2.0f;

    Octree octree(parent_length, 4, center, mat_vertices, triangle_indices, 8);
    EXPECT_FALSE(octree.isRootLeaf());
    // The root does not count as a node
    EXPECT_EQ(octree.getNodeNumber(), 8);

    bool expected_intersections[6] = {false, false, false, true, true, true};

    for(int i=0; i < edges_indices.rows(); ++i){
        Eigen::Vector3f edge_origin = mat_vertices_edges.row(edges_indices(i, 0));
        Eigen::Vector3f edge_end = mat_vertices_edges.row(edges_indices(i, 1));
        EXPECT_EQ(octree.isEdgeIntersecting(edge_origin, edge_end), expected_intersections[i]);
    }
}

TEST_F(OctreeTest, singleEdgeOutsideOfBBDoesNotIntersect) {
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

    // Now load the edges
    std::vector<Eigen::Vector3f> edge_vertices;
    std::vector<Eigen::Vector2i> edge_assignments;
    read_and_populate_edges_from_obj_file(&edge_vertices, &edge_assignments,
                                          "/home/guibertf/CLionProjects/graph_analysis/bgg22-voxBG/test/Python/single_edge_out.obj");

    Eigen::MatrixX3f mat_vertices_edges(edge_vertices.size(), 3);
    Eigen::MatrixX2i edges_indices(edge_assignments.size(), 2);

    for(int i=0;  i < edge_vertices.size(); ++i){
        mat_vertices_edges.row(i) = edge_vertices.at(i);
    }

    for(int i=0;  i < edge_assignments.size(); ++i){
        edges_indices.row(i) = edge_assignments.at(i);
    }


    // Now let's create a sub node. From Blender, we can have some expectation as to how many triangles will fall within!
    Eigen::Vector3f center;
    center << 0, 0, 0;
    float parent_length = 2.0f;

    Octree octree(parent_length, 4, center, mat_vertices, triangle_indices, 0);
    EXPECT_TRUE(octree.isRootLeaf());
    // The root does not count as a node
    EXPECT_EQ(octree.getNodeNumber(), 0);

    for(int i=0; i < edges_indices.rows(); ++i){
        Eigen::Vector3f edge_origin = mat_vertices_edges.row(edges_indices(i, 0));
        Eigen::Vector3f edge_end = mat_vertices_edges.row(edges_indices(i, 1));
        EXPECT_EQ(false, octree.isEdgeIntersecting(edge_origin, edge_end));
    }
}

TEST_F(OctreeTest, singleEdgeInsideBBAndTriangleIntersects) {
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

    // Now load the edges
    std::vector<Eigen::Vector3f> edge_vertices;
    std::vector<Eigen::Vector2i> edge_assignments;
    read_and_populate_edges_from_obj_file(&edge_vertices, &edge_assignments,
                                          "/home/guibertf/CLionProjects/graph_analysis/bgg22-voxBG/test/Python/single_edge_in_triangle.obj");

    Eigen::MatrixX3f mat_vertices_edges(edge_vertices.size(), 3);
    Eigen::MatrixX2i edges_indices(edge_assignments.size(), 2);

    for(int i=0;  i < edge_vertices.size(); ++i){
        mat_vertices_edges.row(i) = edge_vertices.at(i);
    }

    for(int i=0;  i < edge_assignments.size(); ++i){
        edges_indices.row(i) = edge_assignments.at(i);
    }


    // Now let's create a sub node. From Blender, we can have some expectation as to how many triangles will fall within!
    Eigen::Vector3f center;
    center << 0, 0, 0;
    float parent_length = 2.0f;

    Octree octree(parent_length, 4, center, mat_vertices, triangle_indices, 0);
    EXPECT_TRUE(octree.isRootLeaf());
    // The root does not count as a node
    EXPECT_EQ(octree.getNodeNumber(), 0);

    for(int i=0; i < edges_indices.rows(); ++i){
        Eigen::Vector3f edge_origin = mat_vertices_edges.row(edges_indices(i, 0));
        Eigen::Vector3f edge_end = mat_vertices_edges.row(edges_indices(i, 1));
        EXPECT_EQ(true, octree.isEdgeIntersecting(edge_origin, edge_end));

    }
}

/**
* INTERSECTION FUNCTIONS TESTS
*/


/**
* TREE TRAVERSAL TESTS
*/

