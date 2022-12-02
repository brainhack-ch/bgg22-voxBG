//
// Created by guibertf on 12/2/22.
//

#ifndef GRAPH_ANALYSIS_OCTREE_H
#define GRAPH_ANALYSIS_OCTREE_H
#include <Eigen/Dense>

struct Node {
    bool is_leaf;
    unsigned int first_child_id;
    unsigned int n_children;
};


class Octree {
private:
    unsigned int n_triangles;
    unsigned int max_triangles_per_leaf;
    float side_length;
    Eigen::Vector3f root_center;
    std::vector<Node*> nodes;
    std::vector<unsigned int> triangle_node_ids; // This array stores contiguously for each terminal leaf the indices of its triangles
    Eigen::MatrixX3f vertices; // Triangle vertices
    Eigen::MatrixX3i triangles; // Triangles of the mesh

    bool is_root_leaf;

public:

    /**
     * Constructor of the octree.
     * For a given depth (which defines space subdivision level), the triangles are to be inserted into the octree.
     * Specifically, a triangle is inserted into the octree node with which its barycenter intersects.
     * If several such nodes exist, it is inserted in all of them
     * @param length
     * @param max_triangles_per_leaf The number of triangles in a leaf before subdivision of a leaf will start
     * @param center
     * @param vertices Vertices of the mesh
     * @param triangles Triangles, defined as triplets of vertices ids (points towards vertices array)
     */
    Octree(float length, unsigned int max_triangles_per_leaf, Eigen::Vector3f center, Eigen::MatrixX3f vertices,
           Eigen::MatrixX3i triangles);

    /**
     * Main method of interest for the octree: the intersection checking of an edge with the octree
     * This method performs an iterative BFS traversal of the octree, checking if the edge is intersecting the
     * non leaf nodes. When it reaches a leaf node, it iterates on all its children (triangles) to check for intersection with the edge.
     * Returns false if and only if the edge has no intersection with any triangle in the octree.
     * @param edge_origin
     * @param edge_end
     * @return
     */
    bool isEdgeIntersecting(Eigen::Vector3f edge_origin, Eigen::Vector3f edge_end);

    unsigned int triangleNumbers();

    Eigen::MatrixX3i getTriangles();

    bool isRootLeaf();

    unsigned int getNodeNumber();

    ~Octree();

    /**
     * Tests if an edge intersects with a triangle.
     * @param edge_origin
     * @param edge_end
     * @return
     */
    static bool checkTriangleIntersect(Eigen::Vector3f edge_origin, Eigen::Vector3f edge_end);

    /**
     * Tests if an edge intersects with a square box.
     * @param edge_origin
     * @param edge_end
     * @param box_center
     * @param box_side_length
     * @return
     */
    static bool checkBoxIntersect(Eigen::Vector3f edge_origin, Eigen::Vector3f edge_end, Eigen::Vector3f box_center, float box_side_length);

};

#endif //GRAPH_ANALYSIS_OCTREE_H
