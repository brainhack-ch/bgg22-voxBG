//
// Created by guibertf on 12/2/22.
//

#ifndef GRAPH_ANALYSIS_CONSTRUCTIONNODE_H
#define GRAPH_ANALYSIS_CONSTRUCTIONNODE_H


#include "Eigen/Core"

class ConstructionNode {
public:
    std::vector<int> triangle_ids;
    Eigen::Vector3f center;
    float side_length;
    unsigned int leafid;
    bool isLeaf;

    ConstructionNode(Eigen::Vector3f parent_center, unsigned int node_space_id, float side_length_parent);
    void insertTriangles(const Eigen::MatrixX3f &vertex_coordinates,
                         const Eigen::MatrixX3i &triangle_vertices,
                         const std::vector<int> &integers_of_interest);

    /*void setFirstLeafChild(int leafid);

    void getFirstLeafChild();*/

    ~ConstructionNode();

    int getNumberChildren();

    void setAsNonLeaf();

    bool checkTriangleInBox(const Eigen::Vector3f &v1, const Eigen::Vector3f &v2, const Eigen::Vector3f &v3);


    /**
     * Tests if an edge intersects with a triangle.
     * @param edge_origin
     * @param edge_end
     * @return
     */
    bool
    checkTrianglesIntersect(Eigen::Vector3f edge_origin, Eigen::Vector3f edge_end, const Eigen::MatrixX3f &vertices,
                            const Eigen::MatrixX3i &triangles);

    /**
     * Tests if an edge intersects with a square box.
     * @param edge_origin
     * @param edge_end
     * @param box_center
     * @param box_side_length
     * @return
     */
    bool checkBoxIntersect(Eigen::Vector3f edge_origin, Eigen::Vector3f edge_end);

};


#endif //GRAPH_ANALYSIS_CONSTRUCTIONNODE_H
