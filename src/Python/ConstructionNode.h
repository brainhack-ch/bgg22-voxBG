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
    unsigned int side_length;
    unsigned int leafid;
    bool isLeaf;
    ConstructionNode(Eigen::Vector3f parent_center, unsigned int node_space_id, unsigned int side_length_parent);

    void insertTriangles(const Eigen::MatrixX3f &vertex_coordinates,
                         const Eigen::MatrixX3i &triangle_vertices,
                         const std::vector<int> &integers_of_interest);

    void setFirstLeafChild(int leafid);

    void getFirstLeafChild();

    ~ConstructionNode();

    int getNumberChildren();

    void setAsNonLeaf();

    static bool checkTriangleInBox(Eigen::Vector3f v1, Eigen::Vector3f v2, Eigen::Vector3f v3, Eigen::Vector3f box_center, float box_side_length);

};


#endif //GRAPH_ANALYSIS_CONSTRUCTIONNODE_H
