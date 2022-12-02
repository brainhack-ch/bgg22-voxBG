//
// Created by guibertf on 12/2/22.
//

#include "ConstructionNode.h"

ConstructionNode::ConstructionNode(Eigen::Vector3f parent_center, unsigned int node_space_id,
                                   unsigned int side_length_parent): isLeaf(true) {
    // Decide the center here

}

ConstructionNode::~ConstructionNode() {

}


void ConstructionNode::insertTriangles(const Eigen::MatrixX3f &vertex_coordinates,
                                       const Eigen::MatrixX3i &triangle_vertices,
                                       const std::vector<int> &integers_of_interest) {
    // For every triangle
    for(auto &t: integers_of_interest){
        // Recover vertex coordinates of the triangle
        auto triangle_v = triangle_vertices.row(t);
        auto v1 = vertex_coordinates.row(triangle_v(0));
        auto v2 = vertex_coordinates.row(triangle_v(1));
        auto v3 = vertex_coordinates.row(triangle_v(2));

        // Check if it is in the box. If so, insert triangle id in triangles
        if(checkTriangleInBox(v1, v2, v3, center, side_length)){
            triangle_ids.push_back(t);
        }
    }
}

bool ConstructionNode::checkTriangleInBox(Eigen::Vector3f v1, Eigen::Vector3f v2, Eigen::Vector3f v3,
                                          Eigen::Vector3f box_center, float box_side_length) {
    return false;
}

int ConstructionNode::getNumberChildren() {
    return triangle_ids.size();
}

void ConstructionNode::setAsNonLeaf() {
    isLeaf = false;
}
