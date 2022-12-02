//
// Created by guibertf on 12/2/22.
//

#include "ConstructionNode.h"

ConstructionNode::ConstructionNode(Eigen::Vector3f parent_center, unsigned int node_space_id,
                                   float side_length_parent): isLeaf(true), leafid(0),
                                   side_length(side_length_parent/2) {
    // Decide the center here
    // We must construct the X, Y and Z displacement
    float half_side = side_length*0.5;
    float x_shift = node_space_id % 4 < 2 ? -half_side: half_side;
    float y_shift = node_space_id % 2 == 0 ? -half_side : half_side;
    float z_shift = node_space_id >=4 ? half_side : -half_side;
    center(0) = parent_center.x() + x_shift;
    center(1) = parent_center.y() + y_shift;
    center(2) = parent_center.z() + z_shift;
    //disp << x_shift, y_shift , z_shift;
}

ConstructionNode::~ConstructionNode() = default;


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

bool ConstructionNode::checkTriangleInBox(const Eigen::Vector3f &v1, const Eigen::Vector3f &v2, const Eigen::Vector3f &v3,
                                          const Eigen::Vector3f &box_center, float box_side_length) {
    return false;
}

int ConstructionNode::getNumberChildren() {
    return triangle_ids.size();
}

void ConstructionNode::setAsNonLeaf() {
    isLeaf = false;
}

bool ConstructionNode::checkTriangleIntersect(Eigen::Vector3f edge_origin, Eigen::Vector3f edge_end) {
    return false;
}

bool
ConstructionNode::checkBoxIntersect(Eigen::Vector3f edge_origin, Eigen::Vector3f edge_end, Eigen::Vector3f box_center,
                                    float box_side_length) {
    return false;
}
