//
// Created by guibertf on 12/2/22.
//

#include "ConstructionNode.h"
#include <Eigen/Geometry>
#include <math.h>
#include <iostream>

ConstructionNode::ConstructionNode(Eigen::Vector3f parent_center, unsigned int node_space_id,
                                   float side_length_parent): isLeaf(true), leafid(0),
                                   side_length(side_length_parent/2) {
    if(node_space_id >= 8){
        throw std::invalid_argument("Node space ID can only be between 0 and 7 for an octree.");
    }
    // Decide the center here
    // We must construct the X, Y and Z displacement
    if(abs(side_length_parent) < 10e-14 || side_length_parent < 0){
        throw std::invalid_argument("Side length cannot be null");
    }
    float half_side = side_length*0.5;
    float x_shift = node_space_id % 4 < 2 ? -half_side: half_side;
    float y_shift = node_space_id % 2 == 0 ? -half_side : half_side;
    float z_shift = node_space_id >=4 ? half_side : -half_side;
    center(0) = parent_center.x() + x_shift;
    center(1) = parent_center.y() + y_shift;
    center(2) = parent_center.z() + z_shift;
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
        if(checkTriangleInBox(v1, v2, v3)){
            triangle_ids.push_back(t);
        }
    }
}

bool ConstructionNode::checkTriangleInBox(const Eigen::Vector3f &v1, const Eigen::Vector3f &v2, const Eigen::Vector3f &v3) {
    /*
     * Credits for the base code (adapted to our API) based on Separating axis theorem:
     * https://gdbooks.gitbooks.io/3dcollisions/content/Chapter4/aabb-triangle.html
     * Based on: Akenine-Möller, Tomas. "Fast 3D triangle-box overlap testing." Acm siggraph 2005 courses. 2005. 8-es.
     */

    // Translate points to the origin; conceptually, we want to detect if a collision occurs with a bounding box centered at origin here
    Eigen::Vector3f v1_c = v1 - center;
    Eigen::Vector3f v2_c = v2 - center;
    Eigen::Vector3f v3_c = v3 - center;

    // Compute edge vectors of the triangle (needed for barycentric coordinates)
    Eigen::Vector3f f0 = v2_c - v1_c;
    Eigen::Vector3f f1 = v3_c - v2_c;
    Eigen::Vector3f f2 = v1_c - v3_c;

    // Face normals of the bounding box (trivially the X, Y, Z axis)
    Eigen::Vector3f u0(1,0,0);
    Eigen::Vector3f u1(0,1,0);
    Eigen::Vector3f u2(0,0,1);

    Eigen::MatrixX3f separating_axis_faces_bb(13, 3);

    separating_axis_faces_bb.row(0) = u0.cross(f0);
    separating_axis_faces_bb.row(1) = u0.cross(f1);
    separating_axis_faces_bb.row(2) = u0.cross(f2);

    separating_axis_faces_bb.row(3) = u1.cross(f0);
    separating_axis_faces_bb.row(4) = u1.cross(f1);
    separating_axis_faces_bb.row(5) = u1.cross(f2);

    separating_axis_faces_bb.row(6) = u2.cross(f0);
    separating_axis_faces_bb.row(7) = u2.cross(f1);
    separating_axis_faces_bb.row(8) = u2.cross(f2);

    // Check if triangle extents are aligned with bounding box along face normals of bounding box
    separating_axis_faces_bb.row(9) = u0;
    separating_axis_faces_bb.row(10) = u1;
    separating_axis_faces_bb.row(11) = u2;

    // Check if triangle extents are aligned with bounding box along triangle normal
    separating_axis_faces_bb.row(12) = f0.cross(f1);

    // Test on all face separating axis
    for(int axis_i = 0; axis_i < 13; axis_i++ ) {
        Eigen::Vector3f sep_axis = separating_axis_faces_bb.row(axis_i);

        // Project all 3 vertices onto separating axis
        float p0 = v1_c.dot(sep_axis);
        float p1 = v2_c.dot(sep_axis);
        float p2 = v3_c.dot(sep_axis);

        // Project bounding box along separating axis
        // We only care about the extents of the box to project it along the axis!
        float extent_length_x= side_length/2;
        float extent_length_y= side_length/2;
        float extent_length_z= side_length/2;

        float r = extent_length_x*std::abs(u0.dot(sep_axis)) +
                extent_length_y*std::abs(u1.dot(sep_axis)) +
                extent_length_z*std::abs(u2.dot(sep_axis));

        // See if most extreme point of the triangle intersects r.
        // If we find an axis where this is not the case, there is a separating axis: intersection impossible.
        if(std::max(-std::max({p0, p1, p2}), std::min({p0, p1, p2})) > r) {
            return false;
        }
    }
    // If all tests passed, we can return true: there exist no separating axis!
    return true;
}

int ConstructionNode::getNumberChildren() {
    return triangle_ids.size();
}

void ConstructionNode::setAsNonLeaf() {
    isLeaf = false;
}

bool ConstructionNode::checkTrianglesIntersect(Eigen::Vector3f edge_origin, Eigen::Vector3f edge_end,
                                               const Eigen::MatrixX3f &vertices,
                                               const Eigen::MatrixX3i &triangles) {
    if(isLeaf){
        for(auto &t: triangle_ids){
            // Loop over triangles to figure out if there exists an intersection
            std::vector<int> ind(triangles.row(t).begin(), triangles.row(t).end());
            // Check with this triangle's vertices coordinates if there is indeed an intersection
            if(checkTriangleIntersect(edge_origin, edge_end, vertices(ind, Eigen::placeholders::all))){
                return true;
            }
        }
    } else {
        throw std::logic_error("Cannot check for intersection on non-leaf node!");
    }
    return false;
}

bool ConstructionNode::boxToEdgeIntersection(Eigen::Vector3f edge_origin, Eigen::Vector3f edge_end, Eigen::Vector3f box_center, float box_side_length){
    // Recenter edge vector first around center
    Eigen::Vector3f p1_c = edge_origin - box_center;
    Eigen::Vector3f p2_c = edge_end - box_center;
    Eigen::Vector3f f0 = p2_c - p1_c;

    // Let's apply again SAT

    // Face normals of the bounding box (trivially the X, Y, Z axis)
    Eigen::Vector3f u0(1,0,0);
    Eigen::Vector3f u1(0,1,0);
    Eigen::Vector3f u2(0,0,1);

    Eigen::MatrixX3f separating_axis_faces_bb(6, 3);

    separating_axis_faces_bb.row(0) = u0.cross(f0);
    separating_axis_faces_bb.row(1) = u1.cross(f0);
    separating_axis_faces_bb.row(2) = u2.cross(f0);

    // Check if vector extents are aligned with bounding box along face normals of bounding box
    separating_axis_faces_bb.row(3) = u0;
    separating_axis_faces_bb.row(4) = u1;
    separating_axis_faces_bb.row(5) = u2;

    float extent_length = box_side_length/2;


    // Test on all face separating axis
    for(int axis_i = 0; axis_i < 6; axis_i++ ) {
        Eigen::Vector3f sep_axis = separating_axis_faces_bb.row(axis_i);

        // Project all 2 vertices onto separating axis
        float p0 = p1_c.dot(sep_axis);
        float p1 = p2_c.dot(sep_axis);

        // Project bounding box along separating axis
        // We only care about the extents of the box to project it along the axis!

        // Note: here we consider the extent length as uniform across all axis. Because it can change, we are explictly multiplying each axis separately to allow easier change later
        float r = (std::abs(extent_length* u0.dot(sep_axis)) + extent_length*std::abs(u1.dot(sep_axis)) + extent_length*std::abs(u2.dot(sep_axis)));

        // See if most extreme point of the triangle intersects r.
        // If we find an axis where this is not the case, there is a separating axis: intersection impossible.
        if(std::max(-std::max({p0, p1}), std::min({p0, p1})) > r) {
            return false;
        }
    }
    // If all tests passed, we can return true: there exist no separating axis!
    return true;
}

bool
ConstructionNode::checkBoxIntersect(Eigen::Vector3f edge_origin, Eigen::Vector3f edge_end) {
    return boxToEdgeIntersection(edge_origin, edge_end, center, side_length);
}

bool ConstructionNode::checkTriangleIntersect(Eigen::Vector3f edge_origin, Eigen::Vector3f edge_end,
                                              const Eigen::Matrix3f &triangle_vertices) {
    /**
     * Based on Möller, Tomas, and Ben Trumbore. "Fast, minimum storage ray/triangle intersection." ACM SIGGRAPH 2005 Courses. 2005. 7-es.
     */

    Eigen::Vector3f v0 = triangle_vertices.row(0);
    Eigen::Vector3f v1 = triangle_vertices.row(1);
    Eigen::Vector3f v2 = triangle_vertices.row(2);

    Eigen::Vector3f e1 = v1 - v0;
    Eigen::Vector3f e2 = v2 - v0;
    Eigen::Vector3f D = (edge_end-edge_origin).normalized();
    Eigen::Vector3f T = edge_origin - v0;

    Eigen::Vector3f q = T.cross(e1);
    Eigen::Vector3f P = D.cross(e2);
    float quotient = P.dot(e1);
    float t = q.dot(e2) / quotient;
    // The segment of intersection is greater than our own segment: intersection cannot happen with this edge!
    if(t > (edge_end-edge_origin).norm()){
        return false;
    }
    // Check that we are on the triangle
    float u = P.dot(T) /quotient;
    // We cannot be on the triangle!
    if(u <0 || u > 1){
        return false;
    }

    float v = q.dot(D) / quotient;
    if(v < 0 || v > 1 || u + v > 1){
        return false;
    }

    return true;
}
