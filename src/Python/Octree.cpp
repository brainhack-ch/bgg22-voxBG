//
// Created by guibertf on 12/2/22.
//

#include "Octree.h"

Octree::Octree(float length, unsigned int depth, Eigen::Vector3f center, Eigen::MatrixX3f vertices,
               Eigen::MatrixX3i triangles): n_triangles(triangles.rows()), depth(depth), side_length(length),
                                            root_center(center), vertices(vertices), triangles(triangles),
                                            is_root_leaf(depth == 0) {
}

unsigned int Octree::triangleNumbers() {
    return n_triangles;
}

Eigen::MatrixX3i Octree::getTriangles() {
    return triangles;
}

bool Octree::isRootLeaf() {
    return is_root_leaf;
}

unsigned int Octree::getNodeNumber() {
    return 0;
}

Octree::~Octree() = default;
