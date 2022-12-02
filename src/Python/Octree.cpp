//
// Created by guibertf on 12/2/22.
//

#include <queue>
#include <numeric>
#include "Octree.h"

Octree::Octree(float length, unsigned int max_triangles_per_leaf, Eigen::Vector3f center, Eigen::MatrixX3f vertices,
               Eigen::MatrixX3i triangles): n_triangles(triangles.rows()), max_triangles_per_leaf(max_triangles_per_leaf), side_length(length),
                                            root_center(center), vertices(vertices), triangles(triangles){

    if(max_triangles_per_leaf == 0){
        throw std::invalid_argument("Cannot have 0 triangles per node");
    }
    // The base case is simple: check if the number of inserted triangles is above max_triangles_per_leaf
    if(triangles.rows() < max_triangles_per_leaf) {
        is_root_leaf = true;
    } else {
        is_root_leaf = false;
        //std::vector<ConstructionNode*> construction_nodes;

        // Start by creating first 8 nodes:
        // - Allocate 8 nodes, each with 0 child
        // - Push back 8 new empty vectors in the tmp triangle_assignments vector
        std::vector<int> starting_values(n_triangles);
        std::iota(starting_values.begin(), starting_values.end(), 0);
        for(int i=0; i < 8; ++i){
            ConstructionNode* node = new ConstructionNode(center, i, side_length);
            nodes.push_back(node);
            node->insertTriangles(vertices, triangles, starting_values);
        }


        int current_i = 0;
        int last_i = nodes.size();
        while(current_i < last_i){
            ConstructionNode* current_node = nodes.at(current_i);
            // Subdivision case
            // Follows the steps:
            // - The parent is no longer a leaf (mark it as non leaf)
            // - Create 8 children and make the parent point to them
            // - Pass the triangles to the children
            // - Add children to overall tree structure
            // - Clear triangles from parent's triangles!
            if(current_node->getNumberChildren() > max_triangles_per_leaf) {
                current_node->setAsNonLeaf();

                // Extract the children and pass the triangle ids to them
                for(int i=0; i < 8; ++i){
                    auto* node = new ConstructionNode(current_node->center, i, current_node->side_length);
                    nodes.push_back(node);
                    node->insertTriangles(vertices, triangles, current_node->triangle_ids);
                }

                // Remove triangle ids
                current_node->triangle_ids.clear();

                // Set where this node's children begin: inserted at current vector's end
                current_node->leafid = last_i;

                // Update last_i
                last_i = nodes.size();
            }

            // Go to next node
            current_i += 1;

        }
    }

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
    return nodes.size();
}

Octree::~Octree(){
    for(int i=0; i < nodes.size(); ++i){
        delete nodes.at(i);
        nodes.at(i) = nullptr;
    }
    nodes.clear();

}

bool Octree::isEdgeIntersecting(Eigen::Vector3f edge_origin, Eigen::Vector3f edge_end) {
    // Perform DFS traversal of the nodes and query everytime
    return false;
}
