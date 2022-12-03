//
// Created by guibertf on 12/2/22.
//

#include <queue>
#include <numeric>
#include "Octree.h"
#include <stack>
#include <iostream>

#if PYTHON_BIND == 1
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#endif

/*
 * TODO: ADD MAX SUBDIVISION LEVEL
 */

Octree::Octree(float length, unsigned int max_triangles_per_leaf, unsigned int max_subdivision_level): max_triangles_per_leaf(max_triangles_per_leaf),
side_length(length), max_depth(max_subdivision_level)
{

}

void Octree::create_octree() {
    if(max_triangles_per_leaf == 0){
        throw std::invalid_argument("Cannot have 0 triangles per node");
    }
    // The base case is simple: check if the number of inserted triangles is above max_triangles_per_leaf
    if(triangles.rows() <= max_triangles_per_leaf || max_depth == 0) {
        is_root_leaf = true;
        for(int i=0; i < n_triangles; ++i){
            triangle_node_ids.push_back(i);
        }
    } else {
        is_root_leaf = false;
        //std::vector<ConstructionNode*> construction_nodes;

        // Start by creating first 8 nodes:
        // - Allocate 8 nodes, each with 0 child
        // - Push back 8 new empty vectors in the tmp triangle_assignments vector
        std::vector<int> starting_values(n_triangles);
        unsigned int curr_depth = 0;
        std::iota(starting_values.begin(), starting_values.end(), 0);
        for(int i=0; i < 8; ++i){
            ConstructionNode* node = new ConstructionNode(root_center, i, side_length);
            nodes.push_back(node);
            node->insertTriangles(vertices, triangles, starting_values);
            node->division_level = 1;
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
            if(current_node->getNumberChildren() > max_triangles_per_leaf && current_node->division_level < max_depth) {
                current_node->setAsNonLeaf();

                // Extract the children and pass the triangle ids to them
                for(int i=0; i < 8; ++i){
                    auto* node = new ConstructionNode(current_node->center, i, current_node->side_length);
                    nodes.push_back(node);
                    node->insertTriangles(vertices, triangles, current_node->triangle_ids);
                    node->division_level = current_node->division_level + 1;
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

void Octree::init_octree(Eigen::Vector3f center, Eigen::MatrixX3f vertices,
                 Eigen::MatrixX3i triangles){
    root_center = center;
    n_triangles = triangles.rows();
    this->vertices = vertices;
    this->triangles = triangles;
    create_octree();
}

Octree::Octree(float length, unsigned int max_triangles_per_leaf, Eigen::Vector3f center, Eigen::MatrixX3f vertices,
               Eigen::MatrixX3i triangles, unsigned int max_subdivision_level)
        : n_triangles(triangles.rows()), max_triangles_per_leaf(max_triangles_per_leaf), side_length(length),
          root_center(center), vertices(vertices), triangles(triangles), max_depth(max_subdivision_level){
    create_octree();
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
    bool intersect = false;

    // Is the edge intersecting with the root ?
    if(ConstructionNode::boxToEdgeIntersection(edge_origin, edge_end, root_center, side_length)){
        std::stack<ConstructionNode*> node_stack;

        // Test against all triangles directly
        if(is_root_leaf){
            for(auto &t: triangle_node_ids){
                // Loop over triangles to figure out if there exists an intersection
                std::vector<int> ind(triangles.row(t).begin(), triangles.row(t).end());
                // Check with this triangle's vertices coordinates if there is indeed an intersection
                if(ConstructionNode::checkTriangleIntersect(edge_origin, edge_end, vertices(ind, Eigen::placeholders::all))){
                    intersect=true;
                    break;
                }
            }
        } else {
            // Test against children: perform BFS
            for(int i=0; i < 8; ++i){
                node_stack.push(nodes.at(i));
            }

            while(!node_stack.empty()  &&  !intersect){
                ConstructionNode *c_node = node_stack.top();
                node_stack.pop();

                if(c_node->isLeaf){
                    intersect = c_node->checkTrianglesIntersect(edge_origin, edge_end, vertices, triangles);
                } else {
                    if(c_node->checkBoxIntersect(edge_origin, edge_end) && c_node->leafid != 0){
                        for(int i=0; i < 8; ++i){
                            node_stack.push(nodes.at(i+c_node->leafid));
                        }
                    }
                }

            }
        }
    }
    // Perform DFS traversal of the nodes and query everytime
    return intersect;
}

#if PYTHON_BIND == 1

namespace py = pybind11;

PYBIND11_MODULE(graphPrunerLib, m) {
    py::class_<Octree>(m, "Octree").def(py::init<float, unsigned int, unsigned int>())
             .def("initialize", &Octree::init_octree)
             .def("checkEdgeIntersection", &Octree::isEdgeIntersecting);
}
#endif