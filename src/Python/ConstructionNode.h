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
    int division_level;

    ConstructionNode(const Eigen::Ref<const Eigen::Vector3f> &parent_center, unsigned int node_space_id,
                     float side_length_parent);
    void insertTriangles(const Eigen::Ref<const Eigen::MatrixX3f> &vertex_coordinates,
                         const Eigen::Ref<const Eigen::MatrixX3i> &triangle_vertices,
                         const std::vector<int> &integers_of_interest);

    /*void setFirstLeafChild(int leafid);

    void getFirstLeafChild();*/

    ~ConstructionNode();

    [[nodiscard]]  int getNumberChildren() const;

    void setAsNonLeaf();

    [[nodiscard]] bool checkTriangleInBox(const Eigen::Ref<const Eigen::Vector3f> &v1,
                                          const Eigen::Ref<const Eigen::Vector3f> &v2,
                                          const Eigen::Ref<const Eigen::Vector3f> &v3) const;


    /**
     * Tests if an edge intersects with any of this node's triangles.
     * @param edge_origin
     * @param edge_end
     * @return
     */
    [[nodiscard]] bool
    checkTrianglesIntersect(const Eigen::Ref<const Eigen::Vector3f> &edge_origin,
                            const Eigen::Ref<const Eigen::Vector3f> &edge_end,
                            const Eigen::Ref<const Eigen::MatrixX3f> &vertices,
                            const Eigen::Ref<const Eigen::MatrixX3i> &triangles) const;

    /**
     * Tests if the edge intersects with given triangle.
     * @param edge_origin
     * @param edge_end
     * @param triangle_vertices
     * @return
     */
    static bool checkTriangleIntersect(const Eigen::Ref<const Eigen::Vector3f> &edge_origin,
                                       const Eigen::Ref<const Eigen::Vector3f> &edge_end,
                                       const Eigen::Ref<const Eigen::Matrix3f> &triangle_vertices);

    /**
     * Tests if an edge intersects with a square box.
     * @param edge_origin
     * @param edge_end
     * @param box_center
     * @param box_side_length
     * @return
     */
    [[nodiscard]] bool checkBoxIntersect(const Eigen::Ref<const Eigen::Vector3f> &edge_origin,
                                         const Eigen::Ref<const Eigen::Vector3f> &edge_end) const;

    static bool boxToEdgeIntersection(const Eigen::Ref<const Eigen::Vector3f> &edge_origin,
                                      const Eigen::Ref<const Eigen::Vector3f> &edge_end,
                                      const Eigen::Ref<const Eigen::Vector3f> &box_center, float box_side_length);

};


#endif //GRAPH_ANALYSIS_CONSTRUCTIONNODE_H
