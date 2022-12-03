//
// Created by guibertf on 12/3/22.
//

#ifndef GRAPH_ANALYSIS_UTILS_H
#define GRAPH_ANALYSIS_UTILS_H

#include <string>
#include <vector>
#include <Eigen/Dense>

int split_face_str_and_convert_int(std::string input_str);

void read_and_populate_from_obj_file(std::vector<Eigen::Vector3f> *vertices, std::vector<Eigen::Vector3i> *faces,
                                     std::string filename);

#endif //GRAPH_ANALYSIS_UTILS_H
