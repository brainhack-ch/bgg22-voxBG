//
// Created by guibertf on 12/3/22.
//
#include "utils.h"
#include <iostream>
#include <fstream>


int split_face_str_and_convert_int(std::string input_str){
    return std::stoi(input_str.substr(0, input_str.find('/')));
}

void read_and_populate_from_obj_file(std::vector<Eigen::Vector3f> *vertices, std::vector<Eigen::Vector3i> *faces,
                                     std::string filename){
    vertices->clear();
    faces->clear();

    std::ifstream file(filename);
    std::string v, valuesX, valuesY, valuesZ;

    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            // using printf() in all tests for consistency
            std::istringstream iss(line);
            iss >> v >> valuesX >> valuesY >> valuesZ;
            if(v == "v"){
                Eigen::Vector3f vec(std::stof(valuesX), std::stof(valuesY), std::stof(valuesZ));
                vertices->push_back(vec);
            } else {
                if(v == "f"){
                    Eigen::Vector3i face(split_face_str_and_convert_int(valuesX)-1,
                                         split_face_str_and_convert_int(valuesY)-1,
                                         split_face_str_and_convert_int(valuesZ)-1);
                    faces->push_back(face);
                }
            }
        }
        file.close();
    } else {
        std::cout << "Could not open file?!" << std::endl;
        throw std::logic_error("Could not open file, check your path");
    }
}
