#ifndef JSON_CONFIG_HPP
#define JSON_CONFIG_HPP


#include <iostream>
#include <fstream>
#include <exception>
#include <cxxtest/TestSuite.h>
#include <string>
#include <cstdlib>

#include "FileFinder.hpp"

#include "json.hpp"


using json = nlohmann::json;



/// location relative to chaste-simulations/test
json json_config(const std::string& location) {

    FileFinder file(std::string("projects/chaste-simulations/test/") + location, RelativeTo::ChasteSourceRoot);

    std::cout << "Loading configuration: " << file.GetAbsolutePath() << std::endl;

    std::ifstream input(file.GetAbsolutePath());
    if (!input) {
        throw std::runtime_error(std::string("No such file: ") + file.GetAbsolutePath());
    }

    json config;
    input >> config;
    return config;
}


#endif /*JSON_CONFIG_HPP*/