#ifndef ECG1D_HPP
#define ECG1D_HPP

#include <iostream>
#include <cxxtest/TestSuite.h>
#include <string>
#include <cstdlib>
#include <thread>
#include <vector>
#include <set>

#include "BidomainProblem.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PseudoEcgCalculator.hpp"
#include "FileFinder.hpp"
#include "OdeSystemInformation.hpp"

#include "json.hpp"
#include "JsonConfig.hpp"

using json = nlohmann::json;


class Ecg1d : public CxxTest::TestSuite
{
private:

    json config = json_config("EcgConfig.json"); 

    // initialze it from config
    double stepSize;
    double length;
    double electrodePos;


    void DoEcg(const std::string& sDirectory) {
        std::cout << "Mesuring ECG in `" << sDirectory << "` at " << electrodePos << std::endl; 

        DistributedTetrahedralMesh<1, 1> mesh;
        mesh.ConstructRegularSlabMesh(stepSize, length, 0.0, 0.0);
        
        FileFinder directory(sDirectory, RelativeTo::ChasteTestOutput);
        std::string h5file = "results";

        ChastePoint<1> electrode(electrodePos);

        PseudoEcgCalculator<1, 1, 1> ecgCalculator(mesh, electrode, directory, h5file);

        try {
            ecgCalculator.WritePseudoEcg();
            std::cout << "WritePseudoEcg" << std::endl << std::endl;
        }
        catch (const std::exception& e) {
            std::cout << "Exception: " << e.what() << std::endl;
        }
    }

public:


    void TestEcg() {

        stepSize = config["step-size"];
        length = config["length"];
        electrodePos = config["electrode-position"];

        auto vDirectories = config["directories"];
        std::set<std::string> directories(vDirectories.begin(), vDirectories.end());

        std::vector<std::thread> threads;

        for (const auto& directory : directories) {
            threads.push_back(std::thread([this, directory]() { DoEcg(directory); }));
        }

        for (auto& thread : threads) {
            thread.join();
        }
    }

};

#endif /*ECG1D_HPP*/