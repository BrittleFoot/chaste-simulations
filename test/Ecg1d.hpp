#ifndef ECG1D_HPP
#define ECG1D_HPP

#include <iostream>
#include <exception>
#include <cxxtest/TestSuite.h>
#include <string>
#include <cstdlib>
#include <vector>
#include <set>

#include "BidomainProblem.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PseudoEcgCalculator.hpp"
#include "FileFinder.hpp"
#include "OdeSystemInformation.hpp"
#include "PetscTools.hpp"

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
            std::cout << "WritePseudoEcg" << std::endl;
            Rename(sDirectory, electrodePos);
            std::cout << "Rename" << std::endl << std::endl;
        }
        catch (const std::exception& e) {
            std::cout << "Exception: " << e.what() << std::endl;
        }
    }

    void Rename(const std::string sDirectory, double electrodePos) {
        std::stringstream outputName, newName;
        outputName << "ChasteResults/output/PseudoEcgFromElectrodeAt_" << electrodePos << "_0_0.dat";
        newName << "ChasteResults/output/PseudoEcgFromElectrodeAt_" << electrodePos << "_0_0." << sDirectory << ".dat";

        FileFinder oldFile(outputName.str(), RelativeTo::ChasteTestOutput);
        FileFinder newFile(newName.str(), RelativeTo::ChasteTestOutput);

        oldFile.CopyTo(newFile);
        oldFile.Remove();

        std::cout << "Renamed `" << oldFile.GetLeafName() << "` to `" << newFile.GetLeafName() << "`" << std::endl;
    }

public:


    void TestEcg() {

        stepSize = config["step-size"];
        length = config["length"];
        electrodePos = config["electrode-position"];

        auto vDirectories = config["directories"];
        std::set<std::string> directories(vDirectories.begin(), vDirectories.end());

        for (const auto& directory : directories) {
            DoEcg(directory);
        }
    }

};

#endif /*ECG1D_HPP*/