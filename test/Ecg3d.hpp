#ifndef ECG3D_HPP
#define ECG3D_HPP

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


class Ecg3d : public CxxTest::TestSuite
{
private:

    json config = json_config("Ecg3dConfig.json"); 


    void DoEcg() {
        std::cout << "\nMeasuring ECG\n" << std::endl; 

        TrianglesMeshReader<3,3> mesh_reader(config["mesh"]);

        AbstractTetrahedralMesh<3,3>* pMesh = new TetrahedralMesh<3,3>;
        pMesh->ConstructFromMeshReader(mesh_reader);

        auto h5 = config["h5"];
        
        FileFinder directory(h5["dir"], RelativeTo::AbsoluteOrCwd);
        ChastePoint<3> electrode(-10.0, -11.0, -12.0);
        std::string h5name = h5["name"];

        PseudoEcgCalculator<3, 3, 1> ecgCalculator(*pMesh, electrode, directory, h5name);

        try {
            std::cout << "Begin WritePseudoEcg" << std::endl;
            ecgCalculator.WritePseudoEcg();
            std::cout << "End WritePseudoEcg" << std::endl;
        }
        catch (const std::exception& e) {
            std::cout << "Exception: " << e.what() << std::endl;
        }
    }


public:

    void Test3dEcg() {

        try {
            DoEcg();
        }
        catch (const std::exception& e) {
            std::cout << "Exception: " << e.what() << std::endl;
        }
    }

};

#endif /*ECG1D_HPP*/