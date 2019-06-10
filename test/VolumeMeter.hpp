#ifndef VOLUME_HPP
#define VOLUME_HPP

#include <iostream>
#include <sstream>
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
#include "VtkMeshReader.hpp"



#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkGeometryFilter.h>
#include <vtkGenericGeometryFilter.h>
#include <vtkDataCompressor.h>
#include <vtkFeatureEdges.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkSmartPointer.h>
#include <vtkDelaunay3D.h>
#include <vtkGeometryFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkMassProperties.h>



class VolumeMeter : public CxxTest::TestSuite
{
private:

    std::string Region(int i) {
        std::stringstream ss;
        ss << "projects/chaste-simulations/ischemic-regions/" << i << ".vtu";
        return ss.str();    
    }


public:

    void TestHeartVolume() {


        try {
            VtkMeshReader<3, 3> meshReader("projects/chaste-simulations/mesh/heart.vtu");

            vtkUnstructuredGrid* pGrid = meshReader.OutputMeshAsVtkUnstructuredGrid();

            auto delaunay = vtkSmartPointer<vtkDelaunay3D>::New();
            auto gFilter = vtkSmartPointer<vtkGeometryFilter>::New();
            auto tFilter = vtkSmartPointer<vtkTriangleFilter>::New();
            auto massProps = vtkSmartPointer<vtkMassProperties>::New();

            delaunay->SetInputData(pGrid);
            gFilter->SetInputConnection(delaunay->GetOutputPort());
            tFilter->SetInputConnection(gFilter->GetOutputPort());
            massProps->SetInputConnection(tFilter->GetOutputPort());

            std::cout << " [---] " << ". Heart volume: " << massProps->GetVolume() << std::endl;
        }
        catch (const std::exception& e) {
            std::cout << " [---] " << ". Exception: " << e.what() << std::endl;
        }
        catch (...) {
            std::cout << " [---] " << ". Undefined exception! " << std::endl;
        }
    }

    void TestPartVolume() {
        // unstructured grid -> vtkDelaunay3D -> vtkGeometryFilter -> vtkTriangleFilter -> vtkMassProperties -> volume.


        for (int i = 1; i < 35; ++i)
        {
            if (i == 9) continue; // 9 volume == 0

            try {
                VtkMeshReader<3, 3> meshReader(Region(i));

                vtkUnstructuredGrid* pGrid = meshReader.OutputMeshAsVtkUnstructuredGrid();

                auto delaunay = vtkSmartPointer<vtkDelaunay3D>::New();
                auto gFilter = vtkSmartPointer<vtkGeometryFilter>::New();
                auto tFilter = vtkSmartPointer<vtkTriangleFilter>::New();
                auto massProps = vtkSmartPointer<vtkMassProperties>::New();

                delaunay->SetInputData(pGrid);
                gFilter->SetInputConnection(delaunay->GetOutputPort());
                tFilter->SetInputConnection(gFilter->GetOutputPort());
                massProps->SetInputConnection(tFilter->GetOutputPort());

                std::cout << " [---] " << i << ". ischemia region volume: " << massProps->GetVolume() << std::endl;
            }
            catch (const std::exception& e) {
                std::cout << " [---] " << i << ". Exception: " << e.what() << std::endl;
            }
            catch (...) {
                std::cout << " [---] " << i << ". Undefined exception! " << std::endl;
            }

        }

    }

};

#endif /*ECG1D_HPP*/