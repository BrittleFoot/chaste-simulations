#ifndef HEARTTISSUE1D_HPP
#define HEARTTISSUE1D_HPP

#include <iostream>
#include <cxxtest/TestSuite.h>
#include <string>
#include <cstdlib>

#include "LuoRudy1991BackwardEuler.hpp"
#include "SimpleStimulus.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "BidomainProblem.hpp"
#include "PseudoEcgCalculator.hpp"
#include "FileFinder.hpp"

#include "ten_tusscher_model_2006_endoBackwardEuler.hpp"
#include "ten_tusscher_model_2006_epiBackwardEuler.hpp"
#include "ten_tusscher_model_2006_MBackwardEuler.hpp"

typedef Cellten_tusscher_model_2006_endoFromCellMLBackwardEuler TenTusser2006Endo_BckwardEuler;
typedef Cellten_tusscher_model_2006_epiFromCellMLBackwardEuler  TenTusser2006Epi_BckwardEuler;
typedef Cellten_tusscher_model_2006_MFromCellMLBackwardEuler    TenTusser2006Mid_BckwardEuler;



class CellFactory_HeartTissue1d : public AbstractCardiacCellFactory<1>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

    double length;
    double stepSize; 

public:
    CellFactory_HeartTissue1d(const double length, const double stepSize, const double stimulus)
        : AbstractCardiacCellFactory<1>(),
          mpStimulus(new SimpleStimulus(stimulus, 2, 0)), //nA/cm3, ms, ms
          length(length), // cm
          stepSize(stepSize) // cm (stepSize * N = length)
    {
    }
    
    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<1>* pNode)
    {
        double x = pNode->rGetLocation()[0];

        // 20% endocardial cells, 30% midwall cells, and 50% epicardial and first cell with stimulus

        if (x <= stepSize * 5)
            return new TenTusser2006Endo_BckwardEuler(mpSolver, mpStimulus);
        

        if (x < length * 0.2)
            return new TenTusser2006Endo_BckwardEuler(mpSolver, mpZeroStimulus); 
        
        if (x < length * 0.5)
            return new TenTusser2006Mid_BckwardEuler(mpSolver, mpZeroStimulus); 

        return new TenTusser2006Epi_BckwardEuler(mpSolver, mpZeroStimulus);
    }
};


class HeartTissue1dTest : public CxxTest::TestSuite
{
private:

    double stepSize = 0.01;
    double length = 1.63;
    double stimulus = -100 * 1000;
    DistributedTetrahedralMesh<1, 1> mesh;

public:

    void NoTestHeartTissue1d()
    {
        mesh.ConstructRegularSlabMesh(stepSize, length, 0.0, 0.0);

        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(true);
        HeartConfig::Instance()->SetSimulationDuration(600);//ms
        HeartConfig::Instance()->SetOutputDirectory("HeartTissue1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetVisualizeWithVtk(true);
        
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(12, 1.3333, 1.3333));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(45, 5, 5));
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.05, 1.0);//ms
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1400);
        HeartConfig::Instance()->SetCapacitance(2.0);
        

        CellFactory_HeartTissue1d cell_factory(length, stepSize, stimulus);

        BidomainProblem<1> bidomain_problem( &cell_factory );
        std::cout << "bidomain_problem(.)" << std::endl;

        bidomain_problem.SetMesh(&mesh);
        std::cout << "SetMesh" << std::endl;

        bidomain_problem.SetWriteInfo();
        std::cout << "SetWriteInfo" << std::endl;

        bidomain_problem.Initialise();
        std::cout << "Initialise" << std::endl;

        try {
            bidomain_problem.Solve();
            std::cout << "Solve" << std::endl;
        }
        catch (const std::exception& e) {
            std::cout << "Exception: " << e.what() << std::endl;
        }
    }
};

#endif /*HEARTTISSUE1D_HPP*/