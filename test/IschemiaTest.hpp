#ifndef ICHEMIATEST_HPP
#define ICHEMIATEST_HPP

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
#include "OdeSystemInformation.hpp"

#include "ten_tusscher_model_2006_endoBackwardEuler.hpp"
#include "ten_tusscher_model_2006_epiBackwardEuler.hpp"
#include "ten_tusscher_model_2006_MBackwardEuler.hpp"

typedef Cellten_tusscher_model_2006_endoFromCellMLBackwardEuler TenTusser2006Endo_BckwardEuler;
typedef Cellten_tusscher_model_2006_epiFromCellMLBackwardEuler  TenTusser2006Epi_BckwardEuler;
typedef Cellten_tusscher_model_2006_MFromCellMLBackwardEuler    TenTusser2006Mid_BckwardEuler;


void InduceIschemia(AbstractCardiacCell* cell) {
    
    // Specifically, the [K+]o was increased from its normal value, 5.4 mmol/L, to 12 mmol/L;15,34 (double)
    cell->SetParameter("rapid_time_dependent_potassium_current_Xr1_gate__Xr1", cell->GetParameter("rapid_time_dependent_potassium_current_Xr1_gate__Xr1") * 2);
    cell->SetParameter("rapid_time_dependent_potassium_current_Xr2_gate__Xr2", cell->GetParameter("rapid_time_dependent_potassium_current_Xr2_gate__Xr2") * 2);
    cell->SetParameter("slow_time_dependent_potassium_current_Xs_gate__Xs", cell->GetParameter("slow_time_dependent_potassium_current_Xs_gate__Xs") * 2);

    // the maximum conductances of the Na+
    cell->SetParameter("membrane_fast_sodium_current_m_gate", cell->GetParameter("membrane_fast_sodium_current_m_gate") * 0.75);
    cell->SetParameter("membrane_fast_sodium_current_h_gate", cell->GetParameter("membrane_fast_sodium_current_h_gate") * 0.75);
    cell->SetParameter("membrane_fast_sodium_current_j_gate", cell->GetParameter("membrane_fast_sodium_current_j_gate") * 0.75);
    // and L-Type Ca2+ channels (currents INa and ICaL, respectively) 
    cell->SetParameter("membrane_L_type_calcium_current_d_gate", cell->GetParameter("membrane_L_type_calcium_current_d_gate") * 0.75);
    cell->SetParameter("membrane_L_type_calcium_current_f_gate", cell->GetParameter("membrane_L_type_calcium_current_f_gate") * 0.75);
    cell->SetParameter("membrane_L_type_calcium_current_f2_gate", cell->GetParameter("membrane_L_type_calcium_current_f2_gate") * 0.75);
    cell->SetParameter("membrane_L_type_calcium_current_fCa_gate", cell->GetParameter("membrane_L_type_calcium_current_fCa_gate") * 0.75);
    // were decreased by 25%, representing inhibition by acidosis
}


class CellFactory_Ichemia : public AbstractCardiacCellFactory<1>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

    double length;
    double stepSize; 

public:
    CellFactory_Ichemia(const double length, const double stepSize, const double stimulus)
        : AbstractCardiacCellFactory<1>(),
          mpStimulus(new SimpleStimulus(stimulus, 2, 0)), //nA/cm3, ms, ms
          length(length), // cm
          stepSize(stepSize) // cm (stepSize * N = length)
    {
    }
    
    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<1>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        AbstractCardiacCell* pCell = NULL;

        // 20% endocardial cells, 30% midwall cells, and 50% epicardial and first cell with stimulus

        if (x <= stepSize * 5)
            pCell = new TenTusser2006Endo_BckwardEuler(mpSolver, mpStimulus);

        else if (x < length * 0.2)
            pCell = new TenTusser2006Endo_BckwardEuler(mpSolver, mpZeroStimulus); 
        else if (x < length * 0.5)
            pCell = new TenTusser2006Mid_BckwardEuler(mpSolver, mpZeroStimulus); 
        else
            pCell = new TenTusser2006Epi_BckwardEuler(mpSolver, mpZeroStimulus);

        return pCell;
    }
};


class IschemiaTest : public CxxTest::TestSuite
{
private:

    double stepSize = 0.01;
    double length = 1.63;
    double stimulus = -100 * 1000;
    DistributedTetrahedralMesh<1, 1> mesh;

public:

    void TestIschemicHeart()
    {
        mesh.ConstructRegularSlabMesh(stepSize, length, 0.0, 0.0);

        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(true);
        HeartConfig::Instance()->SetSimulationDuration(600);//ms
        HeartConfig::Instance()->SetOutputDirectory("IschemiaTest");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetVisualizeWithVtk(true);
        
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(12, 1.3333, 1.3333));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(45, 5, 5));
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.05, 1.0);//ms
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1400);
        HeartConfig::Instance()->SetCapacitance(2.0);
        

        CellFactory_Ichemia cell_factory(length, stepSize, stimulus);

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

    void NoTestPseudoEcgMeasurement() {
        mesh.ConstructRegularSlabMesh(stepSize, length, 0.0, 0.0);
        
        FileFinder directory("IschemiaTest", RelativeTo::ChasteTestOutput);
        std::string h5file = "results";

        // Electrode is here |
        //                   V
        //  oooooooooooo - - E
        ChastePoint<1> electrode(length * 3 / 2);

        PseudoEcgCalculator<1, 1, 1> ecgCalculator(mesh, electrode, directory, h5file);

        try {
            ecgCalculator.WritePseudoEcg();
            std::cout << "WritePseudoEcg" << std::endl;
        }
        catch (const std::exception& e) {
            std::cout << "Exception: " << e.what() << std::endl;
        }

    }
};

#endif /*ICHEMIATEST_HPP*/