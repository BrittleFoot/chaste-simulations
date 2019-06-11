#ifndef ONEDIMENSIONEXPERIMENTS
#define ONEDIMENSIONEXPERIMENTS

#include <iostream>
#include <sstream>
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

#include "json.hpp"
#include "JsonConfig.hpp"

#include "ten_tusscher_model_2006_endoBackwardEuler.hpp"
#include "ten_tusscher_model_2006_epiBackwardEuler.hpp"
#include "ten_tusscher_model_2006_MBackwardEuler.hpp"

#include "ten_tusscher_ischemic_model_endoBackwardEuler.hpp"
#include "ten_tusscher_ischemic_model_epiBackwardEuler.hpp"
#include "ten_tusscher_ischemic_model_MBackwardEuler.hpp"


using json = nlohmann::json;

typedef Cellten_tusscher_ischemic_model_endoFromCellMLBackwardEuler IschemicTenTusser2006Endo_BckwardEuler;
typedef Cellten_tusscher_ischemic_model_epiFromCellMLBackwardEuler  IschemicTenTusser2006Epi_BckwardEuler;
typedef Cellten_tusscher_ischemic_model_MFromCellMLBackwardEuler    IschemicTenTusser2006Mid_BckwardEuler;

typedef Cellten_tusscher_model_2006_endoFromCellMLBackwardEuler TenTusser2006Endo_BckwardEuler;
typedef Cellten_tusscher_model_2006_epiFromCellMLBackwardEuler  TenTusser2006Epi_BckwardEuler;
typedef Cellten_tusscher_model_2006_MFromCellMLBackwardEuler    TenTusser2006Mid_BckwardEuler;


void InduceIschemia(AbstractCardiacCell* cell) {
    
    // Specifically, the [K+]o was increased from its normal value, 5.4 mmol/L, to 12 mmol/L;15,34 (double)
    cell->SetParameter("extracellular_potassium_concentration", cell->GetParameter("extracellular_potassium_concentration") * 2);

    // the maximum conductances of the Na+  and L-Type Ca2+ channels (currents INa and ICaL, respectively)
    cell->SetParameter("membrane_fast_sodium_current_conductance", cell->GetParameter("membrane_fast_sodium_current_conductance") * 0.75);
    cell->SetParameter("membrane_calcium_pump_current_conductance", cell->GetParameter("membrane_calcium_pump_current_conductance") * 0.75);
    // were decreased by 25%, representing inhibition by acidosis

    // IK(ATP) was activated by decreasing the intracellular ATP concentration from 6.8 mmol/L to 4.6 mmol/L
    // ATP-sensitive inward-rectifying potassium current (IK(ATP)) 
    // cell->SetParameter("membrane_inward_rectifier_potassium_current_conductance", cell->GetParameter("membrane_inward_rectifier_potassium_current_conductance") * 0.75);

    // and increasing the intracellular ADP concentration from 15 μmol/L to 99 μmol/L
}


class ExperimentalCellFactory : public AbstractCardiacCellFactory<1>
{
private:

    int ncells;

    int id;
    int position;
    double stimulus;
    double length;
    double stepSize; 

    json config = json_config("OneDimensionExperiments.json");


    bool isIschemic(double x) {
        double proc = config["experiment"][id];
        proc = proc / 10;

        // end of stick
        if (position == 1) {
            double lasts = length - x;
            double lastsProc = lasts / length;
            return lastsProc <= proc;
        }
        // begin of stick
        else if (position == 2) {
            double currProc = x / length;
            return currProc <= proc;
        } 
        // middle of stick
        else if (position == 3) {
            double currProc = x / length;
            double midProc = 0.5 - proc / 2;
            double endProc = 0.5 + proc / 2;

            return currProc >= midProc && currProc <= endProc;

        } else {
            return false;
        }


    }

public:
    ExperimentalCellFactory(int id, int position, double length, double stepSize, double stimulus)
        : AbstractCardiacCellFactory<1>(),
          id(id),
          position(position),
          stimulus(stimulus), 
          length(length), // cm
          stepSize(stepSize) // cm (stepSize * N = length)
    {
        ncells = 0;
    }
    
    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<1>* pNode)
    {
        double x = pNode->rGetLocation()[0];

        boost::shared_ptr<AbstractStimulusFunction> mpStimulus(new SimpleStimulus(stimulus, 2, 0));

        if (x > length * 0.1) {
            mpStimulus = mpZeroStimulus;
        }


        // 20% endocardial cells, 30% midwall cells, and 50% epicardial and first cell with stimulus
        AbstractCardiacCell* pCell = NULL;

        if (isIschemic(x)) {

            if (x < length * 0.2)
                pCell = new IschemicTenTusser2006Endo_BckwardEuler(mpSolver, mpStimulus); 
            else if (x < length * 0.5)
                pCell = new IschemicTenTusser2006Epi_BckwardEuler(mpSolver, mpStimulus); 
            else
                pCell = new IschemicTenTusser2006Mid_BckwardEuler(mpSolver, mpStimulus);

            std::cout << "Add another ischemic cell to model (total=" << ++ncells << ")\n";
            InduceIschemia(pCell);
        } else {

            if (x < length * 0.2)
                pCell = new TenTusser2006Endo_BckwardEuler(mpSolver, mpStimulus); 
            else if (x < length * 0.5)
                pCell = new TenTusser2006Epi_BckwardEuler(mpSolver, mpStimulus); 
            else
                pCell = new TenTusser2006Mid_BckwardEuler(mpSolver, mpStimulus);
        }

        if (pCell == NULL) {
            pCell = new TenTusser2006Epi_BckwardEuler(mpSolver, mpZeroStimulus);
        }


        return pCell;
    }
};


class OneDimensionExperiments : public CxxTest::TestSuite
{
private:
    json config = json_config("OneDimensionExperiments.json");

    double stepSize = 0.005;
    double length = 1.63;
    double stimulus = -100 * 1000;
    double electrodePos = 2.4;
    int duration = 600;

    DistributedTetrahedralMesh<1, 1> mesh;


    std::string Name(int i, int position) {
        std::stringstream ss;
        ss << "experiment_" << i << ".position_" << position;
        return ss.str();    
    }


    void TestInstance(int number, int position) {
        std::string sDirectory = Name(number, position);


        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(true);
        HeartConfig::Instance()->SetSimulationDuration(duration);//ms
        HeartConfig::Instance()->SetOutputDirectory(sDirectory);
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetVisualizeWithVtk(true);
        
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(12, 1.3333, 1.3333));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(45, 5, 5));
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.05, 1.0);//ms
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1400);
        HeartConfig::Instance()->SetCapacitance(2.0);
        

        ExperimentalCellFactory cell_factory(number, position, length, stepSize, stimulus);

        BidomainProblem<1> bidomain_problem( &cell_factory );
        std::cout << "bidomain_problem(.)" << std::endl;

        bidomain_problem.SetMesh(&mesh);
        std::cout << "SetMesh" << std::endl;

        bidomain_problem.SetWriteInfo();
        std::cout << "SetWriteInfo" << std::endl;

        bidomain_problem.Initialise();
        std::cout << "Initialise" << std::endl;

        bidomain_problem.Solve();
        std::cout << "Solve" << std::endl;


        std::cout << "Mesuring ECG in `" << sDirectory << "` at " << electrodePos << std::endl; 

        
        FileFinder directory(sDirectory, RelativeTo::ChasteTestOutput);
        std::string h5file = "results";

        ChastePoint<1> electrode(electrodePos);

        PseudoEcgCalculator<1, 1, 1> ecgCalculator(mesh, electrode, directory, h5file);

        ecgCalculator.WritePseudoEcg();
        std::cout << "WritePseudoEcg" << std::endl;
    }

public:

    void TestIschemicHeart()
    {
        mesh.ConstructRegularSlabMesh(stepSize, length, 0.0, 0.0);

        stepSize = config["step-size"];
        length = config["length"];
        electrodePos = config["electrode-position"];
        stimulus = config["stimulus"];
        duration = config["duration"];

        for (int i = 1; i < 35; ++i)
        {
            for (int position = 1; position <= 3; ++position) {
                std::cout << "\n=== TEST #" << i << "\n";
                try {
                    TestInstance(i, position);
                }
                catch (const std::exception& e) {
                    std::cout << "Exception: " << e.what() << std::endl;
                }
                catch (...) {
                    std::cout << "Unexpected error" << std::endl;
                }
                std::cout << "\n===\n";
            }
        }
    }
};

#endif /*ICHEMIATEST_HPP*/
