#ifndef ICHEMIATEST_HPP
#define ICHEMIATEST_HPP

#include <iostream>
#include <iomanip>
#include <cxxtest/TestSuite.h>
#include <string>
#include <vector>
#include <set>
#include <cstdlib>

#include "LuoRudy1991BackwardEuler.hpp"
#include "SimpleStimulus.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "BidomainProblem.hpp"
#include "PseudoEcgCalculator.hpp"
#include "FileFinder.hpp"
#include "OdeSystemInformation.hpp"
#include "VtkMeshWriter.hpp"


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


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

typedef Cellten_tusscher_model_2006_endoFromCellMLBackwardEuler     TenTusser2006Endo_BckwardEuler;
typedef Cellten_tusscher_model_2006_epiFromCellMLBackwardEuler      TenTusser2006Epi_BckwardEuler;
typedef Cellten_tusscher_model_2006_MFromCellMLBackwardEuler        TenTusser2006Mid_BckwardEuler;


const int CSV_LINENO = 44515;


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



class CellFactory_Ichemia : public AbstractCardiacCellFactory<3>
{
private:

    json config = json_config("IschemiaProblemHeartConfig.json");

    boost::numeric::ublas::matrix<double> cellMapping;
    std::vector<c_vector<double, 3>> ischemicCenters;


    double get_idx(int line_num)        { return cellMapping(line_num, 0); }
    double get_x(int line_num)          { return cellMapping(line_num, 1); }
    double get_y(int line_num)          { return cellMapping(line_num, 2); }
    double get_z(int line_num)          { return cellMapping(line_num, 3); }
    double get_d_Ks(int line_num)       { return cellMapping(line_num, 4); }
    double get_amplitude(int line_num)  { return cellMapping(line_num, 5); }
    double get_duration(int line_num)   { return cellMapping(line_num, 6); }
    double get_time(int line_num)       { return cellMapping(line_num, 7); }
    double get_model(int line_num)      { return cellMapping(line_num, 8); }


    bool isIschemicCell(double x, double y, double z) {

        for (auto center : ischemicCenters) {
            auto cx = center[0];
            auto cy = center[1];
            auto cz = center[2];

            if ((x - cx)*(x - cx) + (y - cy)*(y - cy) + (z - cz)*(z - cz) < 15 * 15) {
                return true;
            }
        }   
        return false;
    }


public:
        
    CellFactory_Ichemia() : AbstractCardiacCellFactory<3>()
    {

        cellMapping = boost::numeric::ublas::matrix<double>(CSV_LINENO, 9);

        std::string csv = config["csv"];
        FileFinder ff(csv, RelativeTo::ChasteSourceRoot);

        std::ifstream infile(ff.GetAbsolutePath());

        double idx, x, y, z, d_Ks, amplitude, duration, time, model;

        int line_num = 0;
        while (infile >> idx >> x >> y >> z >> d_Ks >> amplitude >> duration >> time >> model)
        {
            cellMapping(line_num, 0) = idx;
            cellMapping(line_num, 1) = x;
            cellMapping(line_num, 2) = y;
            cellMapping(line_num, 3) = z;
            cellMapping(line_num, 4) = d_Ks;
            cellMapping(line_num, 5) = amplitude;
            cellMapping(line_num, 6) = duration;
            cellMapping(line_num, 7) = time;
            cellMapping(line_num, 8) = model;
            line_num++;
        }

        std::cout << std::setw(10) << "idx " << get_idx(5) << std::endl;
        std::cout << std::setw(10) << "x " << get_x(5) << std::endl;
        std::cout << std::setw(10) << "y " << get_y(5) << std::endl;
        std::cout << std::setw(10) << "z " << get_z(5) << std::endl;
        std::cout << std::setw(10) << "d_Ks " << get_d_Ks(5) << std::endl;
        std::cout << std::setw(10) << "amplitude " << get_amplitude(5) << std::endl;
        std::cout << std::setw(10) << "duration " << get_duration(5) << std::endl;
        std::cout << std::setw(10) << "time " << get_time(5) << std::endl;
        std::cout << std::setw(10) << "model " << get_model(5) << std::endl;


        for (auto region: config["ischemia"]) {
            ischemicCenters.push_back(Create_c_vector(region[0], region[1], region[2]));
        }
    }
    
    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {

        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
        double z = pNode->rGetLocation()[2];

        bool ischemiaApplied = config["apply-ischemia"];

        AbstractCardiacCell* cell = NULL;

        
        for(int idx=0; idx < CSV_LINENO; idx++){
            // find this cell in csv and apply supplied config 
            if ((x - get_x(idx))*(x - get_x(idx)) + (y - get_y(idx))*(y - get_y(idx)) + (z - get_z(idx))*(z - get_z(idx)) >= 0.001) {
                continue;
            }

            boost::shared_ptr<AbstractStimulusFunction> mpStimulus(new SimpleStimulus(
                get_amplitude(idx),
                get_duration(idx),
                get_time(idx)
            ));

            if (abs(get_amplitude(idx) <= 1)) {
                mpStimulus = mpZeroStimulus;
            }

            double cellModel = get_model(idx);

            if (ischemiaApplied && isIschemicCell(x, y, z)) {
                
                if (cellModel <= 0.2){
                    cell = new IschemicTenTusser2006Endo_BckwardEuler(mpSolver, mpStimulus);
                }
                if ((cellModel >= 0.2) && (cellModel <= 0.5)){
                    cell = new IschemicTenTusser2006Mid_BckwardEuler(mpSolver, mpStimulus);
                }
                if (cellModel >= 0.5){
                    cell = new IschemicTenTusser2006Epi_BckwardEuler(mpSolver, mpStimulus);
                }
                InduceIschemia(cell);
            }
            else {

                if (cellModel <= 0.2){
                    cell = new TenTusser2006Endo_BckwardEuler(mpSolver, mpStimulus);
                }
                if ((cellModel >= 0.2) && (cellModel <= 0.5)){
                    cell = new TenTusser2006Mid_BckwardEuler(mpSolver, mpStimulus);
                }
                if (cellModel >= 0.5){
                    cell = new TenTusser2006Epi_BckwardEuler(mpSolver, mpStimulus);
                }
                
            }
            
           break;
        
        }

        if (cell == NULL) {
            cell = new TenTusser2006Epi_BckwardEuler(mpSolver, mpZeroStimulus);
        }

        return cell;
    }
};


class IschemiaTest : public CxxTest::TestSuite
{
private:

    json config = json_config("IschemiaProblemHeartConfig.json");

public:
    
    void TestIschemicHeart()
    {
        HeartConfig::Instance()->SetMeshFileName(config["mesh"], cp::media_type::Axisymmetric);
        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(true);
        HeartConfig::Instance()->SetSimulationDuration(600);
        HeartConfig::Instance()->SetOutputDirectory("IschemiaOnHeart");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetVisualizeWithVtk(false);

        
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(12, 1.3333, 1.3333));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(45, 5, 5));
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.05, 1.0);
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1400);
        HeartConfig::Instance()->SetCapacitance(2.0);
        

        AbstractCardiacCellFactory<3> *pCellFactory = new CellFactory_Ichemia();

        BidomainProblem<3> bidomain_problem( pCellFactory );
        std::cout << "bidomain_problem(.)" << std::endl;


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

        AbstractTetrahedralMesh<3,3>* pMesh = &(bidomain_problem.rGetMesh());
        VtkMeshWriter<3,3> vtkWriter(HeartConfig::Instance()->GetOutputDirectory(), "init_mesh_1", false);
        vtkWriter.WriteFilesUsingMesh(*pMesh);

        VtkMeshWriter<3,3> vtkWriter2(HeartConfig::Instance()->GetOutputDirectory(), "init_mesh_2", false);
        std::string originalFile = pMesh->GetMeshFileBaseName();
        std::shared_ptr<AbstractMeshReader<3, 3> > pOriginalMeshWriter = GenericMeshReader<3, 3>(originalFile);
        vtkWriter2.WriteFilesUsingMeshReader(*pOriginalMeshWriter);
    }
};

#endif /*ICHEMIATEST_HPP*/