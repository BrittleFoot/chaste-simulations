#ifndef TESTSTICK5CM_HPP
#define TESTSTICK5CM_HPP


#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "LuoRudy1991BackwardEuler.hpp"
#include "SimpleStimulus.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "BidomainProblem.hpp"


class CellFactory_measurement_of_speed : public AbstractCardiacCellFactory<1>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;
public:
    CellFactory_measurement_of_speed()
    : AbstractCardiacCellFactory<1>(),
      mpStimulus(new SimpleStimulus(-19000, 2, 0 )) //nA/cm3, ms, ms 
    {
    }
    
    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<1>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        if (x < 1){
            return new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, mpStimulus);            
        }
        else{
            return new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, mpZeroStimulus); 
        }
        
        
    }

};



class Stick5cm : public CxxTest::TestSuite
{
public:
    void TestStick5cm()
    {
        DistributedTetrahedralMesh<1, 1> mesh;
        double h = 0.025; //размерноcть шага по сетке в см
        mesh.ConstructRegularSlabMesh(h, 5, 0.0, 0.0);
        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(true);
        HeartConfig::Instance()->SetSimulationDuration(800);//ms
        HeartConfig::Instance()->SetOutputDirectory("measurement_of_speed_output_luo_rudy2");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetVisualizeWithVtk(true);
        //HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(2.81, 2.81,2.81));//mS/cm
        //HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.61, 0.61, 0.61));//mS/cm
        
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.5, 0.5, 0.5));//mS/cm
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.4, 0.4, 0.4));//mS/cm
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.05, 1.0);//ms
        HeartConfig::Instance()->SetCapacitance(0.7); // uF/cm^2
        

        CellFactory_measurement_of_speed cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );
        bidomain_problem.SetMesh(&mesh);
        bidomain_problem.SetWriteInfo();
        bidomain_problem.Initialise();
        bidomain_problem.Solve();
    }
};

#endif /*TESTSTICK5CM_HPP*/