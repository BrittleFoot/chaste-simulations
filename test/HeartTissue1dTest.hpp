#ifndef HEARTTISSUE1D_HPP
#define HEARTTISSUE1D_HPP


#include <iostream>
#include <cxxtest/TestSuite.h>

#include "LuoRudy1991BackwardEuler.hpp"
#include "SimpleStimulus.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "BidomainProblem.hpp"

#include "ten_tusscher_model_2006_IK1Ko_epi_units.hpp"
#include "ten_tusscher_model_2006_IK1Ko_endo_units.hpp"
#include "ten_tusscher_model_2006_IK1Ko_M_units.hpp"


const double LENGTH = 1.63;


class CellFactory_HeartTissue1d : public AbstractCardiacCellFactory<1>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    CellFactory_HeartTissue1d()
    : AbstractCardiacCellFactory<1>(),
      mpStimulus(new SimpleStimulus(-19000, 2, 0 )) //nA/cm3, ms, ms 
    {
    }
    
    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<1>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        if (x < 1){
            /* Cellten_tusscher_model_2006_IK1Ko_endo_unitsFromCellML */  
            return new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, mpStimulus);            
        }
        else{
            return new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, mpZeroStimulus); 
        }
    }
};



class HeartTissue1dTest : public CxxTest::TestSuite
{
public:
    void TestHeartTissue1d()
    {
        DistributedTetrahedralMesh<1, 1> mesh;
        double h = 0.01; //размерноcть шага по сетке в см
        double length = LENGTH; // length должна делиться на h

        mesh.ConstructRegularSlabMesh(h, length, 0.0, 0.0);
        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(true);
        HeartConfig::Instance()->SetSimulationDuration(600);//ms
        HeartConfig::Instance()->SetOutputDirectory("HeartTissue1d.v2");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetVisualizeWithVtk(true);
        
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.5, 0.5, 0.5));//mS/cm
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.4, 0.4, 0.4));//mS/cm
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.05, 1.0);//ms
        HeartConfig::Instance()->SetCapacitance(0.7); // uF/cm^2
        

        CellFactory_HeartTissue1d cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );
        bidomain_problem.SetMesh(&mesh);
        bidomain_problem.SetWriteInfo();
        bidomain_problem.Initialise();
        bidomain_problem.Solve();
    }
};

#endif /*HEARTTISSUE1D_HPP*/