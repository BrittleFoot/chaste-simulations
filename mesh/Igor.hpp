#include <cxxtest/TestSuite.h>
//#include "BidomainWithBathProblem.hpp"
#include "BidomainProblem.hpp"
#include "ten_tusscher_model_2006_IK1Ko_epi_units.hpp"
#include "ten_tusscher_model_2006_IK1Ko_endo_units.hpp"
#include "ten_tusscher_model_2006_IK1Ko_M_units.hpp"
#include "GenericMeshReader.hpp"
#include "SimpleStimulus.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "FixedModifier.hpp"
#include <fstream>
#include <string>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <unistd.h>
#include <sys/resource.h>
#include "Debug.hpp"
#include <mpi.h>
#include "VtkMeshWriter.hpp"
#include "LuoRudy1991BackwardEuler.hpp"
#include "MonodomainProblem.hpp"

const int MAGIC_NUMBER = 44515;

class HeartCellFactory_CASE_R4_epi : public AbstractCardiacCellFactory<3> // <3> here
{
private:
    boost::numeric::ublas::matrix<double> m;

public:
    HeartCellFactory_CASE_R4_epi ()
            : AbstractCardiacCellFactory<3>()// <3> here as well!)
    {
        m = boost::numeric::ublas::matrix<double>(MAGIC_NUMBER, 9); //TODO: Memory leak

        std::ifstream infile("heart/test/autorun/IGOR/data.csv");
        double idx, x, y, z, d_Ks, amplitude, duration, time, model;
        int line_num = 0;
        while (infile >> idx >> x >> y >> z >> d_Ks >> amplitude >> duration >> time >> model)
        {
            m(line_num, 0) = idx;
            m(line_num, 1) = x;
            m(line_num, 2) = y;
            m(line_num, 3) = z;
            m(line_num, 4) = d_Ks;
            m(line_num, 5) = amplitude;
            m(line_num, 6) = duration;
            m(line_num, 7) = time;
            m(line_num, 8) = model;
            line_num++;
        }

        PRINT_VARIABLE( m(5, 0));
        PRINT_VARIABLE( m(5, 1));
        PRINT_VARIABLE( m(5, 2));
        PRINT_VARIABLE( m(5, 3));
        PRINT_VARIABLE( m(5, 4));
        PRINT_VARIABLE( m(5, 5));
        PRINT_VARIABLE( m(5, 6));
        PRINT_VARIABLE( m(5, 7));
        PRINT_VARIABLE( m(5, 8));

    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
        double z = pNode->rGetLocation()[2];
        unsigned int node_index = pNode->GetIndex();

        AbstractCardiacCell* cell = NULL;

        for(int idx=0; idx < MAGIC_NUMBER; idx++){
            if ((x - m(idx,1))*(x - m(idx,1)) + (y - m(idx,2))*(y - m(idx,2)) + (z - m(idx,3))*(z - m(idx,3)) < 0.001) {

                if (abs(m(idx, 5)) > 1){
                    
                    boost::shared_ptr <SimpleStimulus> mpStimulus(new SimpleStimulus(m(idx, 5), m(idx, 6), m(idx, 7)));
                    //TODO: Memory leak  
                    if (m(idx, 8)<= 0.2){
                        cell = new CML_tentusscher_panfilov_2006_endo_cell_pe_lut_be(mpSolver, mpStimulus);
                        cell->SetParameter("ScaleFactorGks", m(idx,4)/0.392);
                    }
                    if ((m(idx, 8) >= 0.2) && (m(idx, 8)<= 0.5)){
                        cell = new CML_tentusscher_panfilov_2006_M_cell_pe_lut_be(mpSolver, mpStimulus);
                        cell->SetParameter("ScaleFactorGks", m(idx,4)/0.098);                        
                    }
                    if (m(idx, 8) >= 0.5){
                        cell = new CML_tentusscher_panfilov_2006_epi_cell_pe_lut_be(mpSolver, mpStimulus);
                        cell->SetParameter("ScaleFactorGks", m(idx,4)/0.392);
                    }
                }
                else {
                    if (m(idx, 8)<= 0.2){
                        cell = new CML_tentusscher_panfilov_2006_endo_cell_pe_lut_be(mpSolver, mpZeroStimulus);
                        cell->SetParameter("ScaleFactorGks", m(idx,4)/0.392);
                    }
                    if ((m(idx, 8) >= 0.2) && (m(idx, 8)<= 0.5)){
                        cell = new CML_tentusscher_panfilov_2006_M_cell_pe_lut_be(mpSolver, mpZeroStimulus);
                        cell->SetParameter("ScaleFactorGks", m(idx,4)/0.098);
                    }
                    if (m(idx, 8) >= 0.5){
                        cell = new CML_tentusscher_panfilov_2006_epi_cell_pe_lut_be(mpSolver, mpZeroStimulus);
                        cell->SetParameter("ScaleFactorGks", m(idx,4)/0.392);
                    }
                    
                }
                
               break; 
                                   
            }
        }

        if (cell == NULL)
        {
            cell = new CML_tentusscher_panfilov_2006_endo_cell_pe_lut_be(mpSolver, mpZeroStimulus);
        }
        PRINT_VARIABLE(node_index);
        return cell;
    }
};

class TestClass_amikard_CASE_R4_epi : public CxxTest::TestSuite
{
public:
    double GetMemoryUsage()
    {
        struct rusage rusage;
        getrusage( RUSAGE_SELF, &rusage );

        return (double)(rusage.ru_maxrss)/(1024);// Convert KB to MB
    }

    void Test_CASE_R4_epi() throw(Exception)
    {
        int rank, size;
        MPI_Comm_rank (PETSC_COMM_WORLD, &rank);	/* get current process id */
        MPI_Comm_size (PETSC_COMM_WORLD, &size);	/* get number of processes */
        PRINT_VARIABLE(rank);
        PRINT_VARIABLE(size);

        //HeartConfig::Instance()->SetDomain(cp::domain_type::BiWithBath);
        HeartConfig::Instance()->SetMeshFileName("heart/test/autorun/IGOR/data", cp::media_type::Axisymmetric);

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(12, 1.3333, 1.3333));   //1.3333, 1.3333
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(45, 5, 5));  //45, 5, 5
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1400);
        HeartConfig::Instance()->SetCapacitance(2.0);

        HeartConfig::Instance()->SetSimulationDuration(100); //ms
        HeartConfig::Instance()->SetOutputDirectory("IGOR");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetVisualizeWithVtk(true);
        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(false);

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.05, 1.);

        HeartCellFactory_CASE_R4_epi cell_factory;
        BidomainProblem<3> bidomain_problem( &cell_factory );
        TRACE("bidomain_problem(.) complete!");

        bidomain_problem.SetWriteInfo();
        TRACE("SetWriteInfo() complete!");
        
        bidomain_problem.Initialise();
        TRACE("Initialise() complete!");
        bidomain_problem.Solve();

        AbstractTetrahedralMesh<3,3>* p_mesh = &(bidomain_problem.rGetMesh());
        VtkMeshWriter<3,3> vtk_writer("./IGOR/", "init_mesh_1", false);
        vtk_writer.WriteFilesUsingMesh(*p_mesh);

        VtkMeshWriter<3,3> vtk_writer2("./IGOR/", "init_mesh_2", false);
        std::string original_file = p_mesh->GetMeshFileBaseName();
        std::auto_ptr<AbstractMeshReader<3, 3> > p_original_mesh_reader = GenericMeshReader<3, 3>(original_file);
        vtk_writer2.WriteFilesUsingMeshReader(*p_original_mesh_reader);

    }
};
