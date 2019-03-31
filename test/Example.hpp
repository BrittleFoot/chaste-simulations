#include <cxxtest/TestSuite.h>
#include "BidomainWithBathProblem.hpp"
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

const int MAGIC_NUMBER = 43794;

class HeartCellFactory_HeartTissue1d : public AbstractCardiacCellFactory<3> // <3> here
{
private:
    boost::numeric::ublas::matrix<double> m;

public:
    HeartCellFactory_HeartTissue1d ()
            : AbstractCardiacCellFactory<3>()// <3> here as well!)
    {
        m = boost::numeric::ublas::matrix<double>(MAGIC_NUMBER, 9); //TODO: Memory leak

        std::ifstream infile("heart/test/autorun/HeartTissue1d/data.csv");
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
            if ((x - m(idx,1))*(x - m(idx,1)) + (y - m(idx,2))*(y - m(idx,2)) + (z - m(idx,3))*(z - m(idx,3)) < 0.00001) {
                if (abs(m(idx, 5)) < 1) { //Zero stimulus
                    //TODO: Memory leak
                    if (abs(m(idx, 8) - 0.0) < 0.0001){
                        cell = new CML_tentusscher_panfilov_2006_epi_cell_pe_lut_be(mpSolver, mpZeroStimulus);
                        cell->SetParameter("ScaleFactorGks", m(idx,4)/0.392);
                    }
                    if (abs(m(idx, 8) - 1.0) < 0.0001){
                        cell = new CML_tentusscher_panfilov_2006_M_cell_pe_lut_be(mpSolver, mpZeroStimulus);
                        cell->SetParameter("ScaleFactorGks", m(idx,4)/0.098);
                    }
                    if (abs(m(idx, 8) - 2.0) < 0.0001){
                        cell = new CML_tentusscher_panfilov_2006_endo_cell_pe_lut_be(mpSolver, mpZeroStimulus);
                        cell->SetParameter("ScaleFactorGks", m(idx,4)/0.392);
                    }
                }
                else    //Normal stimulus
                {
                    boost::shared_ptr <SimpleStimulus> mpStimulus(new SimpleStimulus(m(idx, 5), m(idx, 6), m(idx, 7)));
                    //TODO: Memory leak
                    if (abs(m(idx, 8) - 0.0) < 0.0001){
                        cell = new CML_tentusscher_panfilov_2006_epi_cell_pe_lut_be(mpSolver, mpStimulus);
                        cell->SetParameter("ScaleFactorGks", m(idx,4)/0.392);
                    }
                    if (abs(m(idx, 8) - 1.0) < 0.0001){
                        cell = new CML_tentusscher_panfilov_2006_M_cell_pe_lut_be(mpSolver, mpStimulus);
                        cell->SetParameter("ScaleFactorGks", m(idx,4)/0.098);
                    }
                    if (abs(m(idx, 8) - 2.0) < 0.0001){
                        cell = new CML_tentusscher_panfilov_2006_endo_cell_pe_lut_be(mpSolver, mpStimulus);
                        cell->SetParameter("ScaleFactorGks", m(idx,4)/0.392);
                    }
                }

                break;
            }
        }

        if (cell == NULL)
        {
            cell = new CML_tentusscher_panfilov_2006_epi_cell_pe_lut_be(mpSolver, mpZeroStimulus);
        }
        PRINT_VARIABLE(node_index);
        return cell;
    }
};

class HeartTissue1d : public CxxTest::TestSuite
{
public:
    double GetMemoryUsage()
    {
        struct rusage rusage;
        getrusage( RUSAGE_SELF, &rusage );

        return (double)(rusage.ru_maxrss)/(1024);// Convert KB to MB
    }

    void Test_HeartTissue1d() throw(Exception)
    {
        int rank, size;
        MPI_Comm_rank (PETSC_COMM_WORLD, &rank);	/* get current process id */
        MPI_Comm_size (PETSC_COMM_WORLD, &size);	/* get number of processes */
        PRINT_VARIABLE(rank);
        PRINT_VARIABLE(size);

        HeartConfig::Instance()->SetDomain(cp::domain_type::BiWithBath);
        HeartConfig::Instance()->SetMeshFileName("heart/test/autorun/HeartTissue1d/data", cp::media_type::Axisymmetric);

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(12, 1.3333, 1.3333));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(45, 5, 5));
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio( 1400 );
        HeartConfig::Instance()->SetCapacitance( 2.0 );

        std::set<unsigned> tissue_ids;
        static unsigned tissue_id=0;
        tissue_ids.insert(tissue_id);

        std::set<unsigned> bath_ids;
        static unsigned bath_id1=1;
        bath_ids.insert(bath_id1);
        static unsigned bath_id2=2;
        bath_ids.insert(bath_id2);
	    static unsigned bath_id3=3;
        bath_ids.insert(bath_id3);

        std::map< unsigned, double > cond_map;
        cond_map[1] = 1.; // Lung
        cond_map[2] = 1.; // Lung
        cond_map[3] = 7.; // Torso

        HeartConfig::Instance()->SetTissueAndBathIdentifiers(tissue_ids, bath_ids);
        HeartConfig::Instance()->SetBathMultipleConductivities(cond_map);

        HeartConfig::Instance()->SetSimulationDuration(600); //ms
        HeartConfig::Instance()->SetOutputDirectory("HeartTissue1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetVisualizeWithVtk(false);
        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(false);

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.05, 1.);

        HeartCellFactory_HeartTissue1d cell_factory;
        BidomainWithBathProblem<3> bidomain_problem( &cell_factory );
        TRACE("bidomain_problem(.) complete!");

        bidomain_problem.SetWriteInfo();
        TRACE("SetWriteInfo() complete!");
        
        bidomain_problem.Initialise();
        TRACE("Initialise() complete!");
        bidomain_problem.Solve();

        AbstractTetrahedralMesh<3,3>* p_mesh = &(bidomain_problem.rGetMesh());
        VtkMeshWriter<3,3> vtk_writer("./HeartTissue1d/", "init_mesh_1", false);
        vtk_writer.WriteFilesUsingMesh(*p_mesh);

        VtkMeshWriter<3,3> vtk_writer2("./HeartTissue1d/", "init_mesh_2", false);
        std::string original_file = p_mesh->GetMeshFileBaseName();
        std::auto_ptr<AbstractMeshReader<3, 3> > p_original_mesh_reader = GenericMeshReader<3, 3>(original_file);
        vtk_writer2.WriteFilesUsingMeshReader(*p_original_mesh_reader);

    }
};
