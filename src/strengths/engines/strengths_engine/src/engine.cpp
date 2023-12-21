#include <iostream>
#include <random>

#include "SimulationAlgorithm3DBase.hpp"
#include "Euler3D.hpp"
#include "TauLeap3D.hpp"
#include "Gillespie3D.hpp"

#include "SimulationAlgorithmGraphBase.hpp"
#include "EulerGraph.hpp"
#include "TauLeapGraph.hpp"
#include "GillespieGraph.hpp"

#include <chrono>

#ifdef CPYEMVER
  #include <Python.h>

  extern "C" PyObject * PyInit_engine()
    {
    return NULL;
    }
#endif

std::vector<double> GenerateStochasticDistribution (std::vector<double> mesh_x, int n_meshes, int n_species, int seed)
  {
  /// generate a poisson distributed stochastic state that respects the floored total quantities of the input floating point state.

  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> uiud(0, 1);

  std::vector<double> mesh_x_sto = std::vector<double>(mesh_x.size(), 0);
  std::vector<double> tot_species(n_species, 0);
  std::vector<double> tot2_species(n_species, 0);
  std::vector<double> dtot_species(n_species, 0);

  //setp 1 : count species

  for(int i=0; i<n_meshes; i++)
    {
    for(int j=0; j<n_species; j++)
      {
      tot_species[j] += mesh_x[i*n_species+j];
      }
    }

  for(int i=0; i<n_species; i++)
    {
    tot_species[i] = std::floor(tot_species[i]);
    }

  // step 2 : poisson_distribution

  for(int i=0; i<n_meshes*n_species; i++)
    {
    if (mesh_x[i]<100)
      {
      mesh_x_sto[i] = std::poisson_distribution<int>(mesh_x[i])(rng);
      }
    else
      {
      mesh_x_sto[i] = std::max(0.0, std::floor(std::normal_distribution<double>(mesh_x[i], sqrt(mesh_x[i]))(rng)));
      }
    }

  // step 3 : count species on the poisson distributed state

  for(int i=0; i<n_meshes; i++)
    {
    for(int j=0; j<n_species; j++)
      {
      tot2_species[j] += mesh_x_sto[i*n_species+j];
      }
    }

  // step 4 : difference of the two totals

  for(int i=0; i<n_species; i++)
    {
    dtot_species[i] = tot2_species[i] - tot_species[i];
    }

  // step 5 : correction of the poisson distributed state

  for(int s=0; s<n_species; s++)
    {

    int delta = static_cast<int>(dtot_species[s]);
    if(delta == 0) continue;

    bool rm_species = false;
    if(delta>0) rm_species = true;

    delta = abs(delta);
    int delta_count = 0;

    for(;;)
      {
      double cumul = 0;
      double target = uiud(rng) * tot_species[s];

      for(int i=0; i<n_meshes; i++)
        {
        cumul += mesh_x[i*n_species+s];
        if(target<cumul)
          {
          if(rm_species) // remoive species
            {
            if(mesh_x_sto[i*n_species+s]>0)
              {
              mesh_x_sto[i*n_species+s]--;
              delta_count++;
              }
            }
          else // add species
            {
            mesh_x_sto[i*n_species+s]++;
            delta_count++;
            }
          break;
          }
        }
      if(delta_count == delta) break;
      }
    }

  return mesh_x_sto;
  }

template<typename T_out, typename T_in> std::vector<T_out> MkVec(T_in * a, int len)
    {
    std::vector<T_out> v(len);
    for(int i=0;i<len;i++)
        {
        v[i] = static_cast<T_out>(a[i]);
        }
    return v;
    }

SimulationAlgorithm3DBase * global_grid_algo;
SimulationAlgorithmGraphBase * global_graph_algo;
int global_space_type; //0 : grid, 1 : graph
bool global_algo_freed = true;

template<typename T> std::vector<T> SpeciesFirstToMeshFirstArray(std::vector<T> species_first_array, int n_species, int n_meshes)
    {
    std::vector<T> mesh_first_array(species_first_array.size());

    for(int s=0;s<n_species;s++)
        for(int i=0;i<n_meshes;i++)
            mesh_first_array[i*n_species+s] = species_first_array[s*n_meshes+i];
    return mesh_first_array;
    }

bool CompareStr(const char * str1, const char * str2)
    {
    return (std::string(str1) == std::string(str2));
    }

extern "C" int Initialize3D (
    int w,               //system width
    int h,               //system height
    int d,               //system depth
    int n_species,       //number of species
    int n_reactions,     //number of reactions
    int n_env,           //number of environments
    double * mesh_state, //species quantities
                         //size N*width*height*depth
                         //species first array : x = [species[depth[height[width]]]] or [species[mesh]]
    int *    mesh_chstt, //species chemostat flag
                         //size N*width*height*depth
                         //species first array : x = [species[depth[height[width]]]] or [species[mesh]]
    int *    mesh_env,   //species chemostat flag //size N*width*height*depth
    double mesh_vol,     //volume of a square mesh
    double * k,          //reaction rates //size M
    int * sub,           //N*M substrate matrix // (rows * columns) //size N*M
    int * sto,           //N*M stoechiometry matrix //size N*M
    int * r_env,         //reactions environmenrs
    double * D,          //diffusion coefficient for each species in each mesh type
    const char * boundary_conditions_x, //boundary consitions to be applied along the x axis
    const char * boundary_conditions_y, //boundary consitions to be applied along the y axis
    const char * boundary_conditions_z, //boundary consitions to be applied along the z axis
    int sample_n,        //number of sample timepoints
    double * sample_t,   //sample timepoints //size sample_n

    const char * sampling_policy, //tells how the sampling should be done
    double sampling_interval,     //time interval at which the system should be sampled, if used
    double t_max,                 //time past which the simulation should be stopped

    double time_step,    //time step
    int seed,            //rng seed
    const char * option  //option
    )
    //return codes :
    //  0 : success
    //  1 : invalid option
    //  2 : invalid boudary condition
    //  3 : invalid sampling policy
    {
    global_space_type = 0;
    int n_meshes = w*h*d;

    std::vector<int> boundary_conditions(3);

    // boundary conditions
    if      (CompareStr(boundary_conditions_x, "reflecting")) boundary_conditions[0] = 0;
    else if (CompareStr(boundary_conditions_x, "periodical")) boundary_conditions[0] = 1;
    else return 2;

    if      (CompareStr(boundary_conditions_y, "reflecting")) boundary_conditions[1] = 0;
    else if (CompareStr(boundary_conditions_y, "periodical")) boundary_conditions[1] = 1;
    else return 2;

    if      (CompareStr(boundary_conditions_z, "reflecting")) boundary_conditions[2] = 0;
    else if (CompareStr(boundary_conditions_z, "periodical")) boundary_conditions[2] = 1;
    else return 2;

    // sampling_policy
    int sampling_policy_code;
    if     (CompareStr(sampling_policy, "on_t_sample" )) sampling_policy_code = 0;
    else if(CompareStr(sampling_policy, "on_iteration")) sampling_policy_code = 1;
    else if(CompareStr(sampling_policy, "on_interval" )) sampling_policy_code = 2;
    else if(CompareStr(sampling_policy, "no_sampling" )) sampling_policy_code = 3;
    else return 3;

    // option
    if      (CompareStr(option, "gillespie"))   {global_grid_algo = new Gillespie3D(); global_algo_freed = false;}
    else if (CompareStr(option, "tauleap"))     {global_grid_algo = new TauLeap3D();   global_algo_freed = false;}
    else if (CompareStr(option, "euler"))       {global_grid_algo = new Euler3D();     global_algo_freed = false;}
    else return 1;

    std::vector<double> mesh_x;
    bool is_stochastic = (CompareStr(option, "tauleap") || CompareStr(option, "gillespie"));

    if(is_stochastic)
      {
      mesh_x = GenerateStochasticDistribution (
        SpeciesFirstToMeshFirstArray(MkVec<double, double>(mesh_state, n_meshes*n_species), n_species, n_meshes),
        n_meshes,
        n_species,
        seed);
      }
    else
      {
      mesh_x = SpeciesFirstToMeshFirstArray(MkVec<double, double>(mesh_state, n_meshes*n_species), n_species, n_meshes);
      }

    global_grid_algo->Init(
          w,
          h,
          d,
          n_species,
          n_reactions,
          n_env,
          mesh_x, //species first to mesh first
          SpeciesFirstToMeshFirstArray(MkVec<int,    int   >(mesh_chstt, n_meshes*n_species),
                                       n_species,
                                       n_meshes), //species first to mesh first
          MkVec<int,    int   >(mesh_env, n_meshes),
          mesh_vol,
          MkVec<double, double>(k, n_reactions),
          MkVec<double, int   >(sub, n_species*n_reactions),
          MkVec<double, int   >(sto, n_species*n_reactions),
          MkVec<double, int   >(r_env, n_reactions*n_env),
          MkVec<double, double>(D, n_species*n_env),
          boundary_conditions,
          sample_n,
          MkVec<double, double>(sample_t, sample_n),

          sampling_policy_code,
          sampling_interval,
          t_max,

          time_step,
          seed
          );

    return 0;
    }

extern "C" int InitializeGraph (
    int n_nodes,         //system width
    int n_species,       //number of species
    int n_reactions,     //number of reactions
    int n_env,           //number of environments

    int n_edges,         //number of edges
    int * edge_i,        //edge first node index
    int * edge_j,        //edge second node index
    double * edge_sfc,   //edge nodes contact surface
    double * edge_dst,   //edge nodes center distances

    double * mesh_state, //species quantities
                         //size N*width*height*depth
                         //species first array : x = [species[depth[height[width]]]] or [species[mesh]]
    int *    mesh_chstt, //species chemostat flag
                         //size N*width*height*depth
                         //species first array : x = [species[depth[height[width]]]] or [species[mesh]]
    int *    mesh_env,   //species chemostat flag //size N*width*height*depth
    double * mesh_vol,     //volume of a square mesh

    double * k,          //reaction rates //size M
    int * sub,           //N*M substrate matrix // (rows * columns) //size N*M
    int * sto,           //N*M stoechiometry matrix //size N*M
    int * r_env,         //reactions environmenrs
    double * D,          //diffusion coefficient for each species in each mesh type

    int sample_n,        //number of sample timepoints
    double * sample_t,   //sample timepoints //size sample_n

    const char * sampling_policy, //tells how the sampling should be done
    double sampling_interval,     //time interval at which the system should be sampled, if used
    double t_max,                 //time past which the simulation should be stopped

    double time_step,    //time step
    int seed,            //rng seed
    const char * option  //option
    )
    //return codes :
    //  0 : success
    //  1 : invalid option
    //  2 : invalid boudary condition
    //  3 : invalid sampling policy
    {
    global_space_type = 1;
    int n_meshes = n_nodes;

    // sampling_policy
    int sampling_policy_code;
    if     (CompareStr(sampling_policy, "on_t_sample" )) sampling_policy_code = 0;
    else if(CompareStr(sampling_policy, "on_iteration")) sampling_policy_code = 1;
    else if(CompareStr(sampling_policy, "on_interval" )) sampling_policy_code = 2;
    else if(CompareStr(sampling_policy, "no_sampling" )) sampling_policy_code = 3;
    else return 3;

    // option
    if      (CompareStr(option, "gillespie"))   {global_graph_algo = new GillespieGraph(); global_algo_freed = false;}
    else if (CompareStr(option, "tauleap"))     {global_graph_algo = new TauLeapGraph();   global_algo_freed = false;}
    else if (CompareStr(option, "euler"))       {global_graph_algo = new EulerGraph();     global_algo_freed = false;}
    else return 1;

    std::vector<double> mesh_x;
    bool is_stochastic = (CompareStr(option, "tauleap") || CompareStr(option, "gillespie"));

    if(is_stochastic)
      {
      mesh_x = GenerateStochasticDistribution (
        SpeciesFirstToMeshFirstArray(MkVec<double, double>(mesh_state, n_meshes*n_species), n_species, n_meshes),
        n_meshes,
        n_species,
        seed);
      }
    else
      {
      mesh_x = SpeciesFirstToMeshFirstArray(MkVec<double, double>(mesh_state, n_meshes*n_species), n_species, n_meshes);
      }

    global_graph_algo->Init(
          n_meshes,
          n_species,
          n_reactions,
          n_env,

          n_edges,
          MkVec<int,    int   >(edge_i, n_edges),
          MkVec<int,    int   >(edge_j, n_edges),
          MkVec<double,    double   >(edge_sfc, n_edges),
          MkVec<double,    double   >(edge_dst, n_edges),


          mesh_x, //species first to mesh first
          SpeciesFirstToMeshFirstArray(MkVec<int,    int   >(mesh_chstt, n_meshes*n_species),
                                       n_species,
                                       n_meshes), //species first to mesh first
          MkVec<int,    int   >(mesh_env, n_meshes),
          MkVec<double,    double>(mesh_vol, n_meshes),


          MkVec<double, double>(k, n_reactions),
          MkVec<double, int   >(sub, n_species*n_reactions),
          MkVec<double, int   >(sto, n_species*n_reactions),
          MkVec<double, int   >(r_env, n_reactions*n_env),
          MkVec<double, double>(D, n_species*n_env),

          sample_n,
          MkVec<double, double>(sample_t, sample_n),

          sampling_policy_code,
          sampling_interval,
          t_max,

          time_step,
          seed
          );

    return 0;
    }

extern "C" int Run(int breathe_dt)
    {
    bool unfinished = true;
    auto t0 = std::chrono::system_clock::now();
    for(;;)
        {
        if      (global_space_type == 0) unfinished = global_grid_algo->Iterate();
        else if (global_space_type == 1) unfinished = global_graph_algo->Iterate();
        int dt = static_cast<int>(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - t0).count());
        if(!unfinished || dt>=breathe_dt)
            break;
        }
    return unfinished;
    }

extern "C" int IterateN(int n_iterations)
    {
    bool unfinished = true;
    for(int i=0; i<n_iterations; i++)
        {
        if      (global_space_type == 0) unfinished = global_grid_algo->Iterate();
        else if (global_space_type == 1) unfinished = global_graph_algo->Iterate();
        if(!unfinished)
            break;
        }
    return unfinished;
    }

extern "C" int Iterate()
    {
    bool unfinished = true;
    if      (global_space_type == 0) unfinished = global_grid_algo->Iterate();
    else if (global_space_type == 1) unfinished = global_graph_algo->Iterate();
    return unfinished;
    }

extern "C" double GetProgress()
    {
    //return t/tmax
    double progress=0;
    if      (global_space_type == 0) progress = global_grid_algo->GetProgress();
    else if (global_space_type == 1) progress = global_graph_algo->GetProgress();
    return progress;
    }

extern "C" int GetOutput(double * trajectory_data)
    {
    if (global_space_type == 0)
      {

      int n_samples = global_grid_algo->NSamples();
      int n_species = global_grid_algo->NSpecies();
      int n_meshes  = global_grid_algo->NMeshes();

      std::vector<std::vector<double>> & trajectory_data_vec = global_grid_algo->GetSampledStates();
      for(int n=0;n<n_samples;n++)
          {
          for(int s=0;s<n_species;s++)
              {
              for(int i=0; i<n_meshes; i++)
                  {
                  //mesh first to species first
                  trajectory_data[n*n_meshes*n_species+ s*n_meshes + i] = trajectory_data_vec[n][i*n_species+s];
                  }
              }
          }
      return 0;
      }
    else
      {

      int n_samples = global_graph_algo->NSamples();
      int n_species = global_graph_algo->NSpecies();
      int n_meshes  = global_graph_algo->NMeshes();

      std::vector<std::vector<double>> & trajectory_data_vec = global_graph_algo->GetSampledStates();
      for(int n=0;n<n_samples;n++)
          {
          for(int s=0;s<n_species;s++)
              {
              for(int i=0; i<n_meshes; i++)
                  {
                  //mesh first to species first
                  trajectory_data[n*n_meshes*n_species+ s*n_meshes + i] = trajectory_data_vec[n][i*n_species+s];
                  }
              }
          }
      return 0;
      }
    }

extern "C" int GetState(double * state_data)
    {
    if (global_space_type == 0)
      {
      int n_species = global_grid_algo->NSpecies();
      int n_meshes  = global_grid_algo->NMeshes();

      std::vector<double> & state_data_vec = global_grid_algo->GetState();

      for(int s=0;s<n_species;s++)
          {
          for(int i=0; i<n_meshes; i++)
              {
              //mesh first to species first
              state_data[s*n_meshes + i] = state_data_vec[i*n_species+s];
              }
          }

      return 0;
      }
    else
      {
      int n_species = global_graph_algo->NSpecies();
      int n_meshes  = global_graph_algo->NMeshes();

      std::vector<double> & state_data_vec = global_graph_algo->GetState();

      for(int s=0;s<n_species;s++)
          {
          for(int i=0; i<n_meshes; i++)
              {
              //mesh first to species first
              state_data[s*n_meshes + i] = state_data_vec[i*n_species+s];
              }
          }

      return 0;
      }
    }

extern "C" double GetT()
    {
    if (global_space_type == 0)
      return global_grid_algo->GetT();
    else
      return global_graph_algo->GetT();
    }

extern "C" int GetTSample(double * t_sample)
    {
    if (global_space_type == 0)
      {
      std::vector<double> & t_sample_vec = global_grid_algo->GetSampledT();
      int n_sample = global_grid_algo->NSamples();

      for(int i=0; i<n_sample ; i++)
          {
          //mesh first to species first
          t_sample[i] = t_sample_vec[i];
          }

      return 0;
      }
    else
      {
      std::vector<double> & t_sample_vec = global_graph_algo->GetSampledT();
      int n_sample = global_graph_algo->NSamples();

      for(int i=0; i<n_sample ; i++)
          {
          //mesh first to species first
          t_sample[i] = t_sample_vec[i];
          }

      return 0;
      }

    }

extern "C" int GetNSamples()
    {
    if (global_space_type == 0)
      return global_grid_algo->NSamples();
    else
      return global_graph_algo->NSamples();
    }

extern "C" int Sample()
    {
    if (global_space_type == 0)
      global_grid_algo->Sample();
    else
      global_graph_algo->Sample();
    return 0;
    }

extern "C" int Finalize ()
    {
    if(global_algo_freed)
      return 0;

    if (global_space_type == 0)
      delete global_grid_algo;
    else
      delete global_graph_algo;

    return 0;
    }
