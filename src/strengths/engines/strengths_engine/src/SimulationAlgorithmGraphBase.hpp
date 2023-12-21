//implements the base for the kinetics simulation algorithms in graph space

class SimulationAlgorithmGraphBase
    {
    protected :

    int n_meshes, n_edges;                          // number of meshes and mumber of edges
    int n_species, n_reactions, n_env;              // number of species, number of reactions// N=n_species, M=n_reactions

    std::vector<double> mesh_x;                     // species quantities
    std::vector<int>    mesh_chstt;                 // meshes chemostates
    std::vector<int>    mesh_env;                   // meshes environments
    std::vector<double> mesh_vol;                   // meshes volumes

    std::vector<int> mesh_neighbor_n;           // number of neighbors for each mesh
    std::vector<std::vector<int>> mesh_neighbor_index;      // neighbor mesh index
    std::vector<std::vector<double>> mesh_neighbor_sfc;      // neighbor contact suface
    std::vector<std::vector<double>> mesh_neighbor_dst;      // neighbor center-center distance

    std::vector<std::vector<double>> mesh_kd_out;      // diffusion rate constant to neighbor meshes
    std::vector<std::vector<double>> mesh_kd_in;      // diffusion rate constant from neighbor meshes

    std::vector<double> sto;                        // reaction species change stoechiometry matrix
    std::vector<double> sub;                        // substrates stroechiometry matrix
    std::vector<double> mesh_kr;                    // reaction kinetic rates (accounting for mesh volumes)//n_meshes * n_reactions

    int n_samples;                                  // size of t_samples
    int sample_pos;                                 // index of the next sample time
    std::vector<double> t_samples;                  // timepoints at which the system state should be sampled

    std::vector<std::vector<double>> sampled_mesh_x; // sampled system states.
    std::vector<double> sampled_t;                   // exact time at which each system states in sampled_mesh_x were sampled.

    int sampling_policy_code;                        // see Init arguments
    double sampling_interval;                        // see Init arguments
    double t_max;                                    // see Init arguments
    bool sampling_done_this_iteration;               // flags that tells if sample() have been called for the current iteration.
    double last_tsi_ratio;                           // floored value of the t/sampling_interval used for on_interval sampling.

    double t;                                        // time
    double dt;                                       // time step
    bool complete;                                   // true if all the sampling is done.
    std::mt19937 rng;                                // pseudo random number generator
    std::uniform_real_distribution<double> uiud;     // floating point uniform distribution in [0,1[


    int Poisson(double lambda)
        {
        return std::poisson_distribution<int>(lambda)(rng);
        }

    void SetNeighbors(
          int n_edges,
          std::vector<int>&edge_i,
          std::vector<int>&edge_j,
          std::vector<double>&edge_sfc,
          std::vector<double>&edge_dst)
      {
        mesh_neighbor_n.resize(n_meshes, 0);
        mesh_neighbor_index.resize(n_meshes);
        mesh_neighbor_sfc.resize(n_meshes);
        mesh_neighbor_dst.resize(n_meshes);

      for(int i=0; i<n_edges; i++)
          {
          mesh_neighbor_n[edge_i[i]]++;
          mesh_neighbor_n[edge_j[i]]++;

          mesh_neighbor_index[edge_i[i]].push_back(edge_j[i]);
          mesh_neighbor_index[edge_j[i]].push_back(edge_i[i]);

          mesh_neighbor_sfc[edge_i[i]].push_back(edge_sfc[i]);
          mesh_neighbor_sfc[edge_j[i]].push_back(edge_sfc[i]);

          mesh_neighbor_dst[edge_i[i]].push_back(edge_dst[i]);
          mesh_neighbor_dst[edge_j[i]].push_back(edge_dst[i]);
          }
      }

    void Build_mesh_kr(const std::vector<double> & k, const std::vector<double> & r_env)
    // builds mesh_kr
        {
        mesh_kr.clear();
        mesh_kr.resize(n_meshes*n_reactions, 0);
        for(int i=0;i<n_meshes;i++)
          {
          for(int r=0; r<n_reactions; r++)
            {
            double q = 0;
            for(int s=0; s<n_species; s++)
                q+=sub[s*n_reactions+r];
            mesh_kr[i*n_reactions+r] = k[r]*pow(mesh_vol[i],1-q)*r_env[r*n_env+mesh_env[i]];
            }
          }
        }

    void Build_mesh_kd(const std::vector<double> & D)
    // builds mesh_kd
        {
        mesh_kd_out.clear();
        mesh_kd_out.resize(n_meshes);

        mesh_kd_in.clear();
        mesh_kd_in.resize(n_meshes);

        for(int i=0;i<n_meshes;i++)
            {
            mesh_kd_out[i].resize(n_species*mesh_neighbor_n[i]);
            mesh_kd_in [i].resize(n_species*mesh_neighbor_n[i]);

            for(int s=0;s<n_species;s++)
                {
                for(int n=0; n<mesh_neighbor_n[i]; n++)
                    {
                    int j = mesh_neighbor_index[i][n];

                    // #########################################################################
                    // diffusion coefficient between cells i and j is calculated according to David Bernstein's method.
                    // reference :
                    // Bernstein, D. (2005). Simulating mesoscopic reaction-diffusion systems using the Gillespie algorithm.
                    // Physical Review E, 71(4), Article 041103. https://doi.org/10.1103/PhysRevE.71.041103
                    double hi = pow(mesh_vol[i], 1.0/3.0);
                    double hj = pow(mesh_vol[j], 1.0/3.0);

                    double Di = D[s*n_env+mesh_env[i]];
                    double Dj = D[s*n_env+mesh_env[j]];
                    double Dij = 0;

                    if(Di!=0 && Dj!=0)
                      {
                      Dij = (hi+hj)/(hi/Di + hj/Dj);
                      }
                    // #########################################################################
                    
                    // the way the diffusion rate is computed from Dij here have been slightly extended from what is presented in 
                    // [Bernstein, D. (2005). Simulating mesoscopic reaction-diffusion systems using the Gillespie algorithm.
                    // Physical Review E, 71(4), Article 041103. https://doi.org/10.1103/PhysRevE.71.041103], 
                    // in order to take into account the specific exchage surface between cells i and j
                    mesh_kd_out[i][s*mesh_neighbor_n[i]+n] = Dij * mesh_neighbor_sfc[i][n] / (mesh_vol[i] * mesh_neighbor_dst[i][n]);
                    mesh_kd_in [i][s*mesh_neighbor_n[i]+n] = Dij * mesh_neighbor_sfc[i][n] / (mesh_vol[j] * mesh_neighbor_dst[i][n]);
                    }
                }
            }
        }

    void CheckTMax()
      {
      if(t_max>=0 && t>t_max)
        {
        FlagAsComplete();
        }
      }

    void SampleOnTSample()
    // sample as close a possible from t_samples.
        {
        while(t>=t_samples[sample_pos] && sample_pos<n_samples)
            {
            Sample();
            sample_pos ++;
            }
        }

    void SampleOnInterval()
    // sample at the given interval sampling_interval.
        {
        double tsi_ratio = floor(t/sampling_interval);
        if(tsi_ratio > last_tsi_ratio)
          {
          Sample();
          last_tsi_ratio = tsi_ratio;
          }
        }

    void SamplingStep()
    // manage the sampling procedure according to the samplng policy.
        {
        switch(sampling_policy_code)
          {
          case 0 : SampleOnTSample(); break;  //sample on t sample
          case 1 : Sample(); break;           //sample on iteration
          case 2 : SampleOnInterval(); break; //sample on interval
          case 3 : break;                     //no sample
          };
        }

    void FlagAsComplete()
    // flag the simulation as complete
    // whatever the completion cause is
      {
      complete = true;
      }

    double ReactionRate(int mesh_index, int reaction_index)
    // computes the deterministic reaction rate
        {
        double r = mesh_kr[mesh_index*n_reactions+reaction_index];
        for(int s=0; s<n_species; s++)
            r *= pow(mesh_x[mesh_index*n_species+s], sub[s*n_reactions+reaction_index]);
        return r;
        }

    double ReactionProp(int mesh_index, int reaction_index)
    // computes the Gillespie reaction propensity
        {
        // #######################################################################################
        // reactions propensities for the Gillespie algorithm.
        // reference :
        // Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical reactions.
        // The Journal of Physical Chemistry, 81(25), 2340-2361. https://doi.org/10.1021/j100540a008
        double a = mesh_kr[mesh_index*n_reactions+reaction_index];
        for(int s = 0; s<n_species; s++)
            {
            if (mesh_x[mesh_index*n_species+s] >= sub[s*n_reactions+reaction_index])
                {
                for (int q=0;q<sub[s*n_reactions+reaction_index];q++)
                    {
                    a *= (mesh_x[mesh_index*n_species+s]-q);
                    }
                }
            else
                {
                a = 0;
                break;
                }
            }
        return a;
        // #######################################################################################
        }

    double DiffusionRate(int mesh_index, int species_index, int direction)
        {
        // #######################################################################################
        // for diffusion events treated as first order reactions, as described by David Bernstein :
        // reference :
        // Bernstein, D. (2005). Simulating mesoscopic reaction-diffusion systems using the Gillespie algorithm.
        // Physical Review E, 71(4), Article 041103. https://doi.org/10.1103/PhysRevE.71.041103
        return mesh_x[mesh_index*n_species+species_index] * mesh_kd_out[mesh_index][species_index*mesh_neighbor_n[mesh_index]+direction];
        // #######################################################################################
        }

    double DiffusionRateDifference(int mesh_index, int species_index, int direction)
        {
        return
        mesh_x[mesh_index*n_species+species_index] * mesh_kd_out[mesh_index][species_index*mesh_neighbor_n[mesh_index]+direction] -
        mesh_x[mesh_neighbor_index[mesh_index][direction]*n_species+species_index] * mesh_kd_in[mesh_index][species_index*mesh_neighbor_n[mesh_index]+direction];
        }

    double DiffusionProp(int mesh_index, int species_index, int direction)
        {
        // #######################################################################################
        // reactions propensities for the Gillespie algorithm :
        // reference :
        // Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical reactions.
        // The Journal of Physical Chemistry, 81(25), 2340-2361. https://doi.org/10.1021/j100540a008
        //
        // for diffusion events treated as first order reactions, as described by David Bernstein :
        // reference :
        // Bernstein, D. (2005). Simulating mesoscopic reaction-diffusion systems using the Gillespie algorithm.
        // Physical Review E, 71(4), Article 041103. https://doi.org/10.1103/PhysRevE.71.041103
        return mesh_x[mesh_index*n_species+species_index] * mesh_kd_out[mesh_index][species_index*mesh_neighbor_n[mesh_index]+direction];
        // #######################################################################################
        }

    virtual void AlgorithmSpecificInit() = 0;
    // should contain algorithm specific initialization steps, in order to avoid overwriting the constructor

    public :

    SimulationAlgorithmGraphBase()
        {
        }

    virtual ~SimulationAlgorithmGraphBase()
        {
        }

    void Init
    // constructor of the class
        (
        int    n_nodes,                  //number of meshes
        int    n_species,               //number of species
        int    n_reactions,             //number of reactions
        int    n_env,                   //number of encironments

        int n_edges,
        std::vector<int> edge_i,
        std::vector<int> edge_j,
        std::vector<double> edge_sfc,
        std::vector<double> edge_dst,


        std::vector<double> mesh_x0,    //initial state //mesh first array : [mesh [species]]
        std::vector<int>    mesh_chstt, //species chemostats //mesh first array : [mesh [species]]
        std::vector<int>    mesh_env,   //meshes environment indices
        std::vector<double> mesh_vol,   //meshes volume

        std::vector<double> k,          //reaction rates
        std::vector<double> sub,        //N*M substrate matrix
        std::vector<double> sto,        //N*M stoechiometry matrix
        std::vector<double> r_env,      //reactions environments
        std::vector<double> D,          //reactions diffusion coefficients

        int sample_n,                   //number of sample timepoints
        std::vector<double> t_samples,  //sample timepoints
        int sampling_policy_code,       //tells when the system state should be sampled
        double sampling_interval,       //interval at which the system state should be sampled (if sampling_policy_code=2)
        double t_max,                   //time past which the simulation should be flagged as complete (if negative, there is no t_max).
        double time_step,               //time step
        int seed                        //rng seed
        )
        {
        this->n_meshes = n_nodes;

        SetNeighbors(n_edges, edge_i, edge_j, edge_sfc, edge_dst);

        this->n_species = n_species;
        this->n_reactions = n_reactions;
        this->n_env = n_env;

        this->mesh_x = mesh_x0;
        this->mesh_chstt = mesh_chstt;
        this->mesh_env = mesh_env;
        this->mesh_vol = mesh_vol;

        this->sub = sub;
        this->sto = sto;
        this->n_samples = sample_n;
        this->t_samples = t_samples;
        this->sample_pos = 0;

        this->sampled_mesh_x.clear();
        this->sampled_t.clear();

        this->sampling_policy_code = sampling_policy_code;
        this->sampling_interval = sampling_interval;
        this->t_max = t_max;
        this->sampling_done_this_iteration = false;
        this->last_tsi_ratio = -1; // rather than 0, to allow for t0 sampling.

        this->t = 0.0;
        this->dt = time_step;
        this->complete = false;
        Build_mesh_kr(k, r_env);
        Build_mesh_kd(D);
        this->rng = std::mt19937(seed);
        this->uiud = std::uniform_real_distribution<double> (0.0, 1.0);
        this->AlgorithmSpecificInit();

        SamplingStep(); //for t0 sampling if necessary
        }

    virtual bool Iterate() = 0;
    // one itetation of the simulation algorithm. Returns true if the simulation should continue. False otherwise.

    double GetProgress()
    // returns 100*t/t_max
        {
        if(t_max>0) return 100.0 * t/t_max;
        else return 0;
        }

    std::vector<std::vector<double>> & GetSampledStates()
    // returns sample_mesh_xO
        {
        return sampled_mesh_x;
        }

    int NSamples()
    // return the number of states currently sampled,
    // not the n_sample value which was given as Init argument
        {
        return static_cast<int>(sampled_t.size());
        }

    int NSpecies()
        {
        return n_species;
        }

    int NMeshes()
        {
        return n_meshes;
        }

    std::vector<double> & GetState()
        {
        return mesh_x;
        }

    double GetT()
        {
        return t;
        }

    std::vector<double> & GetSampledT()
        {
        return sampled_t;
        }

    void Sample()
    // sample the current system state and time.
        {
        if(!sampling_done_this_iteration)
          {
          sampled_mesh_x.push_back(mesh_x);
          sampled_t.push_back(t);
          sampling_done_this_iteration = true;
          }
        }
    };
