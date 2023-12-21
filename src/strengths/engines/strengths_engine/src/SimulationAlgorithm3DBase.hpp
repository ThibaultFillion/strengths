//implements the base for the 3D kinetics simulation algorithms

class SimulationAlgorithm3DBase
    {
    protected :

    int w, h, d, n_meshes;                          // system dimensions and number of meshes
    int n_species, n_reactions, n_env;              // number of species, number of reactions// N=n_species, M=n_reactions
    std::vector<int> delta_i;                       // mesh index offset associated with each of the 6 directions
    std::vector<int> opposed_direction;             // opposed direction index
    std::vector<double> mesh_x;                     // species quantities
    std::vector<int> mesh_neighbors;                // index of mesh neighbor in each direction
    std::vector<int> mesh_chstt;                    // chemostates
    std::vector<int> mesh_env;                      // chemostates
    double mesh_vol;                                // mesh volume
    double mesh_edge;                               // mesh edge
    std::vector<double> sto;                        // reaction species change stoechiometry matrix
    std::vector<double> sub;                        // substrates stroechiometry matrix
    std::vector<double> mesh_kr;                    // reaction kinetic rates (accounting for mesh volumes)//n_meshes * n_reactions
    std::vector<double> mesh_kd;                    // diffusion kinetic rates //n_species * n_meshes * n_meshes

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
    std::vector<int> boundary_conditions;            // see Init arguments

    int Poisson(double lambda)
        {
        return std::poisson_distribution<int>(lambda)(rng);
        }

    bool AreNeighbors(int i, int j)
    // returns true if meshes i and j are neighbors. false ohterwise
        {
        int xi = i%w;
        int yi = i%(w*h)/w;
        int zi = i/(w*h);

        int xj = j%w;
        int yj = j%(w*h)/w;
        int zj = j/(w*h);

        return ((abs(xi-xj)+abs(yi-yj)+abs(zi-zj)) == 1);
        }

    int GetNeighborIndex(int x, int y, int z, int direction)
        {
        int xn = x;
        int yn = y;
        int zn = z;
        switch(direction)
            {
            case 0 : xn+=1; break;
            case 1 : xn-=1; break;
            case 2 : yn+=1; break;
            case 3 : yn-=1; break;
            case 4 : zn+=1; break;
            case 5 : zn-=1; break;
            };

        if (boundary_conditions[0] == 1) xn = (w+xn)%w;
        if (boundary_conditions[1] == 1) yn = (h+yn)%h;
        if (boundary_conditions[2] == 1) zn = (d+zn)%d;

        if (xn>=0 && xn<w &&
            yn>=0 && yn<h &&
            zn>=0 && zn<d)
              return w*h*zn + w*yn + xn;
        else return -1;
        }

    void BuildMeshNeighbors()
        {
        this->mesh_neighbors = std::vector<int>(w*h*d*6);
        for(int i=0; i<w*h*d; i++)
            {
            int xcoord = i%w;
            int ycoord = i%(w*h)/w;
            int zcoord = i/(w*h);

            for(int n=0; n<6; n++)
                {
                this->mesh_neighbors[i*6+n] = GetNeighborIndex(xcoord, ycoord, zcoord, n);
                }
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
            mesh_kr[i*n_reactions+r] = k[r]*pow(mesh_vol,1-q)*r_env[r*n_env+mesh_env[i]];
            }
          }
        }

    void Build_mesh_kd(const std::vector<double> & D)
    // builds mesh_kd
        {
        mesh_kd.clear();
        mesh_kd.resize(n_species*n_meshes*6, 0);
        for(int s=0;s<n_species;s++)
            {
            for(int i=0;i<n_meshes;i++)
                {
                for(int n=0; n<6; n++)
                    {
                    int j = mesh_neighbors[i*6+n];
                    if(j==-1)
                        {
                        mesh_kd[i*n_species*6 + s*6+ n] = 0;
                        continue;
                        }

                    // #########################################################################
                    // diffusion reaction rate constants are calculated according to David Bernstein's method.
                    // reference :
                    // Bernstein, D. (2005). Simulating mesoscopic reaction-diffusion systems using the Gillespie algorithm.
                    // Physical Review E, 71(4), Article 041103. https://doi.org/10.1103/PhysRevE.71.041103
                    double Di = D[s*n_env+mesh_env[i]];
                    double Dj = D[s*n_env+mesh_env[j]];
                    double Dij = 0;
                    if(Di!=0 && Dj!=0)
                      {
                      Dij = (2*mesh_edge)/(mesh_edge/Di+mesh_edge/Dj);
                      }

                    mesh_kd[i*n_species*6 + s*6+ n] = Dij/(mesh_edge*mesh_edge);
                    // #########################################################################
                    }
                }
            }
        }

    /*void CompleteSampling()
    // samples the current state for the remaining sample times.
    // The simulation is then flagged as complete, in case it is not already.
        {
        for(;;)
            {
            sampled_mesh_x.push_back(mesh_x);
            sample_pos ++;
            if (sample_pos == n_samples)
                {
                FlagAsComplete();
                break;
                }
            }
        }*/
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
        return mesh_x[mesh_index*n_species+species_index] * mesh_kd[mesh_index*n_species*6+species_index*6+direction];
        // #######################################################################################
        }

    double DiffusionRateDifference(int src_mesh_index, int species_index, int direction)
        {
        return DiffusionRate(src_mesh_index, species_index, direction) -
               DiffusionRate(mesh_neighbors[src_mesh_index*6+direction], species_index, opposed_direction[direction]);
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
        return mesh_x[mesh_index*n_species+species_index] * mesh_kd[mesh_index*n_species*6+species_index*6+direction];
        // #######################################################################################
        }

    virtual void AlgorithmSpecificInit() = 0;
    // should contain algorithm specific initialization steps, in order to avoid overwriting the constructor

    public :

    SimulationAlgorithm3DBase()
        {
        }

    virtual ~SimulationAlgorithm3DBase()
        {
        }

    void Init
    // constructor of the class
        (
        int w,                          //system width
        int h,                          //system height
        int d,                          //system depth
        int    n_species,               //number of species
        int    n_reactions,             //number of reactions
        int    n_env,                   //number of encironments
        std::vector<double> mesh_x0,    //initial state //mesh first array : [mesh [species]]
        std::vector<int>    mesh_chstt, //species chemostats //mesh first array : [mesh [species]]
        std::vector<int>    mesh_env,   //meshes environment indices
        double mesh_vol,                //volume of a mesh
        std::vector<double> k,          //reaction rates
        std::vector<double> sub,        //N*M substrate matrix
        std::vector<double> sto,        //N*M stoechiometry matrix
        std::vector<double> r_env,      //reactions environments
        std::vector<double> D,          //reactions diffusion coefficients
        std::vector<int> boundary_conditions,
        int sample_n,                   //number of sample timepoints
        std::vector<double> t_samples,  //sample timepoints
        int sampling_policy_code,       //tells when the system state should be sampled
        double sampling_interval,       //interval at which the system state should be sampled (if sampling_policy_code=2)
        double t_max,                   //time past which the simulation should be flagged as complete (if negative, there is no t_max).
        double time_step,               //time step
        int seed                        //rng seed
        )
        {
        this->boundary_conditions = boundary_conditions;
        this->w = w;
        this->h = h;
        this->d = d;
        this->delta_i = std::vector<int>{+1, -1, +w, -w, +(w*h), -(w*h)};
        BuildMeshNeighbors();
        this->opposed_direction = std::vector<int>{1, 0, 3, 2, 5, 4};
        this->n_meshes = w*h*d;
        this->n_species = n_species;
        this->n_reactions = n_reactions;
        this->n_env = n_env;
        this->mesh_x = mesh_x0;
        this->mesh_chstt = mesh_chstt;
        this->mesh_env = mesh_env;
        this->mesh_vol = mesh_vol;
        this->mesh_edge = pow(mesh_vol, 1.0/3.0);
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
