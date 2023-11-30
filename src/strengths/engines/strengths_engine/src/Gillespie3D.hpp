//implementation using the Gillespie algorithm

// #######################################################################################
// the Gillespie algorithm.
// reference :
// Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical reactions.
// The Journal of Physical Chemistry, 81(25), 2340-2361. https://doi.org/10.1021/j100540a008
// #######################################################################################

class Gillespie3D : public SimulationAlgorithm3DBase
    {
    private :

    std::vector<double> mesh_ar; //reaction propensities
    std::vector<double> mesh_ad; //diffusion porpensities
    std::vector<double> mesh_a0r; //a0
    std::vector<double> mesh_a0d; //a0d
    double a0;

    void ComputePropensities()
        {
        a0 = 0; //
        for(int i=0; i<n_meshes; i++)
            {
            //mesh reinitialize propensities sums
            mesh_a0d[i] = 0;
            mesh_a0r[i] = 0;

            //reaction rates
            for(int r=0; r<n_reactions; r++)
              {
              mesh_ar[i*n_reactions+r] = ReactionProp(i, r);
              mesh_a0r[i] += mesh_ar[i*n_reactions+r];
              a0 += mesh_ar[i*n_reactions+r];
              }

            for(int s=0; s<n_species; s++)
              {
              //diffusion
              for (int n=0; n<6; n++)
                {
                if(mesh_neighbors[i*6+n] != -1)
                  mesh_ad[i*6*n_species+s*6+n] = DiffusionProp(i, s, n);
                else
                  mesh_ad[i*6*n_species+s*6+n] = 0;

                mesh_a0d[i] += mesh_ad[i*6*n_species+s*6+n];
                a0 += mesh_ad[i*6*n_species+s*6+n];
                }
              }
            }
        }

    void ApplyReaction(int mesh_index, int reaction_index)
        {
        for(int s=0; s<n_species; s++)
            {
            mesh_x[mesh_index*n_species+s] += sto[s*n_reactions+reaction_index];
            }
        }

    void ApplyDiffusion(int mesh_index, int species_index, int direction)
        {
        int j = mesh_neighbors[mesh_index*6+direction];
        mesh_x[mesh_index*n_species+species_index] -= 1;
        mesh_x[j*n_species+species_index] += 1;
        }

    void DrawAndApplyEvent()
        {
        double r = uiud(rng)*a0;
        double a0_cumul = 0;
        for(int i=0; i<n_meshes; i++)
            {
            if(r < a0_cumul + mesh_a0r[i])
                {
                //reaction
                double r2 = r - a0_cumul;
                double a_cumul = 0;
                for(int j=0; j<n_reactions; j++)
                    {
                    a_cumul += mesh_ar[i*n_reactions+j];
                    if(r2<a_cumul)
                        {
                        ApplyReaction(i, j);
                        break;
                        }
                    }
                break;
                }
            a0_cumul += mesh_a0r[i];

            if(r < a0_cumul + mesh_a0d[i])
                {
                //diffusion
                double r2 = r - a0_cumul;
                double a_cumul = 0;
                bool diff_is_done = false;
                for(int j=0; j<n_species; j++)
                    {
                    for(int n=0; n<6; n++)
                        {
                        a_cumul += mesh_ad[i*n_species*6+j*6+n];
                        if(r2<a_cumul)
                            {
                            ApplyDiffusion(i, j, n);
                            diff_is_done = true;
                            break;
                            }
                        }
                      if(diff_is_done) break;
                      }
                break;
                }
            a0_cumul += mesh_a0d[i];
            }
        }

    virtual void AlgorithmSpecificInit()
        {
        this->mesh_ar.resize(n_reactions*n_meshes);
        this->mesh_ad.resize(6*n_species*n_meshes);
        this->mesh_a0r.resize(n_meshes);
        this->mesh_a0d.resize(n_meshes);
        }

    public :

    Gillespie3D()
        {
        }

    virtual ~Gillespie3D()
        {
        }

    virtual bool Iterate()
        {
        sampling_done_this_iteration = false; // reset the flag

        if(complete)
          return false;

        ComputePropensities();
        if(a0 == 0)
            {
            FlagAsComplete();
            //CompleteSampling();
            }
        else
            {
            DrawAndApplyEvent();
            dt = log(1/uiud(rng))/a0;
            t += dt;
            SamplingStep();
            CheckTMax();
            }
        return !complete;
        }
    };
