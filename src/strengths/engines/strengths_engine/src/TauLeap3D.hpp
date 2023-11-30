//implementation using the Gillespie algorithm with the tau-leap approximation

// #######################################################################################
// the tau leap algorithm.
// reference :
// Gillespie, D. T. (2001). Approximate accelerated stochastic simulation of chemically reacting systems.
// The Journal of Chemical Physics, 115(4), 1716-1733. https://doi.org/10.1063/1.1378322
// #######################################################################################

class TauLeap3D : public SimulationAlgorithm3DBase
    {
    private :

    std::vector<int> mesh_nr; //species quantities
    std::vector<int> mesh_nd; //species quantities

    void Compute_nevt()
        {
        for(int i=0; i<n_meshes; i++)
            {
            //reaction rates
            for(int r=0; r<n_reactions; r++)
                mesh_nr[i*n_reactions+r] = Poisson(ReactionProp(i, r)*dt);

            for(int s=0; s<n_species; s++)
                {
                //diffusion
                for (int n=0; n<6; n++)
                  {
                  if(mesh_neighbors[i*6+n] != -1)
                    mesh_nd[i*6*n_species+s*6+n] = Poisson(DiffusionProp(i, s, n)*dt);
                  else
                    mesh_nd[i*6*n_species+s*6+n] = 0;
                  }
                }
            }
        }

    void Apply_nevt()
        {
        for(int i=0; i<n_meshes; i++)
            {
            for(int r=0; r<n_reactions; r++)
              {
              for(int j=0; j<n_species; j++)
                {
                if(mesh_chstt[i*n_species+j]) continue;
                mesh_x[i*n_species+j] += sto[j*n_reactions+r]*mesh_nr[i*n_reactions+r];
                }
              }

            for(int s=0; s<n_species; s++)
                {
                for (int n=0; n<6; n++)
                    {
                    if(mesh_nd[i*6*n_species+s*6+n]==0) continue;

                    if(! mesh_chstt[i*n_species+s])
                        {
                        mesh_x[i*n_species+s] -= mesh_nd[i*6*n_species+s*6+n];
                        }
                    int j = mesh_neighbors[i*6+n];
                    if(! mesh_chstt[j*n_species+s])
                        {
                        mesh_x[j*n_species+s] += mesh_nd[i*6*n_species+s*6+n];
                        }
                    }
                }
            }
        }

    virtual void AlgorithmSpecificInit()
        {
        this->mesh_nr.resize(n_reactions*n_meshes);
        this->mesh_nd.resize(6*n_species*n_meshes);
        }

    public :

    TauLeap3D()
        {
        }

    virtual ~TauLeap3D()
        {
        }

    virtual bool Iterate()
        {
        sampling_done_this_iteration = false; // reset the flag
        
        if(complete)
          return false;

        Compute_nevt();
        Apply_nevt();
        t += dt;
        SamplingStep();
        CheckTMax();
        return !complete;
        }

    };
