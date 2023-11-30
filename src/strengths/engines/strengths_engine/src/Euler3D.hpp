//implementation using the Euler method

class Euler3D : public SimulationAlgorithm3DBase
    {
    private :

    std::vector<double> mesh_dxdt; //species quantities

    void Compute_dxdt()
        {
        for(int i=0; i<n_meshes; i++)
            {
            //reaction rates
            std::vector<double> rr(n_reactions);

            for(int r=0; r<n_reactions; r++)
                rr[r] = ReactionRate(i, r);

            for(int s=0; s<n_species; s++)
              {
              mesh_dxdt[i*n_species+s] = 0;
              if(mesh_chstt[i*n_species+s]) continue;

              //reaction
              for(int r=0; r<n_reactions; r++)
                {
                mesh_dxdt[i*n_species+s] += sto[s*n_reactions+r]*rr[r];
                }

              //diffusion
              for (int n=0; n<6; n++)
                {
                if(mesh_neighbors[i*6+n] != -1)
                  mesh_dxdt[i*n_species+s] -= DiffusionRateDifference(i, s, n);
                }
              }
            }
        }

    void Apply_dxdt()
        {
        for(int i=0; i<n_meshes; i++)
            {
            for(int j=0; j<n_species; j++)
                {
                mesh_x[i*n_species+j] += mesh_dxdt[i*n_species+j]*dt;
                }
            }
        }

    virtual void AlgorithmSpecificInit()
        {
        this->mesh_dxdt.resize(n_species*n_meshes);
        }

    public :

    Euler3D()
        {
        }

    virtual ~Euler3D()
        {
        }

    virtual bool Iterate()
        {
        sampling_done_this_iteration = false; // reset the flag

        if(complete)
          return false;

        Compute_dxdt();
        Apply_dxdt();
        t += dt;
        SamplingStep();
        CheckTMax();
        return !complete;
        }

    };
