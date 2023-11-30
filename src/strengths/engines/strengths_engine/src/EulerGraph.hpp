//implementation using the Euler method
//in a graph space

class EulerGraph : public SimulationAlgorithmGraphBase
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
              for (int n=0; n<mesh_neighbor_n[i]; n++)
                {
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

    EulerGraph()
        {
        }

    virtual ~EulerGraph()
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
