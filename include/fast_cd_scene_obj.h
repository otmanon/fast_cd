#pragma once
#include "fast_cd_arap_sim.h"
#include "fast_cd_viewer_custom_shader.h"
#include "fast_cd_subspace.h"
#include "read_rig_from_json.h"
#include "surface_to_volume_weights.h"
#include "fit_rig_to_mesh.h"
#include "read_rig_anim_from_json.h"
#include "rig_parameters.h"

/*
Helper class for one object of many in the fast cd scene
*/
struct fast_cd_scene_obj
{
	fast_cd_arap_sim sim;
    fast_cd_subspace sub;
    cd_sim_state st;

	bool do_cd;
    bool loop;
    int anim_length;
	MatrixXd anim_P;
    MatrixXd W; //control rig weights
    MatrixXd Ws; //secondary rig weights

    MatrixXd P0; //initial world space rig positions (for rendering/vis)
    VectorXd pL; //rig bone lengths (for rendering/vis)

    MatrixXd V;
    MatrixXi T, F;
    int timestep; //which timestep of the simulation is this obejct on.
    std::string name;
    

    //Texutre info if applicable
    bool do_texture; 
    MatrixXd V_tex, TC, N;
    MatrixXi F_tex, FTC, FN;
    SparseMatrix<double> P;

    fast_cd_scene_obj(std::string name, const MatrixXd& V, const MatrixXi& T, std::string rig_path, std::string rig_anim_path, fast_cd_subspace& sub, fast_cd_arap_sim& sim)
	{
        this->name = name;
		this->sim = sim;
        this->sub = sub;
		do_cd = true;
        loop = true;
        timestep = 0;
        this->V = V;
        this->T = T;
        this->Ws = sub.W;
        boundary_facets(T, F);

        //read rig from .json file. Only care about weights, P0, and PL
        VectorXi pI;
        MatrixXd Vs;
        MatrixXi Fs;
        std::string rig_type;
        read_rig_from_json(rig_path, W, P0, pI, pL, Vs, Fs, rig_type);
        int d = P0.cols();
        int b = P0.rows() / (d + 1);

        // map surface weights to volume
        MatrixXd rA(3, 4);
        if (rig_type == "surface")  
        {
         //   printf("Made it before surface_to_volume_weights call!\n");
            W = surface_to_volume_weights(W, Vs, V, T);
           // printf("Made it past surface_to_volume_weights call!\n");
            fit_rig_to_mesh(V, T, Vs, P0, rA);
        }
        else if (rig_type == "volume")
        {
            fit_rig_to_mesh_vertices(V, Vs, P0, rA);
        }
        else
        {
            printf("Rig type is not supported \n");
            exit(1);
        }

        MatrixXd anim_P_w;
        read_rig_anim_from_json(rig_anim_path, anim_P_w);
        transform_rig_parameters_anim(anim_P_w, rA);
        VectorXd p0 = (Map<VectorXd>(P0.data(), b * (d + 1) * d, 1)).cast<double>();
        world_to_rel_rig_anim(anim_P_w, p0, anim_P);

        VectorXd z; VectorXd p;
        p = anim_P.col(0);
        z = VectorXd::Zero(sub.W.cols() * 12);
        st.init(z, z, p, p);

        anim_length = anim_P.cols();

        do_texture = false;
	}
	

    virtual void transform_animation(MatrixXd& A)
    {
        transform_rig_parameters_anim(anim_P, A);
    }
    /*
    Steps the fast cd animation for this object forward in time
    */
    virtual void step(VectorXd& p, VectorXd& z)
    {   //get animation rig parameters
        p = anim_P.col(timestep % anim_length);
        z = VectorXd::Zero(sub.W.cols()*12);
        if (loop)
        {
            if (timestep % anim_length == 0)
            {
                st.init(z, z, p, p);
            }
        }
        if (do_cd)
        {
            VectorXd f_ext = VectorXd::Zero(sub.W.cols()*12); //  ii.momentum_leaking_force* reduced_sim->params.invh2* BMDJ* (2.0 * state.p_curr - state.p_prev - sol_p);           
            VectorXd bc; //no boundary conditions
            z = st.z_curr;
            z = sim.step(z, p, st, f_ext, bc);
            st.update(z, p);
        }
        timestep += 1;
    }



};