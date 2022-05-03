#include "InteractiveCDHook.h"

#include "arap_hessian.h"
#include "kmeans.h"
#include "augment_with_linear_constraints.h"
#include "momentum_leaking_matrix.h"
#include "sparse_diag.h"
#include "covariance_scatter_matrix.h"
#include "interweaving_matrix.h"
#include "grouping_matrix_from_clusters.h"

#include <chrono>
#include <cassert>

#include <filesystem>

#include "igl/repdiag.h"
#include "igl/cat.h"
#include <igl/slice.h>
#include <igl/average_onto_faces.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include <igl/slice_into.h>
#include <igl/find.h>

#include <write_sparse_ijv_DMAT.h>

#include <igl/get_seconds.h>
bool InteractiveCDHook::simulateOneStep()
{
    //Q += dt * Eigen::MatrixXd::Ones(Q.rows(), Q.cols());
    auto start = std::chrono::high_resolution_clock::now();

    //check if any simulation change events have occured
    if (refresh)
    { 
        poll_sim_changes();
        refresh = !refresh;
    }
    if (as.animation_mode == ANIMATION_MODE::INTERACTIVE_ANIMATION)
    {
       
        if (as.constraint_type == CONSTRAINT_TYPE::PINNING)
        {
            if (as.do_reduction)
            {//
                reduced_sim_step_pinning_control();
            }
            else
            {
                full_sim_step_pinning_control();
            }
        }
        else
        {
         //   as.constraint_controller->get_interactive_motion(bc);
  
            if (as.do_reduction)
            {//
                reduced_sim_step_cd_control();
            }
            else
            {
                full_sim_step_cd_control();
            }
        }
        step += 1;

        as.rig_controller->set_scripted_motion(step);  // change our rig controller parameters if they're scripted
    }
    if (as.animation_mode == ANIMATION_MODE::EIGENMODES_ANIMATION)
    {
        sim_step_modal_animation();
    }

    auto end = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1e6;

    int timing_i =step % timing_range;
    timings[timing_i] = duration;

    if (as.record_metrics)
    {
        if (step == as.record_timelength)
        {
            save_results();
        }
    }
    //  std::cout << "Sim FPS : " + std::to_string(duration) << std::endl;
    //  std::cout << "Sim FPS : " + std::to_string(1./timings.mean()) << std::endl;
    return false;
}



void InteractiveCDHook::full_sim_step_pinning_control()
{
    Eigen::VectorXd bc;
    as.rig_controller->get_pinned_motion(V0, rig->bI, bc);
    Eigen::VectorXd u_next;
    if (as.record_metrics)
    {
        double energy, grad_norm, start, time;
        Eigen::VectorXd u_null = sim.full_step(u_curr, u_prev, bc, energy, grad_norm);  //get the energy
        start = igl::get_seconds(); //call this again to get timing, method above takes more time to record.
        u_next = sim.full_step(u_curr, u_prev, bc);
        time = igl::get_seconds() - start;
        time_energy.conservativeResize(time_energy.rows() + 1, 2);
        time_energy.row(time_energy.rows() - 1) = Eigen::RowVector2d(time, energy);
    }
    else
    {
        u_next = sim.full_step(u_curr, u_prev, bc);
    }

    u_prev = u_curr;
    u_curr = u_next;

 //   Eigen::MatrixXd U_next = Eigen::Map<Eigen::MatrixXd>(u_next.data(), u_next.rows() / 3, 3);

 //   V = V0;
 //   if (v_state.vis_cd)
 //   {
 //       V += U_next;
 //   }
 //   if (v_state.vis_mode == TEXTURES) //visualize high res embedded mesh.
//    {
 //       Eigen::VectorXd v_high_res = W_low_to_high * u_next;
 //       V_high_res = V_high_res0 + Eigen::Map<Eigen::MatrixXd>(v_high_res.data(), v_high_res.rows() / 3, 3);
 //   }

    igl::slice(V, ext_ind, 1, V_ext);
}

void InteractiveCDHook::reduced_sim_step_pinning_control()
{
    Eigen::VectorXd bc;
    as.rig_controller->get_pinned_motion(V0, rig->bI, bc);
    Eigen::VectorXd z_next;
    if (as.record_metrics)
    {
        double energy, grad_norm, start, time;
        Eigen::VectorXd z_null = sim.reduced_step(z_curr, z_prev, bc, energy, grad_norm);  //get the energy
        start = igl::get_seconds(); //call this again to get timing, method above takes more time to record.
        z_next = sim.reduced_step(z_curr, z_prev, bc);
        Eigen::VectorXd u_null = pinned_B_ext * z_next;   //gotta include this step so that it's fair
        time = igl::get_seconds() - start;
        time_energy.conservativeResize(time_energy.rows() + 1, 2);
        time_energy.row(time_energy.rows() - 1) = Eigen::RowVector2d(time, energy);
    }
    else
    {
        z_next = sim.reduced_step(z_curr, z_prev, bc);
    }

    z_prev = z_curr;
    z_curr = z_next;
    /*

    V_ext = V0_ext;
    if (v_state.vis_cd)
    {
        //get exterior surface vertex quantities only.
        const Eigen::VectorXd u_ext = pinned_B_ext * z_next;
        const Eigen::MatrixXd U_ext = Eigen::Map<const Eigen::MatrixXd>(u_ext.data(), u_ext.rows() / 3, 3);

        V_ext += U_ext;
    }

    if (v_state.vis_mode == TEXTURES) //visualize high res embedded mesh.
    {
        V_high_res = V_high_res0;
        if (v_state.vis_cd)
        {
            Eigen::VectorXd u_high_res = WB * z_next;
            V_high_res += Eigen::Map<Eigen::MatrixXd>(u_high_res.data(), u_high_res.rows() / 3, 3);
        }
    }*/
}


void InteractiveCDHook::full_sim_step_cd_control()
{
    Eigen::VectorXd  p_next = as.rig_controller->p_rel;
    Eigen::VectorXd uc_next = cd_sim.full_step(p_next, p_curr, p_prev, uc_curr, uc_prev);
    uc_prev = uc_curr;
    uc_curr = uc_next;

    p_prev = p_curr;
    p_curr = p_next;

    Eigen::MatrixXd Uc_next = Eigen::Map<Eigen::MatrixXd>(uc_next.data(), uc_next.rows() / 3, 3);

  //  Eigen::VectorXd r = rig->J * p_next;
  //  Eigen::MatrixXd R = Eigen::Map<Eigen::MatrixXd>(r.data(), r.rows() / 3, 3);
  //  V = R;
    if (v_state.vis_cd)
    {
  //      V += Uc_next;
    }
    igl::slice(V, ext_ind, 1, V_ext);
    if (v_state.vis_mode == TEXTURES) //visualize high res embedded mesh.
    {
 //       Eigen::VectorXd v = W_low_to_high * (uc_next + r);
 //       V_high_res = Eigen::Map<Eigen::MatrixXd>(v.data(), v.rows() / 3, 3);
    }
}

void InteractiveCDHook::reduced_sim_step_cd_control()
{
   //nice and simple
    //local_global_solver_reducedh(z, z_next);
    p_next = as.rig_controller->p_rel;
    z_next = cd_sim.reduced_step(p_next, p_curr, p_prev, z_curr, z_prev);

    z_prev = z_curr;
    z_curr = z_next;

    p_prev = p_curr;
    p_curr = p_next;

   // if (v_state.vis_mode == TEXTURES)
   // {
   //     Eigen::VectorXd v_fine, v_coarse;
   //     if (v_state.vis_cd)
   //     {
   //         v_coarse = (cd_sim.B * z_next + rig->J * p_curr);
   //         v_fine = W_low_to_high * v_coarse; // a little slow... some things can be sped up, whatevs
   //     }
   //     else
   //     {
   //         v_coarse = ( rig->J * p_curr);
   //         v_fine = W_low_to_high * v_coarse;
   //     }
   //     V_high_res = Eigen::Map<Eigen::MatrixXd>(v_fine.data(), v_fine.rows() / 3, 3);
   //     V = Eigen::Map<Eigen::MatrixXd>(v_coarse.data(), v_coarse.rows() / 3, 3);
   //     igl::slice(V, ext_ind, 1, V_ext); // dont forget to update cage
   // }
   // else
   // {
  //      const Eigen::VectorXd r_ext = cd_J_ext * p_next;
  //      Eigen::MatrixXd R_ext = Eigen::Map<const Eigen::MatrixXd>(r_ext.data(), r_ext.rows() / 3, 3);

  //      V_ext = R_ext; 
        if (v_state.vis_cd)  //don't slice into it if we arent visualizing cd
        {
            //get exterior surface vertex quantities only.
  //          const Eigen::VectorXd uc_ext = cd_B_ext * z_next;
  //          const Eigen::MatrixXd Uc_ext = Eigen::Map<const Eigen::MatrixXd>(uc_ext.data(), uc_ext.rows() / 3, 3);

  //          V_ext += Uc_ext;
            //igl::slice_into(Uc_ext, ext_ind, 1, V_ext);
        }
    //}
}

void InteractiveCDHook::sim_step_modal_animation()
{
    step += 1;
    if (step % mas.period == 0)
    {
        mas.mode += 1;
        mas.mode = mas.mode % as.r;

        V = V0;
    }

    int progress = step - mas.period * (step / mas.period);
    if (progress == mas.half_period)
    {
        mas.scale *= -1;
    }
    Eigen::VectorXd B_flat;
    if (as.constraint_type == COMPLEMENTARY_DYNAMICS)
    {
        B_flat = cd_sim.B.col(mas.mode);
    }
    else if (as.constraint_type == PINNING)
    {
        B_flat = sim.B.col(mas.mode);
    }


   // V += mas.scale * Eigen::Map<Eigen::MatrixXd>(B_flat.data(), V.rows(), 3);
   // igl::slice(V, ext_ind, 1, V_ext);

    if (v_state.vis_mode == TEXTURES) //visualize high res embedded mesh.
    {
   //     Eigen::VectorXd v = W_low_to_high * Eigen::Map<Eigen::VectorXd>(V.data(), V.rows() * 3);
   //     V_high_res = Eigen::Map<Eigen::MatrixXd>(v.data(), v.rows() / 3, 3);
    }
}