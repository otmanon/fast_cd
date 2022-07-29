#pragma once
#include "fast_cd_init.h"
#include "fast_complementary_dynamics_sim.h"
#include <Eigen/Core>
#include "sim_params.h"


#include "ScalarRecorder.h"
#include "VectorRecorder.h"
#include "MeshRecorder.h"
#include "screenshotter.h"

class Recorder
{
public:
	FastCDSim* sim;
	FastCDInit* init;
	
	ScalarRecorder m_kinetic_energy, m_bending_energy, m_volume_energy, m_kinetic_grad_norm, m_bending_grad_norm, m_volume_grad_norm, m_timings_proj, m_timings_sim, m_timings_basis_rotate;
	VectorRecorder m_z, m_kinetic_grad, m_bending_grad, m_volume_grad, rig_recorder;
	MeshRecorder mesh_recorder;
	sim_params s;
	//class that holds and records any scene data, metrics and is in charge of writing them to memory
	Recorder() {};
	Recorder(FastCDInit * init, FastCDSim* sim, int num_modes, int num_rig_params)
	{
		this->init = init;
		this->sim = sim;
		m_z = VectorRecorder(num_modes);
		m_kinetic_grad = VectorRecorder(num_modes);
		m_bending_grad = VectorRecorder(num_modes);
		m_volume_grad = VectorRecorder(num_modes);
		rig_recorder = VectorRecorder(num_rig_params);
	};

	void update_quantities(Eigen::VectorXd& z_next, Eigen::VectorXd& p_next, Eigen::MatrixXd& B, Eigen::SparseMatrix<double>& J, Eigen::MatrixXd& U, Eigen::MatrixXi& F, double dt_proj, double dt_sim, double dt_basis_rotate)
	{

        FastCDInit& ii = *init;
        FastCDSim& fcd_sim = *sim;
        //////////////////////////////////// Record Quantities ///////////////////////////////////////
        if (ii.record_energy)
        {   //calculating metrics we care about
            double inertia, bending, volume;
            if (ii.use_complementary_recordings)
            {
                if (ii.do_reduction)
                    fcd_sim.energy_complementary_reduced(z_next, s.z_curr, s.z_prev, bending, volume, inertia);
                else
                    fcd_sim.energy(z_next, s.z_curr, s.z_prev, bending, volume, inertia);
            }
            else
            {
                if (ii.do_reduction)
                    fcd_sim.energy_reduced(z_next, s.z_curr, s.z_prev, p_next, s.p_curr, s.p_prev, bending, volume, inertia);
                else
                {
                    Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(fcd_sim.X.data(), fcd_sim.X.rows() * fcd_sim.X.cols());
                    Eigen::VectorXd u_next = J * p_next + z_next - x;
                    Eigen::VectorXd u_curr = J * s.p_curr + s.z_curr - x;
                    Eigen::VectorXd u_prev = J * s.p_prev + s.z_prev - x;
                    fcd_sim.energy(u_next, u_curr, u_prev, bending, volume, inertia);
                }
            }
            m_kinetic_energy.record_frame(inertia);
            m_bending_energy.record_frame(bending);
            m_volume_energy.record_frame(volume);
        }
        if (ii.record_gradient_norm || ii.record_gradient)
        {
            Eigen::VectorXd inertia_grad, bending_grad, volume_grad;
            if (ii.use_complementary_recordings)
            {
                if (ii.do_reduction)
                    fcd_sim.grad_complementary_reduced(z_next, s.z_curr, s.z_prev, bending_grad, volume_grad, inertia_grad);
                else
                    fcd_sim.grad(z_next, s.z_curr, s.z_prev, bending_grad, volume_grad, inertia_grad);
            }
            else
            {
                if (ii.do_reduction)
                    fcd_sim.grad_reduced(z_next, s.z_curr, s.z_prev, p_next, s.p_curr, s.p_prev, bending_grad, volume_grad, inertia_grad);
                else
                {
                    Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(fcd_sim.X.data(), fcd_sim.X.rows() * fcd_sim.X.cols());
                    Eigen::VectorXd u_next = J * p_next + z_next - x;
                    Eigen::VectorXd u_curr = J * s.p_curr + s.z_curr - x;
                    Eigen::VectorXd u_prev = J * s.p_prev + s.z_prev - x;
                    fcd_sim.grad(u_next, u_curr, u_prev, bending_grad, volume_grad, inertia_grad);
                }
            }
            if (ii.record_gradient_norm)
            {
                m_kinetic_grad_norm.record_frame(inertia_grad.norm());
                m_bending_grad_norm.record_frame(bending_grad.norm());
                m_volume_grad_norm.record_frame(volume_grad.norm());
            }
            if (ii.record_gradient)
            {
                m_kinetic_grad.record_frame(inertia_grad);
                m_bending_grad.record_frame(bending_grad);
                m_volume_grad.record_frame(volume_grad);
            }
        }
        if (ii.record_modal_activations)
        {
            m_z.record_frame(z_next);
        }
        if (ii.record_timings)
        {
            m_timings_proj.record_frame(dt_proj);
            m_timings_sim.record_frame(dt_sim);
            m_timings_basis_rotate.record_frame(dt_basis_rotate);
        }
        if (ii.record_mesh)
        {
            mesh_recorder.record_frame(U, F);
        }
        if (ii.record_rig)
        {
            rig_recorder.record_frame(p_next);
        }

	}

	void save()
	{
      
        FastCDInit ii = *init;
        //////////////////////////////////// Record Quantities ///////////////////////////////////////
        if (ii.record_energy)
        {
            m_kinetic_energy.save(ii.results + "kinetic_energy.DMAT");
            m_bending_energy.save(ii.results + "bending_energy.DMAT");
            m_volume_energy.save(ii.results + "volume_energy.DMAT");
        }

        if (ii.record_gradient_norm)
        {
            m_kinetic_grad_norm.save(ii.results + "kinetic_grad_norm.DMAT");
            m_bending_grad_norm.save(ii.results + "bending_grad_norm.DMAT");
            m_volume_grad_norm.save(ii.results + "volume_grad_norm.DMAT");
        }
        if (ii.record_gradient)
        {
            m_kinetic_grad.save(ii.results + "kinetic_grad.DMAT");
            m_bending_grad.save(ii.results + "bending_grad.DMAT");
            m_volume_grad.save(ii.results + "volume_grad.DMAT");
        }
        if (ii.record_modal_activations)
        {
            m_z.save(ii.results + "modal_activations.DMAT");
        }
        if (ii.record_timings)
        {
            m_timings_proj.save(ii.results + "timings_proj.DMAT");
            m_timings_sim.save(ii.results + "timings_sim.DMAT");
            m_timings_basis_rotate.save(ii.results + "timings_basis_rotate.DMAT");
        }
        if (ii.record_mesh)
        {
            mesh_recorder.save(ii.results + "mesh_recordings/");
            //TODO add this
        }
        if (ii.record_rig)
        {
            rig_recorder.save(ii.results + "rig.DMAT");//ignore;
        }


	}
	

};