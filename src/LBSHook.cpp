#include "LBSHook.h"
#include <filesystem>
#include <stdio.h>
#include <cassert>
#include <chrono>

#include "arap_hessian.h"
#include "kmeans.h"
#include "compute_modes.h"
#include "augment_with_linear_constraints.h"
#include "momentum_leaking_matrix.h"
#include "sparse_diag.h"
#include "covariance_scatter_matrix.h"
#include "interweaving_matrix.h"
#include <normalize_geometry_unit_height.h>
#include "grouping_matrix_from_clusters.h"

#include "igl/repdiag.h"
#include "igl/cat.h"
#include <igl/min_quad_with_fixed.h>
#include "igl/polar_svd3x3.h"
#include <igl/average_onto_faces.h>
#include "igl/readMSH.h"
#include <igl/boundary_facets.h>
#include <igl/colon.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include "igl/slice.h"
#include <igl/unique.h>
#include <igl/png/readPNG.h>
#include <igl/unproject.h>
#include <igl/project.h>

LBSHook::LBSHook(std::string mesh_name) :PhysicsHook()
{
    r = 250;
    l = 250;

    dt = 1.0 / 60;        //60 fps, each timestep is 1/60 seconds
    k = 330000;     // young'as modulus... the units are definitely wrong.
    use_cd = true;
    use_inertia = true;
    do_reduction = true;
    do_clustering = true;

    this->mesh_name = mesh_name;
    load_file_paths();

    Eigen::VectorXi tritag, tettag;
    //read tet mesh. F is filled with junk after this
    igl::readMSH(data_dir + "/" + mesh_name + ".msh", V0, F, T, tritag, tettag);
    l = do_clustering ? l : T.rows();           //if we don't do clustering, then set num clusters to num tets
    
   // V0 = normalize_geometry_unit_height(V0);
    V = V0;
    Eigen::VectorXi _n; //throwaway var
    igl::boundary_facets(T, F, FiT, _n);

    //get list of exterior vertex indices
    igl::unique(F, ext_ind);

    //build our affine rig
    rig = LBSRig(rig_file_dir, V0, T, ext_ind);

   // Eigen::MatrixXd surfaceV;
    //igl::slice(V, ext_ind, 1, surfaceV);
    V = V0;
    //only update V_ext each timestep instead of all of V
   // igl::slice(V, ext_ind, 2, V_ext);

    //mouse control parameterts
    changed_constraint = false;
    moved_handle = false;
    mouse_dragging = false;

    //gui simulation change parameters
    new_stiffness = k;
    new_do_clustering = do_clustering;
    new_num_modes = r;
    new_num_clusters = l;

    //simulation parameters... make sure these are all initialized as they may persist when changing init_file_names
    init_simulation_state();



}

bool LBSHook::simulateOneStep()
{

    poll_sim_changes();

    rig.get_rig_motion(ur);
  //  x = Eigen::Map<Eigen::VectorXd>(rig.X.data(), X.rows() * X.cols());
    Eigen::MatrixXd U_next;
    Eigen::VectorXd z_next;

    local_global_solver_reduced(z, z_next);
    z = z_next;
    uc = B * z_next;
    u = B * z_next + ur;
    //update inertia history
    u_prev = u_curr;
    u_curr = u;


    U_next = Eigen::Map<Eigen::MatrixXd>(u.data(), X.rows(), X.cols());
    V = U_next + X;
    U = U_next;

  //  u = x + ur;
   // V = Eigen::Map<Eigen::MatrixXd>(u.data(), X.rows(), X.cols());
    
    return false;
}

void LBSHook::initSimulation()
{

    init_simulation_state();
   // rig.reset();

   
    //arap hessian
    arap_hessian(X, T, H);
    H *= -1;

    //cotan laplacian
    igl::cotmatrix(X, T, C);
    C = -igl::repdiag(C, 3);

    //mass matrix
    igl::massmatrix(X, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
    M = igl::repdiag(M, 3);


    //reset up point constraint info
    ui.resize(0); ub.resize(0);
    handleI.resize(0); handleV.resize(0, 3);
    currentHandle = -1;

    //constructs A matrix, and does preomputation if we will use it directly
    build_A_matrix();
    compute_cd_equality_constraints();

    //this line takes forever so only run it when we arent using a reduced system:
    if (!do_reduction)
        igl::min_quad_with_fixed_precompute(A, ui, Aeq, false, precomp);        //precompute_with_constraints reduced_solver with complementarity constraints for full sim

     //builds reduced matrices we will need for our reduced simulation. assumes we already have the A matrix pre-built.
    build_reduced_matrices();

    //construct CSM matrix, taking clusters into account. 
    build_CSM();

    //reduced_solver parameters
    max_iters = 10;
    tol = 1e-4;

    mode = 0;
    step = 0;
    //initialize simulation variables:


    timing_range = 120;
    timings = Eigen::VectorXd::Zero(timing_range);
   
}

void LBSHook::init_simulation_state()
{
    X = V0;
    x = Eigen::Map<Eigen::VectorXd>(X.data(), X.rows() * X.cols());
    V = V0;
    U = MatrixXd::Zero(V.rows(), V.cols()); Ur = U; Uc = U;
    u = Eigen::Map<Eigen::VectorXd>(U.data(), U.rows() * U.cols());
    u_prev = u; u_curr = u;
    uc = u; ur = u;
}
void LBSHook::fillModes(const Eigen::SparseMatrix<double>& H, const Eigen::SparseMatrix<double> M, Rig& rig, int num_modes)
{
    //append tiny diagonal entries to sparse matrix...
    Eigen::VectorXd m = M.diagonal();
    Eigen::VectorXd z = 1e-12 * Eigen::VectorXd::Ones(rig.get_J().cols());
    Eigen::VectorXd diag = igl::cat(1, m, z);
    Eigen::SparseMatrix<double> M_exp = sparse_diag(diag);

    //get momentum leaking matrix
    Eigen::SparseMatrix<double> D;
    momentum_leaking_matrix(X, T, fast_cd::MOMENTUM_LEAK_DIFFUSION, D);
    D = igl::repdiag(D, 3);


    //augment hessian with complementary constraint jacobians.
    Eigen::SparseMatrix<double> H_exp;
    Eigen::SparseMatrix<double> A_eq = (rig.get_J().transpose() * (M * D));
    augment_with_linear_constraints(H, A_eq, H_exp);

    Eigen::MatrixXd S_full_mat;
    bool found_modes = igl::readDMAT(B_file_dir, B_full);
    bool found_evals = igl::readDMAT(S_file_dir, S_full_mat);
    if (found_modes && found_evals)
    {
        std::cout << "Found " + std::to_string(B_full.cols()) + "/" + std::to_string(num_modes) + " Cached Modes." << std::endl;
        S_full = Eigen::Map<Eigen::VectorXd>(S_full_mat.data(), S_full_mat.rows());
    }
    bool enough_modes = B_full.cols() >= num_modes;
    if (found_modes && found_evals && enough_modes)
    {
        std::cout << "Using  " + std::to_string(num_modes) + " Cached Modes." << std::endl;
    }
    else
    {
        std::cout << "Computing " + std::to_string(num_modes) + " Modes." << std::endl;
        compute_modes(H_exp, M_exp, num_modes, B_full, S_full);
        if (B_full.cols() == num_modes)
        {
            std::cout << "Computation of Modes Successful, saving to:" << std::endl;
            std::cout << B_file_dir << std::endl;
            //  write_binary(B_file_path, B_full);
            igl::writeDMAT(B_file_dir, B_full);
            igl::writeDMAT(S_file_dir, S_full);
        }

    }
    B = B_full.block(0, 0, B_full.rows() - rig.get_p().rows(), num_modes);
    S = S_full.topRows(num_modes);
}

void LBSHook::makeClusters()
{
    Eigen::MatrixXi labels_mat;
    bool found_clusters = igl::readDMAT(labels_file_dir, labels_mat);
    if (!found_clusters)
    {
        //100 modes is usually enough for a fine clustering
        if (B_full.rows() < 100)
            fillModes(H, M, rig, 100);
        compute_grouping_matrix(100, l);
        Eigen::MatrixXi labels_mat = Eigen::Map<Eigen::MatrixXi>(labels.data(), labels.rows(), 1);
        igl::writeDMAT(labels_file_dir, labels_mat);
    }
    else
    {
        std::cout << "Found " << std::to_string(l) << " Cached Clusters" << std::endl;
        labels = Eigen::Map<Eigen::VectorXi>(labels_mat.data(), labels_mat.rows());
        //build up grouping matrix dimesnion by kron Id
        Eigen::SparseMatrix<double> G_small, G_large, S_cols, S_rows;
        grouping_matrix_from_clusters(labels, G_small);

        //I wish there was a kron function
        G_large = igl::repdiag(G_small, 3);
        interweaving_matrix(G_small.cols(), 3, S_cols);
        interweaving_matrix(G_small.rows(), 3, S_rows);
        G = S_rows.transpose() * G_large * S_cols;
    }

    //now to premultiply CSM and CSMB
    CSM = CSM * G.transpose();
}

void LBSHook::compute_grouping_matrix(int num_modes, int num_clusters)
{
    Eigen::MatrixXd B, B_faces;
    B = B_full.block(0, 0, B_full.rows() - rig.get_p().rows(), num_modes);
    Eigen::MatrixXd S;  //S contains eigenvals
    S.resize(1, num_modes);
    //   S.block(0, 0, 1, num_modes) = S_full.topRows(num_modes).transpose();
     //  B.array().rowwise() *= S.row(0).array();
     //  B.array().rowwise() /= S.row(0).array();
    igl::average_onto_faces(T, B, B_faces);

    Eigen::MatrixXd C; Eigen::VectorXi I;

    //clustering!
    std::cout << "Clustering " << std::to_string(num_clusters) << " Clusters" << std::endl;
    igl::kmeans(B_faces, num_clusters, C, I);
    std::cout << "Clustering done!";
    labels = I;


    Eigen::SparseMatrix<double> G_small, G_large, S_cols, S_rows;
    grouping_matrix_from_clusters(labels, G_small);

    //I wish there was a kron function
    G_large = igl::repdiag(G_small, 3);
    interweaving_matrix(G_small.cols(), 3, S_cols);
    interweaving_matrix(G_small.rows(), 3, S_rows);
    G = S_rows.transpose() * G_large * S_cols;

}


void LBSHook::local_global_solver_reduced(Eigen::VectorXd& z, Eigen::VectorXd& z_next)
{

    Eigen::VectorXd z0 = z;
    double diff = 1;
    Eigen::VectorXd r = ur + x;
    CSMr = CSM.transpose() * r;

    //precompute_with_constraints some linear terms, multiplications so we don't have to do it each iteration
    linear_term_reduced = BTA * r;        //positional linear term that pops out when we optimize for displacement
    Eigen::MatrixXd R_stack;
    if (use_inertia)
    {
        Eigen::VectorXd y = 2 * (u_curr + x) - (u_prev + x);            //u_curr and u_prev are total displacements, eg u_prev = uc_prev + ur_prev;
        linear_term_reduced -= (1.0 / (dt * dt)) * BTM * y;       //inertial linear term if it applies
    }
    for (int i = 0; i < max_iters; i++)
    {
        z0 = z;
        //ger best fit rotations with current reduced coefficients z
        local_step_reduced(z, R_stack);
        //global step will return displacement but will be split up into different cases if we are using complementarity or not
        global_step_reduced(R_stack, z);
        diff = (z - z0).squaredNorm();
        if (diff < tol)
            break;
    }
    z_next = z;
    //need to fille these out outside the function
  //  u_prev = u_curr;
   // u_curr = u;
   // U_next = Eigen::Map<Eigen::MatrixXd>(u.data(), U.rows(), U.cols());
}

void LBSHook::local_step_reduced(Eigen::VectorXd& z, Eigen::MatrixXd& R)
{
    //get covariance stack:
    Eigen::VectorXd cov_stack_flat = CSMB.transpose() * z;
    cov_stack_flat += CSMr;
    cov_stack_flat *= k;

    int num_clusters = do_clustering ? l : T.rows();
    //   Eigen::VectorXd cov_stack_flat_check = k * (CSMB.transpose() * z + CSMr);
    Eigen::MatrixXd cov_stack = Eigen::Map<Eigen::MatrixXd>(cov_stack_flat.data(), num_clusters * V.cols(), V.cols());
    R.resizeLike(cov_stack);
    //fit rotations
    Eigen::Matrix3d rot, cov;
    for (int i = 0; i < num_clusters; i++)
    {
        cov = cov_stack.block(3 * i, 0, 3, 3);
        igl::polar_svd3x3(cov, rot);
        R.block(3 * i, 0, 3, 3) = rot;
    }
}

void LBSHook::global_step_reduced(Eigen::MatrixXd& R, Eigen::VectorXd& z_next)
{
    Eigen::VectorXd R_flat = Eigen::Map<Eigen::VectorXd>(R.data(), R.rows() * R.cols());
    Eigen::VectorXd linear_term_rot = k * CSMB * R_flat;


    // Eigen::VectorXd linear_term_rot_check =  CSMB * R_flat;
     //fill linear_term_reduced in local_global_solver so we don't have to recompute it each timestep
    Eigen::VectorXd linear_term = linear_term_reduced - linear_term_rot;

    //do a solve_with_constraints, using precomp as a prefactorization. Make sure this prefactorization is updated everytime the system changes:
    //ie for removing inertia, reducing system, adding point constraints.
    z_next = precomp_reduced.solve(-linear_term);
}

void LBSHook::local_global_solver(Eigen::MatrixXd& U, Eigen::MatrixXd& U_next)
{
    u = Eigen::Map<Eigen::VectorXd>(U.data(), U.rows() * U.cols());
    Eigen::VectorXd u0 = u;
    double diff = 1;
    Eigen::MatrixXd R_stack;
    for (int i = 0; i < max_iters; i++)
    {
        u0 = u;
        //keep local step in terms of total displacement for now... it'as easier this way
        local_step(u, R_stack);
        //global step will return displacement but will be split up into different cases if we are using complementarity or not
        global_step(R_stack, u);
        diff = (u - u0).squaredNorm();
        if (diff < tol)
            break;
    }
    u_prev = u_curr;
    u_curr = u;
    U_next = Eigen::Map<Eigen::MatrixXd>(u.data(), U.rows(), U.cols());

}

void LBSHook::local_step(Eigen::VectorXd& u, Eigen::MatrixXd& R)
{
    int num_clusters = do_clustering ? l : T.rows();
    //get best fit rotations:
    Eigen::VectorXd cov_stack_flat = k * CSM.transpose() * (u + x);
    Eigen::MatrixXd cov_stack = Eigen::Map<Eigen::MatrixXd>(cov_stack_flat.data(), num_clusters * V.cols(), V.cols());
    R.resizeLike(cov_stack);
    //fit rotations
    Eigen::Matrix3d rot, cov;
    for (int i = 0; i < num_clusters; i++)
    {
        cov = cov_stack.block(3 * i, 0, 3, 3);
        igl::polar_svd3x3(cov, rot);
        R.block(3 * i, 0, 3, 3) = rot;
    }

}

void LBSHook::global_step(Eigen::MatrixXd& R, Eigen::VectorXd& u)
{
    //update constraints, precomputation if we added on in render thread.

    Eigen::VectorXd R_flat = Eigen::Map<Eigen::VectorXd>(R.data(), R.rows() * R.cols());
    Eigen::VectorXd linear_term_rot = k * CSM * R_flat;
    Eigen::VectorXd linear_term = -linear_term_rot;

    //optimize for displacement, additional linear term for positions

    linear_term += A * x;

    if (use_inertia)
    {
        Eigen::VectorXd y = 2 * (u_curr + x) - (u_prev + x);
        linear_term -= (1.0 / (dt * dt)) * M * y;
        //optimize for displacement, additional linear term for positions
    }

    if (use_cd)
    {
        linear_term += A * ur;
        igl::min_quad_with_fixed_solve<double>(precomp, linear_term, ub, Beq, uc);
        u = ur + uc;
    }
    else
    {
        //do a solve_with_constraints, using precomp as a prefactorization. Make sure this prefactorization is updated everytime the system changes:
        //ie for removing inertia, reducing system, adding point constraints.
        igl::min_quad_with_fixed_solve<double>(precomp, linear_term, ub, Beq, u);
    }

}

void LBSHook::load_file_paths()
{
    data_dir = "../data/" + mesh_name + "/";
    cache_dir = data_dir + "cache/";

    //if cache dir doesn't exist, make it dammit!
    bool cache_exists = std::filesystem::exists(cache_dir);
    if (!cache_exists)
        std::filesystem::create_directory(cache_dir);
    //eigenvectors
    B_file_dir = cache_dir + "B_lbs_rig.DMAT";
    //eigenvalues
    S_file_dir = cache_dir + "S_lbs_rig.DMAT";

    //rig information
    rig_file_dir = cache_dir + "rig_data.json";

    //clustered labels, use this to build grouping matrix G
    labels_file_dir = cache_dir + "labels_" + std::to_string(l) + "_lbs_rig.DMAT";
}

void LBSHook::compute_cd_equality_constraints()
{
    //get momentum leaking matrix... should I make this global??  don't need to yet
    Eigen::SparseMatrix<double> D;
    momentum_leaking_matrix(V0, T, fast_cd::MOMENTUM_LEAK_DIFFUSION, D);
    D = igl::repdiag(D, 3);
    Aeq = (rig.get_J().transpose() * (M * D));
    Beq = Eigen::VectorXd::Zero(Aeq.rows());
}


void LBSHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu& menu)
{

}

void LBSHook::renderRenderGeometry(igl::opengl::glfw::Viewer& viewer)
{
    rig.render(viewer);
    viewer.data_list[0].set_vertices(V);

    viewer.data_list[0].compute_normals();
    if (handleV.rows() > 0)
    {
        viewer.data_list[0].set_points(handleV, Eigen::RowVector3d(1.0, 0.2, 0.2));

    }
    else if (viewer.data().points.rows() > 0)
    {
        viewer.data_list[0].clear_points();
    }

    if (do_clustering)
    {
        Eigen::VectorXi drawing_labels;
        igl::slice(labels, FiT, drawing_labels);
        //      viewer.data().set_data(drawing_labels.cast<double>());
    }

 
}


void LBSHook::initViewer(igl::opengl::glfw::Viewer& viewer)
{
    this->viewer = &viewer;
    viewer.data_list[0].set_mesh(V, F);

    viewer.data_list[0].point_size = 10;
  
    //TODO: should keep matcap separate for each mesh, and load it up once we pick our mesh
    Eigen::Matrix<unsigned char, -1, -1> R, G, B, A;
    std::string matcap_filepath = "../data/jade.png";
    igl::png::readPNG(matcap_filepath, R, G, B, A);

    viewer.data_list[0].invert_normals = true;
    viewer.data_list[0].double_sided = false;
    viewer.data_list[0].set_face_based(true);
  
    //turn data_list[0] on/off
    viewer.data_list[0].set_texture(R, G, B, A);
    viewer.data_list[0].use_matcap = true;
    viewer.data_list[0].show_lines = true;
    viewer.data_list[0].show_faces = true;

    rig.init_viewer(viewer);
}


void LBSHook::change_stiffness()
{
    build_A_matrix();
    // handle reduced matrices that depend on A, and then precompute_with_constraints as needed

    BTA = B.transpose() * A;
    Q = BTA * B;

    //recompute both reduced and unreduced factorizations
    if (do_reduction)
    {
        precomp_reduced.compute(Q);
    }
    else
    {
        igl::min_quad_with_fixed_precompute(A, ui, Aeq, false, precomp);
    }
    
}

void LBSHook::build_A_matrix()
{
    if (use_inertia)
    {
        A = k * C + (1.0 / (dt * dt)) * M;
    }
    else
    {
        Eigen::SparseMatrix<double> speye = Eigen::SparseMatrix<double>(C.rows(), C.cols());
        speye.setIdentity(); // add  tiny tikhnov reg, ill posed when there are no constraints
        A = k * C + 1e-12 * speye;
    }
}

void LBSHook::build_CSM()
{
    covariance_scatter_matrix(X, T, CSM);
    if (do_clustering)
    {
        makeClusters();
    }
    CSM = igl::repdiag(CSM, 3);

    if (do_reduction)
    {
        assert((B.cols() == r) && "Eigenmodes not properly filled out before computing CSMB");
        //get a reduced covariance scatter matrix ready to go. if you want to make use of this, needs to have B matrices already computed
        CSMB = B.transpose() * CSM;
    }

}


void LBSHook::build_reduced_matrices()
{
    if (do_reduction)
    {
        //fill modes... writes out the modes internally and stores the eigenvectors in B, and values in S. First tries to read themf rom a cache
        fillModes(H, M, rig, r);
        //quadratic coefficients
        Q = B.transpose() * A * B;

        //need this matrix because a linear term pops out if we optimize for complementary displacements (instead of positions)
        BTA = B.transpose() * A;

        //need this matrix for inertial term
        BTM = B.transpose() * M;


        //precompute_with_constraints system factorization to make solve_with_constraints extra fast
        precomp_reduced.compute(Q);

        //get vector of reduced coefficients ready to go
        z = Eigen::VectorXd::Zero(r);
    }
}

void LBSHook::poll_sim_changes()
{
    if (do_clustering != new_do_clustering)
    {
        do_clustering = new_do_clustering;
        if (do_clustering)
            std::cout << "Turned on clustering, clustering with " + std::to_string(l) + " clusters" << std::endl;
        else
            std::cout << "Turned off clustering" << std::endl;
        build_CSM();
    }
    if (new_stiffness != k)
    {
        k = new_stiffness;
        std::cout << "Changed stiffness to " + std::to_string(k) + ", recomputing matrices.." << std::endl;

        change_stiffness();
    }

    if (new_num_clusters != l)
    {
        l = new_num_clusters;
        std::cout << "Changed #Clusters to " + std::to_string(l) + " checking cache..." << std::endl;
        load_file_paths();          //need to reload these if we're going to try to find a cached file
        build_CSM();
    }

    if (new_num_modes != r)
    {
        r = new_num_modes;
        std::cout << "Changed #Modes to " + std::to_string(r) + " checking cache..." << std::endl;
        load_file_paths();              //need to reload these to find appropriate caches
        //builds reduced matrices we will need for our reduced simulation. assumes we already have the A matrix pre-built.
        build_reduced_matrices();
        //rebuild CSM, though really we only need one line
        build_CSM();
    }




    //if (changed_constraint) set_constraints(); not handling point constraints for CPP code y
}



bool LBSHook::mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier) {

    return rig.mouse_down(viewer,  button, modifier);
    /*

    if (button == 0) //left click means we are making a new vertex
    {
        Eigen::Vector2f mouse_win_2d = Eigen::Vector2f(viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y);
        mouse_win = Eigen::Vector3d(mouse_win_2d.x(), mouse_win_2d.y(), 0);
        int fid;
        Eigen::Vector3d bc;
        bool hit = false;
        //bool hit = igl::unproject_onto_mesh(mouse_win_2d, viewer.core().view, viewer.core().proj, viewer.core().viewport, hook->V, hook->F, fid, bc);
        //if you click on the mesh select the vertex, otherwise do nothing
        if (hit)
        {
            std::cout << "we have a hit at face id " + std::to_string(fid) << std::endl;
            igl::unproject(
                mouse_win,
                viewer.core().view,
                viewer.core().proj,
                viewer.core().viewport,
                mouse_world);
            //  hook->addConstraintHandle(fid, bc);
            return true;
        }
    
    }

    if (button == 1)  //middle click means nothing
    {
    }

    if (button == 2) //right click means we are moving a vertex
    {   //copying geometry processing assignment, pba assignment ui doesn't properly move vertices
        if (handleI.rows() > 0)
        {
            mouse_win = Eigen::Vector3d(viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, 0);
            Eigen::MatrixXd CP; //closest point
            igl::project(
                handleV,
                viewer.core().view,
                viewer.core().proj, viewer.core().viewport, CP);
            Eigen::VectorXd D = (CP.rowwise() - mouse_win.transpose()).rowwise().norm();
            currentHandle = (D.minCoeff(&currentHandle) < 30) ? currentHandle : -1;
            if (currentHandle != -1)
            {
                mouse_dragging = true;
                mouse_win(2) = CP(currentHandle, 2);
            }
            //find closest vertex
            //set selected index to it
            //wait for move callback toa ctually move the vertex
        }

    }
        */

    return false;

}

bool LBSHook::mouse_move(igl::opengl::glfw::Viewer& viewer, int x, int y) {

    return rig.mouse_move(viewer, x, y);
    /*
    if (mouse_dragging)
    {
        mouse_drag_win = Eigen::Vector3d(viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, mouse_win(2));
        //   mouse_win = Eigen::Vector3d(viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, 0.);

        igl::unproject(
            mouse_drag_win,
            viewer.core().view,
            viewer.core().proj,
            viewer.core().viewport,
            mouse_drag_world);

        //   hook->moveConstraintHandle(mouse_drag_world);
        return true;
    }
    */

   
}

bool LBSHook::mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
    return rig.mouse_up(viewer, button, modifier);
   // mouse_dragging = false;
   // return false;
}
