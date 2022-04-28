#ifdef IGL_STATIC_LIBRARY
#include <igl/repdiag.cpp>
template void igl::repdiag<double>(Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, int, Eigen::Matrix<double, -1, -1, 0, -1, -1>&);
#include <igl/average_onto_faces.cpp>
template void igl::average_onto_faces<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#include <igl/slice_into.cpp>
template void igl::slice_into<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase< Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::slice_into< Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, 1, 3, 1, 1, 3> >(Eigen::MatrixBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase< Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase< Eigen::Matrix<int, 1, 3, 1, 1, 3> > const&, Eigen::PlainObjectBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#include <igl/writeDMAT.cpp>
template bool igl::writeDMAT<Eigen::Matrix<int, -1, 1, 0, -1, 1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> >, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, bool);

#include <igl/unproject_in_mesh.cpp>
template int igl::unproject_in_mesh<class Eigen::Matrix<double, -1, -1, 0, -1, -1>, class Eigen::Matrix<int, -1, -1, 0, -1, -1>, class Eigen::Matrix<double, -1, -1, 0, -1, -1> >(class Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, class Eigen::Matrix<float, 4, 4, 0, 4, 4> const&, class Eigen::Matrix<float, 4, 4, 0, 4, 4> const&, class Eigen::Matrix<float, 4, 1, 0, 4, 1> const&, class Eigen::MatrixBase<class Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, class Eigen::MatrixBase<class Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, class Eigen::PlainObjectBase<class Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, class std::vector<struct igl::Hit, class std::allocator<struct igl::Hit> >&);
template int igl::unproject_in_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, 3, 1, 0, 3, 1> >(Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 4, 4, 0, 4, 4> const&, Eigen::Matrix<float, 4, 4, 0, 4, 4> const&, Eigen::Matrix<float, 4, 1, 0, 4, 1> const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> >&, std::vector<igl::Hit, std::allocator<igl::Hit> >&);

#include <igl/procrustes.cpp>
template void  igl::procrustes< Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, double, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<class Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, bool, bool, double&, Eigen::PlainObjectBase< Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<class Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
#include <igl/polar_dec.cpp>
template void igl::polar_dec<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);

#endif
