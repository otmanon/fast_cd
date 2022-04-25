cmake_minimum_required(VERSION 3.1)

file(GLOB SRC_FILES ./src/*)
file(GLOB HEADER_FILES ./include/*)


option(LIBIGL_PNG  "Build target igl::png"  ON)
option(LIBIGL_GLFW "Build target igl::glfw" ON)
option(LIBIGL_IMGUI   "Build target igl::imgui" ON)
# libigl
add_subdirectory(./deps/libigl/)

#add spectra
add_library(spectra INTERFACE)
target_include_directories(spectra INTERFACE ./deps/spectra/include/)
target_link_libraries(spectra INTERFACE Eigen3::Eigen)


add_library(fast_cd ${SRC_FILES} ${HEADER_FILES})
target_include_directories(fast_cd PRIVATE ./include/)
target_include_directories(fast_cd PUBLIC ./deps/json/)
set_target_properties(fast_cd PROPERTIES LINKER_LANGUAGE CXX)
set_property(TARGET fast_cd PROPERTY CXX_STANDARD 17)

#find suitesparse (for windows, follow https://github.com/jlblancoc/suitesparse-metis-for-windows)
option(FAST_CD_USE_SUITESPARSE "FAST_CD_USE_SUITESPARSE" ON)
add_definitions( -DFAST_CD_USE_SUITESPARSE=${FAST_CD_USE_SUITESPARSE} )

option(LIBIGL_RESTRICTED_MATLAB "USE_MATLAB" ON)
add_definitions( -DFAST_CD_USE_MATLAB=${FAST_CD_USE_MATLAB} )


if (FAST_CD_USE_SUITESPARSE)
	find_package(SuiteSparse CONFIG REQUIRED)
	if (FAST_CD_USE_MATLAB)
	target_link_libraries(fast_cd spectra SuiteSparse::umfpack SuiteSparse::cholmod igl::core igl::glfw igl::imgui igl::png igl_restricted::matlab)
	else()
	target_link_libraries(fast_cd spectra SuiteSparse::umfpack SuiteSparse::cholmod igl::core igl::glfw igl::imgui igl::png igl_restricted::matlab)
	endif()
	
else()
	target_link_libraries(fast_cd spectra igl::core igl::glfw igl::imgui igl::png igl_restricted::matlab)
endif()

