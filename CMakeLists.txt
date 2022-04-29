cmake_minimum_required(VERSION 3.1)

file(GLOB SRC_FILES ./src/*)
file(GLOB HEADER_FILES ./include/*)


option(FAST_CD_USE_MATLAB "use matlab" ON)
option(FAST_CD_USE_SUITESPARSE "use suitesparse" ON)
option(LIBIGL_PNG  "Build target igl::png"  ON)
option(LIBIGL_GLFW "Build target igl::glfw" ON)
option(LIBIGL_IMGUI   "Build target igl::imgui" ON)
option(LIBIGL_RESTRICTED_MATLAB   "Build target igl_restricted::matlab" ${FAST_CD_USE_MATLAB})
# libigl
add_subdirectory(./deps/libigl/)

#add spectra
add_library(spectra INTERFACE)
target_include_directories(spectra INTERFACE ./deps/spectra/include/)
target_link_libraries(spectra INTERFACE Eigen3::Eigen)


add_library(fast_cd ${SRC_FILES} ${HEADER_FILES})
set(FAST_CD_INCLUDE_DIRS ./include ./deps/json)
target_include_directories(fast_cd PUBLIC ${FAST_CD_INCLUDE_DIRS})
set_target_properties(fast_cd PROPERTIES LINKER_LANGUAGE CXX)
set_property(TARGET fast_cd PROPERTY CXX_STANDARD 17)

#find suitesparse (for windows, follow https://github.com/jlblancoc/suitesparse-metis-for-windows)
option(FAST_CD_USE_SUITESPARSE "FAST_CD_USE_SUITESPARSE" OFF)
if(FAST_CD_USE_SUITESPARSE)
  target_compile_definitions(fast_cd PUBLIC -DFAST_CD_USE_SUITESPARSE)
endif()

option(LIBIGL_RESTRICTED_MATLAB "USE_MATLAB" ON)
if(FAST_CD_USE_MATLAB)
  target_compile_definitions(fast_cd PUBLIC -DFAST_CD_USE_MATLAB)
endif()


if (FAST_CD_USE_SUITESPARSE)
    set(WITH_FORTRAN OFF CACHE BOOL "")
    set(WITH_CUDA    OFF CACHE BOOL "")
    set(WITH_OPENMP  OFF CACHE BOOL "")
    set(WITH_PARTITION OFF CACHE BOOL "")
    set(WITH_DEMOS OFF CACHE BOOL "")
    set(BUILD_CXSPARSE OFF CACHE BOOL "")
    add_subdirectory("SuiteSparse/")
	if (FAST_CD_USE_MATLAB)
          set(FAST_CD_LIBRARIES spectra umfpack cholmod igl::core igl::glfw igl::imgui igl::png igl_restricted::matlab)
	else()
          set(FAST_CD_LIBRARIES spectra umfpack cholmod igl::core igl::glfw igl::imgui igl::png )
	endif()
	
else()
  set(FAST_CD_LIBRARIES spectra igl::core igl::glfw igl::imgui igl::png igl_restricted::matlab)
endif()
target_link_libraries(fast_cd ${FAST_CD_LIBRARIES})

add_library(fast_cd::fast_cd INTERFACE IMPORTED GLOBAL)
target_include_directories(fast_cd::fast_cd INTERFACE ${FAST_CD_INCLUDE_DIRS})
target_link_libraries(fast_cd::fast_cd INTERFACE ${FAST_CD_LIBRARIES} fast_cd)
add_dependencies(fast_cd::fast_cd fast_cd)
