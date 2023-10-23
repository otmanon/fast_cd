cmake_minimum_required(VERSION 3.1)

file(GLOB SRC_FILES ./src/*)
file(GLOB HEADER_FILES ./include/*)

#add_subdirectory(./deps/OpenBLAS/)

add_compile_definitions(IGL_VIEWER_VIEWER_QUIET)
option(LIBIGL_PNG  "Build target igl::png"  ON)
option(LIBIGL_GLFW "Build target igl::glfw" ON)
option(LIBIGL_IMGUI   "Build target igl::imgui" ON)
option(LIBIGL_RESTRICTED_MATLAB   "Build target igl_restricted::matlab" OFF)

add_subdirectory(./deps/libigl/)



add_library(fast_cd ${SRC_FILES} ${HEADER_FILES})
set(FAST_CD_INCLUDE_DIRS ./include/ ./deps/json/)
target_include_directories(fast_cd PUBLIC ${FAST_CD_INCLUDE_DIRS})
set_target_properties(fast_cd PROPERTIES LINKER_LANGUAGE CXX)
set_property(TARGET fast_cd PROPERTY CXX_STANDARD 17)

#find_package(SuiteSparse)

set(FAST_CD_LIBRARIES igl::core igl::glfw igl::imgui igl::png)# SuiteSparse::umfpack)

target_link_libraries(fast_cd ${FAST_CD_LIBRARIES})

add_library(fast_cd::fast_cd INTERFACE IMPORTED GLOBAL)
target_include_directories(fast_cd::fast_cd INTERFACE ${FAST_CD_INCLUDE_DIRS})
target_link_libraries(fast_cd::fast_cd INTERFACE ${FAST_CD_LIBRARIES} fast_cd)
add_dependencies(fast_cd::fast_cd fast_cd)