cmake_minimum_required(VERSION 3.1)

file(GLOB SRC_FILES ./src/*)
file(GLOB HEADER_FILES ./include/*)


option(LIBIGL_PNG  "Build target igl::png"  ON)
option(LIBIGL_GLFW "Build target igl::glfw" ON)
option(LIBIGL_IMGUI   "Build target igl::imgui" ON)
option(LIBIGL_RESTRICTED_MATLAB   "Build target igl_restricted::matlab" ON)
# libigl
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")



add_library(fast_cd ${SRC_FILES} ${HEADER_FILES})
set(FAST_CD_INCLUDE_DIRS ./include/ ./deps/json/)
target_include_directories(fast_cd PUBLIC ${FAST_CD_INCLUDE_DIRS})
set_target_properties(fast_cd PROPERTIES LINKER_LANGUAGE CXX)
set_property(TARGET fast_cd PROPERTY CXX_STANDARD 17)



set(FAST_CD_LIBRARIES igl::core igl::glfw igl::imgui igl::png igl_restricted::matlab)

target_link_libraries(fast_cd ${FAST_CD_LIBRARIES})

add_library(fast_cd::fast_cd INTERFACE IMPORTED GLOBAL)
target_include_directories(fast_cd::fast_cd INTERFACE ${FAST_CD_INCLUDE_DIRS})
target_link_libraries(fast_cd::fast_cd INTERFACE ${FAST_CD_LIBRARIES} fast_cd)
add_dependencies(fast_cd::fast_cd fast_cd)
