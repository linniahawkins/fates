cmake_minimum_required(VERSION 3.12)
project(fatesMLexample LANGUAGES Fortran)
enable_language(Fortran)
set(example_sources exampleTest.F90)

#set(USE_FTORCH ON CACHE BOOL "Enable internal FTorch build")
#add_subdirectory(${SRC_ROOT}/libraries/FTorch ${CMAKE_BINARY_DIR}/external/ftorch)
add_executable(ml_example_exe ${example_sources})
target_link_directories(ml_example_exe PRIVATE ${PFUNIT_TOP_DIR}/lib)
target_link_libraries(ml_example_exe PRIVATE
                  fates
                  ftorch_wrapper
                  ftorch
                  csm_share
                  ${Torch_LIB_PATH}
                  funit)


set_target_properties(ml_example_exe PROPERTIES
    LINKER_LANGUAGE Fortran
)
