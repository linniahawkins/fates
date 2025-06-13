cmake_minimum_required(VERSION 3.12)
project(fatesMLBtran LANGUAGES Fortran)
enable_language(Fortran)
set(btran_sources MLBtranTest.F90)

set(USE_FTORCH ON CACHE BOOL "Enable internal FTorch build")
add_subdirectory(${SRC_ROOT}/libraries/FTorch ${CMAKE_BINARY_DIR}/external/ftorch)
add_executable(ml_btran_exe ${phen_sources})
target_link_directories(ml_btran_exe PRIVATE ${PFUNIT_TOP_DIR}/lib)
target_link_libraries(ml_btran_exe PRIVATE
                  fates
                  ftorch_wrapper
                  ftorch
                  csm_share
                  ${Torch_LIB_PATH}
                  funit)


set_target_properties(ml_btran_exe PROPERTIES
    LINKER_LANGUAGE Fortran
)
