set(math_sources FatesTestMathUtils.F90)

link_directories(${PFUNIT_TOP_DIR}/lib)

add_executable(ml_phen_exe ${math_sources})

target_link_libraries(FATES_math_exe
                  fates
                  csm_share
                  funit)