set(phen_sources MLPhenologyTest.F90)

link_directories(${PFUNIT_TOP_DIR}/lib)

add_executable(ml_phen_exe ${phen_sources})

target_link_libraries(ml_phen_exe
                  fates
                  csm_share
                  funit)