if (Pism_BUILD_EXTRA_EXECS)

  foreach (TEST IN ITEMS
      ssa_testi ssa_testj ssa_test_const ssa_test_linear ssa_test_plug ssa_test_cfbc)
    add_executable (${TEST} ${TEST}.cc)
    target_link_libraries (${TEST} pism)
    install (TARGETS ${TEST} DESTINATION ${Pism_BIN_DIR})
  endforeach()

endif ()
