#  The tests defined here will only be run if "benchmarks" configuration is specified to ctest,
#  for example: "ctest -C benchmarks" will run all tests including the benchmark tests, and
#               "ctest -C benchmarks -L benchmark" will run only the benchmark tests, while
#               "ctest" will run all tests except the "benchmark" configuration tests.
#               ("ctest -L benchmark" will run no tests.)

macro(dirac_benchmark _name _label)
    add_test(NAME ${_name}
             COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/test/${_name}/test --binary-dir=${PROJECT_BINARY_DIR} --work-dir=${PROJECT_BINARY_DIR}/test/${_name} --verbose
             CONFIGURATIONS benchmarks)
    if(NOT "${_label}" STREQUAL "")
        set_tests_properties(${_name} PROPERTIES LABELS "${_label}")
    endif()
    # do not include benchmark tests in default ctest test suite
endmacro()

dirac_benchmark(benchmark_cc "benchmark;cc" "10000")
dirac_benchmark(benchmark_cc_linear "benchmark;cc" "3600")
dirac_benchmark(benchmark_gaunt "benchmark;gaunt" "3600")
dirac_benchmark(benchmark_mcscf_energy "benchmark;mcscf" "3600")
dirac_benchmark(benchmark_moltra_4or6 "benchmark;moltra" "3600")
dirac_benchmark(benchmark_eomip_energy_symmetry "benchmark;cc" "")
dirac_benchmark(benchmark_eomea_energy_symmetry "benchmark;cc" "")
dirac_benchmark(benchmark_eomee_energy_symmetry "benchmark;cc" "")
