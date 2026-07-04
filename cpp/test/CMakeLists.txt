add_executable(run_tests test_static_prefix_sum.cpp)
target_link_libraries(run_tests PRIVATE Catch2::Catch2WithMain)
add_test(NAME run_tests COMMAND run_tests)
