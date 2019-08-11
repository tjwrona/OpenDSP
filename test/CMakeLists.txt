include(GoogleTest)

add_executable(opendsp_test
    main.cpp
    dft_test.cpp
)

target_compile_options(opendsp_test
  PRIVATE -D_USE_MATH_DEFINES
)

target_link_libraries(opendsp_test
    opendsp::opendsp
    CONAN_PKG::gtest
)

gtest_discover_tests(opendsp_test)