add_library(opendsp
    analysis.cpp
    dft.cpp
)

add_library(opendsp::opendsp
  ALIAS opendsp
)

target_compile_options(opendsp
  PRIVATE -D_USE_MATH_DEFINES
)

target_include_directories(opendsp
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:include>
  PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/opendsp>
          $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
)
