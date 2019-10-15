
include_directories(/usr/local/include/)

file( GLOB SRCS *.cpp *.hpp *.c *.h)
add_executable(run ${SRCS}) # create the executable file
add_library(libTumorCoal STATIC ${SRCS})

find_package (PkgConfig REQUIRED)
find_package(Boost 1.67 COMPONENTS regex system filesystem unit_test_framework program_options REQUIRED)
find_package (GSL REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

target_link_libraries(run GSL::gsl GSL::gslcblas /usr/local/lib/libboost_system.a /usr/local/lib/libboost_filesystem.a /usr/local/lib/libboost_program_options.a /usr/local/lib/libpll.a /usr/local/lib/libpll_algorithm.a /usr/local/lib/libpll_tree.a /usr/local/lib/libpll_util.a /usr/local/lib/libpll_optimize.a /usr/local/lib/libpll_msa.a    /usr/local/lib/libpll_binary.a /usr/local/lib/libpll_tree.0.dylib /usr/local/lib/libpll_binary.a /usr/local/lib/libpll_tree.dylib /usr/local/lib/libpll_binary.0.dylib /usr/local/lib/libpll_binary.dylib /usr/local/lib/libpll_algorithm.0.dylib /usr/local/lib/libpll_algorithm.dylib /usr/local/lib/libpll_msa.0.dylib /usr/local/lib/libpll_msa.dylib /usr/local/lib/libpll_optimize.0.dylib /usr/local/lib/libpll_util.0.dylib /usr/local/lib/libpll.0.dylib /usr/local/lib/libpll.dylib)

target_link_libraries(libTumorCoal GSL::gsl GSL::gslcblas /usr/local/lib/libboost_system.a /usr/local/lib/libboost_filesystem.a /usr/local/lib/libboost_program_options.a /usr/local/lib/libpll.a /usr/local/lib/libpll_algorithm.a /usr/local/lib/libpll_tree.a /usr/local/lib/libpll_util.a /usr/local/lib/libpll_optimize.a /usr/local/lib/libpll_msa.a    /usr/local/lib/libpll_binary.a /usr/local/lib/libpll_tree.0.dylib /usr/local/lib/libpll_binary.a /usr/local/lib/libpll_tree.dylib /usr/local/lib/libpll_binary.0.dylib /usr/local/lib/libpll_binary.dylib /usr/local/lib/libpll_algorithm.0.dylib /usr/local/lib/libpll_algorithm.dylib /usr/local/lib/libpll_msa.0.dylib /usr/local/lib/libpll_msa.dylib  /usr/local/lib/libpll_optimize.0.dylib  /usr/local/lib/libpll_util.0.dylib /usr/local/lib/libpll.0.dylib /usr/local/lib/libpll.dylib)

INSTALL(TARGETS run DESTINATION ${CMAKE_SOURCE_DIR}/)
install(TARGETS libTumorCoal DESTINATION bin)

