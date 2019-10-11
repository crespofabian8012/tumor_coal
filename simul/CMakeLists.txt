# code relating to tree structured stick breaking process
include_directories(/usr/local/include/ ${CMAKE_SOURCE_DIR}/inference/)

file( GLOB SRCS *.cpp *.hpp *.c *.h)
add_executable(simul ${SRCS}) # create the executable file

find_package (PkgConfig REQUIRED)
find_package(Boost 1.67 COMPONENTS regex system filesystem unit_test_framework program_options REQUIRED)
find_package (GSL REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

target_link_libraries(simul GSL::gsl GSL::gslcblas /usr/local/lib/libboost_system.a /usr/local/lib/libboost_filesystem.a /usr/local/lib/libboost_program_options.a /usr/local/lib/libpll.a /usr/local/lib/libpll_algorithm.a /usr/local/lib/libpll_tree.a /usr/local/lib/libpll_util.a /usr/local/lib/libpll_optimize.a /usr/local/lib/libpll_msa.a	/usr/local/lib/libpll_binary.a libTumorCoal)

INSTALL(TARGETS simul DESTINATION ${CMAKE_SOURCE_DIR}/)
