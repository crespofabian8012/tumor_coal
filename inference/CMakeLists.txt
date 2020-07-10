set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lgsl -fopenmp  -std=c++17  -stdlib=libc++")
include_directories(/usr/local/opt/llvm/include ${CMAKE_SOURCE_DIR}/install/include/ /usr/local/include/
/usr/local/include/gsl )

message(${CMAKE_SOURCE_DIR})


file( GLOB SRCS *.cpp *.hpp *.c *.h)
add_executable(run ${SRCS}) # create the executable file
add_library(libTumorCoal STATIC ${SRCS})



find_package (PkgConfig REQUIRED)
find_package (GSL REQUIRED)
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


find_package(Boost 1.67 COMPONENTS regex system filesystem unit_test_framework program_options REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()



target_link_libraries(run ${CMAKE_SOURCE_DIR}/install/lib/libpll_algorithm.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_optimize.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_tree.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_util.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_msa.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_binary.a ${CMAKE_SOURCE_DIR}/install/lib/libpll.a ${Boost_LIBRARIES}  GSL::gsl OpenMP::OpenMP_CXX)

target_link_libraries(libTumorCoal
    ${Boost_LIBRARIES}  ${CMAKE_SOURCE_DIR}/install/lib/libpll_algorithm.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_optimize.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_tree.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_util.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_msa.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_binary.a ${CMAKE_SOURCE_DIR}/install/lib/libpll.a  GSL::gsl OpenMP::OpenMP_CXX)
message(${CMAKE_SOURCE_DIR})
INSTALL(TARGETS run DESTINATION ${CMAKE_SOURCE_DIR}/)
INSTALL(TARGETS libTumorCoal DESTINATION bin)

include_directories(${Boost_INCLUDE_DIRS})


