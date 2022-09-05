set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++17 -stdlib=libc++ -D_REENTRANT -Denable-glibcxx-debug=ON -DBOOST_DISABLE_ASSERTS")

set(HEADER_FILES ${CMAKE_SOURCE_DIR}/install/include/stan/lib/stan_math/lib/stan⁩/lib⁩/stan_math⁩/stan/math/prim/fun/Eigen.hpp )
include_directories(${CMAKE_SOURCE_DIR}/install/include/stan/lib/stan_math/lib/stan⁩/lib⁩/stan_math⁩/stan/math/prim/fun/Eigen.hpp)

include_directories(/usr/local/opt/llvm/include ${CMAKE_SOURCE_DIR}/install/include/
/usr/local/include/
/usr/local/include/gsl
${CMAKE_SOURCE_DIR}/install/include/stan/src/stan/
${CMAKE_SOURCE_DIR}/install/include/stan/lib/stan_math/
${CMAKE_SOURCE_DIR}/install/include/stan/lib/stan_math/lib/sundials_5.7.0/include/
${CMAKE_SOURCE_DIR}/install/include/stan/lib/stan_math/lib/tbb/tbb_2020.3/include/
${CMAKE_SOURCE_DIR}/install/include/stan/lib/stan_math/lib/eigen_3.3.9/include/
${CMAKE_SOURCE_DIR}/install/include/stan/lib/stan_math/lib/stan⁩/lib⁩/stan_math⁩/stan/math/prim/fun/
/usr/local/include/tbb/
 )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake/FINDTBB.cmake")

include(FeatureSummary)
set(CMAKE_VERBOSE_MAKEFILE ON)

link_directories(/usr/local/lib)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}")


file( GLOB SRCS *.cpp *.hpp *.c *.h)
add_executable(run ${SRCS}) # create the executable file
add_library(libTumorCoal STATIC ${SRCS}  )

set ( PROJECT_LINK_LIBS libSPF.a )
link_directories(/usr/local/bin)
include_directories(/usr/local/include)
include_directories(/usr/local/include/spf)

file(GLOB SRC *.cpp *.hpp)
find_package(GSL REQUIRED)


find_package (PkgConfig REQUIRED)
find_package (GSL REQUIRED)
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#find_package (Eigen3 3.3.9  REQUIRED NO_MODULE)
find_package (Eigen3 3.3.9 REQUIRED NO_MODULE)

find_path(EIGEN3_INCLUDE_DIR NAMES signature_of_eigen3_matrix_library
    PATHS
    ${CMAKE_INSTALL_PREFIX}/include
    ${KDE4_INCLUDE_DIR}
    PATH_SUFFIXES eigen3 eigen
)

if (TARGET Eigen3::Eigen)
  message( "--Found EIGEN:" )
  message(${EIGEN3_INCLUDE_DIR})
endif (TARGET Eigen3::Eigen)

find_package(GTest  QUIET)
set_package_properties(GTest PROPERTIES
        DESCRIPTION "C++ unit test framework"
        URL "https://github.com/google/googletest"
        TYPE OPTIONAL)
        
include_directories(${GTEST_INCLUDE_DIRS})

find_package(Sundials 5.8.0 QUIET)
set_package_properties(Sundials PROPERTIES
        DESCRIPTION "SUite of Nonlinear and DIfferential/ALgebraic equation Solvers."
        URL "https://github.com/LLNL/sundials"
        TYPE REQUIRED)

find_package(OpenCL  QUIET)
set_package_properties(OpenCL PROPERTIES
        DESCRIPTION "Framework for writing programs that execute across heterogeneous platforms"
        PURPOSE "Necessary for running OpenCL-enabled unit tests."
        TYPE OPTIONAL)

find_package(MPI QUIET)
set_package_properties(MPI PROPERTIES
        DESCRIPTION "Message passing interface implementation"
        PURPOSE "Neccessary in order to run MPI-enabled unit tests."
        TYPE RUNTIME)

find_package(Doxygen QUIET)
set_package_properties(Doxygen PROPERTIES
        DESCRIPTION "Tool for generating code documentation"
        URL "https://github.com/doxygen/doxygen"
        TYPE OPTIONAL)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads QUIET)
set_package_properties(Threads PROPERTIES
        DESCRIPTION "POSIX-compatible system threading library."
        PURPOSE "Required to run multithreaded Stan tests."
        TYPE OPTIONAL)
    
OPTION(USE_STAN "Use Stan Math" ON)
if(USE_STAN)
    message( "--Found STAN" )
    add_definitions(-pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS -DUSE_STAN )
    #set(STAN_LIBRARIES ${CMAKE_CURRENT_LIST_DIR}/install/include/stan/lib/stan_math/lib/tbb/libtbb.dylib)
    include_directories(${CMAKE_CURRENT_LIST_DIR}/install/include/stan/lib/stan_math)
    include_directories(${CMAKE_CURRENT_LIST_DIR}/install/include/stan/src/stan/io)
    include_directories(${CMAKE_CURRENT_LIST_DIR}/install/include/stan/src/stan/model)
    set(Eigen3_FOUND true)
    # use Eigen distribution from Stan Math
    include_directories(${CMAKE_CURRENT_LIST_DIR}/install/include/stan/lib/stan_math/lib/boost_1.75.0)
    include_directories(${CMAKE_CURRENT_LIST_DIR}/install/include/stan/lib/stan_math/lib/eigen_3.3.9)
    include_directories(${CMAKE_CURRENT_LIST_DIR}/install/include/stan/lib/stan_math/lib/cpplint_1.4.5)
    include_directories(${CMAKE_CURRENT_LIST_DIR}/install/include/stan/lib/stan_math/lib/opencl_3.0.0)
    include_directories(${CMAKE_CURRENT_LIST_DIR}/install/include/stan/lib/stan_math/lib/sundials_5.7.0/include)
    include_directories(${CMAKE_CURRENT_LIST_DIR}/install/include/stan/lib/stan_math/lib/tbb/tbb_2020.3/include)
endif(USE_STAN)


#find_package(Boost 1.67 COMPONENTS regex system filesystem unit_test_framework program_options REQUIRED)
#find_package(Boost 1.67 COMPONENTS  system filesystem unit_test_framework program_options date_time filesystem iostreams REQUIRED)

find_package(Boost 1.67 QUIET COMPONENTS  serialization   regex system filesystem unit_test_framework program_options iostreams date_time)
set_package_properties(Boost PROPERTIES
        DESCRIPTION "Free peer-reviewed portable C++ source libraries"
        URL "https://www.boost.org/"
        TYPE REQUIRED)

if(Boost_FOUND)
    #target_include_directories(run PRIVATE ${Boost_INCLUDE_DIRS})
    include_directories( PRIVATE ${Boost_INCLUDE_DIRS} ${Boost_LIBRARIES})
    target_compile_definitions(Boost::headers INTERFACE BOOST_DISABLE_ASSERTS)
    add_feature_info("MPI Support" Boost_mpi_FOUND "This requires Boost::MPI.")
endif()

find_package(TBB REQUIRED)



feature_summary(INCLUDE_QUIET_PACKAGES WHAT ALL
FATAL_ON_MISSING_REQUIRED_PACKAGES)

target_link_libraries(run  ${CMAKE_SOURCE_DIR}/install/lib/libpll_algorithm.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_optimize.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_tree.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_util.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_msa.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_binary.a ${CMAKE_SOURCE_DIR}/install/lib/libpll.a  ${Boost_LIBRARIES} Boost::date_time Boost::filesystem
    Boost::iostreams GSL::gsl  OpenMP::OpenMP_CXX ${CMAKE_SOURCE_DIR}/install/lib/libSPF.a Eigen3::Eigen
    ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES}
    ${CMAKE_SOURCE_DIR}/install/lib/libsundials_nvecserial.a
    ${CMAKE_SOURCE_DIR}/install/lib/libsundials_kinsol.a
    ${CMAKE_SOURCE_DIR}/install/lib/libsundials_idas.a
    ${CMAKE_SOURCE_DIR}/install/lib/libsundials_cvodes.a
    TBB::tbb
    #${STAN_LIBRARIES}
    )

target_link_libraries(libTumorCoal
    ${Boost_LIBRARIES}  ${CMAKE_SOURCE_DIR}/install/lib/libpll_algorithm.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_optimize.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_tree.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_util.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_msa.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_binary.a ${CMAKE_SOURCE_DIR}/install/lib/libpll.a ${Boost_LIBRARIES} Boost::date_time Boost::filesystem
    Boost::iostreams  GSL::gsl OpenMP::OpenMP_CXX ${CMAKE_SOURCE_DIR}/install/lib/libSPF.a Eigen3::Eigen
    ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES}
    ${CMAKE_SOURCE_DIR}/install/lib/libsundials_nvecserial.a
    ${CMAKE_SOURCE_DIR}/install/lib/libsundials_kinsol.a
    ${CMAKE_SOURCE_DIR}/install/lib/libsundials_idas.a
    ${CMAKE_SOURCE_DIR}/install/lib/libsundials_cvodes.a
    TBB::tbb
    #${STAN_LIBRARIES}
    )

INSTALL(TARGETS run DESTINATION ${CMAKE_SOURCE_DIR}/)
INSTALL(TARGETS libTumorCoal DESTINATION bin)

include_directories(${Boost_INCLUDE_DIRS} )


