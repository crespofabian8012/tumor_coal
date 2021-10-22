
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++17 -stdlib=libc++")
include_directories(/usr/local/opt/llvm/include ${CMAKE_SOURCE_DIR}/install/include/ /usr/local/include/
/usr/local/include/gsl )

link_directories(${CMAKE_SOURCE_DIR}/install/lib/)

link_directories(/usr/local/lib)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}")

find_package(Boost 1.67 COMPONENTS regex unit_test_framework REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
endif()

file( GLOB TESTS *.cpp *.hpp *.c *.h)



set ( PROJECT_LINK_LIBS libSPF.a )
link_directories(/usr/local/bin)
include_directories(/usr/local/include)
include_directories(/usr/local/include/spf)


find_package(GSL REQUIRED)

find_package (PkgConfig REQUIRED)
find_package (GSL REQUIRED)
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package (Eigen3 3.4 REQUIRED NO_MODULE)

#find_package(Boost 1.67 COMPONENTS regex system filesystem unit_test_framework program_options REQUIRED)
find_package(Boost 1.67 COMPONENTS  system filesystem unit_test_framework program_options date_time filesystem iostreams REQUIRED)
if(Boost_FOUND)
    include_directories( PRIVATE ${Boost_INCLUDE_DIRS} ${Boost_LIBRARIES})
endif()


include_directories(simul /usr/local/include/ ${CMAKE_SOURCE_DIR}/inference/ ${CMAKE_SOURCE_DIR}/install/include/ /usr/local/opt/llvm/include)

add_executable(tests ${TESTS})


target_link_libraries(tests gtest gtest_main  libTumorCoal )

add_test(SMCTests tests)
