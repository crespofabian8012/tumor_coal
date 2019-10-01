# code relating to tree structured stick breaking process
include_directories(/usr/local/include/)

file( GLOB SRCS *.cpp *.hpp)
add_executable(simul ${SRCS}) # create the executable file

find_package (PkgConfig REQUIRED)
find_package(Boost 1.67 COMPONENTS regex system filesystem unit_test_framework program_options REQUIRED)
find_package (GSL REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

target_link_libraries(simul GSL::gsl GSL::gslcblas /usr/local/lib/libboost_system.a /usr/local/lib/libboost_filesystem.a /usr/local/lib/libboost_program_options.a)

INSTALL(TARGETS simul DESTINATION ${CMAKE_SOURCE_DIR}/)
