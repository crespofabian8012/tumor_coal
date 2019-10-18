include_directories(/usr/local/opt/llvm/include ${CMAKE_SOURCE_DIR}/install/include/ /usr/local/include/)

file( GLOB SRCS *.cpp *.hpp *.c *.h)
add_executable(run ${SRCS}) # create the executable file
add_library(libTumorCoal STATIC ${SRCS})

find_package (PkgConfig REQUIRED)
find_package(Boost 1.67 COMPONENTS regex system filesystem unit_test_framework program_options REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

target_link_libraries(run ${CMAKE_SOURCE_DIR}/install/lib/libpll_algorithm.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_optimize.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_tree.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_util.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_msa.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_binary.a ${CMAKE_SOURCE_DIR}/install/lib/libpll.a /usr/local/lib/libboost_regex.a)

target_link_libraries(libTumorCoal ${CMAKE_SOURCE_DIR}/install/lib/libpll_algorithm.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_optimize.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_tree.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_util.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_msa.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_binary.a ${CMAKE_SOURCE_DIR}/install/lib/libpll.a /usr/local/lib/libboost_regex.a)

INSTALL(TARGETS run DESTINATION ${CMAKE_SOURCE_DIR}/)
INSTALL(TARGETS libTumorCoal DESTINATION bin)

