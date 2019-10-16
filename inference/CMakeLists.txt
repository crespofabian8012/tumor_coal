include_directories(/usr/local/opt/llvm/include ${CMAKE_SOURCE_DIR}/install/include/)

file( GLOB SRCS *.cpp *.hpp *.c *.h)
add_executable(run ${SRCS}) # create the executable file
add_library(libTumorCoal STATIC ${SRCS})

find_package (PkgConfig REQUIRED)

target_link_libraries(run ${CMAKE_SOURCE_DIR}/install/lib/libpll_algorithm.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_optimize.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_tree.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_util.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_msa.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_binary.a ${CMAKE_SOURCE_DIR}/install/lib/libpll.a)

target_link_libraries(libTumorCoal ${CMAKE_SOURCE_DIR}/install/lib/libpll_algorithm.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_optimize.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_tree.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_util.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_msa.a ${CMAKE_SOURCE_DIR}/install/lib/libpll_binary.a ${CMAKE_SOURCE_DIR}/install/lib/libpll.a)

INSTALL(TARGETS run DESTINATION ${CMAKE_SOURCE_DIR}/)
INSTALL(TARGETS libTumorCoal DESTINATION bin)

