# code relating to tree structured stick breaking process
link_directories(${CMAKE_SOURCE_DIR}/install/lib/)

file( GLOB SRCS *.cpp *.hpp *.c *.h)
add_executable(simul ${SRCS}) # create the executable file

find_package (PkgConfig REQUIRED)

include_directories(simul /usr/local/include/ ${CMAKE_SOURCE_DIR}/inference/ ${CMAKE_SOURCE_DIR}/install/include/ /usr/local/opt/llvm/include)
target_link_libraries(simul libTumorCoal)

INSTALL(TARGETS simul DESTINATION ${CMAKE_SOURCE_DIR}/)
