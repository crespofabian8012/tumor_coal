# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.15.3/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.15.3/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal

# Include any dependencies generated for this target.
include inference/CMakeFiles/run.dir/depend.make

# Include the progress variables for this target.
include inference/CMakeFiles/run.dir/progress.make

# Include the compile flags for this target's objects.
include inference/CMakeFiles/run.dir/flags.make

inference/CMakeFiles/run.dir/data_utils.cpp.o: inference/CMakeFiles/run.dir/flags.make
inference/CMakeFiles/run.dir/data_utils.cpp.o: inference/data_utils.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object inference/CMakeFiles/run.dir/data_utils.cpp.o"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/run.dir/data_utils.cpp.o -c /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/data_utils.cpp

inference/CMakeFiles/run.dir/data_utils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/run.dir/data_utils.cpp.i"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/data_utils.cpp > CMakeFiles/run.dir/data_utils.cpp.i

inference/CMakeFiles/run.dir/data_utils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/run.dir/data_utils.cpp.s"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/data_utils.cpp -o CMakeFiles/run.dir/data_utils.cpp.s

inference/CMakeFiles/run.dir/eigen.cpp.o: inference/CMakeFiles/run.dir/flags.make
inference/CMakeFiles/run.dir/eigen.cpp.o: inference/eigen.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object inference/CMakeFiles/run.dir/eigen.cpp.o"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/run.dir/eigen.cpp.o -c /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/eigen.cpp

inference/CMakeFiles/run.dir/eigen.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/run.dir/eigen.cpp.i"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/eigen.cpp > CMakeFiles/run.dir/eigen.cpp.i

inference/CMakeFiles/run.dir/eigen.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/run.dir/eigen.cpp.s"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/eigen.cpp -o CMakeFiles/run.dir/eigen.cpp.s

inference/CMakeFiles/run.dir/main.cpp.o: inference/CMakeFiles/run.dir/flags.make
inference/CMakeFiles/run.dir/main.cpp.o: inference/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object inference/CMakeFiles/run.dir/main.cpp.o"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/run.dir/main.cpp.o -c /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/main.cpp

inference/CMakeFiles/run.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/run.dir/main.cpp.i"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/main.cpp > CMakeFiles/run.dir/main.cpp.i

inference/CMakeFiles/run.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/run.dir/main.cpp.s"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/main.cpp -o CMakeFiles/run.dir/main.cpp.s

inference/CMakeFiles/run.dir/mcmc_chain.cpp.o: inference/CMakeFiles/run.dir/flags.make
inference/CMakeFiles/run.dir/mcmc_chain.cpp.o: inference/mcmc_chain.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object inference/CMakeFiles/run.dir/mcmc_chain.cpp.o"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/run.dir/mcmc_chain.cpp.o -c /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/mcmc_chain.cpp

inference/CMakeFiles/run.dir/mcmc_chain.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/run.dir/mcmc_chain.cpp.i"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/mcmc_chain.cpp > CMakeFiles/run.dir/mcmc_chain.cpp.i

inference/CMakeFiles/run.dir/mcmc_chain.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/run.dir/mcmc_chain.cpp.s"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/mcmc_chain.cpp -o CMakeFiles/run.dir/mcmc_chain.cpp.s

inference/CMakeFiles/run.dir/mcmc_move.cpp.o: inference/CMakeFiles/run.dir/flags.make
inference/CMakeFiles/run.dir/mcmc_move.cpp.o: inference/mcmc_move.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object inference/CMakeFiles/run.dir/mcmc_move.cpp.o"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/run.dir/mcmc_move.cpp.o -c /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/mcmc_move.cpp

inference/CMakeFiles/run.dir/mcmc_move.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/run.dir/mcmc_move.cpp.i"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/mcmc_move.cpp > CMakeFiles/run.dir/mcmc_move.cpp.i

inference/CMakeFiles/run.dir/mcmc_move.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/run.dir/mcmc_move.cpp.s"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/mcmc_move.cpp -o CMakeFiles/run.dir/mcmc_move.cpp.s

inference/CMakeFiles/run.dir/mutationModel.cpp.o: inference/CMakeFiles/run.dir/flags.make
inference/CMakeFiles/run.dir/mutationModel.cpp.o: inference/mutationModel.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object inference/CMakeFiles/run.dir/mutationModel.cpp.o"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/run.dir/mutationModel.cpp.o -c /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/mutationModel.cpp

inference/CMakeFiles/run.dir/mutationModel.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/run.dir/mutationModel.cpp.i"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/mutationModel.cpp > CMakeFiles/run.dir/mutationModel.cpp.i

inference/CMakeFiles/run.dir/mutationModel.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/run.dir/mutationModel.cpp.s"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/mutationModel.cpp -o CMakeFiles/run.dir/mutationModel.cpp.s

inference/CMakeFiles/run.dir/output_functions.cpp.o: inference/CMakeFiles/run.dir/flags.make
inference/CMakeFiles/run.dir/output_functions.cpp.o: inference/output_functions.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object inference/CMakeFiles/run.dir/output_functions.cpp.o"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/run.dir/output_functions.cpp.o -c /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/output_functions.cpp

inference/CMakeFiles/run.dir/output_functions.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/run.dir/output_functions.cpp.i"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/output_functions.cpp > CMakeFiles/run.dir/output_functions.cpp.i

inference/CMakeFiles/run.dir/output_functions.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/run.dir/output_functions.cpp.s"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/output_functions.cpp -o CMakeFiles/run.dir/output_functions.cpp.s

inference/CMakeFiles/run.dir/population.cpp.o: inference/CMakeFiles/run.dir/flags.make
inference/CMakeFiles/run.dir/population.cpp.o: inference/population.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object inference/CMakeFiles/run.dir/population.cpp.o"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/run.dir/population.cpp.o -c /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/population.cpp

inference/CMakeFiles/run.dir/population.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/run.dir/population.cpp.i"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/population.cpp > CMakeFiles/run.dir/population.cpp.i

inference/CMakeFiles/run.dir/population.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/run.dir/population.cpp.s"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/population.cpp -o CMakeFiles/run.dir/population.cpp.s

inference/CMakeFiles/run.dir/random.cpp.o: inference/CMakeFiles/run.dir/flags.make
inference/CMakeFiles/run.dir/random.cpp.o: inference/random.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object inference/CMakeFiles/run.dir/random.cpp.o"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/run.dir/random.cpp.o -c /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/random.cpp

inference/CMakeFiles/run.dir/random.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/run.dir/random.cpp.i"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/random.cpp > CMakeFiles/run.dir/random.cpp.i

inference/CMakeFiles/run.dir/random.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/run.dir/random.cpp.s"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/random.cpp -o CMakeFiles/run.dir/random.cpp.s

inference/CMakeFiles/run.dir/tree_node.cpp.o: inference/CMakeFiles/run.dir/flags.make
inference/CMakeFiles/run.dir/tree_node.cpp.o: inference/tree_node.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object inference/CMakeFiles/run.dir/tree_node.cpp.o"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/run.dir/tree_node.cpp.o -c /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/tree_node.cpp

inference/CMakeFiles/run.dir/tree_node.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/run.dir/tree_node.cpp.i"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/tree_node.cpp > CMakeFiles/run.dir/tree_node.cpp.i

inference/CMakeFiles/run.dir/tree_node.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/run.dir/tree_node.cpp.s"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/tree_node.cpp -o CMakeFiles/run.dir/tree_node.cpp.s

inference/CMakeFiles/run.dir/utils.cpp.o: inference/CMakeFiles/run.dir/flags.make
inference/CMakeFiles/run.dir/utils.cpp.o: inference/utils.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object inference/CMakeFiles/run.dir/utils.cpp.o"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/run.dir/utils.cpp.o -c /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/utils.cpp

inference/CMakeFiles/run.dir/utils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/run.dir/utils.cpp.i"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/utils.cpp > CMakeFiles/run.dir/utils.cpp.i

inference/CMakeFiles/run.dir/utils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/run.dir/utils.cpp.s"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/utils.cpp -o CMakeFiles/run.dir/utils.cpp.s

# Object files for target run
run_OBJECTS = \
"CMakeFiles/run.dir/data_utils.cpp.o" \
"CMakeFiles/run.dir/eigen.cpp.o" \
"CMakeFiles/run.dir/main.cpp.o" \
"CMakeFiles/run.dir/mcmc_chain.cpp.o" \
"CMakeFiles/run.dir/mcmc_move.cpp.o" \
"CMakeFiles/run.dir/mutationModel.cpp.o" \
"CMakeFiles/run.dir/output_functions.cpp.o" \
"CMakeFiles/run.dir/population.cpp.o" \
"CMakeFiles/run.dir/random.cpp.o" \
"CMakeFiles/run.dir/tree_node.cpp.o" \
"CMakeFiles/run.dir/utils.cpp.o"

# External object files for target run
run_EXTERNAL_OBJECTS =

inference/run: inference/CMakeFiles/run.dir/data_utils.cpp.o
inference/run: inference/CMakeFiles/run.dir/eigen.cpp.o
inference/run: inference/CMakeFiles/run.dir/main.cpp.o
inference/run: inference/CMakeFiles/run.dir/mcmc_chain.cpp.o
inference/run: inference/CMakeFiles/run.dir/mcmc_move.cpp.o
inference/run: inference/CMakeFiles/run.dir/mutationModel.cpp.o
inference/run: inference/CMakeFiles/run.dir/output_functions.cpp.o
inference/run: inference/CMakeFiles/run.dir/population.cpp.o
inference/run: inference/CMakeFiles/run.dir/random.cpp.o
inference/run: inference/CMakeFiles/run.dir/tree_node.cpp.o
inference/run: inference/CMakeFiles/run.dir/utils.cpp.o
inference/run: inference/CMakeFiles/run.dir/build.make
inference/run: install/lib/libpll_algorithm.a
inference/run: install/lib/libpll_optimize.a
inference/run: install/lib/libpll_tree.a
inference/run: install/lib/libpll_util.a
inference/run: install/lib/libpll_msa.a
inference/run: install/lib/libpll_binary.a
inference/run: install/lib/libpll.a
inference/run: /usr/local/lib/libboost_regex.a
inference/run: /usr/local/lib/libgsl.a
inference/run: inference/CMakeFiles/run.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Linking CXX executable run"
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/run.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
inference/CMakeFiles/run.dir/build: inference/run

.PHONY : inference/CMakeFiles/run.dir/build

inference/CMakeFiles/run.dir/clean:
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference && $(CMAKE_COMMAND) -P CMakeFiles/run.dir/cmake_clean.cmake
.PHONY : inference/CMakeFiles/run.dir/clean

inference/CMakeFiles/run.dir/depend:
	cd /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference /Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/inference/CMakeFiles/run.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : inference/CMakeFiles/run.dir/depend

