# DO NOT EDIT
# This makefile makes sure all linkable targets are
# up-to-date with anything they link to
default:
	echo "Do not invoke directly"

# Rules to remove targets that are older than anything to which they
# link.  This forces Xcode to relink the targets from scratch.  It
# does not seem to check these dependencies itself.
PostBuild.simul.Debug:
/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/tumor_coal/xcode/simul/Debug/simul:\
	/usr/local/lib/libgsl.dylib\
	/usr/local/lib/libgslcblas.dylib\
	/usr/local/lib/libboost_system.a\
	/usr/local/lib/libboost_filesystem.a\
	/usr/local/lib/libboost_program_options.a\
	/usr/local/lib/libpll.a\
	/usr/local/lib/libpll_algorithm.a\
	/usr/local/lib/libpll_tree.a\
	/usr/local/lib/libpll_util.a\
	/usr/local/lib/libpll_optimize.a\
	/usr/local/lib/libpll_msa.a\
	/usr/local/lib/libpll_binary.a
	/bin/rm -f /Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/tumor_coal/xcode/simul/Debug/simul


PostBuild.simul.Release:
/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/tumor_coal/xcode/simul/Release/simul:\
	/usr/local/lib/libgsl.dylib\
	/usr/local/lib/libgslcblas.dylib\
	/usr/local/lib/libboost_system.a\
	/usr/local/lib/libboost_filesystem.a\
	/usr/local/lib/libboost_program_options.a\
	/usr/local/lib/libpll.a\
	/usr/local/lib/libpll_algorithm.a\
	/usr/local/lib/libpll_tree.a\
	/usr/local/lib/libpll_util.a\
	/usr/local/lib/libpll_optimize.a\
	/usr/local/lib/libpll_msa.a\
	/usr/local/lib/libpll_binary.a
	/bin/rm -f /Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/tumor_coal/xcode/simul/Release/simul


PostBuild.simul.MinSizeRel:
/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/tumor_coal/xcode/simul/MinSizeRel/simul:\
	/usr/local/lib/libgsl.dylib\
	/usr/local/lib/libgslcblas.dylib\
	/usr/local/lib/libboost_system.a\
	/usr/local/lib/libboost_filesystem.a\
	/usr/local/lib/libboost_program_options.a\
	/usr/local/lib/libpll.a\
	/usr/local/lib/libpll_algorithm.a\
	/usr/local/lib/libpll_tree.a\
	/usr/local/lib/libpll_util.a\
	/usr/local/lib/libpll_optimize.a\
	/usr/local/lib/libpll_msa.a\
	/usr/local/lib/libpll_binary.a
	/bin/rm -f /Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/tumor_coal/xcode/simul/MinSizeRel/simul


PostBuild.simul.RelWithDebInfo:
/Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/tumor_coal/xcode/simul/RelWithDebInfo/simul:\
	/usr/local/lib/libgsl.dylib\
	/usr/local/lib/libgslcblas.dylib\
	/usr/local/lib/libboost_system.a\
	/usr/local/lib/libboost_filesystem.a\
	/usr/local/lib/libboost_program_options.a\
	/usr/local/lib/libpll.a\
	/usr/local/lib/libpll_algorithm.a\
	/usr/local/lib/libpll_tree.a\
	/usr/local/lib/libpll_util.a\
	/usr/local/lib/libpll_optimize.a\
	/usr/local/lib/libpll_msa.a\
	/usr/local/lib/libpll_binary.a
	/bin/rm -f /Users/seonghwanjun/Dropbox/Research/single-cell-research/repos/tumor_coal/xcode/simul/RelWithDebInfo/simul




# For each target create a dummy ruleso the target does not have to exist
/usr/local/lib/libboost_filesystem.a:
/usr/local/lib/libboost_program_options.a:
/usr/local/lib/libboost_system.a:
/usr/local/lib/libgsl.dylib:
/usr/local/lib/libgslcblas.dylib:
/usr/local/lib/libpll.a:
/usr/local/lib/libpll_algorithm.a:
/usr/local/lib/libpll_binary.a:
/usr/local/lib/libpll_msa.a:
/usr/local/lib/libpll_optimize.a:
/usr/local/lib/libpll_tree.a:
/usr/local/lib/libpll_util.a:
