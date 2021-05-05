#include "SESync/SESync.h"
#include "SESync/SESync_utils.h"

#include "read_initial_guess.h"
#include "srrg2_parser.h"

#ifdef GPERFTOOLS
#include <gperftools/profiler.h>
#endif

using namespace std;
using namespace SESync;


int main(int argc, char **argv) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " [input .g2o file]" << endl;
    exit(1);
  }

  string outfile = "";
  if (argc == 3){
    cout << "output file selected: " << argv[2] << endl;
    outfile = argv[2];
  }

  size_t num_poses;
  measurements_t measurements = read_g2o_file(argv[1], num_poses);
  cout << "Loaded " << measurements.size() << " measurements between "
       << num_poses << " poses from file " << argv[1] << endl
       << endl;

  if (measurements.size() == 0){
    cout << "Error: No measurements were read!"
         << " Are you sure the file exists?" << endl;
    exit(1);
  }

  SESyncOpts opts;
  opts.verbose = true; // Print output to stdout
  opts.num_threads = 4;
  opts.r0             = 3;
  opts.formulation    = Formulation::Explicit;
  opts.initialization = Initialization::Chordal;
  size_t n = 0, m = 0;
  Matrix Y0 = readInitialGuess(argv[1], n, m);
  // Matrix Y0;

#ifdef GPERFTOOLS
    ProfilerStart("SE-Sync.prof"); 
#endif

  /// RUN SE-SYNC!
    SESyncResult results = SESync::SESync(measurements, opts, Y0);

#ifdef GPERFTOOLS
  ProfilerStop();
#endif
  std::ofstream stats("sesync_stats.txt");
  int total_iterations = 0;
  for (const auto& level_iterations : results.function_values) {
    total_iterations += level_iterations.size();
  }
  float normalizer = m - n;
  stats << results.SDPval / normalizer << " " << total_iterations << " " << results.total_computation_time << std::endl;

  if (outfile != ""){
    srrg2::ReaderWriterSRRG2 io;
    // ldg populate containers by parsing the original file, ie covariances have to be retrieved
    io.read_g2o_file(argv[1]);
    const size_t rot_offset = (opts.formulation == Formulation::Simplified ? 0 : num_poses);
    io.write_g2o_file(results.xhat, outfile, rot_offset);
  }
}
