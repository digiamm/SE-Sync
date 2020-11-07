#include "SESync/SESync.h"
#include "SESync/SESync_utils.h"

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
  opts.formulation = Formulation::Explicit;
  opts.initialization = Initialization::Chordal;

#ifdef GPERFTOOLS
  ProfilerStart("SE-Sync.prof"); 
#endif

  /// RUN SE-SYNC!
  SESyncResult results = SESync::SESync(measurements, opts);

#ifdef GPERFTOOLS
  ProfilerStop();
#endif

  if (outfile != ""){
    srrg2::ReaderWriterSRRG2 io;
    // ldg populate containers by parsing the original file, ie covariances have to be retrieved
    io.read_g2o_file(argv[1]);
    const size_t rot_offset = (opts.formulation == Formulation::Simplified ? 0 : num_poses);
    io.write_g2o_file(results.xhat, outfile, rot_offset);
  }
}
