#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

#include <Eigen/Geometry>
#include <Eigen/StdVector>

#include "SESync/SESync_utils.h"

namespace SESync {
  using Quaterniond      = Eigen::Quaterniond;
  using Vector3d         = Eigen::Vector3d;
  using QuaternionVector = std::vector<Quaterniond, Eigen::aligned_allocator<Quaterniond>>;
  using Vector3dVector   = std::vector<Vector3d, Eigen::aligned_allocator<Vector3d>>;
  Matrix readInitialGuess(const std::string& filename, size_t& num_variables, size_t& num_measurements) {
    QuaternionVector q;
    Vector3dVector t;
    std::string line;
    std::string token;
    size_t id;
    Scalar qx, qy, qz, qw, tx, ty, tz;
    // Open the file for reading
    std::ifstream infile(filename);

    num_variables    = 0;
    num_measurements = 0;

    while (std::getline(infile, line)) {
      // Construct a stream from the string
      std::stringstream strstrm(line);

      // Extract the first token from the string
      strstrm >> token;
      if (token == "VERTEX_SE3:QUAT") {
        strstrm >> id >> tx >> ty >> tz >> qx >> qy >> qz >> qw;
        Quaterniond qi(qw, qx, qy, qz);
        Vector3d ti(tx, ty, tz);
        q.push_back(qi);
        t.push_back(ti);
        ++num_variables;
      }
      if (token == "EDGE_SE3:QUAT") {
        ++num_measurements;
      }
    }
    infile.close();
    Matrix Y = Matrix::Zero(3, 4 * num_variables);
    for (size_t i = 0; i < num_variables; ++i) {
      Eigen::Matrix3d Ri      = q[i].toRotationMatrix();
      Vector3d ti             = t[i];
      Y.block(0, i, 3, 1)     = ti;
      Y.block(0, num_variables + i * 3, 3, 3) = Ri;
    }
    num_variables *= 6;
    num_measurements *= 6;
    return Y;
  }
} // namespace SESync
