#include "SESync/SESync_utils.h"
#include <fstream>

const std::string prefix = "srrg2_parser";
#define LOG std::cerr << prefix + "|"

namespace srrg2{
    using namespace SESync;

    struct RelativePoseMeasurementSRRG2{
        size_t i = 0;
        size_t j = 0;
        Matrix R;
        Vector t;
        Matrix trans_cov;
        Matrix rot_cov;
        float kappa = 0.f;

        /** A utility function for streaming Nodes to cout */
        inline friend std::ostream &
        operator<<(std::ostream &os, const RelativePoseMeasurementSRRG2 &measurement){
            os << "i: " << measurement.i << std::endl;
            os << "j: " << measurement.j << std::endl;
            os << "R: " << std::endl
               << measurement.R << std::endl;
            os << "t: " << std::endl
               << measurement.t << std::endl;
            os << "trans cov: " << measurement.trans_cov << std::endl;
            os << "rot cov: " << measurement.rot_cov << std::endl;
            os << "kappa: " << measurement.kappa << std::endl;
            return os;
        }
    };

    class ReaderWriterSRRG2{

    public:
        ReaderWriterSRRG2(){
            LOG << "prepare to IO g2o file!" << std::endl;
        }

        inline void write_g2o_file(const Matrix &solution_, const std::string outfile_, const size_t rot_offset_){
            // ldg some retarded sanity checks
            if (_num_poses == 0 && _num_measures == 0 && _dim == 0){
                throw std::runtime_error("fields not set, first read g2o file in order to populate");
            }

            if (outfile_ == ""){
                throw std::runtime_error("set name of .g2o output file!");
            }

            // ldg prepare output buffer
            std::ofstream g2o_writer(outfile_);

            // ldg get edge type
            std::string vertex_type = (_dim == 3 ? "VERTEX_SE3:QUAT" : "VERTEX_SE2");

            for (size_t i = 0; i < _num_poses; ++i){
                // ldg write invariant parts
                g2o_writer << vertex_type << " " << i << " ";

                const Vector &translation = solution_.block(0, i, _dim, 1);
                const Matrix &rot = solution_.block(0, rot_offset_ + i * _dim, _dim, _dim);

                g2o_writer << translation(0) << " " << translation(1) << " ";

                // ldg depending on edge dump shit
                if (_dim == 3){
                    // casting to so3 type
                    const Eigen::Matrix<Scalar, 3, 3> rot3 = rot;
                    // ldg convert current rotation to quaternion qw, qx, qy, qz
                    Eigen::Quaternion<Scalar> quat(rot3);
                    g2o_writer << translation(2) << " " << quat.coeffs().transpose() << "\n";
                }
                else{
                    const Scalar theta = std::atan2(rot(1, 0), rot(0, 0));
                    g2o_writer << theta << "\n";
                }
            }

            // ldg now edges
            std::string edge_type = (_dim == 3 ? "EDGE_SE3:QUAT" : "EDGE_SE2");
            const size_t measurements_size = _measurements.size();
            for (size_t z = 0; z < measurements_size; ++z){
                // get current meas
                const RelativePoseMeasurementSRRG2 &measurement = _measurements[z];
                // ldg write invariant parts
                g2o_writer << edge_type << " " << measurement.i << " " << measurement.j << " "
                           << measurement.t(0) << " " << measurement.t(1) << " ";

                if (_dim == 3){
                    // ldg convert current rotation to quaternion qw, qx, qy, qz
                    // casting to so3 type
                    const Eigen::Matrix<Scalar, 3, 3> z_rot3 = measurement.R;
                    Eigen::Quaternion<Scalar> quat(z_rot3);

                    // g2o output dqx, dqy, dqz, dqw
                    g2o_writer << measurement.t(2) << " " << quat.coeffs().transpose() << " ";
                    // ldg add covariance
                    const Eigen::Matrix<Scalar, 3, 3> tcov = measurement.trans_cov;
                    const Eigen::Matrix<Scalar, 3, 3> rcov = measurement.rot_cov;

                    // clang-format off
                    // covariance diagonal matrix in our case
                    g2o_writer << tcov(0, 0) << " " << tcov(0, 1) << " " << tcov(0, 2) << " " << 0   << " " << 0 << " " << 0 << " "
                                                    << tcov(1, 1) << " " << tcov(1, 2) << " " << 0   << " " << 0 << " " << 0 << " "
                                                                         << tcov(2, 2) << " " << 0   << " " << 0 << " " << 0 << " "
                                                                                       << rcov(0, 0) << " " << rcov(0, 1) << " " << rcov(0, 2) << " "
                                                                                                            << rcov(1, 1) << " " << rcov(1, 2) << " "
                                                                                                                                 << rcov(2, 2) << "\n";
                    // clang-format on
                }else{
                    const Scalar theta = std::atan2(measurement.R(1, 0), measurement.R(0, 0));
                    const Matrix &tcov = measurement.trans_cov;
                    g2o_writer << theta << " " << tcov(0, 0) << " " << tcov(0, 1) << " " << tcov(0, 2)
                                        << " " << tcov(1, 1) << " " << tcov(1, 2) << " " << measurement.kappa << "\n";
                }
            }
            g2o_writer.close();
            LOG << "pose-graph optimized output written on .g2o file: " << outfile_ << std::endl;
        }

        // this is needed in order to retrieve covariances without
        // afftecting original author code!
        inline void read_g2o_file(const std::string infile_){

            _measurements.clear();
            RelativePoseMeasurementSRRG2 measurement;

            // A string used to contain the contents of a single line
            std::string line;

            // A string used to extract tokens from each line one-by-one
            std::string token;

            // Preallocate various useful quantities
            Scalar dtheta = 0.0;
            Scalar dx, dy, dz, dqx, dqy, dqz, dqw, I11, I12, I13, I14, I15, I16,
                I22, I23, I24, I25, I26, I33, I34, I35, I36, I44, I45, I46, I55, I56, I66;

            size_t i, j;

            // Open the file for reading
            std::ifstream infile(infile_);

            _num_poses = 0;

            while (std::getline(infile, line)){
                // Construct a stream from the string
                std::stringstream strstrm(line);

                // Extract the first token from the string
                strstrm >> token;

                if (token == "EDGE_SE2"){

                    // This is a 2D pose measurement
                    // Extract formatted output
                    strstrm >> i >> j >> dx >> dy >> dtheta >> I11 >> I12 >> I13 >> I22 >> I23 >> I33;

                    // Fill in elements of this measurement

                    // Pose ids
                    measurement.i = i;
                    measurement.j = j;

                    // Raw measurements
                    measurement.t = Eigen::Matrix<Scalar, 2, 1>(dx, dy);
                    measurement.R = Eigen::Rotation2D<Scalar>(dtheta).toRotationMatrix();

                    measurement.trans_cov = Eigen::Matrix<Scalar, 2, 2>();
                    measurement.trans_cov << I11, I12, I12, I22;

                    measurement.kappa = I33;
                }else if (token == "EDGE_SE3:QUAT"){
                    // This is a 3D pose measurement

                    // Extract formatted output
                    strstrm >> i >> j >> dx >> dy >> dz >> dqx >> dqy >> dqz >> dqw >> I11 >>
                        I12 >> I13 >> I14 >> I15 >> I16 >> I22 >> I23 >> I24 >> I25 >> I26 >>
                        I33 >> I34 >> I35 >> I36 >> I44 >> I45 >> I46 >> I55 >> I56 >> I66;

                    // Fill in elements of the measurement

                    // Pose ids
                    measurement.i = i;
                    measurement.j = j;

                    // Raw measurements
                    measurement.t = Eigen::Matrix<Scalar, 3, 1>(dx, dy, dz);
                    measurement.R = Eigen::Quaternion<Scalar>(dqw, dqx, dqy, dqz).toRotationMatrix();

                    // store covariances
                    measurement.trans_cov = Eigen::Matrix<Scalar, 3, 3>();
                    measurement.trans_cov.setZero();
                    measurement.trans_cov << I11, I12, I13, I12, I22, I23, I13, I23, I33;

                    measurement.rot_cov = Eigen::Matrix<Scalar, 3, 3>();
                    measurement.rot_cov.setZero();
                    measurement.rot_cov << I44, I45, I46, I45, I55, I56, I46, I56, I66;
                }else if ((token == "VERTEX_SE2") || (token == "VERTEX_SE3:QUAT")){
                    // This is just initialization information, so do nothing
                    continue;
                }else{
                    LOG << "error: unrecognized type: " << token << "!" << std::endl;
                    assert(false);
                }

                // Update maximum value of poses found so far
                size_t max_pair = std::max<size_t>(measurement.i, measurement.j);

                _num_poses = ((max_pair > _num_poses) ? max_pair : _num_poses);
                _measurements.emplace_back(measurement);
            } // while

            infile.close();

            _num_poses++; // account for the use of zero-based indexing
            _dim = (dtheta == Scalar(0) ? 3 : 2);
            // store some final values
            _num_measures = _measurements.size();

            LOG << "reparsed .g2o input to get original covariances, input file: " << infile_ << std::endl;
        }

    private:
        size_t _num_poses = 0;
        size_t _num_measures = 0;
        size_t _dim = 0;
        std::vector<RelativePoseMeasurementSRRG2> _measurements;
    };

} // namespace srrg2
