#ifndef PERSPECTIVE_CAMERA_H
#define PERSPECTIVE_CAMERA_H

#include <generic_camera.h>

/*
 * This class inherits from GenericCamera and implements the classical perspective camera model
 * Methods to be filled up are:
 * - project
 * - computeJacobianIntrinsic
 * - computeJacobianExtrinsic
 */


namespace covis
{


// xi = (px, py, u0, v0)
class PerspectiveCamera : public GenericCamera
{
public:

    PerspectiveCamera(const double &_px, const double &_py, const double &_u0, const double &_v0, const bool &_calibrated=false)
    {
        xi_.resize(4);
        xi_[0] = _px;
        xi_[1] = _py;
        xi_[2] = _u0;
        xi_[3] = _v0;
    }


    // check that the parameters stay meaningfull in case of wrong update
    void updateIntrinsic(const vpColVector &_dxi)
    {
        // update
        xi_ += _dxi;
        // all parameters should be positive
        for(unsigned int i=0;i<4;++i){
            if(xi_[i] < 0){
                xi_[i] = 0;
            }
        }

    }


    // compute pixel coordinates of a 3D point
    // we assume the point is already projected in the camera frame
    void project(const vpPoint &_P, double &_u, double &_v)
    {
        _u = xi_[0]*_P.get_X()/_P.get_Z() + xi_[2];
        _v = xi_[1]*_P.get_Y()/_P.get_Z() + xi_[3];

    }


    // write the Jacobian corresponding to the intrinsic parameters
    // Jacobian should be 2x4
    // we assume the point is already projected in the camera frame
    void computeJacobianIntrinsic(const vpPoint &_P, vpMatrix &_J)
    {
        /* Ji = |X/Z 0   1 0|
                |0   Y/Z 0 1|
        */

        _J[0][0] = _P.get_X()/_P.get_Z();
        _J[0][1] = 0;
        _J[0][2] = 1;
        _J[0][3] = 0;

        _J[1][0] = 0;
        _J[1][1] = _P.get_Y()/_P.get_Z();
        _J[1][2] = 0;
        _J[1][3] = 1;

    }


    // write the Jacobian wrt extrinsic parameters
    // J should be 2x6
    // we assume the point is already projected in the camera frame
    void computeJacobianExtrinsic(const vpPoint &_P, vpMatrix &_J)
    {
        

    double X = _P.get_X();
    double Y = _P.get_Y();
    double Z = _P.get_Z();
    double x_ = xi_[0];
    double y_ = xi_[0];
    double xx = X/Z;
    double yy = Y/Z;
        
    _J.resize(2,6);
    _J[0][0] = -x_/Z ;
    _J[0][1] = 0;
    _J[0][2] = x_ * xx/Z;
    _J[0][3] = x_ * xx * yy;
    _J[0][4] = - x_ * (1 + (xx*xx)) ;
    _J[0][5] = x_ * yy;
    
    _J[1][0] = 0;
    _J[1][1] = - y_/Z;
    _J[1][2] = y_ * yy/Z;
    _J[1][3] = y_ * (1+ (yy*yy));
    _J[1][4] = - y_ * xx * yy;
    _J[1][5] = - y_ * xx;

    }
};
}


#endif // PERSPECTIVE_CAMERA_H
