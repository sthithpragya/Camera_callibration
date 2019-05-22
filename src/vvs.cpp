
#include <vvs.h>
#include <visp/vpSubColVector.h>
#include <visp/vpExponentialMap.h>

using namespace covis;
using std::cout;
using std::endl;
using std::vector;
using std::string;

/* Calibrate the camera from sequence of images and extracted points
    _pat is a list of image patterns with:
        - _pat[i].im is the i-eth image
        - _pat[i].point are the pixel points extracted from this image
        - _pat[i].window is the name of the window to display the results

*/
void VVS::calibrate(std::vector<Pattern> &_pat)
{
    // number of images
    const int n = _pat.size();
    // number of points
    const int m = r_*c_;
    // number of parameters
    const int n_xi = cam_->nbParam();

    // the initial guess for M
    vector<vpHomogeneousMatrix> M(n);
    for(unsigned int k=0;k<n;++k)
        M[k].buildFrom(
            0,0.,0.5,                                                                                   // translation
            0,0,atan2(_pat[k].point[5].y-_pat[k].point[0].y, _pat[k].point[5].x-_pat[k].point[0].x));   // rotation

    // current and desired features: 2 coord. for each m points for each n images
    vpColVector s(2*m*n), sd(2*m*n);

    // write measured positions as desired features
    unsigned int row = 0;

    unsigned int col = 0; // Declared for subsequent use

    unsigned int k = 0;                     // image number
    for(auto &pat: _pat)                    // loop through all images
    {
        for(unsigned int i=0;i<m;++i)       // loop through all points
        {
            row = 2*i+2*m*k;                // corresponding row in the global error vector
            sd[row] = pat.point[i].x;       // x value
            sd[row+1] =  pat.point[i].y;    // y value
        }
        k++;
    }

    // update vector
    vpColVector dx(n_xi + 6*n);      // intrinsic parameters + 6 velocities for each image
    vpSubColVector dxi(dx, 0, n_xi); // part of dx that represents the intrinsic parameters

    // Jacobians
    vpMatrix J(2*m*n, n_xi+6*n);     // global Jacobian
    vpMatrix Ji(2, n_xi);            // Jacobian for 1 point for intrinsic
    vpMatrix Li(2, 6);               // Jacobian for 1 point for extrinsic

    unsigned int iter=0;

    dx = 1; // set to 1 so that it is higher than minimum error at first
    while(dx.euclideanNorm() > 0.00001 && iter++ < 100)
    {
        /* first we have to compute for all points from all images:
         * - the pixel coordinates corresponding to the current estimation of intrinsic and extrinsic
         * - the intrinsic Jacobian
         * - the extrinsic Jacobian
         */

        unsigned int k = 0;                     // image number
        for(auto &pat: _pat)                    // loop through all images with index k
        {
            for(unsigned int i=0;i<m;++i)       // loop through all points with index i
            {
                // corresponding row in the global error vector
                unsigned int ROW = 2*i+2*m*k;

                // track the point X_[i] from the current estimation of M_k


                vpPoint world_vpPoint = X_[i]; 
                vpPoint camera_vpPoint = world_vpPoint;
                camera_vpPoint.track(M[k]);

                double Ui;
                double Vi;

                // use the current intrinsic estimation to project into pixels

                cam_->project(camera_vpPoint, Ui, Vi);

                s[ROW] = Ui;
                s[ROW+1] = Vi;

                cam_->computeJacobianIntrinsic(camera_vpPoint, Ji);

                cam_->computeJacobianExtrinsic(camera_vpPoint, Li);

                //col = 4 + 6*i+6*m*k;
                putAt(J, Ji, ROW, 0);
                putAt(J, Li, ROW, (4+6*k));

            }
            k++;
        }


        /* Here we have the error (s-sd) and the corresponding Jacobian J
         * We can compute an update for the unknown with the pseudo inverse
         */


        // dx = ??

        dx = -0.1*J.pseudoInverse()*(s - sd);

        // we now update the intrinsic parameters, common to all images
        cam_->updateIntrinsic(dxi);

        // and the extrinsic parameters for each image
        for(unsigned int k=0;k<n;++k)
        {
            vpSubColVector v(dx, n_xi+6*k, 6);          // extract velocity vector from global update dx
            M[k] = vpExponentialMap::direct(v).inverse() * M[k];
        }

        // if we want to see what's happening during the loop
        display(_pat, M, s, 5);

    }

    // here the camera should be calibrated and all poses M_k should be estimated, display the result
    display(_pat, M, s, 0);

}


/* Compute the pose of the camera from a single image and extracted points and a guess on the current pose M
    _pat is an image patterns with:
        - _pat.im is the image
        - _pat.point are the pixel points extracted from this image
        - _pat.window is the name of the window to display the results

    _reset can be used as an external way to tell the VVS not to use the last M found but to reinitialize one from a wild guess
            (used when the pose was badly calibrated)
*/

void VVS::computePose(Pattern &_pat, vpHomogeneousMatrix &_M, const bool &_reset)
{   
    vpHomogeneousMatrix M_guess; 
    M_guess.buildFrom( //guess estimate from the image
            0,0.,0.5,                                                                                   // translation
            0,0,atan2(_pat.point[5].y-_pat.point[0].y, _pat.point[5].x-_pat.point[0].x));   // rotation

    if(_reset == true){ //if we want to reinitialize M from a new guess
        _M = M_guess;
    }

    // number of points
    const int m = r_*c_;
    vpColVector s(2*m), sd(2*m);

    unsigned int row = 0;

    //making sd vector of known pixel coordinates
    for(unsigned int i=0;i<m;++i){       // loop through all points
        row = 2*i;                // corresponding row in the global error vector
        sd[row] = _pat.point[i].x;       // x value
        sd[row+1] =  _pat.point[i].y;    // y value
    }    


    vpColVector dx(6);      // screw of pose
    vpMatrix L(2*m, 6); // extrinsic jacobian matrix
    vpMatrix Li(2,6);

    unsigned int iter=0;

    dx = 1; // set to 1 so that it is higher than minimum error at first

    while(dx.euclideanNorm() > 0.00001 && iter++ < 100)
    {
        /* first we have to compute for all points from all images:
         * - the pixel coordinates corresponding to the current estimation of intrinsic and extrinsic
         * - the intrinsic Jacobian
         * - the extrinsic Jacobian
         */

        for(unsigned int i=0;i<m;++i)       // loop through all points with index i
        {
            // corresponding row in the global error vector
            unsigned int ROW = 2*i;

            // track the point X_[i] from the current estimation of M_k


            vpPoint world_vpPoint = X_[i]; 
            vpPoint camera_vpPoint = world_vpPoint;
            camera_vpPoint.track(_M); 

            double Ui;
            double Vi;

            cam_->project(camera_vpPoint, Ui, Vi);

            s[ROW] = Ui;
            s[ROW+1] = Vi;

            cam_->computeJacobianExtrinsic(camera_vpPoint, Li);
            putAt(L, Li, ROW,0); 

        }

        dx = -0.1*L.pseudoInverse()*(s - sd);

        // update the extrinsic parameters for the image

        _M = vpExponentialMap::direct(dx).inverse()*_M;
    }
}
