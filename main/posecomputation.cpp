#include <visp/vpHomogeneousMatrix.h>

#include <visp/vpPoint.h>
#include <visp/vpSubColVector.h>
#include <visp/vpSubMatrix.h>
#include <visp/vpFeaturePoint.h>
#include <visp/vpFeatureBuilder.h>
#include <visp/vpExponentialMap.h>
#include <visp/vpAdaptiveGain.h>
#include <visp/vpIoTools.h>
#include <fstream>

#include <opencv2/calib3d/calib3d.hpp>

#include <vvs.h>
#include <grid_tracker.h>
#include <perspective_camera.h>
#include <distortion_camera.h>
#include <cb_tracker.h>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::stringstream;
using cv::waitKey;
using namespace covis;

int main()
{
// load calibration images from hard drive
    const string base = "../images/";
    const string prefix = "img";

    // init empty vector of detected patterns
    
    vector<Pattern> patterns;
    patterns.clear();
    patterns.reserve(36);


    //GridTracker tracker;      // this tracker detects a 6x6 grid of points
    CBTracker tracker(8,6);     // this one is to be given the chessboard dimension (8x6)

    // read images while the corresponding file exists
    // images are displayed to ensure the detection was performed
    while(true)
    {
        stringstream ss;
        ss << prefix << patterns.size() << ".jpg";
        std::ifstream testfile(base + ss.str());
        if(testfile.good())
        {
            testfile.close();
            Pattern pat;
            pat.im =  cv::imread(base + ss.str());
            tracker.detect(pat.im, pat.point);
            pat.window = ss.str();
            // draw extraction results
            drawSeq(pat.window, pat.im, pat.point);
            patterns.push_back(pat);
            waitKey(0);
        }
        else
            break;
    }

    PerspectiveCamera cam(544.6583785,546.1634166,319.780979,235.3760346);   // not a very good guess

    // initiate virtual visual servoing with inter-point distance and pattern dimensions
    VVS vvs(cam, 0.03, 8, 6);

    // calibrate from all images
    vvs.calibrate(patterns);   
    
    vector<Pattern> patterns_pc;
    patterns_pc.clear();
    patterns_pc.reserve(36);

    vpHomogeneousMatrix M_pose; //extrinsic matrix
    
    // using image 0 for pose computation
    stringstream ss_test;
    ss_test << prefix << 0 << ".jpg";
    Pattern pat_test;
    pat_test.im =  cv::imread(base + ss_test.str());
    tracker.detect(pat_test.im, pat_test.point);
    pat_test.window = ss_test.str();
    drawSeq(pat_test.window, pat_test.im, pat_test.point);
    vvs.computePose(pat_test, M_pose, true);    

    // print results
    cout << "Final pose matrix for image0: " << M_pose << endl;

    // this will wait for a key pressed to stop the program
    waitKey(0);
}