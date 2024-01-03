## List of Known Dependencies
### ORB-SLAM3 v1.1 (Modernized 2025)

In this document we list all the pieces of code included by ORB-SLAM3 and linked libraries which are not property of the authors of ORB-SLAM3.

**Note**: This version has been modernized to use C++23 standard and latest library versions. All dependencies have been updated accordingly while maintaining compatibility and original licensing.

##### Code in **src** and **include** folders

* *ORBextractor.cc*.
This is a modified version of orb.cpp of OpenCV library. The original code is BSD licensed.

* *PnPsolver.h, PnPsolver.cc*.
This is a modified version of the epnp.h and epnp.cc of Vincent Lepetit. 
This code can be found in popular BSD licensed computer vision libraries as [OpenCV](https://github.com/Itseez/opencv/blob/master/modules/calib3d/src/epnp.cpp) and [OpenGV](https://github.com/laurentkneip/opengv/blob/master/src/absolute_pose/modules/Epnp.cpp). The original code is FreeBSD.

* *MLPnPsolver.h, MLPnPsolver.cc*.
This is a modified version of the MLPnP of Steffen Urban from [here](https://github.com/urbste/opengv). 
The original code is BSD licensed.

* Function *ORBmatcher::DescriptorDistance* in *ORBmatcher.cc*.
The code is from: http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel.
The code is in the public domain.

##### Code in Thirdparty folder

* All code in **DBoW2** folder.
This is a modified version of [DBoW2](https://github.com/dorian3d/DBoW2) and [DLib](https://github.com/dorian3d/DLib) library. All files included are BSD licensed.
**Modernization**: Updated for C++23 compatibility with modern CMake build system.

* All code in **g2o** folder.
This is a modified version of [g2o](https://github.com/RainerKuemmerle/g2o). All files included are BSD licensed.
**Modernization**: Updated for C++23 standard with modern compiler flags and dependencies.

* All code in **Sophus** folder.
This is a modified version of [Sophus](https://github.com/strasdat/Sophus). [MIT license](https://en.wikipedia.org/wiki/MIT_License).
**Modernization**: Updated for C++23 compatibility with warning suppressions for legacy template usage.

##### Library dependencies 

* **Pangolin (visualization and user interface)**.
[MIT license](https://en.wikipedia.org/wiki/MIT_License). 
**Modernized version**: Built from latest source (2025) for C++23 compatibility.

* **OpenCV**.
BSD license.
**Modernized version**: 4.6.0+ with modern header includes (opencv2/opencv.hpp). Updated from legacy OpenCV 3.x headers for improved performance and C++23 compatibility.

* **Eigen3**.
For versions greater than 3.1.1 is MPL2, earlier versions are LGPLv3.
**Modernized version**: 3.4.0+ (tested with system packages) for C++23 template compatibility.

* **CMake**.
**Modernized requirement**: 3.16+ for C++23 standard support and modern build features.

* **C++ Compiler**.
**Modernized requirement**: GCC 11+ or Clang 14+ with full C++23 standard support. Upgraded from original C++11 requirement.

* **Boost (Serialization)**.
**Modernized addition**: Latest Boost serialization libraries required for updated DBoW2 dependencies.

* **ROS (Optional, only if you build Examples/ROS)**.
BSD license. In the manifest.xml the only declared package dependencies are roscpp, tf, sensor_msgs, image_transport, cv_bridge, which are all BSD licensed.
**Note**: ROS examples updated for compatibility with modern dependencies.

#####Modernization Summary (2025)

This version represents a comprehensive modernization effort:

**Standards Upgrade**: 
- C++11 → C++23 throughout the entire codebase
- All CMakeLists.txt files updated for modern CMake practices
- Compiler flags and build system modernized

**Library Updates**:
- OpenCV: Legacy headers → Modern OpenCV 4.6.0+ with opencv2/opencv.hpp
- Eigen3: Updated to 3.4.0+ for template compatibility
- Pangolin: Latest source build for C++23 support
- Boost: Added serialization dependencies

**Compatibility**:
- All third-party libraries (DBoW2, g2o, Sophus) updated for C++23
- Cross-platform testing on Ubuntu 24.04 LTS
- Backward compatibility maintained where possible

**Performance**:
- Benefits from modern OpenCV optimizations
- Improved build times with modern CMake
- Enhanced type safety with C++23 features

All licensing remains unchanged from original dependencies.




