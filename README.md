
### Building for Sphinx

This should be similar to the Geoflow build, the only difference is that you need to 
use the flag `-DENABLE_PBAM_SPHINX=ON` when you run cmake. Or you can ignore this and
read the following steps:

1. `cd prototype/plugins/PB_S_AM/src`
2. `mkdir build`
3. `cd build`
4. `cmake -DENABLE_PBAM_SPHINX=on ../`
5. `make`
6. `cp pbam_sph.so ../../`


### General XCode set-up: ###

#### Adding a new target ####
1. Editor --> Add target
2. OS X --> Application --> Command Line Tool 

#### Steps ####
1. Add path to user-made header files to Build Settings --> Search Paths --> User Header Search Paths for all Targets
2. Set Apple LLVM Language Dialects
    C Language Dialect    GNU99
    C++ Language Dialect    GNU++11
    C++ Standard Library    libc++
3. For Target pbsolvers
    Add source .cpp and .h Files to pbsolvers Project
    Under Build Phases --> Compile Sources : only need .cpp files
    


#### For the simple schemes: ####

Need to change the Search Paths : User Header Search Paths
to where your code is on the local machine!

### GTest Setup ###
For gtest, it's a little more complicated...


1. Compile gtest.framework for version 1.7!!
  1. The first step is to download the latest version of gtest
  2. Go to the XCode directory and I had to do the following hacks
    1. Comment out the following lines in gtest-port.h:
         #include "gtest/internal/gtest-port-arch.h"
         #include "gtest/internal/custom/gtest-port.h"
    2. Comment out the following options in xcode/Config/General.xconfig:
        `SDKROOT`, `MACOS_DEPLOYMENT_TARGET`, and `GCC_VERSION`
  3. With those changes the framework should compile
  4. To find it, right click on the gtest.framework icon and select "Show in Finder"
  5. Copy this to a directory of your choice

2. Add Target gtests

3. Link to gtest.framework
    Add path to gtest.framework to Target gtests --> Build Settings --> Search Paths --> Framework Search Paths
    Under Build Phases --> Link with Binary Libraries : add gtest.framwork 
    Under Product --> Scheme --> Edit Scheme --> Run --> Environment Variables : Add DYLD_FRAMEWORK_PATH and path to gtest.framework

4. Files with Google Test format are called ClassUnitTest.h for each class.
    The Target Membership for these is to Test only.

5. The .cpp files of the main source must have Target Membership to Test and to pbsolvers.
    The .cpp files of the main source should also appear in Build Phases --> Compile Sources for Target gtests



#### Other XCode hints : ####
- To see page width:
    XCode --> Preferences --> Text Editing --> Show : Page guide at column 80

- To increase font size :
    XCode --> Preferences --> Fonts and Colors 
    Duplicate one of the styles. Select all text in example box. Change font size. 

- To access Instruments :
    XCode --> Product --> Profile ; select profiling template
        For timings, select Counters
    make sure path to executable is correct -- XCode makes a Debug and Release version in separate folders
    Product --> Clean ; Product --> Build for --> Profiling ; Product --> Profile
        should display data in Instruments window 
