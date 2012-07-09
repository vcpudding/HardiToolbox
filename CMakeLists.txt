cmake_minimum_required(VERSION 2.8)
PROJECT(hiermrf)

#if you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE OFF)

ADD_DEFINITIONS(-std=gnu++0x)

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

INCLUDE_DIRECTORIES (/home/collab/xhzhu/Libraries/arma-bin/usr/local/include)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR})

LINK_DIRECTORIES (/home/collab/xhzhu/Libraries/arma-bin/usr/local/lib64)
#LINK_DIRECTORIES (/home/sci/weiliu/packages/lemon_build/lemon)



ADD_EXECUTABLE(HardiToolbox
  main.cpp HardiToolbox.cpp testProc.cpp utils.cpp
  )

ADD_EXECUTABLE(RunMultiTensorTest
  testMultiTensor.cpp HardiToolbox.cpp
  )

ADD_EXECUTABLE(TestWeightEst
  tests/testWeightEst.cpp HardiToolbox.cpp utils.cpp
  )

ADD_EXECUTABLE(TestBrainData
  tests/testBrainData.cpp HardiToolbox.cpp utils.cpp testProc.cpp
  )

ADD_EXECUTABLE(TestPhantomData
  tests/testPhantomData.cpp HardiToolbox.cpp utils.cpp
  )

#need to link to some other libraries ? just add them here

TARGET_LINK_LIBRARIES(HardiToolbox libarmadillo.so ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES(TestWeightEst libarmadillo.so)
TARGET_LINK_LIBRARIES(TestBrainData libarmadillo.so ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES(TestPhantomData libarmadillo.so ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES(RunMultiTensorTest libarmadillo.so)
#TARGET_LINK_LIBRARIES(singlemrf emon ${ITK_LIBRARIES} boost_program_options boost_filesystem)
