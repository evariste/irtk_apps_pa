SET(PA_INCLUDES
../include/abcdPointDistanceMap.h
../include/abcdRandomRotationSource.h
../include/abcdUtils.h
../include/vtkPolyDataSmoothCustom.h
)

SET(PA_SRCS
abcdPointDistanceMap.cc
abcdRandomRotationSource.cc
abcdUtils.cc
vtkPolyDataSmoothCustom.cc)

ADD_LIBRARY(pa++ ${PA_SRCS} ${PA_INCLUDES})

INSTALL_FILES(/include FILES ${PA_INCLUDES})
