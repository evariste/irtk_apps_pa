SET(PA_INCLUDES
../include/abcdRandomRotationSource.h
../include/vtkPolyDataSmoothCustom.h
)

SET(PA_SRCS
abcdRandomRotationSource.cc
vtkPolyDataSmoothCustom.cc)

ADD_LIBRARY(pa++ ${PA_SRCS} ${PA_INCLUDES})

INSTALL_FILES(/include FILES ${PA_INCLUDES})
