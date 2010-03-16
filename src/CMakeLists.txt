SET(PA_INCLUDES
../include/vtkPolyDataSmoothCustom.h)

SET(PA_SRCS
vtkPolyDataSmoothCustom.cc)

ADD_LIBRARY(pa++ ${PA_SRCS} ${PA_INCLUDES})

INSTALL_FILES(/include FILES ${PA_INCLUDES})
