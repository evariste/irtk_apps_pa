

#ifndef ABCDVTKUTILS_H
#define ABCDVTKUTILS_H

#include <string.h>
#include <sys/stat.h>

#include <irtkImage.h>


bool is_numeric(const char *str);

#ifdef HAS_VTK
bool is_vtkPolyDataFile(const char* filename);
#endif




#endif

