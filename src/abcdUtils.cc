
#include <abcdUtils.h>



bool is_numeric(const char *str)
{
  char *pEnd;
  double dummy;
  dummy = strtod(str, &pEnd);

  if (*pEnd != '\0')
    return false;

  return true;
}




bool is_vtkPolyDataFile(const char* filename) {
  // Rough and ready check to see if a file contains polydata.
  char *c = NULL;
  char buff[1000];
  struct stat statStruct;

  // Does file name end in .vtk?
  sprintf(buff, "%s", filename);
  c = strstr(buff, ".vtk\0");

  if (c == NULL) {
    return false;
  }

  // Does file exist?
  int i = stat(buff, &statStruct);
  if (i != 0) {
    return false;
  }

  // Have a go at reading it.
  vtkSmartPointer<vtkPolyDataReader> reader =
      vtkSmartPointer<vtkPolyDataReader>::New();
  reader->SetFileName(filename);
  if (!reader->IsFilePolyData()) {
    return false;
  }

  return true;
}
