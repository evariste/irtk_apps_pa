#ifdef HAS_VTK

#include <irtkImage.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>

#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

void usage(){
  cerr << "polydatadeletearray [input] [output] <-name arrayName | -all>" << endl;
  cerr << "" << endl;
  exit(1);
}

char *input_name = NULL;
char *output_name = NULL;
char *scalar_name = NULL;

int main(int argc, char **argv)
{

	int ok, deleteAll = False;

	if (argc < 4){
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;
  output_name  = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-name") == 0)){
      argc--;
      argv++;
      scalar_name  = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-all") == 0)){
      argc--;
      argv++;
      deleteAll = True;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  if (deleteAll == False && scalar_name == NULL){
  	cerr << "No array(s) specified." << endl;
  	usage();
  }

  if (deleteAll == True){
  	cout << "Deleting all arrays." << endl;
  	scalar_name = NULL;
  }


  // Read surface
  vtkPolyDataReader *surface_reader = vtkPolyDataReader::New();
  surface_reader->SetFileName(input_name);
  surface_reader->Modified();
  surface_reader->Update();

  vtkPolyData *surface = surface_reader->GetOutput();

  int i, noOfArrays;
  int success = False;

  noOfArrays = surface->GetPointData()->GetNumberOfArrays();

  if (noOfArrays < 1){
  	cout << "Surface has no arrays." << endl;
  	exit(1);
  }

  if (deleteAll == True){
  	while(surface->GetPointData()->GetNumberOfArrays() > 0){
  		vtkFloatArray *currArray;
  		currArray = (vtkFloatArray*) surface->GetPointData()->GetArray(0);
			cout << "Deleting array : " << currArray->GetName() << endl;
  		surface->GetPointData()->RemoveArray(currArray->GetName());
  		success = True;
  	}
  } else {
  	for (i = 0; i < noOfArrays; i++){
  		vtkFloatArray *currArray;
  		currArray = (vtkFloatArray*) surface->GetPointData()->GetArray(i);

  		if (strcmp(currArray->GetName(), scalar_name) == 0){
  			cout << "Deleting array : " << currArray->GetName() << endl;
  			surface->GetPointData()->RemoveArray(scalar_name);

  			surface->Update();
  			success = True;
  			break;
  		}
  	}
  }


  if (success == True){
    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    cout << "Writing surface to file: " << output_name << endl;
    writer->SetInput(surface);
    writer->SetFileName(output_name);
    writer->SetFileTypeToBinary();
    writer->Write();
  } else {
  	if (deleteAll == False){
  		cerr << "No such scalars : " << scalar_name << endl;
  	}
  }



}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
