#include <vtkVersion.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVRMLImporter.h>
#include <vtkDataSet.h>
#include <vtkActorCollection.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkActor.h>
#include <vtkSmartPointer.h>
#include <vtkRendererCollection.h>
#include <vtkLight.h>

int main ( int argc, char *argv[])
{
  if(argc != 2)
  {
    std::cout << "Required arguments: Filename.wrl" << std::endl;
    std::cout << "File should be a VRML 2.0 format file" << std::endl;
    return EXIT_FAILURE;
  }

  std::string filename = argv[1];
  std::cout << "Reading " << filename << std::endl;

  // TODO: Check that file starts with the string  '#VRML V2.0'


  vtkSmartPointer<vtkRenderer> ren1= vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(ren1);

  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);


  // VRML Import
  vtkSmartPointer<vtkVRMLImporter> importer = vtkSmartPointer<vtkVRMLImporter>::New();
  importer->SetFileName ( filename.c_str() );
  cout << "bla 1" << endl;
  importer->Read();
  cout << "bla" << endl;
  importer->SetRenderWindow(renderWindow);
  importer->Update();



  ren1=importer->GetRenderer();
  vtkRendererCollection *renCollection=vtkRendererCollection::New();
  renCollection = renderWindow->GetRenderers();
  renCollection->InitTraversal();
  vtkRenderer *ren= vtkRenderer::New();
  ren = renCollection->GetNextItem();

  ren->ResetCamera();
  vtkLight *light1=vtkLight::New();
  light1->SetLightTypeToCameraLight ();
  ren->AddLight(light1);

  vtkActorCollection *actorcol=vtkActorCollection::New();
  actorcol=ren->GetActors();

  vtkActor *actor=vtkActor::New();
  actor=actorcol->GetLastActor();

  vtkMapper *map=actor->GetMapper();
  vtkDataSet *PolyData=map->GetInput();

  ren->AddActor(actor);
  renderWindow->Render();
  renderWindowInteractor->Start();

  return 0;
}


