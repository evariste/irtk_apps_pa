set -vx

# Run like this: bash TESTS.sh  &> test-output.txt

# Then can compare test output with previously saved output in
# test-output-archive.txt


tar xf outputFilesArchive.tar  

for f in `ls outputFilesArchive`; do rm -f $f; done


polydatasphere bla.vtk -solid 0

polydatasphere bla2.vtk -solid 1

polydataappend  2 bla.vtk bla2.vtk  bla3.vtk

polydatabumps  bumps.vtk  -l 4  -randomRotation  -seed 23455

polydatacurvatureindices  bumps.vtk 


polydatadecimate bumps.vtk bumps-d.vtk

polydatalistscalars bumps.vtk 

polydatadeletescalars bumps.vtk bumps.vtk -name theta

polydatalistscalars bumps.vtk 



binarize tr-aal-to-neo-nr-40-crop.nii.gz  temp.nii.gz 1 1000  1000 0

mcubes temp.nii.gz temp.vtk 500

polydata2maskimage  temp.vtk temp.nii.gz  temp-recovered.nii.gz  -value 10

polydataassignlabels  tr-aal-to-neo-nr-40-crop.nii.gz  temp.vtk    temp-label.vtk -name aal



extract_roi tr-aal-to-neo-nr-40-crop.nii.gz  rois.nii.gz  10 20

rescale rois.nii.gz  rois.nii.gz  0 1000

dilation  rois.nii.gz rois.nii.gz 

mcubes  rois.nii.gz rois.vtk 500 

polydatalcc  rois.vtk  rois-lcc.vtk   


polydatamaths bumps.vtk -array_name phi -mul 5 -add 2 bumps-maths.vtk


polydatascalarstats  bumps.vtk -name phi


polydatascalarstats  bumps-maths.vtk -name phi



polydatarecalculatenormals bumps.vtk  bumps-recalcN.vtk



polydatareflect bumps.vtk  bumps-r.vtk -x 



polydatadecimate bumps.vtk bumps-dec.vtk  -reduction 0.7

polydataremesh bumps-dec.vtk bumps.vtk bumps-remesh.vtk  

# Two runs of random scalars starting with the same polydata.
polydatarandomscalars  bumps.vtk bumps-rs.vtk -seed 234

polydatarandomscalars  bumps.vtk bumps-rs-2.vtk -seed 456

# The stats should be different for each run of the above.
polydatascalarstats  bumps-rs.vtk -name Random

polydatascalarstats  bumps-rs-2.vtk -name Random


polydatascalarsmooth bumps-rs.vtk bumps-ss.vtk 20  2   -name Random



polydatascalarstats  bumps-rs.vtk -name Random


polydatascalarstats  bumps-ss.vtk -name Random



polydatasmooth  bumps.vtk  bumps-smooth.vtk   20 0.7


polydatarandomsphericalpoints ranSph-1.vtk -seed 234
polydatarandomsphericalpoints ranSph-2.vtk -seed 567

# Following should contain 1000 points at random on a unit sphere
# each. They should be different slightly, check by looking at centre
# of gravity.
polydatacentreofgravity ranSph-1.vtk 
polydatacentreofgravity ranSph-2.vtk 


# Following should give different outputs:
polydatageodesic bumps.vtk -reps 10 -seed 234
polydatageodesic bumps.vtk -reps 10 -seed 456


# Following should give different outputs:
polydataeuclidean bumps.vtk -reps 10 -seed 234
polydataeuclidean bumps.vtk -reps 10 -seed 456


# File comparisons

for f in `ls outputFilesArchive`
do
    ls -l $f | awk '{print $5, $9}'
    ls -l outputFilesArchive/$f | awk '{print $5, $9}'
    echo
done

