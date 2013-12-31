
# Script to deal with the absence of quoting around library
# specifications for Mac OS. These use the style '-framework XXX'
# which is produced unquoted by cmake. If we want to use them in a
# second project, then the linker misinterprets them as corrupt
# versions of the -lYYY style.

irtkBuildDir=/Users/paulaljabar/work/packages/irtk/build


if [[ -f IRTKInclude.cmake ]]; then
    cp IRTKInclude.cmake IRTKInclude.cmake.BAK
fi

cp $irtkBuildDir/lib/IRTKInclude.cmake temp.cmake

sed -i -e 's/;\-framework/ -framework/' temp.cmake

sed -i  -e 's/\-framework \([^ )]*\)/"-framework \1"/g' temp.cmake


cp temp.cmake IRTKInclude.cmake

