#rootcling -f PROfit_dict.cxx -c sbnanaobj/sbnanaobj/StandardRecord/*.h LinkDef.h
genreflex build/_deps/sbnanaobj-src/sbnanaobj/StandardRecord/classes.h \
    -s build/_deps/sbnanaobj-src/sbnanaobj/StandardRecord/classes_def.xml \
    -o PROfit_dict.cxx \
    --noIncludePaths \
    --interpreteronly \
    -Ibuild/_deps/sbnanaobj-src/
g++ PROfit_dict.cxx -o libPROfit_dict.so -shared -fPIC -I${ROOTSYS}/include -Ibuild/_deps/sbnanaobj-src/ `root-config --cflags --libs`
