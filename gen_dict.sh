#rootcling -f PROfit_dict.cxx -c sbnanaobj/sbnanaobj/StandardRecord/*.h LinkDef.h
genreflex sbnanaobj/sbnanaobj/StandardRecord/classes.h \
    -s sbnanaobj/sbnanaobj/StandardRecord/classes_def.xml \
    -o PROfit_dict.cxx \
    --noIncludePaths \
    --interpreteronly \
    -Isbnanaobj/
g++ PROfit_dict.cxx -o libPROfit_dict.so -shared -fPIC -I${ROOTSYS}/include -Isbnanaobj `root-config --cflags --libs`
