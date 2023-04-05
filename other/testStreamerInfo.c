
#include <TFile.h>
#include <TStreamerInfo.h>

void PrintVersionNumbers(const char* filename, const char* classname)
{
    TFile* file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: failed to open file " << filename << std::endl;
        return;
    }

    TList* list = file->GetStreamerInfoList();
    TStreamerInfo* info = static_cast<TStreamerInfo*>(list->FindObject(classname));
    if (!info) {
        std::cerr << "Error: failed to find streamer info for class " << classname << std::endl;
        return;
    }

    std::cout << "Version number for class " << classname << " is " << info->GetClassVersion() << std::endl;

    file->Close();
    delete file;
}

void testStreamerInfo(){

    PrintVersionNumbers("~/Downloads/test.flat.caf.root", "caf::SRWeightParam");
    

}
