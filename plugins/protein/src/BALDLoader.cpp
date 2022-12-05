/**
 * @author Nils Lichtenberg
 */

#include "BALDLoader.h"
#include "mmcore/Call.h"
#include "mmcore/param/FilePathParam.h"
#include "protein_calls/BALDCall.h"
#include <fstream>

using namespace megamol;
using namespace megamol::core;
using namespace megamol::protein;
using namespace megamol::protein_calls;

BALDLoader::BALDLoader()
        : AnimDataModule()
        , baldFilenameSlot("pdbFilename", "The path to the PDB data file to be loaded")
        , getDataSlot("dataOut", "Gets the data from the data source") {

    this->baldFilenameSlot << new param::FilePathParam(
        "", param::FilePathParam::FilePathFlags_::Flag_File_RestrictExtension, {"bald"});
    this->MakeSlotAvailable(&this->baldFilenameSlot);

    this->getDataSlot.SetCallback(
        BALDCall::ClassName(), BALDCall::FunctionName(BALDCall::CallForGetData), &BALDLoader::getData);
    this->getDataSlot.SetCallback(
        BALDCall::ClassName(), BALDCall::FunctionName(BALDCall::CallForGetExtent), &BALDLoader::getExtent);
    this->MakeSlotAvailable(&this->getDataSlot);
}

BALDLoader::~BALDLoader() = default;

bool BALDLoader::create() {
    // intentionally empty
    return true;
}

bool BALDLoader::getData(core::Call& call) {
    using megamol::core::utility::log::Log;
    auto* baldCall = dynamic_cast<BALDCall*>(&call);

    if (this->baldFilenameSlot.IsDirty()) {
        this->baldFilenameSlot.ResetDirty();

        // read header data
        if (!readHeaderData()) {
            return false;
        }

        baldCall->SetFrameCount(headerData.numFrames);

        AnimDataModule::setFrameCount(headerData.numFrames);
        this->initFrameCache(headerData.numFrames);
    }
    return true;
}

bool BALDLoader::getExtent(core::Call& call) {
    auto* baldCall = dynamic_cast<BALDCall*>(&call);
    if (baldCall == nullptr)
        return false;

    if (this->baldFilenameSlot.IsDirty()) {
        // read header data
        if (!readHeaderData()) {
            return false;
        }
    }

    baldCall->SetDataHash(1);
    baldCall->SetExtent(headerData.numFrames, -1.0f, -1.0f, -1.0f, 1.0f, 1.0f, 1.0f);

    return true;
}

void BALDLoader::release() {}

view::AnimDataModule::Frame* BALDLoader::constructFrame() const {
    auto* f = new Frame(*const_cast<BALDLoader*>(this));
    return f;
}

void BALDLoader::loadFrame(view::AnimDataModule::Frame* frame, unsigned int idx) {

    auto* fr = dynamic_cast<BALDLoader::Frame*>(frame);
    if (fr == nullptr)
        return;

    std::filesystem::path path = this->baldFilenameSlot.Param<core::param::FilePathParam>()->Value();
    std::ifstream baldFile;
    baldFile.open(path, std::ios::in | std::ios::binary);

    // skip header data
    baldFile.ignore(sizeof(int) * 3 + // magicNumber, numAtoms, numFrames
                    sizeof(float));   // sasRadius

    // skip data to the desired frame
    unsigned long long frameByteSize = sizeof(int)                           // maxLayer
                                       + sizeof(char) * headerData.numAtoms; // layer info per atom
    baldFile.ignore(static_cast<std::make_signed_t<long long>>(frameByteSize) * idx);

    unsigned int maxLayer;
    baldFile.read(reinterpret_cast<char*>(&maxLayer), sizeof(maxLayer));
    fr->setMaxLayer(maxLayer);
    fr->setFrameIdx(idx);

    if (fr->getMaxLayer() != -1) {
        std::vector<char> layers(headerData.numAtoms);
        baldFile.read((char*)layers.data(), headerData.numAtoms);
        fr->setLayers(std::vector<float>(layers.begin(), layers.end()));
    }
}

bool BALDLoader::readHeaderData() {
    std::filesystem::path path = this->baldFilenameSlot.Param<core::param::FilePathParam>()->Value();
    std::ifstream baldFile;
    baldFile.open(path, std::ios::in | std::ios::binary);

    baldFile.read(reinterpret_cast<char*>(&headerData.magicNumber), sizeof(headerData.magicNumber));
    baldFile.read(reinterpret_cast<char*>(&headerData.numAtoms), sizeof(headerData.numAtoms));
    baldFile.read(reinterpret_cast<char*>(&headerData.numFrames), sizeof(headerData.numFrames));
    baldFile.read(reinterpret_cast<char*>(&headerData.sasRadius), sizeof(headerData.sasRadius));

    if (headerData.magicNumber != 1111575620 // := 42 41 4C 44 := "BALD" in Big Endian
        || headerData.numAtoms <= 0 || headerData.numFrames <= 0 || headerData.sasRadius <= 0) {
        megamol::core::utility::log::Log::DefaultLog.WriteError("[importLayerData:]: Header invalid!");
        return false;
    }

    return true;
}

megamol::protein::BALDLoader::Frame::Frame(megamol::core::view::AnimDataModule& owner)
        : view::AnimDataModule::Frame(owner)
        , maxLayer(0) {
    // Intentionally empty
}

megamol::protein::BALDLoader::Frame::~Frame() = default;
