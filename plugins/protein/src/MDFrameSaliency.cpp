/**
 * @author Nils Lichtenberg
 */
#include "MDFrameSaliency.h"
#include "mmcore/param/FloatParam.h"
#include "protein_calls/BALDCall.h"
#include "protein_calls/SalientFrameCall.h"

using namespace megamol;
using namespace megamol::core;
using namespace megamol::protein;
using namespace megamol::protein_calls;

MDFrameSaliency::MDFrameSaliency()
        : Module()
        , inMolDataSlot("molDataIn", "Molecular data source (usually PDBLoader)")
        , inBALDSlot("BALDIn", "Additional data provided by BALDLoader (optional)")
        , outDataSlot("dataOut", "The slot providing saliency information for each data frame")
        , outMolDataSlot("molDataOut", "Pass through of the molecular data") {

    this->inMolDataSlot.SetCompatibleCall<MolecularDataCallDescription>();
    this->inMolDataSlot.SetNecessity(AbstractCallSlotPresentation::Necessity::SLOT_REQUIRED);
    this->MakeSlotAvailable(&this->inMolDataSlot);

    this->inBALDSlot.SetCompatibleCall<BALDCallDescription>();
    this->inBALDSlot.SetNecessity(AbstractCallSlotPresentation::Necessity::SLOT_OPTIONAL);
    this->MakeSlotAvailable(&this->inBALDSlot);

    this->outMolDataSlot.SetCallback(MolecularDataCall::ClassName(),
        MolecularDataCall::FunctionName(MolecularDataCall::CallForGetData),
        &MDFrameSaliency::getDataPassthroughMolecularData);
    this->outMolDataSlot.SetCallback(MolecularDataCall::ClassName(),
        MolecularDataCall::FunctionName(MolecularDataCall::CallForGetExtent),
        &MDFrameSaliency::getExtentPassthroughMolecularData);
    this->MakeSlotAvailable(&this->outMolDataSlot);

    this->outDataSlot.SetCallback(
        SalientFrameCall::ClassName(), SalientFrameCall::FunctionName(0), &MDFrameSaliency::getData);
    this->outDataSlot.SetCallback(
        SalientFrameCall::ClassName(), SalientFrameCall::FunctionName(1), &MDFrameSaliency::getExtent);
    this->MakeSlotAvailable(&this->outDataSlot);
}

MDFrameSaliency::~MDFrameSaliency() {
    this->Release();
}

bool MDFrameSaliency::create() {
    return true;
}

bool MDFrameSaliency::getData(Call& call) {
    auto* baldCall = this->inBALDSlot.CallAs<protein_calls::BALDCall>();
    auto* molDataCall = this->inMolDataSlot.CallAs<MolecularDataCall>();

    if (molDataCall == nullptr) {
        return false;
    }

    if (!(*molDataCall)(MolecularDataCall::CallForGetData)) {
        return false;
    }

    if (!(*molDataCall)(MolecularDataCall::CallForGetExtent)) {
        return false;
    }

    if (baldCall != nullptr) {
        // execute the call
        if (!(*baldCall)(protein_calls::BALDCall::CallForGetData))
            return false;
    }

    if (baldCall != nullptr && molDataCall->FrameCount() + 1 != baldCall->FrameCount()) {
        // The number of frames provided by the molDataCall are not equal to the number of frames.
        // See Comment in PDBLoader.cpp:
        /*
         * //frames in xtc-file - 1 (without the last frame)
         * this->setFrameCount(this->numXTCFrames);
         */
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "Frame count of trajectory and BALD file do not match.");
        return false;
    }

    return true;
}

bool MDFrameSaliency::getExtent(Call& call) {
    auto* mol = this->inMolDataSlot.CallAs<MolecularDataCall>();
    if (mol == nullptr)
        return false;
    if (!(*mol)(MolecularDataCall::CallForGetExtent))
        return false;


    auto* baldCall = this->inBALDSlot.CallAs<protein_calls::BALDCall>();

    if (baldCall != nullptr) {
        // execute the call
        if (!(*baldCall)(protein_calls::BALDCall::CallForGetExtent))
            return false;
    }

    if (baldCall != nullptr && mol->FrameCount() + 1 != baldCall->FrameCount()) {
        // The number of frames provided by the molDataCall are not equal to the number of frames.
        // See Comment in PDBLoader.cpp:
        /*
         * //frames in xtc-file - 1 (without the last frame)
         * this->setFrameCount(this->numXTCFrames);
         */
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "Frame count of trajectory and BALD file do not match.");
        return false;
    }

    return true;
}

bool MDFrameSaliency::getDataPassthroughMolecularData(Call& call) {
    auto* outCall = dynamic_cast<MolecularDataCall*>(&call);
    if (outCall == nullptr) {
        return false;
    }

    auto* inCall = this->inMolDataSlot.CallAs<MolecularDataCall>();
    if (inCall == nullptr) {
        return false;
    }

    if (!(*inCall)(MolecularDataCall::CallForGetData)) {
        return false;
    }

    inCall->SetCalltime(outCall->Calltime());
    inCall->SetUnlocker(outCall->GetUnlocker());
    inCall->SetDataHash(outCall->DataHash());
    inCall->SetFrameID(outCall->FrameID());
    *outCall = *inCall;

    return true;
}

bool MDFrameSaliency::getExtentPassthroughMolecularData(Call& call) {
    auto* outCall = dynamic_cast<MolecularDataCall*>(&call);
    if (outCall == nullptr) {
        return false;
    }

    auto* inCall = this->inMolDataSlot.CallAs<MolecularDataCall>();
    if (inCall == nullptr) {
        return false;
    }

    if (!(*inCall)(MolecularDataCall::CallForGetExtent)) {
        return false;
    }

    inCall->SetCalltime(outCall->Calltime());
    inCall->SetUnlocker(outCall->GetUnlocker());
    inCall->SetDataHash(outCall->DataHash());
    inCall->SetFrameID(outCall->FrameID());
    *outCall = *inCall;

    return true;
}

void MDFrameSaliency::release() {}