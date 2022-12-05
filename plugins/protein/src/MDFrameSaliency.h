/**
 * @author Nils Lichtenberg
 */
#ifndef MMPROTEINPLUGIN_SALIENTFRAME_H_INCLUDED
#define MMPROTEINPLUGIN_SALIENTFRAME_H_INCLUDED
#if (defined(_MSC_VER) && (_MSC_VER > 1000))
#pragma once
#endif /* (defined(_MSC_VER) && (_MSC_VER > 1000)) */

#include "mmcore/Call.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"
#include "protein_calls/MolecularDataCall.h"

namespace megamol::protein {

class MDFrameSaliency : public core::Module {
public:
    /** Ctor */
    MDFrameSaliency();

    /** Dtor */
    ~MDFrameSaliency() override;

    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName() {
        return "MDFrameSaliency";
    }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description() {
        return "Analyses MD data an assigns one or more saliency measures to each frame. A subsequent filter should be used to find the most interesting frames in an MD simulation.";
    }

    /**
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable() {
        return true;
    }

protected:
    /**
     * Implementation of 'Create'.
     *
     * @return 'true' on success, 'false' otherwise.
     */
    bool create() override;

    /**
     * Implementation of 'Release'.
     */
    void release() override;

    /**
     * Call callback to get the data
     *
     * @param c The calling call
     * @return True on success
     */
    bool getData(core::Call& call);

    /**
     * Call callback to get the extents
     *
     * @param c The calling call
     * @return True on success
     */
    bool getExtent(core::Call& call);


    /**
     * Call callback to get the molecular data
     * that is being passed through.
     *
     * @param c The calling call
     *
     * @return True on success
     */
    bool getDataPassthroughMolecularData(core::Call& call);

    /**
     * Call callback to get the extent of the molecular data
     * that is being passed through.
     *
     * @param c The calling call
     *
     * @return True on success
     */
    bool getExtentPassthroughMolecularData(core::Call& call);

private:
    /** caller slot */
    core::CallerSlot inMolDataSlot;

    /** caller slot */
    core::CallerSlot inBALDSlot;

    /** callee slot */
    core::CalleeSlot outDataSlot;

    /** callee slot */
    core::CalleeSlot outMolDataSlot;
};

} // namespace megamol::protein

#endif /* MMPROTEINPLUGIN_SALIENTFRAME_H_INCLUDED */
