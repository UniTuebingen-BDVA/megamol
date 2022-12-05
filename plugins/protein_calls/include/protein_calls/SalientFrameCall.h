/**
 * SalientFrameCall.h
 *
 * Copyright (C) 2010 by VISUS (University of Stuttgart)
 * Alle Rechte vorbehalten.
 *
 * @author Nils Lichtenberg
 */
#ifndef MEGAMOL_PROTEIN_SALIENTFRAMECALL_H_INCLUDED
#define MEGAMOL_PROTEIN_SALIENTFRAMECALL_H_INCLUDED
#if (defined(_MSC_VER) && (_MSC_VER > 1000))
#pragma once
#endif /* (defined(_MSC_VER) && (_MSC_VER > 1000)) */

#include "mmcore/Call.h"
#include "mmcore/factories/CallAutoDescription.h"
#include "mmstd/data/AbstractGetData3DCall.h"
#include "mmstd/data/AbstractGetDataCall.h"

namespace megamol {
namespace protein_calls {

/**
 * Get data describing MD frame saliency.
 */
class SalientFrameCall : public megamol::core::AbstractGetData3DCall {
public:
    /**
     * Answer the name of the objects of this description.
     *
     * @return The name of the objects of this description.
     */
    static const char* ClassName(void) {
        return "SalientFrameCall";
    }

    /**
     * Gets a human readable description of the module.
     *
     * @return A human readable description of the module.
     */
    static const char* Description(void) {
        return "Call to get information about MD frame saliency";
    }

    /**
     * Answer the number of functions used for this call.
     *
     * @return The number of functions used for this call.
     */
    static unsigned int FunctionCount(void) {
        return megamol::core::AbstractGetData3DCall::FunctionCount();
    }

    /**
     * Answer the name of the function used for this call.
     *
     * @param idx The index of the function to return it's name.
     *
     * @return The name of the requested function.
     */
    static const char* FunctionName(unsigned int idx) {
        return megamol::core::AbstractGetData3DCall::FunctionName(idx);
    }

    /** Ctor */
    SalientFrameCall(void);

    /** Dtor. */
    virtual ~SalientFrameCall(void);

    /** Index of the 'GetData' function */
    static const unsigned int CallForGetData = 0;
};

/** Description class typedef */
typedef megamol::core::factories::CallAutoDescription<SalientFrameCall> SalientFrameCallDescription;

} // namespace protein_calls
} /* end namespace megamol */

#endif /*  MEGAMOL_PROTEIN_SALIENTFRAMECALL_H_INCLUDED */
